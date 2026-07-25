#include "Resampling.h"
#include "ResamplingFilterbank.h"
#include "../python-bindings/include/Python_ResamplingFilterbank.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr unsigned TestAllpassNumber{3};
constexpr double TestTransitionBandwidth{0.1};

void require(const bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void requireNear(
    const std::vector<double> &actual,
    const std::vector<double> &expected,
    const double tolerance,
    const std::string &message)
{
    require(actual.size() == expected.size(), message + ": vector sizes differ");
    double maximumError = 0.;
    for (std::size_t index = 0; index < actual.size(); ++index)
    {
        maximumError = std::max(maximumError, std::abs(actual[index] - expected[index]));
    }
    require(maximumError <= tolerance, message + ": maximum error is " + std::to_string(maximumError));
}

std::vector<double> makeNoise(const std::size_t size)
{
    std::mt19937 generator(0x5eedU);
    std::uniform_real_distribution<double> distribution(-1., 1.);
    std::vector<double> result(size);
    for (double &sample : result)
    {
        sample = distribution(generator);
    }
    return result;
}

void testHalfBandPartitionInvariance()
{
    const std::vector<double> input = makeNoise(1024);

    Cqt::HalfBandLowpass<double, TestAllpassNumber> wholeDownsampler;
    wholeDownsampler.init(1024, true, TestTransitionBandwidth);
    double *wholeDownData = wholeDownsampler.processDown(input.data());
    const std::vector<double> wholeDown(wholeDownData, wholeDownData + 512);

    Cqt::HalfBandLowpass<double, TestAllpassNumber> partitionedDownsampler;
    partitionedDownsampler.init(256, true, TestTransitionBandwidth);
    std::vector<double> partitionedDown;
    partitionedDown.reserve(512);
    for (std::size_t offset = 0; offset < input.size(); offset += 256)
    {
        double *block = partitionedDownsampler.processDown(input.data() + offset);
        partitionedDown.insert(partitionedDown.end(), block, block + 128);
    }
    requireNear(partitionedDown, wholeDown, 1.e-14, "Half-band downsampling depends on block partition");

    Cqt::HalfBandLowpass<double, TestAllpassNumber> wholeUpsampler;
    wholeUpsampler.init(512, false, TestTransitionBandwidth);
    double *wholeUpData = wholeUpsampler.processUp(wholeDown.data());
    const std::vector<double> wholeUp(wholeUpData, wholeUpData + 1024);

    Cqt::HalfBandLowpass<double, TestAllpassNumber> partitionedUpsampler;
    partitionedUpsampler.init(128, false, TestTransitionBandwidth);
    std::vector<double> partitionedUp;
    partitionedUp.reserve(1024);
    for (std::size_t offset = 0; offset < wholeDown.size(); offset += 128)
    {
        double *block = partitionedUpsampler.processUp(wholeDown.data() + offset);
        partitionedUp.insert(partitionedUp.end(), block, block + 256);
    }
    requireNear(partitionedUp, wholeUp, 1.e-14, "Half-band upsampling depends on block partition");

    partitionedDownsampler.init(256, true, TestTransitionBandwidth);
    double *reinitializedData = partitionedDownsampler.processDown(input.data());
    Cqt::HalfBandLowpass<double, TestAllpassNumber> freshDownsampler;
    freshDownsampler.init(256, true, TestTransitionBandwidth);
    double *freshData = freshDownsampler.processDown(input.data());
    requireNear(
        std::vector<double>(reinitializedData, reinitializedData + 128),
        std::vector<double>(freshData, freshData + 128),
        1.e-14,
        "Reinitializing a half-band filter did not reset its state");
}

void testResamplingHandlerPartitionInvariance()
{
    const std::vector<double> input = makeNoise(1024);

    Cqt::ResamplingHandler<double, TestAllpassNumber> whole;
    whole.init(3, 1024, Cqt::DirectionConfig::Down);
    double *wholeData = whole.processDown(input.data());
    const std::vector<double> expected(wholeData, wholeData + 128);

    Cqt::ResamplingHandler<double, TestAllpassNumber> partitioned;
    partitioned.init(3, 256, Cqt::DirectionConfig::Down);
    std::vector<double> actual;
    actual.reserve(128);
    for (std::size_t offset = 0; offset < input.size(); offset += 256)
    {
        double *block = partitioned.processDown(input.data() + offset);
        actual.insert(actual.end(), block, block + 32);
    }

    requireNear(actual, expected, 1.e-14, "Power-of-two resampling depends on block partition");

    Cqt::ResamplingHandler<double, TestAllpassNumber> wholeRoundTrip;
    wholeRoundTrip.init(3, 1024, Cqt::DirectionConfig::DownUp);
    double *wholeDown = wholeRoundTrip.processDown(input.data());
    double *wholeUp = wholeRoundTrip.processUp(wholeDown);
    const std::vector<double> expectedRoundTrip(wholeUp, wholeUp + 1024);

    Cqt::ResamplingHandler<double, TestAllpassNumber> partitionedRoundTrip;
    partitionedRoundTrip.init(3, 256, Cqt::DirectionConfig::DownUp);
    std::vector<double> actualRoundTrip;
    actualRoundTrip.reserve(1024);
    for (std::size_t offset = 0; offset < input.size(); offset += 256)
    {
        double *down = partitionedRoundTrip.processDown(input.data() + offset);
        double *up = partitionedRoundTrip.processUp(down);
        actualRoundTrip.insert(actualRoundTrip.end(), up, up + 256);
    }
    requireNear(
        actualRoundTrip,
        expectedRoundTrip,
        1.e-14,
        "Power-of-two round-trip resampling depends on block partition");
}

template <int StageNumber>
struct FilterbankRun
{
    std::vector<std::vector<double>> stages{
        static_cast<std::size_t>(StageNumber)};
    std::vector<double> output;
};

template <int StageNumber>
FilterbankRun<StageNumber> runFilterbank(
    const double samplerate,
    const int maximumCallbackSize,
    const std::vector<double> &input,
    const std::vector<int> &partitions,
    const bool synthesize)
{
    Cqt::ResamplingFilterbank<StageNumber> filterbank;
    filterbank.init(samplerate, maximumCallbackSize, static_cast<int>(input.size() * 4));

    FilterbankRun<StageNumber> result;
    result.output.reserve(input.size());

    std::size_t inputPosition = 0;
    std::size_t partitionIndex = 0;
    while (inputPosition < input.size())
    {
        const int requestedSize = partitions[partitionIndex % partitions.size()];
        const int blockSize = std::min<int>(
            requestedSize, static_cast<int>(input.size() - inputPosition));
        require(blockSize <= maximumCallbackSize, "Test partition exceeds configured callback size");

        filterbank.inputBlock(
            input.data() + inputPosition, blockSize);

        for (int stage = 0; stage < StageNumber; ++stage)
        {
            Cqt::BufferPtr inputBuffer = filterbank.getStageInputBuffer(stage);
            const int stageBlockSize = static_cast<int>(inputBuffer->getWriteReadDistance());
            if (stageBlockSize == 0)
            {
                continue;
            }

            std::vector<double> stageData(static_cast<std::size_t>(stageBlockSize));
            inputBuffer->pullBlock(stageData.data(), stageBlockSize);
            result.stages[static_cast<std::size_t>(stage)].insert(
                result.stages[static_cast<std::size_t>(stage)].end(),
                stageData.begin(),
                stageData.end());

            if (synthesize)
            {
                filterbank.getStageOutputBuffer(stage)->pushBlock(
                    stageData.data(), stageBlockSize);
            }
        }

        if (synthesize)
        {
            double *outputBlock = filterbank.outputBlock(blockSize);
            result.output.insert(result.output.end(), outputBlock, outputBlock + blockSize);
        }

        inputPosition += static_cast<std::size_t>(blockSize);
        ++partitionIndex;
    }

    return result;
}

void testFilterbankBlockSizing()
{
    Cqt::ResamplingFilterbank<1> singleStage;
    singleStage.init(48000., 3, 16);
    require(singleStage.getProcessingBlockSize() == 3, "Single-stage filterbank was unnecessarily aligned");
    require(singleStage.getOriginBlockSize() == 3, "Unexpected single-stage origin block");
    const auto singleStageResult = runFilterbank<1>(
        48000., 3, makeNoise(12), {1, 2, 3}, true);
    require(singleStageResult.stages[0].size() == 12, "Single-stage analysis lost samples");
    require(singleStageResult.output.size() == 12, "Single-stage synthesis lost samples");

    Cqt::ResamplingFilterbank<9> filterbank48;
    filterbank48.init(48000., 64, 1024);
    require(filterbank48.getProcessingBlockSize() == 256, "Unexpected 48 kHz processing block");
    require(filterbank48.getOriginBlockSize() == 256, "Unexpected 48 kHz origin block");
    require(filterbank48.getLatencySamples() == 256, "Unexpected 48 kHz latency");

    Cqt::ResamplingFilterbank<9> filterbank44;
    filterbank44.init(44100., 64, 1024);
    require(filterbank44.getProcessingBlockSize() == 256, "Unexpected 44.1 kHz processing block");
    require(filterbank44.getOriginBlockSize() == 256, "Unexpected 44.1 kHz origin block");

    Cqt::ResamplingFilterbank<9> filterbank96;
    filterbank96.init(96000., 64, 1024);
    require(filterbank96.getProcessingBlockSize() == 512, "Unexpected 96 kHz processing block");
    require(filterbank96.getOriginBlockSize() == 256, "Unexpected 96 kHz origin block");

    Cqt::ResamplingFilterbank<9> filterbank88;
    filterbank88.init(88200., 64, 1024);
    require(filterbank88.getProcessingBlockSize() == 512, "Unexpected 88.2 kHz processing block");
    require(filterbank88.getOriginBlockSize() == 256, "Unexpected 88.2 kHz origin block");

    Cqt::ResamplingFilterbank<9> largeCallbackFilterbank;
    largeCallbackFilterbank.init(48000., 1000, 4096);
    require(largeCallbackFilterbank.getProcessingBlockSize() == 1024, "Callback was not aligned upward");

    bool rejectedUnsupportedRate = false;
    try
    {
        Cqt::ResamplingFilterbank<9> unsupported;
        unsupported.init(88201., 64, 1024);
    }
    catch (const std::invalid_argument &)
    {
        rejectedUnsupportedRate = true;
    }
    require(rejectedUnsupportedRate, "Unsupported sample rate was accepted");
}

void testFilterbankBatching()
{
    Cqt::ResamplingFilterbank<9> filterbank;
    filterbank.init(48000., 64, 1024);
    std::vector<double> block(64, 1.);

    for (int callback = 0; callback < 3; ++callback)
    {
        filterbank.inputBlock(block.data(), static_cast<int>(block.size()));
        require(filterbank.getLastProcessedInputSize() == 0, "Filterbank released a partial internal block");
    }

    filterbank.inputBlock(block.data(), static_cast<int>(block.size()));
    require(filterbank.getLastProcessedInputSize() == 256, "Filterbank did not release a complete internal block");
    for (int stage = 0; stage < 9; ++stage)
    {
        const std::size_t expectedSize = static_cast<std::size_t>(256 >> stage);
        require(
            filterbank.getStageInputBuffer(stage)->getWriteReadDistance() == expectedSize,
            "Unexpected stage block size at stage " + std::to_string(stage));
    }
}

void testCallbackSizeMatrix()
{
    const std::vector<std::pair<int, int>> callbackAndProcessingSizes{
        {1, 256},
        {3, 256},
        {64, 256},
        {257, 512},
        {1024, 1024},
    };

    for (const auto &[callbackSize, processingSize] : callbackAndProcessingSizes)
    {
        Cqt::ResamplingFilterbank<9> filterbank;
        filterbank.init(48000., callbackSize, processingSize * 4);
        require(
            filterbank.getProcessingBlockSize() == processingSize,
            "Unexpected processing size for callback size " + std::to_string(callbackSize));

        const std::vector<double> input = makeNoise(
            static_cast<std::size_t>(processingSize * 2));
        const auto result = runFilterbank<9>(
            48000., callbackSize, input, {callbackSize}, true);
        for (int stage = 0; stage < 9; ++stage)
        {
            require(
                result.stages[static_cast<std::size_t>(stage)].size() ==
                    (input.size() >> stage),
                "Callback-size matrix lost stage samples");
        }
        require(
            std::all_of(
                result.output.begin(),
                result.output.end(),
                [](const double sample) { return std::isfinite(sample); }),
            "Callback-size matrix produced a non-finite output");
    }
}

void testFilterbankPartitionInvariance()
{
    const std::vector<double> input = makeNoise(512);
    const auto regular = runFilterbank<5>(48000., 13, input, {13}, true);
    const auto irregular = runFilterbank<5>(48000., 13, input, {1, 7, 3, 11, 5, 2, 13}, true);

    for (int stage = 0; stage < 5; ++stage)
    {
        const std::size_t expectedSize = input.size() >> stage;
        require(
            regular.stages[static_cast<std::size_t>(stage)].size() == expectedSize,
            "Unexpected accumulated stage size");
        requireNear(
            irregular.stages[static_cast<std::size_t>(stage)],
            regular.stages[static_cast<std::size_t>(stage)],
            1.e-14,
            "Filterbank analysis depends on callback partition");
    }

    requireNear(irregular.output, regular.output, 1.e-14, "Filterbank synthesis depends on callback partition");
    require(
        std::all_of(
            regular.output.begin(),
            regular.output.begin() + 16,
            [](const double sample) { return sample == 0.; }),
        "The declared one-block filterbank latency is missing");
    require(
        std::any_of(
            regular.output.begin() + 16,
            regular.output.end(),
            [](const double sample) { return std::abs(sample) > 1.e-12; }),
        "Filterbank synthesis produced only silence");
}

void testFilterbankAt96kHz()
{
    const std::vector<double> input = makeNoise(512);
    const auto regular = runFilterbank<5>(96000., 13, input, {13}, true);
    const auto irregular = runFilterbank<5>(96000., 13, input, {13, 5, 1, 9}, true);
    for (int stage = 0; stage < 5; ++stage)
    {
        const std::size_t expectedSize = input.size() >> (stage + 1);
        require(
            regular.stages[static_cast<std::size_t>(stage)].size() == expectedSize,
            "Unexpected 96 kHz stage size");
        requireNear(
            irregular.stages[static_cast<std::size_t>(stage)],
            regular.stages[static_cast<std::size_t>(stage)],
            1.e-14,
            "96 kHz filterbank analysis depends on callback partition");
    }
    requireNear(
        irregular.output,
        regular.output,
        1.e-14,
        "96 kHz filterbank synthesis depends on callback partition");
}

void testPythonPlotAdapter()
{
    Cqt::Python_ResamplingFilterbank<9> filterbank;
    filterbank.init(48000., 64);
    std::vector<double> block(64, 0.);
    block[0] = 1.;

    for (int callback = 0; callback < 4; ++callback)
    {
        const auto result = filterbank.process(block);
        require(result.second.size() == block.size(), "Python adapter returned the wrong output size");
        const bool processingBoundary = callback == 3;
        require(
            result.first[0].size() == (processingBoundary ? 256U : 0U),
            "Python adapter exposed the wrong stage block size");
        std::fill(block.begin(), block.end(), 0.);
    }
}

void run(const std::string &name, const std::function<void()> &test)
{
    test();
    std::cout << "[pass] " << name << '\n';
}

}

int main()
{
    try
    {
        run("half-band partition invariance", testHalfBandPartitionInvariance);
        run("resampling-handler partition invariance", testResamplingHandlerPartitionInvariance);
        run("filterbank block sizing", testFilterbankBlockSizing);
        run("filterbank batching", testFilterbankBatching);
        run("callback-size matrix", testCallbackSizeMatrix);
        run("filterbank partition invariance", testFilterbankPartitionInvariance);
        run("filterbank at 96 kHz", testFilterbankAt96kHz);
        run("Python plot adapter", testPythonPlotAdapter);
    }
    catch (const std::exception &error)
    {
        std::cerr << "[fail] " << error.what() << '\n';
        return 1;
    }

    return 0;
}
