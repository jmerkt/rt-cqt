#include "ConstantQTransform.h"
#include "SlidingCqt.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void require(const bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void requireFinite(const double *const data, const int size, const std::string &message)
{
    for (int sample = 0; sample < size; ++sample)
    {
        require(std::isfinite(data[sample]), message);
    }
}

std::vector<double> makeSineBlock(
    const int blockSize,
    const int blockIndex,
    const double samplerate,
    const double frequency)
{
    std::vector<double> block(static_cast<std::size_t>(blockSize));
    for (int sample = 0; sample < blockSize; ++sample)
    {
        const int position = blockIndex * blockSize + sample;
        block[static_cast<std::size_t>(sample)] =
            std::sin(2. * audio_utils::Pi<double>() * frequency * position / samplerate);
    }
    return block;
}

double measureToneAmplitude(
    const std::vector<double> &data,
    const std::size_t start,
    const double samplerate,
    const double frequency)
{
    double real = 0.;
    double imaginary = 0.;
    for (std::size_t sample = start; sample < data.size(); ++sample)
    {
        const double phase =
            2. * audio_utils::Pi<double>() * frequency *
            static_cast<double>(sample) / samplerate;
        real += data[sample] * std::cos(phase);
        imaginary -= data[sample] * std::sin(phase);
    }
    const double sampleCount = static_cast<double>(data.size() - start);
    return 2. * std::hypot(real, imaginary) / sampleCount;
}

void testSlidingCqtBatchedInput()
{
    constexpr int BlockSize = 64;
    Cqt::SlidingCqt<12, 9, false> cqt;
    cqt.init(48000., BlockSize);

    for (int blockIndex = 0; blockIndex < 12; ++blockIndex)
    {
        std::vector<double> input = makeSineBlock(BlockSize, blockIndex, 48000., 440.);
        cqt.inputBlock(input.data(), BlockSize);

        const bool processingBoundary = ((blockIndex + 1) % 4) == 0;
        for (int stage = 0; stage < 9; ++stage)
        {
            const std::size_t expectedSamples =
                processingBoundary ? static_cast<std::size_t>(256 >> stage) : 0U;
            require(
                cqt.getSamplesToProcess(stage) == expectedSamples,
                "Sliding CQT received an unexpected stage block size");
        }

        double *output = cqt.outputBlock(BlockSize);
        requireFinite(output, BlockSize, "Sliding CQT produced a non-finite sample");
    }
}

template <bool Windowing>
void testSlidingCqtAmplitudeNormalization()
{
    constexpr int BlockSize = 64;
    constexpr int BlockCount = 512;
    constexpr double Samplerate = 48000.;
    constexpr int Tone = 8;
    const double frequency = Cqt::computeBinFrequency(
        Cqt::computeReferenceFrequency(440.), 24, 0, Tone);

    Cqt::SlidingCqt<24, 9, Windowing> cqt;
    cqt.init(Samplerate, BlockSize);

    std::vector<double> output(
        static_cast<std::size_t>(BlockSize * BlockCount), 0.);
    double magnitudeSum = 0.;
    int magnitudeCount = 0;
    for (int blockIndex = 0; blockIndex < BlockCount; ++blockIndex)
    {
        std::vector<double> input =
            makeSineBlock(BlockSize, blockIndex, Samplerate, frequency);
        cqt.inputBlock(input.data(), BlockSize);

        if (blockIndex >= (BlockCount / 2) && cqt.getSamplesToProcess(0) > 0)
        {
            magnitudeSum += std::abs(
                cqt.getOctaveCqtBuffer(0)[Tone].pullDelaySample(0));
            ++magnitudeCount;
        }

        const double *const outputBlock = cqt.outputBlock(BlockSize);
        std::copy_n(
            outputBlock,
            BlockSize,
            output.data() + static_cast<std::size_t>(blockIndex * BlockSize));
    }

    const double coefficientMagnitude =
        magnitudeSum / static_cast<double>(magnitudeCount);
    require(
        std::abs(coefficientMagnitude - 0.5) < 1.e-3,
        "Sliding CQT coefficient normalization is incorrect: " +
            std::to_string(coefficientMagnitude));

    const double reconstructedAmplitude =
        measureToneAmplitude(output, output.size() / 2, Samplerate, frequency);
    const double minimumAmplitude = Windowing ? 0.9 : 0.75;
    require(
        reconstructedAmplitude > minimumAmplitude &&
            reconstructedAmplitude < 1.1,
        "Sliding CQT resynthesis normalization is incorrect: " +
            std::to_string(reconstructedAmplitude));
}

void testConstantCqtBatchedSchedule()
{
    constexpr int BlockSize = 64;
    Cqt::ConstantQTransform<12, 9> cqt;
    cqt.init(64);
    cqt.initFs(48000., BlockSize);

    for (int blockIndex = 0; blockIndex < 12; ++blockIndex)
    {
        std::vector<double> input = makeSineBlock(BlockSize, blockIndex, 48000., 440.);
        cqt.inputBlock(input.data(), BlockSize);

        const auto &schedule = cqt.getCqtSchedule();
        const bool processingBoundary = ((blockIndex + 1) % 4) == 0;
        require(
            processingBoundary ? !schedule.empty() : schedule.empty(),
            "Constant CQT schedule is not aligned with filterbank processing");

        std::array<int, 9> firstDelay;
        std::array<int, 9> previousSynthesisOffset;
        firstDelay.fill(-1);
        previousSynthesisOffset.fill(-1);
        for (const Cqt::ScheduleElement &element : schedule)
        {
            require(
                element.sample() >= 0 && element.sample() < 256,
                "Schedule position is outside the internal processing block");
            const int octave = element.octave();
            if (firstDelay[octave] < 0)
            {
                firstDelay[octave] = element.delayOctaveRate();
            }
            require(
                element.synthesisOffset() ==
                    firstDelay[octave] - element.delayOctaveRate(),
                "Synthesis offset is inconsistent with the analysis delay");
            require(
                element.synthesisOffset() > previousSynthesisOffset[octave],
                "Synthesis offsets must advance within an internal block");
            previousSynthesisOffset[octave] = element.synthesisOffset();
            cqt.cqt(element);
            cqt.icqt(element);
        }

        double *output = cqt.outputBlock(BlockSize);
        requireFinite(output, BlockSize, "Constant CQT produced a non-finite sample");
    }
}

void testConstantCqtHighFrequencyResynthesis()
{
    constexpr int BlockSize = 64;
    constexpr int BlockCount = 512;
    constexpr double Samplerate = 48000.;
    const double frequency = Cqt::computeBinFrequency(
        Cqt::computeReferenceFrequency(440.), 24, 0, 8);

    Cqt::ConstantQTransform<24, 9> cqt;
    cqt.init(64);
    cqt.initFs(Samplerate, BlockSize);

    std::vector<double> output(
        static_cast<std::size_t>(BlockSize * BlockCount), 0.);
    for (int blockIndex = 0; blockIndex < BlockCount; ++blockIndex)
    {
        std::vector<double> input =
            makeSineBlock(BlockSize, blockIndex, Samplerate, frequency);
        cqt.inputBlock(input.data(), BlockSize);
        for (const Cqt::ScheduleElement &element : cqt.getCqtSchedule())
        {
            cqt.cqt(element);
            cqt.icqt(element);
        }

        const double *const outputBlock = cqt.outputBlock(BlockSize);
        std::copy_n(
            outputBlock,
            BlockSize,
            output.data() + static_cast<std::size_t>(blockIndex * BlockSize));
    }

    const double reconstructedAmplitude =
        measureToneAmplitude(output, output.size() / 2, Samplerate, frequency);
    require(
        reconstructedAmplitude > 0.5,
        "Constant CQT lost high-frequency resynthesis amplitude: " +
            std::to_string(reconstructedAmplitude));
}

}

int main()
{
    try
    {
        testSlidingCqtBatchedInput();
        std::cout << "[pass] sliding CQT batched input\n";
        testSlidingCqtAmplitudeNormalization<false>();
        std::cout << "[pass] rectangular sliding CQT amplitude normalization\n";
        testSlidingCqtAmplitudeNormalization<true>();
        std::cout << "[pass] windowed sliding CQT amplitude normalization\n";
        testConstantCqtBatchedSchedule();
        std::cout << "[pass] constant CQT batched schedule\n";
        testConstantCqtHighFrequencyResynthesis();
        std::cout << "[pass] constant CQT high-frequency resynthesis\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << "[fail] " << error.what() << '\n';
        return 1;
    }

    return 0;
}
