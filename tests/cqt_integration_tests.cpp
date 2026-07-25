#include "ConstantQTransform.h"
#include "SlidingCqt.h"

#include <algorithm>
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

        for (const Cqt::ScheduleElement &element : schedule)
        {
            require(
                element.sample() >= 0 && element.sample() < 256,
                "Schedule position is outside the internal processing block");
            cqt.cqt(element);
            cqt.icqt(element);
        }

        double *output = cqt.outputBlock(BlockSize);
        requireFinite(output, BlockSize, "Constant CQT produced a non-finite sample");
    }
}

}

int main()
{
    try
    {
        testSlidingCqtBatchedInput();
        std::cout << "[pass] sliding CQT batched input\n";
        testConstantCqtBatchedSchedule();
        std::cout << "[pass] constant CQT batched schedule\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << "[fail] " << error.what() << '\n';
        return 1;
    }

    return 0;
}
