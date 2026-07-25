/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

#include "../../include/ResamplingFilterbank.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Cqt
{

    template <int StageNumber>
    class Python_ResamplingFilterbank
    {
    public:
        void init(const double samplerate, const int maximumCallbackBlockSize)
        {
            if (maximumCallbackBlockSize <= 0)
            {
                throw std::invalid_argument("The callback block size must be positive");
            }
            mMaximumCallbackBlockSize = maximumCallbackBlockSize;
            mFilterbank.init(
                samplerate,
                maximumCallbackBlockSize,
                std::max(4096, maximumCallbackBlockSize * 4));
        }

        std::pair<std::vector<std::vector<double>>, std::vector<double>> processBlock(
            const std::vector<double> &input)
        {
            if (input.size() > static_cast<std::size_t>(mMaximumCallbackBlockSize))
            {
                throw std::invalid_argument("The input exceeds the configured callback block size");
            }

            mFilterbank.inputBlock(
                input.data(), static_cast<int>(input.size()));

            std::vector<std::vector<double>> stages(static_cast<std::size_t>(StageNumber));
            for (int stage = 0; stage < StageNumber; ++stage)
            {
                BufferPtr stageInput = mFilterbank.getStageInputBuffer(stage);
                const int stageBlockSize =
                    static_cast<int>(stageInput->getWriteReadDistance());
                stages[static_cast<std::size_t>(stage)].resize(
                    static_cast<std::size_t>(stageBlockSize));
                if (stageBlockSize == 0)
                {
                    continue;
                }

                stageInput->pullBlock(
                    stages[static_cast<std::size_t>(stage)].data(),
                    stageBlockSize);
                // Echo every stage into the synthesis side. This exposes the
                // complete filterbank path without involving a CQT.
                mFilterbank.getStageOutputBuffer(stage)->pushBlock(
                    stages[static_cast<std::size_t>(stage)].data(),
                    stageBlockSize);
            }

            std::vector<double> output(input.size(), 0.);
            const double *const outputBlock =
                mFilterbank.outputBlock(static_cast<int>(input.size()));
            std::copy_n(outputBlock, output.size(), output.data());
            return {std::move(stages), std::move(output)};
        }

        int getProcessingBlockSize() const
        {
            return mFilterbank.getProcessingBlockSize();
        }

        int getLatencySamples() const
        {
            return mFilterbank.getLatencySamples();
        }

    private:
        int mMaximumCallbackBlockSize{0};
        ResamplingFilterbank<StageNumber> mFilterbank;
    };

}
