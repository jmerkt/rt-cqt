/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

#include "../../include/resampling_filterbank.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rt_cqt
{

    template <int StageCount>
    class PythonResamplingFilterbank
    {
    public:
        void init(const double sample_rate, const int maximum_callback_block_size)
        {
            if (maximum_callback_block_size <= 0)
            {
                throw std::invalid_argument("The callback block size must be positive");
            }
            maximum_callback_block_size_ = maximum_callback_block_size;
            filterbank_.init(sample_rate, maximum_callback_block_size, std::max(4096, maximum_callback_block_size * 4));
        }

        std::pair<std::vector<std::vector<double>>, std::vector<double>> process(const std::vector<double> &input)
        {
            if (input.size() > static_cast<std::size_t>(maximum_callback_block_size_))
            {
                throw std::invalid_argument("The input exceeds the configured callback block size");
            }

            filterbank_.input_block(input.data(), static_cast<int>(input.size()));

            std::vector<std::vector<double>> stages(static_cast<std::size_t>(StageCount));
            for (int stage = 0; stage < StageCount; ++stage)
            {
                BufferPtr stage_input = filterbank_.get_stage_input_buffer(stage);
                const int stage_block_size = static_cast<int>(stage_input->get_write_read_distance());
                stages[static_cast<std::size_t>(stage)].resize(static_cast<std::size_t>(stage_block_size));
                if (stage_block_size == 0)
                {
                    continue;
                }

                stage_input->pull_block(stages[static_cast<std::size_t>(stage)].data(), stage_block_size);
                // Echo every stage into the synthesis side. This exposes the
                // complete filterbank path without involving a CQT.
                filterbank_.get_stage_output_buffer(stage)->push_block(stages[static_cast<std::size_t>(stage)].data(),
                                                                       stage_block_size);
            }

            std::vector<double> output(input.size(), 0.);
            const double *const output_data = filterbank_.output_block(static_cast<int>(input.size()));
            std::copy_n(output_data, output.size(), output.data());
            return {std::move(stages), std::move(output)};
        }

        int get_processing_block_size() const { return filterbank_.get_processing_block_size(); }

        int get_latency_samples() const { return filterbank_.get_latency_samples(); }

    private:
        int maximum_callback_block_size_{0};
        ResamplingFilterbank<StageCount> filterbank_;
    };

}
