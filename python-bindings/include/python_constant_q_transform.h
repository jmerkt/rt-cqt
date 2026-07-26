/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

#include "../../include/constant_q_transform.h"

#include <cstring>
#include <vector>

namespace Cqt
{
    template <int BinsPerOctave, int OctaveCount>
    class PythonConstantQTransform
    {
    public:
        void init(const int hop_size) { transform_.init(hop_size); }

        void init_sample_rate(const double sample_rate, const int block_size)
        {
            transform_.init_sample_rate(sample_rate, block_size);
        }

        void input_block(std::vector<double> &data)
        {
            transform_.input_block(data.data(), static_cast<int>(data.size()));
        }

        std::vector<double> output_block(const int block_size)
        {
            std::vector<double> output(block_size, 0.);
            const double *const output_data = transform_.output_block(block_size);
            std::memcpy(output.data(), output_data, block_size * sizeof(double));
            return output;
        }

        std::vector<ScheduleElement> &get_cqt_schedule() { return transform_.get_cqt_schedule(); }

        void cqt(const ScheduleElement schedule) { transform_.cqt(schedule); }

        void icqt(const ScheduleElement schedule) { transform_.icqt(schedule); }

        CqtBufferType *get_octave_cqt_buffer(const int octave) { return transform_.get_octave_cqt_buffer(octave); }

    private:
        ConstantQTransform<BinsPerOctave, OctaveCount> transform_;
    };
}
