/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

#include "../../include/sliding_cqt.h"

#include <complex>
#include <cstring>
#include <vector>

namespace Cqt
{
    template <int BinsPerOctave, int OctaveCount, bool Windowing>
    class PythonSlidingCqt
    {
    public:
        void init(const double sample_rate, const int block_size) { transform_.init(sample_rate, block_size); }

        void input_block(std::vector<double> &data, const int block_size)
        {
            transform_.input_block(data.data(), block_size);
        }

        std::vector<double> output_block(const int block_size)
        {
            std::vector<double> output(block_size, 0.);
            const double *const output_data = transform_.output_block(block_size);
            std::memcpy(output.data(), output_data, block_size * sizeof(double));
            return output;
        }

        std::vector<std::complex<double>> get_octave_values(const int octave)
        {
            std::vector<std::complex<double>> values(BinsPerOctave, {0., 0.});
            audio_utils::CircularBuffer<std::complex<double>> *octave_cqt_buffer =
                transform_.get_octave_cqt_buffer(octave);
            for (int tone = 0; tone < BinsPerOctave; tone++)
            {
                values[tone] = octave_cqt_buffer[tone].pull_delay_sample(0);
            }
            return values;
        }

        std::vector<double> get_octave_bin_frequencies(const int octave)
        {
            std::vector<double> values(BinsPerOctave, 0.);
            const double *octave_bin_frequencies = transform_.get_octave_bin_frequencies(octave);
            for (int tone = 0; tone < BinsPerOctave; tone++)
            {
                values[tone] = octave_bin_frequencies[tone];
            }
            return values;
        }

    private:
        SlidingCqt<BinsPerOctave, OctaveCount, Windowing> transform_;
    };
}
