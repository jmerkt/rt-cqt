/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

#include <algorithm>
#include <atomic>
#include <complex>
#include <cstddef>

#include "../submodules/audio-utils/include/utils.h"
#include "resampling_filterbank.h"
#include "util.h"

namespace rt_cqt
{

    using namespace std::complex_literals;

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing = false>
    class SlidingCqt
    {
    public:
        SlidingCqt();
        ~SlidingCqt() = default;

        void init(const double sample_rate, const int block_size);

        void set_concert_pitch(double concert_pitch);
        inline void recalculate_kernels() { kernels_dirty_.store(true); };

        void input_block(double *const data, const int block_size);
        double *output_block(const int block_size);

        inline audio_utils::CircularBuffer<std::complex<double>> *get_octave_cqt_buffer(const int octave)
        {
            return &cqt_data_[octave][0];
        };
        inline std::size_t get_samples_to_process(const int octave) { return samples_to_process_[octave]; };
        inline void pull_bin_cqt_data(const int octave, const int tone, std::complex<double> *const data);
        inline void push_bin_cqt_data(const int octave, const int tone, std::complex<double> *const data);

        inline double get_octave_sample_rate(const int octave) { return octave_sample_rates_[octave]; };
        inline int get_octave_block_size(const int octave) { return octave_block_sizes_[octave]; };

        inline double *get_octave_bin_frequencies(const int octave) { return &bin_frequencies_[octave][0]; };

    private:
        void compute_kernels();

        double sample_rate_{48000.};
        double octave_sample_rates_[OctaveCount];
        int octave_block_sizes_[OctaveCount];
        double concert_pitch_{440.};

        std::atomic<bool> kernels_dirty_{true};

        ResamplingFilterbank<OctaveCount> filterbank_;

        audio_utils::CircularBuffer<double> delay_lines_[OctaveCount];

        // Pre-calculated exp stuff
        std::complex<double> exp_q_[OctaveCount][BinsPerOctave][3];
        std::complex<double> exp_q_nk_[OctaveCount][BinsPerOctave][3];
        std::complex<double> previous_transform_[OctaveCount][BinsPerOctave][3];

        double normalized_frequencies_[OctaveCount][BinsPerOctave][3];
        double window_lengths_[OctaveCount][BinsPerOctave];
        double inverse_window_lengths_[OctaveCount][BinsPerOctave];
        double bin_frequencies_[OctaveCount][BinsPerOctave];

        std::size_t samples_to_process_[OctaveCount];

        // CQT data
        audio_utils::CircularBuffer<std::complex<double>> cqt_data_[OctaveCount][BinsPerOctave];

        // Windowing
        // A periodic Hann window has a coherent gain of 1/2. The analysis
        // normalization therefore keeps its bin-centred amplitudes consistent
        // with the rectangular transform. The remaining 4/3 synthesis factor
        // complements it to the Hann energy normalization 1 / mean(w^2) = 8/3.
        static constexpr double WINDOW_ANALYSIS_NORMALIZATION{2.};
        static constexpr double WINDOW_SYNTHESIS_NORMALIZATION{4. / 3.};
        // Only the positive-frequency coefficients are stored. Restore the
        // conjugate half when producing a real signal.
        static constexpr double REAL_SYNTHESIS_NORMALIZATION{2.};
        static constexpr double WINDOW_COEFFICIENTS[3] = {0.5, -0.25, -0.25};
        static constexpr double Q_OFFSETS[3] = {0., -1., 1};

        // Buffers for block processing
        std::vector<double> input_samples_[OctaveCount];
        std::vector<double> delayed_input_samples_[OctaveCount][BinsPerOctave];
        std::vector<std::complex<double>> input_transforms_[OctaveCount][BinsPerOctave];

        std::vector<std::complex<double>> output_transforms_[OctaveCount][BinsPerOctave];
        std::vector<double> tone_outputs_[OctaveCount][BinsPerOctave];
        std::vector<double> octave_outputs_[OctaveCount];
    };

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing>
    SlidingCqt<BinsPerOctave, OctaveCount, Windowing>::SlidingCqt()
    {
        for (std::size_t octave = 0; octave < OctaveCount; octave++)
        {
            octave_sample_rates_[octave] = 48000.;
            octave_block_sizes_[octave] = 0;
            samples_to_process_[octave] = 0;
            for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
            {
                bin_frequencies_[octave][tone] = 0.;
                window_lengths_[octave][tone] = 0.;
                inverse_window_lengths_[octave][tone] = 0.;

                for (std::size_t window_index = 0; window_index < 3; window_index++)
                {
                    normalized_frequencies_[octave][tone][window_index] = 0.;
                    exp_q_nk_[octave][tone][window_index] = 0. + 0i;
                    previous_transform_[octave][tone][window_index] = 0. + 0.i;
                }
            }
        }
    }

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing>
    inline void SlidingCqt<BinsPerOctave, OctaveCount, Windowing>::init(const double sample_rate, const int block_size)
    {
        filterbank_.init(sample_rate, block_size, block_size * 2);
        sample_rate_ = filterbank_.get_origin_sample_rate();
        const int origin_block_size = filterbank_.get_origin_block_size();
        for (std::size_t octave = 0; octave < OctaveCount; octave++)
        {
            octave_sample_rates_[octave] = sample_rate_ / std::pow(2., octave);
            octave_block_sizes_[octave] =
                (origin_block_size / std::pow(2, octave) >= 1) ? origin_block_size / std::pow(2, octave) : 1;
        }
        compute_kernels();
        // initialize delay lines
        for (std::size_t octave = 0; octave < OctaveCount; octave++)
        {
            int maximum_delay_size = -1;
            for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
            {
                const int delay_size =
                    static_cast<int>(std::ceil(window_lengths_[octave][tone])) + octave_block_sizes_[octave] + 1;
                if (delay_size > maximum_delay_size)
                    maximum_delay_size = delay_size;
            }
            delay_lines_[octave].change_size(static_cast<std::size_t>(maximum_delay_size));
        }

        for (std::size_t octave = 0; octave < OctaveCount; octave++)
        {
            for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
            {
                const std::size_t octave_block_size = static_cast<std::size_t>(octave_block_sizes_[octave]);
                const std::size_t octave_buffer_size = std::max<std::size_t>(2, octave_block_size * 2);
                cqt_data_[octave][tone].change_size(octave_buffer_size);

                for (std::size_t window_index = 0; window_index < 3u; window_index++)
                {
                    previous_transform_[octave][tone][window_index] = 0. + 0.i;
                }
            }
        }

        // Buffers for block processing
        for (std::size_t octave = 0; octave < OctaveCount; octave++)
        {
            const std::size_t octave_block_size = octave_block_sizes_[octave];

            input_samples_[octave].resize(octave_block_size, 0.);
            octave_outputs_[octave].resize(octave_block_size, 0.);
            for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
            {
                delayed_input_samples_[octave][tone].resize(octave_block_size, 0.);
                input_transforms_[octave][tone].resize(octave_block_size, {0., 0.});

                output_transforms_[octave][tone].resize(octave_block_size, {0., 0.});
                tone_outputs_[octave][tone].resize(octave_block_size, 0.);
            }
        }
    };

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing>
    inline void SlidingCqt<BinsPerOctave, OctaveCount, Windowing>::set_concert_pitch(double concert_pitch)
    {
        concert_pitch_ = concert_pitch;
        recalculate_kernels();
    };

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing>
    inline void SlidingCqt<BinsPerOctave, OctaveCount, Windowing>::input_block(double *const data, const int block_size)
    {
        // check for new kernels
        if (kernels_dirty_.load())
        {
            kernels_dirty_.store(false);
            // calc the windows and give them to handlers
            compute_kernels();
        }

        // push data into multirate resampling
        filterbank_.input_block(data, block_size);
        // Process all CQT samples pushed into the stage buffers.
        for (std::size_t octave = 0; octave < OctaveCount; octave++)
        {
            BufferPtr input_buffer = filterbank_.get_stage_input_buffer(octave);
            const int octave_sample_count = input_buffer->get_available_sample_count();
            samples_to_process_[octave] = octave_sample_count;
            if (octave_sample_count <= 0)
            {
                continue;
            }

            input_buffer->pull_block(input_samples_[octave].data(), octave_sample_count);
            delay_lines_[octave].push_block(input_samples_[octave].data(), octave_sample_count);
            for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
            {
                const double window_length = window_lengths_[octave][tone];
                delay_lines_[octave].pull_delay_block(delayed_input_samples_[octave][tone].data(),
                                                      static_cast<int>(window_length) + octave_sample_count - 1,
                                                      octave_sample_count);
            }
            for (std::size_t sample = 0; sample < octave_sample_count; sample++)
            {
                // #pragma omp simd
                for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
                {
                    const double inverse_window_length = inverse_window_lengths_[octave][tone];
                    const double delayed_sample = delayed_input_samples_[octave][tone][sample];

                    const std::complex<double> complex_delayed_sample{delayed_sample, 0.};

                    if constexpr (Windowing == false)
                    {
                        const std::complex<double> exp_q = exp_q_[octave][tone][0];
                        const std::complex<double> exp_q_nk = exp_q_nk_[octave][tone][0];
                        const std::complex<double> previous_transform = previous_transform_[octave][tone][0];
                        const std::complex<double> modulated_input = exp_q * input_samples_[octave][sample];

                        const std::complex<double> transform =
                            exp_q_nk *
                            (previous_transform + (modulated_input - complex_delayed_sample) * inverse_window_length);

                        previous_transform_[octave][tone][0] = transform;
                        input_transforms_[octave][tone][sample] = transform;
                    }
                    else
                    {
                        std::complex<double> transform_sum = 0. + 0.i;
                        for (std::size_t window_index = 0; window_index < 3u; window_index++)
                        {
                            const std::complex<double> exp_q = exp_q_[octave][tone][window_index];
                            const std::complex<double> exp_q_nk = exp_q_nk_[octave][tone][window_index];
                            const std::complex<double> previous_transform =
                                previous_transform_[octave][tone][window_index];
                            const std::complex<double> modulated_input = exp_q * input_samples_[octave][sample];

                            const std::complex<double> transform =
                                exp_q_nk * (previous_transform +
                                            (modulated_input - complex_delayed_sample) * inverse_window_length);

                            previous_transform_[octave][tone][window_index] = transform;

                            transform_sum += WINDOW_COEFFICIENTS[window_index] * transform;
                        }
                        input_transforms_[octave][tone][sample] = transform_sum * WINDOW_ANALYSIS_NORMALIZATION;
                    }
                }
            }
            for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
            {
                cqt_data_[octave][tone].push_block(input_transforms_[octave][tone].data(), octave_sample_count);
            }
        }
    };

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing>
    inline double *SlidingCqt<BinsPerOctave, OctaveCount, Windowing>::output_block(const int block_size)
    {
        for (std::size_t octave = 0; octave < OctaveCount; octave++)
        {
            const std::size_t octave_sample_count = samples_to_process_[octave];
            if (octave_sample_count <= 0)
            {
                continue;
            }

            for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
            {
                cqt_data_[octave][tone].pull_block(output_transforms_[octave][tone].data(), octave_sample_count);
            }
            for (std::size_t sample = 0; sample < octave_sample_count; sample++)
            {
                // #pragma omp simd
                for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
                {
                    const std::complex<double> exp_q_nk = exp_q_nk_[octave][tone][0];
                    const std::complex<double> transform = output_transforms_[octave][tone][sample];
                    tone_outputs_[octave][tone][sample] = (transform * exp_q_nk).real();
                }
                octave_outputs_[octave][sample] = 0.;
                for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
                {
                    octave_outputs_[octave][sample] += tone_outputs_[octave][tone][sample];
                }
                octave_outputs_[octave][sample] *= REAL_SYNTHESIS_NORMALIZATION;
                if constexpr (Windowing)
                {
                    octave_outputs_[octave][sample] *= WINDOW_SYNTHESIS_NORMALIZATION;
                }
            }

            BufferPtr output_buffer = filterbank_.get_stage_output_buffer(octave);
            output_buffer->push_block(octave_outputs_[octave].data(), octave_sample_count);
        }
        return filterbank_.output_block(block_size);
    };

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing>
    inline void SlidingCqt<BinsPerOctave, OctaveCount, Windowing>::compute_kernels()
    {
        const double q_initial = 1. / (std::pow(2., 1. / static_cast<double>(BinsPerOctave)) - 1.);

        const double reference_frequency = compute_reference_frequency(concert_pitch_);
        for (std::size_t octave = 0; octave < OctaveCount; octave++)
        {
            // sample_rate
            const double sample_rate = octave_sample_rates_[octave];
            for (std::size_t tone = 0; tone < BinsPerOctave; tone++)
            {
                // bin_frequency
                const double bin_frequency = compute_bin_frequency(reference_frequency, BinsPerOctave, octave, tone);
                bin_frequencies_[octave][tone] = bin_frequency;
                // Nk
                window_lengths_[octave][tone] = std::floor((sample_rate / bin_frequency) * q_initial);
                inverse_window_lengths_[octave][tone] = 1. / window_lengths_[octave][tone];

                for (std::size_t window_index = 0; window_index < 3u; window_index++)
                {
                    // Q
                    normalized_frequencies_[octave][tone][window_index] =
                        window_lengths_[octave][tone] * bin_frequency / sample_rate;
                    normalized_frequencies_[octave][tone][window_index] += Q_OFFSETS[window_index];

                    // exp multiplication
                    exp_q_[octave][tone][window_index] = std::exp(-1i * audio_utils::two_pi<double>() *
                                                                  normalized_frequencies_[octave][tone][window_index]);
                    exp_q_nk_[octave][tone][window_index] = std::exp(
                        1i * audio_utils::two_pi<double>() * normalized_frequencies_[octave][tone][window_index] *
                        inverse_window_lengths_[octave][tone]);
                }
            }
        }
    };

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing>
    inline void SlidingCqt<BinsPerOctave, OctaveCount, Windowing>::pull_bin_cqt_data(const int octave,
                                                                                     const int tone,
                                                                                     std::complex<double> *const data)
    {
        cqt_data_[octave][tone].pull_block(data, samples_to_process_[octave]);
    };

    template <std::size_t BinsPerOctave, std::size_t OctaveCount, bool Windowing>
    inline void SlidingCqt<BinsPerOctave, OctaveCount, Windowing>::push_bin_cqt_data(const int octave,
                                                                                     const int tone,
                                                                                     std::complex<double> *const data)
    {
        cqt_data_[octave][tone].push_block(data, samples_to_process_[octave]);
    }

}
