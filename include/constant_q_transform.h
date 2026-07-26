/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

/*
 * Limitations:
 * - Fixed Hann window function - for resynthesis a HopSize >= FFT_SIZE / 2 should be used - HopSize == FFT_SIZE will
 * produce NaNs
 * - No Zero Padding as of now
 */

#include "../submodules/audio-utils/include/utils.h"
#include "resampling_filterbank.h"
#include "util.h"
#include <array>
#include <atomic>
#include <memory>

#define SIMD_SZ 1
#define PFFFT_ENABLE_DOUBLE
#include "../submodules/pffft/pffft.hpp"
#include <complex>

namespace Cqt
{

    using namespace std::complex_literals;
    constexpr int FFT_SIZE{512};
    constexpr int FFT_DOMAIN_SIZE{FFT_SIZE / 2};
    constexpr double KERNEL_THRESHOLD{1.e-4};
    constexpr double WINDOW_ENERGY_LOSS_COMPENSATION{1.63};   // Fixed for Hanning window as of now
    constexpr double WINDOW_AMPLITUDE_LOSS_COMPENSATION{2.0}; // Fixed for Hanning window as of now

    typedef pffft::AlignedVector<pffft::Types<double>::Complex> CplxVector;
    typedef pffft::AlignedVector<double> RealVector;
    typedef RealVector TimeDataType;
    typedef CplxVector CqtBufferType;

    /*
    Structure to schedule transformation timings.
    */
    class ScheduleElement
    {
    public:
        ScheduleElement(const int sample,
                        const int octave,
                        const int delay_at_octave_rate,
                        const int synthesis_offset = 0) :
            sample_(sample),
            octave_(octave),
            delay_at_octave_rate_(delay_at_octave_rate),
            synthesis_offset_(synthesis_offset) {};
        ~ScheduleElement() = default;
        int sample() const { return sample_; };
        int octave() const { return octave_; };
        int delay_at_octave_rate() const { return delay_at_octave_rate_; };
        int synthesis_offset() const { return synthesis_offset_; };

    private:
        int sample_{0}; // Position in the filterbank's internal processing block
        int octave_{0};
        int delay_at_octave_rate_{0}; // Delay in samples in the corresponding octave sample buffer
        int synthesis_offset_{0};     // Forward offset in the current octave output block
    };

    /*
    Handles the transformations of the single octaves.
    */
    template <int BinsPerOctave>
    class TransformationHandler
    {
    public:
        TransformationHandler();
        ~TransformationHandler() = default;

        void init(const int hop_size);
        void init_buffers(BufferPtr input_buffer = nullptr, BufferPtr output_buffer = nullptr);
        void init_kernels(const CplxVector *const kernel_array,
                          const CplxVector *const inverse_kernel_array,
                          const std::vector<int> *const kernel_mask,
                          const std::vector<int> *const inverse_kernel_mask);
        void resize_output_buffer(const int block_size);

        void cqt_transform(const ScheduleElement schedule);
        void icqt_transform(const ScheduleElement schedule);

        static void calculate_window(double *const window_data, const int size);
        static void calculate_inverse_window(double *const window_data,
                                             double *const inverse_window_data,
                                             const int size,
                                             const int hop_size);

        inline CqtBufferType *get_cqt_buffer() { return &cqt_buffer_; };
        inline BufferPtr get_output_buffer() { return stage_output_; };

    private:
        double window_[FFT_SIZE];
        double inverse_window_[FFT_SIZE];
        // kernel storage
        CplxVector kernels_[BinsPerOctave];
        CplxVector inverse_kernels_[BinsPerOctave];
        std::vector<int> kernel_masks_[BinsPerOctave];
        std::vector<int> inverse_kernel_masks_[BinsPerOctave];

        // scaling
        const double fft_scaling_factor_{1. / std::sqrt(static_cast<double>(FFT_SIZE))};
        // pffft
        pffft::Fft<double> fft_;
        RealVector fft_input_;
        CplxVector spectrum_;
        CplxVector inverse_spectrum_;
        CqtBufferType cqt_buffer_;
        RealVector output_;
        RealVector inverse_fft_output_;

        // input / output buffers
        BufferPtr stage_input_;
        BufferPtr stage_output_;
    };

    template <int BinsPerOctave>
    TransformationHandler<BinsPerOctave>::TransformationHandler() :
        cqt_buffer_(BinsPerOctave),
        fft_(FFT_SIZE)
    {
        // hard-coded Hann window as of now
        calculate_window(window_, FFT_SIZE);
        // fft and buffers
        fft_input_ = fft_.valueVector();
        spectrum_ = fft_.spectrumVector();
        inverse_spectrum_ = fft_.spectrumVector();
        inverse_fft_output_ = fft_.valueVector();
        std::fill(inverse_fft_output_.begin(), inverse_fft_output_.end(), 0.);
        output_ = fft_.valueVector();
        std::fill(output_.begin(), output_.end(), 0.);
        for (int tone = 0; tone < BinsPerOctave; tone++)
        {
            cqt_buffer_[tone] = 0. + 0.i;
        }
        // kernels
        for (int tone = 0; tone < BinsPerOctave; tone++)
        {
            kernels_[tone] = fft_.spectrumVector();
            inverse_kernels_[tone] = fft_.spectrumVector();
        }
    }

    template <int BinsPerOctave>
    inline void TransformationHandler<BinsPerOctave>::init(const int hop_size)
    {
        calculate_inverse_window(window_, inverse_window_, FFT_SIZE, hop_size);
    }

    template <int BinsPerOctave>
    inline void TransformationHandler<BinsPerOctave>::init_buffers(BufferPtr input_buffer, BufferPtr output_buffer)
    {
        stage_input_ = input_buffer;
        stage_output_ = output_buffer;
    }

    template <int BinsPerOctave>
    inline void TransformationHandler<BinsPerOctave>::init_kernels(const CplxVector *const kernel_array,
                                                                   const CplxVector *const inverse_kernel_array,
                                                                   const std::vector<int> *const kernel_mask,
                                                                   const std::vector<int> *const inverse_kernel_mask)
    {
        for (int tone = 0; tone < BinsPerOctave; tone++)
        {
            for (int i = 0; i < FFT_DOMAIN_SIZE; i++)
            {
                kernels_[tone][i] = kernel_array[tone].at(i);
                inverse_kernels_[tone][i] = inverse_kernel_array[tone].at(i);
            }
            kernel_masks_[tone].resize(kernel_mask[tone].size(), 0);
            for (size_t i = 0; i < kernel_mask[tone].size(); i++)
            {
                kernel_masks_[tone][i] = kernel_mask[tone][i];
            }
            inverse_kernel_masks_[tone].resize(inverse_kernel_mask[tone].size(), 0);
            for (size_t i = 0; i < inverse_kernel_mask[tone].size(); i++)
            {
                inverse_kernel_masks_[tone][i] = inverse_kernel_mask[tone][i];
            }
        }
    };

    template <int BinsPerOctave>
    inline void TransformationHandler<BinsPerOctave>::resize_output_buffer(const int block_size)
    {
        output_.resize(FFT_SIZE + block_size);
    }

    template <int BinsPerOctave>
    inline void TransformationHandler<BinsPerOctave>::cqt_transform(const ScheduleElement schedule)
    {
        // collect the fft input data
        stage_input_->pull_delay_block(fft_input_.data(),
                                       static_cast<int>(FFT_SIZE) + schedule.delay_at_octave_rate() - 1,
                                       static_cast<int>(FFT_SIZE));
        // apply time window
        for (int i = 0; i < FFT_SIZE; i++)
        {
            fft_input_[i] *= window_[i];
        }
        // fft
        fft_.forward(fft_input_, spectrum_);
        // scale data
        for (int i = 0; i < FFT_DOMAIN_SIZE; i++)
        {
            spectrum_[i] *= fft_scaling_factor_;
        }
        // kernel multipications
        for (int tone = 0; tone < BinsPerOctave; tone++)
        {
            cqt_buffer_[tone] = 0. + 0.i;
            for (size_t i = 0; i < kernel_masks_[tone].size(); i++)
            {
                const int index = kernel_masks_[tone][i];
                cqt_buffer_[tone] += spectrum_[index] * kernels_[tone][index];
            }
        }
    };

    template <int BinsPerOctave>
    inline void TransformationHandler<BinsPerOctave>::icqt_transform(const ScheduleElement schedule)
    {
        // kernel multiplications
        for (int i = 0; i < FFT_DOMAIN_SIZE; i++)
        {
            inverse_spectrum_[i] = 0. + 0.i;
        }
        for (int tone = 0; tone < BinsPerOctave; tone++)
        {
            for (size_t i = 0; i < inverse_kernel_masks_[tone].size(); i++)
            {
                const int index = inverse_kernel_masks_[tone][i];
                inverse_spectrum_[index] += cqt_buffer_[tone] * inverse_kernels_[tone][index];
            }
        }
        // ifft
        fft_.inverse(inverse_spectrum_, inverse_fft_output_);
        // scale and window data
        for (int i = 0; i < FFT_SIZE; i++)
        {
            inverse_fft_output_[i] *= fft_scaling_factor_;
        }
        for (int i = 0; i < FFT_SIZE; i++)
        {
            inverse_fft_output_[i] *= inverse_window_[i];
        }

        // overlap-add
        std::fill(output_.begin(), output_.end(), 0.);
        //// pull whats left from the previous transform
        stage_output_->pull_block(output_.data(), stage_output_->get_write_read_distance());
        //// add new data
        int output_index = 0;
        for (int i = schedule.synthesis_offset(); i < (schedule.synthesis_offset() + FFT_SIZE); i++)
        {
            output_[i] += inverse_fft_output_[output_index];
            output_index++;
        }
        //// push new data
        stage_output_->push_block(output_.data(), FFT_SIZE + schedule.synthesis_offset());
    };

    template <int BinsPerOctave>
    inline void TransformationHandler<BinsPerOctave>::calculate_window(double *const window_data, const int size)
    {
        for (int i = 0; i < size; i++)
        {
            window_data[i] = (1. / 2.) * (1. - std::cos((2. * audio_utils::pi<double>() * static_cast<double>(i)) /
                                                        static_cast<double>(size - 1)));
        }
    }

    template <int BinsPerOctave>
    inline void TransformationHandler<BinsPerOctave>::calculate_inverse_window(double *const window_data,
                                                                               double *const inverse_window_data,
                                                                               const int size,
                                                                               const int hop_size)
    {
        std::vector<double> window_sum(size, 0.);
        for (int i = 0; i < size; i += hop_size)
        {
            for (int j = 0; j < size; j++)
            {
                window_sum[(i + j) % size] += window_data[j] * window_data[j];
            }
        }
        for (int i = 0; i < size; i++)
        {
            inverse_window_data[i] = window_data[i] / window_sum[i] * std::pow(WINDOW_ENERGY_LOSS_COMPENSATION, 2);
        }
    }

    /*
    Main CQT class
    */
    template <int BinsPerOctave, int OctaveCount>
    class ConstantQTransform
    {
    public:
        ConstantQTransform();
        ~ConstantQTransform() = default;

        void init(int hop_size);
        void init(std::vector<int> octave_hop_sizes);
        void init_sample_rate(double sample_rate, const int block_size);
        void set_concert_pitch(double concert_pitch);
        inline void recalculate_kernels() { kernels_dirty_.store(true); };

        void input_block(double *const data, const int block_size);
        double *output_block(const int block_size);
        void cqt(const ScheduleElement schedule);
        void icqt(const ScheduleElement schedule);

        inline std::vector<ScheduleElement> &get_cqt_schedule() { return schedule_; };
        inline CqtBufferType *get_octave_cqt_buffer(const int octave)
        {
            return transform_handlers_[octave].get_cqt_buffer();
        };
        inline BufferPtr get_octave_output_buffer(const int octave)
        {
            return transform_handlers_[octave].get_output_buffer();
        };
        inline int get_hop_size(const int octave) { return hop_sizes_[octave]; };
        inline size_t get_latency_samples(const int octave) { return latency_samples_[octave]; };
        inline double get_latency_ms(const int octave) { return latency_ms_[octave]; };
        inline double get_octave_sample_rate(const int octave) { return octave_sample_rates_[octave]; };
        inline std::vector<std::vector<double>> &get_kernel_frequencies() { return kernel_frequencies_; };
        inline std::vector<std::vector<double>> &get_inverse_kernel_frequencies()
        {
            return inverse_kernel_frequencies_;
        };
        inline void reset_kernel_frequencies() { init_kernel_frequencies(); };

    protected:
        void init_kernel_frequencies();
        void calculate_kernels();

        double concert_pitch_{440.};
        int bin_count_;
        int overlaps_[OctaveCount];
        double latency_ms_[OctaveCount];
        size_t latency_samples_[OctaveCount];
        int hop_sizes_[OctaveCount];
        double sample_rate_;
        double octave_sample_rates_[OctaveCount];
        double octave_rate_ratios_[OctaveCount];

        TransformationHandler<BinsPerOctave> transform_handlers_[OctaveCount];
        ResamplingFilterbank<OctaveCount> filterbank_;
        size_t sample_counters_[OctaveCount];

        std::vector<ScheduleElement> schedule_;

        pffft::Fft<std::complex<double>> fft_;
        pffft::Fft<std::complex<double>> inverse_kernel_fft_;
        pffft::Fft<double> fft_allocator_;
        CplxVector kernel_spectrum_;
        CplxVector inverse_kernel_spectrum_;

        CplxVector kernels_[BinsPerOctave];
        CplxVector inverse_kernels_[BinsPerOctave];
        std::vector<int> kernel_masks_[BinsPerOctave];
        std::vector<int> inverse_kernel_masks_[BinsPerOctave];
        CplxVector time_kernels_[BinsPerOctave];
        CplxVector inverse_time_kernels_[BinsPerOctave];
        RealVector analysis_window_;
        std::atomic<bool> kernels_dirty_{true};

        std::vector<std::vector<double>> kernel_frequencies_;
        std::vector<std::vector<double>> inverse_kernel_frequencies_;
    };

    template <int BinsPerOctave, int OctaveCount>
    ConstantQTransform<BinsPerOctave, OctaveCount>::ConstantQTransform() :
        fft_(FFT_SIZE),
        inverse_kernel_fft_(FFT_SIZE),
        fft_allocator_(FFT_SIZE)
    {
        bin_count_ = BinsPerOctave * OctaveCount;

        kernel_spectrum_ = fft_.spectrumVector();
        inverse_kernel_spectrum_ = fft_.spectrumVector();
        // configure all the buffer sizes
        for (int tone = 0; tone < BinsPerOctave; tone++)
        {
            time_kernels_[tone] = fft_.valueVector();
            inverse_time_kernels_[tone] = fft_.valueVector();
            kernels_[tone] = fft_allocator_.spectrumVector();
            inverse_kernels_[tone] = fft_allocator_.spectrumVector();
        }
        // generate window function
        analysis_window_.resize(FFT_SIZE);
        TransformationHandler<BinsPerOctave>::calculate_window(analysis_window_.data(), FFT_SIZE);
        // transformation in/out buffers
        kernel_frequencies_.resize(OctaveCount);
        inverse_kernel_frequencies_.resize(OctaveCount);
        for (int octave = 0; octave < OctaveCount; octave++)
        {
            transform_handlers_[octave].init_buffers(filterbank_.get_stage_input_buffer(octave),
                                                     filterbank_.get_stage_output_buffer(octave));
            kernel_frequencies_[octave].resize(BinsPerOctave, 0.);
            inverse_kernel_frequencies_[octave].resize(BinsPerOctave, 0.);
            sample_counters_[octave] = 0;
        }
    };

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::init(int hop_size)
    {
        init_kernel_frequencies();

        hop_size = audio_utils::clip<int>(hop_size, 1, FFT_SIZE);
        for (int octave = 0; octave < OctaveCount; octave++)
        {
            hop_sizes_[octave] = hop_size;
            overlaps_[octave] = FFT_SIZE - hop_size;
        }
        for (int octave = 0; octave < OctaveCount; octave++)
        {
            transform_handlers_[octave].init(hop_sizes_[octave]);
        }
    }

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::init(std::vector<int> octave_hop_sizes)
    {
        init_kernel_frequencies();

        assert(octave_hop_sizes.size() == OctaveCount);
        for (int octave = 0; octave < OctaveCount; octave++)
        {
            int hop_size = audio_utils::clip<int>(octave_hop_sizes.at(octave), 1, FFT_SIZE);
            hop_sizes_[octave] = hop_size;
            overlaps_[octave] = FFT_SIZE - hop_size;
        }
        for (int octave = 0; octave < OctaveCount; octave++)
        {
            transform_handlers_[octave].init(hop_sizes_[octave]);
        }
    }

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::init_sample_rate(double sample_rate,
                                                                                 const int block_size)
    {
        filterbank_.init(sample_rate, block_size, block_size + FFT_SIZE);
        sample_rate_ = filterbank_.get_origin_sample_rate();

        for (int octave = 0; octave < OctaveCount; octave++)
        {
            // latency per octave
            latency_samples_[octave] = static_cast<size_t>(hop_sizes_[octave]) *
                                       static_cast<size_t>(std::pow(2, octave)) *
                                       static_cast<size_t>(std::pow(2, filterbank_.get_origin_downsampling()));
            sample_counters_[octave] = latency_samples_[octave];
            // samplerates
            octave_sample_rates_[octave] = sample_rate_ / std::pow(2., octave);
            latency_ms_[octave] = static_cast<double>(hop_sizes_[octave]) / octave_sample_rates_[octave] * 1000.;
            octave_rate_ratios_[octave] = octave_sample_rates_[octave] / sample_rate;
        }
        for (int octave = 0; octave < OctaveCount; octave++)
        {
            transform_handlers_[octave].resize_output_buffer(filterbank_.get_origin_block_size());
        }
        // calc the windows and give em to handlers
        recalculate_kernels();
    };

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::init_kernel_frequencies()
    {
        const double reference_frequency = compute_reference_frequency(concert_pitch_);
        for (int octave = 0; octave < OctaveCount; octave++)
        {
            for (int tone = 0; tone < BinsPerOctave; tone++)
            {
                kernel_frequencies_[octave][tone] =
                    compute_bin_frequency(reference_frequency, BinsPerOctave, octave, tone);
                inverse_kernel_frequencies_[octave][tone] = kernel_frequencies_[octave][tone];
            }
        }
    };

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::calculate_kernels()
    {
        // calculate the time domain kernels
        for (int k = 0; k < BinsPerOctave; k++)
        {
            const double kernel_frequency = kernel_frequencies_[0][k];
            const double inverse_kernel_frequency = inverse_kernel_frequencies_[0][k];
            for (int n = 0; n < FFT_SIZE; n++)
            {
                time_kernels_[k][n] = std::conj((1. / static_cast<double>(FFT_SIZE)) * analysis_window_[n] *
                                                std::exp(-1i * 2. * audio_utils::pi<double>() * static_cast<double>(n) *
                                                         (kernel_frequency / octave_sample_rates_[0])));
                inverse_time_kernels_[k][n] =
                    std::conj((1. / static_cast<double>(FFT_SIZE)) * analysis_window_[n] *
                              std::exp(-1i * 2. * audio_utils::pi<double>() * static_cast<double>(n) *
                                       (inverse_kernel_frequency / octave_sample_rates_[0])));
            }
        }
        // fft transform kernels and extract necessary (right side of the spectrum) parts
        for (int k = 0; k < BinsPerOctave; k++)
        {
            fft_.forward(time_kernels_[k], kernel_spectrum_);
            inverse_kernel_fft_.forward(inverse_time_kernels_[k], inverse_kernel_spectrum_);
            // extract real part
            for (int n = 0; n < FFT_DOMAIN_SIZE; n++)
            {
                kernels_[k][n] = kernel_spectrum_[n];
                inverse_kernels_[k][n] = std::conj(inverse_kernel_spectrum_[n]);
            }
        }
        // mark relevant kernel values
        for (int k = 0; k < BinsPerOctave; k++)
        {
            kernel_masks_[k].clear();
            inverse_kernel_masks_[k].clear();
            for (int n = 0; n < FFT_DOMAIN_SIZE; n++)
            {
                const double kernel_magnitude = std::abs(kernels_[k][n]);
                if (kernel_magnitude > KERNEL_THRESHOLD)
                {
                    kernel_masks_[k].push_back(n);
                }
                const double inverse_kernel_magnitude = std::abs(inverse_kernels_[k][n]);
                if (inverse_kernel_magnitude > KERNEL_THRESHOLD)
                {
                    inverse_kernel_masks_[k].push_back(n);
                }
            }
        }
        // pass kernels to handlers
        for (int octave = 0; octave < OctaveCount; octave++)
        {
            transform_handlers_[octave].init_kernels(kernels_, inverse_kernels_, kernel_masks_, inverse_kernel_masks_);
        }
    };

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::set_concert_pitch(double concert_pitch)
    {
        concert_pitch_ = concert_pitch;
        init_kernel_frequencies();
        recalculate_kernels();
    };

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::input_block(double *const data, const int block_size)
    {
        // check for new kernels
        if (kernels_dirty_.load())
        {
            kernels_dirty_.store(false);
            // calc the windows and give them to handlers
            calculate_kernels();
        }
        // process Filterbank and create Schedule
        filterbank_.input_block(data, block_size);
        const int processed_input_size = filterbank_.get_last_processed_input_size();
        // determine cqt positions and schedule them
        schedule_.clear();
        std::array<int, OctaveCount> first_delay_at_octave_rate;
        first_delay_at_octave_rate.fill(-1);
        for (int i = 0; i < processed_input_size; i++)
        {
            for (int octave = (OctaveCount - 1); octave >= 0;
                 octave--) // starting with lowest pitched octave for historical reasons
            {
                sample_counters_[octave]++;
                if (sample_counters_[octave] >= latency_samples_[octave])
                {
                    sample_counters_[octave] = 0;
                    const int delay_at_octave_rate = static_cast<int>(
                        static_cast<double>(processed_input_size - i - 1) * octave_rate_ratios_[octave]);
                    if (first_delay_at_octave_rate[octave] < 0)
                    {
                        first_delay_at_octave_rate[octave] = delay_at_octave_rate;
                    }
                    const int synthesis_offset = first_delay_at_octave_rate[octave] - delay_at_octave_rate;
                    schedule_.push_back({i, octave, delay_at_octave_rate, synthesis_offset});
                }
            }
        }
    };

    template <int BinsPerOctave, int OctaveCount>
    inline double *ConstantQTransform<BinsPerOctave, OctaveCount>::output_block(const int block_size)
    {
        return filterbank_.output_block(block_size);
    };

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::cqt(const ScheduleElement schedule)
    {
        transform_handlers_[schedule.octave()].cqt_transform(schedule);
    };

    template <int BinsPerOctave, int OctaveCount>
    inline void ConstantQTransform<BinsPerOctave, OctaveCount>::icqt(const ScheduleElement schedule)
    {
        transform_handlers_[schedule.octave()].icqt_transform(schedule);
    };

};
