/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "../submodules/audio-utils/include/circular_buffer.h"
#include "resampling.h"

/*
This class handles multirate resampling of the input data. Input callbacks may
have any size up to the size supplied to init(). Internally, callbacks are
collected into a block that is divisible by every power-of-two resampling stage.
Consequently all resampling is block based and every internal block produces at
least one sample for the lowest stage.

The supported input sample rates are power-of-two multiples of either 44.1 kHz
or 48 kHz.
*/

namespace rt_cqt
{

    using BufferPtr = audio_utils::CircularBuffer<double> *;

    constexpr double FILTER_TRANSITION_BANDWIDTH{0.1};
    constexpr unsigned ALLPASS_COUNT{3};

    template <int StageCount>
    class ResamplingFilterbank
    {
        static_assert(StageCount > 0, "A resampling filterbank needs at least one stage");
        static_assert(StageCount < 30, "StageCount is too large for 32-bit block sizes");

    public:
        ResamplingFilterbank() = default;
        ~ResamplingFilterbank() = default;

        void init(const double sample_rate, const int block_size, const int buffer_size);

        void input_block(const double *const data, const int block_size);
        double *output_block(const int block_size);

        double get_origin_sample_rate() const { return origin_sample_rate_; }
        int get_origin_block_size() const { return origin_block_size_; }
        int get_processing_block_size() const { return processing_block_size_; }
        int get_latency_samples() const { return processing_block_size_; }
        int get_last_processed_input_size() const { return last_processed_input_size_; }
        BufferPtr get_stage_input_buffer(const int stage)
        {
            return &stage_input_buffers_.at(static_cast<std::size_t>(stage));
        }
        BufferPtr get_stage_output_buffer(const int stage)
        {
            return &stage_output_buffers_.at(static_cast<std::size_t>(stage));
        }
        int get_origin_downsampling() const { return origin_downsampling_; }

    private:
        static constexpr int RESAMPLING_STAGE_COUNT = StageCount - 1;

        static bool is_power_of_two(const int value);
        static int power_of_two_exponent(const int value);
        static int round_up_to_multiple(const int value, const int multiple);

        void process_input();
        void process_output();
        void push_output_samples(const double *const data, const int block_size);
        void pull_output_samples(double *const data, const int block_size);

        ResamplingHandler<double, ALLPASS_COUNT> input_resampler_;
        std::array<HalfBandLowpass<double, ALLPASS_COUNT>, RESAMPLING_STAGE_COUNT> downsampling_filters_;
        std::array<HalfBandLowpass<double, ALLPASS_COUNT>, RESAMPLING_STAGE_COUNT> upsampling_filters_;
        std::array<audio_utils::CircularBuffer<double>, StageCount> stage_input_buffers_;
        std::array<audio_utils::CircularBuffer<double>, StageCount> stage_output_buffers_;

        double origin_sample_rate_{48000.};
        int origin_block_size_{0};
        int origin_downsampling_{0};
        int processing_block_size_{0};
        int maximum_callback_block_size_{0};
        int lowest_stage_block_size_{0};

        std::vector<double> input_data_;
        std::vector<double> lowest_stage_output_;
        std::vector<double> output_data_;
        int input_data_size_{0};
        int last_processed_input_size_{0};
        int pending_output_blocks_{0};

        audio_utils::CircularBuffer<double> output_queue_;
    };

    template <int StageCount>
    inline bool ResamplingFilterbank<StageCount>::is_power_of_two(const int value)
    {
        return value > 0 && (value & (value - 1)) == 0;
    }

    template <int StageCount>
    inline int ResamplingFilterbank<StageCount>::power_of_two_exponent(const int value)
    {
        assert(is_power_of_two(value));
        int exponent = 0;
        for (int remaining = value; remaining > 1; remaining >>= 1)
        {
            ++exponent;
        }
        return exponent;
    }

    template <int StageCount>
    inline int ResamplingFilterbank<StageCount>::round_up_to_multiple(const int value, const int multiple)
    {
        assert(value > 0);
        assert(multiple > 0);
        const long long result = ((static_cast<long long>(value) + multiple - 1) / multiple) * multiple;
        if (result > std::numeric_limits<int>::max())
        {
            throw std::overflow_error("The requested resampling block size is too large");
        }
        return static_cast<int>(result);
    }

    template <int StageCount>
    inline void
    ResamplingFilterbank<StageCount>::init(const double sample_rate, const int block_size, const int buffer_size)
    {
        if (!std::isfinite(sample_rate) || sample_rate <= 0.)
        {
            throw std::invalid_argument("The sample rate must be finite and positive");
        }
        if (block_size <= 0 || buffer_size <= 0)
        {
            throw std::invalid_argument("Block and buffer sizes must be positive");
        }

        const long long rounded_sample_rate = std::llround(sample_rate);
        if (std::abs(sample_rate - static_cast<double>(rounded_sample_rate)) > 1.e-6 ||
            rounded_sample_rate > std::numeric_limits<int>::max())
        {
            throw std::invalid_argument("Only integer sample rates are supported");
        }

        const int sample_rate_integer = static_cast<int>(rounded_sample_rate);
        int origin_factor = 0;
        if ((sample_rate_integer % 44100) == 0 && is_power_of_two(sample_rate_integer / 44100))
        {
            origin_sample_rate_ = 44100.;
            origin_factor = sample_rate_integer / 44100;
        }
        else if ((sample_rate_integer % 48000) == 0 && is_power_of_two(sample_rate_integer / 48000))
        {
            origin_sample_rate_ = 48000.;
            origin_factor = sample_rate_integer / 48000;
        }
        else
        {
            throw std::invalid_argument("The sample rate must be a power-of-two multiple of 44.1 kHz or 48 kHz");
        }

        origin_downsampling_ = power_of_two_exponent(origin_factor);
        const int filterbank_factor = 1 << RESAMPLING_STAGE_COUNT;
        if (origin_factor > (std::numeric_limits<int>::max() / filterbank_factor))
        {
            throw std::overflow_error("The resampling factor is too large");
        }
        const int input_alignment = origin_factor * filterbank_factor;

        maximum_callback_block_size_ = block_size;
        processing_block_size_ = round_up_to_multiple(block_size, input_alignment);
        origin_block_size_ = processing_block_size_ / origin_factor;
        lowest_stage_block_size_ = origin_block_size_ / filterbank_factor;

        input_data_.assign(static_cast<std::size_t>(processing_block_size_), 0.);
        lowest_stage_output_.assign(static_cast<std::size_t>(lowest_stage_block_size_), 0.);
        output_data_.assign(static_cast<std::size_t>(maximum_callback_block_size_), 0.);
        input_data_size_ = 0;
        last_processed_input_size_ = 0;
        pending_output_blocks_ = 0;

        input_resampler_.init(origin_downsampling_, processing_block_size_, DirectionConfig::DownUp);

        for (int stage = 0; stage < RESAMPLING_STAGE_COUNT; ++stage)
        {
            const int stage_input_size = origin_block_size_ / (1 << stage);
            const int stage_output_size = stage_input_size / 2;
            downsampling_filters_[static_cast<std::size_t>(stage)].init(
                stage_input_size, true, FILTER_TRANSITION_BANDWIDTH);
            upsampling_filters_[static_cast<std::size_t>(stage)].init(
                stage_output_size, false, FILTER_TRANSITION_BANDWIDTH);
        }

        for (int stage = 0; stage < StageCount; ++stage)
        {
            const int stage_block_size = origin_block_size_ / (1 << stage);
            const int required_buffer_size = std::max(buffer_size, stage_block_size * 2);
            stage_input_buffers_[static_cast<std::size_t>(stage)].change_size(required_buffer_size);
            stage_output_buffers_[static_cast<std::size_t>(stage)].change_size(required_buffer_size);
        }

        const std::size_t output_queue_capacity = static_cast<std::size_t>(processing_block_size_) * 2U;
        output_queue_.change_size(output_queue_capacity);
        for (int sample = 0; sample < processing_block_size_; ++sample)
        {
            output_queue_.push_sample(0.);
        }
    }

    template <int StageCount>
    inline void ResamplingFilterbank<StageCount>::process_input()
    {
        double *input = input_resampler_.process_down(input_data_.data());
        int data_size = origin_block_size_;
        stage_input_buffers_[0].push_block(input, data_size);

        for (int stage = 0; stage < RESAMPLING_STAGE_COUNT; ++stage)
        {
            input = downsampling_filters_[static_cast<std::size_t>(stage)].process_down(input);
            data_size /= 2;
            stage_input_buffers_[static_cast<std::size_t>(stage + 1)].push_block(input, data_size);
        }

        ++pending_output_blocks_;
        last_processed_input_size_ += processing_block_size_;
    }

    template <int StageCount>
    inline void ResamplingFilterbank<StageCount>::input_block(const double *const data, const int block_size)
    {
        if (block_size < 0 || block_size > maximum_callback_block_size_ || (data == nullptr && block_size > 0))
        {
            throw std::invalid_argument("The input block must not exceed the callback size supplied to init()");
        }

        last_processed_input_size_ = 0;
        int input_position = 0;
        while (input_position < block_size)
        {
            const int samples_to_copy =
                std::min(block_size - input_position, processing_block_size_ - input_data_size_);
            std::copy_n(data + input_position, samples_to_copy, input_data_.data() + input_data_size_);
            input_position += samples_to_copy;
            input_data_size_ += samples_to_copy;

            if (input_data_size_ == processing_block_size_)
            {
                process_input();
                input_data_size_ = 0;
            }
        }
    }

    template <int StageCount>
    inline void ResamplingFilterbank<StageCount>::push_output_samples(const double *const data, const int block_size)
    {
        if (static_cast<std::size_t>(block_size) > output_queue_.get_free_sample_count())
        {
            throw std::logic_error("Too many input blocks were submitted without consuming output");
        }
        output_queue_.push_block(data, block_size);
    }

    template <int StageCount>
    inline void ResamplingFilterbank<StageCount>::pull_output_samples(double *const data, const int block_size)
    {
        if (static_cast<std::size_t>(block_size) > output_queue_.get_available_sample_count())
        {
            throw std::logic_error("More output was requested than the filterbank has buffered");
        }
        output_queue_.pull_block(data, block_size);
    }

    template <int StageCount>
    inline void ResamplingFilterbank<StageCount>::process_output()
    {
        stage_output_buffers_[StageCount - 1].pull_block(lowest_stage_output_.data(), lowest_stage_block_size_);
        double *output = lowest_stage_output_.data();

        for (int stage = RESAMPLING_STAGE_COUNT - 1; stage >= 0; --stage)
        {
            output = upsampling_filters_[static_cast<std::size_t>(stage)].process_up(output);
            const int stage_block_size = origin_block_size_ / (1 << stage);
            stage_output_buffers_[static_cast<std::size_t>(stage)].pull_block_add(output, stage_block_size);
        }

        output = input_resampler_.process_up(output);
        push_output_samples(output, processing_block_size_);
    }

    template <int StageCount>
    inline double *ResamplingFilterbank<StageCount>::output_block(const int block_size)
    {
        if (block_size < 0 || block_size > maximum_callback_block_size_)
        {
            throw std::invalid_argument("The output block must not exceed the callback size supplied to init()");
        }

        while (pending_output_blocks_ > 0)
        {
            process_output();
            --pending_output_blocks_;
        }

        pull_output_samples(output_data_.data(), block_size);
        return output_data_.data();
    }

}
