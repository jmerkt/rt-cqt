/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

#include "../submodules/audio-utils/include/utils.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace rt_cqt
{

    template <typename FloatType>
    class Delay
    {
    public:
        Delay() = default;
        ~Delay() = default;
        inline void reset() { storage_ = static_cast<FloatType>(0.); }

        inline void process(FloatType *const data, const int block_size)
        {
            for (int i = 0; i < block_size; i++)
            {
                const FloatType input = data[i];
                data[i] = storage_;
                storage_ = input;
            }
        }

    private:
        FloatType storage_{0.};
    };

    template <typename FloatType>
    class FirstOrderAllpass
    {
    public:
        FirstOrderAllpass() = default;
        ~FirstOrderAllpass() = default;

        inline void init_coefficient(FloatType coefficient)
        {
            coefficient_ = static_cast<FloatType>(coefficient);
            previous_input_ = static_cast<FloatType>(0.);
            previous_output_ = static_cast<FloatType>(0.);
        };

        inline void process(FloatType *const samples, const int block_size)
        {
            for (int i = 0; i < block_size; i++)
            {
                const FloatType sample = samples[i];
                samples[i] = coefficient_ * (sample - previous_output_) + previous_input_;
                previous_input_ = sample;
                previous_output_ = samples[i];
            }
        };

    private:
        FloatType coefficient_{0.};
        FloatType previous_input_{0.};
        FloatType previous_output_{0.};
    };

    /*
     * Polyphase IIR lowpass for resampling blocks by a factor of two.
     *
     * See "Digital signal processing schemes for efficient interpolation and
     * decimation" by Reinaldo A. Valenzuela and A. G. Constantinides. This
     * class handles allocation and downsampling or upsampling for incoming
     * blocks. A single instance must only process in one direction.
     */
    template <typename FloatType, std::size_t AllpassCount>
    class HalfBandLowpass
    {
    public:
        HalfBandLowpass();
        ~HalfBandLowpass() = default;

        /*
        Initialize the filters and allocate memory. Has to get called before processing starts.
        */
        bool init(const int expected_block_size = 128, bool downsampling = true, double transition_bandwidth = 0.02);

        FloatType *process_down(const FloatType *const input_block);
        FloatType *process_up(const FloatType *const input_block);

        int get_output_block_size() { return output_block_size_; };
        int get_input_block_size() { return input_block_size_; };

    private:
        double transition_bandwidth_;
        std::size_t allpass_count_;
        std::size_t filter_order_;
        std::vector<double> coefficients_;

        FirstOrderAllpass<FloatType> direct_path_filters_[AllpassCount];
        FirstOrderAllpass<FloatType> delay_path_filters_[AllpassCount];
        Delay<FloatType> delay_;

        int output_block_size_{0};
        int input_block_size_{1};
        int filter_block_size_{0};

        std::vector<FloatType> direct_path_buffer_;
        std::vector<FloatType> delay_path_buffer_;
        std::vector<FloatType> output_block_;

        std::vector<double> design_filter();
    };

    template <typename FloatType, std::size_t AllpassCount>
    HalfBandLowpass<FloatType, AllpassCount>::HalfBandLowpass()
    {
        allpass_count_ = AllpassCount * 2;
        filter_order_ = 2 * allpass_count_ + 1;
    }

    template <typename FloatType, std::size_t AllpassCount>
    inline bool HalfBandLowpass<FloatType, AllpassCount>::init(const int expected_block_size,
                                                               bool downsampling,
                                                               double transition_bandwidth)
    {
        if (expected_block_size <= 0 || (downsampling && ((expected_block_size % 2) != 0)))
        {
            throw std::invalid_argument(
                "Resampling blocks must be positive; downsampling blocks must also be divisible by two");
        }

        // init filters
        transition_bandwidth_ = transition_bandwidth * 2. * audio_utils::pi<double>();
        coefficients_.clear();
        coefficients_ = design_filter();
        int filter_index = 0;
        for (std::size_t i = 0; i < coefficients_.size(); i += 2)
        {
            direct_path_filters_[filter_index].init_coefficient(coefficients_[i]);
            filter_index++;
        }
        filter_index = 0;
        for (std::size_t i = 1; i < coefficients_.size(); i += 2)
        {
            delay_path_filters_[filter_index].init_coefficient(coefficients_[i]);
            filter_index++;
        }
        // init buffers
        direct_path_buffer_.clear();
        delay_path_buffer_.clear();
        output_block_.clear();
        delay_.reset();
        input_block_size_ = expected_block_size;
        if (downsampling)
        {
            filter_block_size_ = expected_block_size / 2;
            output_block_size_ = expected_block_size / 2;
        }
        else
        {
            filter_block_size_ = expected_block_size;
            output_block_size_ = expected_block_size * 2;
        }
        direct_path_buffer_.resize(filter_block_size_, static_cast<FloatType>(0.));
        delay_path_buffer_.resize(filter_block_size_, static_cast<FloatType>(0.));
        output_block_.resize(output_block_size_, static_cast<FloatType>(0.));
        return true;
    }

    template <typename FloatType, std::size_t AllpassCount>
    inline FloatType *HalfBandLowpass<FloatType, AllpassCount>::process_down(const FloatType *const input_block)
    {
        int direct_output_index = 0;
        for (int i = 0; i < input_block_size_; i += 2)
        {
            direct_path_buffer_[direct_output_index] = input_block[i];
            direct_output_index++;
        }
        int delay_output_index = 0;
        for (int i = 1; i < input_block_size_; i += 2)
        {
            delay_path_buffer_[delay_output_index] = input_block[i];
            delay_output_index++;
        }
        delay_.process(delay_path_buffer_.data(), filter_block_size_);
        for (std::size_t i = 0; i < AllpassCount; i++)
        {
            direct_path_filters_[i].process(direct_path_buffer_.data(), filter_block_size_);
            delay_path_filters_[i].process(delay_path_buffer_.data(), filter_block_size_);
        }
        for (int i = 0; i < output_block_size_; i++)
        {
            output_block_[i] = static_cast<FloatType>(0.5) * (direct_path_buffer_[i] + delay_path_buffer_[i]);
        }
        return output_block_.data();
    };

    template <typename FloatType, std::size_t AllpassCount>
    inline FloatType *HalfBandLowpass<FloatType, AllpassCount>::process_up(const FloatType *const input_block)
    {
        for (int i = 0; i < input_block_size_; i++)
        {
            direct_path_buffer_[i] = input_block[i];
            delay_path_buffer_[i] = input_block[i];
        }
        for (std::size_t i = 0; i < AllpassCount; i++)
        {
            direct_path_filters_[i].process(direct_path_buffer_.data(), filter_block_size_);
            delay_path_filters_[i].process(delay_path_buffer_.data(), filter_block_size_);
        }
        int direct_input_index = 0;
        for (int i = 0; i < output_block_size_; i += 2)
        {
            output_block_[i] = direct_path_buffer_[direct_input_index];
            direct_input_index += 1;
        }
        int delay_input_index = 0;
        for (int i = 1; i < output_block_size_; i += 2)
        {
            output_block_[i] = delay_path_buffer_[delay_input_index];
            delay_input_index += 1;
        }
        return output_block_.data();
    };

    template <typename FloatType, std::size_t AllpassCount>
    inline std::vector<double> HalfBandLowpass<FloatType, AllpassCount>::design_filter()
    {
        // step 1
        const double k = std::pow(std::tan((audio_utils::pi<double>() - transition_bandwidth_) / 4.), 2);
        const double k_dash = std::sqrt(1. - std::pow(k, 2));
        const double e = (1. / 2.) * ((1. - std::sqrt(k_dash)) / (1. + std::sqrt(k_dash)));
        const double q = e + 2. * std::pow(e, 5) + 15. * std::pow(e, 9.) + 150. * std::pow(e, 13.);
        // step 2
        const std::size_t n = filter_order_;
        // step 3
        std::vector<double> w;
        std::vector<double> a_dash;
        for (std::size_t i = 1; i <= ((n - 1) / 2); i++)
        {
            // w_i
            double delta = 1.;
            double num = 0.;
            double m = 0.;
            while (delta > 1.e-100)
            {
                delta = std::pow((-1.), m) * std::pow(q, (m * (m + 1.))) *
                        std::sin((2. * m + 1.) * audio_utils::pi<double>() * static_cast<double>(i) /
                                 static_cast<double>(n));
                num += delta;
                m += 1.;
            }
            num = 2. * std::pow(q, (1. / 4.)) * num;
            delta = 1.;
            double den = 0.;
            m = 1.;
            while (delta > 1.e-100)
            {
                delta = std::pow((-1.), m) * std::pow(q, std::pow(m, 2.)) *
                        std::cos(2. * m * audio_utils::pi<double>() * static_cast<double>(i) / static_cast<double>(n));
                den += delta;
                m += 1.;
            }
            den = den * 2. + 1.;
            double w_i = num / den;
            w.push_back(w_i);
            // a'_i
            num = std::pow(((1. - std::pow(w_i, 2.) * k) * (1. - std::pow(w_i, 2.) / k)), (1. / 2.));
            den = 1. + std::pow(w_i, 2.);
            double a_dash_i = num / den;
            a_dash.push_back(a_dash_i);
        }
        // step 4
        std::vector<double> a;
        for (double a_dash_i : a_dash)
        {
            double a_i = (1. - a_dash_i) / (1. + a_dash_i);
            a.push_back(a_i);
        }
        return a;
    };

    enum class DirectionConfig
    {
        Up = 0,
        Down,
        UpDown,
        DownUp
    };

    /**
    Resampling of audio blocks by a given power of 2.
    */
    template <typename FloatType, std::size_t AllpassCount>
    class ResamplingHandler
    {
    public:
        ResamplingHandler(double transition_bandwidth = 0.02);
        ~ResamplingHandler() = default;
        /**
        Initialization of the ResamplingHandler
        The configuration gives the amount of memory that will be allocated.
        DirectionConfig has to match the order the processing functions will be called from outside.
        */
        void init(const int power_of_two_exponent = 0,
                  const int expected_block_size = 128,
                  DirectionConfig direction = DirectionConfig::UpDown);

        FloatType *process_down(const FloatType *const input_block);
        FloatType *process_up(const FloatType *const input_block);

        int get_upsampled_block_size() { return upsampled_block_size_; };
        int get_downsampled_block_size() { return downsampled_block_size_; };

    private:
        double transition_bandwidth_;
        int power_of_two_exponent_{0};
        int upsampled_block_size_{0};
        int downsampled_block_size_{0};
        std::vector<HalfBandLowpass<FloatType, AllpassCount>> downsampling_filters_;
        std::vector<HalfBandLowpass<FloatType, AllpassCount>> upsampling_filters_;
    };

    template <typename FloatType, std::size_t AllpassCount>
    ResamplingHandler<FloatType, AllpassCount>::ResamplingHandler(double transition_bandwidth)
    {
        transition_bandwidth_ = transition_bandwidth;
    }

    template <typename FloatType, std::size_t AllpassCount>
    inline void ResamplingHandler<FloatType, AllpassCount>::init(const int power_of_two_exponent,
                                                                 const int expected_block_size,
                                                                 DirectionConfig direction)
    {
        if (power_of_two_exponent < 0 || power_of_two_exponent >= 30)
        {
            throw std::invalid_argument("The power-of-two exponent must be between 0 and 29");
        }
        if (expected_block_size <= 0)
        {
            throw std::invalid_argument("The resampling block size must be positive");
        }

        downsampling_filters_.clear();
        upsampling_filters_.clear();
        downsampling_filters_.reserve(static_cast<std::size_t>(power_of_two_exponent));
        upsampling_filters_.reserve(static_cast<std::size_t>(power_of_two_exponent));
        power_of_two_exponent_ = power_of_two_exponent;
        downsampled_block_size_ = expected_block_size;
        upsampled_block_size_ = expected_block_size;

        const int resampling_factor = 1 << power_of_two_exponent_;
        if ((direction == DirectionConfig::Down) || (direction == DirectionConfig::DownUp))
        {
            if ((expected_block_size % resampling_factor) != 0)
            {
                throw std::invalid_argument("The downsampling block size must be divisible by the resampling factor");
            }
        }
        if ((direction == DirectionConfig::Up || direction == DirectionConfig::UpDown) &&
            expected_block_size > (std::numeric_limits<int>::max() / resampling_factor))
        {
            throw std::overflow_error("The upsampling block size is too large");
        }

        for (int i = 0; i < power_of_two_exponent_; i++)
        {
            const int stage_factor = 1 << i;
            if (direction == DirectionConfig::Down)
            {
                downsampling_filters_.emplace_back();
                downsampling_filters_.back().init(expected_block_size / stage_factor, true, transition_bandwidth_);
            }
            else if (direction == DirectionConfig::Up)
            {
                upsampling_filters_.emplace_back();
                upsampling_filters_.back().init(expected_block_size * stage_factor, false, transition_bandwidth_);
            }
            else if (direction == DirectionConfig::DownUp)
            {
                downsampling_filters_.emplace_back();
                upsampling_filters_.emplace_back();
                downsampling_filters_.back().init(expected_block_size / stage_factor, true, transition_bandwidth_);
                upsampling_filters_.back().init(
                    expected_block_size / (1 << (power_of_two_exponent_ - i)), false, transition_bandwidth_);
            }
            else if (direction == DirectionConfig::UpDown)
            {
                upsampling_filters_.emplace_back();
                downsampling_filters_.emplace_back();
                upsampling_filters_.back().init(expected_block_size * stage_factor, false, transition_bandwidth_);
                downsampling_filters_.back().init(
                    expected_block_size * (1 << (power_of_two_exponent_ - i)), true, transition_bandwidth_);
            }
        }

        if (!downsampling_filters_.empty())
        {
            downsampled_block_size_ = downsampling_filters_.back().get_output_block_size();
        }
        if (!upsampling_filters_.empty())
        {
            upsampled_block_size_ = upsampling_filters_.back().get_output_block_size();
        }
    };

    template <typename FloatType, std::size_t AllpassCount>
    inline FloatType *ResamplingHandler<FloatType, AllpassCount>::process_down(const FloatType *const input_block)
    {
        assert(power_of_two_exponent_ == 0 || static_cast<int>(downsampling_filters_.size()) == power_of_two_exponent_);
        const FloatType *input = input_block;
        FloatType *output = const_cast<FloatType *>(input_block);
        for (int i = 0; i < power_of_two_exponent_; i++)
        {
            output = downsampling_filters_[i].process_down(input);
            input = output;
        }
        return output;
    };

    template <typename FloatType, std::size_t AllpassCount>
    inline FloatType *ResamplingHandler<FloatType, AllpassCount>::process_up(const FloatType *const input_block)
    {
        assert(power_of_two_exponent_ == 0 || static_cast<int>(upsampling_filters_.size()) == power_of_two_exponent_);
        const FloatType *input = input_block;
        FloatType *output = const_cast<FloatType *>(input_block);
        for (int i = 0; i < power_of_two_exponent_; i++)
        {
            output = upsampling_filters_[i].process_up(input);
            input = output;
        }
        return output;
    };

}
