#include "constant_q_transform.h"
#include "sliding_cqt.h"

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

    void require_finite(const double *const data, const int size, const std::string &message)
    {
        for (int sample = 0; sample < size; ++sample)
        {
            require(std::isfinite(data[sample]), message);
        }
    }

    std::vector<double>
    make_sine_block(const int block_size, const int block_index, const double sample_rate, const double frequency)
    {
        std::vector<double> block(static_cast<std::size_t>(block_size));
        for (int sample = 0; sample < block_size; ++sample)
        {
            const int position = block_index * block_size + sample;
            block[static_cast<std::size_t>(sample)] =
                std::sin(2. * audio_utils::pi<double>() * frequency * position / sample_rate);
        }
        return block;
    }

    double measure_tone_amplitude(const std::vector<double> &data,
                                  const std::size_t start,
                                  const double sample_rate,
                                  const double frequency)
    {
        double real = 0.;
        double imaginary = 0.;
        for (std::size_t sample = start; sample < data.size(); ++sample)
        {
            const double phase = 2. * audio_utils::pi<double>() * frequency * static_cast<double>(sample) / sample_rate;
            real += data[sample] * std::cos(phase);
            imaginary -= data[sample] * std::sin(phase);
        }
        const double sample_count = static_cast<double>(data.size() - start);
        return 2. * std::hypot(real, imaginary) / sample_count;
    }

    void test_sliding_cqt_batched_input()
    {
        constexpr int BLOCK_SIZE = 64;
        rt_cqt::SlidingCqt<12, 9, false> cqt;
        cqt.init(48000., BLOCK_SIZE);

        for (int block_index = 0; block_index < 12; ++block_index)
        {
            std::vector<double> input = make_sine_block(BLOCK_SIZE, block_index, 48000., 440.);
            cqt.input_block(input.data(), BLOCK_SIZE);

            const bool processing_boundary = ((block_index + 1) % 4) == 0;
            for (int stage = 0; stage < 9; ++stage)
            {
                const std::size_t expected_samples = processing_boundary ? static_cast<std::size_t>(256 >> stage) : 0U;
                require(cqt.get_samples_to_process(stage) == expected_samples,
                        "Sliding CQT received an unexpected stage block size");
            }

            double *output = cqt.output_block(BLOCK_SIZE);
            require_finite(output, BLOCK_SIZE, "Sliding CQT produced a non-finite sample");
        }
    }

    template <bool Windowing>
    void test_sliding_cqt_amplitude_normalization()
    {
        constexpr int BLOCK_SIZE = 64;
        constexpr int BLOCK_COUNT = 512;
        constexpr double SAMPLE_RATE = 48000.;
        constexpr int TONE = 8;
        const double frequency = rt_cqt::compute_bin_frequency(rt_cqt::compute_reference_frequency(440.), 24, 0, TONE);

        rt_cqt::SlidingCqt<24, 9, Windowing> cqt;
        cqt.init(SAMPLE_RATE, BLOCK_SIZE);

        std::vector<double> output(static_cast<std::size_t>(BLOCK_SIZE * BLOCK_COUNT), 0.);
        double magnitude_sum = 0.;
        int magnitude_count = 0;
        for (int block_index = 0; block_index < BLOCK_COUNT; ++block_index)
        {
            std::vector<double> input = make_sine_block(BLOCK_SIZE, block_index, SAMPLE_RATE, frequency);
            cqt.input_block(input.data(), BLOCK_SIZE);

            if (block_index >= (BLOCK_COUNT / 2) && cqt.get_samples_to_process(0) > 0)
            {
                magnitude_sum += std::abs(cqt.get_octave_cqt_buffer(0)[TONE].pull_delay_sample(0));
                ++magnitude_count;
            }

            const double *const output_block = cqt.output_block(BLOCK_SIZE);
            std::copy_n(output_block, BLOCK_SIZE, output.data() + static_cast<std::size_t>(block_index * BLOCK_SIZE));
        }

        const double coefficient_magnitude = magnitude_sum / static_cast<double>(magnitude_count);
        require(std::abs(coefficient_magnitude - 0.5) < 1.e-3,
                "Sliding CQT coefficient normalization is incorrect: " + std::to_string(coefficient_magnitude));

        const double reconstructed_amplitude =
            measure_tone_amplitude(output, output.size() / 2, SAMPLE_RATE, frequency);
        const double minimum_amplitude = Windowing ? 0.9 : 0.75;
        require(reconstructed_amplitude > minimum_amplitude && reconstructed_amplitude < 1.1,
                "Sliding CQT resynthesis normalization is incorrect: " + std::to_string(reconstructed_amplitude));
    }

    void test_constant_cqt_batched_schedule()
    {
        constexpr int BLOCK_SIZE = 64;
        rt_cqt::ConstantQTransform<12, 9> cqt;
        cqt.init(64);
        cqt.init_sample_rate(48000., BLOCK_SIZE);

        for (int block_index = 0; block_index < 12; ++block_index)
        {
            std::vector<double> input = make_sine_block(BLOCK_SIZE, block_index, 48000., 440.);
            cqt.input_block(input.data(), BLOCK_SIZE);

            const auto &schedule = cqt.get_cqt_schedule();
            const bool processing_boundary = ((block_index + 1) % 4) == 0;
            require(processing_boundary ? !schedule.empty() : schedule.empty(),
                    "Constant CQT schedule is not aligned with filterbank processing");

            std::array<int, 9> first_delay;
            std::array<int, 9> previous_synthesis_offset;
            first_delay.fill(-1);
            previous_synthesis_offset.fill(-1);
            for (const rt_cqt::ScheduleElement &element : schedule)
            {
                require(element.sample() >= 0 && element.sample() < 256,
                        "Schedule position is outside the internal processing block");
                const int octave = element.octave();
                if (first_delay[octave] < 0)
                {
                    first_delay[octave] = element.delay_at_octave_rate();
                }
                require(element.synthesis_offset() == first_delay[octave] - element.delay_at_octave_rate(),
                        "Synthesis offset is inconsistent with the analysis delay");
                require(element.synthesis_offset() > previous_synthesis_offset[octave],
                        "Synthesis offsets must advance within an internal block");
                previous_synthesis_offset[octave] = element.synthesis_offset();
                cqt.cqt(element);
                cqt.icqt(element);
            }

            double *output = cqt.output_block(BLOCK_SIZE);
            require_finite(output, BLOCK_SIZE, "Constant CQT produced a non-finite sample");
        }
    }

    void test_constant_cqt_high_frequency_resynthesis()
    {
        constexpr int BLOCK_SIZE = 64;
        constexpr int BLOCK_COUNT = 512;
        constexpr double SAMPLE_RATE = 48000.;
        const double frequency = rt_cqt::compute_bin_frequency(rt_cqt::compute_reference_frequency(440.), 24, 0, 8);

        rt_cqt::ConstantQTransform<24, 9> cqt;
        cqt.init(64);
        cqt.init_sample_rate(SAMPLE_RATE, BLOCK_SIZE);

        std::vector<double> output(static_cast<std::size_t>(BLOCK_SIZE * BLOCK_COUNT), 0.);
        for (int block_index = 0; block_index < BLOCK_COUNT; ++block_index)
        {
            std::vector<double> input = make_sine_block(BLOCK_SIZE, block_index, SAMPLE_RATE, frequency);
            cqt.input_block(input.data(), BLOCK_SIZE);
            for (const rt_cqt::ScheduleElement &element : cqt.get_cqt_schedule())
            {
                cqt.cqt(element);
                cqt.icqt(element);
            }

            const double *const output_block = cqt.output_block(BLOCK_SIZE);
            std::copy_n(output_block, BLOCK_SIZE, output.data() + static_cast<std::size_t>(block_index * BLOCK_SIZE));
        }

        const double reconstructed_amplitude =
            measure_tone_amplitude(output, output.size() / 2, SAMPLE_RATE, frequency);
        require(reconstructed_amplitude > 0.5,
                "Constant CQT lost high-frequency resynthesis amplitude: " + std::to_string(reconstructed_amplitude));
    }

}

int main()
{
    try
    {
        test_sliding_cqt_batched_input();
        std::cout << "[pass] sliding CQT batched input\n";
        test_sliding_cqt_amplitude_normalization<false>();
        std::cout << "[pass] rectangular sliding CQT amplitude normalization\n";
        test_sliding_cqt_amplitude_normalization<true>();
        std::cout << "[pass] windowed sliding CQT amplitude normalization\n";
        test_constant_cqt_batched_schedule();
        std::cout << "[pass] constant CQT batched schedule\n";
        test_constant_cqt_high_frequency_resynthesis();
        std::cout << "[pass] constant CQT high-frequency resynthesis\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << "[fail] " << error.what() << '\n';
        return 1;
    }

    return 0;
}
