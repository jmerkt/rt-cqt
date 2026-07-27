#include "../python-bindings/include/python_resampling_filterbank.h"
#include "../submodules/audio-utils/include/circular_buffer.h"
#include "resampling.h"
#include "resampling_filterbank.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

    constexpr unsigned TEST_ALLPASS_COUNT{3};
    constexpr double TEST_TRANSITION_BANDWIDTH{0.1};

    void require(const bool condition, const std::string &message)
    {
        if (!condition)
        {
            throw std::runtime_error(message);
        }
    }

    void require_near(const std::vector<double> &actual,
                      const std::vector<double> &expected,
                      const double tolerance,
                      const std::string &message)
    {
        require(actual.size() == expected.size(), message + ": vector sizes differ");
        double maximum_error = 0.;
        for (std::size_t index = 0; index < actual.size(); ++index)
        {
            maximum_error = std::max(maximum_error, std::abs(actual[index] - expected[index]));
        }
        require(maximum_error <= tolerance, message + ": maximum error is " + std::to_string(maximum_error));
    }

    std::vector<double> make_noise(const std::size_t size)
    {
        std::mt19937 generator(0x5eedU);
        std::uniform_real_distribution<double> distribution(-1., 1.);
        std::vector<double> result(size);
        for (double &sample : result)
        {
            sample = distribution(generator);
        }
        return result;
    }

    void test_circular_buffer_states()
    {
        audio_utils::CircularBuffer<int> buffer(4);
        require(buffer.get_buffer_size() == 4, "Circular buffer has an unexpected capacity");
        require(buffer.is_empty(), "A new circular buffer is not empty");
        require(!buffer.is_full(), "A new circular buffer is full");
        require(buffer.get_available_sample_count() == 0, "A new circular buffer contains samples");
        require(buffer.get_free_sample_count() == 4, "A new circular buffer reports the wrong free space");

        const std::array<int, 4> initial{1, 2, 3, 4};
        buffer.push_block(initial.data(), static_cast<int>(initial.size()));
        require(buffer.is_full(), "A capacity-sized write did not make the circular buffer full");
        require(!buffer.is_empty(), "A full circular buffer is also reported as empty");
        require(buffer.get_available_sample_count() == 4, "A full circular buffer reports zero available samples");
        require(buffer.get_free_sample_count() == 0, "A full circular buffer reports free space");

        buffer.push_sample(5);
        require(buffer.is_full(), "Overwriting the oldest sample changed the full state");
        require(buffer.pull_delay_sample(0) == 5, "The newest delayed sample was not retained");
        require(buffer.pull_delay_sample(3) == 2, "The oldest retained delayed sample is incorrect");

        std::array<int, 4> overwritten{};
        buffer.pull_block(overwritten.data(), static_cast<int>(overwritten.size()));
        require(overwritten == std::array<int, 4>{2, 3, 4, 5},
                "A full circular buffer did not discard exactly its oldest sample");
        require(buffer.is_empty(), "Consuming every sample did not empty the circular buffer");

        const std::array<int, 3> wrapped_input{10, 11, 12};
        buffer.push_block(wrapped_input.data(), static_cast<int>(wrapped_input.size()));
        require(buffer.pull_sample() == 10, "Circular buffer returned the wrong sample before wrapping");
        const std::array<int, 2> wrapped_tail{13, 14};
        buffer.push_block(wrapped_tail.data(), static_cast<int>(wrapped_tail.size()));

        std::array<int, 4> wrapped_output{};
        buffer.pull_block(wrapped_output.data(), static_cast<int>(wrapped_output.size()));
        require(wrapped_output == std::array<int, 4>{11, 12, 13, 14},
                "Circular buffer order changed across pointer wraparound");

        buffer.push_sample(21);
        std::array<int, 2> unchanged{7, 8};
        bool block_underflow_rejected = false;
        try
        {
            buffer.pull_block(unchanged.data(), static_cast<int>(unchanged.size()));
        }
        catch (const std::underflow_error &)
        {
            block_underflow_rejected = true;
        }
        require(block_underflow_rejected, "Circular buffer accepted a block underflow");
        require(buffer.get_available_sample_count() == 1, "A rejected block underflow consumed data");
        require(unchanged == std::array<int, 2>{7, 8}, "A rejected block underflow modified its destination");
        require(buffer.pull_sample() == 21, "A rejected block underflow corrupted the buffered sample");

        bool sample_underflow_rejected = false;
        try
        {
            static_cast<void>(buffer.pull_sample());
        }
        catch (const std::underflow_error &)
        {
            sample_underflow_rejected = true;
        }
        require(sample_underflow_rejected, "Circular buffer accepted a sample underflow");

        buffer.push_sample(42);
        buffer.clear();
        require(buffer.is_empty(), "Clearing the circular buffer did not reset its logical state");

        bool zero_size_rejected = false;
        try
        {
            buffer.change_size(0);
        }
        catch (const std::invalid_argument &)
        {
            zero_size_rejected = true;
        }
        require(zero_size_rejected, "Circular buffer accepted a zero-sized allocation");
    }

    void test_half_band_partition_invariance()
    {
        const std::vector<double> input = make_noise(1024);

        rt_cqt::HalfBandLowpass<double, TEST_ALLPASS_COUNT> whole_downsampler;
        whole_downsampler.init(1024, true, TEST_TRANSITION_BANDWIDTH);
        double *whole_down_data = whole_downsampler.process_down(input.data());
        const std::vector<double> whole_down(whole_down_data, whole_down_data + 512);

        rt_cqt::HalfBandLowpass<double, TEST_ALLPASS_COUNT> partitioned_downsampler;
        partitioned_downsampler.init(256, true, TEST_TRANSITION_BANDWIDTH);
        std::vector<double> partitioned_down;
        partitioned_down.reserve(512);
        for (std::size_t offset = 0; offset < input.size(); offset += 256)
        {
            double *block = partitioned_downsampler.process_down(input.data() + offset);
            partitioned_down.insert(partitioned_down.end(), block, block + 128);
        }
        require_near(partitioned_down, whole_down, 1.e-14, "Half-band downsampling depends on block partition");

        rt_cqt::HalfBandLowpass<double, TEST_ALLPASS_COUNT> whole_upsampler;
        whole_upsampler.init(512, false, TEST_TRANSITION_BANDWIDTH);
        double *whole_up_data = whole_upsampler.process_up(whole_down.data());
        const std::vector<double> whole_up(whole_up_data, whole_up_data + 1024);

        rt_cqt::HalfBandLowpass<double, TEST_ALLPASS_COUNT> partitioned_upsampler;
        partitioned_upsampler.init(128, false, TEST_TRANSITION_BANDWIDTH);
        std::vector<double> partitioned_up;
        partitioned_up.reserve(1024);
        for (std::size_t offset = 0; offset < whole_down.size(); offset += 128)
        {
            double *block = partitioned_upsampler.process_up(whole_down.data() + offset);
            partitioned_up.insert(partitioned_up.end(), block, block + 256);
        }
        require_near(partitioned_up, whole_up, 1.e-14, "Half-band upsampling depends on block partition");

        partitioned_downsampler.init(256, true, TEST_TRANSITION_BANDWIDTH);
        double *reinitialized_data = partitioned_downsampler.process_down(input.data());
        rt_cqt::HalfBandLowpass<double, TEST_ALLPASS_COUNT> fresh_downsampler;
        fresh_downsampler.init(256, true, TEST_TRANSITION_BANDWIDTH);
        double *fresh_data = fresh_downsampler.process_down(input.data());
        require_near(std::vector<double>(reinitialized_data, reinitialized_data + 128),
                     std::vector<double>(fresh_data, fresh_data + 128),
                     1.e-14,
                     "Reinitializing a half-band filter did not reset its state");
    }

    void test_resampling_handler_partition_invariance()
    {
        const std::vector<double> input = make_noise(1024);

        rt_cqt::ResamplingHandler<double, TEST_ALLPASS_COUNT> whole;
        whole.init(3, 1024, rt_cqt::DirectionConfig::Down);
        double *whole_data = whole.process_down(input.data());
        const std::vector<double> expected(whole_data, whole_data + 128);

        rt_cqt::ResamplingHandler<double, TEST_ALLPASS_COUNT> partitioned;
        partitioned.init(3, 256, rt_cqt::DirectionConfig::Down);
        std::vector<double> actual;
        actual.reserve(128);
        for (std::size_t offset = 0; offset < input.size(); offset += 256)
        {
            double *block = partitioned.process_down(input.data() + offset);
            actual.insert(actual.end(), block, block + 32);
        }

        require_near(actual, expected, 1.e-14, "Power-of-two resampling depends on block partition");

        rt_cqt::ResamplingHandler<double, TEST_ALLPASS_COUNT> whole_round_trip;
        whole_round_trip.init(3, 1024, rt_cqt::DirectionConfig::DownUp);
        double *whole_down = whole_round_trip.process_down(input.data());
        double *whole_up = whole_round_trip.process_up(whole_down);
        const std::vector<double> expected_round_trip(whole_up, whole_up + 1024);

        rt_cqt::ResamplingHandler<double, TEST_ALLPASS_COUNT> partitioned_round_trip;
        partitioned_round_trip.init(3, 256, rt_cqt::DirectionConfig::DownUp);
        std::vector<double> actual_round_trip;
        actual_round_trip.reserve(1024);
        for (std::size_t offset = 0; offset < input.size(); offset += 256)
        {
            double *down = partitioned_round_trip.process_down(input.data() + offset);
            double *up = partitioned_round_trip.process_up(down);
            actual_round_trip.insert(actual_round_trip.end(), up, up + 256);
        }
        require_near(actual_round_trip,
                     expected_round_trip,
                     1.e-14,
                     "Power-of-two round-trip resampling depends on block partition");
    }

    template <int StageCount>
    struct FilterbankRun
    {
        std::vector<std::vector<double>> stages_{static_cast<std::size_t>(StageCount)};
        std::vector<double> output_;
    };

    template <int StageCount>
    FilterbankRun<StageCount> run_filterbank(const double sample_rate,
                                             const int maximum_callback_size,
                                             const std::vector<double> &input,
                                             const std::vector<int> &partitions,
                                             const bool synthesize)
    {
        rt_cqt::ResamplingFilterbank<StageCount> filterbank;
        filterbank.init(sample_rate, maximum_callback_size, static_cast<int>(input.size() * 4));

        FilterbankRun<StageCount> result;
        result.output_.reserve(input.size());

        std::size_t input_position = 0;
        std::size_t partition_index = 0;
        while (input_position < input.size())
        {
            const int requested_size = partitions[partition_index % partitions.size()];
            const int block_size = std::min<int>(requested_size, static_cast<int>(input.size() - input_position));
            require(block_size <= maximum_callback_size, "Test partition exceeds configured callback size");

            filterbank.input_block(input.data() + input_position, block_size);

            for (int stage = 0; stage < StageCount; ++stage)
            {
                rt_cqt::BufferPtr input_buffer = filterbank.get_stage_input_buffer(stage);
                const int stage_block_size = static_cast<int>(input_buffer->get_available_sample_count());
                if (stage_block_size == 0)
                {
                    continue;
                }

                std::vector<double> stage_data(static_cast<std::size_t>(stage_block_size));
                input_buffer->pull_block(stage_data.data(), stage_block_size);
                result.stages_[static_cast<std::size_t>(stage)].insert(
                    result.stages_[static_cast<std::size_t>(stage)].end(), stage_data.begin(), stage_data.end());

                if (synthesize)
                {
                    filterbank.get_stage_output_buffer(stage)->push_block(stage_data.data(), stage_block_size);
                }
            }

            if (synthesize)
            {
                double *output_block = filterbank.output_block(block_size);
                result.output_.insert(result.output_.end(), output_block, output_block + block_size);
            }

            input_position += static_cast<std::size_t>(block_size);
            ++partition_index;
        }

        return result;
    }

    void test_filterbank_block_sizing()
    {
        rt_cqt::ResamplingFilterbank<1> single_stage;
        single_stage.init(48000., 3, 16);
        require(single_stage.get_processing_block_size() == 3, "Single-stage filterbank was unnecessarily aligned");
        require(single_stage.get_origin_block_size() == 3, "Unexpected single-stage origin block");
        const auto single_stage_result = run_filterbank<1>(48000., 3, make_noise(12), {1, 2, 3}, true);
        require(single_stage_result.stages_[0].size() == 12, "Single-stage analysis lost samples");
        require(single_stage_result.output_.size() == 12, "Single-stage synthesis lost samples");

        rt_cqt::ResamplingFilterbank<9> filterbank_48;
        filterbank_48.init(48000., 64, 1024);
        require(filterbank_48.get_processing_block_size() == 256, "Unexpected 48 kHz processing block");
        require(filterbank_48.get_origin_block_size() == 256, "Unexpected 48 kHz origin block");
        require(filterbank_48.get_latency_samples() == 256, "Unexpected 48 kHz latency");

        rt_cqt::ResamplingFilterbank<9> filterbank_44;
        filterbank_44.init(44100., 64, 1024);
        require(filterbank_44.get_processing_block_size() == 256, "Unexpected 44.1 kHz processing block");
        require(filterbank_44.get_origin_block_size() == 256, "Unexpected 44.1 kHz origin block");

        rt_cqt::ResamplingFilterbank<9> filterbank_96;
        filterbank_96.init(96000., 64, 1024);
        require(filterbank_96.get_processing_block_size() == 512, "Unexpected 96 kHz processing block");
        require(filterbank_96.get_origin_block_size() == 256, "Unexpected 96 kHz origin block");

        rt_cqt::ResamplingFilterbank<9> filterbank_88;
        filterbank_88.init(88200., 64, 1024);
        require(filterbank_88.get_processing_block_size() == 512, "Unexpected 88.2 kHz processing block");
        require(filterbank_88.get_origin_block_size() == 256, "Unexpected 88.2 kHz origin block");

        rt_cqt::ResamplingFilterbank<9> large_callback_filterbank;
        large_callback_filterbank.init(48000., 1000, 4096);
        require(large_callback_filterbank.get_processing_block_size() == 1024, "Callback was not aligned upward");

        bool rejected_unsupported_rate = false;
        try
        {
            rt_cqt::ResamplingFilterbank<9> unsupported;
            unsupported.init(88201., 64, 1024);
        }
        catch (const std::invalid_argument &)
        {
            rejected_unsupported_rate = true;
        }
        require(rejected_unsupported_rate, "Unsupported sample rate was accepted");
    }

    void test_filterbank_batching()
    {
        rt_cqt::ResamplingFilterbank<9> filterbank;
        filterbank.init(48000., 64, 1024);
        std::vector<double> block(64, 1.);

        for (int callback = 0; callback < 3; ++callback)
        {
            filterbank.input_block(block.data(), static_cast<int>(block.size()));
            require(filterbank.get_last_processed_input_size() == 0, "Filterbank released a partial internal block");
        }

        filterbank.input_block(block.data(), static_cast<int>(block.size()));
        require(filterbank.get_last_processed_input_size() == 256,
                "Filterbank did not release a complete internal block");
        for (int stage = 0; stage < 9; ++stage)
        {
            const std::size_t expected_size = static_cast<std::size_t>(256 >> stage);
            require(filterbank.get_stage_input_buffer(stage)->get_available_sample_count() == expected_size,
                    "Unexpected stage block size at stage " + std::to_string(stage));
        }
    }

    void test_callback_size_matrix()
    {
        const std::vector<std::pair<int, int>> callback_and_processing_sizes{
            {1, 256},
            {3, 256},
            {64, 256},
            {257, 512},
            {1024, 1024},
        };

        for (const auto &[callback_size, processing_size] : callback_and_processing_sizes)
        {
            rt_cqt::ResamplingFilterbank<9> filterbank;
            filterbank.init(48000., callback_size, processing_size * 4);
            require(filterbank.get_processing_block_size() == processing_size,
                    "Unexpected processing size for callback size " + std::to_string(callback_size));

            const std::vector<double> input = make_noise(static_cast<std::size_t>(processing_size * 2));
            const auto result = run_filterbank<9>(48000., callback_size, input, {callback_size}, true);
            for (int stage = 0; stage < 9; ++stage)
            {
                require(result.stages_[static_cast<std::size_t>(stage)].size() == (input.size() >> stage),
                        "Callback-size matrix lost stage samples");
            }
            require(std::all_of(result.output_.begin(),
                                result.output_.end(),
                                [](const double sample) { return std::isfinite(sample); }),
                    "Callback-size matrix produced a non-finite output");
        }
    }

    void test_filterbank_partition_invariance()
    {
        const std::vector<double> input = make_noise(512);
        const auto regular = run_filterbank<5>(48000., 13, input, {13}, true);
        const auto irregular = run_filterbank<5>(48000., 13, input, {1, 7, 3, 11, 5, 2, 13}, true);

        for (int stage = 0; stage < 5; ++stage)
        {
            const std::size_t expected_size = input.size() >> stage;
            require(regular.stages_[static_cast<std::size_t>(stage)].size() == expected_size,
                    "Unexpected accumulated stage size");
            require_near(irregular.stages_[static_cast<std::size_t>(stage)],
                         regular.stages_[static_cast<std::size_t>(stage)],
                         1.e-14,
                         "Filterbank analysis depends on callback partition");
        }

        require_near(irregular.output_, regular.output_, 1.e-14, "Filterbank synthesis depends on callback partition");
        require(std::all_of(regular.output_.begin(),
                            regular.output_.begin() + 16,
                            [](const double sample) { return sample == 0.; }),
                "The declared one-block filterbank latency is missing");
        require(std::any_of(regular.output_.begin() + 16,
                            regular.output_.end(),
                            [](const double sample) { return std::abs(sample) > 1.e-12; }),
                "Filterbank synthesis produced only silence");
    }

    void test_filterbank_at_96_khz()
    {
        const std::vector<double> input = make_noise(512);
        const auto regular = run_filterbank<5>(96000., 13, input, {13}, true);
        const auto irregular = run_filterbank<5>(96000., 13, input, {13, 5, 1, 9}, true);
        for (int stage = 0; stage < 5; ++stage)
        {
            const std::size_t expected_size = input.size() >> (stage + 1);
            require(regular.stages_[static_cast<std::size_t>(stage)].size() == expected_size,
                    "Unexpected 96 kHz stage size");
            require_near(irregular.stages_[static_cast<std::size_t>(stage)],
                         regular.stages_[static_cast<std::size_t>(stage)],
                         1.e-14,
                         "96 kHz filterbank analysis depends on callback partition");
        }
        require_near(
            irregular.output_, regular.output_, 1.e-14, "96 kHz filterbank synthesis depends on callback partition");
    }

    void test_python_plot_adapter()
    {
        rt_cqt::PythonResamplingFilterbank<9> filterbank;
        filterbank.init(48000., 64);
        std::vector<double> block(64, 0.);
        block[0] = 1.;

        for (int callback = 0; callback < 4; ++callback)
        {
            const auto result = filterbank.process(block);
            require(result.second.size() == block.size(), "Python adapter returned the wrong output size");
            const bool processing_boundary = callback == 3;
            require(result.first[0].size() == (processing_boundary ? 256U : 0U),
                    "Python adapter exposed the wrong stage block size");
            std::fill(block.begin(), block.end(), 0.);
        }
    }

    void run(const std::string &name, const std::function<void()> &test)
    {
        test();
        std::cout << "[pass] " << name << '\n';
    }

}

int main()
{
    try
    {
        run("circular-buffer full and empty states", test_circular_buffer_states);
        run("half-band partition invariance", test_half_band_partition_invariance);
        run("resampling-handler partition invariance", test_resampling_handler_partition_invariance);
        run("filterbank block sizing", test_filterbank_block_sizing);
        run("filterbank batching", test_filterbank_batching);
        run("callback-size matrix", test_callback_size_matrix);
        run("filterbank partition invariance", test_filterbank_partition_invariance);
        run("filterbank at 96 kHz", test_filterbank_at_96_khz);
        run("Python plot adapter", test_python_plot_adapter);
    }
    catch (const std::exception &error)
    {
        std::cerr << "[fail] " << error.what() << '\n';
        return 1;
    }

    return 0;
}
