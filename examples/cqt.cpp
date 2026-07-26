#include "constant_q_transform.h"

#include <cstddef>

int main(int argc, char *argv[])
{
    const int hop_size = 256;
    const int octave_number = 9;
    const int bins_per_octave = 12;

    const int block_size = 1024;
    const double sample_rate = 48000.;

    std::vector<double> audio_input_block(block_size, 0.);
    std::vector<double> audio_output_block(block_size, 0.);

    rt_cqt::ConstantQTransform<bins_per_octave, octave_number> cqt;
    cqt.init(hop_size); // Separate hop-sizes for each octave can be initialized using the .init(std::vector<int>
                        // octave_hop_sizes) overload
    cqt.init_sample_rate(sample_rate, block_size);

    cqt.input_block(audio_input_block.data(), block_size);
    const auto &schedule = cqt.get_cqt_schedule();
    for (const auto &element : schedule)
    {
        cqt.cqt(element);
        auto cqt_domain_buffer =
            cqt.get_octave_cqt_buffer(element.octave()); // The data could now be manipulated in cqt domain
        cqt.icqt(element);
    }
    auto cqt_audio_block = cqt.output_block(audio_input_block.size());
    for (std::size_t i = 0U; i < audio_input_block.size(); ++i)
    {
        audio_output_block[i] = cqt_audio_block[i];
    }

    return 0;
}
