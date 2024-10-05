#include "ConstantQTransform.h"

int main(int argc, char *argv[])
{
    const int hop_size = 256;
    const int octave_number = 9;
    const int bins_per_octave = 12;

    const int block_size = 1024;
    const double sample_rate = 48000.;

    std::vector<double> audio_input_block(block_size, 0.);
    std::vector<double> audio_output_block(block_size, 0.);

    Cqt::ConstantQTransform<bins_per_octave, octave_number> cqt;
    cqt.init(hop_size); // Separate hop-sizes for each octave can be initialized using the .init(std::vector<int> octaveHopSizes) overload
    cqt.initFs(sample_rate, block_size);

    cqt.inputBlock(audio_input_block.data(), block_size);
    const auto &schedule = cqt.getCqtSchedule();
    for (const auto &s : schedule)
    {
        cqt.cqt(s);
        auto cqt_domain_buffer = cqt.getOctaveCqtBuffer(s.octave()); // The data could now be manipulated in cqt domain
        cqt.icqt(s);
    }
    auto cqt_audio_block = cqt.outputBlock(audio_input_block.size());
    for (size_t i = 0U; i < audio_input_block.size(); ++i)
    {
        audio_output_block[i] = cqt_audio_block[i];
    }

    return 0;
}