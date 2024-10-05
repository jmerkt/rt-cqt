#include "ConstantQTransform.h"

int main(int argc, char *argv[])
{
    const int hop_size = 256;
    const int octave_number = 9;
    const int bins_per_octave = 12;

    const int block_size = 1024;
    const double sample_rate = 48000.;

    std::vector<double> audioInputBlock(block_size, 0.);
    std::vector<double> audioOutputBlock(block_size, 0.);

    Cqt::ConstantQTransform<bins_per_octave, octave_number> cqt;
    cqt.init(hop_size); // Separate hop-sizes for each octave can be initialized using the .init(std::vector<int> octaveHopSizes) overload
    cqt.initFs(sample_rate, block_size);

    cqt.inputBlock(audioInputBlock.data(), block_size);
    auto schedule = cqt.getCqtSchedule();
    for (const auto &s : schedule)
    {
        cqt.cqt(s);
        auto cqtDomainBuffer = cqt.getOctaveCqtBuffer(s.octave()); // The data could now be manipulated in cqt domain
        cqt.icqt(s);
    }
    auto cqtAudioBlock = cqt.outputBlock(audioInputBlock.size());
    for (int i = 0; i < audioInputBlock.size(); i++)
    {
        audioOutputBlock[i] = cqtAudioBlock[i];
    }

    return 0;
}