#include "SlidingCqt.h"

int main(int argc, char *argv[])
{
    constexpr int octave_number = 9;
    constexpr int bins_per_octave = 12;
    constexpr bool use_windowing = false;

    const int blockSize = 1024;
    const double samplerate = 48000.;

    std::vector<double> audioInputBlock(blockSize, 0.);
    std::vector<double> audioOutputBlock(blockSize, 0.);
    std::vector<std::complex<double>> cqtDomainBuffers[octave_number][bins_per_octave];

    Cqt::SlidingCqt<bins_per_octave, octave_number, use_windowing> cqt;
    cqt.init(samplerate, blockSize);

    for (unsigned i_octave = 0u; i_octave < octave_number; i_octave++)
    {
        // The sample rates and block sizes of each downsampled octave can be accessed
        const double octaveRate = cqt.getOctaveSampleRate(i_octave);
        const int octaveSize = cqt.getOctaveBlockSize(i_octave);
        for (unsigned i_tone = 0u; i_tone < bins_per_octave; i_tone++)
        {
            cqtDomainBuffers[i_octave][i_tone].resize(octaveSize, {0., 0.});
        }
    }

    cqt.inputBlock(audioInputBlock.data(), blockSize);
    for (unsigned i_octave = 0u; i_octave < octave_number; i_octave++)
    {
        // All complex values for each octave and bin are stored in circular buffers and can be accessed by getting a pointer to that buffer.
        // Because the data for each octave is downsampled, the number of samples per octave and block varies.
        const size_t nSamplesOctave = cqt.getSamplesToProcess(i_octave);
        audio_utils::CircularBuffer<std::complex<double>> *octaveCqtBuffer = cqt.getOctaveCqtBuffer(i_octave);
        for (unsigned i_tone = 0u; i_tone < bins_per_octave; i_tone++)
        {
            octaveCqtBuffer[i_tone].pullBlock(cqtDomainBuffers[i_octave][i_tone].data(), nSamplesOctave);
            for (size_t i_sample = 0u; i_sample < nSamplesOctave; i_sample++)
            {
                // Here we can manipulate the complex values for each sample, bin and octave
                cqtDomainBuffers[i_octave][i_tone][i_sample] *= 2.0;
            }
            octaveCqtBuffer[i_tone].pushBlock(cqtDomainBuffers[i_octave][i_tone].data(), nSamplesOctave);
        }
    }
    auto cqtAudioBlock = cqt.outputBlock(audioInputBlock.size());
    for (int i_sample = 0; i_sample < audioInputBlock.size(); i_sample++)
    {
        audioOutputBlock[i_sample] = cqtAudioBlock[i_sample];
    }

    return 0;
}
