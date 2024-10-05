#include "SlidingCqt.h"

int main(int argc, char *argv[])
{
    constexpr int octave_number = 9;
    constexpr int bins_per_octave = 12;
    constexpr bool use_windowing = false;

    const int block_size = 1024;
    const double sample_rate = 48000.;

    std::vector<double> audio_input_block(block_size, 0.);
    std::vector<double> audio_output_block(block_size, 0.);
    std::vector<std::complex<double>> cqt_domain_buffers[octave_number][bins_per_octave];

    Cqt::SlidingCqt<bins_per_octave, octave_number, use_windowing> cqt;
    cqt.init(sample_rate, block_size);

    for (int octave = 0; octave < octave_number; ++octave)
    {
        // The sample rates and block sizes of each downsampled octave can be accessed
        const double octave_rate = cqt.getOctaveSampleRate(octave);
        const int octave_size = cqt.getOctaveBlockSize(octave);
        for (int tone = 0; tone < bins_per_octave; ++tone)
        {
            cqt_domain_buffers[octave][tone].resize(octave_size, {0., 0.});
        }
    }

    cqt.inputBlock(audio_input_block.data(), block_size);
    for (int octave = 0; octave < octave_number; ++octave)
    {
        // Because the data for each octave is downsampled, the number of samples per octave and block varies.
        const size_t number_octave_samples = cqt.getSamplesToProcess(octave);
        for (int tone = 0; tone < bins_per_octave; ++tone)
        {
            // Pull the cqt domain data for this octave and tone
            cqt.pullBinCqtData(octave, tone, cqt_domain_buffers[octave][tone].data());
            for (size_t sample = 0U; sample < number_octave_samples; ++sample)
            {
                // Here, we can manipulate the complex values for each sample, bin and octave
                cqt_domain_buffers[octave][tone][sample] *= 2.0;
            }
            // Push back the manipulated cqt domain data for this octave and tone
            // The number of pulled and pushed samples must match number_octave_samples
            cqt.pushBinCqtData(octave, tone, cqt_domain_buffers[octave][tone].data());
        }
    }
    auto cqt_audio_block = cqt.outputBlock(audio_input_block.size());
    for (size_t sample = 0U; sample < audio_input_block.size(); ++sample)
    {
        audio_output_block[sample] = cqt_audio_block[sample];
    }

    return 0;
}
