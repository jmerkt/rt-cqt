# rt-cqt: Real-Time Constant-Q Transform
rt-cqt aims to be a reasonable fast header-only C++11 implementation of the Constant-Q transform. It is especially designed for easy usage in real-time audio applications, handling various block-sizes and samplerates. There are two different implementation currently: Constant-Q Transform (`include/ConstantQTransform.h`), doing a FFT per octave and Sliding Constant-Q Transform (`include/SlidingCqt.h`), which minimizes latency by updating the frequency bins for every sample. Polyphase IIR lowpasses are used to perform upsampling / downsampling in both cases. This reduces the computational costs per octave, because lower octaves don't require high sample-rates. 

## Constant-Q Transform
The implementation is roughly based on the [Judith C. Brown, Miller S. Puckette: An efficient algorithm  for the calculation  of a constant Q transform](http://academics.wellesley.edu/Physics/brown/pubs/effalgV92P2698-P2701.pdf) paper.
[pffft](https://github.com/marton78/pffft) is used to handle the ffts.

### Example Usage CQT
```cpp
#include "include/ConstantQTransform.h"

const int hopSize = 256;
const int octaveNumber = 9;
const int binsPerOctave = 12;

const int blockSize = 1024;
const double samplerate = 48000.;

std::vector<double> audioInputBlock(blockSize, 0.);
std::vector<double> audioOutputBlock(blockSize, 0.);

Cqt::ConstantQTransform<binsPerOctave, octaveNumber> cqt;
cqt.init(hopSize); // separate hop-sizes for each octave can be initialized using the .init(std::vector<int> octaveHopSizes) overload 
cqt.initFs(samplerate, blockSize);

cqt.inputBlock(audioInputBlock.data());
auto schedule = cqt.getCqtSchedule();
for(const auto& s : schedule)
{
  cqt.cqt(s);
  auto cqtDomainBuffer = cqt.getOctaveCqtBuffer(s.octave()); // the data could now be manipulated in cqt domain
  cqt.icqt(s);
}
auto cqtAudioBlock = cqt.outputBlock(audioInputBlock.size());
for(int i = 0; i < audioInputBlock.size(); i++)
{
  audioOutputBlock[i] = cqtAudioBlock[i];
}
```
Even though this library is header-only, the pffft implementation is not. Hence, the following files have to be linked into your project:
```cpp
submodules/pffft/pffft.c
submodules/pffft/pffft_common.c
submodules/pffft/pffft_double.c
```

## Sliding Constant-Q Transform
While the regular Constant-Q Transform above is doing a FFT per octave and hop, the CQT can also be implemented in a sliding fashion, which corresponds to a hop size of 1. This is still under development. The forward transform should work, but the backward transform is currently not working 100% as expected. 
The implementation roughly follows [Russell Bradford, John ffitch, Richard Dobson: SLIDING WITH A CONSTANT Q](https://purehost.bath.ac.uk/ws/portalfiles/portal/377255/constQ.pdf). Downsampling filterbanks are used to optimize processing time for lower octaves. 
Because there are no FFTs, the library is header-only without dependecies on pffft.


### Example Usage Sliding CQT
```cpp
#include "include/SlidingCqt.h"

constexpr int OctaveNumber = 9;
constexpr int BinsPerOctave = 12;
constexpr bool Windowing = false;

const int blockSize = 1024;
const double samplerate = 48000.;

std::vector<double> audioInputBlock(blockSize, 0.);
std::vector<double> audioOutputBlock(blockSize, 0.);
std::vector<std::complex<double>> cqtDomainBuffers[OctaveNumber][BinsPerOctave];

Cqt::SlidingCqt<BinsPerOctave, OctaveNumber, Windowing> cqt;
cqt.init(samplerate, blockSize);

for(unsigned i_octave = 0u; i_octave < OctaveNumber; i_octave++)
{
  // The sample rates and block sizes of each downsampled octave can be accessed
  const double octaveRate = cqt.getOctaveSampleRate(i_octave);
  const int octaveSize = cqt.getOctaveBlockSize(i_octave);
  for(unsigned i_tone = 0u; i_tone < BinsPerOctave; i_tone++)
  {
    cqtDomainBuffers[i_octave][i_tone].resize(octaveSize, {0., 0.});
  }
}

cqt.inputBlock(audioInputBlock.data(), blockSize);
for(unsigned i_octave = 0u; i_octave < OctaveNumber; i_octave++)
{
  // All complex values for each octave and bin are stored in circular buffers and can be accessed by getting a pointer to that buffer.
  // Because the data for each octave is downsampled, the number of samples per octave and block varies.
  const size_t nSamplesOctave = cqt.getSamplesToProcess(i_octave);
  CircularBuffer<std::complex<double>>* octaveCqtBuffer = cqt.getOctaveCqtBuffer(i_octave);
  for(unsigned i_tone = 0u; i_tone < BinsPerOctave; i_tone++)
  {
    octaveCqtBuffer[i_tone].pullBlock(cqtDomainBuffers[i_octave][i_tone].data(), nSamplesOctave);
    for(size_t i_sample = 0u; i_sample < nSamplesOctave; i_sample++)
    {
      // Here we can manipulate the complex values for each sample, bin and octave
      cqtDomainBuffers[i_octave][i_tone][i_sample] *= 2.0;
    }
    octaveCqtBuffer[i_tone].pushBlock(cqtDomainBuffers[i_octave][i_tone].data(), nSamplesOctave);
  }
}
auto cqtAudioBlock = cqt.outputBlock(audioInputBlock.size());
for(int i_sample = 0; i_sample < audioInputBlock.size(); i_sample++)
{
  audioOutputBlock[i_sample] = cqtAudioBlock[i_sample];
}
```

## Python Bindings
This is WIP. Python bindings can be found in the `python-bindings` folder, using [pybind11](https://github.com/pybind/pybind11).

## Example Projects
[CQT Analyzer Audio-Plugin based on iPlug2](https://github.com/jmerkt/cqt-analyzer)

[WIP Reverb plugin operating in Sliding CQT domain](https://github.com/jmerkt/harmonic-reverb)

## Limitations and Future Work
* While the Polyphase IIR Lowpasses are cheap, they distort the phase of the signal. For future a linear phase approach could be added.


