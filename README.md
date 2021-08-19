# rt-cqt: Real-Time Constant-Q Transform
rt-cqt aims to be a reasonable fast header-only C++11 implementation of the Constant-Q transform. It is especially designed for easy usage in real-time audio applications, handling various block-sizes and samplerates.
The implementation is roughly based on the [Judith C. Brown, Miller S. Puckette: An efficient algorithm  for the calculation  of a constant Q transform](http://academics.wellesley.edu/Physics/brown/pubs/effalgV92P2698-P2701.pdf) paper.
[pffft](https://github.com/marton78/pffft) is used to handle the ffts, Polyphase IIR lowpasses to perform upsampling / downsampling.

## Example Usage
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
  auto cqtDomainBuffer = cqt.getOctaveCqtBuffer(s.octave); // the data could now be manipulated in cqt domain
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

## Example Project
[CQT Analyzer Audio-Plugin based on iPlug2](https://github.com/jmerkt/cqt-analyzer)

## Limitations and Future Work
* While the Polyphase IIR Lowpasses are cheap, they distort the phase of the signal. For future a linear phase approach could be added.
