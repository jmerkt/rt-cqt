# rt-cqt: Real-Time Constant-Q-Transform
rt-cqt is a reasonable fast header-only C++11 implementation of the Constant-Q-Transform. It is especially designed for the usage in real-time audio applications handling various block-sizes and samplerates.
The implementation is roughly based on the [Judith C. Brown, Miller S. Puckette: An efficient algorithm  for the calculation  of a constant Q transform](http://academics.wellesley.edu/Physics/brown/pubs/effalgV92P2698-P2701.pdf) paper.
As fft library pffft is used.

## Example Usage
```cpp
#include <include/ConstantQTransform>

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
auto cqtAudioBlock = cqt.outputBlock(audioBlock.size());
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

## Limitations and Future Work

