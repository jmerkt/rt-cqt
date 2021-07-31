# rt-cqt: Real-Time Constant-Q-Transform
rt-cqt is a reasonable fast header-only C++11 implementation of the Constant-Q-Transform. It is especially designed for the usage in real-time audio applications handling various block-sizes and samplerates.
The implementation is roughly based on the [Judith C. Brown, Miller S. Puckette: An efficient algorithm  for the calculation  of a constant Q transform](http://academics.wellesley.edu/Physics/brown/pubs/effalgV92P2698-P2701.pdf) paper.
As fft library pffft is used.

## Example Usage
```cpp
const int hopSize = 256;
const int octaveNumber = 9;
const int binsPerOctave = 12;

const int blockSize = 1024;
const double samplerate = 48000.;

std::vector<double> audioBlock(blockSize, 0.);

Cqt::ConstantQTransform<binsPerOctave, octaveNumber> cqt;
cqt.init(hopSize);
cqt.initFs(samplerate, blockSize);

cqt.inputBlock(audioBlock.data());
auto schedule = cqt.getCqtSchedule();
for(const auto& s : schedule)
{
  cqt.cqt(s);
  auto cqtBuffer = cqt.getOctaveCqtBuffer(s.octave);
  // manipulate data here
  cqt.icqt(s);
}
```
