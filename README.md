# rt-cqt: Real-Time Constant-Q Transform
rt-cqt aims to be a reasonable fast, header-only C++11 library for performing the Constant-Q Transform (CQT), optimized for real-time audio applications. It supports dynamic handling of different block sizes and sample rates, making it ideal for scenarios with varying audio processing demands. The library offers two distinct implementations:

1. Constant-Q Transform `(include/ConstantQTransform.h)`: This version performs a Fast Fourier Transform (FFT) for each octave.
2. Sliding Constant-Q Transform `(include/SlidingCqt.h)`: This version minimizes latency by continuously updating frequency bins with every new audio sample, making it particularly suitable low-latency applications.

Both implementations utilize polyphase IIR lowpass filters for efficient upsampling and downsampling, which reduces computational overhead by processing lower octaves at reduced sample rates. 

## Constant-Q Transform
The implementation is roughly based on the [Judith C. Brown, Miller S. Puckette: An efficient algorithm  for the calculation  of a constant Q transform](http://academics.wellesley.edu/Physics/brown/pubs/effalgV92P2698-P2701.pdf) paper.
[pffft](https://github.com/marton78/pffft) is used to handle the ffts.

An example can be found in `examples/cqt.cpp`.

Even though this library is header-only, the pffft implementation is not. Hence, the following files have to be linked into your project:
```cpp
submodules/pffft/pffft.c
submodules/pffft/pffft_common.c
submodules/pffft/pffft_double.c
```

## Sliding Constant-Q Transform
While the regular Constant-Q Transform above is doing a FFT per octave and hop, the CQT can also be implemented in a sliding fashion, which corresponds to a hop size of 1.
The implementation roughly follows [Russell Bradford, John ffitch, Richard Dobson: SLIDING WITH A CONSTANT Q](https://purehost.bath.ac.uk/ws/portalfiles/portal/377255/constQ.pdf). Because there are no FFTs, the library is header-only without dependecies on pffft.

An example can be found in `examples/scqt.cpp`.

## Python Bindings
This is WIP. Python bindings can be found in the `python-bindings` folder, using [pybind11](https://github.com/pybind/pybind11).

## Example Projects
[CQT Analyzer Audio-Plugin based on iPlug2](https://github.com/jmerkt/cqt-analyzer)

[WIP Reverb plugin operating in Sliding CQT domain](https://github.com/jmerkt/harmonic-reverb)

## Limitations and Future Work
* While the Polyphase IIR Lowpasses are cheap, they distort the phase of the signal. For future, a linear phase approach could be added.


