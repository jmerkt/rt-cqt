# rt-cqt: Real-Time Constant-Q Transform
rt-cqt aims to be a reasonable fast header-only C++11 implementation of the Constant-Q transform. It is especially designed for easy usage in real-time audio applications, handling various block-sizes and samplerates. There are two different implementation currently: Constant-Q Transform (`include/ConstantQTransform.h`), doing a FFT per octave and Sliding Constant-Q Transform (`include/SlidingCqt.h`), which minimizes latency by updating the frequency bins for every sample. Polyphase IIR lowpasses are used to perform upsampling / downsampling in both cases. This reduces the computational costs per octave, because lower octaves don't require high sample-rates. 

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
While the regular Constant-Q Transform above is doing a FFT per octave and hop, the CQT can also be implemented in a sliding fashion, which corresponds to a hop size of 1. This is still under development. The forward transform should work, but the backward transform is currently not working 100% as expected. 
The implementation roughly follows [Russell Bradford, John ffitch, Richard Dobson: SLIDING WITH A CONSTANT Q](https://purehost.bath.ac.uk/ws/portalfiles/portal/377255/constQ.pdf). Downsampling filterbanks are used to optimize processing time for lower octaves. 
Because there are no FFTs, the library is header-only without dependecies on pffft.

An example can be found in `examples/scqt.cpp`.

## Python Bindings
This is WIP. Python bindings can be found in the `python-bindings` folder, using [pybind11](https://github.com/pybind/pybind11).

## Example Projects
[CQT Analyzer Audio-Plugin based on iPlug2](https://github.com/jmerkt/cqt-analyzer)

[WIP Reverb plugin operating in Sliding CQT domain](https://github.com/jmerkt/harmonic-reverb)

## Limitations and Future Work
* While the Polyphase IIR Lowpasses are cheap, they distort the phase of the signal. For future a linear phase approach could be added.


