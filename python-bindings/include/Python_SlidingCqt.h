/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

#include "../../include/SlidingCqt.h"

#include <complex>
#include <cstring>
#include <vector>

namespace Cqt
{
    template <int B, int OctaveNumber, bool Windowing>
    class Python_SlidingCqt
    {
    public:
        void init(const double samplerate, const int blockSize)
        {
            mTransform.init(samplerate, blockSize);
        }

        void inputBlock(std::vector<double> &data, const int blockSize)
        {
            mTransform.inputBlock(data.data(), blockSize);
        }

        std::vector<double> outputBlock(const int blockSize)
        {
            std::vector<double> outputVector(blockSize, 0.);
            const double *const outputBlock = mTransform.outputBlock(blockSize);
            std::memcpy(outputVector.data(), outputBlock, blockSize * sizeof(double));
            return outputVector;
        }

        std::vector<std::complex<double>> getOctaveValues(const int octave)
        {
            std::vector<std::complex<double>> valueVector(B, {0., 0.});
            audio_utils::CircularBuffer<std::complex<double>> *octaveCqtBuffer =
                mTransform.getOctaveCqtBuffer(octave);
            for (int i_tone = 0; i_tone < B; i_tone++)
            {
                valueVector[i_tone] = octaveCqtBuffer[i_tone].pullDelaySample(0);
            }
            return valueVector;
        }

        std::vector<double> getOctaveBinFreqs(const int octave)
        {
            std::vector<double> valueVector(B, 0.);
            const double *octaveBinFreqs = mTransform.getOctaveBinFreqs(octave);
            for (int i_tone = 0; i_tone < B; i_tone++)
            {
                valueVector[i_tone] = octaveBinFreqs[i_tone];
            }
            return valueVector;
        }

    private:
        SlidingCqt<B, OctaveNumber, Windowing> mTransform;
    };
}
