/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

#include "../../include/ConstantQTransform.h"

#include <cstring>
#include <vector>

namespace Cqt
{
    template <int B, int OctaveNumber>
    class Python_ConstantQTransform
    {
    public:
        void init(const int hopSize) { mTransform.init(hopSize); }

        void initFs(const double samplerate, const int blockSize) { mTransform.initFs(samplerate, blockSize); }

        void inputBlock(std::vector<double> &data)
        {
            mTransform.inputBlock(data.data(), static_cast<int>(data.size()));
        }

        std::vector<double> outputBlock(const int blockSize)
        {
            std::vector<double> outputVector(blockSize, 0.);
            const double *const outputBlock = mTransform.outputBlock(blockSize);
            std::memcpy(outputVector.data(), outputBlock, blockSize * sizeof(double));
            return outputVector;
        }

        std::vector<ScheduleElement> &getCqtSchedule() { return mTransform.getCqtSchedule(); }

        void cqt(const ScheduleElement schedule) { mTransform.cqt(schedule); }

        void icqt(const ScheduleElement schedule) { mTransform.icqt(schedule); }

        CqtBufferType *getOctaveCqtBuffer(const int octave) { return mTransform.getOctaveCqtBuffer(octave); }

    private:
        ConstantQTransform<B, OctaveNumber> mTransform;
    };
}
