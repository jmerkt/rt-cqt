/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

#include <pybind11/numpy.h>

#include "../../include/ConstantQTransform.h"

namespace Cqt
{
    template <int B, int OctaveNumber>
    class Python_ConstantQTransform : public ConstantQTransform<B, OctaveNumber>
    {
    public:
        void Python_inputBlock(std::vector<double> &data)
        {
            this->inputBlock(data.data(), data.size());
        };

        std::vector<double> Python_outputBlock(const int blockSize)
        {
            std::vector<double> outputVector(blockSize, 0.);
            const auto outputBlock = this->outputBlock(blockSize);
            std::memcpy(outputVector.data(), outputBlock, blockSize * sizeof(double));
            return outputVector;
        };
    };
}