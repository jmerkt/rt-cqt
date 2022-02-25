/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once


#include <pybind11/numpy.h>

#include "../../include/SlidingCqt.h"

namespace Cqt
{
    template <int B, int OctaveNumber>
    class Python_SlidingCqt : public SlidingCqt<B, OctaveNumber>
    {
    public:

        void Python_inputBlock(std::vector<double>& data, const int blockSize) 
        { 
            this->inputBlock(data.data(), blockSize); 
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