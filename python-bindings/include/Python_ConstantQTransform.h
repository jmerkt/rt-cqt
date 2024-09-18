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
        void Python_inputBlock(std::vector<double> &data, const int blockSize)
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

        pybind11::array_t<ScheduleElement> Python_getCqtSchedule()
        {
            auto result = pybind11::array_t<ScheduleElement>(this->mCqtSchedule.size());
            pybind11::buffer_info buf = result.request();
            ScheduleElement *ptr = static_cast<ScheduleElement *>(buf.ptr);
            for (int i = 0; i < this->mCqtSchedule.size(); i++)
            {
                ptr[i] = this->mCqtSchedule.at(i);
            }
            return result;
        };

        void Python_cqt(pybind11::array_t<ScheduleElement> schedule)
        {
            pybind11::buffer_info buf = schedule.request();
            ScheduleElement *ptr = static_cast<ScheduleElement *>(buf.ptr);
            this->cqt(ptr[0]);
        };

        void Python_icqt(pybind11::array_t<ScheduleElement> schedule)
        {
            pybind11::buffer_info buf = schedule.request();
            ScheduleElement *ptr = static_cast<ScheduleElement *>(buf.ptr);
            this->icqt(ptr[0]);
        };
    };
}