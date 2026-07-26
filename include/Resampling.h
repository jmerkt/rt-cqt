/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

#include "../submodules/audio-utils/include/Utils.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Cqt
{

    template <typename FloatType>
    class Delay
    {
    public:
        Delay() = default;
        ~Delay() = default;
        inline void reset() { mStorage = static_cast<FloatType>(0.); }

        inline void process(FloatType *const data, const int blockSize)
        {
            for (int i = 0; i < blockSize; i++)
            {
                const FloatType input = data[i];
                data[i] = mStorage;
                mStorage = input;
            }
        }

    private:
        FloatType mStorage{0.};
    };

    template <typename FloatType>
    class FirstOrderAllpass
    {
    public:
        FirstOrderAllpass() = default;
        ~FirstOrderAllpass() = default;

        inline void initCoeff(FloatType ak)
        {
            mAk = static_cast<FloatType>(ak);
            mXm1 = static_cast<FloatType>(0.);
            mYm1 = static_cast<FloatType>(0.);
        };

        inline void process(FloatType *const samples, const int blocksize)
        {
            for (int i = 0; i < blocksize; i++)
            {
                const FloatType sample = samples[i];
                samples[i] = mAk * (sample - mYm1) + mXm1;
                mXm1 = sample;
                mYm1 = samples[i];
            }
        };

    private:
        FloatType mAk{0.};
        FloatType mXm1{0.};
        FloatType mYm1{0.};
    };

    /*
     * Polyphase IIR lowpass for resampling blocks by a factor of two.
     *
     * See "Digital signal processing schemes for efficient interpolation and
     * decimation" by Reinaldo A. Valenzuela and A. G. Constantinides. This
     * class handles allocation and downsampling or upsampling for incoming
     * blocks. A single instance must only process in one direction.
     */
    template <typename FloatType, size_t AllpassNumber>
    class HalfBandLowpass
    {
    public:
        HalfBandLowpass();
        ~HalfBandLowpass() = default;

        /*
        Initialize the filters and allocate memory. Has to get called before processing starts.
        */
        bool init(const int expectedBlockSize = 128, bool isDownsampling = true, double transitionBandwidth = 0.02);

        FloatType *processDown(const FloatType *const inputBlock);
        FloatType *processUp(const FloatType *const inputBlock);

        int getOutputBlockSize() { return mTargetBlockSize; };
        int getInputBlockSize() { return mInputBlockSize; };

    private:
        double mTransitionBandwidth;
        size_t mAllpassNumberTotal;
        size_t mFilterOrder;
        std::vector<double> mCoefficients;

        FirstOrderAllpass<FloatType> mDirectPathFilters[AllpassNumber];
        FirstOrderAllpass<FloatType> mDelayPathFilters[AllpassNumber];
        Delay<FloatType> mDelay;

        int mTargetBlockSize{0};
        int mInputBlockSize{1};
        int mFilterBufferSize{0};

        std::vector<FloatType> mDirectPathBuffer;
        std::vector<FloatType> mDelayPathBuffer;
        std::vector<FloatType> mOutputBlock;

        std::vector<double> filterDesign();
    };

    template <typename FloatType, size_t AllpassNumber>
    HalfBandLowpass<FloatType, AllpassNumber>::HalfBandLowpass()
    {
        mAllpassNumberTotal = AllpassNumber * 2;
        mFilterOrder = 2 * mAllpassNumberTotal + 1;
    }

    template <typename FloatType, size_t AllpassNumber>
    inline bool HalfBandLowpass<FloatType, AllpassNumber>::init(const int expectedBlockSize,
                                                                bool isDownsampling,
                                                                double transitionBandwidth)
    {
        if (expectedBlockSize <= 0 || (isDownsampling && ((expectedBlockSize % 2) != 0)))
        {
            throw std::invalid_argument(
                "Resampling blocks must be positive; downsampling blocks must also be divisible by two");
        }

        // init filters
        mTransitionBandwidth = transitionBandwidth * 2. * audio_utils::Pi<double>();
        mCoefficients.clear();
        mCoefficients = filterDesign();
        int filterCount = 0;
        for (size_t i = 0; i < mCoefficients.size(); i += 2)
        {
            mDirectPathFilters[filterCount].initCoeff(mCoefficients[i]);
            filterCount++;
        }
        filterCount = 0;
        for (size_t i = 1; i < mCoefficients.size(); i += 2)
        {
            mDelayPathFilters[filterCount].initCoeff(mCoefficients[i]);
            filterCount++;
        }
        // init buffers
        mDirectPathBuffer.clear();
        mDelayPathBuffer.clear();
        mOutputBlock.clear();
        mDelay.reset();
        mInputBlockSize = expectedBlockSize;
        if (isDownsampling)
        {
            mFilterBufferSize = expectedBlockSize / 2;
            mTargetBlockSize = expectedBlockSize / 2;
        }
        else
        {
            mFilterBufferSize = expectedBlockSize;
            mTargetBlockSize = expectedBlockSize * 2;
        }
        mDirectPathBuffer.resize(mFilterBufferSize, static_cast<FloatType>(0.));
        mDelayPathBuffer.resize(mFilterBufferSize, static_cast<FloatType>(0.));
        mOutputBlock.resize(mTargetBlockSize, static_cast<FloatType>(0.));
        return true;
    }

    template <typename FloatType, size_t AllpassNumber>
    inline FloatType *HalfBandLowpass<FloatType, AllpassNumber>::processDown(const FloatType *const inputBlock)
    {
        int outCountDirect = 0;
        for (int i = 0; i < mInputBlockSize; i += 2)
        {
            mDirectPathBuffer[outCountDirect] = inputBlock[i];
            outCountDirect++;
        }
        int outCountDelay = 0;
        for (int i = 1; i < mInputBlockSize; i += 2)
        {
            mDelayPathBuffer[outCountDelay] = inputBlock[i];
            outCountDelay++;
        }
        mDelay.process(mDelayPathBuffer.data(), mFilterBufferSize);
        for (size_t i = 0; i < AllpassNumber; i++)
        {
            mDirectPathFilters[i].process(mDirectPathBuffer.data(), mFilterBufferSize);
            mDelayPathFilters[i].process(mDelayPathBuffer.data(), mFilterBufferSize);
        }
        for (int i = 0; i < mTargetBlockSize; i++)
        {
            mOutputBlock[i] = static_cast<FloatType>(0.5) * (mDirectPathBuffer[i] + mDelayPathBuffer[i]);
        }
        return mOutputBlock.data();
    };

    template <typename FloatType, size_t AllpassNumber>
    inline FloatType *HalfBandLowpass<FloatType, AllpassNumber>::processUp(const FloatType *const inputBlock)
    {
        for (int i = 0; i < mInputBlockSize; i++)
        {
            mDirectPathBuffer[i] = inputBlock[i];
            mDelayPathBuffer[i] = inputBlock[i];
        }
        for (size_t i = 0; i < AllpassNumber; i++)
        {
            mDirectPathFilters[i].process(mDirectPathBuffer.data(), mFilterBufferSize);
            mDelayPathFilters[i].process(mDelayPathBuffer.data(), mFilterBufferSize);
        }
        int inCountDirect = 0;
        for (int i = 0; i < mTargetBlockSize; i += 2)
        {
            mOutputBlock[i] = mDirectPathBuffer[inCountDirect];
            inCountDirect += 1;
        }
        int inCountDelay = 0;
        for (int i = 1; i < mTargetBlockSize; i += 2)
        {
            mOutputBlock[i] = mDelayPathBuffer[inCountDelay];
            inCountDelay += 1;
        }
        return mOutputBlock.data();
    };

    template <typename FloatType, size_t AllpassNumber>
    inline std::vector<double> HalfBandLowpass<FloatType, AllpassNumber>::filterDesign()
    {
        // step 1
        const double k = std::pow(std::tan((audio_utils::Pi<double>() - mTransitionBandwidth) / 4.), 2);
        const double k_dash = std::sqrt(1. - std::pow(k, 2));
        const double e = (1. / 2.) * ((1. - std::sqrt(k_dash)) / (1. + std::sqrt(k_dash)));
        const double q = e + 2. * std::pow(e, 5) + 15. * std::pow(e, 9.) + 150. * std::pow(e, 13.);
        // step 2
        const size_t n = mFilterOrder;
        // step 3
        std::vector<double> w;
        std::vector<double> a_dash;
        for (size_t i = 1; i <= ((n - 1) / 2); i++)
        {
            // w_i
            double delta = 1.;
            double num = 0.;
            double m = 0.;
            while (delta > 1.e-100)
            {
                delta = std::pow((-1.), m) * std::pow(q, (m * (m + 1.))) *
                        std::sin((2. * m + 1.) * audio_utils::Pi<double>() * static_cast<double>(i) /
                                 static_cast<double>(n));
                num += delta;
                m += 1.;
            }
            num = 2. * std::pow(q, (1. / 4.)) * num;
            delta = 1.;
            double den = 0.;
            m = 1.;
            while (delta > 1.e-100)
            {
                delta = std::pow((-1.), m) * std::pow(q, std::pow(m, 2.)) *
                        std::cos(2. * m * audio_utils::Pi<double>() * static_cast<double>(i) / static_cast<double>(n));
                den += delta;
                m += 1.;
            }
            den = den * 2. + 1.;
            double w_i = num / den;
            w.push_back(w_i);
            // a'_i
            num = std::pow(((1. - std::pow(w_i, 2.) * k) * (1. - std::pow(w_i, 2.) / k)), (1. / 2.));
            den = 1. + std::pow(w_i, 2.);
            double a_dash_i = num / den;
            a_dash.push_back(a_dash_i);
        }
        // step 4
        std::vector<double> a;
        for (double a_dash_i : a_dash)
        {
            double a_i = (1. - a_dash_i) / (1. + a_dash_i);
            a.push_back(a_i);
        }
        return a;
    };

    enum class DirectionConfig
    {
        Up = 0,
        Down,
        UpDown,
        DownUp
    };

    /**
    Resampling of audio blocks by a given power of 2.
    */
    template <typename FloatType, size_t AllpassNumber>
    class ResamplingHandler
    {
    public:
        ResamplingHandler(double transitionBandwidth = 0.02);
        ~ResamplingHandler() = default;
        /**
        Initialization of the ResamplingHandler
        The configuration gives the amount of memory that will be allocated.
        DirectionConfig has to match the order the processing functions will be called from outside.
        */
        void init(const int powToExponent = 0,
                  const int expectedBlockSize = 128,
                  DirectionConfig directionConfig = DirectionConfig::UpDown);

        FloatType *processDown(const FloatType *const inputBlock);
        FloatType *processUp(const FloatType *const inputBlock);

        int getOutputBlockSizeUp() { return mTargetBlockSizeUp; };
        int getOutputBlockSizeDown() { return mTargetBlockSizeDown; };

    private:
        double mTransitionBandwidth;
        int mPowTwoFactor{0};
        int mTargetBlockSizeUp{0};
        int mTargetBlockSizeDown{0};
        std::vector<HalfBandLowpass<FloatType, AllpassNumber>> mDownFilters;
        std::vector<HalfBandLowpass<FloatType, AllpassNumber>> mUpFilters;
    };

    template <typename FloatType, size_t AllpassNumber>
    ResamplingHandler<FloatType, AllpassNumber>::ResamplingHandler(double transitionBandwidth)
    {
        mTransitionBandwidth = transitionBandwidth;
    }

    template <typename FloatType, size_t AllpassNumber>
    inline void ResamplingHandler<FloatType, AllpassNumber>::init(const int powToExponent,
                                                                  const int expectedBlockSize,
                                                                  DirectionConfig directionConfig)
    {
        if (powToExponent < 0 || powToExponent >= 30)
        {
            throw std::invalid_argument("The power-of-two exponent must be between 0 and 29");
        }
        if (expectedBlockSize <= 0)
        {
            throw std::invalid_argument("The resampling block size must be positive");
        }

        mDownFilters.clear();
        mUpFilters.clear();
        mDownFilters.reserve(static_cast<std::size_t>(powToExponent));
        mUpFilters.reserve(static_cast<std::size_t>(powToExponent));
        mPowTwoFactor = powToExponent;
        mTargetBlockSizeDown = expectedBlockSize;
        mTargetBlockSizeUp = expectedBlockSize;

        const int resamplingFactor = 1 << mPowTwoFactor;
        if ((directionConfig == DirectionConfig::Down) || (directionConfig == DirectionConfig::DownUp))
        {
            if ((expectedBlockSize % resamplingFactor) != 0)
            {
                throw std::invalid_argument("The downsampling block size must be divisible by the resampling factor");
            }
        }
        if ((directionConfig == DirectionConfig::Up || directionConfig == DirectionConfig::UpDown) &&
            expectedBlockSize > (std::numeric_limits<int>::max() / resamplingFactor))
        {
            throw std::overflow_error("The upsampling block size is too large");
        }

        for (int i = 0; i < mPowTwoFactor; i++)
        {
            const int stageFactor = 1 << i;
            if (directionConfig == DirectionConfig::Down)
            {
                mDownFilters.emplace_back();
                mDownFilters.back().init(expectedBlockSize / stageFactor, true, mTransitionBandwidth);
            }
            else if (directionConfig == DirectionConfig::Up)
            {
                mUpFilters.emplace_back();
                mUpFilters.back().init(expectedBlockSize * stageFactor, false, mTransitionBandwidth);
            }
            else if (directionConfig == DirectionConfig::DownUp)
            {
                mDownFilters.emplace_back();
                mUpFilters.emplace_back();
                mDownFilters.back().init(expectedBlockSize / stageFactor, true, mTransitionBandwidth);
                mUpFilters.back().init(expectedBlockSize / (1 << (mPowTwoFactor - i)), false, mTransitionBandwidth);
            }
            else if (directionConfig == DirectionConfig::UpDown)
            {
                mUpFilters.emplace_back();
                mDownFilters.emplace_back();
                mUpFilters.back().init(expectedBlockSize * stageFactor, false, mTransitionBandwidth);
                mDownFilters.back().init(expectedBlockSize * (1 << (mPowTwoFactor - i)), true, mTransitionBandwidth);
            }
        }

        if (!mDownFilters.empty())
        {
            mTargetBlockSizeDown = mDownFilters.back().getOutputBlockSize();
        }
        if (!mUpFilters.empty())
        {
            mTargetBlockSizeUp = mUpFilters.back().getOutputBlockSize();
        }
    };

    template <typename FloatType, size_t AllpassNumber>
    inline FloatType *ResamplingHandler<FloatType, AllpassNumber>::processDown(const FloatType *const inputBlock)
    {
        assert(mPowTwoFactor == 0 || static_cast<int>(mDownFilters.size()) == mPowTwoFactor);
        const FloatType *inBlock = inputBlock;
        FloatType *outBlock = const_cast<FloatType *>(inputBlock);
        for (int i = 0; i < mPowTwoFactor; i++)
        {
            outBlock = mDownFilters[i].processDown(inBlock);
            inBlock = outBlock;
        }
        return outBlock;
    };

    template <typename FloatType, size_t AllpassNumber>
    inline FloatType *ResamplingHandler<FloatType, AllpassNumber>::processUp(const FloatType *const inputBlock)
    {
        assert(mPowTwoFactor == 0 || static_cast<int>(mUpFilters.size()) == mPowTwoFactor);
        const FloatType *inBlock = inputBlock;
        FloatType *outBlock = const_cast<FloatType *>(inputBlock);
        for (int i = 0; i < mPowTwoFactor; i++)
        {
            outBlock = mUpFilters[i].processUp(inBlock);
            inBlock = outBlock;
        }
        return outBlock;
    };

}
