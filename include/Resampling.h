/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

#include <cmath>
#include <vector>
#include "Utils.h"

namespace Cqt
{ 

template<typename FloatType>
class Delay
{
public:
	Delay() = default;
	~Delay() = default;

	inline void processBlock(FloatType* const data, const int blockSize)
	{
		for (int i = 0; i < blockSize; i++)
		{
			const FloatType input = data[i];
			data[i] = mStorage;
			mStorage = input;
		}
	}
private:
	FloatType mStorage{ 0. };
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
	};

	inline FloatType processSample(const FloatType sample)
	{
		const FloatType y = mAk * (sample - mYm1) + mXm1;
		mXm1 = sample;
		mYm1 = y;
		return y;
	};

	inline void processBlock(FloatType* const samples, const int blocksize)
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
	FloatType mAk{ 0. };
	FloatType mXm1{ 0. };
	FloatType mYm1{ 0. };
};


/*
*	Polyphase IIR Lowpass to sample up/down by 2.
	See paper: "Reinaldo A. Valenzuela, A. G. Constantinides: Digital signal processing schemes for efficient interpolation and decimation".
	This class handles memory allocation and downsampling/upsampling filtering for incoming samples/blocks by 2.
	One instance of this class can either be used for Downsampling OR for Upsampling. Not for both at the same time.
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
	bool init(const int expectedBlockSize = 128, bool isDownsampling = true, bool isSampleBased = false, double transitionBandwidth = 0.02);

	/*
	Call with 2 * fs_target.
	*/
	FloatType processSampleDown(const FloatType sample, bool& isSampleRet);
	/*
	Call with 2 * fs_original.
	*/
	FloatType processSampleUp(const FloatType sample);

	FloatType* processBlockDown(FloatType* const inputBlock, const int inputBlockSize);
	FloatType* processBlockUp(FloatType* const inputBlock, const int inputBlockSize);

	int getOutputBlockSize() { return mTargetBlockSize; };
	int getInputBlockSize() { return mInputBlockSize; };
private:
	bool mIsSample = false;
	double mTransitionBandwidth;
	size_t mAllpassNumberTotal;
	size_t mFilterOrder;
	std::vector<double> mCoefficients;

	FirstOrderAllpass<FloatType> mDirectPathFilters[AllpassNumber];
	FirstOrderAllpass<FloatType> mDelayPathFilters[AllpassNumber];
	Delay<FloatType> mDelay;

	int mTargetBlockSize{ 0 };
	int mInputBlockSize{ 1 };
	int mFilterBufferSize{ 0 };

	FloatType mDelayPathStorage = 0.;
	FloatType mDelayPathInput = 0.;
	FloatType mDirectPathInput = 0.;
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
inline bool HalfBandLowpass<FloatType, AllpassNumber>::init(const int expectedBlockSize, bool isDownsampling, bool isSampleBased, double transitionBandwidth)
{
	// init filters
	mTransitionBandwidth = transitionBandwidth * 2. * Pi<double>();
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
	if (!isSampleBased)
	{
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
	}
	return true;
}

template <typename FloatType, size_t AllpassNumber>
inline FloatType HalfBandLowpass<FloatType, AllpassNumber>::processSampleDown(const FloatType sample, bool& isSampleRet) {
	FloatType output = 0.;
	mDelayPathInput = sample;
	mDirectPathInput = sample;
	if (mIsSample) 
	{
		for (size_t i = 0; i < AllpassNumber; i++) 
		{
			mDirectPathInput = mDirectPathFilters[i].processSample(mDirectPathInput);
			mDelayPathStorage = mDelayPathFilters[i].processSample(mDelayPathStorage);
		}
		output = (mDirectPathInput + mDelayPathStorage) * 0.5;
	}
	mDelayPathStorage = mDelayPathInput;
	isSampleRet = mIsSample;
	mIsSample = !mIsSample;
	return output;
};
	
template <typename FloatType, size_t AllpassNumber>
inline FloatType HalfBandLowpass<FloatType, AllpassNumber>::processSampleUp(const FloatType sample) {
	FloatType output = 0;
	if (mIsSample) 
	{
		mDirectPathInput = sample;
		mDelayPathInput = sample;
		for (size_t i = 0; i < AllpassNumber; i++) 
		{
			mDirectPathInput = mDirectPathFilters[i].processSample(mDirectPathInput);
		}
		output = mDirectPathInput;
	}
	else 
	{
		for (size_t i = 0; i < AllpassNumber; i++) 
		{
			mDelayPathInput = mDelayPathFilters[i].processSample(mDelayPathInput);
		}
		output = mDelayPathInput;
	}
	mIsSample = !mIsSample;
	return output;
};

template <typename FloatType, size_t AllpassNumber>
inline FloatType* HalfBandLowpass<FloatType, AllpassNumber>::processBlockDown(FloatType* const inputBlock, const int inputBlockSize)
{
	int outCountDirect = 0;
	for (int i = 0; i < inputBlockSize; i += 2)
	{
		mDirectPathBuffer[outCountDirect] = inputBlock[i];
		outCountDirect++;
	}
	int outCountDelay = 0;
	for (int i = 1; i < inputBlockSize; i += 2)
	{
		mDelayPathBuffer[outCountDelay] = inputBlock[i];
		outCountDelay++;
	}
	mDelay.processBlock(mDelayPathBuffer.data(), mFilterBufferSize);
	for (size_t i = 0; i < AllpassNumber; i++) 
	{
		mDirectPathFilters[i].processBlock(mDirectPathBuffer.data(), mFilterBufferSize);
		mDelayPathFilters[i].processBlock(mDelayPathBuffer.data(), mFilterBufferSize);
	}
	for (int i = 0; i < mTargetBlockSize;i++)
	{
		mOutputBlock[i] = static_cast<FloatType>(0.5) * (mDirectPathBuffer[i] + mDelayPathBuffer[i]);
	}
	return mOutputBlock.data();
};

template <typename FloatType, size_t AllpassNumber>
inline FloatType* HalfBandLowpass<FloatType, AllpassNumber>::processBlockUp(FloatType* const inputBlock, const int inputBlockSize)
{
	for (int i = 0; i < inputBlockSize; i++)
	{
		mDirectPathBuffer[i] = inputBlock[i];
		mDelayPathBuffer[i] = inputBlock[i];
	}
	for (size_t i = 0; i < AllpassNumber; i++)
	{
		mDirectPathFilters[i].processBlock(mDirectPathBuffer.data(), mFilterBufferSize);
		mDelayPathFilters[i].processBlock(mDelayPathBuffer.data(), mFilterBufferSize);
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
	const double k = std::pow(std::tan((Pi<double>() - mTransitionBandwidth) / 4.), 2);
	const double k_dash = std::sqrt(1. - std::pow(k, 2));
	const double e = (1. / 2.) * ((1. - std::sqrt(k_dash)) / (1. + std::sqrt(k_dash)));
	const double q = e + 2. * std::pow(e, 5) + 15. * std::pow(e, 9.) + 150. * std::pow(e, 13.);
	// step 2
	const size_t n = mFilterOrder;
	// step 3
	const double q_1 = std::pow(q, n);
	const double k_1 = 4. * std::sqrt(q_1);
	const double d_s = std::sqrt((k_1) / (1. + k_1));
	const double d_p = 1. - std::sqrt(1. - std::pow(d_s, 2));
	// step 4
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
			delta = std::pow((-1.), m) * std::pow(q, (m * (m + 1.))) * std::sin((2. * m + 1.) * Pi<double>() * static_cast<double>(i) / static_cast<double>(n));
			num += delta;
			m += 1.;
		}
		num = 2. * std::pow(q, (1. / 4.)) * num;
		delta = 1.;
		double den = 0.;
		m = 1.;
		while (delta > 1.e-100) 
		{
			delta = std::pow((-1.), m) * std::pow(q, std::pow(m, 2.)) * std::cos(2. * m * Pi<double>() * static_cast<double>(i) / static_cast<double>(n));
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
	// step 5
	std::vector<double> a;
	for (double a_dash_i : a_dash) 
	{
		double a_i = (1. - a_dash_i) / (1. + a_dash_i);
		a.push_back(a_i);
	}
	return a;
};


/*
Config of ResamplingHandler
*/
enum class ProcessConfig
{
	Sample = 0,
	Block
};
enum class DirectionConfig
{
	Up = 0,
	Down,
	UpDown,
	DownUp
};

/**
Resampling of audio blocks or samples by a given power of 2.	
*/
template <typename FloatType, size_t AllpassNumber>
class ResamplingHandler
{
public:
	ResamplingHandler(double transitionBandwidth = 0.02);
	~ResamplingHandler();
	/**
	Initialization of the ResamplingHandler
	The configuration gives the amount of memory that will be allocated.
	DirectionConfig has to match the order the processing functions will be called from outside.
	*/
	void init(const int powToExponent = 0, const int expectedBlockSize = 128, ProcessConfig sampleConfig = ProcessConfig::Sample, DirectionConfig directionConfig = DirectionConfig::UpDown);
	/*
	Called with 2*fs target. isSampleRet indicates whether sample is relevant.
	*/
	FloatType processSampleDown(FloatType sample, bool& isSampleRet);
	/*
	Called with target samplerate.
	*/
	FloatType processSampleUp(FloatType sample);

	FloatType* processBlockDown(FloatType* const  inputBlock);
	FloatType* processBlockUp(FloatType* const inputBlock);

	int getOutputBlockSizeUp() { return mTargetBlockSizeUp; };
	int getOutputBlockSizeDown() { return mTargetBlockSizeDown; };
private:
	double mTransitionBandwidth;
	int mPowTwoFactor;
	int mResamplingFactor;
	int mInputBlockSizeDown{ 0 };
	int mInputBlockSizeUp{ 0 };
	int mTargetBlockSizeUp{ 0 };
	int mTargetBlockSizeDown{ 0 };
	std::vector<HalfBandLowpass<FloatType, AllpassNumber>*> mDownFilters;
	std::vector<HalfBandLowpass<FloatType, AllpassNumber>*> mUpFilters;
	// sample based storage
	std::vector<bool> mIsSample;
	std::vector<FloatType> mUpsamplingStorage;
};

template <typename FloatType, size_t AllpassNumber>
ResamplingHandler<FloatType, AllpassNumber>::ResamplingHandler(double transitionBandwidth)
{
	mTransitionBandwidth = transitionBandwidth;
}

template <typename FloatType, size_t AllpassNumber>
ResamplingHandler<FloatType, AllpassNumber>::~ResamplingHandler()
{
	for (auto filter : mDownFilters)
	{
		delete filter;
	}
	for (auto filter : mUpFilters)
	{
		delete filter;
	}
}

template <typename FloatType, size_t AllpassNumber>
inline void ResamplingHandler<FloatType, AllpassNumber>::init(const int powToExponent, const int expectedBlockSize, ProcessConfig sampleConfig, DirectionConfig directionConfig)
{
	for (auto filter : mDownFilters)
	{
		delete filter;
	}
	for (auto filter : mUpFilters)
	{
		delete filter;
	}
	mDownFilters.clear();
	mUpFilters.clear();
	mPowTwoFactor = powToExponent;
	mResamplingFactor = std::pow(2, powToExponent);
	for (int i = 0; i < mPowTwoFactor;i++)
	{
		if (sampleConfig == ProcessConfig::Sample)
		{
			mIsSample.push_back(false);
			mUpsamplingStorage.push_back(0.);
			mDownFilters.push_back(new HalfBandLowpass<FloatType, AllpassNumber>());
			mUpFilters.push_back(new HalfBandLowpass<FloatType, AllpassNumber>());
			mDownFilters[i]->init(1, true, true, mTransitionBandwidth);
			mUpFilters[i]->init(1, false, true, mTransitionBandwidth);
		}
		else if (sampleConfig == ProcessConfig::Block)
		{
			if (directionConfig == DirectionConfig::Down)
			{
				mDownFilters.push_back(new HalfBandLowpass<FloatType, AllpassNumber>());
				mDownFilters[i]->init(expectedBlockSize / std::pow(2, i), true, false, mTransitionBandwidth);
				mInputBlockSizeDown = expectedBlockSize;
			}
			else if (directionConfig == DirectionConfig::Up)
			{
				mUpFilters.push_back(new HalfBandLowpass<FloatType, AllpassNumber>());
				mUpFilters[i]->init(expectedBlockSize * std::pow(2, i), false, false, mTransitionBandwidth);
				mInputBlockSizeUp = expectedBlockSize;
			}
			else if (directionConfig == DirectionConfig::DownUp)
			{
				mDownFilters.push_back(new HalfBandLowpass<FloatType, AllpassNumber>());
				mUpFilters.push_back(new HalfBandLowpass<FloatType, AllpassNumber>());
				mDownFilters[i]->init(expectedBlockSize / std::pow(2, i), true, false, mTransitionBandwidth);
				mUpFilters[i]->init(expectedBlockSize / std::pow(2, mPowTwoFactor - i), false, false, mTransitionBandwidth);
				mInputBlockSizeDown = expectedBlockSize;
				mInputBlockSizeUp = expectedBlockSize / std::pow(2, mPowTwoFactor);
			}
			else if (directionConfig == DirectionConfig::UpDown)
			{
				mDownFilters.push_back(new HalfBandLowpass<FloatType, AllpassNumber>());
				mUpFilters.push_back(new HalfBandLowpass<FloatType, AllpassNumber>());
				mUpFilters[i]->init(expectedBlockSize * std::pow(2, i), false, false, mTransitionBandwidth);
				mDownFilters[i]->init(expectedBlockSize * std::pow(2, mPowTwoFactor - i), true, false, mTransitionBandwidth);
				mInputBlockSizeDown = expectedBlockSize * std::pow(2, mPowTwoFactor);
				mInputBlockSizeUp = expectedBlockSize;
			}
		}
	}
	if (sampleConfig == ProcessConfig::Block)
	{
		if (mDownFilters.size() > 0)
		{
			mTargetBlockSizeDown = mDownFilters[mPowTwoFactor - 1]->getOutputBlockSize();

		}
		else 
		{
			mTargetBlockSizeDown = expectedBlockSize;
		}
		if (mUpFilters.size() > 0)
		{
			mTargetBlockSizeUp = mUpFilters[mPowTwoFactor - 1]->getOutputBlockSize();
		}
		else
		{
			mTargetBlockSizeUp = expectedBlockSize;
		}
	}
};

template <typename FloatType, size_t AllpassNumber>
inline FloatType ResamplingHandler<FloatType, AllpassNumber>::processSampleDown(FloatType sample, bool& isSampleRet) {
	if (mPowTwoFactor > 0) {
		for (int i = 0; i < mPowTwoFactor; i++) 
		{
			mIsSample[i] = false;
		}
		for (int i = 0; i < mPowTwoFactor; i++) 
		{
			bool isSampleTmp = mIsSample[i];
			sample = mDownFilters[i]->processSampleDown(sample, isSampleTmp);
			mIsSample[i] = isSampleTmp;
			if (!mIsSample[i]) 
			{
				break;
			}
		}
		isSampleRet = mIsSample[mPowTwoFactor - 1];
		return sample;
	}
	else 
	{
		isSampleRet = true;
		return sample;
	}
};

template <typename FloatType, size_t AllpassNumber>
inline FloatType ResamplingHandler<FloatType, AllpassNumber>::processSampleUp(FloatType sample) {
	if (mPowTwoFactor > 0) {
		if (mIsSample[mPowTwoFactor - 1]) 
		{
			mUpsamplingStorage[mPowTwoFactor - 1] = sample;
		}
		for (int i = (mPowTwoFactor - 1); i > 0; i -= 1) 
		{
			if (mIsSample[i - 1]) 
			{
				mUpsamplingStorage[i - 1] = mUpFilters[i]->processSampleUp(mUpsamplingStorage[i]);
			}
		}
		sample = mUpFilters[0]->processSampleUp(mUpsamplingStorage[0]);
		return sample;
	}
	else 
	{
		return sample;
	}
};

template <typename FloatType, size_t AllpassNumber>
inline FloatType* ResamplingHandler<FloatType, AllpassNumber>::processBlockDown(FloatType* const  inputBlock)
{
	FloatType* inBlock = inputBlock;
	FloatType* outBlock = inputBlock;
	int inputBlockSize = mInputBlockSizeDown;
	for (int i = 0; i < mPowTwoFactor;i++)
	{
		outBlock = mDownFilters[i]->processBlockDown(inBlock, inputBlockSize);
		inputBlockSize = mDownFilters[i]->getOutputBlockSize();
		inBlock = outBlock;
	}
	return outBlock;
};

template <typename FloatType, size_t AllpassNumber>
inline FloatType* ResamplingHandler<FloatType, AllpassNumber>::processBlockUp(FloatType* const  inputBlock)
{
	FloatType* inBlock = inputBlock;
	FloatType* outBlock = inputBlock;
	int inputBlockSize = mInputBlockSizeUp;
	for (int i = 0; i < mPowTwoFactor;i++)
	{
		outBlock = mUpFilters[i]->processBlockUp(inBlock, inputBlockSize);
		inputBlockSize = mUpFilters[i]->getOutputBlockSize();
		inBlock = outBlock;
	}
	return outBlock;
};

}

