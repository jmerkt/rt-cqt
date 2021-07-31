/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

#include <vector>
#include "Resampling.h"
#include "CircularBuffer.h"

/*
This class handles multirate resampling of the input data. First to 44.1/48 kHz and from there to the rates of the different stages.
The incoming blocks are downsampled by 2 until the remaining sample number cannot be divided by 2 anymore. From there sample-based downsampling is used. 
This method should cover all possible input block sizes.

Note: at the moment only (samplerates > 44.1 kHz) and (((fsInt % 44100) == 0) || ((fsInt % 48000) == 0)) are supported. 
*/

namespace Cqt
{

typedef CircularBuffer<double>* BufferPtr;

template<int StageNumber>
class ResamplingFilterbank
{
public:
	ResamplingFilterbank() = default;
	~ResamplingFilterbank() = default;

	void init(const double samplerate, const int blockSize, const int bufferSize);

	void inputBlock(double* const data, const int blockSize);
	double* outputBlock(const int blockSize);

	double getOriginSamplerate() { return mOriginSamplerate; };
	BufferPtr getStageInputBuffer(const int stage) { return &mStageInputBuffers[stage]; };
	BufferPtr getStageOutputBuffer(const int stage) { return &mStageOutputBuffers[stage]; };
	int getOriginDownsampling() { return mOriginDownsampling; };
protected:
	ResamplingHandler<double, 3> mInputResamplingHandler;
	HalfBandLowpass<double, 3> mDownsamplingFilters[StageNumber - 1];
	HalfBandLowpass<double, 3> mUpsamplingFilters[StageNumber - 1];
	CircularBuffer<double> mStageInputBuffers[StageNumber];
	CircularBuffer<double> mStageOutputBuffers[StageNumber];

	double mOriginSamplerate;
	int mOriginBlockSize;
	int mOriginDownsampling;
	double mStageSamplerates[StageNumber];
	int mDownsamplingBlockSizes[StageNumber - 1];
	int mUpsamplingBlockSizes[StageNumber - 1];

	// block based
	int mBlockFilterNumber;

	// sample based
	int mSampleInputSize{ 0 };
	int mSampleFilterNumber;
	std::vector<std::vector<bool>> mIsSample;
	std::vector<double> mUpsamplingSampleBuffer;
	std::vector<double> mUpsamplingSampleOutputBuffer;

	// handling of input / output buffering
	CircularBuffer<double> mInputBuffer;
	CircularBuffer<double> mOutputBuffer;
	std::vector<double> mInputData;
	std::vector<double> mOutputData;
	size_t mInputDataCounter{ 0 };
	size_t mOutputDataCounter{ 0 };
	int mExpectedBlockSize{ 0 };
};

template<int StageNumber>
inline void ResamplingFilterbank<StageNumber>::init(const double samplerate, const int blockSize, const int bufferSize)
{
	// handling of input / output buffering
	mInputBuffer.changeSize(blockSize * 2);
	mOutputBuffer.changeSize(blockSize * 2);
	mInputData.resize(blockSize, 0.);
	mOutputData.resize(blockSize, 0.);
	mInputDataCounter = 0;
	mOutputDataCounter = blockSize;
	mExpectedBlockSize = blockSize;

	mOriginDownsampling = 0;
	int fsInt = static_cast<int>(samplerate);

	//assert(((fsInt % 44100) == 0) || ((fsInt % 48000) == 0));
	if ((fsInt % 44100) == 0)
	{
		mOriginDownsampling = std::log2(fsInt / 44100);
		mOriginSamplerate = 44100.;
	}
	else if ((fsInt % 48000) == 0)
	{
		mOriginDownsampling = std::log2(fsInt / 48000);
		mOriginSamplerate = 48000.;
	}
	else
	{
		mOriginDownsampling = std::log2(fsInt / 44100);
		mOriginSamplerate = 44100.;
	}

	// init origin downsampling
	mInputResamplingHandler.init(mOriginDownsampling, blockSize, ProcessConfig::Block, DirectionConfig::DownUp);
	mOriginBlockSize = mInputResamplingHandler.getOutputBlockSizeDown();
	for (int i = 0; i < StageNumber; i++)
	{
		mStageSamplerates[i] = mOriginSamplerate / pow(2., i);
	}

	// configure internal stages
	for (int i = 0; i < (StageNumber - 1); i++)
	{
		mDownsamplingBlockSizes[i] = 0;
		mUpsamplingBlockSizes[i] = 0;
	}
	int originBlockSize = mOriginBlockSize;
	mBlockFilterNumber = 0;
	mSampleFilterNumber = 0;
	while ((originBlockSize > 1) && (mBlockFilterNumber < (StageNumber - 1)) && ((originBlockSize % 2) == 0))
	{
		mDownsamplingBlockSizes[mBlockFilterNumber] = originBlockSize;
		originBlockSize /= 2;
		mUpsamplingBlockSizes[mBlockFilterNumber] = originBlockSize;
		mBlockFilterNumber++;
	}
	// init Filters
	for (int stage = 0; stage < (StageNumber - 1); stage++)
	{
		if (stage < mBlockFilterNumber)
		{
			mDownsamplingFilters[stage].init(mOriginBlockSize / std::pow(2, stage), true, false, 0.02);
			mUpsamplingFilters[stage].init(mOriginBlockSize / std::pow(2, stage + 1), false, false, 0.02);
		}
		else
		{
			mDownsamplingFilters[stage].init(mOriginBlockSize / std::pow(2, stage), true, true, 0.02);
			mUpsamplingFilters[stage].init(mOriginBlockSize / std::pow(2, stage + 1), false, true, 0.02);
		}
	}
	// buffers for samplebased / block based intermediate processing
	mIsSample.clear();
	mUpsamplingSampleBuffer.clear();
	if ((StageNumber - 1) > mBlockFilterNumber)
	{
		mSampleInputSize = originBlockSize;
		mUpsamplingSampleOutputBuffer.resize(originBlockSize, 0.);
		mSampleFilterNumber = (StageNumber - 1) - mBlockFilterNumber;
		mUpsamplingSampleBuffer.resize(mSampleFilterNumber, 0.);
		mIsSample.resize(originBlockSize);
		for (int i = 0; i < originBlockSize; i++)
		{
			mIsSample[i].resize(mSampleFilterNumber, false);
		}
	}
	else
	{
		mUpsamplingSampleOutputBuffer.resize(mUpsamplingFilters[mBlockFilterNumber - 1].getInputBlockSize(), 0.);
		mSampleInputSize = 0;
	}
	// init Buffers
	for (int stage = 0; stage < StageNumber; stage++)
	{
		mStageInputBuffers[stage].changeSize(bufferSize);
		mStageOutputBuffers[stage].changeSize(bufferSize);
	}
}

template<int StageNumber>
inline void ResamplingFilterbank<StageNumber>::inputBlock(double* const data, const int blockSize)
{
	mInputBuffer.pushBlock(data, blockSize);
	mInputDataCounter += blockSize;
	if (mInputDataCounter >= mExpectedBlockSize)
	{
		mInputDataCounter -= mExpectedBlockSize;
		mInputBuffer.pullDelayBlock(mInputData.data(), mExpectedBlockSize + mInputDataCounter - 1, mExpectedBlockSize);
		// Filterbank
		int stageIdx = 0;
		int filterIdx = 0;
		// input Resampler to FsOrigin
		double* dataIn = mInputResamplingHandler.processBlockDown(mInputData.data());
		int blockSize = mInputResamplingHandler.getOutputBlockSizeDown();
		// save origin data
		mStageInputBuffers[stageIdx].pushBlock(dataIn, blockSize);
		stageIdx++;
		// process block based stages
		for (int i = 0; i < mBlockFilterNumber; i++)
		{
			dataIn = mDownsamplingFilters[filterIdx].processBlockDown(dataIn, blockSize);
			blockSize = mDownsamplingFilters[filterIdx].getOutputBlockSize();
			// save stage data
			mStageInputBuffers[stageIdx].pushBlock(dataIn, blockSize);
			stageIdx++;
			filterIdx++;
		}

		// handle sample based stages, where block based is not possible anymore
		int stageIdxSave = stageIdx;
		int filterIdxSave = filterIdx;
		for (int s = 0; s < mSampleInputSize; s++)
		{
			stageIdx = stageIdxSave;
			filterIdx = filterIdxSave;
			double sampleIn = dataIn[s];
			std::fill(mIsSample[s].begin(), mIsSample[s].end(), false);
			for (int i = 0; i < mSampleFilterNumber; i++)
			{
				bool isSampleTmp = mIsSample[s][i];
				sampleIn = mDownsamplingFilters[filterIdx].processSampleDown(sampleIn, isSampleTmp);
				mIsSample[s][i] = isSampleTmp;
				if (mIsSample[s][i])
				{
					mStageInputBuffers[stageIdx].pushSample(sampleIn);
				}
				if (!mIsSample[s][i])
				{
					break;
				}
				filterIdx++;
				stageIdx++;
			}
		}
	}
}

template<int StageNumber>
inline double* ResamplingFilterbank<StageNumber>::outputBlock(const int blockSize)
{
	mOutputDataCounter += blockSize;
	if (mOutputDataCounter >= mExpectedBlockSize)
	{
		mOutputDataCounter -= mExpectedBlockSize;
		//Filterbank
		int stageIdx = StageNumber - 1;
		int filterIdx = StageNumber - 2;
		// handle sample based stages
		for (int s = 0; s < mSampleInputSize; s++)
		{
			stageIdx = StageNumber - 1;
			filterIdx = StageNumber - 2;
			// lowest Stage
			if (mSampleFilterNumber > 0)
			{
				if (mIsSample[s][(mSampleFilterNumber - 1)])
				{
					mUpsamplingSampleBuffer[(mSampleFilterNumber - 1)] = mStageOutputBuffers[stageIdx].pullSample();
				}
				stageIdx--;
			}
			// rest of stages
			for (int i = (mSampleFilterNumber - 2); i >= 0; i -= 1)
			{
				if (mIsSample[s][i])
				{
					// pull own sample and store - add together with upsampled lower level
					mUpsamplingSampleBuffer[i] = mStageOutputBuffers[stageIdx].pullSample() + mUpsamplingFilters[filterIdx].processSampleUp(mUpsamplingSampleBuffer[i + 1]);
				}
				stageIdx--;
				filterIdx--;
			}
			// highest sample based stage
			if (mSampleFilterNumber > 0)
			{
				mUpsamplingSampleOutputBuffer[s] = mStageOutputBuffers[stageIdx].pullSample() + mUpsamplingFilters[filterIdx].processSampleUp(mUpsamplingSampleBuffer[0]);
				stageIdx--;
				filterIdx--;
			}
		}

		// get output of last sample based stage
		int inputBlockSize = mUpsamplingFilters[filterIdx].getInputBlockSize();
		double* inBlock = mUpsamplingSampleOutputBuffer.data();
		if (mSampleFilterNumber <= 0)
		{
			mStageOutputBuffers[stageIdx].pullBlock(mUpsamplingSampleOutputBuffer.data(), inputBlockSize);
			stageIdx--;
		}
		// handle block based stages 
		for (int i = 0; i < mBlockFilterNumber; i++)
		{
			inBlock = mUpsamplingFilters[filterIdx].processBlockUp(inBlock, inputBlockSize);
			inputBlockSize = mUpsamplingFilters[filterIdx].getOutputBlockSize();
			filterIdx--;
			// combine upsampled data
			mStageOutputBuffers[stageIdx].pullBlockAdd(inBlock, inputBlockSize);
			stageIdx--;
		}

		// input Resampler to FsOrigin
		double* dataOut = mInputResamplingHandler.processBlockUp(inBlock);
		// store output
		mOutputBuffer.pushBlock(dataOut, mExpectedBlockSize);
	}
	mOutputBuffer.pullBlock(mOutputData.data(), blockSize);
	return mOutputData.data();
}

}










