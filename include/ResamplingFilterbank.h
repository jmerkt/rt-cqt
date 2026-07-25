/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for more info.

 ==============================================================================
*/

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Resampling.h"
#include "../submodules/audio-utils/include/CircularBuffer.h"

/*
This class handles multirate resampling of the input data. Input callbacks may
have any size up to the size supplied to init(). Internally, callbacks are
collected into a block that is divisible by every power-of-two resampling stage.
Consequently all resampling is block based and every internal block produces at
least one sample for the lowest stage.

The supported input sample rates are power-of-two multiples of either 44.1 kHz
or 48 kHz.
*/

namespace Cqt
{

	using BufferPtr = audio_utils::CircularBuffer<double> *;

	constexpr double FilterTransitionBandwidth{0.1};
	constexpr unsigned AllpassNumber{3};

	template <int StageNumber>
	class ResamplingFilterbank
	{
		static_assert(StageNumber > 0, "A resampling filterbank needs at least one stage");
		static_assert(StageNumber < 30, "StageNumber is too large for 32-bit block sizes");

	public:
		ResamplingFilterbank() = default;
		~ResamplingFilterbank() = default;

		void init(const double samplerate, const int blockSize, const int bufferSize);

		void inputBlock(const double *const data, const int blockSize);
		double *outputBlock(const int blockSize);

		double getOriginSamplerate() const { return mOriginSamplerate; }
		int getOriginBlockSize() const { return mOriginBlockSize; }
		int getProcessingBlockSize() const { return mProcessingBlockSize; }
		int getLatencySamples() const { return mProcessingBlockSize; }
		int getLastProcessedInputSize() const { return mLastProcessedInputSize; }
		BufferPtr getStageInputBuffer(const int stage) { return &mStageInputBuffers.at(static_cast<std::size_t>(stage)); }
		BufferPtr getStageOutputBuffer(const int stage) { return &mStageOutputBuffers.at(static_cast<std::size_t>(stage)); }
		int getOriginDownsampling() const { return mOriginDownsampling; }

	private:
		static constexpr int ResamplingStageNumber = StageNumber - 1;

		static bool isPowerOfTwo(const int value);
		static int powerOfTwoExponent(const int value);
		static int roundUpToMultiple(const int value, const int multiple);

		void processInputBlock();
		void processOutputBlock();
		void pushOutputSamples(const double *const data, const int blockSize);
		void pullOutputSamples(double *const data, const int blockSize);

		ResamplingHandler<double, AllpassNumber> mInputResamplingHandler;
		std::array<HalfBandLowpass<double, AllpassNumber>, ResamplingStageNumber> mDownsamplingFilters;
		std::array<HalfBandLowpass<double, AllpassNumber>, ResamplingStageNumber> mUpsamplingFilters;
		std::array<audio_utils::CircularBuffer<double>, StageNumber> mStageInputBuffers;
		std::array<audio_utils::CircularBuffer<double>, StageNumber> mStageOutputBuffers;

		double mOriginSamplerate{48000.};
		int mOriginBlockSize{0};
		int mOriginDownsampling{0};
		int mProcessingBlockSize{0};
		int mMaximumCallbackBlockSize{0};
		int mLowestStageBlockSize{0};

		std::vector<double> mInputData;
		std::vector<double> mLowestStageOutput;
		std::vector<double> mOutputData;
		int mInputDataSize{0};
		int mLastProcessedInputSize{0};
		int mPendingOutputBlocks{0};

		// A counted ring buffer is used here because CircularBuffer cannot
		// distinguish an empty buffer from a completely full one.
		std::vector<double> mOutputQueue;
		std::size_t mOutputReadPosition{0};
		std::size_t mOutputWritePosition{0};
		std::size_t mOutputSampleCount{0};
	};

	template <int StageNumber>
	inline bool ResamplingFilterbank<StageNumber>::isPowerOfTwo(const int value)
	{
		return value > 0 && (value & (value - 1)) == 0;
	}

	template <int StageNumber>
	inline int ResamplingFilterbank<StageNumber>::powerOfTwoExponent(const int value)
	{
		assert(isPowerOfTwo(value));
		int exponent = 0;
		for (int remaining = value; remaining > 1; remaining >>= 1)
		{
			++exponent;
		}
		return exponent;
	}

	template <int StageNumber>
	inline int ResamplingFilterbank<StageNumber>::roundUpToMultiple(const int value, const int multiple)
	{
		assert(value > 0);
		assert(multiple > 0);
		const long long result =
			((static_cast<long long>(value) + multiple - 1) / multiple) * multiple;
		if (result > std::numeric_limits<int>::max())
		{
			throw std::overflow_error("The requested resampling block size is too large");
		}
		return static_cast<int>(result);
	}

	template <int StageNumber>
	inline void ResamplingFilterbank<StageNumber>::init(const double samplerate, const int blockSize, const int bufferSize)
	{
		if (!std::isfinite(samplerate) || samplerate <= 0.)
		{
			throw std::invalid_argument("The sample rate must be finite and positive");
		}
		if (blockSize <= 0 || bufferSize <= 0)
		{
			throw std::invalid_argument("Block and buffer sizes must be positive");
		}

		const long long roundedSamplerate = std::llround(samplerate);
		if (std::abs(samplerate - static_cast<double>(roundedSamplerate)) > 1.e-6 ||
			roundedSamplerate > std::numeric_limits<int>::max())
		{
			throw std::invalid_argument("Only integer sample rates are supported");
		}

		const int samplerateInt = static_cast<int>(roundedSamplerate);
		int originFactor = 0;
		if ((samplerateInt % 44100) == 0 && isPowerOfTwo(samplerateInt / 44100))
		{
			mOriginSamplerate = 44100.;
			originFactor = samplerateInt / 44100;
		}
		else if ((samplerateInt % 48000) == 0 && isPowerOfTwo(samplerateInt / 48000))
		{
			mOriginSamplerate = 48000.;
			originFactor = samplerateInt / 48000;
		}
		else
		{
			throw std::invalid_argument(
				"The sample rate must be a power-of-two multiple of 44.1 kHz or 48 kHz");
		}

		mOriginDownsampling = powerOfTwoExponent(originFactor);
		const int filterbankFactor = 1 << ResamplingStageNumber;
		if (originFactor > (std::numeric_limits<int>::max() / filterbankFactor))
		{
			throw std::overflow_error("The resampling factor is too large");
		}
		const int inputAlignment = originFactor * filterbankFactor;

		mMaximumCallbackBlockSize = blockSize;
		mProcessingBlockSize = roundUpToMultiple(blockSize, inputAlignment);
		mOriginBlockSize = mProcessingBlockSize / originFactor;
		mLowestStageBlockSize = mOriginBlockSize / filterbankFactor;

		mInputData.assign(static_cast<std::size_t>(mProcessingBlockSize), 0.);
		mLowestStageOutput.assign(static_cast<std::size_t>(mLowestStageBlockSize), 0.);
		mOutputData.assign(static_cast<std::size_t>(mMaximumCallbackBlockSize), 0.);
		mInputDataSize = 0;
		mLastProcessedInputSize = 0;
		mPendingOutputBlocks = 0;

		mInputResamplingHandler.init(
			mOriginDownsampling,
			mProcessingBlockSize,
			DirectionConfig::DownUp);

		for (int stage = 0; stage < ResamplingStageNumber; ++stage)
		{
			const int stageInputSize = mOriginBlockSize / (1 << stage);
			const int stageOutputSize = stageInputSize / 2;
			mDownsamplingFilters[static_cast<std::size_t>(stage)].init(
				stageInputSize, true, FilterTransitionBandwidth);
			mUpsamplingFilters[static_cast<std::size_t>(stage)].init(
				stageOutputSize, false, FilterTransitionBandwidth);
		}

		for (int stage = 0; stage < StageNumber; ++stage)
		{
			const int stageBlockSize = mOriginBlockSize / (1 << stage);
			const int requiredBufferSize = std::max(bufferSize, stageBlockSize * 2);
			mStageInputBuffers[static_cast<std::size_t>(stage)].changeSize(requiredBufferSize);
			mStageOutputBuffers[static_cast<std::size_t>(stage)].changeSize(requiredBufferSize);
		}

		const std::size_t outputQueueCapacity =
			static_cast<std::size_t>(mProcessingBlockSize) * 2U;
		mOutputQueue.assign(outputQueueCapacity, 0.);
		mOutputReadPosition = 0;
		mOutputWritePosition = static_cast<std::size_t>(mProcessingBlockSize);
		mOutputSampleCount = static_cast<std::size_t>(mProcessingBlockSize);
	}

	template <int StageNumber>
	inline void ResamplingFilterbank<StageNumber>::processInputBlock()
	{
		double *dataIn = mInputResamplingHandler.processBlockDown(mInputData.data());
		int dataSize = mOriginBlockSize;
		mStageInputBuffers[0].pushBlock(dataIn, dataSize);

		for (int stage = 0; stage < ResamplingStageNumber; ++stage)
		{
			dataIn = mDownsamplingFilters[static_cast<std::size_t>(stage)].processBlockDown(dataIn);
			dataSize /= 2;
			mStageInputBuffers[static_cast<std::size_t>(stage + 1)].pushBlock(dataIn, dataSize);
		}

		++mPendingOutputBlocks;
		mLastProcessedInputSize += mProcessingBlockSize;
	}

	template <int StageNumber>
	inline void ResamplingFilterbank<StageNumber>::inputBlock(const double *const data, const int blockSize)
	{
		if (blockSize < 0 || blockSize > mMaximumCallbackBlockSize ||
			(data == nullptr && blockSize > 0))
		{
			throw std::invalid_argument(
				"The input block must not exceed the callback size supplied to init()");
		}

		mLastProcessedInputSize = 0;
		int inputPosition = 0;
		while (inputPosition < blockSize)
		{
			const int samplesToCopy =
				std::min(blockSize - inputPosition, mProcessingBlockSize - mInputDataSize);
			std::copy_n(
				data + inputPosition,
				samplesToCopy,
				mInputData.data() + mInputDataSize);
			inputPosition += samplesToCopy;
			mInputDataSize += samplesToCopy;

			if (mInputDataSize == mProcessingBlockSize)
			{
				processInputBlock();
				mInputDataSize = 0;
			}
		}
	}

	template <int StageNumber>
	inline void ResamplingFilterbank<StageNumber>::pushOutputSamples(const double *const data, const int blockSize)
	{
		if (mOutputSampleCount + static_cast<std::size_t>(blockSize) > mOutputQueue.size())
		{
			throw std::logic_error(
				"Too many input blocks were submitted without consuming output");
		}
		for (int sample = 0; sample < blockSize; ++sample)
		{
			mOutputQueue[mOutputWritePosition] = data[sample];
			mOutputWritePosition = (mOutputWritePosition + 1U) % mOutputQueue.size();
		}
		mOutputSampleCount += static_cast<std::size_t>(blockSize);
	}

	template <int StageNumber>
	inline void ResamplingFilterbank<StageNumber>::pullOutputSamples(double *const data, const int blockSize)
	{
		if (mOutputSampleCount < static_cast<std::size_t>(blockSize))
		{
			throw std::logic_error(
				"More output was requested than the filterbank has buffered");
		}
		for (int sample = 0; sample < blockSize; ++sample)
		{
			data[sample] = mOutputQueue[mOutputReadPosition];
			mOutputReadPosition = (mOutputReadPosition + 1U) % mOutputQueue.size();
		}
		mOutputSampleCount -= static_cast<std::size_t>(blockSize);
	}

	template <int StageNumber>
	inline void ResamplingFilterbank<StageNumber>::processOutputBlock()
	{
		mStageOutputBuffers[StageNumber - 1].pullBlock(
			mLowestStageOutput.data(), mLowestStageBlockSize);
		double *dataOut = mLowestStageOutput.data();

		for (int stage = ResamplingStageNumber - 1; stage >= 0; --stage)
		{
			dataOut = mUpsamplingFilters[static_cast<std::size_t>(stage)].processBlockUp(dataOut);
			const int stageBlockSize = mOriginBlockSize / (1 << stage);
			mStageOutputBuffers[static_cast<std::size_t>(stage)].pullBlockAdd(
				dataOut, stageBlockSize);
		}

		dataOut = mInputResamplingHandler.processBlockUp(dataOut);
		pushOutputSamples(dataOut, mProcessingBlockSize);
	}

	template <int StageNumber>
	inline double *ResamplingFilterbank<StageNumber>::outputBlock(const int blockSize)
	{
		if (blockSize < 0 || blockSize > mMaximumCallbackBlockSize)
		{
			throw std::invalid_argument(
				"The output block must not exceed the callback size supplied to init()");
		}

		while (mPendingOutputBlocks > 0)
		{
			processOutputBlock();
			--mPendingOutputBlocks;
		}

		pullOutputSamples(mOutputData.data(), blockSize);
		return mOutputData.data();
	}

}
