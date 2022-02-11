/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

/*
* Limitations:
* - Fixed Hann window function - for resynthesis a HopSize >= Fft_Size / 2 should be used - HopSize == Fft_Size will produce NaNs
* - No Zero Padding as of now
*/

#include "CircularBuffer.h"
#include "ResamplingFilterbank.h"
#include "Utils.h"
#include <atomic>
#include <memory>

#define SIMD_SZ 1
#define PFFFT_ENABLE_DOUBLE
#include "../submodules/pffft/pffft.hpp"
#include <complex>

namespace Cqt 
{

using namespace std::complex_literals;
constexpr int Fft_Size{ 512 };
constexpr int Fft_Domain_Size{ Fft_Size / 2 };
constexpr double KernelThreshold{ 1.e-4 };
constexpr double WindowEnergyLossCompensation{ 1.63 }; // Fixed for Hanning window as of now
constexpr double WindowAmplitudeLossCompensation{ 2.0 }; // Fixed for Hanning window as of now

typedef pffft::AlignedVector<pffft::Types<double>::Complex> CplxVector;
typedef pffft::AlignedVector<double> RealVector;
typedef RealVector TimeDataType;
typedef CplxVector CqtBufferType;


/*
Structure to schedule transformation timings.
*/
struct ScheduleElement
{
	int sample{ 0 }; // Position in actual block
	int octave{ 0 };
	int delayOctaveRate{ 0 }; // Delay of the transformation in samplesin the corresponding octave sample buffer
};

/*
Handles the transformations of the single octaves.
*/
template <int B>
class TransformationHandler  
{
public:
	TransformationHandler ();
	~TransformationHandler() = default;

	void init(const int hopSize);
	void initBuffers(BufferPtr inputBuffer = nullptr, BufferPtr outputBuffer = nullptr);
	void initKernels(const CplxVector* const kernelArray, const CplxVector* const kernelArrayInverse, const std::vector<int>* const kernelMask, const std::vector<int>* const kernelMaskInv);
	void initFs(const int blockSize);

	void cqtTransform(const ScheduleElement schedule);
	void icqtTransform(const ScheduleElement schedule);

	static void calculateWindow(double* const windowData, const int size);
	static void calculateInverseWindow(double* const windowData, double* const invWindowData, const int size, const int hopSize);

	inline CqtBufferType* getCqtBuffer() { return &mCqtBuffer; }; 
	inline BufferPtr getOutputBuffer() { return mStageOutputBuffer; };
private:
	double mWindow[Fft_Size];
	double mInvWindow[Fft_Size];
	// kernel storage
	CplxVector mKernelArray[B]; // mB x Fft_Domain_Size
	CplxVector mKernelArrayInverse[B]; // mB x Fft_Domain_Size
	std::vector<int> mKernelMask[B]; // mB x Fft_Domain_Size
	std::vector<int> mKernelMaskInv[B]; // mB x Fft_Domain_Size

	// scaling 
	const double mFftScalingFactor{ 1. / std::sqrt(static_cast<double>(Fft_Size)) };
	// pffft
	pffft::Fft<double> mFft;
	RealVector mFftInputBuffer;
	CplxVector mFftBuffer;
	CplxVector mIfftInputBuffer;
	CqtBufferType mCqtBuffer; // mB
	RealVector mOutputBuffer;
	RealVector mIfftOutputBuffer;

	// input / output buffers
	BufferPtr mStageInputBuffer;
	BufferPtr mStageOutputBuffer;
};

template <int B>
TransformationHandler<B>::TransformationHandler() :
	mCqtBuffer(B),
	mFft(Fft_Size)
{
	// hard-coded Hann window as of now
	calculateWindow(mWindow, Fft_Size);
	// fft and buffers
	mFftInputBuffer = mFft.valueVector();
	mFftBuffer = mFft.spectrumVector();
	mIfftInputBuffer = mFft.spectrumVector();
	mIfftOutputBuffer = mFft.valueVector();
	std::fill(mIfftOutputBuffer.begin(), mIfftOutputBuffer.end(), 0.);
	mOutputBuffer = mFft.valueVector();
	std::fill(mOutputBuffer.begin(), mOutputBuffer.end(), 0.);
	for (int tone = 0; tone < B; tone++) 
	{
		mCqtBuffer[tone] = 0. + 0.i;
	}
	// kernels
	for (int tone = 0; tone < B; tone++) 
	{
		mKernelArray[tone] = mFft.spectrumVector();
		mKernelArrayInverse[tone] = mFft.spectrumVector();
	}
}

template <int B>
inline void TransformationHandler<B>::init(const int hopSize)
{
	calculateInverseWindow(mWindow, mInvWindow, Fft_Size, hopSize);
}

template <int B>
inline void TransformationHandler<B>::initBuffers(BufferPtr inputBuffer, BufferPtr outputBuffer)
{
	mStageInputBuffer = inputBuffer;
	mStageOutputBuffer = outputBuffer;
}

template <int B>
inline void TransformationHandler<B>::initKernels(const CplxVector* const kernelArray, const CplxVector* const kernelArrayInverse, const std::vector<int>* const kernelMask, const std::vector<int>* const kernelMaskInv)
{
	for (int tone = 0; tone < B; tone++) 
	{
		for (int i = 0; i < Fft_Domain_Size; i++) 
		{
			mKernelArray[tone][i] = kernelArray[tone].at(i);
			mKernelArrayInverse[tone][i] = kernelArrayInverse[tone].at(i);
		}
		mKernelMask[tone].resize(kernelMask[tone].size(), 0);
		for (size_t i = 0; i < kernelMask[tone].size(); i++)
		{
			mKernelMask[tone][i] = kernelMask[tone][i];
		}
		mKernelMaskInv[tone].resize(kernelMaskInv[tone].size(), 0);
		for (size_t i = 0; i < kernelMaskInv[tone].size(); i++)
		{
			mKernelMaskInv[tone][i] = kernelMaskInv[tone][i];
		}
	}
};

template <int B>
inline void TransformationHandler<B>::initFs(const int blockSize)
{
	mOutputBuffer.resize(Fft_Size + blockSize);
}

template <int B>
inline void TransformationHandler<B>::cqtTransform(const ScheduleElement schedule)
{
	// collect the fft input data
	mStageInputBuffer->pullDelayBlock(mFftInputBuffer.data(), static_cast<int>(Fft_Size) + schedule.delayOctaveRate - 1, static_cast<int>(Fft_Size)); 
	// apply time window
	for (int i = 0; i < Fft_Size; i++)
	{
		mFftInputBuffer[i] *= mWindow[i];
	}
	// fft
	mFft.forward(mFftInputBuffer, mFftBuffer);
	// scale data
	for (int i = 0; i < Fft_Domain_Size; i++) 
	{
		mFftBuffer[i] *= mFftScalingFactor;
	}
	// kernel multipications
	for (int tone = 0; tone < B; tone++) 
	{
		mCqtBuffer[tone] = 0. + 0.i;
		for (size_t i = 0; i < mKernelMask[tone].size(); i++)
		{
			const int idx = mKernelMask[tone][i];
			mCqtBuffer[tone] += mFftBuffer[idx] * mKernelArray[tone][idx];
		}
	}
};

template <int B>
inline void TransformationHandler<B>::icqtTransform(const ScheduleElement schedule)
{
	// kernel multiplications
	for (int i = 0; i < Fft_Domain_Size; i++)
	{
		mIfftInputBuffer[i] = 0. + 0.i;
	}
	for (int tone = 0; tone < B; tone++)
	{
		for (size_t i = 0; i < mKernelMaskInv[tone].size(); i++)
		{
			const int idx = mKernelMaskInv[tone][i];
			mIfftInputBuffer[idx] += mCqtBuffer[tone] * mKernelArrayInverse[tone][idx];
		}
	}
	// ifft
	mFft.inverse(mIfftInputBuffer, mIfftOutputBuffer);
	// scale and window data
	for (int i = 0; i < Fft_Size; i++) 
	{
		mIfftOutputBuffer[i] *= mFftScalingFactor;
	}
	for (int i = 0; i < Fft_Size; i++)
	{
		mIfftOutputBuffer[i] *= mInvWindow[i];
	}

	// overlap-add
	std::fill(mOutputBuffer.begin(), mOutputBuffer.end(), 0.);
	//// pull whats left from the previous transform
	mStageOutputBuffer->pullBlock(mOutputBuffer.data(), mStageOutputBuffer->getWriteReadDistance() ); 
	//// add new data
	int count = 0; 
	for (int i = schedule.delayOctaveRate; i < (schedule.delayOctaveRate + Fft_Size); i++)
	{
		mOutputBuffer[i] += mIfftOutputBuffer[count];
		count++;
	}
	//// push new data 
	mStageOutputBuffer->pushBlock(mOutputBuffer.data(), Fft_Size + schedule.delayOctaveRate);
};

template <int B>
inline void TransformationHandler<B>::calculateWindow(double* const windowData, const int size)
{
	for (int i = 0; i < size; i++)
	{
		windowData[i] = (1. / 2.) * (1. - std::cos((2. * Pi<double>() * static_cast<double>(i)) / static_cast<double>(size - 1)));
	}
}

template <int B>
inline void TransformationHandler<B>::calculateInverseWindow(double* const windowData, double* const invWindowData, const int size, const int hopSize)
{
	std::vector<double> windowTmp(size, 0.);
	for (int i = 0; i < size; i += hopSize)
	{
		for (int j = 0; j < size; j++)
		{
			windowTmp[(i + j) % size] += windowData[j] * windowData[j];
		}
	}
	for (int i = 0; i < size; i++)
	{
		invWindowData[i] = windowData[i] / windowTmp[i] * std::pow(WindowEnergyLossCompensation, 2);
	}
}

/*
Main CQT class
*/
template <int B, int OctaveNumber>
class ConstantQTransform 
{
public:
	ConstantQTransform();
	~ConstantQTransform() = default;

	void init(int hopSize);
	void init(std::vector<int> octaveHopSizes);
	void initFs(double fs, const int blockSize);
	void setConcertPitch(double concertPitch);
	inline void recalculateKernels() { mNewKernels.store(true); };

	void inputBlock(double* const data, const int blockSize);
	double* outputBlock(const int blockSize);
	void cqt(const ScheduleElement schedule);
	void icqt(const ScheduleElement schedule);

	inline std::vector<ScheduleElement>& getCqtSchedule() { return mCqtSchedule; };
	inline CqtBufferType* getOctaveCqtBuffer(const int octave) { return mTransformationHandlers[octave].getCqtBuffer(); };
	inline BufferPtr getOctaveOutputBuffer(const int octave) { return mTransformationHandlers[octave].getOutputBuffer(); };
	inline int getHopSize(const int octave) { return mHopSizes[octave]; };
	inline size_t getLatencySamples(const int octave) { return mLatencySamples[octave]; };
	inline double getLatencyMs(const int octave) { return mLatencyMs[octave]; };
	inline double getOctaveSampleRate(const int octave) { return mSampleRates[octave]; };
	inline std::vector<std::vector<double>>& getKernelFreqs() { return mKernelFreqs; };
	inline std::vector<std::vector<double>>& getKernelFreqsInv() { return mKernelFreqsInv; };
	inline void resetKernelFreqs() { initKernelFreqs(); };

protected:
	void initKernelFreqs();
	void calculateKernels();

	double mConcertPitch{ 440. };
	int mBinNumber;
	int mOverlaps[OctaveNumber];
	double mLatencyMs[OctaveNumber];
	size_t mLatencySamples[OctaveNumber];
	int mHopSizes[OctaveNumber];
	double mFs;
	double mSampleRates[OctaveNumber];
	double mSampleRatesByOriginRate[OctaveNumber];

	TransformationHandler<B> mTransformationHandlers[OctaveNumber];
	ResamplingFilterbank<OctaveNumber> mFilterbank;
	size_t mSampleCounters[OctaveNumber];

	std::vector<ScheduleElement> mCqtSchedule;

	pffft::Fft<std::complex<double>>  mFft;
	pffft::Fft<std::complex<double>>  mFftInv;
	pffft::Fft<double>  mFftAllocation;
	CplxVector mFftTmpStorage;
	CplxVector mFftTmpStorageInv;

	CplxVector mKernelStorage[B];
	CplxVector mKernelStorageInv[B];
	std::vector<int> mKernelMask[B];
	std::vector<int> mKernelMaskInv[B];
	CplxVector mKernelStorageTime[B];
	CplxVector mKernelStorageTimeInv[B];
	RealVector window;
	std::atomic<bool> mNewKernels{ true };

	std::vector<std::vector<double>> mKernelFreqs;
	std::vector<std::vector<double>> mKernelFreqsInv;
	
};

template <int B, int OctaveNumber>
ConstantQTransform<B, OctaveNumber>::ConstantQTransform():
	mFft(Fft_Size), mFftInv(Fft_Size), mFftAllocation(Fft_Size)
{
	mBinNumber = B * OctaveNumber;

	mFftTmpStorage = mFft.spectrumVector();
	mFftTmpStorageInv = mFft.spectrumVector();
	// configure all the buffer sizes
	for (int tone = 0; tone < B; tone++) 
	{
		mKernelStorageTime[tone] = mFft.valueVector();
		mKernelStorageTimeInv[tone] = mFft.valueVector();
		mKernelStorage[tone] = mFftAllocation.spectrumVector();
		mKernelStorageInv[tone] = mFftAllocation.spectrumVector();
	}
	// generate window function
	window.resize(Fft_Size);
	TransformationHandler<B>::calculateWindow(window.data(), Fft_Size);
	// transformation in/out buffers
	mKernelFreqs.resize(OctaveNumber);
	mKernelFreqsInv.resize(OctaveNumber);
	for (int octave = 0; octave < OctaveNumber; octave++)
	{
		mTransformationHandlers[octave].initBuffers(mFilterbank.getStageInputBuffer(octave), mFilterbank.getStageOutputBuffer(octave));
		mKernelFreqs[octave].resize(B, 0.);
		mKernelFreqsInv[octave].resize(B, 0.);
		mSampleCounters[octave] = 0;
	}
};

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::init(int hopSize)
{
	initKernelFreqs();

	hopSize = Clip<int>(hopSize, 1, Fft_Size);
	for (int octave = 0; octave < OctaveNumber; octave++) 
	{
		mHopSizes[octave] = hopSize;
		mOverlaps[octave] = Fft_Size - hopSize;
	}
	for (int octave = 0; octave < OctaveNumber; octave++)
	{
		mTransformationHandlers[octave].init(mHopSizes[octave]);
	}
}

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::init(std::vector<int> octaveHopSizes)
{
	initKernelFreqs();

	assert(octaveHopSizes.size() == OctaveNumber);
	for (int octave = 0; octave < OctaveNumber; octave++)
	{
		int hopSize = Clip<int>(octaveHopSizes.at(octave), 1, Fft_Size);
		mHopSizes[octave] = hopSize;
		mOverlaps[octave] = Fft_Size - hopSize;
	}
	for (int octave = 0; octave < OctaveNumber; octave++)
	{
		mTransformationHandlers[octave].init(mHopSizes[octave]);
	}
}

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::initFs(double fs, const int blockSize)
{
	mFilterbank.init(fs, blockSize, blockSize + Fft_Size);
	mFs = mFilterbank.getOriginSamplerate();

	for (int octave = 0; octave < OctaveNumber; octave++) 
	{
		// latency per octave
		mLatencySamples[octave] = static_cast<size_t>(mHopSizes[octave]) * static_cast<size_t>(std::pow(2, octave)) * static_cast<size_t>(std::pow(2, mFilterbank.getOriginDownsampling()));
		mSampleCounters[octave] = mLatencySamples[octave];
		// samplerates
		mSampleRates[octave] = mFs / std::pow(2., octave);
		mLatencyMs[octave] = static_cast<double>(mHopSizes[octave]) / mSampleRates[octave] * 1000.;
		mSampleRatesByOriginRate[octave] = mSampleRates[octave] / mSampleRates[0];
	}
	for (int octave = 0; octave < OctaveNumber; octave++) 
	{
		mTransformationHandlers[octave].initFs(blockSize);
	}
	// calc the windows and give em to handlers
	recalculateKernels();
};

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::initKernelFreqs() 
{
	/*
	https://en.wikipedia.org/wiki/Piano_key_frequencies
	f(n) = 2^((n-49)/12) * mConcertPitch
	Range in highest octave from n = [100, 111] ~ [8.37 kHz, 15.804 kHz]
	*/
	const double fRef = std::pow(2., ((100. - 49.) / 12.)) * mConcertPitch;
	for (int octave = 0; octave < OctaveNumber; octave++) 
	{
		for (int tone = 0; tone < B; tone++) 
		{
			mKernelFreqs[octave][tone] = (fRef / std::pow(2., (octave + 1))) * std::pow(2., static_cast<double>(B + tone) / static_cast<double>(B));
			mKernelFreqsInv[octave][tone] = mKernelFreqs[octave][tone];
		}
	}
};

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::calculateKernels()
{
	// calculate the time domain kernels
	for (int k = 0; k < B; k++) 
	{
		const double fk = mKernelFreqs[0][k];
		const double fkInv = mKernelFreqsInv[0][k];
		for (int n = 0; n < Fft_Size; n++) 
		{
			mKernelStorageTime[k][n] = std::conj((1. / static_cast<double>(Fft_Size)) * window[n] * std::exp(-1i * 2. * Pi<double>() * static_cast<double>(n) * (fk / mSampleRates[0])));
			mKernelStorageTimeInv[k][n] = std::conj((1. / static_cast<double>(Fft_Size)) * window[n] * std::exp(-1i * 2. * Pi<double>() * static_cast<double>(n) * (fkInv / mSampleRates[0])));
		}
	}
	// fft transform kernels and extract necessary (right side of the spectrum) parts
	for (int k = 0; k < B; k++) 
	{
		mFft.forward(mKernelStorageTime[k], mFftTmpStorage);
		mFftInv.forward(mKernelStorageTimeInv[k], mFftTmpStorageInv);
		// extract real part
		for (int n = 0; n < Fft_Domain_Size; n++) 
		{
			mKernelStorage[k][n] = mFftTmpStorage[n];
			mKernelStorageInv[k][n] = std::conj(mFftTmpStorageInv[n]);
		}
	}
	// mark relevant kernel values
	for (int k = 0; k < B; k++) 
	{
		mKernelMask[k].clear();
		mKernelMaskInv[k].clear();
		for (int n = 0; n < Fft_Domain_Size; n++) 
		{
			const double kernelAbs = std::abs(mKernelStorage[k][n]);
			if (kernelAbs > KernelThreshold)
			{ 
				mKernelMask[k].push_back(n);
			}
			const double kernelAbsInv = std::abs(mKernelStorageInv[k][n]);
			if (kernelAbsInv > KernelThreshold)
			{
				mKernelMaskInv[k].push_back(n);
			}
		}
	}
	// pass kernels to handlers
	for (int octave = 0; octave < OctaveNumber; octave++)
	{
		mTransformationHandlers[octave].initKernels(mKernelStorage, mKernelStorageInv, mKernelMask, mKernelMaskInv);
	}
};

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::setConcertPitch(double concertPitch)
{
	mConcertPitch = concertPitch;
	initKernelFreqs();
	recalculateKernels();
};

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::inputBlock(double* const data, const int blockSize)
{
	// check for new kernels
	if (mNewKernels.load())
	{
		mNewKernels.store(false);
		// calc the windows and give them to handlers
		calculateKernels();
	}
	// process Filterbank and create Schedule
	mFilterbank.inputBlock(data, blockSize);
	// determine cqt positions and schedule them
	mCqtSchedule.clear();
	for (int i = 0; i < blockSize; i++)
	{
		for (int octave = (OctaveNumber - 1); octave >= 0; octave--) // starting with lowest pitched octave for historical reasons
		{
			mSampleCounters[octave]++;
			if (mSampleCounters[octave] >= mLatencySamples[octave])
			{
				mSampleCounters[octave] = 0;
				const int delayOctaveRate = static_cast<int>(static_cast<double>(blockSize - i - 1) * mSampleRatesByOriginRate[octave]);
				mCqtSchedule.push_back({ i, octave, delayOctaveRate });
			}
		}
	}
};

template <int B, int OctaveNumber>
inline double* ConstantQTransform<B, OctaveNumber>::outputBlock(const int blockSize)
{
	return mFilterbank.outputBlock(blockSize);
};

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::cqt(const ScheduleElement schedule)
{
	mTransformationHandlers[schedule.octave].cqtTransform(schedule);
};

template <int B, int OctaveNumber>
inline void ConstantQTransform<B, OctaveNumber>::icqt(const ScheduleElement schedule)
{
	mTransformationHandlers[schedule.octave].icqtTransform(schedule);
};


};