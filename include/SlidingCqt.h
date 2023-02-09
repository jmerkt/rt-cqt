/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

#include "ResamplingFilterbank.h"
#include "Utils.h"
#include <complex>

/*
TODO: 
- Optimization:
    * Parallelization / Vectorization
*/


namespace Cqt 
{

using namespace std::complex_literals;

template <size_t B, size_t OctaveNumber, bool Windowing = true>
class SlidingCqt
{
public:
    SlidingCqt();
    ~SlidingCqt() = default;

    void init(const double samplerate, const int blockSize);

    void setConcertPitch(double concertPitch);
	inline void recalculateKernels() { mNewKernels.store(true); };

    void inputBlock(double* const data, const int blockSize);
	double* outputBlock(const int blockSize);

    inline CircularBuffer<std::complex<double>>* getOctaveCqtBuffer(const int octave) { return &mCqtData[octave][0]; };
    inline size_t getSamplesToProcess(const int octave){return mSamplesToProcess[octave];};

    inline double getOctaveSampleRate(const int octave) { return mSampleRates[octave]; };
    inline int getOctaveBlockSize(const int octave) { return mBlockSizes[octave]; };
    
    inline double* getOctaveBinFreqs(const int octave){ return &mBinFreqs[octave][0]; };

private:
    void computeKernels();

    double mFs{ 48000. };
    double mSampleRates[OctaveNumber];
    int mBlockSizes[OctaveNumber];
    double mConcertPitch{ 440. };

    std::atomic<bool> mNewKernels{ true };

    ResamplingFilterbank<OctaveNumber> mFilterbank;

    // Do we need interpolating lines to match Nk more precicely?
    CircularBuffer<double> mDelayLines[OctaveNumber];

    // Pre-calculated exp stuff
    double mQ[3]{0., 0., 0.};
    std::complex<double> mExpQ[3]{0.+0.i, 0.+0.i, 0.+0.i};
    std::complex<double> mExpQNk[OctaveNumber][B][3];
    std::complex<double> mFtPrev[OctaveNumber][B][3];

    int mNk[OctaveNumber][B];
    double mNkDouble[OctaveNumber][B];
    double mOneDivNkDouble[OctaveNumber][B];
    double mBinFreqs[OctaveNumber][B];
    

    size_t mSamplesToProcess[OctaveNumber];

    // Cqt data
    CircularBuffer<std::complex<double>> mCqtData[OctaveNumber][B];

    // Windowing
    static constexpr double mWindowCoeffs[3] = {0.5, -0.25, -0.25};
    static constexpr double mQAdd[3] = {0., -1., 1};

    // Buffers for block processing
    std::vector<double> mInputSamplesBuffer[OctaveNumber];
    std::vector<double> mInputDelaySamplesBuffer[OctaveNumber][B];
    std::vector<std::complex<double>> mInputFtBuffer[OctaveNumber][B];

    std::vector<std::complex<double>> mOutputFtBuffer[OctaveNumber][B];
    std::vector<double> mOutputSamplesTonesBuffer[OctaveNumber][B];
    std::vector<double> mOutputSamplesBuffer[OctaveNumber];
};

template <size_t B, size_t OctaveNumber, bool Windowing>
SlidingCqt<B, OctaveNumber, Windowing>::SlidingCqt()
{
    for (size_t i_octave = 0; i_octave < OctaveNumber; i_octave++) 
	{
        mSampleRates[i_octave] = 48000.;
        mBlockSizes[i_octave] = 0;
        mSamplesToProcess[i_octave] = 0;
        for(size_t i_tone = 0; i_tone < B; i_tone++)
        {
            mBinFreqs[i_octave][i_tone] = 0.; 
            mNkDouble[i_octave][i_tone] = 0.;
            mOneDivNkDouble[i_octave][i_tone] = 0.;

            for(size_t i_window = 0; i_window < 3; i_window++)
            {
                mExpQNk[i_octave][i_tone][i_window] = 0. + 0i;
                mFtPrev[i_octave][i_tone][i_window] = 0. + 0.i;
            }    
        }
    }
}

template <size_t B, size_t OctaveNumber, bool Windowing>
inline void SlidingCqt<B, OctaveNumber, Windowing>::init(const double samplerate, const int blockSize)
{
	mFilterbank.init(samplerate, blockSize, blockSize * 2);
	mFs = mFilterbank.getOriginSamplerate();
    const int originBlockSize = mFilterbank.getOriginBlockSize();
	for (size_t i_octave = 0; i_octave < OctaveNumber; i_octave++) 
	{
		mSampleRates[i_octave] = mFs / std::pow(2., i_octave);
        mBlockSizes[i_octave] = (originBlockSize / std::pow(2, i_octave) >= 1) ? originBlockSize / std::pow(2, i_octave) : 1;
	}
	computeKernels();
    // initialize delay lines
    for(size_t i_octave = 0; i_octave < OctaveNumber; i_octave++)
    {
        int maxNk = -1;
        for(size_t i_tone = 0; i_tone < B; i_tone++)
        {
            const int delaySize = mNk[i_octave][i_tone] + mBlockSizes[i_octave] + 1u;
            if(delaySize > maxNk)
                maxNk = delaySize;
        }
        mDelayLines[i_octave].changeSize(static_cast<size_t>(maxNk));
    }

    for (size_t i_octave = 0; i_octave < OctaveNumber; i_octave++) 
	{
        for(size_t i_tone = 0; i_tone < B; i_tone++)
        {
            const size_t octaveBlockSize = blockSize / std::pow(2, i_octave);
            const size_t octaveBlockSizeClipped = octaveBlockSize > 2 ? octaveBlockSize : 2; 
            mCqtData[i_octave][i_tone].changeSize(octaveBlockSizeClipped);

            for(size_t i_window = 0; i_window < 3u; i_window++)
            {
                mFtPrev[i_octave][i_tone][i_window] = 0. + 0.i;
            }    
        }
    }

    // Buffers for block processing
    for (size_t i_octave = 0; i_octave < OctaveNumber; i_octave++) 
	{
        const size_t octaveBlockSize = mBlockSizes[i_octave];

        mInputSamplesBuffer[i_octave].resize(octaveBlockSize, 0.);
        mOutputSamplesBuffer[i_octave].resize(octaveBlockSize, 0.);
        for(size_t i_tone = 0; i_tone < B; i_tone++)
        {
            mInputDelaySamplesBuffer[i_octave][i_tone].resize(octaveBlockSize, 0.);
            mInputFtBuffer[i_octave][i_tone].resize(octaveBlockSize, {0., 0.});

            mOutputFtBuffer[i_octave][i_tone].resize(octaveBlockSize, {0., 0.});  
            mOutputSamplesTonesBuffer[i_octave][i_tone].resize(octaveBlockSize, 0.);
        }
    }
};

template <size_t B, size_t OctaveNumber, bool Windowing>
inline void SlidingCqt<B, OctaveNumber, Windowing>::setConcertPitch(double concertPitch)
{
	mConcertPitch = concertPitch;
	recalculateKernels();
};

template <size_t B, size_t OctaveNumber, bool Windowing>
inline void SlidingCqt<B, OctaveNumber, Windowing>::inputBlock(double* const data, const int blockSize) 
{
    // check for new kernels
	if (mNewKernels.load())
	{
		mNewKernels.store(false);
		// calc the windows and give them to handlers
		computeKernels();
	}

    // push data into multirate resampling
    mFilterbank.inputBlock(data, blockSize);
    // process all cqt sample based on numbers pushed into stage buffers
    for(size_t i_octave = 0; i_octave < OctaveNumber; i_octave++)
    {
        const BufferPtr inputBuffer = mFilterbank.getStageInputBuffer(i_octave);
        const size_t nOctaveSamples = inputBuffer->getWriteReadDistance();
        mSamplesToProcess[i_octave] = nOctaveSamples;

        inputBuffer->pullBlock(mInputSamplesBuffer[i_octave].data(), nOctaveSamples);
        mDelayLines[i_octave].pushBlock(mInputSamplesBuffer[i_octave].data(), nOctaveSamples);
        for (size_t i_tone = 0; i_tone < B; i_tone++) 
        {
            const int nk = mNk[i_octave][i_tone];
            mDelayLines[i_octave].pullDelayBlock(mInputDelaySamplesBuffer[i_octave][i_tone].data(), nk + nOctaveSamples - 1, nOctaveSamples);
        }
        for(size_t i_sample = 0; i_sample < nOctaveSamples; i_sample++)
        {
            // #pragma omp simd
            for (size_t i_tone = 0; i_tone < B; i_tone++) 
            {
                const double oneDivNkDouble = mOneDivNkDouble[i_octave][i_tone];
                const double delaySample = mInputDelaySamplesBuffer[i_octave][i_tone][i_sample];

                const std::complex<double> delayCplx{delaySample, 0.};

                if constexpr (Windowing == false)
                {
                    const std::complex<double> expQ = mExpQ[0];
                    const std::complex<double> expQNk = mExpQNk[i_octave][i_tone][0];
                    const std::complex<double> FtPrev = mFtPrev[i_octave][i_tone][0];
                    const std::complex<double> expMInput = expQ * mInputSamplesBuffer[i_octave][i_sample];

                    const std::complex<double> Ft = expQNk * (FtPrev + (expMInput - delayCplx) * oneDivNkDouble);

                    mFtPrev[i_octave][i_tone][0] = Ft;
                    mInputFtBuffer[i_octave][i_tone][i_sample] = Ft;
                }
                else
                {
                    std::complex<double> FtSum = 0. + 0.i;
                    for(size_t i_window = 0; i_window < 3u; i_window++)
                    {
                        const std::complex<double> expQ = mExpQ[i_window];
                        const std::complex<double> expQNk = mExpQNk[i_octave][i_tone][i_window];
                        const std::complex<double> FtPrev = mFtPrev[i_octave][i_tone][i_window];
                        const std::complex<double> expMInput = expQ * mInputSamplesBuffer[i_octave][i_sample];

                        const std::complex<double> Ft = expQNk * (FtPrev + (expMInput - delayCplx) * oneDivNkDouble);

                        mFtPrev[i_octave][i_tone][i_window] = Ft;

                        FtSum += mWindowCoeffs[i_window] * Ft;
                    }
                    mInputFtBuffer[i_octave][i_tone][i_sample] = FtSum;
                }
            }
        }
        for (size_t i_tone = 0; i_tone < B; i_tone++) 
        {
            mCqtData[i_octave][i_tone].pushBlock(mInputFtBuffer[i_octave][i_tone].data(), nOctaveSamples);
        }
    }
};

template <size_t B, size_t OctaveNumber, bool Windowing>
inline double* SlidingCqt<B, OctaveNumber, Windowing>::outputBlock(const int blockSize)
{
    for(size_t i_octave = 0; i_octave < OctaveNumber; i_octave++)
    {
        const BufferPtr outputBuffer = mFilterbank.getStageOutputBuffer(i_octave);
        const size_t nOctaveSamples = mSamplesToProcess[i_octave];

        for (size_t i_tone = 0; i_tone < B; i_tone++) 
        {
            mCqtData[i_octave][i_tone].pullBlock(mOutputFtBuffer[i_octave][i_tone].data(), nOctaveSamples);
        }
        for(size_t i_sample = 0; i_sample < nOctaveSamples; i_sample++)
        {
            // #pragma omp simd
            for (size_t i_tone = 0; i_tone < B; i_tone++) 
            {
                const std::complex<double> expQNk = mExpQNk[i_octave][i_tone][0];
                const std::complex<double> Ft = mOutputFtBuffer[i_octave][i_tone][i_sample];
                mOutputSamplesTonesBuffer[i_octave][i_tone][i_sample] = (Ft * expQNk).real();
            }
            mOutputSamplesBuffer[i_octave][i_sample] = 0.;
            for (size_t i_tone = 0; i_tone < B; i_tone++) 
            {
                mOutputSamplesBuffer[i_octave][i_sample] += mOutputSamplesTonesBuffer[i_octave][i_tone][i_sample];
            }
        }
        outputBuffer->pushBlock(mOutputSamplesBuffer[i_octave].data(), nOctaveSamples);
    }
	return mFilterbank.outputBlock(blockSize);
};

template <size_t B, size_t OctaveNumber, bool Windowing>
inline void SlidingCqt<B, OctaveNumber, Windowing>::computeKernels()
{
    
    for(size_t i_window = 0; i_window < 3u; i_window++)
    {
        // Q
        mQ[i_window] = 1. / (std::pow(2., 1. / static_cast<double>(B)) - 1.);
        mQ[i_window] += mQAdd[i_window];

        // exp multiplication 
        mExpQ[i_window] = std::exp(-1i * TwoPi<double>() * mQ[i_window]);
    }

    /*
	https://en.wikipedia.org/wiki/Piano_key_frequencies
	f(n) = 2^((n-49)/12) * mConcertPitch
	Range in highest octave from n = [100, 111] ~ [8.37 kHz, 15.804 kHz]
	*/
    const double fRef = std::pow(2., ((100. - 49.) / 12.)) * mConcertPitch;

    for(size_t i_octave = 0; i_octave < OctaveNumber; i_octave++)
    {
        // fs
        const double fs = mSampleRates[i_octave];
        for (size_t i_tone = 0; i_tone < B; i_tone++) 
		{
            // fk
            const double fk = (fRef / std::pow(2., (i_octave + 1))) * std::pow(2., static_cast<double>(B + i_tone) / static_cast<double>(B));
            mBinFreqs[i_octave][i_tone] = fk;
            // Nk
            mNk[i_octave][i_tone] = static_cast<int>((fs / fk) * mQ[0]);
            const double nkDouble = static_cast<double>(mNk[i_octave][i_tone]);
            mNkDouble[i_octave][i_tone] = nkDouble;
            mOneDivNkDouble[i_octave][i_tone] = 1. / nkDouble;

            for(size_t i_window = 0; i_window < 3u; i_window++)
            {
                mExpQNk[i_octave][i_tone][i_window] = std::exp(1i * TwoPi<double>() * mQ[i_window] * mOneDivNkDouble[i_octave][i_tone]);
            }
        }      
    }
};


}