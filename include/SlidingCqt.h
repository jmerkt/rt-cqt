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
- Windowing
- Parallelization
- Block processing instead of sample
*/


namespace Cqt 
{

using namespace std::complex_literals;

template <int B, int OctaveNumber, bool Windowing = true>
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
    void calculateKernels();

    double mFs{ 48000. };
    double mSampleRates[OctaveNumber];
    int mBlockSizes[OctaveNumber];
    double mConcertPitch{ 440. };

    std::atomic<bool> mNewKernels{ true };

    ResamplingFilterbank<OctaveNumber> mFilterbank;

    // delay lines / circular buffers for cqt stuff - do they need to be interpolating to match Nk more precicely?
    CircularBuffer<double> mDelayLines[OctaveNumber][B];

    // precalculated exp stuff
    double mQ[3]{0., 0., 0.};
    std::complex<double> mExpQ[3]{0.+0.i, 0.+0.i, 0.+0.i};
    std::complex<double> mExpQNk[OctaveNumber][B][3];
    std::complex<double> mFtPrev[OctaveNumber][B][3];

    int mNk[OctaveNumber][B];
    double mNkDouble[OctaveNumber][B];
    double mOneDivNkDouble[OctaveNumber][B];
    double mBinFreqs[OctaveNumber][B];
    

    size_t mSamplesToProcess[OctaveNumber];

    // cqt data
    CircularBuffer<std::complex<double>> mCqtData[OctaveNumber][B];

    // windowing
    static constexpr double mWindowCoeffs[3] = {0.5, -0.25, -0.25};
    static constexpr double mQAdd[3] = {0., -1., 1};
};

template <int B, int OctaveNumber, bool Windowing>
SlidingCqt<B, OctaveNumber, Windowing>::SlidingCqt()
{
    for (int o = 0; o < OctaveNumber; o++) 
	{
        mSampleRates[o] = 48000.;
        mBlockSizes[o] = 0;
        mSamplesToProcess[o] = 0;
        for(int tone = 0; tone < B; tone++)
        {
            mBinFreqs[o][tone] = 0.; 
            mNkDouble[o][tone] = 0.;
            mOneDivNkDouble[o][tone] = 0.;

            for(int w = 0; w < 3; w++)
            {
                mExpQNk[o][tone][w] = 0. + 0i;
                mFtPrev[o][tone][w] = 0. + 0.i;
            }    
        }
    }
}

template <int B, int OctaveNumber, bool Windowing>
inline void SlidingCqt<B, OctaveNumber, Windowing>::init(const double samplerate, const int blockSize)
{
	mFilterbank.init(samplerate, blockSize, blockSize * 2);
	mFs = mFilterbank.getOriginSamplerate();
    const int originBlockSize = mFilterbank.getOriginBlockSize();
	for (int octave = 0; octave < OctaveNumber; octave++) 
	{
		mSampleRates[octave] = mFs / std::pow(2., octave);
        mBlockSizes[octave] = (originBlockSize / std::pow(2, octave) >= 1) ? originBlockSize / std::pow(2, octave) : 1;
	}
	calculateKernels();
    // initialize delay lines
    for(int o = 0; o < OctaveNumber; o++)
    {
        for(int tone = 0; tone < B; tone++)
        {
            mDelayLines[o][tone].changeSize(mNk[o][tone]);
        }
    }

    for (int o = 0; o < OctaveNumber; o++) 
	{
        for(int tone = 0; tone < B; tone++)
        {
            const size_t octaveBlockSize = blockSize / std::pow(2, o);
            const size_t octaveBlockSizeClipped = octaveBlockSize > 2 ? octaveBlockSize : 2; 
            mCqtData[o][tone].changeSize(octaveBlockSizeClipped);

            for(int w = 0; w < 3; w++)
            {
                mFtPrev[o][tone][w] = 0. + 0.i;
            }    
        }
    }
};

template <int B, int OctaveNumber, bool Windowing>
inline void SlidingCqt<B, OctaveNumber, Windowing>::setConcertPitch(double concertPitch)
{
	mConcertPitch = concertPitch;
	recalculateKernels();
};

template <int B, int OctaveNumber, bool Windowing>
inline void SlidingCqt<B, OctaveNumber, Windowing>::inputBlock(double* const data, const int blockSize)
{
    // check for new kernels
	if (mNewKernels.load())
	{
		mNewKernels.store(false);
		// calc the windows and give them to handlers
		calculateKernels();
	}

    // push data into multirate resampling
    mFilterbank.inputBlock(data, blockSize);
    // process all cqt sample based on numbers pushed into stage buffers
    for(int o = 0; o < OctaveNumber; o++)
    {
        const BufferPtr inputBuffer = mFilterbank.getStageInputBuffer(o);
        const size_t nSamples = inputBuffer->getWriteReadDistance();
        mSamplesToProcess[o] = nSamples;

        for(size_t i = 0; i < nSamples; i++)
        {
            const double inputSample = inputBuffer->pullSample();
            for (int tone = 0; tone < B; tone++) 
            {
                const double oneDivNkDouble = mOneDivNkDouble[o][tone];
                const int nk = mNk[o][tone];
                const double delaySample = mDelayLines[o][tone].pullDelaySample(nk - 1);

                const std::complex<double> delayCplx{delaySample, 0.};

                if constexpr (Windowing == false)
                {
                    const std::complex<double> expQ = mExpQ[0];
                    const std::complex<double> expQNk = mExpQNk[o][tone][0];
                    const std::complex<double> FtPrev = mFtPrev[o][tone][0];
                    const std::complex<double> expMInput = expQ * inputSample;

                    const std::complex<double> Ft = expQNk * (FtPrev + (expMInput - delayCplx) * oneDivNkDouble);

                    mFtPrev[o][tone][0] = Ft;

                    mCqtData[o][tone].pushSample(Ft); 
                }
                else
                {
                    std::complex<double> FtSum = 0. + 0.i;
                    for(int w = 0; w < 3; w++)
                    {
                        const std::complex<double> expQ = mExpQ[w];
                        const std::complex<double> expQNk = mExpQNk[o][tone][w];
                        const std::complex<double> FtPrev = mFtPrev[o][tone][w];
                        const std::complex<double> expMInput = expQ * inputSample;

                        const std::complex<double> Ft = expQNk * (FtPrev + (expMInput - delayCplx) * oneDivNkDouble);

                        mFtPrev[o][tone][w] = Ft;

                        FtSum += mWindowCoeffs[w] * Ft;
                    }
                    mCqtData[o][tone].pushSample(FtSum);
                }
             
                mDelayLines[o][tone].pushSample(inputSample);

                //const std::complex<double> inputCplx{inputSample, 0.};
                //const std::complex<double> expMDelay = expQ * delaySample;
                //const std::complex<double> Ft = expQNk * (FtPrev + (expMDelay - inputCplx) * oneDivNkDouble);
            }
        }
    }
};

template <int B, int OctaveNumber, bool Windowing>
inline double* SlidingCqt<B, OctaveNumber, Windowing>::outputBlock(const int blockSize)
{
    for(int o = 0; o < OctaveNumber; o++)
    {
        const BufferPtr outputBuffer = mFilterbank.getStageOutputBuffer(o);
        const size_t nSamples = mSamplesToProcess[o];
        for(size_t i = 0; i < nSamples; i++)
        {
            double outputSample = 0.;
            for (int tone = 0; tone < B; tone++) 
            {
                const std::complex<double> expQNk = mExpQNk[o][tone][0];
                const std::complex<double> Ft = mCqtData[o][tone].pullSample();
                outputSample += (Ft * expQNk).real();
            }
            outputBuffer->pushSample(outputSample);
        }
    }
	return mFilterbank.outputBlock(blockSize);
};

template <int B, int OctaveNumber, bool Windowing>
inline void SlidingCqt<B, OctaveNumber, Windowing>::calculateKernels()
{
    
    for(int w = 0; w < 3; w++)
    {
        // Q
        mQ[w] = 1. / (std::pow(2., 1. / static_cast<double>(B)) - 1.);
        mQ[w] += mQAdd[w];

        // exp multiplication 
        mExpQ[w] = std::exp(-1i * TwoPi<double>() * mQ[w]);
    }

    /*
	https://en.wikipedia.org/wiki/Piano_key_frequencies
	f(n) = 2^((n-49)/12) * mConcertPitch
	Range in highest octave from n = [100, 111] ~ [8.37 kHz, 15.804 kHz]
	*/
    const double fRef = std::pow(2., ((100. - 49.) / 12.)) * mConcertPitch;

    for(int o = 0; o < OctaveNumber; o++)
    {
        // fs
        const double fs = mSampleRates[o];
        for (int tone = 0; tone < B; tone++) 
		{
            // fk
            const double fk = (fRef / std::pow(2., (o + 1))) * std::pow(2., static_cast<double>(B + tone) / static_cast<double>(B));
            mBinFreqs[o][tone] = fk;
            // Nk
            mNk[o][tone] = static_cast<int>((fs / fk) * mQ[0]);
            const double nkDouble = static_cast<double>(mNk[o][tone]);
            mNkDouble[o][tone] = nkDouble;
            mOneDivNkDouble[o][tone] = 1. / nkDouble;

            for(int w = 0; w < 3; w++)
            {
                mExpQNk[o][tone][w] = std::exp(1i * TwoPi<double>() * mQ[w] * mOneDivNkDouble[o][tone]);
            }
        }      
    }
};


}