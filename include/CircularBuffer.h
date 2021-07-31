/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once
#include <vector>
#include <cmath>

template <typename T>
class CircularBuffer
{
public:
	CircularBuffer(size_t bufferSize = 128) { changeSize(bufferSize); };
	~CircularBuffer() = default;
	void changeSize(size_t bufferSize);
	size_t getBufferSize() { return mBufferSize; };
	void pushSample(const T value);
	void pushBlock(const T* const data, const int blockSize);
	T pullSample();
	void pullBlock(T* const data, const int blockSize);
	void pullBlockAdd(T* const data, const int blockSize);
	T pullDelaySample(const int delay);
	void pullDelayBlock(T* const data, const int delay, const int blockSize);
	void modulateDelayBlock(const T* const data, const int delay, const int blockSize); 
	inline size_t nextPowOfTwo(size_t size) { return std::pow(2, std::ceil(std::log(size) / std::log(2))); };
	size_t getWriteReadDistance();
	inline void resetReadPointer() { mReadPointer = mWritePointer; };
protected:
	std::vector<T> mBuffer;
	size_t mBufferSize{ 0 };
	size_t mBufferSizeMinOne{ 0 };
	size_t mWritePointer{ 0 };
	size_t mReadPointer{ 0 };
};

template <typename T>
inline void CircularBuffer<T>::changeSize(size_t bufferSize)
{
	mBufferSize = nextPowOfTwo(bufferSize);
	mBufferSizeMinOne = mBufferSize - 1;
	mBuffer.resize(mBufferSize, static_cast<T>(0.));
	mWritePointer = 0;
	mReadPointer = 0;
}

template <typename T>
inline void CircularBuffer<T>::pushSample(const T value) 
{
	mWritePointer += 1;
	mWritePointer = mWritePointer & mBufferSizeMinOne;
	mBuffer[mWritePointer] = value;
}

template <typename T>
inline void CircularBuffer<T>::pushBlock(const T* const data, const int blockSize)
{
	for (int i = 0; i < blockSize; i++)
	{
		mWritePointer += 1;
		mWritePointer = mWritePointer & mBufferSizeMinOne;
		mBuffer[mWritePointer] = data[i];
	}
}

template <typename T>
inline T CircularBuffer<T>::pullDelaySample(const int delay)
{
	int position = (static_cast<int>(mWritePointer) - delay);
	if (position < 0)
	{
		position += mBufferSize;
	}
	return mBuffer[position];
}

template <typename T>
inline void CircularBuffer<T>::pullDelayBlock(T* const data, const int delay, const int blockSize)
{
	int position = (static_cast<int>(mWritePointer) - delay);
	if (position < 0)
	{
		position += mBufferSize;
	}
	size_t readPointer = static_cast<size_t>(position);
	for (int i = 0; i < blockSize; i++)
	{
		data[i] = mBuffer[readPointer];
		readPointer += 1;
		readPointer = readPointer & mBufferSizeMinOne;
	}
}

template <typename T>
inline void CircularBuffer<T>::modulateDelayBlock(const T* const data, const int delay, const int blockSize)
{
	int position = (static_cast<int>(mWritePointer) - delay);
	if (position < 0)
	{
		position += mBufferSize;
	}
	size_t readPointer = static_cast<size_t>(position);
	for (int i = 0; i < blockSize; i++)
	{
		mBuffer[readPointer] *= data[i];
		readPointer += 1;
		readPointer = readPointer & mBufferSizeMinOne;
	}
}

template <typename T>
inline T CircularBuffer<T>::pullSample()
{
	mReadPointer += 1;
	mReadPointer = mReadPointer & mBufferSizeMinOne;
	return mBuffer[mReadPointer];
}

template <typename T>
inline void CircularBuffer<T>::pullBlock(T* const data, const int blockSize)
{
	for (int i = 0; i < blockSize; i++)
	{
		mReadPointer += 1;
		mReadPointer = mReadPointer & mBufferSizeMinOne;
		data[i] = mBuffer[mReadPointer];
	}
}

template <typename T>
inline void CircularBuffer<T>::pullBlockAdd(T* const data, const int blockSize)
{
	for (int i = 0; i < blockSize; i++)
	{
		mReadPointer += 1;
		mReadPointer = mReadPointer & mBufferSizeMinOne;
		data[i] += mBuffer[mReadPointer];
	}
}

template <typename T>
inline size_t CircularBuffer<T>::getWriteReadDistance()
{
	return  mWritePointer >= mReadPointer ? mWritePointer - mReadPointer : mBufferSize - mReadPointer + mWritePointer;
}

