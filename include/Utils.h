/*
 ==============================================================================

 This file is part of the rt-cqt library. Copyright (C) the rt-cqt developers.

 See LICENSE.txt for  more info.

 ==============================================================================
*/

#pragma once

namespace Cqt
{

template<typename FloatType>
static inline FloatType Pi()
{
	return static_cast<FloatType>(3.14159265358979323846);
}

template<typename FloatType>
static inline FloatType TwoPi()
{
	return static_cast<FloatType>(6.283185307179586476925286766);
}

template<typename T>
static inline T Clip(const T input, const T lowerBound, const T upperBound)
{
	T y = input > upperBound ? upperBound : input;
	y = y < lowerBound ? lowerBound : y;
	return y;
}

template<typename T>
static inline T ClipDelete(const T input, const T lowerBound, const T upperBound)
{
	T y = input > upperBound ? static_cast<T>(-1) : input;
	y = y < lowerBound ? static_cast<T>(-1) : y;
	return y;
}

}



