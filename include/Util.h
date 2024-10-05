#pragma once

#include <cmath>

namespace Cqt
{
    /** 
     * Computes the frequency of the first bin of the highest octave for a given concert pitch.
     * https://en.wikipedia.org/wiki/Piano_key_frequencies
     * For a concert pitch of 440 Hz this results in 8.37 kHz.
    */
    inline double computeReferenceFrequency(const double concert_pitch)
    {
        return std::pow(2., ((100. - 49.) / 12.)) * concert_pitch;
    }

    inline double computeBinFrequency(const double reference_frequency, const int bins_per_octave, const int octave, const int tone)
    {
        return (reference_frequency / std::pow(2., (octave + 1))) * std::pow(2., static_cast<double>(bins_per_octave + tone) / static_cast<double>(bins_per_octave));
    }
}