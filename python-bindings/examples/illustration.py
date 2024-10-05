#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 20:32:31 2023

@author: jmerkt
"""

import sys
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import chirp

sys.path.insert(1, './../build')
import prtcqt as cqt


bins_per_octave: int = 24
number_octaves: int = 9
number_bins: int = bins_per_octave * number_octaves
block_size: int = 64
sample_rate: float = 48000.
number_blocks: int = 5000
cqt_hop_size = 64


def synth_sine_mixture(time: np.ndarray, bins: np.ndarray, tones: List[int]) -> np.ndarray:
    data = np.zeros(time.shape[-1], dtype=np.float64)
    for octave in range(bins.shape[0]):
        for tone in range(bins.shape[1]):
            if tone in tones:
                f0 = bins[octave, tone]
                data += np.sin(2.0 * np.pi * time_data * f0, dtype=np.float64)
    data /= float(len(tones) * bins.shape[0])
    return data


def synth_chirp(time: np.ndarray):
    f0 = 20.0
    f1 = 17000.0
    data = chirp(time_data, f0=f0, t1=time[-1], f1=f1, method='logarithmic')
    return data


def process_block_cqt(cqt_instance: cqt.Cqt24, data: np.ndarray):
    cqt_instance.inputBlock(data)
    schedule = cqt_instance.getCqtSchedule()
    for s in schedule:
        cqt_instance.cqt(s)
        cqt_instance.icqt(s)
    magnitudes = np.zeros(number_bins)
    for i_octave in range(number_octaves):
        octave_data = cqt_instance.getOctaveCqtBuffer(i_octave)
        magnitudes[i_octave * bins_per_octave: (i_octave + 1) * bins_per_octave] = np.flip(np.abs(octave_data))
    output_block = cqt_instance.outputBlock(data.shape[0])
    return magnitudes, output_block


def process_block_scqt(cqt_instance: cqt.SlidingCqt24, data: np.ndarray):
    cqt_instance.inputBlock(data, data.shape[0])
    magnitudes = np.zeros(number_bins)
    for i_octave in range(number_octaves):
        octave_data = cqt_instance_sliding.getOctaveValues(i_octave)
        magnitudes[i_octave * bins_per_octave: (i_octave + 1) * bins_per_octave] = np.flip(np.abs(octave_data))
    output_block = cqt_instance_sliding.outputBlock(block_size)
    return magnitudes, output_block


def process_signal(signal: np.ndarray, cqt_instance: cqt.Cqt24, sliding_cqt_instance: cqt.SlidingCqt24):
    cqt_magnitudes = np.zeros((number_blocks, number_bins), np.float32)
    cqt_magnitudes_sliding = np.zeros((number_blocks, number_bins), np.float32)

    audio_output = np.zeros(signal.shape[0], np.float64)
    audio_output_sliding = np.zeros(signal.shape[0], np.float64)

    for i_block in range(number_blocks):
        input_data = signal[i_block * block_size: (i_block + 1) * block_size]
        # Cqt
        magnitudes, output_data = process_block_cqt(cqt_instance, input_data)
        cqt_magnitudes[i_block, :] = magnitudes
        audio_output[i_block * block_size: (i_block + 1) * block_size] = output_data
        # SlidingCqt
        magnitudes, output_data = process_block_scqt(cqt_instance_sliding, input_data)
        cqt_magnitudes_sliding[i_block, :] = magnitudes
        audio_output_sliding[i_block * block_size: (i_block + 1) * block_size] = output_data

    # Empty buffers
    for i_block in range(number_blocks):
        foo1, foo2 = process_block_cqt(cqt_instance, np.zeros(block_size))
        foo1, foo2 = process_block_scqt(cqt_instance_sliding, np.zeros(block_size))

    return cqt_magnitudes, cqt_magnitudes_sliding, audio_output, audio_output_sliding


def plot_signal_on_axis(ax, time_data, input_data, output_data, cqt_magnitudes):
    cqt_magnitudes = cqt_magnitudes.transpose()

    input_data_f = np.fft.rfft(input_data)
    output_data_f = np.fft.rfft(output_data)
    freq = np.fft.rfftfreq(input_data.shape[-1]) * sample_rate

    bins = np.arange(0, number_bins)
    bin_magnitude_means = np.mean(cqt_magnitudes, axis=1)

    ax[0].plot(time_data, input_data, label='input')
    ax[0].plot(time_data, output_data, label='output')
    ax[0].set_xlabel('time in s')
    ax[0].set_ylabel('amplitude')
    ax[0].legend()

    ax[1].plot(freq, np.abs(input_data_f), label='input')
    ax[1].plot(freq, np.abs(output_data_f), label='output')
    ax[1].set_xlabel('frequency in Hz')
    ax[1].set_ylabel('magnitude')
    ax[1].set_xscale('log')
    ax[1].set_xlim(10, sample_rate / 2)
    ax[1].legend()

    ax[2].imshow(cqt_magnitudes, interpolation='bilinear')
    ax[2].set_aspect('auto')
    ax[2].set_title('Cqt magnitudes')
    ax[2].set_xlabel('i_block')
    ax[2].set_ylabel('i_bin')

    ax[3].plot(bin_magnitude_means)
    ax[3].set_xlabel('i_bin')
    ax[3].set_ylabel('mean magnitude')
    ax[3].invert_xaxis()


cqt_instance = cqt.Cqt24()
cqt_instance.init(cqt_hop_size)
cqt_instance.initFs(sample_rate, block_size)

cqt_instance_sliding = cqt.SlidingCqt24()
cqt_instance_sliding.init(sample_rate, block_size)

# Get bin frequencies and print
bin_freqs = np.zeros((number_octaves, bins_per_octave), np.float64)
for i_octave in range(number_octaves):
    bin_freqs[i_octave, :] = np.flip(cqt_instance_sliding.getOctaveBinFreqs(i_octave))
# print(bin_freqs)

# Synthesize sine tones
time_data = np.linspace(0.0, block_size / sample_rate * number_blocks, block_size * number_blocks)
sine_mixture_signal = synth_sine_mixture(time_data, bin_freqs, [0, 5])
chirp_signal = synth_chirp(time_data)

cqt_magnitudes, cqt_magnitudes_sliding, audio_output, audio_output_sliding \
    = process_signal(sine_mixture_signal, cqt_instance, cqt_instance_sliding)

fig_cqt, ax_cqt = plt.subplots(4, 2)
fig_cqt.suptitle("CQT")
fig_sliding_cqt, ax_sliding_cqt = plt.subplots(4, 2)
fig_sliding_cqt.suptitle("Sliding CQT")

plot_signal_on_axis(ax_cqt[:, 0], time_data, sine_mixture_signal, audio_output, cqt_magnitudes)
plot_signal_on_axis(ax_sliding_cqt[:, 0], time_data, sine_mixture_signal, audio_output_sliding, cqt_magnitudes_sliding)

cqt_magnitudes, cqt_magnitudes_sliding, audio_output, audio_output_sliding \
    = process_signal(chirp_signal, cqt_instance, cqt_instance_sliding)

plot_signal_on_axis(ax_cqt[:, 1], time_data, chirp_signal, audio_output, cqt_magnitudes)
plot_signal_on_axis(ax_sliding_cqt[:, 1], time_data, chirp_signal, audio_output_sliding, cqt_magnitudes_sliding)
plt.show()





