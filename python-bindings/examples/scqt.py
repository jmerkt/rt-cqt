#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 20:32:31 2023

@author: jmerkt
"""

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import chirp

sys.path.insert(1, './../build')
import rtcqt as cqt

# matplotlib.use('Qt5Agg')

bins_per_octave: int = 24
n_octaves: int = 9
n_bins: int = bins_per_octave * n_octaves
block_size: int = 64
sample_rate: float = 48000.
n_blocks: int = 5000

cqt_instance = cqt.SlidingCqt24()
cqt_instance.init(sample_rate, block_size)

# Get bin frequencies and print
bin_freqs = np.zeros((n_octaves, bins_per_octave), np.float64)
for i_octave in range(n_octaves):
    bin_freqs[i_octave, :] = np.flip(cqt_instance.getOctaveBinFreqs(i_octave))
# print(bin_freqs)

# Synthesize sine tones
time_data = np.linspace(0.0, block_size / sample_rate * n_blocks, block_size * n_blocks)

audio_input_data = np.zeros(time_data.shape[-1], dtype=np.float64)
for octave in range(n_octaves):
    for bin in range(bins_per_octave):
        f0 = bin_freqs[octave, bin]
        scaling = 1. / float(n_bins)
        audio_input_data += np.sin(2.0 * np.pi * time_data * f0, dtype=np.float64) * scaling

# Evaluate using cqt
cqt_magnitudes = np.zeros((n_blocks, n_bins), np.float32)
audio_output_data = np.zeros(time_data.shape[0], np.float64)
for i_block in range(n_blocks):
    cqt_instance.inputBlock(audio_input_data[i_block * block_size: (i_block + 1) * block_size], block_size)
    for i_octave in range(n_octaves):
        octave_data = cqt_instance.getOctaveValues(i_octave) 
        cqt_magnitudes[i_block, i_octave * bins_per_octave: (i_octave + 1) * bins_per_octave] = np.flip(np.abs(octave_data))
    output_block = cqt_instance.outputBlock(block_size)
    audio_output_data[i_block * block_size: (i_block + 1) * block_size] = output_block

# Empty buffers
for i_block in range(n_blocks):
    cqt_instance.inputBlock(np.zeros(block_size), block_size)
    output_block = cqt_instance.outputBlock(block_size)
    
# Synthesize sine sweep
f0 = 20.0
f1 = 17000.0
audio_input_data_chirp = chirp(time_data, f0=f0, t1=time_data[-1], f1=f1, method='logarithmic')

cqt_magnitudes_chirp = np.zeros((n_blocks, n_bins), np.float32)
audio_output_data_chirp = np.zeros(time_data.shape[0], np.float64)
for i_block in range(n_blocks):
    cqt_instance.inputBlock(audio_input_data_chirp[i_block * block_size: (i_block + 1) * block_size], block_size)
    for i_octave in range(n_octaves):
        octave_data = cqt_instance.getOctaveValues(i_octave) 
        cqt_magnitudes_chirp[i_block, i_octave * bins_per_octave: (i_octave + 1) * bins_per_octave] = np.flip(np.abs(octave_data))
    output_block = cqt_instance.outputBlock(block_size)
    audio_output_data_chirp[i_block * block_size: (i_block + 1) * block_size] = output_block

    
# Plot result
cqt_magnitudes = cqt_magnitudes.transpose()
cqt_magnitudes_chirp = cqt_magnitudes_chirp.transpose()

# fft of input and output signals
audio_input_f = np.fft.rfft(audio_input_data)
audio_input_chirp_f = np.fft.rfft(audio_input_data_chirp)
audio_output_f = np.fft.rfft(audio_output_data)
audio_output_chirp_f = np.fft.rfft(audio_output_data_chirp)
freq = np.fft.rfftfreq(audio_input_data.shape[-1]) * sample_rate

# averaged bin magnitudes
bins = np.arange(0, n_bins)
bin_magnitude_means = np.mean(cqt_magnitudes, axis=1)
bin_magnitude_means_chirp = np.mean(cqt_magnitudes_chirp, axis=1)

fig, ax = plt.subplots(4, 2)
ax[0, 0].imshow(cqt_magnitudes, interpolation='bilinear')
ax[0, 0].set_aspect('auto')
ax[0, 0].set_title('Cqt magnitudes')
ax[0, 0].set_xlabel('i_block')
ax[0, 0].set_ylabel('i_bin')

ax[1, 0].plot(time_data, audio_input_data, label='input')
ax[1, 0].plot(time_data, audio_output_data, label='output')
ax[1, 0].set_xlabel('time in s')
ax[1, 0].set_ylabel('amplitude')
ax[1, 0].legend()

ax[0, 1].imshow(cqt_magnitudes_chirp, interpolation='bilinear')
ax[0, 1].set_aspect('auto')
ax[0, 1].set_title('Cqt magnitudes')
ax[0, 1].set_xlabel('i_block')
ax[0, 1].set_ylabel('i_bin')

ax[1, 1].plot(time_data, audio_input_data_chirp, label='input')
ax[1, 1].plot(time_data, audio_output_data_chirp, label='output')
ax[1, 1].set_xlabel('time in s')
ax[1, 1].set_ylabel('amplitude')
ax[1, 1].legend()

ax[2, 0].plot(freq, np.abs(audio_input_f))
ax[2, 0].plot(freq, np.abs(audio_output_f))
ax[2, 1].plot(freq, np.abs(audio_input_chirp_f))
ax[2, 1].plot(freq, np.abs(audio_output_chirp_f))
ax[2, 0].set_xlabel('frequency in Hz')
ax[2, 0].set_ylabel('magnitude')
ax[2, 1].set_xlabel('frequency in Hz')
ax[2, 1].set_ylabel('magnitude')
ax[2, 0].set_xscale('log')
ax[2, 0].set_xlim(10, sample_rate/2)
ax[2, 1].set_xscale('log')
ax[2, 1].set_xlim(10, sample_rate/2)

ax[3, 0].plot(bin_magnitude_means)
ax[3, 1].plot(bin_magnitude_means_chirp)
ax[3, 0].set_xlabel('i_bin')
ax[3, 1].set_xlabel('i_bin')
ax[3, 0].set_ylabel('mean magnitude')
ax[3, 1].set_ylabel('mean magnitude')
ax[3, 0].invert_xaxis()
ax[3, 1].invert_xaxis()

plt.show()



