#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 20:32:31 2023

@author: jmerkt
"""

import sys
sys.path.insert(1, './../build')

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import chirp

import rtcqt as cqt

bins_per_octave: int = 24
n_octaves: int = 9
n_bins: int = bins_per_octave * n_octaves
block_size: int = 512
sample_rate: float = 48000.
n_blocks: int = 2000

cqt_instance = cqt.SlidingCqt24()
cqt_instance.init(sample_rate, block_size)

# Synthesize single sine tone
time_data = np.linspace(0.0, block_size / sample_rate * n_blocks, block_size * n_blocks)
f0: float = 440.0
audio_input_data = np.sin(2.0 * np.pi * time_data * f0, dtype=np.float64)

# Evaluate using cqt
cqt_magnitudes = np.zeros((n_blocks, n_bins), np.float32)
audio_output_data = np.zeros(time_data.shape[0], np.float64)
for i_block in range(n_blocks):
    cqt_instance.inputBlock(audio_input_data[i_block * block_size: (i_block + 1) * block_size], block_size)
    for i_octave in range(n_octaves):
        octave_data = cqt_instance.getOctaveValues(i_octave) 
        cqt_magnitudes[i_block, i_octave * bins_per_octave : (i_octave + 1) * bins_per_octave] = np.flip(np.abs(octave_data))
    output_block = cqt_instance.outputBlock(block_size)
    audio_output_data[i_block * block_size: (i_block + 1) * block_size] = output_block
    
    
# Synthesize exponential sine sweep
f0 = 30.0
f1 = 16000.0
audio_input_data_chirp = chirp(time_data, f0=f0, t1=time_data[-1], f1=f1, method='logarithmic') 

cqt_magnitudes_chirp = np.zeros((n_blocks, n_bins), np.float32)
audio_output_data_chirp = np.zeros(time_data.shape[0], np.float64)
for i_block in range(n_blocks):
    cqt_instance.inputBlock(audio_input_data_chirp[i_block * block_size: (i_block + 1) * block_size], block_size)
    for i_octave in range(n_octaves):
        octave_data = cqt_instance.getOctaveValues(i_octave) 
        cqt_magnitudes_chirp[i_block, i_octave * bins_per_octave : (i_octave + 1) * bins_per_octave] = np.flip(np.abs(octave_data))
    output_block = cqt_instance.outputBlock(block_size)
    audio_output_data_chirp[i_block * block_size: (i_block + 1) * block_size] = output_block
    

    
# Plot result
cqt_magnitudes  = cqt_magnitudes.transpose()
cqt_magnitudes_chirp = cqt_magnitudes_chirp.transpose()

fig, ax = plt.subplots(2, 2)
ax[0, 0].imshow(cqt_magnitudes, interpolation='bilinear')
ax[0, 0].set_aspect('auto')
ax[0, 0].set_title('Cqt magnitudes')
ax[0, 0].set_xlabel('i_block')
ax[0, 0].set_ylabel('i_bin')

ax[1, 0].plot(time_data, audio_input_data, label='input')
ax[1, 0].plot(time_data, audio_output_data, label='output')
ax[1, 0].set_xlabel('time')
ax[1, 0].set_ylabel('amplitude')
ax[1, 0].legend()

ax[0, 1].imshow(cqt_magnitudes_chirp, interpolation='bilinear')
ax[0, 1].set_aspect('auto')
ax[0, 1].set_title('Cqt magnitudes')
ax[0, 1].set_xlabel('i_block')
ax[0, 1].set_ylabel('i_bin')

ax[1, 1].plot(time_data, audio_input_data_chirp, label='input')
ax[1, 1].plot(time_data, audio_output_data_chirp, label='output')
ax[1, 1].set_xlabel('time')
ax[1, 1].set_ylabel('amplitude')
ax[1, 1].legend()


# TODO: further evaluation with not just block-resampled scqt results, but actual samples on octave rates

