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

import rtcqt as cqt

bins_per_octave: int = 12
n_octaves: int = 9
n_bins: int = bins_per_octave * n_octaves
block_size: int = 256
sample_rate: float = 48000.
n_blocks: int = 100

cqt_instance = cqt.SlidingCqt12()
cqt_instance.init(sample_rate, block_size)

# Synthesize single sine tone
time_data = np.linspace(0.0, block_size / sample_rate * n_blocks, block_size * n_blocks)
tone_frequency: float = 440.0
sine_data = np.sin(2.0 * np.pi * time_data * tone_frequency)

# Evaluate using cqt
cqt_magnitudes = np.zeros((n_blocks, n_bins), np.float32)
for i_block in range(n_blocks):
    cqt_instance.inputBlock(sine_data[i_block * block_size: (i_block + 1) * block_size], block_size)
    for i_octave in range(n_octaves):
        octave_data = cqt_instance.getOctaveValues(i_octave) 
        cqt_magnitudes[i_block, i_octave * bins_per_octave : (i_octave + 1) * bins_per_octave] = np.abs(octave_data)

# Plot result
cqt_magnitudes  = cqt_magnitudes.transpose()

fig, ax = plt.subplots(1)
ax.imshow(cqt_magnitudes, interpolation='bilinear')
ax.set_aspect('auto')
ax.set_title('Cqt magnitudes')
ax.set_xlabel('i_block')
ax.set_ylabel('i_bin')


# TODO: further evaluation with not just block-resampled scqt results, but actual samples on octave rates
# TODO: forward / backward transform and check input vs. output

