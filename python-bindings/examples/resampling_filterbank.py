#!/usr/bin/env python3
"""Plot the block-only resampling filterbank's impulse response.

Build the Python bindings, then run this file from its directory:

    cd ..
    uv sync
    uv run cmake -S . -B build -DPython_EXECUTABLE=.venv/bin/python
    uv run cmake --build build --target prtcqt
    uv run python examples/resampling_filterbank.py

The filterbank wrapper echoes every analysis stage into the corresponding
synthesis stage. Consequently the reconstructed signal visualizes the complete
multirate path, but it is not expected to equal the input signal.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


THIS_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIRECTORY.parent / "build"))
import prtcqt


SAMPLE_RATE = 48_000.0
CALLBACK_SIZE = 64
SAMPLE_COUNT = 4096
STAGE_COUNT = 9


def process(partitions):
    filterbank = prtcqt.ResamplingFilterbank9()
    filterbank.init(SAMPLE_RATE, CALLBACK_SIZE)

    input_signal = np.zeros(SAMPLE_COUNT, dtype=np.float64)
    input_signal[0] = 1.0
    stages = [[] for _ in range(STAGE_COUNT)]
    output = []

    position = 0
    partition_index = 0
    while position < input_signal.size:
        block_size = min(
            partitions[partition_index % len(partitions)],
            input_signal.size - position,
        )
        stage_blocks, output_block = filterbank.process(
            input_signal[position : position + block_size].tolist()
        )
        for stage, stage_block in enumerate(stage_blocks):
            stages[stage].extend(stage_block)
        output.extend(output_block)
        position += block_size
        partition_index += 1

    return (
        input_signal,
        [np.asarray(stage) for stage in stages],
        np.asarray(output),
        filterbank.get_processing_block_size(),
        filterbank.get_latency_samples(),
    )


input_signal, stages, output, processing_size, latency = process([CALLBACK_SIZE])
_, irregular_stages, irregular_output, _, _ = process([1, 7, 31, 5, 64, 13])

stage_error = max(
    np.max(np.abs(regular - irregular))
    for regular, irregular in zip(stages, irregular_stages)
)
output_error = np.max(np.abs(output - irregular_output))
print(f"processing block: {processing_size} samples")
print(f"reported latency: {latency} samples ({latency / SAMPLE_RATE * 1000:.2f} ms)")
print(f"maximum stage partition error: {stage_error:.3e}")
print(f"maximum output partition error: {output_error:.3e}")

figure, axes = plt.subplots(6, 2, figsize=(12, 15))
figure.suptitle("Block-only resampling filterbank impulse response")

axes[0, 0].plot(input_signal[: processing_size * 2])
axes[0, 0].set_title("Input impulse")
axes[0, 1].plot(output[: processing_size * 4])
axes[0, 1].axvline(latency, color="tab:red", linestyle="--", label="reported latency")
axes[0, 1].set_title("Echoed-stage synthesis")
axes[0, 1].legend()

for stage, axis in enumerate(axes.flat[2 : 2 + STAGE_COUNT]):
    stage_data = stages[stage]
    shown_samples = min(256, stage_data.size)
    axis.plot(stage_data[:shown_samples])
    axis.set_title(
        f"Stage {stage}: {SAMPLE_RATE / (2 ** stage):g} Hz, "
        f"{stage_data.size} samples"
    )

for axis in axes.flat:
    axis.set_xlabel("sample")
    axis.grid(alpha=0.2)

for axis in axes.flat[2 + STAGE_COUNT :]:
    axis.set_visible(False)

figure.tight_layout()
output_path = THIS_DIRECTORY / "resampling_filterbank.png"
figure.savefig(output_path, dpi=150)
print(f"wrote {output_path}")
