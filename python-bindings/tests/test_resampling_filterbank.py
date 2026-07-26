import numpy as np
import pytest

import prtcqt


SAMPLE_RATE = 48_000.0
CALLBACK_SIZE = 64
SAMPLE_COUNT = 2048
STAGE_COUNT = 9


def process(partitions):
    filterbank = prtcqt.ResamplingFilterbank9()
    filterbank.init(SAMPLE_RATE, CALLBACK_SIZE)

    generator = np.random.default_rng(0x5EED)
    input_signal = generator.uniform(-1.0, 1.0, SAMPLE_COUNT)
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
        [np.asarray(stage) for stage in stages],
        np.asarray(output),
        filterbank.get_processing_block_size(),
        filterbank.get_latency_samples(),
    )


def test_callback_partition_invariance():
    regular_stages, regular_output, processing_size, latency = process(
        [CALLBACK_SIZE]
    )
    irregular_stages, irregular_output, _, _ = process([1, 7, 31, 5, 64, 13])

    assert processing_size == 256
    assert latency == 256
    for stage, (regular, irregular) in enumerate(
        zip(regular_stages, irregular_stages)
    ):
        assert regular.size == SAMPLE_COUNT >> stage
        np.testing.assert_allclose(irregular, regular, rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(
        irregular_output, regular_output, rtol=0.0, atol=1.0e-14
    )


def test_reported_latency_is_zero_filled():
    _, output, _, latency = process([CALLBACK_SIZE])

    np.testing.assert_array_equal(output[:latency], np.zeros(latency))
    assert np.any(np.abs(output[latency:]) > 1.0e-12)


def test_rejects_oversized_callback():
    filterbank = prtcqt.ResamplingFilterbank9()
    filterbank.init(SAMPLE_RATE, CALLBACK_SIZE)

    with pytest.raises(ValueError):
        filterbank.process([0.0] * (CALLBACK_SIZE + 1))
