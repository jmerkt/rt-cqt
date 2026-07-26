import math

import prtcqt
import pytest


BLOCK_SIZE = 256
OCTAVE_NUMBER = 9


@pytest.mark.parametrize(
    ("class_name", "bin_count"),
    [("Cqt12", 12), ("Cqt24", 24)],
)
def test_constant_q_binding_uses_composed_transform(class_name, bin_count):
    transform = getattr(prtcqt, class_name)()
    transform.init(128)
    transform.init_sample_rate(48_000, BLOCK_SIZE)
    transform.input_block([0.0] * BLOCK_SIZE)

    schedule = transform.get_cqt_schedule()
    for element in schedule:
        transform.cqt(element)
        assert len(transform.get_octave_cqt_buffer(element.octave())) == bin_count
        transform.icqt(element)

    output = transform.output_block(BLOCK_SIZE)
    assert len(output) == BLOCK_SIZE
    assert all(math.isfinite(sample) for sample in output)


@pytest.mark.parametrize(
    ("class_name", "bin_count"),
    [("SlidingCqt12", 12), ("SlidingCqt24", 24)],
)
def test_sliding_cqt_binding_uses_composed_transform(class_name, bin_count):
    transform = getattr(prtcqt, class_name)()
    transform.init(48_000, BLOCK_SIZE)
    transform.input_block([0.0] * BLOCK_SIZE, BLOCK_SIZE)

    values = transform.get_octave_values(OCTAVE_NUMBER - 1)
    frequencies = transform.get_octave_bin_frequencies(OCTAVE_NUMBER - 1)
    output = transform.output_block(BLOCK_SIZE)

    assert len(values) == bin_count
    assert len(frequencies) == bin_count
    assert len(output) == BLOCK_SIZE
    assert all(math.isfinite(value.real) and math.isfinite(value.imag) for value in values)
    assert all(math.isfinite(frequency) for frequency in frequencies)
    assert all(math.isfinite(sample) for sample in output)
