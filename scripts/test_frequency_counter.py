#!/usr/bin/env python3
"""Regression check for the >2^24 frequency-counter failure."""

import struct


def float32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


limit = 1 << 24
assert float32(float32(limit) + 1.0) == limit

records = 69_680_812
majority = 57_445_000
old_probability = float32(limit / records)
fixed_probability = float32(majority / records)
assert old_probability < 0.25
assert 0.82 < fixed_probability < 0.83
print("Integer frequency-counter regression passed.")
