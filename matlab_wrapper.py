from operator import mul

import array
import struct
import subprocess

def _write_randoms(randoms_n, out_path):
    out_file = open(out_path, "wb")

    bin_floats = []

    for randoms in randoms_n:
        bin_floats.append(struct.pack("f" * len(randoms), *randoms))

    out_file.write(struct.pack("I", len(bin_floats)))
    for bin_float in bin_floats:
        out_file.write(bin_float)

    out_file.close()

def read_luminances(in_path):
    f = open(in_path, "rb")

    luminances = []

    bytes = f.read(4)
    while bytes:
        luminance, = struct.unpack("f", bytes)
        luminances.append(luminance)
        bytes = f.read(4)

    return luminances

def f_1(randoms_1):
    _write_randoms([randoms_1], "randoms.bin")
    return _f()[0]

def f_n(randoms_memoryview):
    randoms_n = randoms_memoryview.tolist()
    _write_randoms(randoms_n, "randoms.bin")
    return _f()

# see order in runner.py or cjh.cpp
X_INDEX = 3
Y_INDEX = 4

def render(cols, rows, spp, randoms_memoryview):
    cols = int(cols)
    rows = int(rows)
    spp = int(spp)

    total_samples_needed = cols * rows * spp

    total_samples_provided, _ = randoms_memoryview.shape
    if total_samples_needed != total_samples_provided:
        raise BaseException(
            "Need (%d x %d x %d) = %d samples. Received %d."
            % (cols, rows, spp, total_samples_needed, total_samples_provided)
        )

    randoms_n = randoms_memoryview.tolist()

    x_delta = 1. / cols
    y_delta = 1. / rows
    for row in range(rows):
        for col in range(cols):
            for sample in range(spp):
                sample_index = cols * row + col + sample
                x_sample = randoms_n[sample_index][X_INDEX]
                y_sample = randoms_n[sample_index][Y_INDEX]
                randoms_n[sample_index][X_INDEX] = (float(col) / cols) + x_sample * x_delta
                randoms_n[sample_index][Y_INDEX] = (float(row) / rows) + y_sample * x_delta

    _write_randoms(randoms_n, "randoms.bin")
    return array.array("f", _f())

def _f():
    return_code = subprocess.call(
        ["dist\\mitsuba", "cjh.xml"],
        env={"PATH": "./dist", "LD_LIBRARY_PATH": "./dist"}
    )
    if return_code != 0:
        raise BaseException("Failed call to Mitsuba")

    return read_luminances("luminance.bin")
