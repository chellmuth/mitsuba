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

def _f_1(randoms_1, scene):
    _write_randoms([randoms_1], "randoms.bin")
    return _f(scene)[0]

def fa_1(randoms_1):
    return _f_1(randoms_1, "cjh.xml")

def fb_1(randoms_1):
    return _f_1(randoms_1, "staircase/cjh-scene.xml")


def _f_n(randoms_memoryview, scene):
    randoms_n = randoms_memoryview.tolist()

    _write_randoms(randoms_n, "randoms.bin")
    return _f(scene)

def fa_n(randoms_memoryview):
    return _f_n(randoms_memoryview, "cjh.xml")

def fb_n(randoms_memoryview):
    return _f_n(randoms_memoryview, "staircase/cjh-scene.xml")


# see order in runner.py or cjh.cpp
X_INDEX = 3
Y_INDEX = 4

def _render(rows, cols, spp, randoms_memoryview, scene):
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
                sample_index = ((cols * spp) * row) + (spp * col) + sample

                x_sample = randoms_n[sample_index][X_INDEX]
                y_sample = randoms_n[sample_index][Y_INDEX]

                randoms_n[sample_index][X_INDEX] = (float(col) / cols) + x_sample * x_delta
                randoms_n[sample_index][Y_INDEX] = (float(row) / rows) + y_sample * y_delta

    _write_randoms(randoms_n, "randoms.bin")
    fs = _f(scene)

    results = [ 0 for _ in range(rows * cols)]
    for row in range(rows):
        for col in range(cols):
            pixel_sum = 0
            for sample in range(spp):
                result_index = ((cols * spp) * row) + (spp * col) + sample
                pixel_sum += fs[result_index]

            pixel_value = pixel_sum / spp
            pixel_value = min(1, pixel_value ** (1/2.2))

            # transpose
            results[rows * col + row] = pixel_value

    return array.array("f", results)

def rendera(rows, cols, spp, randoms_memoryview):
    return _render(rows, cols, spp, randoms_memoryview, "cjh.xml")

def renderb(rows, cols, spp, randoms_memoryview):
    return _render(rows, cols, spp, randoms_memoryview, "staircase/cjh-scene.xml")

def _f(scene):
    return_code = subprocess.call(
        ["dist\\mitsuba", scene],
        env={"PATH": "./dist", "LD_LIBRARY_PATH": "./dist"}
    )
    if return_code != 0:
        raise BaseException("Failed call to Mitsuba")

    return read_luminances("luminance.bin")
