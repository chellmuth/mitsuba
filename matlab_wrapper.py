import subprocess
import struct

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

def _f():
    return_code = subprocess.call(
        ["mitsuba", "staircase/cjh-scene.xml"],
        env={"PATH": "./dist", "LD_LIBRARY_PATH": "./dist"}
    )
    if return_code != 0:
        raise BaseException("Failed call to Mitsuba")

    return read_luminances("luminance.bin")
