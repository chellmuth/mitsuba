import subprocess
import struct

def _write_randoms(randoms, out_path):
    out_file = open(out_path, "wb")

    bin_floats = []

    bin_floats.append(struct.pack("f" * len(randoms), *randoms))

    out_file.write(struct.pack("I", len(bin_floats)))
    for bin_float in bin_floats:
        out_file.write(bin_float)

    out_file.close()

def read_luminance(in_path):
    f = open(in_path, "rb")

    luminance, = struct.unpack("f", f.read(4))
    return luminance

def f(args):
    _write_randoms(args, "randoms.bin")
    subprocess.call(
        ["mitsuba", "cjh.xml"],
        env={"PATH": "./dist", "LD_LIBRARY_PATH": "./dist"}
    )
    return read_luminance("luminance.bin")
