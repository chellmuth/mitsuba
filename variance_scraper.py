import glob
import re

import numpy
import pyexr

WIDTH = 640
HEIGHT = 360
SAMPLES = 32

FEATURES = [
    {
        "name": "color",
        "channels": 3,
    },
    {
        "name": "albedo",
        "channels": 3,
    },
    {
        "name": "normal",
        "channels": 3,
    },
    {
        "name": "depth",
        "channels": 1,
    },
]

CHANNEL_COUNT = 10

def parse(log, buffer):
    pixel_re = re.compile("pixel: \\((\\d+), (\\d+)\\)\n")
    spectrum_re = re.compile("\\[(.*)\\]")

    while log:
        pixel_line = log.readline()
        if pixel_line == "": break

        pixel_match = pixel_re.match(pixel_line)
        x, y = [ int(coord) for coord in pixel_match.group(1, 2) ]

        for sample_index in range(SAMPLES):
            log.readline() # sample: \d+

            channel_base = 0
            for feature in FEATURES:
                spectrum_match = spectrum_re.match(log.readline())
                spectrum = [ float(value) for value in spectrum_match.group(1).split(",") ]

                for channel_index in range(feature["channels"]):
                    value = spectrum[channel_index]
                    buffer[x][y][channel_base + channel_index][sample_index] = value

                channel_base += feature["channels"]

def gather_log_paths(root_dir):
    return glob.glob(root_dir + "/multichannel_log*.txt")

def channel_name(feature):
    return feature["name"] + "Variance"

if __name__ == "__main__":
    buffer = numpy.empty((WIDTH, HEIGHT, CHANNEL_COUNT, SAMPLES), dtype=numpy.float)

    log_paths = gather_log_paths("runtime")
    print("Found %d logs" % len(log_paths))

    for i, log_path in enumerate(log_paths):
        parse(open(log_path), buffer)
        if i > 0 and i % 50 == 0:
            print("Parsed %d" % (i))

    print("Calculating variances")

    channels = {}

    for feature in FEATURES:
        channels[channel_name(feature)] = numpy.empty((HEIGHT, WIDTH, feature["channels"]))

    for x in range(HEIGHT):
        for y in range(WIDTH):
            channel_base = 0

            for feature in FEATURES:
                for channel_index in range(feature["channels"]):
                    value = numpy.var(buffer[y][x][channel_base + channel_index])
                    channels[channel_name(feature)][x][y][channel_base + channel_index] = value

    print("Writing EXR")
    pyexr.write("variance.exr", channels)
