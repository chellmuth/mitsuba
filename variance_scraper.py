import glob
import re

import numpy
import pyexr

WIDTH = 1280
HEIGHT = 720
SAMPLES = 16

ALBEDO_R = 0
ALBEDO_G = 1
ALBEDO_B = 2
NORMAL_R = 3
NORMAL_G = 4
NORMAL_B = 5
DISTANCE = 6
CHANNEL_COUNT = 7

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
            log.readline() # path traced radiance

            albedo_match = spectrum_re.match(log.readline())
            r, g, b = [ float(value) for value in albedo_match.group(1).split(",") ]

            buffer[x][y][ALBEDO_R][sample_index] = r
            buffer[x][y][ALBEDO_G][sample_index] = g
            buffer[x][y][ALBEDO_B][sample_index] = b

            normal_match = spectrum_re.match(log.readline())
            r, g, b = [ float(value) for value in normal_match.group(1).split(",") ]

            buffer[x][y][NORMAL_R][sample_index] = r
            buffer[x][y][NORMAL_G][sample_index] = g
            buffer[x][y][NORMAL_B][sample_index] = b

            distance_match = spectrum_re.match(log.readline())
            r, g, b = [ float(value) for value in distance_match.group(1).split(",") ]

            buffer[x][y][DISTANCE][sample_index] = r

def gather_log_paths(root_dir):
    return glob.glob(root_dir + "/multichannel_log*.txt")

if __name__ == "__main__":
    buffer = numpy.empty((WIDTH, HEIGHT, CHANNEL_COUNT, SAMPLES), dtype=numpy.float)

    log_paths = gather_log_paths("runtime")
    print("Found %d logs" % len(log_paths))

    for i, log_path in enumerate(log_paths):
        parse(open(log_path), buffer)
        if i > 0 and i % 50 == 0:
            print("Parsed %d" % (i))

    print("Calculating variances")
    channels = {
        "albedo": numpy.empty((HEIGHT, WIDTH, 3)),
        "normal": numpy.empty((HEIGHT, WIDTH, 3)),
        "depth": numpy.empty((HEIGHT, WIDTH, 1)),
    }

    for x in range(HEIGHT):
        for y in range(WIDTH):
            for channel_index in range(3):
                channels["albedo"][x][y][channel_index] = \
                    numpy.var(buffer[y][x][ALBEDO_R + channel_index])

            for channel_index in range(3):
                channels["normal"][x][y][channel_index] = \
                    numpy.var(buffer[y][x][NORMAL_R + channel_index])

            for channel_index in range(1):
                channels["depth"][x][y][channel_index] = \
                    numpy.var(buffer[y][x][DISTANCE + channel_index])

    print("Writing EXR")
    pyexr.write("variance.exr", channels)
