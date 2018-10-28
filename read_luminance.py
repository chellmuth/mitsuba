import struct

import numpy as np
from PIL import Image

cols = 500
rows = 500
samples_per_coordinate = 20

def run():
    values = []
    f = open("luminance.bin", "rb")

    for _ in range(cols * rows):
        local_luminance = []
        for _ in range(samples_per_coordinate):
            luminance, = struct.unpack("f", f.read(4))
            local_luminance.append(luminance)

        values.append(np.mean(local_luminance))

    byte_values = (np.array(values) * 255).astype("uint8")

    image = Image.frombytes("L", (cols, rows), byte_values)
    image.save("output.png")

if __name__ == "__main__":
    run()
