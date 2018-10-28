import random as r
import re
import subprocess

import numpy as np
from PIL import Image

def run():
    cols = 500
    rows = 500
    samples_per_coordinate = 20

    args = []
    values = []

    f = open("randoms.txt", "w")

    x_delta = 1. / cols
    y_delta = 1. / rows

    for y in range(rows):
        for x in range(cols):
            for _ in range(samples_per_coordinate):
                args = [
                    ("-Du=%f", 0.2),
                    ("-Dv=%f", 0.8),
                    ("-Dx=%f", (float(x) / cols) + r.random() * x_delta),
                    ("-Dy=%f", (float(y) / rows) + r.random() * y_delta),
                    ("-Ddirect1_1=%f", r.random()),
                    ("-Ddirect1_2=%f", r.random()),
                    ("-Dbsdf1_1=%f", r.random()),
                    ("-Dbsdf1_2=%f", r.random()),
                    ("-Ddirect2_1=%f", r.random()),
                    ("-Ddirect2_2=%f", r.random()),
                    ("-Dbsdf2_1=%f", r.random()),
                    ("-Dbsdf2_2=%f", r.random()),
                ]

                values = [ str(val) for arg, val in args ]
                f.write(" ".join(values))
                f.write("\n")

    f.close()

    print("done!")

if __name__ == "__main__":
    run()
