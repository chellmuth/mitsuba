import re
import subprocess

import numpy as np
from PIL import Image

def run():
    cols = 80
    rows = 80

    args = []
    values = []
    for y in range(rows):
        args.append("-Dy=%f" % (float(y) / (rows-1)))
        for x in range(cols):
            args.append("-Dx=%f" % (float(x) / (cols-1)))

            stdout = subprocess.check_output(["mitsuba", "cjh.xml"] + args).decode("UTF-8")
            match = re.search("luminance = (\d\.?\d*)", stdout)
            if not match:
                import pdb; pdb.set_trace()

            luminance = float(match.group(1))
            values.append(luminance)

    byte_values = (np.array(values) * 255).astype("uint8")
    
    image = Image.frombytes("L", (cols, rows), byte_values)
    image.save("output.png")

    print("done!")

if __name__ == "__main__":
    run()
