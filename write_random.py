import struct

def run(in_path, out_path):
    in_file = open("randoms.txt", "r")
    out_file = open("randoms.bin", "wb")

    bin_floats = []

    for line in in_file:
        floats = [ float(f) for f in line.strip().split(" ") ]
        bin_floats.append(struct.pack("f" * len(floats), *floats))

    out_file.write(struct.pack("I", len(bin_floats)))
    for bin_float in bin_floats:
        out_file.write(bin_float)

    out_file.close()

if __name__ == "__main__":
    run("randoms.txt", "randoms.bin")

