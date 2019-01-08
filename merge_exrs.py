import pyexr

import numpy as np

def merge(mitsuba_exr, variance_exr):
    merged = {}
    for channel in variance_exr.root_channels:
        merged[channel] = variance_exr.get(channel)

    for channel in mitsuba_exr.root_channels:
        output_channel = channel
        if channel == "color":
            output_channel = "default"

        merged[output_channel] = mitsuba_exr.get(channel)

    return merged

def mock_diffuse_and_specular(data):
    shape = data["default"].shape
    data["diffuse"] = data["default"]
    data["diffuseVariance"] = data["colorVariance"]

    data["specular"] = np.empty(shape)
    data["specularVariance"] = np.empty(shape)

    return data

if __name__ == "__main__":
    data = merge(
        pyexr.open("/home/cjh/scenes/mitsuba/bathroom2/scene.exr"),
        pyexr.open("variance.exr")
    )

    data = mock_diffuse_and_specular(data)

    pyexr.write("merged.exr", data)
