import numpy as np
import pandas as pd
import math


def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (
        idx == len(array)
        or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])
    ):
        return array[idx - 1]
    else:
        return array[idx]  # 3d plotting fun


def load_br_file(fname):
    return np.array(pd.read_csv(fname, sep=",", header=None))


# Load BR file and return it as a dictionary
def load_BR_dict(fname, sec_offset=10, channels=["Left", "Right"]):
    txtdata = load_br_file(fname)[:, [0, 2]]

    return {
        chann: txtdata[-(422 * sec_offset) : -1, cidx]
        for cidx, chann in enumerate(channels)
    }
