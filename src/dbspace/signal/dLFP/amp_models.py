import numpy as np


def hard_amp(xin, clip=1):
    xhc = np.copy(xin)
    xhc[np.where(xin > clip)] = clip
    xhc[np.where(xin < -clip)] = -clip

    return xhc
