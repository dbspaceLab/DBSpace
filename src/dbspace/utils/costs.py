import numpy as np


def l2_pow(x: np.ndarray):
    """
    Return the L2 norm (power) of the input vector
    """

    return np.sqrt(np.sum(x**2))
