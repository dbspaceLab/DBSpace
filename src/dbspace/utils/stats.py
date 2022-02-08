import numpy as np
import random


def jk_median(input_signal: np.ndarray, niter: int = 100):
    """
    Using Jack-knife for the median calculation


    Input
    -----
        input_signal : np.ndarray
            assumed to be M x ..., with M is the axis for the median calculation.

    """
    # assume first dim is SEGMENTS/OBSERVATIONS
    med_vec = np.full((niter,), np.nan)
    for ii in range(niter):
        choose_idxs = random.sample(
            range(0, input_signal.shape[0]),
            np.floor(input_signal.shape[0] / 2).astype(np.int),
        )
        med_vec[ii, ...] = np.median(input_signal[choose_idxs, :], axis=0)

    return np.array(med_vec)


def complex_median(inArray: np.ndarray, axis=-1):
    return np.median(np.real(inArray), axis=axis) + 1j * np.median(
        np.imag(inArray), axis=axis
    )


def pca(data, num_comps: int = None, mean_correct: bool = False):
    """
    Very simple PCA implementation to quickly get numComps components


    """
    if mean_correct:
        data -= data.mean(axis=0)

    R = np.cov(data, rowvar=False)
    evals, evecs = np.linalg.eigh(R)
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:, idx]
    evals = evals[idx]

    if num_comps is not None:
        evecs = evecs[:, :num_comps]

    return np.dot(evecs.T, data.T).T, evals, evecs
