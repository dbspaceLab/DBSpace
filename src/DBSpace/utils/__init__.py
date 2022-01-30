from collections import defaultdict

""" Powerful structure that literally changed my life"""


def nestdict():
    return defaultdict(nestdict)


""" Using Jack-knife for the median calculation"""


def jk_median(inX, niter=100):
    # assume first dim is SEGMENTS/OBSERVATIONS
    med_vec = []
    for ii in range(niter):
        choose_idxs = random.sample(
            range(0, inX.shape[0]), np.floor(inX.shape[0] / 2).astype(np.int)
        )
        med_vec.append(np.median(inX[choose_idxs, :], axis=0))

    return np.array(med_vec)


def cmedian(inArray, axis=-1):
    return np.median(np.real(inArray), axis=axis) + 1j * np.median(
        np.imag(inArray), axis=axis
    )


def l2_pow(x):
    return np.sqrt(np.sum(x ** 2))
