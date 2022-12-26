import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
from dbspace.utils.costs import l2_pow
from dbspace.utils.structures import nestdict

#%%
# Basic functions for rotating recordings into particular frames


def gen_T(inpX, Fs=422, nfft=2**10):
    outT = nestdict(dict)
    for chann in inpX.keys():
        outT[chann] = {
            "T": np.linspace(0, inpX[chann].shape[0] / Fs, inpX[chann].shape[0]),
            "V": inpX[chann],
        }

    return outT


""" gen_psd outputs a PSD, not a LogPSD """
""" This function WRAPS F_Domain"""


def gen_psd(inpX, Fs=422, nfft=2**10, polyord=0):
    # inp X is going to be assumed to be a dictionary with different keys for different channels
    outPSD = nestdict()
    outPoly = nestdict()
    # assume input is time x seg
    for chann in inpX.keys():
        # The return here is a dictionary with two keys: F and PSD
        # check the size of the matrix now; it could be that we need to do this many times for each "segment"
        fmatr = np.zeros((inpX[chann].shape[-1], int(nfft / 2) + 1))
        polysub = np.zeros((inpX[chann].shape[-1], polyord + 1))

        if inpX[chann].ndim > 1:
            for seg in range(inpX[chann].shape[-1]):

                psd = np.abs(
                    F_Domain(inpX[chann][:, seg].squeeze(), Fs=Fs, nfft=nfft)["Pxx"]
                )  # Just enveloped this with np.abs 12/15/2020

                fmatr[seg, :] = psd
        else:
            psd = F_Domain(inpX[chann], Fs=Fs, nfft=nfft)["Pxx"]
            fmatr = psd

        outPSD[chann] = fmatr.squeeze()

        # do polysub here
    if polyord != 0:
        print("Polynomial Correcting Stack")
        outPSD = poly_subtr(outPSD, np.linspace(0, Fs / 2, nfft / 2 + 1))

    # Return here is a dictionary with Nchann keys
    return outPSD


#%%
"""Below used to be called poly_subtrLFP, unclear whether it was being used, now renamed and look for errors elsewhere"""


def poly_SG(inSG, fVect, order=4):
    out_sg = np.zeros_like(inSG)

    for seg in range(inSG.shape[1]):
        inpsd = 10 * np.log10(inSG[chann][seg, :])
        polyCoeff = np.polyfit(fVect, inpsd, order)
        polyfunc = np.poly1d(polyCoeff)
        polyitself = polyfunc(fVect)
        out_sg[:, seg] = 10 ** ((curr_psd - polyitself) / 10)

    return out_sg


def gen_SG(inpX, Fs=422, nfft=2**10, plot=False, overlap=True):
    outSG = nestdict()
    for chann in inpX.keys():
        if overlap == True:
            outSG[chann] = TF_Domain(inpX[chann])
        else:
            outSG[chann] = TF_Domain(inpX[chann], noverlap=0, nperseg=422 * 2)

    if plot:
        plot_TF(outSG, chs=inpX.keys())

    return outSG


#%%
# Function to go through and find all the features from the PSD structure of dbo
def calc_feats(psdIn, yvect, dofeats="", modality="eeg", compute_method="median"):
    # psdIn is a VECTOR, yvect is the basis vector
    if dofeats == "":
        dofeats = DEFAULT_FEAT_ORDER

    if psdIn.ndim < 2:
        psdIn = psdIn[..., np.newaxis]

    if modality == "eeg":
        ch_list = np.arange(0, 257)
    elif modality == "lfp":
        ch_list = ["Left", "Right"]

    feat_vect = []
    for feat in dofeats:
        # print(feat_dict[feat]['param'])
        # dofunc = feat_dict[feat]['fn']
        if compute_method == "median":
            computed_featinspace = FEAT_DICT[feat]["fn"](
                psdIn, yvect, FEAT_DICT[feat]["param"]
            )
        elif compute_method == "mean":
            computed_featinspace = FEAT_DICT[feat]["fn"](
                psdIn, yvect, FEAT_DICT[feat]["param"], cmode=np.mean
            )

        cfis_matrix = [computed_featinspace[ch] for ch in ch_list]
        feat_vect.append(cfis_matrix)
        # feat_dict[feat] = dofunc['fn'](datacontainer,yvect,dofunc['param'])[0]

    feat_vect = np.array(feat_vect).squeeze()

    return feat_vect, dofeats


# Convert a feat dict that comes from a get feature function (WHERE IS IT?!)
def featDict_to_Matr(featDict):
    # structure of feat dict is featDict[FEATURE][CHANNEL] = VALUE
    ret_matr = np.array(
        [
            (featDict[feat]["Left"], featDict[feat]["Right"])
            for feat in DEFAULT_FEAT_ORDER
        ]
    )

    # assert that the size is as expected?
    # should be number of feats x number of channels!
    assert ret_matr.shape == (len(DEFAULT_FEAT_ORDER), 2)

    return ret_matr


def get_pow(Pxx, F, frange, cmode=np.median):
    # Pxx is a dictionary where the keys are the channels, the values are the [Pxx desired]
    # Pxx is assumed to NOT be log transformed, so "positive semi-def"

    # check if Pxx is NOT a dict
    if isinstance(Pxx, np.ndarray):
        # Pxx = Pxx.reshape(-1,1)
        # JUST ADDED THIS
        chann_order = range(Pxx.shape[0])
        Pxx = {ch: Pxx[ch, :] for ch in chann_order}

        # ThIS WAS WORKING BEFORE
        # Pxx = {0:Pxx}
    elif len(Pxx.keys()) > 2:
        chann_order = np.arange(0, 257)
    else:
        chann_order = ["Left", "Right"]

    # find the power in the range of the PSD
    # Always assume PSD is a dictionary of channels, and each value is a dictionary with Pxx and F

    # frange should just be a tuple with the low and high bounds of the band
    out_feats = {keys: 0 for keys in Pxx.keys()}

    Fidxs = np.where(np.logical_and(F > frange[0], F < frange[1]))[0]

    # for chans,psd in Pxx.items():
    for cc, chann in enumerate(chann_order):
        # let's make sure the Pxx we're dealing with is as expected and a true PSD
        assert (Pxx[chann] > 0).all()

        # if we want the sum
        # out_feats[chans] = np.sum(psd[Fidxs])
        # if we want the MEDIAN instead

        # log transforming this makes sense, since we find the median of the POLYNOMIAL CORRECTED Pxx, which is still ALWAYS positive
        out_feats[chann] = 10 * np.log10(cmode(Pxx[chann][Fidxs]))

    # return is going to be a dictionary with same elements

    return out_feats  # This returns the out_feats which are 10*log(Pxx)


def get_slope(Pxx, F, params):
    # method to get the fitted polynomial for the range desired
    frange = params["frange"]
    linorder = params["linorder"]

    if isinstance(Pxx, np.ndarray):
        Pxx = {0: Pxx}

    out_feats = {keys: 0 for keys in Pxx.keys()}

    Fidxs = np.where(np.logical_and(F > frange[0], F < frange[1]))

    for chans, psd in Pxx.items():
        logpsd = np.log10(psd[Fidxs])
        logF = np.log10(F[Fidxs])

        fitcoeffs = np.polyfit(logF, logpsd, linorder)

        out_feats[chans] = fitcoeffs[-linorder]

    # return is going to be a dictionary with same channel keys
    return out_feats


def get_ratio(Pxx, F, f_r_set, cmode=np.median):
    bandpow = [None] * len(f_r_set)
    # first get the power for each of the individual bands
    for bb, frange in enumerate(f_r_set):
        bandpow[bb] = get_pow(Pxx, F, frange, cmode=cmode)

    ret_ratio = {ch: bandpow[1][ch] / bandpow[0][ch] for ch in bandpow[0].keys()}
    return ret_ratio


def F_Domain(timeseries, nperseg=512, noverlap=128, nfft=2**10, Fs=422):

    # assert isinstance(timeser,dbs.timeseries)
    # Window size is about 1 second (512 samples is over 1 sec)

    # what are the dimensions of the timeser we're dealing with?

    Fvect, Pxx = sig.welch(
        timeseries,
        Fs,
        window="blackmanharris",
        nperseg=nperseg,
        noverlap=noverlap,
        nfft=nfft,
    )

    FreqReturn = {"F": Fvect, "Pxx": Pxx}

    return FreqReturn


def TF_Domain(
    timeseries: np.ndarray,
    fs: int = 422,
    nperseg: int = 2**10,
    noverlap: int = 2**10 - 50,
):
    # raise Exception
    # assert isinstance(timeser,dbs.timeseries)
    F, T, SG = sig.spectrogram(
        timeseries,
        nperseg=nperseg,
        noverlap=noverlap,
        window=sig.get_window("blackmanharris", nperseg),
        fs=fs,
    )

    TFreqReturn = {"T": T, "F": F, "SG": SG}

    return TFreqReturn


def poly_subtr(input_psd: np.ndarray, fvect: np.ndarray = None, polyord: int = 4):
    # This function takes in a raw PSD, Log transforms it, poly subtracts, and then returns the unloged version.
    # log10 in_psd first
    if fvect is None:
        fvect = np.linspace(0, 1, input_psd.shape[0])

    log_psd = 10 * np.log10(input_psd)
    pfit = np.polyfit(fvect, log_psd, polyord)
    pchann = np.poly1d(pfit)

    bl_correction = pchann(fvect)

    return 10 ** ((log_psd - bl_correction) / 10), pfit


def grab_median_psd(
    TFcont,
    bigmed,
    osc_feat,
    tlim=(880, 900),
    title="",
    do_corr=True,
    band_compute="median",
    band_scheme="Adjusted",
):
    # Plot some PSDs
    # plt.figure()
    chann_label = ["Left", "Right"]
    pf_lPSD = nestdict()

    if do_corr:
        psd_lim = (-20, 50)
    else:
        psd_lim = (-220, -70)

    # Make the big figure that will have both channels
    plt.figure(bigmed.number)
    for cc in range(2):
        chann = chann_label[cc]
        plt.subplot(2, 2, cc + 1)
        T = TFcont["TF"]["T"]
        F = TFcont["TF"]["F"]
        SG = TFcont["TF"]["SG"]

        t_idxs = np.where(np.logical_and(T > tlim[0], T < tlim[1]))

        med_psd = np.median(10 * np.log10(SG[chann][:, t_idxs]).squeeze(), axis=1)
        var_psd = np.var(10 * np.log10(SG[chann][:, t_idxs]).squeeze(), axis=1).reshape(
            -1, 1
        )
        corr_psd = {chann_label[cc]: 10 ** (med_psd / 10)}

        if do_corr:
            # do polynomial subtraction
            fixed_psd, polyitself = poly_subtr(corr_psd, F)
            pf_lPSD[chann_label[cc]] = fixed_psd[chann_label[cc]].reshape(-1, 1)
        else:
            correct_psd, polyitself = poly_subtr(corr_psd, F)

            pf_lPSD[chann_label[cc]] = 10 ** (med_psd / 10).reshape(-1, 1)
            plt.plot(F, polyitself, label="Polynomial Fit", color="black")

        plt.plot(F, 10 * np.log10(pf_lPSD[chann_label[cc]]), label=title)
        plt.title("Channel " + chann_label[cc] + " psd")
        # try: plt.fill_between(F,(10*np.log10(pf_lPSD[chann_label[cc]]))+var_psd,(10*np.log10(pf_lPSD[chann_label[cc]]))-var_psd)

        plt.ylim(psd_lim)

        plt.subplot(2, 2, 2 + (cc + 1))
        plt.plot(F, 10 * np.log10(var_psd), label=title)
        plt.title("Variance in PSD across time: " + chann_label[cc])
    plt.subplot(2, 2, 4)
    plt.legend()

    if band_scheme == "Standard":
        band_wins = ["Delta", "Theta", "Alpha", "Beta", "Gamma"]
    elif band_scheme == "Adjusted":
        band_wins = ["Delta", "Theta", "Alpha", "Beta*", "Gamma1"]

    fcalced, bands = calc_feats(
        pf_lPSD, F, dofeats=band_wins, modality="lfp", compute_method=band_compute
    )

    plt.figure(osc_feat.number)
    plt.subplot(1, 2, 1)
    plt.plot(fcalced[:, 0], label=title)

    plt.title("Left")
    plt.subplot(1, 2, 2)
    plt.plot(fcalced[:, 1], label=title)
    plt.title("Right")
    plt.suptitle("Features " + band_compute + " " + band_scheme)


# Function to go through and find all the features from the PSD structure of dbo
def calc_feats(psdIn, yvect, dofeats="", modality="eeg", compute_method="median"):
    # psdIn is a VECTOR, yvect is the basis vector
    if dofeats == "":
        dofeats = DEFAULT_FEAT_ORDER

    if modality == "eeg":
        ch_list = np.arange(0, 257)
    elif modality == "lfp":
        ch_list = ["Left", "Right"]

    feat_vect = []
    for feat in dofeats:
        # print(feat_dict[feat]['param'])
        # dofunc = feat_dict[feat]['fn']
        if compute_method == "median":
            computed_featinspace = FEAT_DICT[feat]["fn"](
                psdIn, yvect, FEAT_DICT[feat]["param"]
            )
        elif compute_method == "mean":
            computed_featinspace = FEAT_DICT[feat]["fn"](
                psdIn, yvect, FEAT_DICT[feat]["param"], cmode=np.mean
            )

        cfis_matrix = [computed_featinspace[ch] for ch in ch_list]
        feat_vect.append(cfis_matrix)
        # feat_dict[feat] = dofunc['fn'](datacontainer,yvect,dofunc['param'])[0]

    feat_vect = np.array(feat_vect).squeeze()

    return feat_vect, dofeats


# Convert a feat dict that comes from a get feature function (WHERE IS IT?!)
def featDict_to_Matr(featDict):
    # structure of feat dict is featDict[FEATURE][CHANNEL] = VALUE
    ret_matr = np.array(
        [
            (featDict[feat]["Left"], featDict[feat]["Right"])
            for feat in DEFAULT_FEAT_ORDER
        ]
    )

    assert ret_matr.shape == (len(DEFAULT_FEAT_ORDER), 2)

    return ret_matr


#%%
# Variables related to what we're soft-coding as our feature library
FEAT_DICT = {
    "Delta": {"fn": get_pow, "param": (1, 4)},
    "Alpha": {"fn": get_pow, "param": (8, 13)},
    "Theta": {"fn": get_pow, "param": (4, 8)},
    "Beta*": {"fn": get_pow, "param": (13, 20)},
    "Beta": {"fn": get_pow, "param": (13, 30)},
    "Gamma1": {"fn": get_pow, "param": (35, 60)},
    "Gamma2": {"fn": get_pow, "param": (60, 100)},
    "Gamma": {"fn": get_pow, "param": (30, 100)},
    "Stim": {"fn": get_pow, "param": (129, 131)},
    "SHarm": {"fn": get_pow, "param": (30, 34)},  # Secondary Harmonic is at 32Hz
    "THarm": {"fn": get_pow, "param": (64, 68)},  # Tertiary Harmonic is at 66Hz!!!
    "Clock": {"fn": get_pow, "param": (104.5, 106.5)},
    "fSlope": {"fn": get_slope, "param": {"frange": (1, 20), "linorder": 1}},
    "nFloor": {"fn": get_slope, "param": {"frange": (50, 200), "linorder": 0}},
    "GCratio": {"fn": get_ratio, "param": ((63, 65), (65, 67))},
}

DEFAULT_FEAT_ORDER = [
    "Delta",
    "Theta",
    "Alpha",
    "Beta*",
    "Gamma1",
]  # ,'fSlope','nFloor']


#%%
# Plotting functions


def plot_TF(TFR, chs=["Left", "Right"]):
    plt.figure()
    for cc, chann in enumerate(chs):
        plt.subplot(1, len(chs), cc + 1)
        aTFR = TFR[chann]

        plt.pcolormesh(aTFR["T"], aTFR["F"], 10 * np.log10(aTFR["SG"]))
        plt.xlabel("Time")
        plt.ylabel("Frequency")


# This function plots the median/mean across time of the TF representation to get the Frequency representation
# Slightly different than doing the Welch directly
def plot_F_fromTF(TFR, chs=["Left", "Right"]):
    plt.figure()
    for cc, chann in enumerate(chs):
        plt.subplot(1, len(chs), cc + 1)
        aTFR = TFR[chann]

        for ss in range(aTFR["SG"].shape[1]):
            plt.plot(aTFR["F"], 10 * np.log10(aTFR["SG"])[:, ss], alpha=0.1)

        plt.plot(aTFR["F"], np.median(10 * np.log10(aTFR["SG"]), axis=1))
        plt.xlabel("Frequency")
        plt.ylabel("Power")


def plot_T(Tser):
    plt.figure()
    for cc, chann in enumerate(Tser.keys()):
        plt.subplot(1, len(Tser.keys()), cc + 1)
        aT = Tser[chann]
        plt.plot(aT["T"], aT["V"])

        plt.xlabel("Time")


def plot_bands(bandM, bandLabels):
    plt.figure()
    for cc in bandM:
        plt.bar(cc)
        plt.xticks(range(len(cc)))
        plt.xticklabels(bandLabels)
