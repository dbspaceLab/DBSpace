import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
from dbspace.utils.costs import l2_pow

# Function to go through and find all the features from the PSD structure of dbo
def calc_feats(psdIn, yvect, dofeats="", modality="eeg", compute_method="median"):
    # psdIn is a VECTOR, yvect is the basis vector
    if dofeats == "":
        dofeats = feat_order

    if modality == "eeg":
        ch_list = np.arange(0, 257)
    elif modality == "lfp":
        ch_list = ["Left", "Right"]

    feat_vect = []
    for feat in dofeats:
        # print(feat_dict[feat]['param'])
        # dofunc = feat_dict[feat]['fn']
        if compute_method == "median":
            computed_featinspace = feat_dict[feat]["fn"](
                psdIn, yvect, feat_dict[feat]["param"]
            )
        elif compute_method == "mean":
            computed_featinspace = feat_dict[feat]["fn"](
                psdIn, yvect, feat_dict[feat]["param"], cmode=np.mean
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
        [(featDict[feat]["Left"], featDict[feat]["Right"]) for feat in feat_order]
    )

    # assert that the size is as expected?
    # should be number of feats x number of channels!
    assert ret_matr.shape == (len(feat_order), 2)

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
            fixed_psd, polyitself = dbo.poly_subtr(corr_psd, F)
            pf_lPSD[chann_label[cc]] = fixed_psd[chann_label[cc]].reshape(-1, 1)
        else:
            correct_psd, polyitself = dbo.poly_subtr(corr_psd, F)

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

    fcalced, bands = dbo.calc_feats(
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
        dofeats = feat_order

    if modality == "eeg":
        ch_list = np.arange(0, 257)
    elif modality == "lfp":
        ch_list = ["Left", "Right"]

    feat_vect = []
    for feat in dofeats:
        # print(feat_dict[feat]['param'])
        # dofunc = feat_dict[feat]['fn']
        if compute_method == "median":
            computed_featinspace = feat_dict[feat]["fn"](
                psdIn, yvect, feat_dict[feat]["param"]
            )
        elif compute_method == "mean":
            computed_featinspace = feat_dict[feat]["fn"](
                psdIn, yvect, feat_dict[feat]["param"], cmode=np.mean
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
        [(featDict[feat]["Left"], featDict[feat]["Right"]) for feat in feat_order]
    )

    # assert that the size is as expected?
    # should be number of feats x number of channels!
    assert ret_matr.shape == (len(feat_order), 2)

    return ret_matr


""" Higher order measures here"""
# Make a coherence generation function
def gen_coher(inpX, Fs=422, nfft=2**10, polyord=0):
    print("Starting a coherence run...")
    outPLV = nestdict()
    outCSD = nestdict()

    fvect = np.linspace(0, Fs / 2, nfft / 2 + 1)

    for chann_i in inpX.keys():
        print(chann_i)
        for chann_j in inpX.keys():
            csd_ensemble = np.zeros(
                (inpX[chann_i].shape[1], len(feat_order)), dtype=complex
            )
            plv = np.zeros((inpX[chann_i].shape[1], len(feat_order)))

            for seg in range(inpX[chann_i].shape[1]):
                # First we get the cross spectral density
                csd_out = sig.csd(
                    inpX[chann_i][:, seg], inpX[chann_j][:, seg], fs=Fs, nperseg=512
                )[1]

                # normalize the entire CSD for the total power in input signals
                norm_ms_csd = np.abs(csd_out) / np.sqrt(
                    l2_pow(inpX[chann_i][:, seg]) * l2_pow(inpX[chann_j][:, seg])
                )

                # Are we focusing on a band or doing the entire CSD?

                for bb, band in enumerate(feat_order):
                    # what are our bounds?
                    band_bounds = feat_dict[band]["param"]
                    band_idxs = np.where(
                        np.logical_and(fvect >= band_bounds[0], fvect <= band_bounds[1])
                    )
                    csd_ensemble[seg, bb] = cmedian(csd_out[band_idxs])
                    plv[seg, bb] = np.max(norm_ms_csd[band_idxs])

                # Below brings in the entire csd, but this is dumb
                # csd_ensemble[seg] = csd_out

                # Compute the PLV

            # Here we find the median across segments
            # outCSD[chann_i][chann_j] = cmedian(csd_ensemble,axis=0)
            outCSD[chann_i][chann_j] = csd_ensemble

            # Compute the normalized coherence/PLV
            outPLV[chann_i][chann_j] = plv
            # outPLV[chann_i][chann_j] = np.median(plv,axis=0)

            ## PLV abs EEG ->
            ## Coherence value

    return outCSD, outPLV


#%%
# Variables related to what we're soft-coding as our feature library
feat_dict = {
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

feat_order = ["Delta", "Theta", "Alpha", "Beta*", "Gamma1"]  # ,'fSlope','nFloor']
