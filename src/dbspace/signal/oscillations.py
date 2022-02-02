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


def get_ratio(Pxx, F, f_r_set, cmode=np.median):
    bandpow = [None] * len(f_r_set)
    # first get the power for each of the individual bands
    for bb, frange in enumerate(f_r_set):
        bandpow[bb] = get_pow(Pxx, F, frange, cmode=cmode)

    ret_ratio = {ch: bandpow[1][ch] / bandpow[0][ch] for ch in bandpow[0].keys()}
    return ret_ratio


def F_Domain(timeser, nperseg=512, noverlap=128, nfft=2**10, Fs=422):

    # assert isinstance(timeser,dbs.timeseries)
    # Window size is about 1 second (512 samples is over 1 sec)

    # what are the dimensions of the timeser we're dealing with?

    Fvect, Pxx = sig.welch(
        timeser,
        Fs,
        window="blackmanharris",
        nperseg=nperseg,
        noverlap=noverlap,
        nfft=nfft,
    )

    FreqReturn = {"F": Fvect, "Pxx": Pxx}

    return FreqReturn


def TF_Domain(timeser, fs=422, nperseg=2**10, noverlap=2**10 - 50):
    # raise Exception
    # assert isinstance(timeser,dbs.timeseries)
    F, T, SG = sig.spectrogram(
        timeser,
        nperseg=nperseg,
        noverlap=noverlap,
        window=sig.get_window("blackmanharris", nperseg),
        fs=fs,
    )

    TFreqReturn = {"T": T, "F": F, "SG": SG}

    return TFreqReturn
