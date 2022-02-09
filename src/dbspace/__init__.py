#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 13:44:38 2018

@author: virati
Main DBSpace module for SCCwm-DBS project
This is the primary file with the methods needed for the rest of the libraries. Used to be called "DBS_Osc" or "DBSOsc".

Copyright (C) 2018 Vineet Ravi Tiruvadi
"""

print("Using DBSpace LATEST")
import numpy as np
import pandas as pd
from collections import defaultdict
import scipy.signal as sig
from .signal.oscillations import get_pow

# IF you want to do OR related analyses, this needs to be uncommented
# from brpylib import NsxFile


import matplotlib.pyplot as plt

plt.rcParams["image.cmap"] = "jet"

np.seterr(divide="raise")

all_pts = ["901", "903", "905", "906", "907", "908"]

# This is our map for the electrodes that each patient has for ONTarget and OFFTarget
Etrode_map = {
    "OnT": {
        "901": (2, 1),
        "903": (2, 2),
        "905": (2, 1),
        "906": (2, 2),
        "907": (1, 1),
        "908": (2, 1),
    },
    "OffT": {
        "901": (1, 2),
        "903": (1, 1),
        "905": (1, 2),
        "906": (1, 1),
        "907": (2, 2),
        "908": (1, 2),
    },
}
# TODO This needs to be a config/yaml/whatever file in the chapter-specific repos

#%%
# BlackRock Methods
def load_or_file(fname, **kwargs):
    # nsx_file = NsxFile(fname)

    arg_list = ["elec_ids", "start_time_s", "data_time_s", "downsample", "plot_chan"]

    for aa in arg_list:
        print(arg_list[aa])

    # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
    cont_data = nsx_file.getdata(elec_ids, start_time_s, data_time_s, downsample)

    # Close the nsx file now that all data is out
    nsx_file.close()

    return cont_data


#%%
# BRAIN RADIO METHODS
# Method to load in brainradio file


def gen_T(inpX, Fs=422, nfft=2**10):
    outT = defaultdict(dict)
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
    outPSD = defaultdict(dict)
    outPoly = defaultdict(dict)
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


# This function takes a PSD and subtracts out the PSD's fourth order polynomial fit


"""Throw an error if we're calling the old name for the vector function below"""


def poly_subtrLFP(**kwargs):
    raise Exception


"""Below used to be called poly_subtrLFP, unclear whether it was being used, now renamed and look for errors elsewhere"""


def poly_SG(inSG, fVect, order=4):
    out_sg = np.zeros_like(inSG)

    for seg in range(inSG.shape[1]):
        inpsd = 10 * np.log10(inpPSD[chann][seg, :])
        polyCoeff = np.polyfit(fVect, inpsd, order)
        polyfunc = np.poly1d(polyCoeff)
        polyitself = polyfunc(fVect)
        out_sg[:, seg] = 10 ** ((curr_psd - polyitself) / 10)

    return out_sg


def gen_SG(inpX, Fs=422, nfft=2**10, plot=False, overlap=True):
    outSG = defaultdict(dict)
    for chann in inpX.keys():
        if overlap == True:
            outSG[chann] = TF_Domain(inpX[chann])
        else:
            outSG[chann] = TF_Domain(inpX[chann], noverlap=0, nperseg=422 * 2)

    if plot:
        plot_TF(outSG, chs=inpX.keys())

    return outSG


""" Return to us the power in an oscillatory feature"""

"""
F and TF domain plotting 
"""


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


""" Return the phases and the labels associated with them"""
all_phases = (
    ["A0" + str(num) for num in range(4, 0, -1)]
    + ["B0" + str(num) for num in range(1, 5)]
    + ["C0" + str(num) for num in range(1, 10)]
    + ["C" + str(num) for num in range(10, 25)]
)


def Phase_List(exprs="all"):
    if exprs == "all":
        return all_phases
    elif exprs == "ephys":
        return all_phases[4:]
    elif exprs == "therapy":
        return all_phases[8:]
    elif exprs == "notherapy":
        return all_phases[0:8]


""" Function to plot bands since this has been annoying every time I've had to recode the thing from scratch"""


def plot_bands(bandM, bandLabels):
    plt.figure()
    for cc in bandM:
        plt.bar(cc)
        plt.xticks(range(len(cc)))
        plt.xticklabels(bandLabels)
