#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 21:02:09 2018

@author: virati
Main Class for Processed dEEG Data

"""

import random
import sys
from collections import defaultdict

import dbspace as dbo
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import mne
import numpy as np
import scipy.signal as sig
import scipy.stats as stats

import sklearn
from dbspace.signal.oscillations import calc_feats, FEAT_DICT, DEFAULT_FEAT_ORDER, gen_psd
from dbspace.utils.structures import nestdict
from dbspace.viz.MM import EEG_Viz
from scipy.io import loadmat
from sklearn import mixture, svm
from sklearn.decomposition import PCA, FastICA
from sklearn.metrics import auc, confusion_matrix, roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, learning_curve
from sklearn.utils import resample
from statsmodels import robust
import json
import logging
sys.path.append("/home/virati/Dropbox/projects/libs/robust-pca/")
import r_pca

#%%
#%%

class network_action_dEEG:
    keys_of_interest = {"OnT": ["Off_3", "BONT"], "OffT": ["Off_3", "BOFT"]}

    def __init__(
        self,
        pts,
        config_file=None,
        procsteps="liberal",
        condits=["OnT", "OffT"],
        pretty_mode=False,
        polyfix=0,
    ):

        if config_file is None:
            raise ValueError(
                "No config file provided, please provide an experiment json..."
            )

        self.load_config(config_file)

        self.chann_dim = 257
        self.ch_order_list = range(self.chann_dim)
        self.procsteps = procsteps

        self.do_pts = pts
        self.condits = condits

        self.polyorder = polyfix
        self.pretty = pretty_mode

        self.label_map = {"OnT": 1, "OffT": 0}

        # Load in the data
        self.ts_data = self.load_data(pts)

        self.eeg_locs = mne.channels.read_montage(
            "/home/virati/Dropbox/GSN-HydroCel-257.sfp"
        )

        self.gen_output_variables()

    """Setup all the output variables we need"""

    def load_config(self, config_file):
        with open(config_file, "r") as config:
            Targeting = json.load(config)

        self.targeting_config = Targeting

    def gen_output_variables(self):
        # CHECK IF we're still using ANY of these

        # sloppy containers for the outputs of our analyses
        self.psd_trans = {
            pt: {condit: {epoch: [] for epoch in self.keys_of_interest} for condit in self.condits}
            for pt in self.do_pts
        }
        self.PSD_diff = {
            pt: {condit: [] for condit in self.condits} for pt in self.do_pts
        }
        self.PSD_var = {
            pt: {condit: [] for condit in self.condits} for pt in self.do_pts
        }

        self.Feat_trans = {
            pt: {condit: {epoch: [] for epoch in self.keys_of_interest} for condit in self.condits}
            for pt in self.do_pts
        }
        self.Feat_diff = {
            pt: {condit: [] for condit in self.condits} for pt in self.do_pts
        }
        self.Feat_var = {
            pt: {condit: [] for condit in self.condits} for pt in self.do_pts
        }

    """Run the standard pipeline to prepare segments"""

    def standard_pipeline(self):
        print("Doing standard init pipeline")
        self.extract_feats(polyorder=self.polyorder)
        self.pool_patients()  # pool all the DBS RESPONSE vectors
        # self.median_responses = self.median_response(pt=self.do_pts)
        self.median_responses = self.median_bootstrap_response(
            pt="POOL", bootstrap=100
        )["mean"]

    """Load in the MAT data for preprocessed EEG recordings"""

    def load_data(self, pts):
        ts_data = defaultdict(dict)
        for pt in pts:
            ts_data[pt] = defaultdict(dict)
            for condit in self.condits:
                ts_data[pt][condit] = defaultdict(dict)

                temp_data = loadmat(self.targeting_config[self.procsteps][pt][condit])

                for epoch in self.keys_of_interest[condit]:
                    ts_data[pt][condit][epoch] = temp_data[epoch]

        self.fs = temp_data["EEGSamplingRate"][0][0]
        self.donfft = 2**11
        self.fvect = np.linspace(
            0,
            np.round(self.fs / 2).astype(int),
            np.round(self.donfft / 2 + 1).astype(int),
        )

        return ts_data

    """Extract features from the EEG datasegments"""

    def extract_feats(self, polyorder=4):
        pts = self.do_pts

        psd_dict = nestdict()
        osc_dict = nestdict()

        for pt in pts:
            # feat_dict[pt] = defaultdict(dict)

            for condit in self.condits:
                # feat_dict[pt][condit] = defaultdict(dict)
                for epoch in self.keys_of_interest[condit]:
                    # find the mean for all segments
                    data_matr = self.ts_data[pt][condit][
                        epoch
                    ]  # should give us a big matrix with all the crap we care about
                    data_dict = {
                        ch: data_matr[ch, :, :].squeeze()
                        for ch in range(data_matr.shape[0])
                    }  # transpose is done to make it seg x time

                    # TODO check if this is in the right units
                    seg_psds = gen_psd(
                        data_dict, Fs=self.fs, nfft=self.donfft, polyord=polyorder
                    )

                    # gotta flatten the DICTIONARY, so have to do it carefully
                    PSD_matr = np.array([seg_psds[ch] for ch in self.ch_order_list])

                    OSC_matr = np.zeros((seg_psds[0].shape[0], 257, len(DEFAULT_FEAT_ORDER)))
                    # middle_osc = {chann:seg_psd for chann,seg_psd in seg_psds.items}
                    middle_osc = np.array([seg_psds[ch] for ch in range(257)])

                    # have to go to each segment due to code
                    for ss in range(seg_psds[0].shape[0]):
                        try:
                            state_return = calc_feats(middle_osc[:, ss, :], self.fvect)[
                                0
                            ].T
                            state_return[
                                :, 4
                            ] = 0  # we know gamma is nothing for this entire analysis
                            # pdb.set_trace()
                            OSC_matr[ss, :, :] = np.array(
                                [state_return[ch] for ch in range(257)]
                            )
                        except Exception as e:
                            raise ValueError("Something going wrong here...")

                    # find the variance for all segments
                    psd_dict[pt][condit][epoch] = PSD_matr
                    osc_dict[pt][condit][epoch] = OSC_matr

                    # need to do OSCILLATIONS here

        # THIS IS THE PSDs RAW, not log transformed
        self.psd_dict = psd_dict
        # Below is the oscillatory power
        # we should go through each osc_dict element and zero out the gamma
        self.osc_dict = osc_dict

    """Generate Distributions across all channels for each band"""

    def band_distrs(self, pt="POOL"):
        print("Plotting Distribution for Bands")

        meds = nestdict()
        mads = nestdict()

        marker = ["o", "s"]
        color = ["b", "g"]
        plt.figure()
        ax2 = plt.subplot(111)
        for cc, condit in enumerate(["OnT", "OffT"]):
            allsegs = np.median(self.osc_bl_norm[pt][condit][:, :], axis=0).squeeze()
            # pdb.set_trace()
            parts = ax2.violinplot(allsegs, positions=np.arange(5) + (cc - 0.5) / 10)
            for pc in parts["bodies"]:
                pc.set_facecolor(color[cc])
                pc.set_edgecolor(color[cc])
                # pc.set_linecolor(color[cc])

            # plt.ylim((-0.5,0.5))

        for bb in range(5):
            pass
            # rsres = stats.ranksums(meds['OnT'][:,bb],meds['OffT'][:,bb])
            # rsres = stats.wilcoxon(meds['OnT'][:,bb],meds['OffT'][:,bb])
            # rsres = stats.ttest_ind(10**(meds['OnT'][:,bb]/10),10**(meds['OffT'][:,bb]/10))
            # print(rsres)

        # plt.suptitle(condit)

    """ This function sets the response vectors to the targets x patient"""

    def compute_response(self, do_pts=[], condits=["OnT", "OffT"]):
        if do_pts == []:
            do_pts = self.do_pts

        BL = {pt: {condit: [] for condit in condits} for pt in do_pts}
        response = {pt: {condit: [] for condit in condits} for pt in do_pts}

        for pt in do_pts:
            for condit in condits:
                # first, compute the median state during baseline
                try:
                    BL[pt][condit] = np.median(
                        self.osc_dict[pt][condit]["Off_3"], axis=0
                    )
                except:
                    pdb.set_trace()
                # Now, go to each segment during stim and subtract the BL for that
                response[pt][condit] = (
                    self.osc_dict[pt][condit][self.keys_of_interest[condit][1]] - BL[pt][condit]
                )

        self.targ_response = response

    def response_stats(self, band="Alpha", plot=False):
        band_idx = DEFAULT_FEAT_ORDER.index(band)
        response_diff_stats = {pt: [] for pt in self.do_pts}

        ## First, check to see if per-channel h-testing rejects the null
        for pt in self.do_pts:
            for cc in range(256):
                response_diff_stats[pt].append(
                    stats.mannwhitneyu(
                        self.targ_response[pt]["OnT"][:, cc, band_idx],
                        self.targ_response[pt]["OffT"][:, cc, band_idx],
                    )[1]
                )

        self.response_diff_stats = response_diff_stats

        ## Now check variances\
        ONT_var = {pt: [] for pt in self.do_pts}
        OFFT_var = {pt: [] for pt in self.do_pts}
        pool_ONT = []
        pool_OFFT = []
        for pt in self.do_pts:
            for cc in range(256):
                ONT_var[pt].append(
                    np.var(self.targ_response[pt]["OnT"][:, cc, band_idx])
                )
                OFFT_var[pt].append(
                    np.var(self.targ_response[pt]["OffT"][:, cc, band_idx])
                )

            pool_ONT.append(self.targ_response[pt]["OnT"][:, :, band_idx])
            pool_OFFT.append(self.targ_response[pt]["OffT"][:, :, band_idx])

        # Now stack across all patients
        pool_ONT_var = np.var(np.concatenate(pool_ONT, axis=0), axis=0)
        pool_OFFT_var = np.var(np.concatenate(pool_OFFT, axis=0), axis=0)

        ch_response_sig = {pt: np.array(response_diff_stats[pt]) for pt in self.do_pts}
        aggr_resp_sig = np.array(
            [(resp < 0.05 / 256).astype(int) for pt, resp in ch_response_sig.items()]
        )
        union_sig = np.sum(aggr_resp_sig, axis=0) >= 2

        if plot:
            for pt in self.do_pts:
                if 0:
                    pass

                    # Look at each patient's ONT and OFFT VARIANCE
                    bins = np.linspace(0, 40, 100)
                    plt.figure()
                    plt.violinplot(ONT_var[pt])  # ,bins=bins)
                    print(np.median(ONT_var[pt]))
                    plt.violinplot(OFFT_var[pt])  # ,bins=bins)
                    print(np.median(OFFT_var[pt]))

                # Stats for ONT vs OFFT within each patient
                plt.figure()
                plt.plot(response_diff_stats[pt])
                plt.hlines(0.05 / 256, 0, 256)
                n_sig = np.sum((ch_response_sig[pt] < 0.05 / 256).astype(int))
                plt.suptitle(pt + " " + str(n_sig))

            plt.figure()
            plt.plot(pool_ONT_var)
            plt.plot(pool_OFFT_var)
            print(np.median(pool_ONT_var))
            print(np.median(pool_OFFT_var))
            plt.suptitle("Pooled stats for ONT/OFFT consistency check")

    """Generate the pooled/ensemble segment response matrices -> dict"""

    def pool_patients(self):
        print("Pooling Patient Observations")
        self.osc_bl_norm = {
            pt: {
                condit: self.osc_dict[pt][condit][self.keys_of_interest[condit][1]]
                - np.median(self.osc_dict[pt][condit][self.keys_of_interest[condit][0]], axis=0)
                for condit in self.condits
            }
            for pt in self.do_pts
        }
        self.osc_bl_norm["POOL"] = {
            condit: np.concatenate(
                [
                    self.osc_dict[pt][condit][self.keys_of_interest[condit][1]]
                    - np.median(self.osc_dict[pt][condit][self.keys_of_interest[condit][0]], axis=0)
                    for pt in self.do_pts
                ]
            )
            for condit in self.condits
        }

        self.osc_stim = nestdict()
        self.osc_stim = {
            pt: {
                condit: 10 ** (self.osc_dict[pt][condit][self.keys_of_interest[condit][1]] / 10)
                for condit in self.condits
            }
            for pt in self.do_pts
        }
        self.osc_stim["POOL"] = {
            condit: np.concatenate(
                [
                    10 ** (self.osc_dict[pt][condit][self.keys_of_interest[condit][1]] / 10)
                    for pt in self.do_pts
                ]
            )
            for condit in self.condits
        }

    def plot_psd(self, pt, condit, epoch):
        plt.figure()
        # reshape
        reshaped = self.psd_dict[pt][condit][epoch].reshape(-1, 1025)
        print(reshaped.shape)
        plt.plot(
            np.linspace(0, 1000 / 2, 1025),
            np.mean(20 * np.log10(reshaped.T), axis=1),
            alpha=0.9,
        )
        plt.plot(
            np.linspace(0, 1000 / 2, 1025),
            20 * np.log10(reshaped.T)[:, np.random.randint(0, 256, size=(30,))],
            alpha=0.2,
        )
        plt.xlim((0, 300))

    # Median dimensionality reduction here; for now rPCA
    def distr_response(self, pt="POOL"):
        return {condit: self.osc_bl_norm[pt][condit] for condit in self.condits}

    """Compute the median response using bootstrap"""

    def median_bootstrap_response(self, pt="POOL", mfunc=np.mean, bootstrap=100):
        print("Computing Bootstrap Median Response for " + pt)

        bs_mean = []
        bs_var = []
        for ii in range(bootstrap):
            rnd_idxs = {
                condit: random.sample(range(self.osc_bl_norm[pt][condit].shape[0]), 100)
                for condit in self.condits
            }
            bs_mean.append(
                {
                    condit: mfunc(
                        self.osc_bl_norm[pt][condit][rnd_idxs[condit], :, :], axis=0
                    )
                    for condit in self.condits
                }
            )
            # bs_var.append({condit:np.var(self.osc_bl_norm[pt][condit][rnd_idxs[condit],:,:],axis=0) for condit in self.condits})

        mean_of_means = {
            condit: np.mean([iteration[condit] for iteration in bs_mean], axis=0)
            for condit in self.condits
        }
        var_of_means = {
            condit: np.var([iteration[condit] for iteration in bs_mean], axis=0)
            for condit in self.condits
        }

        return {"mean": mean_of_means, "var": var_of_means}

    # Do per-channel, standard stats. Compare pre-stim to stim condition
    def per_chann_stats(self, condit="OnT", band="Alpha"):
        band_idx = DEFAULT_FEAT_ORDER.index(band)

        for pt in self.do_pts:
            ch_stat = np.zeros((257,))
            ch_bl_mean = []
            ch_stim_mean = []
            for ch in range(256):
                # distribution for pre-stimulation period
                # pdb.set_trace()
                baseline_distr = []
                stim_distr = []
                for ii in range(100):
                    bl_rand_idx = random.sample(
                        range(
                            0, self.osc_dict[pt][condit][self.keys_of_interest[condit][0]].shape[0]
                        ),
                        10,
                    )
                    stim_rand_idx = random.sample(
                        range(
                            0, self.osc_dict[pt][condit][self.keys_of_interest[condit][1]].shape[0]
                        ),
                        10,
                    )

                    baseline_distr.append(
                        np.mean(
                            self.osc_dict[pt][condit][self.keys_of_interest[condit][0]][
                                bl_rand_idx, ch, band_idx
                            ]
                        )
                    )
                    stim_distr.append(
                        np.mean(
                            self.osc_dict[pt][condit][self.keys_of_interest[condit][1]][
                                stim_rand_idx, ch, band_idx
                            ]
                        )
                    )

                # baseline_distr = self.osc_dict[pt][condit][self.keys_of_interest[condit][0]][0:20,ch,band_idx]#should be segments x bands
                # stim_distr = self.osc_dict[pt][condit][self.keys_of_interest[condit][1]][0:20,ch,band_idx]
                diff_stat = stats.mannwhitneyu(baseline_distr, stim_distr)
                # diff_stat = stats.f_oneway(baseline_distr,stim_distr)
                print(str(ch) + ":" + str(diff_stat))
                ch_stat[ch] = diff_stat[1]

                # plt.violinplot(baseline_distr)
                # plt.violinplot(stim_distr)
                ch_bl_mean.append(np.mean(baseline_distr))
                ch_stim_mean.append(np.mean(stim_distr))

            plt.figure()
            plt.plot(ch_stat)
            plt.axhline(0.05 / 256, 0, 256)

            plt.figure()
            plt.violinplot(ch_bl_mean)
            plt.violinplot(ch_stim_mean)

    def topo_median_variability(
        self, pt="POOL", band="Alpha", do_condits=[], use_maya=False
    ):
        band_i = DEFAULT_FEAT_ORDER.index(band)

        # medians = self.median_response(pt=pt)

        for condit in do_condits:
            response_dict = np.median(
                self.osc_bl_norm[pt][condit][:, :, :], axis=0
            ).squeeze()

            var_dict = robust.mad(
                self.osc_bl_norm[pt][condit][:, :, :], axis=0
            ).squeeze()
            var_mask = ((var_dict) > 2.5).astype(int)
            # The old scatterplot approach
            bins = np.linspace(0, 5, 50)
            plt.figure()
            plt.hist(var_dict[:, band_i], bins=bins)
            plt.ylim((0, 70))
            if use_maya:
                EEG_Viz.maya_band_display(var_dict[:, band_i])
            else:
                # EEG_Viz.plot_3d_scalp(response_mask[:,band_i]*var_dict[:,band_i],plt.figure(),label=condit + ' Response Var ' + band + ' | ' + pt,unwrap=True,scale=100,clims=(0,4),alpha=0.3,marker_scale=10)
                EEG_Viz.plot_3d_scalp(
                    var_mask[:, band_i],
                    plt.figure(),
                    label=condit + " Response Var " + band + " | " + pt,
                    unwrap=True,
                    scale=100,
                    clims=(0, 4),
                    alpha=0.3,
                    marker_scale=10,
                    anno_top=False,
                    binary_mask=True,
                )
                plt.suptitle(pt)

    def topo_median_response(
        self,
        pt="POOL",
        band="Alpha",
        do_condits=[],
        use_maya=False,
        scale_w_mad=False,
        avg_func=np.median,
    ):
        band_i = DEFAULT_FEAT_ORDER.index(band)
        # medians = self.median_response(pt=pt)

        for condit in do_condits:
            response_dict = avg_func(
                self.osc_bl_norm[pt][condit][:, :, :], axis=0
            ).squeeze()
            # The old scatterplot approach
            if use_maya:
                EEG_Viz.maya_band_display(response_dict[:, band_i])
            else:
                if scale_w_mad:
                    mad_scale = robust.mad(
                        self.osc_bl_norm[pt][condit][:, :, :], axis=0
                    ).squeeze()

                    EEG_Viz.plot_3d_scalp(
                        response_dict[:, band_i],
                        plt.figure(),
                        label=condit + " Mean Response " + band + " | " + pt,
                        unwrap=True,
                        scale=100,
                        clims=(-1, 1),
                        alpha=0.3,
                        marker_scale=10 * np.tanh(mad_scale[:, band_i] - 5),
                    )

                else:
                    EEG_Viz.plot_3d_scalp(
                        response_dict[:, band_i],
                        plt.figure(),
                        label=condit + " Mean Response " + band + " | " + pt,
                        unwrap=True,
                        scale=100,
                        clims=(-1, 1),
                        alpha=0.3,
                        marker_scale=5,
                    )
                plt.suptitle(pt)

    def calc_median_response(self, pt="POOL", mfunc=np.median):
        print("Computing Median Response for " + pt)
        print("Doing " + str(mfunc))
        return {
            condit: mfunc(self.osc_bl_norm[pt][condit], axis=0)
            for condit in self.condits
        }

    def support_analysis(
        self, support_struct, pt="POOL", condit="OnT", voltage="3", band="Alpha"
    ):
        # support_struct = pickle.load(open('/tmp/'+ pt + '_' + condit + '_' + voltage,'rb'))
        if band == "rP0":
            medians = self.dyn_L.swapaxes(
                0, 1
            )  # if we want to use the 0th component of the dyn_rPCA eigenvector
            band_i = 0
        else:
            medians = self.calc_median_response(pt=pt)[
                condit
            ]  # if we want to use the standard median Alpha change
            band_i = DEFAULT_FEAT_ORDER.index(band)

        # medians = np.median(self.targ_response[pt][condit],axis=0)
        fig = plt.figure()
        # First, we'll plot what the medians actually are

        EEG_Viz.plot_3d_scalp(
            medians[:, band_i],
            fig,
            label=condit + " Mean Response " + band,
            unwrap=True,
            scale=10,
        )
        plt.suptitle(pt)

        full_distr = medians[
            :, band_i
        ]  # - np.mean(medians[:,band_i]) #this zeros the means of the distribution

        primary_distr = full_distr[support_struct["primary"] == 1]
        # now we'll circle where the primary nodes are

        print(np.sum((support_struct["primary"] == 1).astype(int)))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        EEG_Viz.plot_3d_scalp(
            support_struct["primary"], ax, scale=10, alpha=0.5, unwrap=True
        )
        plt.title("Primary Channels")

        secondary_distr = full_distr[support_struct["secondary"] == 1]
        print(np.sum((support_struct["secondary"] == 1).astype(int)))
        fig = plt.figure()
        EEG_Viz.plot_3d_scalp(
            support_struct["secondary"], fig, scale=10, alpha=0.5, unwrap=True
        )
        plt.title("Secondary Channels")

        labels = []

        def add_label(violin, label):
            color = violin["bodies"][0].get_facecolor().flatten()
            labels.append((mpatches.Patch(color=color), label))

        plt.figure()
        bins = np.linspace(-2, 2, 20)
        # plt.hist(primary_distr,bins=bins,alpha=0.5,label='Primary')
        print("Primary mean: " + str(np.median(primary_distr)))
        add_label(plt.violinplot(primary_distr), "Primary Nodes")
        # pdb.set_trace()

        # plt.hist(secondary_distr,bins=bins,alpha=0.5,label='Secondary')
        print("Secondary mean: " + str(np.median(secondary_distr)))
        add_label(plt.violinplot(secondary_distr), "Secondary Nodes")
        plt.legend(*zip(*labels), loc=2)

        print(stats.ks_2samp(primary_distr, secondary_distr))

        # plt.hist(full_distr,bins=bins,alpha=0.5,label='FULL')
        # plt.legend(['Primary','','','Secondary'])
        plt.title(pt + " " + condit + " " + band)

    """I guess this is about developing a rPCA approach to *dynamic* response without oscillations?"""

    def OnT_ctrl_dyn(self, pt="POOL", condit="OnT", do_plot=False):
        source_label = "Dyn PCA"

        response_stack = self.osc_bl_norm["POOL"][condit][:, :, 2]
        # Focusing just on alpha
        # response_stack = np.dot(response_stack.T,response_stack)

        # pdb.set_trace()
        rpca = r_pca.R_pca(response_stack)
        L, S = rpca.fit()

        svm_pca = PCA()
        svm_pca.fit(L)

        # pdb.set_trace()
        svm_pca_coeffs = svm_pca.components_
        # ALL PLOTTING BELOW
        if do_plot:
            for comp in range(5):
                fig = plt.figure()
                EEG_Viz.plot_3d_scalp(
                    (L[comp, :]),
                    fig,
                    label="OnT Mean Response",
                    unwrap=True,
                    scale=100,
                    alpha=0.3,
                    marker_scale=5,
                )
                plt.title("rPCA Component " + str(comp))

            plt.figure()
            plt.subplot(221)
            plt.plot(svm_pca.explained_variance_ratio_)
            plt.ylim((0, 1))
            plt.subplot(222)
            plt.plot(np.array(svm_pca_coeffs)[0:4])
            plt.legend(["PC1", "PC2", "PC3", "PC4", "PC5"])
            plt.title("rPCA Components " + source_label)

        self.dyn_pca = svm_pca
        self.dyn_L = L

    def OnT_ctrl_modes_segs_ICA(self, pt="POOL", do_plot=False):

        seg_responses = self.osc_bl_norm[pt]["OnT"][:, :, 0:4]
        source_label = "Segment Responses"

        svm_ica_coeffs = []
        for ii in range(seg_responses.shape[0]):

            # pdb.set_trace()
            rpca = r_pca.R_pca(seg_responses[ii, :, :])
            L, S = rpca.fit()

            # L = seg_responses[ii,:,:]
            # S = 0
            svm_ica = FastICA(tol=0.1, max_iter=1000)

            svm_ica.fit(L)
            rotated_L = svm_ica.fit_transform(L)

            svm_ica_coeffs.append(svm_ica.components_)

        mode_model = {
            "L": L,
            "S": S,
            "Vectors": svm_ica_coeffs,
            "Model": svm_ica,
            "RotatedL": rotated_L,
        }
        return mode_model

    def OnT_alpha_modes_segs(
        self, pt="POOL", data_source=[], do_plot=False, band="Alpha"
    ):
        seg_responses = self.osc_bl_norm[pt]["OnT"][
            :, :, DEFAULT_FEAT_ORDER.index(band)
        ].squeeze()
        source_label = "Alpha Segmental Response"

        lr_pca_coeffs = []
        S_pca_coeffs = []
        rot_L = []
        rot_S = []

        # pdb.set_trace()
        rpca = r_pca.R_pca(seg_responses[:, :].T)
        L, S = rpca.fit()

        # L = seg_responses[ii,:,:]
        # S = 0
        lr_pca = PCA()

        lr_pca.fit(L)
        rotated_L = lr_pca.fit_transform(L)
        rot_L.append(rotated_L)

        lr_pca_coeffs.append(lr_pca.components_)

        S_pca = PCA()
        S_pca.fit(S)
        rotated_S = S_pca.fit_transform(S)
        rot_S.append(rotated_S)
        S_pca_coeffs.append(S_pca.components_)

        mode_model = {
            "L": L,
            "S": S,
            "SModel": S_pca,
            "RotatedS": np.array(rot_S),
            "SVectors": S_pca_coeffs,
            "Vectors": lr_pca_coeffs,
            "Model": lr_pca,
            "RotatedL": np.array(rot_L),
        }
        return mode_model

    def topo_OnT_alpha_ctrl(self, **kwargs):
        self.alpha_ctrl_model = self.OnT_alpha_modes_segs(**kwargs)
        # model = self.alpha_ctrl_model

    def plot_alpha_ctrl_L(self, top_comp=5):
        model = self.alpha_ctrl_model

        L = model["RotatedL"]
        expl_var = model["Model"].explained_variance_ratio_
        coeffs = np.array(model["Vectors"])

        # ALL PLOTTING BELOW
        # Plot the topo for our low-rank component
        for comp in range(0, top_comp):
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(
                np.median(L, axis=0)[:, comp],
                fig,
                label="OnT Mean Response",
                unwrap=True,
                scale=100,
                alpha=0.3,
                marker_scale=5,
            )
            plt.title("rPCA Component " + str(comp) + " " + str(expl_var[comp]))

        plt.figure()
        plt.subplot(221)
        plt.plot(expl_var)
        plt.ylim((0, 1))
        plt.subplot(222)
        for ii in range(4):  # this loops through our COMPONENTS to find the end
            plt.plot(
                np.mean(coeffs, axis=0)[ii, :].T, linewidth=5 - ii, alpha=expl_var[ii]
            )
        plt.ylim((-0.3, 0.3))
        # plt.hlines(0,0,5)
        plt.legend(["PC0", "PC1", "PC2", "PC3", "PC4"])
        plt.title("rPCA Components")

    def plot_alpha_ctrl_S(self, top_comp=10):
        model = self.alpha_ctrl_model

        S = model["RotatedS"]
        s_coeffs = np.array(model["SVectors"])
        s_expl_var = model["SModel"].explained_variance_ratio_

        for comp in range(0, top_comp):
            # Plot sparse component next
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(
                np.median(S, axis=0)[:, comp],
                fig,
                label="OnT Mean Response",
                unwrap=True,
                scale=100,
                alpha=0.3,
                marker_scale=5,
            )
            plt.title("Sparse Component " + str(comp))

        plt.figure()
        plt.subplot(221)
        plt.plot(s_expl_var)
        plt.ylim((0, 1))
        plt.subplot(222)

        for ii in range(4):  # this loops through our COMPONENTS to find the end
            plt.plot(
                np.mean(s_coeffs, axis=0)[ii, :].T, linewidth=5 - ii
            )  # ,alpha=expl_var[ii])
        plt.ylim((-0.3, 0.3))
        plt.hlines(0, 0, 5)
        plt.legend(["PC1", "PC2", "PC3", "PC4", "PC5"])
        plt.title("rPCA Components")

    # def OnT_ctrl_modes_segs(self,pt='POOL',data_source=[],do_plot=False):

    """ We do naive rPCA on each segment and then average things together: this has been superceded by the tensor decomposition"""

    def OnT_ctrl_modes_segs(self, **kwargs):
        pt = kwargs["pt"]
        do_plot = kwargs["do_plot"]

        print("Using BL Norm Segments - RAW")
        seg_responses = self.osc_bl_norm[pt]["OnT"][:, :, 0:4]
        source_label = "Segment Responses"

        lr_pca_coeffs = []
        S_pca_coeffs = []
        rot_L = []
        rot_S = []
        for ii in range(seg_responses.shape[0]):

            # pdb.set_trace()
            rpca = r_pca.R_pca(seg_responses[ii, :, :])
            L, S = rpca.fit()

            # L = seg_responses[ii,:,:]
            # S = 0
            lr_pca = PCA()

            lr_pca.fit(L)
            rotated_L = lr_pca.fit_transform(L)
            rot_L.append(rotated_L)

            lr_pca_coeffs.append(lr_pca.components_)

            S_pca = PCA()
            S_pca.fit(S)
            rotated_S = S_pca.fit_transform(S)
            rot_S.append(rotated_S)
            S_pca_coeffs.append(S_pca.components_)

        mode_model = {
            "L": L,
            "S": S,
            "SModel": S_pca,
            "RotatedS": np.array(rot_S),
            "SVectors": S_pca_coeffs,
            "Vectors": lr_pca_coeffs,
            "Model": lr_pca,
            "RotatedL": np.array(rot_L),
        }
        return mode_model

    """THIS IS OBSOLETE FOR TENSOR ANALYSES"""

    def topo_OnT_ctrl_segs(self, **kwargs):
        model = self.OnT_ctrl_modes_segs(**kwargs)

        L = model["RotatedL"]
        expl_var = model["Model"].explained_variance_ratio_
        coeffs = np.array(model["Vectors"])
        S = model["RotatedS"]
        s_coeffs = np.array(model["SVectors"])
        # s_expl_var = model['SModel'].explained_variance_ratio_

        # ALL PLOTTING BELOW
        # Plot the topo for our low-rank component
        for comp in range(1):
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(
                np.mean(L, axis=0)[:, comp],
                fig,
                label="OnT Mean Response",
                unwrap=True,
                scale=100,
                alpha=0.3,
                marker_scale=5,
            )
            # plt.title('rPCA Component ' + str(comp))
            plt.title("Median Low-Rank Component " + str(comp))

        # Plot our Components
        plt.figure()
        plt.subplot(221)
        plt.plot(expl_var)
        plt.ylim((0, 1))
        plt.subplot(222)
        for ii in range(4):  # this loops through our COMPONENTS to find the end
            plt.plot(
                coeffs[ii, :].T, linewidth=5 - ii, alpha=expl_var[ii]
            )  # ,alpha=expl_var[ii])
        plt.ylim((-1, 1))
        plt.hlines(0, 0, 3)
        plt.legend(["PC0", "PC1", "PC2", "PC3", "PC4"])
        plt.title("Low-Rank Components")

        plot_sparse = False
        if plot_sparse:
            for comp in range(4):
                # Plot sparse component next
                fig = plt.figure()
                EEG_Viz.plot_3d_scalp(
                    np.median(S, axis=0)[:, comp],
                    fig,
                    label="OnT Mean Response",
                    unwrap=True,
                    scale=100,
                    alpha=0.3,
                    marker_scale=5,
                )
                plt.title("Sparse Component " + str(comp))

            plt.figure()
            plt.subplot(221)
            plt.plot(s_expl_var)
            plt.ylim((0, 1))
            plt.subplot(222)

            for ii in range(4):  # this loops through our COMPONENTS to find the end
                plt.plot(
                    np.mean(s_coeffs, axis=0)[ii, :].T, linewidth=5 - ii
                )  # ,alpha=expl_var[ii])
            plt.ylim((-0.3, 0.3))
            plt.hlines(0, 0, 3)
            plt.legend(["PC0", "PC1", "PC2", "PC3", "PC4"])
            plt.title("Sparse Components")

        if kwargs["plot_maya"]:
            response_dict = np.median(L, axis=0)[:, comp].squeeze()
            EEG_Viz.maya_band_display(response_dict)

    # Dimensionality reduction of ONTarget response; for now rPCA
    def topo_OnT_ctrl_tensor(self, **kwargs):
        pt = kwargs["pt"]
        seg_responses = self.osc_bl_norm[pt]["OnT"][:, :, 0:4]

        # factors = parafac(seg_responses,rank=4)
        print(seg_responses.shape)
        # core, factors = tucker(seg_responses)
        factors = parafac(
            seg_responses.swapaxes(0, 1), rank=4
        )  # factors gives us weight, factors
        # plt.plot(core[0])
        # print(len(factors))
        print(factors)
        for ii in range(4):
            fig = plt.figure()
            # EEG_Viz.plot_3d_scalp(seg_responses[20,:,2],fig,label='Raw Segment',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
            EEG_Viz.plot_3d_scalp(
                factors[1][0][:, ii],
                fig,
                label="Tensor Decomp " + str(ii),
                unwrap=True,
                scale=100,
                alpha=0.3,
                marker_scale=5,
            )
        print(factors[0])
        # print(core.shape)
        # print((factors))

    def OnT_ctrl_modes(self, pt="POOL", data_source=[], do_plot=False, plot_maya=True):

        print("Using BL Norm Segments - RAW")
        med_response = np.median(self.osc_bl_norm[pt]["OnT"], axis=0).squeeze()
        source_label = "BL Normed Segments"

        svm_pca_coeffs = []
        rpca = r_pca.R_pca(med_response)
        L, S = rpca.fit()

        # L = med
        svm_pca = PCA()
        svm_pca.fit(L)
        rotated_L = svm_pca.fit_transform(L)

        svm_pca_coeffs = svm_pca.components_

        # plt.scatter(med_response[:,2],med_response[:,3])

        mode_model = {
            "L": L,
            "S": S,
            "SModel": [],
            "RotatedS": [],
            "SVectors": [],
            "Vectors": svm_pca_coeffs,
            "Model": svm_pca,
            "RotatedL": rotated_L,
        }
        return mode_model

    def topo_OnT_ctrl(self, **kwargs):
        model = self.OnT_ctrl_modes(**kwargs)

        L = model["RotatedL"]
        expl_var = model["Model"].explained_variance_ratio_
        coeffs = np.array(model["Vectors"])
        # s_expl_var = model['SModel'].explained_variance_ratio_

        # ALL PLOTTING BELOW
        # Plot the topo for our low-rank component
        for comp in range(2):
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(
                L[:, comp],
                fig,
                label="OnT Mean Response",
                unwrap=True,
                scale=100,
                alpha=0.3,
                marker_scale=5,
            )
            # plt.title('rPCA Component ' + str(comp))
            plt.title("Median Low-Rank Component " + str(comp))

        # Plot our Components
        plt.figure()
        plt.subplot(221)
        plt.plot(expl_var)
        plt.ylim((0, 1))
        plt.subplot(222)
        for cc in range(4):  # this loops through our COMPONENTS to find the end
            plt.plot(coeffs[cc, :], linewidth=5 - cc, alpha=0.2)  # ,alpha=expl_var[ii])
        plt.ylim((-1, 1))
        plt.hlines(0, 0, 3)
        plt.legend(["PC0", "PC1", "PC2", "PC3", "PC4"])
        plt.title("Low-Rank Components")

        plot_sparse = False
        if plot_sparse:
            for comp in range(4):
                # Plot sparse component next
                fig = plt.figure()
                EEG_Viz.plot_3d_scalp(
                    np.median(S, axis=0)[:, comp],
                    fig,
                    label="OnT Mean Response",
                    unwrap=True,
                    scale=100,
                    alpha=0.3,
                    marker_scale=5,
                )
                plt.title("Sparse Component " + str(comp))

            plt.figure()
            plt.subplot(221)
            plt.plot(s_expl_var)
            plt.ylim((0, 1))
            plt.subplot(222)

            for ii in range(4):  # this loops through our COMPONENTS to find the end
                plt.plot(
                    np.mean(s_coeffs, axis=0)[ii, :].T, linewidth=5 - ii
                )  # ,alpha=expl_var[ii])
            plt.ylim((-0.3, 0.3))
            plt.hlines(0, 0, 3)
            plt.legend(["PC0", "PC1", "PC2", "PC3", "PC4"])
            plt.title("Sparse Components")

        if kwargs["plot_maya"]:
            # response_dict = np.median(L,axis=0)#[:,comp].squeeze()
            response = L[:, 1].squeeze()
            EEG_Viz.maya_band_display(response)
            # EEG_Viz.plot_3d_scalp(response)

    def dict_all_obs(self, condits=["OnT"]):
        full_stack = nestdict()
        label_map = self.label_map

        for condit in condits:
            full_stack[condit]["x"] = self.osc_bl_norm["POOL"][condit]
            full_stack[condit]["g"] = [
                label_map[condit] for seg in self.osc_bl_norm["POOL"][condit]
            ]

        # List comprehend version
        # full_stack = {condit:}

        return full_stack

    def control_rotate(self, condits=["OnT", "OffT"]):
        # get our bases
        control_bases = self.control_bases
        control_modes = self.control_modes
        # Get our observations, ONT and OFFT alike
        obs_dict = self.dict_all_obs(condits=condits)

        trajectories = nestdict()
        rotated_stack = nestdict()
        for condit in condits:
            rotated_stack[condit] = np.dot(obs_dict[condit]["x"], control_bases)

            for ii in range(2):
                trajectories[condit][ii] = np.dot(
                    rotated_stack[condit][:, :, ii].squeeze(),
                    control_modes[:, ii].squeeze(),
                )

        # PLOTTING HERE
        plt.figure()
        # plt.subplot(121)
        for ii in [0, -1]:
            plt.scatter(
                trajectories["OnT"][0][ii],
                trajectories["OnT"][1][ii],
                cmap="jet",
                marker="s",
            )
        plt.scatter(trajectories["OnT"][0], trajectories["OnT"][1], cmap="jet")
        # plt.plot(trajectories['OnT'][0],trajectories['OnT'][1],alpha=0.8)
        plt.xlim([-100, 100])
        plt.ylim([-100, 100])

        plt.figure()
        # plt.subplot(122)
        plt.scatter(
            trajectories["OnT"][0],
            trajectories["OnT"][1],
            c=np.linspace(0, 1, trajectories["OnT"][0].shape[0]),
            cmap="jet",
            alpha=0.1,
        )
        # plt.plot(trajectories['OnT'][0],trajectories['OnT'][1],alpha=0.1)

        for ii in [0, -1]:
            plt.scatter(
                trajectories["OffT"][0][ii],
                trajectories["OffT"][1][ii],
                cmap="jet",
                marker="s",
            )
        plt.scatter(trajectories["OffT"][0], trajectories["OffT"][1], cmap="jet")
        # plt.plot(trajectories['OffT'][0],trajectories['OffT'][1],alpha=0.8)
        plt.xlim([-100, 100])
        plt.ylim([-100, 100])

    # in this method, we're going to do per-channel statistics for each patient, channel, band

    def band_stats(self, do_band="Alpha"):
        self.pop_meds()

    def plot_band_stats(self, do_band="Alpha"):
        self.plot_meds(band=do_band, flatten=not self.pretty)

    def find_seg_covar(self):
        covar_matrix = nestdict()

        for condit in self.condits:
            seg_stack = sig.detrend(self.GMM_Osc_stack[condit], axis=1, type="constant")
            seg_num = seg_stack.shape[1]

            for bb, band in enumerate(DEFAULT_FEAT_ORDER):
                covar_matrix[condit][band] = []
                for seg in range(seg_num):

                    net_vect = seg_stack[:, seg, bb].reshape(-1, 1)

                    cov_matr = np.dot(net_vect, net_vect.T)
                    covar_matrix[condit][band].append(cov_matr)

                covar_matrix[condit][band] = np.array(covar_matrix[condit][band])

        self.cov_segs = covar_matrix

    def plot_seg_covar(self, band="Alpha"):
        for condit in self.condits:
            plt.figure()
            plt.subplot(211)
            plt.imshow(np.median(self.cov_segs[condit][band], axis=0), vmin=0, vmax=0.5)
            plt.colorbar()
            plt.title("Mean Covar")

            plt.subplot(212)
            plt.imshow(np.var(self.cov_segs[condit][band], axis=0), vmin=0, vmax=2)
            plt.colorbar()
            plt.title("Var Covar")

            plt.suptitle(condit)

    def OnT_v_OffT_MAD(self, band="Alpha"):
        Xdsgn = self.SVM_stack
        # THIS DOES MAD across the segments!!!!!!
        X_onT = Xdsgn[self.SVM_labels == "OnTON", :, :].squeeze()
        X_offT = Xdsgn[self.SVM_labels == "OffTON", :, :].squeeze()
        X_NONE = Xdsgn[self.SVM_labels == "OFF", :, :].squeeze()

        OnT_MAD = robust.mad(X_onT, axis=0)
        OffT_MAD = robust.mad(X_offT, axis=0)
        NONE_MAD = robust.mad(X_NONE, axis=0)

        print("OnT segs: " + str(X_onT.shape[0]))
        print("OffT segs: " + str(X_offT.shape[0]))
        print("OFF segs: " + str(X_NONE.shape[0]))

        self.Var_Meas = {
            "OnT": {"Med": np.median(X_onT, axis=0), "MAD": OnT_MAD},
            "OffT": {"Med": np.median(X_offT, axis=0), "MAD": OffT_MAD},
            "OFF": {"Med": np.median(X_NONE, axis=0), "MAD": NONE_MAD},
        }

    def plot_pca_decomp(self, pca_condit="OnT", approach="rpca"):
        self.pca_decomp(
            direction="channels",
            condit=pca_condit,
            bl_correct=True,
            pca_type=approach,
            plot_distr=True,
        )

        plt.figure()
        plt.subplot(221)
        plt.imshow(self.PCA_d.components_, cmap=plt.cm.jet, vmax=1, vmin=-1)
        plt.colorbar()
        plt.subplot(222)
        plt.plot(self.PCA_d.components_)
        plt.ylim((-1, 1))
        plt.legend(["PC0", "PC1", "PC2", "PC3", "PC4"])
        plt.xticks(np.arange(0, 5), ["Delta", "Theta", "Alpha", "Beta", "Gamma1"])
        plt.subplot(223)

        plt.plot(self.PCA_d.explained_variance_ratio_)
        plt.ylim((0, 1))

        for cc in range(2):

            # plot the boring views first
            plt.figure()
            plt.subplot(211)
            plt.plot(self.PCA_x[:, cc])
            plt.subplot(212)
            plt.hist(self.PCA_x[:, cc], bins=np.linspace(-1, 1, 50))
            # find the top mode
            chann_high = np.where(np.abs(self.PCA_x[:, cc]) > 0.7)
            print(chann_high)

            # Plot the 3d scalp distribution
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(
                self.PCA_x[:, cc], fig, animate=False, unwrap=True, highlight=chann_high
            )
            plt.title("Plotting component " + str(cc))
            plt.suptitle(approach + " rotated results for " + pca_condit)

    def shape_GMM_dsgn(self, inStack_dict, band="Alpha", mask_channs=False):
        segs_feats = nestdict()

        inStack = inStack_dict["Stack"]

        assert inStack["OnT"].shape[2] < 10
        for condit in self.condits:
            num_segs = inStack[condit].shape[1]

            # I think this makes it nsegs x nchann x nfeats?
            segs_chann_feats = np.swapaxes(inStack[condit], 1, 0)

            # this stacks all the ALPHAS for all channels together

            if mask_channs:
                chann_mask = self.median_mask
            else:
                chann_mask = np.array([True] * 257)

            # CHOOSE ALPHA HERE

            if band == "Alpha":
                band_idx = DEFAULT_FEAT_ORDER.index(band)
                segs_feats[condit] = segs_chann_feats[:, chann_mask, band_idx]

            elif band == "All":
                segs_feats[condit] = segs_chann_feats[:, chann_mask, :]

            # segs_feats = np.reshape(segs_chann_feats,(num_segs,-1),order='F')

            # We want a 257 dimensional vector with num_segs observations

        return segs_feats

    """Plot the population medians"""

    def pop_meds(self, response=True, pt="POOL"):
        print("Doing Population Meds/Mads on Oscillatory RESPONSES")

        # THIS IS THE OLD WAY: #dsgn_X = self.shape_GMM_dsgn(self.gen_GMM_Osc(self.gen_GMM_stack(stack_bl='normalize')['Stack']),band='All')
        if response:
            dsgn_X = self.osc_bl_norm[pt]
        else:
            dsgn_X = self.osc_stim[pt]

        X_med = nestdict()
        X_mad = nestdict()
        X_segnum = nestdict()
        # do some small simple crap here

        # Here we're averaging across axis zero which corresponds to 'averaging' across SEGMENTS
        for condit in self.condits:
            # Old version just does one shot median
            X_med[condit] = 10 * np.median(dsgn_X[condit], axis=0)
            X_med[condit] = 10 * np.mean(dsgn_X[condit], axis=0)

            # VARIANCE HERE

            X_mad[condit] = robust.mad(dsgn_X[condit], axis=0)
            # X_mad[condit] = np.var(dsgn_X[condit],axis=0)
            X_segnum[condit] = dsgn_X[condit].shape[0]

        self.Seg_Med = (X_med, X_mad, X_segnum)

        weigh_mad = 0.3
        try:
            self.median_mask = (
                np.abs(self.Seg_Med[0]["OnT"][:, 2])
                - weigh_mad * self.Seg_Med[1]["OnT"][:, 2]
                >= 0
            )
        except:
            pdb.set_trace()

        # Do a quick zscore to zero out the problem channels
        chann_patt_zs = stats.zscore(X_med["OnT"], axis=0)
        outlier_channs = np.where(chann_patt_zs > 3)

    def pop_meds_jk(self, response=True):
        print("Doing Population Meds/Mads on Oscillatory RESPONSES - JackKnife Version")

        # THIS IS THE OLD WAY: #dsgn_X = self.shape_GMM_dsgn(self.gen_GMM_Osc(self.gen_GMM_stack(stack_bl='normalize')['Stack']),band='All')
        if response:
            dsgn_X = self.osc_bl_norm["POOL"]
        else:
            dsgn_X = self.osc_stim["POOL"]

        X_med = nestdict()
        X_mad = nestdict()
        X_segnum = nestdict()
        # do some small simple crap here

        # Here we're averaging across axis zero which corresponds to 'averaging' across SEGMENTS
        for condit in self.condits:
            # this version does jackknifing of the median estimate
            ensemble_med = dbo.jk_median(dsgn_X[condit])
            X_med[condit] = np.median(ensemble_med, axis=0)

            # VARIANCE HERE
            X_mad[condit] = np.std(ensemble_med, axis=0)
            # X_mad[condit] = robust.mad(dsgn_X[condit],axis=0)
            # X_mad[condit] = np.var(dsgn_X[condit],axis=0)
            X_segnum[condit] = dsgn_X[condit].shape[0]

        self.Seg_Med = (X_med, X_mad, X_segnum)

        weigh_mad = 0.3
        try:
            self.median_mask = (
                np.abs(self.Seg_Med[0]["OnT"][:, 2])
                - weigh_mad * self.Seg_Med[1]["OnT"][:, 2]
                >= 0
            )
        except:
            pdb.set_trace()

        # Do a quick zscore to zero out the problem channels
        chann_patt_zs = stats.zscore(X_med["OnT"], axis=0)
        outlier_channs = np.where(chann_patt_zs > 3)

    def do_ICA_fullstack(self):
        rem_channs = False
        print("ICA Time")
        ICA_inX = X_med["OnT"]
        if rem_channs:
            ICA_inX[outlier_channs, :] = np.zeros_like(ICA_inX[outlier_channs, :])

        PCA_inX = np.copy(ICA_inX)

        # Go ahead and do PCA here since the variables are already here

        # PCA SECTION
        # pca = PCA()

        # ICA
        ica = FastICA(n_components=5)
        ica.fit(ICA_inX)
        self.ICA_d = ica
        self.ICA_inX = ICA_inX
        self.ICA_x = ica.fit_transform(ICA_inX)

    def plot_PCA_stuff(self):
        plt.figure()
        plt.subplot(221)
        plt.imshow(self.PCA_d.components_, cmap=plt.cm.jet, vmax=1, vmin=-1)
        plt.colorbar()

        plt.subplot(222)
        plt.plot(self.PCA_d.components_)
        plt.ylim((-1, 1))
        plt.legend(["PC0", "PC1", "PC2", "PC3", "PC4"])
        plt.xticks(np.arange(0, 5), ["Delta", "Theta", "Alpha", "Beta", "Gamma1"])

        plt.subplot(223)
        plt.plot(self.PCA_d.explained_variance_ratio_)
        plt.ylim((0, 1))

        for cc in range(2):
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(self.PCA_x[:, cc], fig, animate=False, unwrap=True)
            plt.title("Plotting component " + str(cc))
            plt.suptitle("PCA rotated results for OnT")

    def plot_ICA_stuff(self):
        plt.figure()
        plt.subplot(221)
        plt.imshow(self.ICA_d.components_[:, :-1], cmap=plt.cm.jet, vmax=1, vmin=-1)
        plt.colorbar()
        plt.subplot(222)
        plt.plot(self.ICA_d.components_[:, :-1])
        plt.legend(["IC0", "IC1", "IC2", "IC3", "IC4"])
        plt.xticks(np.arange(0, 5), ["Delta", "Theta", "Alpha", "Beta", "Gamma1"])
        plt.subplot(223)

        plt.plot(self.ICA_d.mixing_)

        for cc in range(2):
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(self.ICA_x[:, cc], fig, animate=False, unwrap=True)
            plt.title("Plotting component " + str(cc))
            plt.suptitle("ICA rotated results for OnT")

    """gets us a histogram of the mads"""

    def band_mads(self):
        pass

    """Plot distribution of change for individual bands"""

    def band_distr(self, do_moment="meds"):
        print("Plotting Distribution for Bands (OLD)")

        meds = nestdict()
        mads = nestdict()

        marker = ["o", "s"]
        color = ["b", "g"]
        plt.figure()
        ax2 = plt.subplot(111)
        for cc, condit in enumerate(["OnT", "OffT"]):
            # for bb in range(5):
            meds[condit] = self.Seg_Med[0][condit][
                :, :
            ]  # result here is 257(chann) x 5(bands)
            mads[condit] = self.Seg_Med[1][condit][:, :]
            # band_segnum[condit] = self.Seg_Med[2][condit]
            # plt.plot(np.arange(0,5)+(cc-0.5)/10,meds[condit][:,:].T,color[cc]+'.',markersize=20,alpha=0.05)
            # There's a way to do MATCHED CHANGES here!! TODO

            # plt.scatter((bb+(cc-0.5)/10)*np.ones_like(meds[condit][:,bb]),meds[condit][:,bb],marker=marker[cc],color=color[cc],s=100,alpha=0.2)
            # plt.boxplot(meds[condit][:,:],positions=np.arange(5)+(cc-0.5)/10,labels=feat_order)
            if do_moment == "meds":
                parts = ax2.violinplot(
                    meds[condit][:, :], positions=np.arange(5) + (cc - 0.5) / 10
                )
            elif do_moment == "mads":
                parts = ax2.violinplot(
                    mads[condit][:, :], positions=np.arange(5) + (cc - 0.5) / 10
                )

            for pc in parts["bodies"]:
                pc.set_facecolor(color[cc])
                pc.set_edgecolor(color[cc])
                # pc.set_linecolor(color[cc])

            # plt.ylim((-0.5,0.5))

        for bb in range(5):
            # rsres = stats.ranksums(meds['OnT'][:,bb],meds['OffT'][:,bb])
            # rsres = stats.wilcoxon(meds['OnT'][:,bb],meds['OffT'][:,bb])
            rsres = stats.ttest_ind(
                10 ** (meds["OnT"][:, bb] / 10), 10 ** (meds["OffT"][:, bb] / 10)
            )
            print(rsres)

        # plt.suptitle(condit)
        plt.ylim((-50, 50))
        plt.hlines(0, -1, 4, linestyle="dotted")
        plt.legend(["OnTarget", "OffTarget"])

    def plot_meds(self, band="Alpha", flatten=True, condits=["OnT", "OffT"]):
        print("Doing Population Level Medians and MADs")

        band_median = {key: 0 for key in self.condits}
        band_mad = {key: 0 for key in self.condits}
        band_segnum = {key: 0 for key in self.condits}

        if band == "DSV":
            # lridge = [-0.00583578, -0.00279751,  0.00131825,  0.01770169,  0.01166687]
            # rridge = [-1.06586005e-02,  2.42700023e-05,  7.31445236e-03,  2.68723035e-03,-3.90440108e-06]
            doridge = np.array(
                [-0.00583578, -0.00279751, 0.00131825, 0.01770169, 0.01166687]
            )
            # doridge = np.array([-1.06586005e-02,  2.42700023e-05,  7.31445236e-03,  2.68723035e-03,-3.90440108e-06])
            doridge = doridge / np.linalg.norm(doridge)
            band_idx = np.array([0, 1, 2, 3, 4])
        else:
            band_idx = DEFAULT_FEAT_ORDER.index(band)
            doridge = [0, 0, 0, 0, 0]
            doridge[band_idx] = 1
            # rridge = [0,0,0,0,0]

        for condit in self.condits:
            band_median[condit] = np.dot(self.Seg_Med[0][condit][:, :], doridge)
            band_mad[condit] = self.Seg_Med[1][condit][:, band_idx]
            band_segnum[condit] = self.Seg_Med[2][condit]

            # band_mad[condit] = self.Seg_Med[1][condit][:,band_idx]
            # band_segnum[condit] = self.Seg_Med[2][condit]

        # Plot the Medians across channels
        plt.figure()
        plt.subplot(211)
        serr_med = {key: 0 for key in self.condits}
        for condit in self.condits:

            plt.plot(band_median[condit], label=condit)
            serr_med[condit] = 1.48 * band_mad[condit] / np.sqrt(band_segnum[condit])

            plt.fill_between(
                np.arange(257),
                band_median[condit] - serr_med[condit],
                band_median[condit] + serr_med[condit],
                alpha=0.4,
            )

        plt.hlines(0, 0, 256)
        plt.title("Medians across Channels")
        plt.legend()
        plt.suptitle(band)

        ##
        # Do Wilcoxon signed rank test
        if "OffT" in band_median.keys():
            WCSRtest = stats.wilcoxon(band_median["OnT"], band_median["OffT"])
            print(WCSRtest)

            # This is the plot of MADs
            plt.subplot(212)
            plt.plot(serr_med["OnT"], label="OnT")
            plt.plot(serr_med["OffT"], label="OffT")
            plt.title("Normed MADs across Channels")
            plt.legend()

        # plot EEG TOPOLOGICAL change for conditions
        for condit in self.condits:
            fig = plt.figure()
            # This is MEDS
            EEG_Viz.plot_3d_scalp(
                band_median[condit],
                fig,
                label=condit + "_med",
                animate=False,
                clims=(-0.2, 0.2),
                unwrap=flatten,
            )
            plt.suptitle(
                "Median of Cortical Response across all "
                + condit
                + " segments | Band is "
                + band
            )

        for condit in self.condits:
            # let's plot the exterior matrix for this
            fig = plt.figure()
            band_corr_matr = (
                band_median[condit].reshape(-1, 1)
                * band_median[condit].reshape(-1, 1).T
            )
            # pdb.set_trace()
            plt.imshow(band_corr_matr, vmin=-0.01, vmax=0.05)
            plt.colorbar()

        # plot the scalp EEG changes
        for condit in self.condits:
            fig = plt.figure()
            # this is MADs
            EEG_Viz.plot_3d_scalp(
                band_mad[condit],
                fig,
                label=condit + "_mad",
                animate=False,
                unwrap=flatten,
                clims=(0, 1.0),
            )
            plt.suptitle(
                "MADs of Cortical Response across all "
                + condit
                + " segments | Band is "
                + band
            )

        plt.suptitle(band)

        # Finally, for qualitative, let's look at the most consistent changes
        # This is the MASKED EEG channels
        if 0:
            for condit in self.condits:
                weigh_mad = 0.4
                fig = plt.figure()
                masked_median = self.Seg_Med[0][condit][:, band_idx] * (
                    np.abs(self.Seg_Med[0][condit][:, band_idx])
                    - weigh_mad * self.Seg_Med[1][condit][:, band_idx]
                    >= 0
                ).astype(int)
                EEG_Viz.plot_3d_scalp(
                    masked_median,
                    fig,
                    label=condit + "_masked_median",
                    animate=False,
                    clims=(-0.1, 0.1),
                )
                plt.suptitle(
                    "Medians with small variances ("
                    + str(weigh_mad)
                    + ") "
                    + condit
                    + " segments | Band is "
                    + band
                )

        olap = {key: 0 for key in self.condits}

        ## Figure out which channels have overlap
        for condit in self.condits:
            olap[condit] = np.array(
                (
                    band_median[condit],
                    band_median[condit]
                    - band_mad[condit] / np.sqrt(band_segnum[condit]),
                    band_median[condit]
                    + band_mad[condit] / np.sqrt(band_segnum[condit]),
                )
            )

        for cc in range(257):
            np.hstack((olap["OnT"][1][cc], olap["OnT"][2][cc]))

    def train_GMM(self):
        # shape our dsgn matrix properly
        intermed_X = self.gen_GMM_Osc(self.gen_GMM_stack(stack_bl="normalize")["Stack"])
        dsgn_X = self.shape_GMM_dsgn(intermed_X, mask_channs=True)

        # let's stack the two together and expect 3 components?
        # this is right shape: (n_samples x n_features)

        # break out our dsgn matrix
        condit_dict = [val for key, val in dsgn_X.items()]

        full_X = np.concatenate(condit_dict, axis=0)

        # setup our covariance prior from OUR OTHER ANALYSES

        covariance_prior = self.gen_GMM_priors(mask_chann=True)

        # when below reg_covar was 1e-1 I got SOMETHING to work
        # gMod = mixture.BayesianGaussianMixture(n_components=self.num_gmm_comps,mean_prior=self.Seg_Med[0]['OnT'],mean_precision_prior=0.1,covariance_type='full',covariance_prior=np.dot(self.Seg_Med[1]['OnT'].reshape(-1,1),self.Seg_Med[1]['OnT'].reshape(-1,1).T),reg_covar=1,tol=1e-2)#,covariance_prior=covariance_prior_alpha)

        # BAYESIAN GMM version
        gMod = mixture.BayesianGaussianMixture(
            n_components=self.num_gmm_comps,
            mean_precision_prior=0.1,
            covariance_type="full",
            reg_covar=1e-6,
            tol=1e-2,
        )  # ,covariance_prior=covariance_prior_alpha)
        condit_mean_priors = [np.median(rr) for rr in condit_dict]

        gMod.means_ = condit_mean_priors
        gMod.covariance_prior_ = covariance_prior

        # STANDARD GMM
        # gMod = mixture.GaussianMixture(n_components=self.num_gmm_comps)

        try:
            gMod.fit(full_X)
        except Exception as e:
            print(e)

        self.GMM = gMod

        self.predictions = gMod.predict(full_X)
        self.posteriors = gMod.predict_proba(full_X)

    def train_newGMM(self, mask=False):
        num_segs = self.SVM_stack.shape[0]

        # generate a mask
        if mask:
            sub_X = self.SVM_stack[:, self.median_mask, :]
            dsgn_X = sub_X.reshape(num_segs, -1, order="C")
        else:
            dsgn_X = self.SVM_stack.reshape(num_segs, -1, order="C")

        # doing a one class SVM
        # clf = svm.OneClassSVM(nu=0.1,kernel="rbf", gamma=0.1)
        # clf = svm.LinearSVC(penalty='l2',dual=False)
        # gMod = mixture.BayesianGaussianMixture(n_components=3,mean_precision_prior=0.1,covariance_type='full',reg_covar=1e-6,tol=1e-2)#,covariance_prior=covariance_prior_alpha)
        gMod = mixture.BayesianGaussianMixture(n_components=3)

        # split out into test and train
        Xtr, Xte, Ytr, Yte = sklearn.model_selection.train_test_split(
            dsgn_X, self.SVM_labels, test_size=0.33
        )

        gMod.fit(Xtr, Ytr)

        # predict IN training set
        predlabels = gMod.predict(Xte)

        plt.figure()
        plt.plot(Yte, label="test")
        plt.plot(predlabels, label="predict")

        print(np.sum(np.array(Yte) == np.array(predlabels)) / len(Yte))

    def assess_dynamics(self, band="Alpha"):
        band_idx = DEFAULT_FEAT_ORDER.index(band)
        self.OnT_v_OffT_MAD()

        # Now, move on to plotting
        for stat in ["Med", "MAD"]:
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(
                self.Var_Meas["OnT"][stat][:, band_idx],
                fig,
                clims=(0, 0),
                label="OnT " + stat,
                unwrap=True,
            )
            plt.suptitle("Non-normalized Power " + stat + " in " + band + " OnT")

            plt.figure()
            plt.bar(np.arange(1, 258), self.Var_Meas["OnT"][stat][:, band_idx])

            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(
                self.Var_Meas["OffT"][stat][:, band_idx],
                fig,
                clims=(0, 0),
                label="OffT " + stat,
                unwrap=True,
            )
            plt.suptitle("Non-normalized Power " + stat + " in " + band + " OffT")
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(
                self.Var_Meas["OFF"][stat][:, band_idx],
                fig,
                clims=(0, 0),
                label="OFF " + stat,
                unwrap=True,
            )
            plt.suptitle("Non-normalized Power " + stat + " in " + band + " OFF")

            plt.figure()
            plt.subplot(211)
            plt.hist(self.Var_Meas["OnT"]["Med"][:, band_idx], label="OnT", bins=30)
            plt.hist(self.Var_Meas["OFF"]["Med"][:, band_idx], label="OFF", bins=30)
            plt.title("Distributions of Medians")

            plt.subplot(212)
            plt.hist(
                [
                    self.Var_Meas["OnT"]["MAD"][:, band_idx],
                    self.Var_Meas["OFF"]["MAD"][:, band_idx],
                ],
                label=["OnT", "OFF"],
                bins=30,
            )
            # plt.hist(self.Var_Meas['OFF']['MAD'][:,band_idx],label='OFF',bins=30)
            plt.title("Distributions of MADs")
            plt.legend()

    def assess_binSVM(self, mask=False):
        num_segs = self.SVM_stack.shape[0]
        print("DOING BINARY")
        # generate a mask
        if mask:
            # what mask do we want?
            # self.SVM_Mask = self.median_mask
            self.SVM_Mask = np.zeros((257,)).astype(bool)
            self.SVM_Mask[np.arange(216, 239)] = True

            sub_X = self.SVM_stack[:, self.SVM_Mask, :]
            dsgn_X = sub_X.reshape(num_segs, -1, order="C")
        else:
            dsgn_X = self.SVM_stack.reshape(num_segs, -1, order="C")

        # doing a one class SVM
        # clf = svm.OneClassSVM(nu=0.1,kernel="rbf", gamma=0.1)

        # get rid of ALL OFF, and only do two labels
        OFFs = np.where(self.SVM_labels == "OFF")
        dsgn_X = np.delete(dsgn_X, OFFs, 0)
        SVM_labels = np.delete(self.SVM_labels, OFFs, 0)

        # Just doing a learning curve on the training data
        tsize, tscore, vscore = learning_curve(
            svm.LinearSVC(penalty="l2", dual=False, C=1),
            dsgn_X,
            SVM_labels,
            train_sizes=np.linspace(0.1, 1, 20),
            shuffle=True,
            cv=5,
            random_state=12342,
        )
        plt.figure()
        plt.plot(tsize, np.mean(tscore, axis=1))
        plt.plot(tsize, np.mean(vscore, axis=1))
        plt.legend(["Training Score", "Cross-validation Score"])

    """ Learning curve for the Binary SVM"""

    def learning_binSVM(self, mask=False):
        label_map = {"OnT": 1, "OffT": 0}

        SVM_stack = np.concatenate(
            [self.osc_bl_norm["POOL"][condit] for condit in self.condits], axis=0
        )
        SVM_labels = np.concatenate(
            [
                [label_map[condit] for seg in self.osc_bl_norm["POOL"][condit]]
                for condit in self.condits
            ],
            axis=0,
        )
        num_segs = SVM_stack.shape[0]

        dsgn_X = SVM_stack.reshape(num_segs, -1, order="C")

        print("DOING BINARY - Learning Curve")
        tsize, tscore, vscore = learning_curve(
            svm.LinearSVC(penalty="l2", dual=False, C=1),
            dsgn_X,
            SVM_labels,
            train_sizes=np.linspace(0.1, 1, 20),
            shuffle=True,
            cv=10,
            random_state=12342,
        )
        plt.figure()
        plt.plot(tsize, np.mean(tscore, axis=1))
        plt.plot(tsize, np.mean(vscore, axis=1))
        plt.legend(["Training Score", "Cross-validation Score"])

    """ The old method of generating a stack from our observations
    May more convoluted, vestiges of code from attempts to do GMM classification
    Need to phase this out completely
    Except it gives results in $\gamma$ that make more sense superficially
    """

    """ This function retrieves a design matrix from the pooled observations """

    def stack_dsgn(self, do_plot=False):
        label_map = self.label_map

        # New method here...
        SVM_stack = np.concatenate(
            [self.osc_bl_norm["POOL"][condit] for condit in self.condits], axis=0
        )
        num_segs = SVM_stack.shape[0]
        flat_dsgn_X = SVM_stack.reshape(
            num_segs, -1, order="C"
        )  # Here we're collapsing all the 5 Osc x 256 Chann features

        dsgn_Y = np.concatenate(
            [
                [label_map[condit] for seg in self.osc_bl_norm["POOL"][condit]]
                for condit in self.condits
            ],
            axis=0,
        )

        if do_plot:
            # collapse along all segments and channels
            plot_stack = dsgn_X["OnT"].swapaxes(0, 2).reshape(5, -1, order="C")
            plt.figure()
            for ii in range(5):
                sns.violinplot(data=plot_stack.T, positions=np.arange(5))

        self.SVM_raw_stack = SVM_stack
        return flat_dsgn_X, dsgn_Y, num_segs

    """ WIP SVM masked classifier where channels can be toggled """

    def masked_SVM(self):
        # generate a mask

        # what mask do we want?
        self.SVM_Mask = np.zeros((257,)).astype(bool)
        self.SVM_Mask[np.arange(216, 239)] = True

        sub_X = self.SVM_stack[:, self.SVM_Mask, :]
        dsgn_X = sub_X.reshape(num_segs, -1, order="C")

    """ Train our Binary SVM """

    def train_binSVM(self, mask=False):
        self.bin_classif = nestdict()

        dsgn_X, SVM_labels, num_segs = self.stack_dsgn()

        # Next, we want to split out a validation set
        Xtr, self.Xva, Ytr, self.Yva = sklearn.model_selection.train_test_split(
            dsgn_X, SVM_labels, test_size=0.7, shuffle=True, random_state=None
        )

        # Next, we want to do CV learning on just the training set
        # Ensemble variables
        big_score = []
        coeffs = []
        models = []

        # Parameters for CV
        nfold = 10
        cv = StratifiedKFold(n_splits=nfold)
        for train, test in cv.split(Xtr, Ytr):
            clf = svm.LinearSVC(penalty="l1", dual=False, C=1)  # l2 works well here...
            mod_score = clf.fit(Xtr[train], Ytr[train]).score(Xtr[test], Ytr[test])
            outpred = clf.predict(Xtr[test])
            coeffs.append(clf.coef_)
            big_score.append(mod_score)
            models.append(clf)
            # Maybe do ROC stuff HERE? TODO

        # Plot the big score for the CVs
        plt.figure()
        plt.plot(big_score)
        plt.title("Plotting the fit scores for the CV training procedure")

        # Find the best model
        best_model_idx = np.argmax(big_score)
        best_model = models[best_model_idx]
        self.bin_classif["Model"] = best_model
        self.bin_classif["Coeffs"] = coeffs
        self.cv_folding = nfold

        # if we want to keep all of our models we can do that here
        self.cv_bin_classif = models

    """This function ASSESSES the best SVM classifier using a bootstrap procedure"""

    def bootstrap_binSVM(self):
        best_model = self.bin_classif

        # randomlt sample the validation set
        validation_accuracy = []
        rocs_auc = []
        total_segments = self.Xva.shape[0]
        for ii in range(100):
            Xva_ss, Yva_ss = resample(
                self.Xva,
                self.Yva,
                n_samples=int(round(total_segments * 0.6)),
                replace=False,
            )
            validation_accuracy.append(best_model["Model"].score(Xva_ss, Yva_ss))
            predicted_Y = best_model["Model"].predict(Xva_ss)
            fpr, tpr, _ = roc_curve(Yva_ss, predicted_Y)
            rocs_auc.append(auc(fpr, tpr))

        plt.figure()
        plt.subplot(311)
        plt.hist(validation_accuracy)
        plt.subplot(312)
        plt.plot(fpr, tpr)
        plt.subplot(313)
        plt.hist(rocs_auc)
        plt.suptitle("Bootstrapped SVM Assessment")

    """One shot assessment of SVM classifier"""

    def oneshot_binSVM(self):
        best_model = self.bin_classif
        # Plotting of confusion matrix and coefficients
        # Validation set assessment now

        validation_accuracy = best_model["Model"].score(self.Xva, self.Yva)
        Ypred = best_model["Model"].predict(self.Xva)
        print(validation_accuracy)

        plt.figure()
        plt.subplot(1, 2, 1)
        # confusion matrix here
        conf_matrix = confusion_matrix(Ypred, self.Yva)
        plt.imshow(conf_matrix)
        plt.yticks(np.arange(0, 2), ["OffT", "OnT"])
        plt.xticks(np.arange(0, 2), ["OffT", "OnT"])
        plt.colorbar()

        plt.subplot(2, 2, 2)
        coeffs = (
            np.array(best_model["Coeffs"])
            .squeeze()
            .reshape(self.cv_folding, 257, 5, order="C")
        )
        # plt.plot(coeffs,alpha=0.2)
        plt.plot(np.median(coeffs, axis=0))
        plt.title("Plotting Median Coefficients for CV-best Model performance")

        plt.subplot(2, 2, 4)
        plt.plot(np.median(np.median(coeffs, axis=0), axis=0))
        plt.suptitle("Oneshot SVM Assessment")

        # self.SVM_coeffs = coeffs

    def assess_binSVM(self):
        best_model = self.bin_classif

        for ii in range(100):
            Xrs, Yrs = resample(self.Xva, self.Yva, 100)
            valid_accuracy = best_model.score(Xva, Yva)

    """Analysis of the binary SVM coefficients should be here"""

    def analyse_binSVM(self, feature_weigh=False):

        # What exactly is this trying to do here??
        # plt.figure()
        # for ii in range(4):
        #     #plt.hist(self.bin_classif['Model'].coef_[:,ii])
        #     plt.scatter(ii,self.bin_classif['Model'].coef_[:,ii])

        # BELOW (order = 'C') IS CORRECT since before, in the features, we collapse to a feature vector that is all 257 deltas, then all 257 thetas, etc...
        # So when we want to reshape that to where we are now, we have to either 'C': (5,257) where C means the last index changes fastest; or 'F': (257,5) where the first index changes fastest.

        if not feature_weigh:
            coeffs = stats.zscore(
                np.sum(
                    np.abs(self.bin_classif["Model"].coef_.reshape(5, 257, order="C")),
                    axis=0,
                )
            )  # what we have here is a reshape where the FEATURE VECTOR is [257 deltas... 257 gammas]
            self.import_mask = coeffs > 0.2
            analysis_title = "Pure Coefficients"
        else:
            avg_feat_value = np.abs(np.median(self.SVM_raw_stack, axis=0))
            # pdb.set_trace()
            coeffs = stats.zscore(
                np.sum(
                    np.multiply(
                        np.abs(
                            self.bin_classif["Model"].coef_.reshape(5, 257, order="C")
                        ),
                        avg_feat_value.T,
                    ),
                    axis=0,
                )
            )
            self.import_mask = coeffs > 0.1
            analysis_title = "Empirical Feature-weighed"

        # WTF is this shit?
        # avg_coeffs = np.mean(np.array(self.bin_classif['Coeffs']),axis=0).reshape(5,257,order='C')
        # coeffs = stats.zscore(np.sum(avg_coeffs**2,axis=0))

        plt.figure()
        plt.hist(coeffs, bins=50)  # ,range=(0,1))
        # plt.figure()

        EEG_Viz.plot_3d_scalp(coeffs, unwrap=True, alpha=0.4)
        EEG_Viz.plot_3d_scalp(self.import_mask.astype(int), unwrap=True, alpha=0.2)
        plt.suptitle("SVM Coefficients " + analysis_title)

    """Do a CV-level analysis here"""

    def analysis_CV_binSVM(self):
        pass

    def NEWcompute_diff(self):
        avg_change = {
            pt: {
                condit: 10
                * (
                    np.log10(avg_psd[pt][condit][self.keys_of_interest[condit][1]])
                    - np.log10(avg_psd[pt][condit]["Off_3"])
                )
                for pt, condit in itertools.product(self.do_pts, self.condits)
            }
        }

    """ Below are functions related to the oscillatory response characterization"""
    # This goes to the psd change average and computed average PSD across all available patients
    def pop_response(self):
        psd_change_matrix = nestdict()
        population_change = nestdict()
        # pop_psds = defaultdict(dict)
        # pop_psds_var = defaultdict(dict)

        # Generate the full PSDs matrix and then find the MEAN along the axes for the CHANGE
        for condit in self.condits:
            # Get us a matrix of the PSD Changes
            psd_change_matrix[condit] = np.array(
                [rr[condit] for pt, rr in self.psd_change.items()]
            )

            population_change[condit]["Mean"] = np.mean(
                psd_change_matrix[condit], axis=0
            )
            population_change[condit]["Var"] = np.var(psd_change_matrix[condit], axis=0)

        self.pop_psd_change = population_change

    def do_pop_stats(self, band="Alpha", threshold=0.4):
        # self.reliablePSD = nestdict()
        avgOsc = nestdict()
        varOsc = nestdict()
        varOsc_mask = nestdict()

        for condit in self.condits:

            # First, let's average across patients

            # band = np.where(np.logical_and(self.fvect > 14, self.fvect < 30))

            avgOsc[condit] = calc_feats(
                10 ** (self.pop_psd_change[condit]["Mean"] / 10), self.fvect
            )
            varOsc[condit] = calc_feats(
                10 ** (self.pop_psd_change[condit]["Var"]), self.fvect
            )
            # varOsc_mask[condit] = np.array(np.sqrt(varOsc[condit]) < threshold).astype(np.int)

        self.pop_osc_change = avgOsc
        self.pop_osc_var = varOsc
        # self.pop_osc_mask = varOsc_mask

    def do_pop_mask(self, threshold):
        cMask = nestdict()

        # weighedPSD = np.divide(self.pop_stats['Mean'][condit].T,np.sqrt(self.pop_stats['Var'][condit].T))
        # do subtract weight
        # THIS SHOULD ONLY BE USED TO DETERMINE THE MOST USEFUL CHANNELS
        # THIS IS A SHITTY THING TO DO
        for condit in self.condits:
            pre_mask = np.abs(self.pop_change["Mean"][condit].T) - threshold * np.sqrt(
                self.pop_change["Var"][condit].T
            )
            # which channels survive?
            cMask[condit] = np.array(pre_mask > 0).astype(int)

        cMask["Threshold"] = threshold

        self.reliability_mask = cMask

    def plot_pop_stats(self):
        for condit in self.condits:
            plt.figure()
            plt.subplot(311)
            plt.plot(self.fvect, self.pop_stats["Mean"][condit].T)
            plt.title("Mean")
            plt.subplot(312)
            plt.plot(self.fvect, np.sqrt(self.pop_stats["Var"][condit].T) / np.sqrt(3))
            plt.title("Standard Error of the Mean; n=3")

            plt.subplot(313)
            # do divide weight

            reliablePSD = self.reliablePSD[condit]["cPSD"]

            plt.plot(self.fvect, reliablePSD)
            plt.title("SubtrPSD by Pop 6*std")

            plt.xlabel("Frequency (Hz)")
            # plt.subplot(313)
            # plt.hist(weighedPSD,bins=50)

            plt.suptitle(condit + " population level")

    def topo_wrap(
        self, band="Alpha", label="", condit="OnT", mask=False, animate=False
    ):
        mainfig = plt.figure()

        if not mask:
            EEG_Viz.plot_3d_scalp(
                self.pop_osc_change[condit][DEFAULT_FEAT_ORDER.index(band)],
                mainfig,
                animate=animate,
                label=label,
                clims=(-1, 1),
            )
        else:
            try:
                EEG_Viz.plot_3d_scalp(
                    self.pop_osc_mask[condit][DEFAULT_FEAT_ORDER.index(band)]
                    * self.pop_osc_change[condit][DEFAULT_FEAT_ORDER.index(band)],
                    mainfig,
                    clims=(-1, 1),
                    animate=animate,
                    label=label,
                )
            except:
                pdb.set_trace()

        plt.suptitle(label)

    def topo_3d_chann_mask(
        self, band="Alpha", animate=True, pt="all", var_fixed=False, label=""
    ):
        for condit in ["OnT", "OffT"]:
            # plot_vect = np.zeros((257))
            # plot_vect[self.reliablePSD[condit]['CMask']] = 1
            # plot_vect = self.reliablePSD[condit]['CMask']

            mainfig = plt.figure()

            band_idxs = np.where(
                np.logical_and(
                    self.fvect > dbo.feat_dict[band]["param"][0],
                    self.fvect < dbo.feat_dict[band]["param"][1],
                )
            )
            if pt == "all":
                if var_fixed:
                    plot_vect = self.reliablePSD[condit]["BandVect"]["Vect"]
                else:

                    plot_vect = np.median(
                        self.avgPSD[condit]["PSD"][:, band_idxs].squeeze(), axis=1
                    )
                    # abov_thresh = np.median(self.reliablePSD[condit]['cPSD'][band_idxs,:],axis=1).T.reshape(-1)

            else:
                # just get the patient's diff

                plot_vect = np.median(
                    self.feat_diff[pt][condit][:, band_idxs].squeeze(), axis=1
                )

            EEG_Viz.plot_3d_scalp(
                plot_vect, mainfig, clims=(-3, 3), label=condit + label, animate=animate
            )
            plt.title(condit + " " + pt + " " + str(var_fixed))
            plt.suptitle(label)

            self.write_sig(plot_vect)

    def write_sig(self, signature):
        for condit in ["OnT", "OffT"]:
            # plot_vect = np.zeros((257))
            # plot_vect[self.reliablePSD[condit]['CMask']] = 1
            # plot_vect = self.reliablePSD[condit]['CMask']

            # write this to a text file
            np.save("/tmp/" + condit + "_sig.npy", signature)

    def plot_diff(self):
        for pt in self.do_pts:
            plt.figure()
            plt.subplot(221)
            plt.plot(self.fvect, self.psd_change[pt]["OnT"].T)
            plt.title("OnT")

            plt.subplot(222)
            plt.plot(self.fvect, self.psd_change[pt]["OffT"].T)
            plt.title("OffT")

            plt.subplot(223)
            plt.plot(self.fvect, 10 * np.log10(self.psd_var[pt]["OnT"]["BONT"].T))
            plt.title("BONT Variance")

            plt.subplot(224)
            plt.plot(self.fvect, 10 * np.log10(self.psd_var[pt]["OffT"]["BOFT"].T))
            plt.title("BOFT Variance")

            plt.suptitle(pt)

    """GMM Classification Modules here - Analysis is obsolete"""

    def GMM_train(self, condit="OnT"):
        # gnerate our big matrix of observations; Should be 256(chann)x4(feats)x(segxpatients)(observations)
        pass

        # this function will generate a big stack of all observations for a given condition across all patients

    def plot_ontvsofft(self, pt="906"):
        if "OffT" not in self.condits:
            raise ValueError

        condit_diff = (self.PSD_diff[pt]["OnT"] - self.PSD_diff[pt]["OffT"]).T
        plt.figure()
        plt.plot(self.fvect, condit_diff, alpha=0.2)
        plt.xlim((0, 150))
        plt.title("Difference from OnT and OffT")
        plt.suptitle(pt)

    def plot_chann_var(self, pt="906", condit="OnT"):
        plt.figure()
        plt.subplot(121)
        plt.plot(self.PSD_var[pt][condit]["Off_3"].T)
        plt.xlim((0, 150))
        plt.subplot(122)
        plt.plot(self.PSD_var[pt][condit][self.keys_of_interest[condit][1]].T)
        plt.xlim((0, 150))

        plt.suptitle(pt + " " + condit)

    # This function quickly gets the power for all channels in each band
    def Full_feat_band(self):
        pass

    def extract_coher_feats(self, do_pts=[], do_condits=["OnT"], epochs="all"):
        if do_pts == []:
            do_pts = self.do_pts

        PLV_dict = nestdict()
        CSD_dict = nestdict()
        if do_condits == []:
            do_condits = self.condits

        for pt in do_pts:
            for condit in do_condits:
                if epochs == "all":
                    do_epochs = self.keys_of_interest[condit]
                else:
                    do_epochs = epochs

                for epoch in do_epochs:
                    print("Doing " + pt + condit + epoch)
                    data_matr = self.ts_data[pt][condit][epoch]
                    data_dict = {
                        ch: data_matr[ch, :, :].squeeze()
                        for ch in range(self.chann_dim)
                    }
                    (
                        CSD_dict[pt][condit][epoch],
                        PLV_dict[pt][condit][epoch],
                    ) = dbo.gen_coher(
                        data_dict, Fs=self.fs, nfft=2**9, polyord=self.polyorder
                    )

        print("Done with coherence... I guess...")
        return CSD_dict, PLV_dict

    # This is the main wrapper function that leads us to the MNE python plot_topomap
    # The topomap is very misleading, so we try not to use it due to its strange interpolation
    def plot_topo(
        self,
        vect,
        vmax=2,
        vmin=-2,
        label="",
    ):
        plt.figure()
        mne.viz.plot_topomap(
            vect,
            pos=self.eeg_locs.pos[:, [0, 1]],
            vmax=vmax,
            vmin=vmin,
            image_interp="none",
        )
        plt.suptitle(label)

    # Here, we'll plot the PSDs for channels of interest for the conditions of interest
    def psd_stats(self, chann_list=[]):
        self.view_PSDs(chann_list=chann_list, zoom_in=True)

    # This function is to just show the raw PSDs of each of the experimental conditions collected
    def view_PSDs(self, zoom_in=True, chann_list=[], plot_var=False):
        print("Showing raw PSDs")
        avg_psd = nestdict()
        var_psd = nestdict()

        if chann_list == []:
            chann_list = np.arange(256)
        else:
            chann_list = np.array(chann_list)

        f_vect = np.linspace(0, 500, 2**10 + 1)

        for pt in self.do_pts:
            # avg_psd[pt] = defaultdict(dict)
            # avg_change[pt] = defaultdict(dict)
            for condit in self.condits:
                # average all the epochs together
                avg_psd[pt][condit] = {
                    epoch: np.median(self.feat_dict[pt][condit][epoch], axis=1)
                    for epoch in self.feat_dict[pt][condit].keys()
                }
                if plot_var:
                    var_psd[pt][condit] = {
                        epoch: robust.mad(self.feat_dict[pt][condit][epoch], axis=1)
                        for epoch in self.feat_dict[pt][condit].keys()
                    }

            psd_fig = plt.figure()
            plt.subplot(2, 2, 1)
            plt.plot(
                f_vect, 10 * np.log10(avg_psd[pt]["OnT"]["Off_3"][chann_list, :].T)
            )
            plt.title("OnT-Pre")

            plt.subplot(2, 2, 2)
            plt.plot(
                f_vect, 10 * np.log10(avg_psd[pt]["OffT"]["Off_3"][chann_list, :].T)
            )
            plt.title("OffT-Pre")

            plt.subplot(2, 2, 3)
            plt.plot(f_vect, 10 * np.log10(avg_psd[pt]["OnT"]["BONT"][chann_list, :].T))
            plt.title("BONT")

            plt.subplot(2, 2, 4)
            plt.plot(
                f_vect, 10 * np.log10(avg_psd[pt]["OffT"]["BOFT"][chann_list, :].T)
            )
            plt.title("BOFT")

            if plot_var:
                plt.figure(psd_fig.number)
                plt.subplot(2, 2, 1)
                plt.fill_between(
                    f_vect, 10 * np.log10(var_psd[pt]["OnT"]["Off_3"][chann_list, :].T)
                )

            if zoom_in:
                for ii in range(1, 5):
                    plt.subplot(2, 2, ii)
                    plt.xlim(0, 40)
                    plt.ylim(-20, 10)

            plt.suptitle(pt)

    """Coherence Statistics here"""

    def coher_stat(self, pt_list=[], chann_list=[]):
        return self.extract_coher_feats(do_pts=pt_list, do_condits=["OnT", "OffT"])
