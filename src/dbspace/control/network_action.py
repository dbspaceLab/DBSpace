#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 17:34:13 2019

@author: virati
"""

import logging
from copy import deepcopy
from itertools import product as cart_prod

import dbspace as dbo
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from dbspace.control.stream_buffers import streamLFP
from dbspace.signal.oscillations import DEFAULT_FEAT_ORDER
from dbspace.utils.structures import nestdict

logging.basicConfig(
    filename="/tmp/network_action.log",
    filemode="w",
    format="%(name)s - %(levelname)s - %(message)s",
)
logging.info("Starting the log...")



class local_response:
    # Setup our main variables for the analysis
    TF_response = nestdict()
    Osc_response = nestdict()
    Osc_prebilat = nestdict()
    Osc_baseline = nestdict()

    # Ancillary analysis variables
    Osc_response_uncorr = nestdict()

    def __init__(self, config_file, do_pts, analysis_windows=["Bilat", "PreBilat"]):
        # Which two epochs are we analysing?

        self.win_list = analysis_windows

        self.do_pts = do_pts

        self.colors = ["b", "g"]

        self.config_file = config_file

    def extract_baselines(self):
        TF_response = self.TF_response
        Osc_response_uncorr = self.Osc_response_uncorr
        Osc_prebilat = self.Osc_prebilat
        Osc_baseline = self.Osc_baseline

        for pt, condit in cart_prod(self.do_pts, ["OnT", "OffT"]):
            eg_rec = streamLFP(config_file=self.config_file, pt=pt, condit=condit)
            rec = eg_rec.time_series(epoch_name="PreBilat")
            TF_response[pt][condit] = eg_rec.tf_transform(epoch_name="Bilat")
            Osc_response_uncorr[pt][condit] = eg_rec.osc_transform(epoch_name="Bilat")

            Osc_prebilat[pt][condit] = eg_rec.osc_transform(epoch_name="PreBilat")
            # Find the mean within the prebilat for both left and right
            Osc_baseline[pt][condit] = [
                np.mean(Osc_prebilat[pt][condit][chann], axis=0)
                for chann in ["Left", "Right"]
            ]

    def extract_response(self):
        Osc_response = deepcopy(self.Osc_response_uncorr)
        Osc_baseline = self.Osc_baseline

        for pt, condit in cart_prod(self.do_pts, ["OnT", "OffT"]):
            for seg in range(self.Osc_response_uncorr[pt][condit]["Left"].shape[0]):
                for cc, chann in enumerate(["Left", "Right"]):
                    Osc_response[pt][condit][chann][seg, :] -= Osc_baseline[pt][condit][
                        cc
                    ]

        self.Osc_response = Osc_response

    def gen_osc_distr(self):
        Osc_response = self.Osc_response
        do_pts = self.do_pts

        self.Osc_indiv_marg = {
            pt: {
                condit: np.array(
                    (
                        Osc_response[pt][condit]["Left"],
                        Osc_response[pt][condit]["Right"],
                    )
                )
                for condit in ["OnT", "OffT"]
            }
            for pt in do_pts
        }
        self.Osc_indiv_med = {
            pt: {
                condit: np.median(self.Osc_indiv_marg[pt][condit], axis=1)
                for condit in ["OnT", "OffT"]
            }
            for pt in do_pts
        }
        self.Osc_indiv_pop = {
            side: {
                condit: np.array(
                    [self.Osc_indiv_med[pt][condit][ss, :] for pt in do_pts]
                )
                for condit in ["OnT", "OffT"]
            }
            for ss, side in enumerate(["Left", "Right"])
        }

        #%%
        # here we'll work with the oscillatory state variables
        self.Osc_pt_marg = {
            condit: np.array(
                [
                    (
                        Osc_response[pt][condit]["Left"],
                        Osc_response[pt][condit]["Right"],
                    )
                    for pt in do_pts
                ]
            )
            for condit in ["OnT", "OffT"]
        }
        self.Osc_pt_marg_bl = {
            condit: np.array(
                [
                    (
                        self.Osc_prebilat[pt][condit]["Left"],
                        self.Osc_prebilat[pt][condit]["Right"],
                    )
                    for pt in do_pts
                ]
            )
            for condit in ["OnT", "OffT"]
        }
        self.Osc_pt_marg_uncorr = {
            condit: np.array(
                [
                    (
                        self.Osc_response_uncorr[pt][condit]["Left"],
                        self.Osc_response_uncorr[pt][condit]["Right"],
                    )
                    for pt in do_pts
                ]
            )
            for condit in ["OnT", "OffT"]
        }

    def plot_responses(self, do_pts, r_type="pt"):
        for cc, chann in enumerate(["Left", "Right"]):
            # do violin plots
            fig = plt.figure()
            ax2 = plt.subplot(111)
            color = ["b", "g"]
            distr = nestdict()

            for co, condit in enumerate(["OnT", "OffT"]):
                # how many segments?
                # Here, we're going to plot ALL segments, marginalized across patients
                segNum = self.Osc_pt_marg[condit].shape[2]
                if r_type == "pt":
                    distr_to_plot = self.Osc_indiv_pop[chann][condit]
                elif r_type == "segs":
                    distr_to_plot = (
                        self.Osc_pt_marg[condit]
                        .swapaxes(1, 2)
                        .reshape(len(do_pts) * segNum, 2, 5)[:, cc, :]
                    )

                plt.plot(
                    np.arange(1, 6) + 0.2 * co,
                    distr_to_plot.T,
                    color[co] + ".",
                    markersize=20,
                    alpha=0.8,
                )
                parts = ax2.violinplot(
                    distr_to_plot,
                    positions=np.array([1, 2, 3, 4, 5]) + 0.2 * co,
                    showmedians=True,
                )

                for partname in ("cbars", "cmins", "cmaxes", "cmedians"):
                    vp = parts[partname]
                    vp.set_edgecolor(color[co])
                    if partname == "cmedians":
                        vp.set_linewidth(5)
                    else:
                        vp.set_linewidth(2)

                for pc in parts["bodies"]:
                    pc.set_facecolor(color[co])
                    pc.set_edgecolor(color[co])
                    # pc.set_linecolor(color[co])

                distr[condit] = distr_to_plot
                # plt.plot([1,2,3,4,5],np.mean(distr_to_plot,axis=0),color=color[co])
                plt.hlines(0, 0, 5, linestyle="dotted", linewidth=5, alpha=0.8)
                plt.suptitle("Looking at patient averages")

            for bb in range(5):
                print(DEFAULT_FEAT_ORDER[bb])
                # rsres = stats.ranksums(distr['OnT'][:,bb],distr['OffT'][:,bb])
                rsres = stats.ks_2samp(distr["OnT"][:, bb], distr["OffT"][:, bb])

                # rsres = stats.wilcoxon(distr['OnT'][:,bb],distr['OffT'][:,bb])
                # rsres = stats.ttest_ind(distr['OnT'][:,bb],distr['OffT'][:,bb])
                print(rsres)

                # ontres = stats.ranksums(distr['OnT'][:,bb])
                # ontres = stats.kstest(distr['OnT'][:,bb],cdf='norm')
                # ontres = stats.mannwhitneyu(distr['OnT'][:,bb])
                ontres = stats.ttest_1samp(distr["OnT"][:, bb], np.zeros((5, 1)))
                print(condit + " " + str(ontres))

            plt.ylim((-30, 50))
            plt.legend()
            plt.title(chann)

    def plot_segment_responses(self, do_pts):
        for cc, chann in enumerate(["Left", "Right"]):
            # do violin plots
            fig = plt.figure()
            ax2 = plt.subplot(111)
            color = ["b", "g"]
            distr = nestdict()

            for co, condit in enumerate(["OnT", "OffT"]):
                # how many segments?
                # Here, we're going to plot ALL segments, marginalized across patients
                segNum = self.Osc_pt_marg[condit].shape[2]
                distr_to_plot = (
                    self.Osc_pt_marg[condit]
                    .swapaxes(1, 2)
                    .reshape(len(do_pts) * segNum, 2, 5)[:, cc, :]
                )

                # plt.plot(np.arange(1,6)+0.2*co,distr_to_plot.T,color[co]+'.',markersize=20,alpha=0.2)
                parts = ax2.violinplot(
                    distr_to_plot,
                    positions=np.array([1, 2, 3, 4, 5]) + 0.2 * co,
                    showmedians=True,
                )

                for partname in ("cbars", "cmins", "cmaxes", "cmedians"):
                    vp = parts[partname]
                    vp.set_edgecolor(color[co])
                    if partname == "cmedians":
                        vp.set_linewidth(5)
                    else:
                        vp.set_linewidth(2)

                for pc in parts["bodies"]:
                    pc.set_facecolor(color[co])
                    pc.set_edgecolor(color[co])
                    # pc.set_linecolor(color[co])

                distr[condit] = distr_to_plot
                # plt.plot([1,2,3,4,5],np.mean(distr_to_plot,axis=0),color=color[co])
                plt.hlines(0, 0, 5, linestyle="dotted", linewidth=5, alpha=0.8)
                plt.suptitle("Looking at all available segments")

            for bb in range(5):
                print(bb)
                # rsres = stats.ranksums(distr['OnT'][:,bb],distr['OffT'][:,bb])
                try:
                    rsres = stats.ks_2samp(distr["OnT"][:, bb], distr["OffT"][:, bb])
                except:
                    raise Exception("Problem with the KS 2 sample test...")
                logging.info(rsres)
                ontres = stats.ttest_1samp(distr["OnT"][:, bb], np.zeros((5, 1)))
                logging.info(condit + " " + str(ontres))

            plt.ylim((-30, 50))
            plt.legend()
            plt.title(chann)

    def plot_patient_responses(self):
        Osc_indiv_pop = self.Osc_indiv_pop
        color = self.colors

        for cc, chann in enumerate(["Left", "Right"]):
            plt.figure()
            ax2 = plt.subplot(111)
            distr = nestdict()
            for co, condit in enumerate(["OnT", "OffT"]):
                distr_to_plot = Osc_indiv_pop[chann][condit]

                plt.plot(
                    np.arange(1, 6) + 0.2 * co,
                    distr_to_plot.T,
                    color[co] + ".",
                    markersize=20,
                )

                parts = ax2.violinplot(
                    distr_to_plot,
                    positions=np.array([1, 2, 3, 4, 5]) + 0.2 * co,
                    showmedians=True,
                )
                for partname in ("cbars", "cmins", "cmaxes", "cmedians"):
                    vp = parts[partname]
                    vp.set_edgecolor(color[co])
                    if partname == "cmedians":
                        vp.set_linewidth(5)
                    else:
                        vp.set_linewidth(2)

                for pc in parts["bodies"]:
                    pc.set_facecolor(color[co])
                    pc.set_edgecolor(color[co])
                    # pc.set_linecolor(color[co])

                plt.ylim((-6, 30))
                distr[condit] = distr_to_plot
            plt.title(chann)
            plt.suptitle("Plotting average response for each patient")
            plt.hlines(0, 0, 5, linestyle="dotted", linewidth=5, alpha=0.8)

            for bb in range(5):
                # rsres = stats.ks_2samp(distr['OnT'][:,bb],distr['OffT'][:,bb])
                rsres = stats.ranksums(distr["OnT"][:, bb], distr["OffT"][:, bb])
                # rsres = stats.ks_2samp(distr['OnT'][:,bb],distr['OffT'][:,bb])
                # rsres = stats.wilcoxon(distr['OnT'][:,bb],distr['OffT'][:,bb])
                # rsres = stats.ttest_ind(distr['OnT'][:,bb],distr['OffT'][:,bb])

                # ontres = stats.ranksums(distr['OnT'][:,bb])
                # ontres = stats.kstest(distr['OnT'][:,bb],cdf='norm')
                # ontres = stats.mannwhitneyu(distr['OnT'][:,bb])
                ontres = stats.ttest_1samp(distr["OnT"][:, bb], 0)
                print(DEFAULT_FEAT_ORDER[bb])
                print(rsres)
                print(ontres)
