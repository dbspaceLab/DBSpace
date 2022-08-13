#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 23:27:33 2018

@author: virati
Clinical Vector Class
"""

import json
from collections import defaultdict
import numpy as np
import sys

from dbspace.utils.r_pca import r_pca
import dbspace as dbo
from dbspace.utils.structures import nestdict
from dbspace.utils.stats import pca

import scipy.stats as stats
import scipy.signal as sig
import scipy.io as sio

import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score, auc


class CStruct:
    all_scales = ["HDRS17", "MADRS", "BDI", "GAF"]
    scale_max = {"HDRS17": 40, "MADRS": 50, "BDI": 60, "GAF": -100, "DSC": 0.01}

    def __init__(self, incl_scales=["HDRS17"]):
        self.phase_list = Phase_List("all")
        ClinVect = json.load(
            open("/home/virati/Dropbox/projects/Research/MDD-DBS/Data/ClinVec.json")
        )
        self.pt_list = [ab["pt"] for ab in ClinVect["HAMDs"]]

        depression_dict = nestdict()
        for pp in ClinVect["HAMDs"]:
            for phph, phase in enumerate(self.phase_list):
                for ss, scale in enumerate(self.all_scales):
                    depression_dict[pp["pt"]][phase][scale] = pp[scale][phph]

        self.depr_dict = depression_dict  # This is patient->phase->scale dictionary

        self.normalize_scales()
        self.load_stim_changes()

    """Wraps self.depr_dict to output a patient->scale->phase ARRAY"""

    def normalize_scales(self):
        baseline_values = nestdict()
        for pt in self.pt_list:
            for scale in self.all_scales:
                # make our timeline
                temp_timeline = np.array(
                    [self.depr_dict[pt][phase][scale] for phase in self.phase_list]
                )
                baseline_values[pt][scale] = np.mean(temp_timeline[0:4])

        for pt in self.pt_list:
            for phase in self.phase_list:
                for scale in self.all_scales:
                    # First, we're going to add absolute normalized scales to our dictionary
                    self.depr_dict[pt][phase]["n" + scale] = (
                        self.depr_dict[pt][phase][scale] / self.scale_max[scale]
                    )

                    # now we're going to generate patient-normalized scales, divided by the average of the first four weeks
                    self.depr_dict[pt][phase]["p" + scale] = (
                        self.depr_dict[pt][phase][scale] / baseline_values[pt][scale]
                    )

        # save our baseline values
        self.pt_baseline_depression = baseline_values

    def gen_DSC(self):
        print("Generating DSC Measure")
        allptX = []

        # get phase lookup table
        ph_lut = dbo.all_phases

        # this is the dictionary of optimal decompositions
        opt_lam_dict = defaultdict(dict)
        pt_ll = defaultdict(dict)

        for pp, pat in enumerate(self.pt_list):
            llscore = np.zeros(50)
            pthd = self.get_pt_scale_timeline(pat, "HDRS17") / 30
            ptgaf = self.get_pt_scale_timeline(pat, "GAF") / 100
            ptmd = self.get_pt_scale_timeline(pat, "MADRS")[np.arange(0, 32, 1)] / 45
            ptbdi = self.get_pt_scale_timeline(pat, "BDI") / 60
            # ptmhd = self.get_pt_scale_timeline(pat,'mHDRS')/25

            sX = np.vstack((pthd, ptmd, ptbdi, ptgaf)).T

            # lump it into a big observation vector AS WELL and do the rPCA on the large one later
            allptX.append(sX)

            min_changes = 100
            for ll, lmbda_s in enumerate(np.linspace(0.3, 0.5, 50)):
                # lmbda = 0.33 did very well here
                RPCA = r_pca.R_pca(sX, lmbda=lmbda_s)
                L, S = RPCA.fit()
                Srcomp, Srevals, Srevecs = pca(S)
                Lrcomp, Lrevals, Lrevecs = pca(L)

                # compare sparse component numbers of nonzero
                # derivative is best bet here
                sdiff = np.diff(Srcomp, axis=0)[
                    :, 0
                ]  # grab just the HDRS sparse deviations

                num_changes = np.sum(np.array(sdiff > 0.006).astype(np.int))

                exp_probs = 3
                nchange_diff = np.abs(num_changes - exp_probs)

                if nchange_diff <= min_changes:
                    opt_sparseness = num_changes
                    min_changes = nchange_diff
                    best_lmbda_s = lmbda_s

                # shift_srcomp = Srcomp - np.median(Srcomp,0)
                # llscore[ll] = num_changes[pp] - len(np.where(np.sum(np.abs(shift_srcomp),1) < 1e-6))
            opt_lam_dict[pat] = {
                "Deviation": min_changes,
                "Lambda": best_lmbda_s,
                "Sparseness": opt_sparseness,
            }

            # We have the "optimal" lambda now and we'll do the final rPCA to generate our components
            RPCA = r_pca.R_pca(sX, lmbda=opt_lam_dict[pat]["Lambda"])
            L, S = RPCA.fit()
            Srcomp, Srevals, Srevecs = pca(S)
            Lrcomp, Lrevals, Lrevecs = pca(L)

            # This generates our DSC scores which are just the negative of the mean of the low rank component
            DSC_scores = -np.mean(Lrcomp[:, :], axis=1)

            # This is our OUTPUT and it goes into DSS_dict
            for phph in range(DSC_scores.shape[0]):
                self.depr_dict[pat][ph_lut[phph]]["DSC"] = DSC_scores[phph]

    """Generate our mHDRS scales"""

    def gen_mHDRS(self):
        print("Generating mHDRS")
        ph_lut = dbo.all_phases

        # Cycle through !! THIS USES DSS DICT
        for pat in self.pt_list:
            mhdrs_tser = sig.medfilt(self.get_pt_scale_timeline(pat, "pHDRS17"), 5)
            for phph in range(mhdrs_tser.shape[0]):
                self.depr_dict[pat][ph_lut[phph]]["mHDRS"] = mhdrs_tser[phph]

    def get_binary_depressed(self, pt, phase, scale="HDRS17"):
        return (
            self.depr_dict[pt][phase][scale]
            > self.pt_baseline_depression[pt][scale] / 2
        )

    def get_pt_binary_timeline(self, pt, scale="HDRS17"):
        return np.array(
            [self.get_binary_depressed(pt, phase) for phase in self.phase_list]
        )

    """Return the specific scale value"""

    def get_depression_measure(self, pt, scale, phase):
        return self.depr_dict[pt][phase][scale]

    """Return an array of the patient's depression measuers over the entire phase list"""

    def get_pt_scale_timeline(self, pt, scale):
        return np.array(
            [self.get_depression_measure(pt, scale, phase) for phase in self.phase_list]
        )

    """PLOTTING FUNCTIONS"""

    def plot_pt_timeline(
        self, pts, scale="HDRS17", overlay_binary=False, plot_stim_changes=True
    ):
        plt.figure()
        for pt in pts:
            if pt[0:3] != "DBS":
                raise Exception

            y = np.array(
                [self.depr_dict[pt][phase][scale] for phase in self.phase_list]
            )

            plt.plot(y, alpha=0.8, linewidth=10, label=pt)
            if overlay_binary:
                plt.plot(self.get_pt_binary_timeline(pt), alpha=0.5, linewidth=10)
            if plot_stim_changes:
                stim_changes = np.array(
                    [
                        self.query_stim_change(pt, self.phase_list[pp])
                        for pp in range(self.get_pt_binary_timeline(pt).shape[0])
                    ]
                ).astype(np.int)
                plt.stem(stim_changes)
                # plt.setp(stemlines, 'color', 'red')

            plt.xticks(np.arange(0, 32), self.phase_list, rotation=90)
            plt.title("Plotting " + scale + " for " + pt)
            plt.ylabel(scale + " Value")
            plt.xlabel("Phase")
        plt.legend()

    def query_stim_change(self, pt, ph, include_init=False):
        if pt[0:3] != "DBS":
            pt = "DBS" + pt

        if include_init:
            stim_change_list = self.Stim_Change_Table()
            return (pt, ph) in stim_change_list
        else:
            stim_change_list = self.Stim_Change_Table()
            stim_change_list_rem_init = [
                (aa, bb) for (aa, bb) in stim_change_list if bb[0] != "B"
            ]
            return (pt, ph) in stim_change_list_rem_init

    def load_stim_changes(self):
        # this is where we'll load in information of when stim changes were done so we can maybe LABEL them in figures
        self.stim_change_mat = sio.loadmat("/home/virati/Dropbox/stim_changes.mat")[
            "StimMatrix"
        ]
        # remove the voltage DECREASES?? from DBS905

    def Stim_Change_Table(self):
        # return stim changes in a meaningful format

        diff_matrix = np.hstack(
            (np.diff(self.stim_change_mat) > 0, np.zeros((6, 1)).astype(np.bool))
        )
        # find the phase corresponding to the stim change
        bump_phases = np.array(
            [np.array(dbo.all_phases)[0:][idxs] for idxs in diff_matrix]
        )

        full_table = [
            [(self.pt_list[rr], ph) for ph in row] for rr, row in enumerate(bump_phases)
        ]

        full_table = [item for sublist in full_table for item in sublist]

        # This returns ALL stim change locations

        # TODO do B-- filtering here

        return full_table


""" Main Class for Clinical Data """


class CFrame:
    do_pts = ["901", "903", "905", "906", "907", "908"]
    scale_max = {"HDRS17": 40, "MADRS": 50, "BDI": 60, "GAF": -100, "DSC": 0.01}
    all_scales = ["HDRS17", "MADRS", "BDI", "GAF"]

    # An easy referencer to find a PATIENT's SCALE at PHASE
    lookup = []
    DSS_dict = []  # PT-SCALE-array
    clin_dict = []  # PT-PHASE-dictionary

    def __init__(self, incl_scales=["HDRS17"], norm_scales=False):
        # load in our JSON file
        # Import the data structure needed for the CVect
        ClinVect = json.load(
            open("/home/virati/Dropbox/projects/Research/MDD-DBS/Data/ClinVec.json")
        )

        # Setup the clinical dictionary structure
        clin_dict = defaultdict(dict)
        # This populates the clinical dictionary structure

        for pp in range(len(ClinVect["HAMDs"])):
            ab = ClinVect["HAMDs"][pp]
            clin_dict[ab["pt"]] = defaultdict(dict)
            for phph, phase in enumerate(ClinVect["HAMDs"][pp]["phases"]):
                for ss, scale in enumerate(self.all_scales):
                    clin_dict[ab["pt"]][phase][scale] = ab[scale][phph]
                    if norm_scales:
                        clin_dict[ab["pt"]][phase]["n" + scale] = (
                            ab[scale][phph] / self.scale_max[scale]
                        )
                    else:
                        clin_dict[ab["pt"]][phase][scale] = ab[scale][phph]

        self.OBS_make_dss(ClinVect)

        add_scales = []
        if norm_scales:
            [add_scales.append("n" + scale) for scale in self.all_scales]
        self.all_scales = add_scales
        self.do_scales = incl_scales
        self.clin_dict = clin_dict

        clin_dict = []

        self.omega_state()
        self.derived_measures()
        self.load_stim_changes()

    def OBS_make_dss(self, ClinVect):
        # Setup derived measures
        # THIS IS JUST A COPY PASTE FROM SCALE DYNAMICS, need to merge this in with above so it's all done properly
        DSS_dict = defaultdict(dict)
        # Here, we cycle through each scale and setup an ARRAY
        for ss, scale in enumerate(self.all_scales):
            for pp in range(len(ClinVect["HAMDs"])):
                ab = ClinVect["HAMDs"][pp]

                DSS_dict[ab["pt"]]["n" + scale] = (
                    np.array(ab[scale]) / self.scale_max[scale]
                )
                DSS_dict[ab["pt"]][scale] = np.array(ab[scale])

        # The DSS Dict is THE MOST IMPORTANT DICT IN THE CLASS
        self.DSS_dict = DSS_dict

    """Omega state is the 'final' state that the clinician cares about: are they depressed or not?"""

    def omega_state(self):
        for pt, ph_dict in self.clin_dict.items():
            # find the average of the first 4 weeks
            baseline = ["A04", "A03", "A02", "A01"]
            bl_obs = [ph_dict]
            for phase in ph_dict:
                print(phase)

    """ this is meant to replace self.DSS_dict as a function that calls and manipulates clin_dict """

    def dss_struct(self):
        # go into clin dict and output its contents in a structure that is consistent with dss_dict
        # dss_dict structure is such that [pt][scale][phase x 1?]
        self.clin_array = {
            pt: {scale: {value} for scale in self.all_scales}
            for pt in self.clin_dict.keys()
        }

    """ Here we go through and generate our derived measures from the established clinical scale measures """

    def derived_measures(self):
        self.mHDRS_gen()
        self.dsc_gen()
        self.fake_gen()

    """ here we generate a random set of clinical scales, uniformly random """

    def fake_gen(self):
        for pat in self.do_pts:
            self.DSS_dict["DBS" + pat]["FAKE_good"] = np.random.uniform(
                0, 50, size=(128, 1)
            )
            self.DSS_dict["DBS" + pat]["FAKE_bad"] = np.random.uniform(
                0, 50, size=(128, 1)
            )

    """ here we compute the median HDRS from the HDRS """

    def mHDRS_gen(self):
        print("Generating mHDRS")
        ph_lut = dbo.all_phases

        # Cycle through !! THIS USES DSS DICT
        for pat in self.DSS_dict.keys():
            mhdrs_tser = sig.medfilt(self.DSS_dict[pat]["HDRS17"], 5)
            self.DSS_dict[pat]["mHDRS"] = mhdrs_tser
            for phph in range(mhdrs_tser.shape[0]):
                self.clin_dict[pat][ph_lut[phph]]["mHDRS"] = (
                    mhdrs_tser[phph] / self.scale_max["HDRS17"]
                )

    """ generate the depression state consensus measure """

    def dsc_gen(self):
        print("Generating DSC Measure")
        allptX = []

        # get phase lookup table
        ph_lut = dbo.all_phases

        # Copy our DSS_Dict reference
        big_dict = self.DSS_dict
        # this is the dictionary of optimal decompositions
        opt_lam_dict = defaultdict(dict)
        pt_ll = defaultdict(dict)

        for pp, pat in enumerate(self.do_pts):
            llscore = np.zeros(50)
            pthd = np.array(big_dict["DBS" + pat]["HDRS17"]) / 30
            ptgaf = np.array(big_dict["DBS" + pat]["GAF"]) / 100
            ptmd = np.array(big_dict["DBS" + pat]["MADRS"])[np.arange(0, 32, 1)] / 45
            ptbdi = np.array(big_dict["DBS" + pat]["BDI"]) / 60
            ptmhd = np.array(big_dict["DBS" + pat]["mHDRS"]) / 25

            sX = np.vstack((ptmhd, ptmd, ptbdi, ptgaf)).T

            # lump it into a big observation vector AS WELL and do the rPCA on the large one later
            allptX.append(sX)

            min_changes = 100
            for ll, lmbda_s in enumerate(np.linspace(0.3, 0.5, 50)):
                # lmbda = 0.33 did very well here
                RPCA = r_pca.R_pca(sX, lmbda=lmbda_s)
                L, S = RPCA.fit()
                Srcomp, Srevals, Srevecs = pca(S)
                Lrcomp, Lrevals, Lrevecs = pca(L)

                # compare sparse component numbers of nonzero
                # derivative is best bet here
                sdiff = np.diff(Srcomp, axis=0)[
                    :, 0
                ]  # grab just the HDRS sparse deviations

                num_changes = np.sum(np.array(sdiff > 0.006).astype(np.int))

                exp_probs = 3
                nchange_diff = np.abs(num_changes - exp_probs)

                if nchange_diff <= min_changes:
                    opt_sparseness = num_changes
                    min_changes = nchange_diff
                    best_lmbda_s = lmbda_s

                # shift_srcomp = Srcomp - np.median(Srcomp,0)
                # llscore[ll] = num_changes[pp] - len(np.where(np.sum(np.abs(shift_srcomp),1) < 1e-6))
            opt_lam_dict[pat] = {
                "Deviation": min_changes,
                "Lambda": best_lmbda_s,
                "Sparseness": opt_sparseness,
            }

            # We have the "optimal" lambda now and we'll do the final rPCA to generate our components
            RPCA = r_pca.R_pca(sX, lmbda=opt_lam_dict[pat]["Lambda"])
            L, S = RPCA.fit()
            Srcomp, Srevals, Srevecs = pca(S)
            Lrcomp, Lrevals, Lrevecs = pca(L)

            # This generates our DSC scores which are just the negative of the mean of the low rank component
            DSC_scores = -np.mean(Lrcomp[:, :], axis=1)

            # This is our OUTPUT and it goes into DSS_dict
            self.DSS_dict["DBS" + pat]["DSC"] = DSC_scores / (self.scale_max["DSC"])

            # WTF does the below do...?
            """
            for phph in range(DSC_scores.shape[0]):
                #self.clin_dict[pt][ph_lut[phph]]['DSC']= new_scores[phph]
                self.clin_dict['DBS'+pat][ph_lut[phph]]['DSC'] = DSC_scores[phph]/3
            """

    """this returns to us a big dictionary with all our scales"""

    def c_dict(self):
        clindict = self.clin_dict
        # This will generate a dictionary with each key being a scale, but each value being a matrix for all patients and timepoints
        big_dict = {
            scale: [
                [clindict[pt][week][scale] for week in week_ordered]
                for pt in self.do_pts
            ]
            for scale in self.do_scales
        }
        self.scale_dict = big_dict

    def c_vect(self):
        # each patient will be a dict key
        c_vects = {el: 0 for el in self.do_pts}
        for pp, pt in enumerate(self.do_pts):
            # vector with all clinical measures in the thing
            # return will be phase x clinscores
            c_vect[pt] = 0

    """Get a patient's timecourse of scale"""

    def pt_scale_tcourse(self, pt):
        # return dictionary with all scales, and each element of that dictionary should be a NP vector
        pt_tcourse = {
            rr: self.clin_dict["DBS" + pt][rr] for rr in self.clin_dict["DBS" + pt]
        }
        # 4/7/2020 - OR we can directly go to DSS_dict
        pt_tcourse = self.DSS_dict["DBS" + pt]

        return pt_tcourse

    """PLOTTING FUNCTIONS BELOW"""

    """Plot scale for all patients"""

    def plot_scale(self, scale="HDRS17", pts="all"):
        if pts == "all":
            pts = dbo.all_pts

        plt.figure()
        for patient in pts:
            # pt_tcourse = {rr:self.clin_dict['DBS'+patient][rr][scale] for rr in self.clin_dict['DBS'+patient]}
            pt_tcourse = self.pt_scale_tcourse(patient)
            # now setup the right order
            prop_order = dbo.Phase_List("all")
            # ordered_tcourse = [pt_tcourse[phase][scale] for phase in prop_order]
            # 4/7/2020
            ordered_tcourse = pt_tcourse[scale]

            plt.plot(ordered_tcourse)
            plt.legend(pts)
        plt.title(scale + " for " + str(pts))

    def pr_curve(self, c1, c2):
        pass

    def c_vs_c_plot(self, c1="HDRS17", c2="HDRS17", plot_v_change=True):
        plt.figure()
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

        # do we want to plot the points when the stim was changed?
        phase_list = dbo.Phase_List("all")
        if plot_v_change:
            stim_change_list = self.Stim_Change_Table()

        big_vchange_list = []

        # This is important for AUC calculations on . Decoupled from *plotting*
        include_bs = False
        if not include_bs:
            start_delay = 8
        else:
            start_delay = 0

        for pat in self.do_pts:
            scale1 = np.array(self.DSS_dict["DBS" + pat][c1][start_delay:32])
            scale2 = np.array(self.DSS_dict["DBS" + pat][c2][start_delay:32])

            ax = plt.subplot(1, 2, 1)
            # plot the A& B periods
            if include_bs:
                plt.scatter(
                    scale1[0:8], scale2[0:8], alpha=0.2, color="black", marker="s"
                )
                # Plot the C periods
                plt.scatter(scale1[8:], scale2[8:], alpha=0.2, color="blue")
            else:
                plt.scatter(scale1, scale2, alpha=0.2, color="blue")

            # plot the changes for the patient
            if include_bs:
                phases_v_changed = [b for a, b in stim_change_list if a == pat]
            else:
                phases_v_changed = [
                    b for a, b in stim_change_list if a == pat and b != "B04"
                ]

            phase_idx_v_changed = (
                np.array([phase_list.index(b) for b in phases_v_changed]) - start_delay
            )

            # Plot the observations with stim changes
            plt.scatter(
                scale1[phase_idx_v_changed],
                scale2[phase_idx_v_changed],
                marker="^",
                s=130,
                alpha=0.3,
                color="red",
            )

            for ii, idx in enumerate(phase_idx_v_changed):
                plt.annotate(
                    phases_v_changed[ii] + " " + pat,
                    (scale1[idx], scale2[idx]),
                    fontsize=8,
                    color="gray",
                )
                # plt.annotate('test',(1,1),fontsize=8,color='gray')

                #

            change_vec = np.zeros_like(scale1)
            change_vec[phase_idx_v_changed] = 1

            big_vchange_list.append((scale1, scale2, change_vec))

        plt.xlabel(c1)
        plt.ylabel(c2)

        # Correlation measures
        corr_matr = np.array(
            [
                (
                    self.DSS_dict["DBS" + pat][c1][0:32],
                    self.DSS_dict["DBS" + pat][c2][0:32],
                )
                for pat in self.do_pts
            ]
        )
        corr_matr = np.swapaxes(corr_matr, 0, 1)
        corr_matr = corr_matr.reshape(2, -1, order="C")

        spearm = stats.spearmanr(corr_matr[0, :], corr_matr[1, :])
        pears = stats.pearsonr(corr_matr[0, :], corr_matr[1, :])

        print("SpearCorr between " + c1 + " and " + c2 + " is: " + str(spearm))
        print("PearsCorr between " + c1 + " and " + c2 + " is: " + str(pears))

        # plt.plot([-1,60],[-1,60])
        # plt.axes().set_aspect('equal')
        # plt.legend(self.do_pts)

        # should be 6x3x32
        self.big_v_change_list = (
            np.array(big_vchange_list).swapaxes(0, 1).reshape(3, -1, order="C")
        )

        scale_labels = (c1, c2, "Min")
        ax2 = plt.subplot(1, 2, 2)
        for ii in range(2):
            # now do the AUC curves and P-R curves
            precision, recall, _ = precision_recall_curve(
                self.big_v_change_list[2, :], self.big_v_change_list[ii, :]
            )
            # Compute AUC directly from pr
            prauc = auc(recall, precision)
            prauc = np.sum(precision) / recall.shape[0]
            # Compute average precision
            avg_precision = average_precision_score(
                self.big_v_change_list[2, :],
                self.big_v_change_list[ii, :],
                average="micro",
            )
            plt.plot(recall, precision)
            # plt.subplot(2,1,2)
            # plt.plot(recall,precision)
            # plt.annotate('Average precision for ' + str(scales[ii]) + ': ' + str(avg_precision)  + ' AUC: ' + str(prauc),(-2,2-(ii/4)),fontsize=8)
            ax.text(
                0.1,
                0.95 - ii / 4,
                "AvgPrec "
                + str(scale_labels[ii])
                + ": "
                + str(avg_precision)
                + " \nAUC: "
                + str(prauc),
                transform=ax.transAxes,
                fontsize=14,
                verticalalignment="top",
                bbox=props,
            )

        plt.ylim((0, 1))
        # This does some other algorithm that runs on just the clinical scales
        # ## do the derived algorithms now
        # ii=2
        # min_algo = np.max(np.vstack((self.big_v_change_list[0,:],self.big_v_change_list[1,:])),axis=0)
        # precision,recall,_ = precision_recall_curve(self.big_v_change_list[2],min_algo)
        # plt.plot(recall,precision)

        # prauc = auc(recall,precision)
        # prauc = np.sum(precision) / recall.shape[0]
        # avg_precision = average_precision_score(self.big_v_change_list[2],min_algo,average="micro")
        # #plt.annotate('Average precision for ' + str(scales[ii]) + ': ' + str(avg_precision) + ' AUC: ' + str(prauc),(-2,1),fontsize=8)
        # ax.text(0.1, 0.95 - 3/4, 'AvgPrec ' + str(scale_labels[ii]) + ': ' + str(avg_precision)  + ' \nAUC: ' + str(prauc), transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)

    def load_stim_changes(self):
        # this is where we'll load in information of when stim changes were done so we can maybe LABEL them in figures
        self.stim_change_mat = sio.loadmat("/home/virati/Dropbox/stim_changes.mat")[
            "StimMatrix"
        ]
        # remove the voltage DECREASES?? from DBS905

    def Stim_Change_Table(self):
        # return stim changes in a meaningful format

        # Diff vector belongs in first part of the diff_matrix
        # Key thing to check for: CHanges are in B04, and DBS907 change is at C15
        # see: https://docs.google.com/spreadsheets/d/1HLZfMoE83ulHm0dc3j8c3ZEDk4LaF-0qQztnavgmAQw/edit#gid=0

        diff_matrix = np.hstack(
            (np.diff(self.stim_change_mat) > 0, np.zeros((6, 1)).astype(np.bool))
        )
        # find the phase corresponding to the stim change
        bump_phases = np.array(
            [np.array(dbo.all_phases)[0:][idxs] for idxs in diff_matrix]
        )

        full_table = [
            [(self.do_pts[rr], ph) for ph in row] for rr, row in enumerate(bump_phases)
        ]

        full_table = [item for sublist in full_table for item in sublist]

        # This returns ALL stim change locations

        # TODO do B-- filtering here

        return full_table

    def week_labels(self):
        week_labels = ["A04", "A03", "A02", "A01", "B01", "B02", "B03", "B04"]
        for ii in range(24):
            if ii < 10:
                ii_label = "0" + str(ii)
            else:
                ii_label = str(ii)

            week_labels.append("C" + ii_label)
        return week_labels

    """ Get the min and max weeks for each patient and the scale associated with that week"""

    def min_max_weeks(self):
        week_labels = self.week_labels()
        hdrs_info = nestdict()

        for pt in self.do_pts:
            pt_hdrs_traj = [a for a in self.DSS_dict["DBS" + pt]["HDRS17"]][8:]

            hdrs_info[pt]["max"]["index"] = np.argmax(pt_hdrs_traj)
            hdrs_info[pt]["min"]["index"] = np.argmin(pt_hdrs_traj)
            hdrs_info[pt]["max"]["week"] = week_labels[np.argmax(pt_hdrs_traj) + 8]
            hdrs_info[pt]["min"]["week"] = week_labels[np.argmin(pt_hdrs_traj) + 8]

            hdrs_info[pt]["max"]["HDRSr"] = pt_hdrs_traj[hdrs_info[pt]["max"]["index"]]
            hdrs_info[pt]["min"]["HDRSr"] = pt_hdrs_traj[hdrs_info[pt]["min"]["index"]]
            hdrs_info[pt]["traj"]["HDRSr"] = pt_hdrs_traj

        return hdrs_info

    """
    TODO
    Here we make a scatter plot that plots MEAS 1 vs MEAS 2, with the stim changes labeled
    """

    def plot_meas_vs_meas_Vchanges(self):
        pass


""" Unit Test for CFrame """
if __name__ == "__main__":
    TestStruct = CStruct()
    TestStruct.gen_DSC()
    TestStruct.gen_mHDRS()
    TestStruct.plot_pt_timeline(
        ["DBS901", "DBS903", "DBS905", "DBS906", "DBS907", "DBS908"],
        "DSC",
        overlay_binary=False,
    )

    # binarization_test = TestStruct.get_pt_binary_timeline('DBS901')
    # plt.plot(binarization_test)


def plot_c_vs_c():
    TestFrame = CFrame(norm_scales=False)
    for c2 in ["mHDRS", "GAF", "BDI", "MADRS", "DSC"]:
        TestFrame.c_vs_c_plot(c1="HDRS17", c2=c2)
    TestFrame.plot_scale(scale="DSC")
    TestFrame.plot_scale(scale="HDRS17")
    TestFrame.c_dict()
    plt.show()


def Phase_List(exprs="all", nmo=-1):
    """
    Phase list for MaybergLab DBS study - patients DBS901-908

    """
    all_phases = ["A04", "A03", "A02", "A01", "B01", "B02", "B03", "B04"]
    for aa in range(1, 25):
        if aa < 10:
            numstr = "0" + str(aa)
        else:
            numstr = str(aa)
        all_phases.append("C" + numstr)

        ephys_phases = all_phases[4:]
    if exprs == "all":
        return all_phases
    elif exprs == "ephys":
        return ephys_phases
    elif exprs == "Nmo_ephys":
        # nmo = 3
        return ephys_phases[0 : 4 * (nmo + 1) - 1]
    elif exprs == "Nmo_onStim":
        # nmo = 5
        return ephys_phases[4 : 4 * (nmo + 1) - 1]
