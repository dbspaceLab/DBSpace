#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:47:25 2020

@author: virati
NEW classes for readout training, testing, and validation
"""
import random

import matplotlib.cm as cm
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy import interp
from sklearn import metrics
from sklearn.linear_model import ElasticNet, RidgeCV, LassoCV
from dbspace.signal.oscillations import poly_subtr
from sklearn.metrics import (
    auc,
    average_precision_score,
    mean_absolute_error,
    mean_squared_error,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle

np.random.seed(seed=2011)
random.seed(2011)

import dbspace as dbo
from dbspace.utils.structures import nestdict
from sklearn import linear_model

default_params = {"CrossValid": 10}

import seaborn as sns

sns.set_context("paper")

sns.set(font_scale=4)
sns.set_style("white")

import copy
import itertools

from dbspace.signal.oscillations import DEFAULT_FEAT_ORDER
from dbspace.readout.ClinVect import Phase_List
import dbspace.signal.oscillations as dbo

def zero_mean(inp):
    return inp - np.mean(inp)


#%%
class base_decoder:
    # Parent readout class
    # Some very fixed constants up here
    circ = "day"
    ch_num = 2
    global_plotting = False

    def __init__(
        self, *args, **kwargs
    ):  # BRFrame,ClinFrame,pts,clin_measure='HDRS17'):
        self.YFrame = kwargs["BRFrame"]
        self.CFrame = kwargs["ClinFrame"]
        self.pts = kwargs["pts"]
        self.c_meas = kwargs["clin_measure"]
        self.fvect = self.YFrame.data_basis["F"]
        self.do_shuffle_null = kwargs["shuffle_null"]

        # here we decide which features we want to do for this analysis
        if kwargs["FeatureSet"] == "stim_check":
            self.do_feats = ["Delta", "Theta", "Alpha", "Beta*", "Gamma1", "THarm"]
        elif kwargs["FeatureSet"] == "main":
            self.do_feats = DEFAULT_FEAT_ORDER

        self.feat_labels = ["L" + feat for feat in self.do_feats] + [
            "R" + feat for feat in self.do_feats
        ]

        self.regression_algo = linear_model.LinearRegression

    """Filter out the recordings we want"""

    def filter_recs(self, rec_class="main_study"):
        if rec_class == "main_study":
            filter_phases = Phase_List(exprs="ephys")
            self.active_rec_list = [
                rec
                for rec in self.YFrame.file_meta
                if rec["Phase"] in filter_phases
                and rec["Patient"] in self.pts
                and rec["Circadian"] in self.circ
            ]

        self.filter_phases = filter_phases

    def y_c_pair(self, rec_list):
        scale_lookup = self.CFrame.clin_dict

        # self.data = [(rec,scale_lookup[pt][phase]['nHDRS'] for rec in rec_list if rec]

    """ Plot things we care about when it comes to how many recordings each patient x phase has, etc."""

    def rec_set_size(self):
        filter_phases = Phase_List(exprs="ephys")
        accounting = np.zeros((len(self.pts), len(filter_phases)))
        detailed_dict = nestdict()

        for pp, pt in enumerate(self.pts):

            print(
                pt
                + " has "
                + str(
                    len(
                        [
                            rec
                            for rec in self.YFrame.file_meta
                            if rec["Phase"] in filter_phases and rec["Patient"] == pt
                        ]
                    )
                )
                + " recordings"
            )
            for ph, phase in enumerate(filter_phases):

                detailed_dict[pt][phase] = [
                    rec
                    for rec in self.YFrame.file_meta
                    if rec["Phase"] == phase and rec["Patient"] == pt
                ]
                print(
                    pt
                    + " has "
                    + str(
                        len(
                            [
                                rec
                                for rec in self.YFrame.file_meta
                                if rec["Phase"] == phase and rec["Patient"] == pt
                            ]
                        )
                    )
                    + " recordings in Phase "
                    + phase
                )

                accounting[pp, ph] = len(detailed_dict[pt][phase])

        # Plot the accounting
        if self.global_plotting:
            plt.figure()
            plt.imshow(accounting)
            plt.figure()
            plt.plot(accounting[:, :].T)

    """Plot PSDs for the first N recordings, sanity check"""

    def plot_psds(self, upper_lim=10):
        plt.figure()
        for ii in range(upper_lim):
            plt.subplot(121)
            plt.plot(
                np.linspace(0, 211, 513), np.log10(self.train_set[ii]["Data"]["Left"])
            )
            plt.subplot(122)
            plt.plot(np.log10(self.train_set[ii]["Data"]["Right"]))

    """split out our training and validation set recordings"""

    def split_train_set(self, train_ratio=0.6):
        self.train_set, self.test_set = train_test_split(
            self.active_rec_list, train_size=train_ratio, shuffle=True
        )

    """Setup our data for training"""

    def train_setup(self):
        self.train_set_y, self.train_set_c = self.calculate_states_in_set(
            self.train_set
        )

    def shuffle_test_c(self):
        np.random.shuffle(self.test_set_c)

    """ Train our model"""

    def train_model(self, do_null=False):
        if do_null:
            shuffled_c = copy.deepcopy(self.train_set_c)
            np.random.shuffle(shuffled_c)
            self.decode_model = self.regression_algo().fit(self.train_set_y, shuffled_c)
        else:
            self.decode_model = self.regression_algo().fit(
                self.train_set_y, self.train_set_c
            )

    """See what the null model generates for stats"""

    def model_analysis(self, do_null=False, n_iter=1, do_plot=False):
        self.train_setup()
        self.test_setup()
        null_stats = []

        for ii in range(n_iter):
            self.train_model(do_null=do_null)

            _, stats = self.test_model()
            null_stats.append(stats)

        slope_results = np.array([a["Slope"] for a in null_stats])
        r2_results = np.array([a["Score"] for a in null_stats])

        if self.global_plotting and do_plot:
            # plot our distribution
            plt.figure()
            plt.hist(slope_results, bins=10)

        return slope_results, r2_results

    """Plot the coefficient path in the regression"""

    def plot_coeff_sig_path(self, do_plot=False):
        # print('Running path')
        offset_train_y = self.train_set_y - np.mean(self.train_set_y)
        offset_train_c = self.train_set_c - np.mean(self.train_set_c)

        coeff_path = self.decode_model.path(
            offset_train_y,
            offset_train_c,
            n_alphas=100,
            cv=False,
            eps=0.001,
            fit_intercept=True,
        )  # Path is STUPIDLY hardcoded to do fit_intercept = False
        if self.global_plotting and do_plot:
            plt.figure()
            plt.subplot(211)
            for ii, label in enumerate(self.feat_labels):
                plt.plot(
                    -np.log(coeff_path[0]),
                    coeff_path[1].squeeze()[ii, :],
                    linewidth=5,
                    label=label,
                )
            # plt.legend(labels = self.feat_labels)
            plt.legend()
            plt.subplot(212)
            plt.plot(self.decode_model.coef_)
            plt.hlines(0, -2, 10, linestyle="dotted")
            plt.xticks(np.arange(10), self.feat_labels)

    """setup our data for the TESTING"""

    def test_setup(self):
        self.test_set_y, self.test_set_c = self.calculate_states_in_set(self.test_set)
        if self.do_shuffle_null:
            self.shuffle_null()

    """Main TESTING method for our model"""

    def test_model(self):
        predicted_c = self.decode_model.predict(self.test_set_y)
        test_stats = self.get_test_stats(self.test_set_y, self.test_set_c, predicted_c)

        return predicted_c, test_stats

    """Plot our predictions here"""

    def plot_test_predictions(self):
        predicted_c, _ = self.test_model()
        plt.figure()
        plt.plot(self.test_set_c, predicted_c, "r.")
        plt.title("Predicted vs Actual")
        # plt.plot(predicted_c,predicted_c - self.test_set_c,'r.');plt.title('Residuals')
        plt.plot([0, 1], [0, 1])
        plt.xlim((0, 1.1))
        plt.ylim((0, 1.1))

    def get_test_stats(self, test_y, true_c, predicted_c):
        # Pearson
        test_y = test_y.squeeze()
        true_c = true_c.squeeze()
        predicted_c = predicted_c.squeeze()
        p_stats = stats.pearsonr(predicted_c, true_c)
        # Spearman
        s_stats = stats.spearmanr(predicted_c, true_c)

        # Linear Regression
        regr_model = linear_model.LinearRegression().fit(
            predicted_c.reshape(-1, 1), true_c.reshape(-1, 1)
        )  # true_c.reshape(-1,1),predicted_c.reshape(-1,1))
        lr_slope = regr_model.coef_[0]
        # Robust regression
        # ransac = linear_model.RANSACRegressor().fit(true_c.reshape(-1,1),predicted_c.reshape(-1,1))
        # ransac_slope = ransac.estimator_.coef_

        stat_dict = {
            "Score": self.decode_model.score(test_y, true_c),
            "Pearson": p_stats,
            "Spearman": s_stats,
            "Slope": lr_slope,
        }  # ,'RANSACm':ransac_slope}
        return stat_dict

    """Plot the test statistics"""

    def plot_test_stats(self):
        predicted_c = self.test_model()

        plt.scatter(self.test_set_c, predicted_c)

        print(corr)

    """Plot the regression visualization of the test procedure"""

    def plot_test_regression(self):
        # do a final test on *all* the data for plotting purposes
        predicted_c = self.decode_model.predict(self.test_set_y)
        r2score = self.decode_model.score(self.test_set_y, self.test_set_c)
        mse = mean_squared_error(self.test_set_c, predicted_c)
        corr = stats.pearsonr(self.test_set_c.squeeze(), predicted_c.squeeze())

        plt.plot([0, 1], [0, 1], color="gray", linestyle="dotted")
        ax = sns.regplot(x=self.test_set_c, y=predicted_c)
        plt.title(
            "R^2:" + str(r2score) + "\n" + " MSE:" + str(mse) + "\n Corr:" + str(corr)
        )
        plt.xlim((0, 1))
        plt.ylim((0, 1))

    """ Calculate oscillatory states for a set of recordings"""

    def calculate_states_in_set(self, data_set):
        state_vector = []
        depr_vector = []
        # pt_ph_vector = []

        for rr in data_set:
            psd_poly_done = {
                ch: poly_subtr(
                    fvect=self.fvect, input_psd=rr["Data"][ch], polyord=5
                )[0]
                for ch in rr["Data"].keys()
            }

            feat_vect = np.zeros(shape=(len(self.do_feats), self.ch_num))
            for ff, featname in enumerate(self.do_feats):
                dofunc = dbo.FEAT_DICT[featname]
                feat_calc = dofunc["fn"](psd_poly_done, self.fvect, dofunc["param"])
                feat_vect[ff, :] = np.array([feat_calc[ch] for ch in ["Left", "Right"]])

            # We need to flatten the state between channels...
            # Then we go ahead and append it to the state vector
            state_vector.append(
                np.reshape(feat_vect, -1, order="F")
            )  # we want our FEATURE index to change quickest so we go (0,0) -> (1,0) -> (2,0) -> ... (4,1)

            # now we need to get a vector of the clinical states
            depr_value = self.CFrame.get_depression_measure(
                "DBS" + rr["Patient"], self.c_meas, rr["Phase"]
            )
            depr_vector.append(depr_value)

            # pt_ph_vector.append(rr['Patient'],rr['Phase'])
        return np.array(state_vector), np.array(depr_vector)  # , pt_ph_vector

    def OBSget_coeffs(self):
        model = self.decode_model
        active_coeffs = self.decode_model.coef_

        return active_coeffs

    """Plot coefficients of our model"""

    def plot_decode_coeffs(self, model):
        active_coeffs = np.array(model.coef_).squeeze()
        # plt.subplot(1,2,side+1)
        plt.figure()
        plt.plot(active_coeffs)
        plt.hlines(0, -2, 10, linestyle="dotted")
        plt.vlines(5, -1, 1, linestyle="solid", color="blue")
        plt.ylim(
            (
                -np.max(np.abs(active_coeffs)) + 0.01,
                np.max(np.abs(active_coeffs)) - 0.01,
            )
        )
        plt.xlim((-1, 10))
        plt.xticks(np.arange(10), self.feat_labels)

    def plot_test_ensemble(self):
        plt.figure()
        plt.subplot(211)
        self.plot_test_regression()
        plt.subplot(212)
        self.plot_decode_coeffs()


#%%
class var_decoder(base_decoder):
    variance_analysis = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def train_setup(self):
        print("Performing Training Setup for Weekly Decoder")

        (
            self.train_set_y,
            self.train_set_c,
            self.train_set_pt,
            self.train_set_ph,
        ) = self.aggregate_weeks(self.train_set)

    def test_setup(self):
        print("Performing TESTING Setup for Weekly Decoder")

        (
            self.test_set_y,
            self.test_set_c,
            self.test_set_pt,
            self.test_set_ph,
        ) = self.aggregate_weeks(self.test_set)
        if self.do_shuffle_null:
            self.shuffle_test_c()

    def train_model(self, do_null):
        self.decode_model = self.regression_algo().fit(
            self.train_set_y, self.train_set_c
        )
        # print('Alpha: ' + str(self.decode_model.alpha_) + ' | L1r: ' + str(self.decode_model.l1_ratio_))
        # self.plot_decode_coeffs(self.decode_model)

    def test_setup(self):
        print("Performing TESTING Setup for Weekly Decoder")

        (
            self.test_set_y,
            self.test_set_c,
            self.test_set_pt,
            self.test_set_ph,
        ) = self.aggregate_weeks(self.test_set)
        if self.do_shuffle_null:
            self.shuffle_test_c()

    def aggregate_weeks(self, dataset):
        # print('Performing Training Setup for Weekly Decoder')
        # go through our training set and aggregate every recording within a given week
        # train_set_y,train_set_c = self.calculate_states_in_set(self.train_set)

        running_list = []
        for pt in self.pts:
            for phase in self.filter_phases:
                block_set = [
                    rr for rr in dataset if rr["Patient"] == pt and rr["Phase"] == phase
                ]
                if block_set != []:
                    y_set, c_set = self.calculate_states_in_set(block_set)
                    weekly_y_set = np.mean(y_set, axis=0)

                    running_list.append(
                        (weekly_y_set, c_set[0], pt, phase)
                    )  # all the c_set values should be the exact same

        y_state = np.array(
            [a for (a, b, c, d) in running_list]
        )  # outputs ~168 observed weeks x 10 features
        c_state = np.array([b for (a, b, c, d) in running_list]).reshape(
            -1, 1
        )  # outputs ~168 observed weeks
        pt_name = np.array([c for (a, b, c, d) in running_list])
        phase_label = np.array([d for (a, b, c, d) in running_list])

        return y_state.reshape(-1, 1), c_state.reshape(-1, 1), pt_name, phase_label


#%%
class weekly_decoder(base_decoder):
    variance_analysis = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if kwargs["algo"] == "ENR":
            self.regression_algo = ElasticNet(alpha=0.2)
        elif kwargs["algo"] == "ENR_all":
            self.regression_algo = ElasticNet
            # self.model_args = {'alpha':np.e **-3.4,'l1_ratio':0.8} #5/28, good stats, visual is small
            self.model_args = {"alpha": [], "l1_ratio": 0.8}
        elif kwargs["algo"] == "Ridge":
            self.regression_algo = RidgeCV()
        elif kwargs["algo"] == "Lasso":
            self.regression_algo = LassoCV()

        if kwargs["variance"] == True:
            self.variance_analysis = True

    def train_model(self):
        self.decode_model = self.regression_algo(**self.model_args).fit(
            self.train_set_y, self.train_set_c
        )
        # print('Alpha: ' + str(self.decode_model.alpha_) + ' | L1r: ' + str(self.decode_model.l1_ratio_))
        # self.plot_decode_coeffs(self.decode_model)

    def aggregate_weeks(self, dataset):
        # print('Performing Training Setup for Weekly Decoder')
        # go through our training set and aggregate every recording within a given week
        # train_set_y,train_set_c = self.calculate_states_in_set(self.train_set)

        running_list = []
        for pt in self.pts:
            for phase in self.filter_phases:
                block_set = [
                    rr for rr in dataset if rr["Patient"] == pt and rr["Phase"] == phase
                ]
                if block_set != []:
                    y_set, c_set = self.calculate_states_in_set(block_set)
                    if not self.variance_analysis:
                        weekly_y_set = np.mean(y_set, axis=0)
                    else:
                        weekly_y_set = np.var(y_set, axis=0)

                    running_list.append(
                        (weekly_y_set, c_set[0], pt, phase)
                    )  # all the c_set values should be the exact same

        y_state = np.array(
            [a for (a, b, c, d) in running_list]
        )  # outputs ~168 observed weeks x 10 features
        c_state = np.array([b for (a, b, c, d) in running_list]).reshape(
            -1, 1
        )  # outputs ~168 observed weeks
        pt_name = np.array([c for (a, b, c, d) in running_list])
        phase_label = np.array([d for (a, b, c, d) in running_list])

        return y_state, c_state, pt_name, phase_label

    def train_setup(self):
        print("Performing Training Setup for Weekly Decoder")

        (
            self.train_set_y,
            self.train_set_c,
            self.train_set_pt,
            self.train_set_ph,
        ) = self.aggregate_weeks(self.train_set)

    def test_setup(self):
        print("Performing TESTING Setup for Weekly Decoder")

        (
            self.test_set_y,
            self.test_set_c,
            self.test_set_pt,
            self.test_set_ph,
        ) = self.aggregate_weeks(self.test_set)
        if self.do_shuffle_null:
            self.shuffle_test_c()

    """This method goes down the regression path and assesses the performance of the model all along the way: Returns the 'optimal' alpha"""

    def _path_slope_regression(
        self, do_plot=False, suppress_vars=0.2, override_alpha=False
    ):
        assess_traj = []

        (
            internal_train_y,
            internal_test_y,
            internal_train_c,
            internal_test_c,
        ) = train_test_split(
            self.train_set_y, self.train_set_c, train_size=0.6, shuffle=True
        )
        path_model = ElasticNet(
            l1_ratio=0.8, fit_intercept=True, normalize=False
        )
        path = path_model.path(
            zero_mean(self.train_set_y),
            zero_mean(self.train_set_c),
            eps=0.0001,
            n_alphas=1000,
        )
        for alpha in path[0]:
            run_model = ElasticNet(
                alpha=alpha, l1_ratio=0.8, fit_intercept=True, normalize=False
            )
            run_model.fit(internal_train_y, internal_train_c)
            # lin regression to identify slope
            predict_c = run_model.predict(internal_test_y)
            score = run_model.score(internal_test_y, internal_test_c)
            slope = stats.linregress(
                internal_test_c.squeeze(), predict_c
            )  # THIS IS BACKWARDS TO AVOID NAN #this used to be linregress(actual,predicted) which is I believe identical to the R^2 of the (predicted,actual) and reflects the percentage of the variance explained
            assess_traj.append({"Alpha": alpha, "Slope": slope, "Score": score})
        # now do the path, this should match up with above
        self._path_slope_results = assess_traj, path

        # Figure out how many coefficients are around
        coeff_present = (np.abs(path[1].squeeze().T) > 0).astype(np.int)
        total_coeffs = np.sum(coeff_present, axis=1)
        slope_traj_vec = np.array([a["Slope"][0] for a in assess_traj])
        score_traj_vec = np.array([a["Score"] for a in assess_traj])

        ## Find the max of the r^2 for our optimal alpha
        optimal_alpha = path[0][
            np.argmax(score_traj_vec)
        ]  # This merely does an R^2 optimal

        # IF YOU WANT TO ALSO SUPPRESS number of variables
        lamb = suppress_vars
        lamb2 = 10
        # optimal_alpha = path[0][np.argmax(lamb2*score_traj_vec - lamb * total_coeffs )] #REMOVED SLOPE FROM THIS

        print("Optimal Alpha: ", optimal_alpha)
        if override_alpha:
            print("But Overriding Alpha with ", override_alpha)
            optimal_alpha = override_alpha

        # optimal_alpha = path[0][np.argmin(np.abs(slope_traj_vec - 1))]

        if self.global_plotting and do_plot:
            fig, ax1 = plt.subplots()

            ax1.plot(-np.log(path[0]), slope_traj_vec, label="Slope")
            plt.title("Slope of readout, Score of readout")
            ax1.plot(-np.log(path[0]), score_traj_vec, label="Score")
            plt.legend()

            plt.vlines(-np.log(optimal_alpha), 0, 0.3, linewidth=10)
            ax2 = ax1.twinx()
            ax2.plot(-np.log(path[0]), total_coeffs)

            fig, ax1 = plt.subplots()
            colors = pl.cm.jet(np.linspace(0, 1, path[1].squeeze().T.shape[1]))

            for ii in range(path[1].squeeze().T.shape[1]):
                ax1.plot(-np.log(path[0]), path[1].squeeze().T[:, ii], color=colors[ii])
                plt.title("Regularization Path")
            ax1.legend(labels=self.feat_labels)
            ax2 = ax1.twinx()

            ax2.plot(-np.log(path[0]), total_coeffs)
            plt.vlines(-np.log(optimal_alpha), 0, 0.3, linewidth=10)

        return optimal_alpha


#%%
class weekly_decoderCV(weekly_decoder):
    def __init__(self, *args, **kwargs):
        print("Initialized the Weekly CV decoder")
        super().__init__(*args, **kwargs)

        if kwargs["algo"] == "ENR":

            self.regression_algo = ElasticNet
            # self.regression_algo = linear_model.Lasso
            # self.model_args = {'alpha':np.e **-3.4,'l1_ratio':0.8} #5/28, good stats, visual is small
            self.model_args = {"alpha": np.e ** kwargs["alpha"], "l1_ratio": 0.8}
            print("Running ENR_CV w/:" + str(self.model_args["l1_ratio"]))

            # BEFORE 5/28
            # self.regression_algo = linear_model.ElasticNetCV
            # self.model_args = {'alphas':np.linspace(0.01,0.04,20),'l1_ratio':np.linspace(0.1,0.3,10),'cv':10}

        self.pt_CV_sets(n=3)

    def pt_CV_sets(self, n=3):
        pt_combos = list(itertools.combinations(self.pts, n))

        self.CV_num_combos = len(pt_combos)
        self.CV_pt_combos = pt_combos

    def train_setup(self):
        print("Performing Training Setup for Weekly Decoder")
        (
            self.train_set_y,
            self.train_set_c,
            self.train_set_pt,
            self.train_set_ph,
        ) = self.aggregate_weeks(self.train_set)

        self.model_args["alpha"] = self._path_slope_regression(do_plot=True)
        print("Set ENR-Alpha at " + str(self.model_args["alpha"]))

    """ Train our model"""

    def train_model(self):
        # Our first goal is to learn a model for each patient combination
        decode_model_combos = [None] * self.CV_num_combos
        model_performance_combos = [None] * self.CV_num_combos
        coeff_path = [None] * self.CV_num_combos

        for run, pt_combo in enumerate(self.CV_pt_combos):
            print(pt_combo)
            combo_train_y = [
                a
                for (a, c) in zip(self.train_set_y, self.train_set_pt)
                if c in pt_combo
            ]
            combo_train_c = [
                b
                for (b, c) in zip(self.train_set_c, self.train_set_pt)
                if c in pt_combo
            ]

            decode_model_combos[run] = self.regression_algo(**self.model_args).fit(
                combo_train_y, combo_train_c
            )

            offset_train_y = combo_train_y - np.mean(combo_train_y)
            offset_train_c = combo_train_c - np.mean(combo_train_c)

            coeff_path[run] = decode_model_combos[run].path(
                offset_train_y,
                offset_train_c,
                n_alphas=100,
            )

            combo_test_y = [
                a
                for (a, c) in zip(self.train_set_y, self.train_set_pt)
                if c not in pt_combo
            ]
            combo_test_c = [
                b
                for (b, c) in zip(self.train_set_c, self.train_set_pt)
                if c not in pt_combo
            ]

            # quick coeff path
            # model_performance_combos[run] = decode_model_combos[run].score(combo_test_y,combo_test_c)
            # pred_c = decode_model_combos[run].predict(combo_test_y)
            # model_performance_combos[run] = mean_absolute_error(combo_test_c,pred_c)

        self.decode_model_combos_ = decode_model_combos
        self.decode_model_combos_paths_ = coeff_path

        average_model_coeffs, _ = self.get_average_model(self.decode_model_combos_)
        self.decode_model = linear_model.LinearRegression()
        self.decode_model.coef_ = average_model_coeffs
        self.decode_model.intercept_ = np.mean(
            [m.intercept_ for m in self.decode_model_combos_]
        )

    """Plot the paths for all CV combos"""

    def plot_combo_paths(self, do_feats=[]):
        if do_feats == []:
            do_feats = self.feat_labels
        plt.figure()

        for ii, label in enumerate(self.feat_labels):
            if label in do_feats:
                path_list = np.array(
                    [path[1].squeeze() for path in self.decode_model_combos_paths_]
                )
                for path in self.decode_model_combos_paths_:
                    plt.plot(
                        -np.log(path[0]),
                        path[1].squeeze()[ii, :],
                        linewidth=5,
                        alpha=0.1,
                    )
                plt.plot(
                    -np.log(path[0]),
                    np.mean(path_list, axis=0)[ii, :],
                    linewidth=10,
                    label=label,
                )
            # plt.legend(labels = self.feat_labels)
            plt.legend()

    def get_average_model(self, model):
        active_coeffs = []
        for ii in self.decode_model_combos_:
            active_coeffs.append([ii.coef_])

        active_coeffs = np.array(active_coeffs).squeeze()
        average_model = np.median(active_coeffs, axis=0)
        # average_model = np.zeros(shape=active_coeffs.shape)

        # return the average model with the stats for each coefficient
        return average_model, stats

    def test_model(self):
        ensemble_score = []
        ensemble_corr = []
        self.test_stats = (
            []
        )  # {'Prediction Score': [], 'Pearson Corr Score': [], 'Spearman Corr Score': []}
        display_test = []

        for tt in range(100):
            test_subset_y, test_subset_c, test_subset_pt, test_subset_ph = zip(
                *random.sample(
                    list(
                        zip(
                            self.test_set_y,
                            self.test_set_c,
                            self.test_set_pt,
                            self.test_set_ph,
                        )
                    ),
                    np.ceil(0.5 * len(self.test_set_y)).astype(np.int),
                )
            )
            test_subset_y = np.array(test_subset_y)
            test_subset_c = np.array(test_subset_c)

            predicted_c = self.decode_model.predict(test_subset_y)
            self.test_stats.append(
                self.get_test_stats(test_subset_y, test_subset_c, predicted_c)
            )

    def one_shot_test(self):

        test_subset_y, test_subset_c, test_subset_pt, test_subset_ph = zip(
            *random.sample(
                list(
                    zip(
                        self.test_set_y,
                        self.test_set_c,
                        self.test_set_pt,
                        self.test_set_ph,
                    )
                ),
                np.ceil(0.5 * len(self.test_set_y)).astype(np.int),
            )
        )
        test_subset_y = np.array(test_subset_y)
        test_subset_c = np.array(test_subset_c)

        predicted_c = self.decode_model.predict(test_subset_y)

        return self.get_test_stats(test_subset_y, test_subset_c, predicted_c)

    """Does a one-shot through the entire testing set to get a single timecourse"""

    def plot_test_timecourse(self):
        # ok... let's think this through...
        pred_c = self.decode_model.predict(self.test_set_y)
        zipped_data = list(
            zip(pred_c, self.test_set_c.squeeze(), self.test_set_pt, self.test_set_ph)
        )

        pred_add = []
        real_add = []

        predicted_dict = nestdict()

        for pt in self.pts:
            print(pt)
            plt.figure()
            sorted_inp = {do_phase: [0, 0] for do_phase in self.filter_phases}
            sorted_inp = {
                do_phase: np.array(
                    [
                        (predc, testc)
                        for predc, testc, patient, phase in zipped_data
                        if phase == do_phase and patient == pt
                    ]
                ).squeeze()
                for do_phase in self.filter_phases
            }
            predicted = []
            actual = []
            for phase in self.filter_phases:
                if sorted_inp[phase].size == 0:
                    predicted.append(0)
                    actual.append(0)
                else:

                    predicted.append(sorted_inp[phase][0])
                    actual.append(sorted_inp[phase][1])

            plt.plot(predicted)
            plt.plot(actual)
            plt.xlabel("Week")
            plt.ylabel("nHDRS")
            plt.title(pt)

    def plot_test_regression_figure(self):
        # do a final test on *all* the data for plotting purposes
        predicted_c = self.decode_model.predict(self.test_set_y)
        slope = stats.linregress(predicted_c.squeeze(), self.test_set_c.squeeze())
        pearson = stats.pearsonr(predicted_c.squeeze(), self.test_set_c.squeeze())
        r2score = self.decode_model.score(self.test_set_y, self.test_set_c)
        mse = mean_squared_error(self.test_set_c, predicted_c)

        #%%

        plt.figure()
        plt.plot([0, 1], [0, 1], color="gray", linestyle="dotted")
        ax = sns.regplot(x=predicted_c, y=self.test_set_c.squeeze())
        for xx, yy, pp, cc in zip(
            predicted_c, self.test_set_c, self.test_set_pt, self.test_set_ph
        ):
            if self.CFrame.query_stim_change(pp, cc, include_init=False):
                ax.scatter(x=xx, y=yy, marker="^", s=100, color="r")
                ax.text(xx, yy, pp + " " + cc)

        plt.xlabel("Predicted")
        plt.ylabel("Actual")
        plt.title(
            "R2:"
            + str(r2score)
            + "\n"
            + " MSE:"
            + str(mse)
            + " Slope:"
            + str(slope[0])
            + " Pearson:"
            + str(pearson)
        )
        # plt.xlim((-0.1,1.5))
        # plt.ylim((-0.1,1.5))

    def plot_test_stats(self):
        plt.figure()
        plt.subplot(311)
        # plt.scatter(test_subset_c,predicted_c)
        plt.hist([a["Score"] for a in self.test_stats])
        plt.vlines(np.mean([a["Score"] for a in self.test_stats]), 0, 10, linewidth=10)
        plt.title("R2 Score")
        plt.subplot(312)
        # plt.hist([a['Pearson'][0] for a in self.test_stats]);
        plt.vlines(
            np.mean([a["Pearson"][0] for a in self.test_stats]), 0, 10, linewidth=10
        )
        plt.title("Pearson")
        # plt.subplot(313)
        # plt.hist([a['Spearman'][0] for a in self.test_stats]);plt.title('Spearman')
        plt.subplot(313)
        plt.hist([a["Slope"][0] for a in self.test_stats])
        plt.vlines(np.mean([a["Slope"] for a in self.test_stats]), 0, 10, linewidth=10)
        plt.title("Slope")

    """PLOTTING--------------------------------------------------------"""

    """Plot the decoding CV coefficients"""

    def plot_decode_CV(self):
        plt.figure()

        active_coeffs = []
        for ii in self.decode_model_combos_:
            active_coeffs.append([ii.coef_[:]])

        active_coeffs = np.array(active_coeffs).squeeze()
        plt.plot(active_coeffs.T, "r.", markersize=20)

        # plt.plot()
        vp_obj = sns.violinplot(data=active_coeffs, scale="width")
        plt.setp(vp_obj.collections, alpha=0.3)

        average_model, _ = self.get_average_model(self.decode_model_combos_)
        plt.plot(average_model)
        plt.hlines(0, -2, 11, linestyle="dotted")
        plt.ylim((-0.2, 0.2))
        plt.xlim((-1, len(self.do_feats) * 2))


class controller_analysis:
    def __init__(self, readout, **kwargs):
        self.readout_model = readout
        # get our binarized disease states
        self.binarized_type = kwargs["bin_type"]

    def gen_binarized_state(self, **kwargs):
        # redo our testing set
        if kwargs["approach"] == "threshold":
            binarized = kwargs["input_c"] > 0.5
        elif kwargs["approach"] == "stim_changes":
            query_array = kwargs["input_ptph"]
            binarized = [
                self.readout_model.CFrame.query_stim_change(pt, ph)
                for pt, ph in query_array
            ]
        else:
            raise Exception

        return binarized

    def pr_classif(self, binarized, predicted):

        precision, recall, thresholds = precision_recall_curve(binarized, predicted)

        # plt.figure()
        # plt.step(recall,precision)
        return precision, recall

    def pr_oracle(self, binarized, level=0.5):
        oracle = np.array(np.copy(binarized)).astype(np.float)
        oracle += np.random.normal(0, level, size=oracle.shape)

        precision, recall, thresholds = precision_recall_curve(binarized, oracle)
        return precision, recall

    def pr_classif_2pred(self, binarized, predicted, empirical):
        empirical = np.array(empirical).squeeze()
        precision, recall, thresholds = precision_recall_curve(
            binarized, empirical - predicted
        )
        return precision, recall

    def bin_classif(self, binarized, predicted):
        fpr, tpr, thresholds = metrics.roc_curve(binarized, predicted)
        roc_curve = (fpr, tpr, thresholds)
        auc = roc_auc_score(binarized, predicted)

        return auc, roc_curve

    def controller_runs(self):
        controller_types = [
            "readout",
            "empirical+readout",
            "oracle",
            "null",
            "empirical",
        ]
        controllers = {key: [] for key in controller_types}
        aucs = {key: [] for key in controller_types}
        pr_curves = {key: [] for key in controller_types}

        for ii in range(100):
            test_subset_y, test_subset_c, test_subset_pt, test_subset_ph = zip(
                *random.sample(
                    list(
                        zip(
                            self.readout_model.test_set_y,
                            self.readout_model.test_set_c,
                            self.readout_model.test_set_pt,
                            self.readout_model.test_set_ph,
                        )
                    ),
                    np.ceil(0.8 * len(self.readout_model.test_set_y)).astype(np.int),
                )
            )
            predicted_c = self.readout_model.decode_model.predict(test_subset_y)

            # test_subset_pt = shuffle(test_subset_pt);print('PR_Classif: Shuffling Data')
            binarized_c = self.gen_binarized_state(
                approach="stim_changes",
                input_ptph=list(zip(test_subset_pt, test_subset_ph)),
            )
            # shuffle?
            # binarized_c = shuffle(binarized_c);print('PR_Classif: Shuffling binarization')
            coinflip = np.random.choice(
                [0, 1], size=(len(test_subset_pt),), p=[0.5, 0.5]
            )

            controllers["readout"].append(self.pr_classif(binarized_c, predicted_c))
            controllers["empirical+readout"].append(
                self.pr_classif_2pred(binarized_c, predicted_c, test_subset_c)
            )
            controllers["oracle"].append(self.pr_oracle(binarized_c, level=0.5))
            controllers["empirical"].append(self.pr_classif(binarized_c, test_subset_c))
            controllers["null"].append(self.pr_classif(binarized_c, coinflip))

        # organize results
        for kk in controller_types:
            for ii in range(100):
                aucs[kk].append(
                    metrics.auc(controllers[kk][ii][1], controllers[kk][ii][0])
                )
                pr_curves[kk].append((controllers[kk][ii][0], controllers[kk][ii][1]))

            self.plot_classif_runs(aucs[kk], pr_curves[kk], title=kk)

    def classif_runs(
        self,
    ):
        aucs = []
        roc_curves = []

        null_aucs = []
        null_roc_curves = []

        for ii in range(100):
            test_subset_y, test_subset_c, test_subset_pt, test_subset_ph = zip(
                *random.sample(
                    list(
                        zip(
                            self.readout_model.test_set_y,
                            self.readout_model.test_set_c,
                            self.readout_model.test_set_pt,
                            self.readout_model.test_set_ph,
                        )
                    ),
                    np.ceil(0.8 * len(self.readout_model.test_set_y)).astype(np.int),
                )
            )
            # THIS IS WHERE WE NEED TO SHUFFLE TO TEST THE READOU
            # test_subset_y, test_subset_c, test_subset_pt, test_subset_ph = shuffle(test_subset_y, test_subset_c, test_subset_pt, test_subset_ph)
            predicted_c = self.readout_model.decode_model.predict(test_subset_y)

            binarized_c = self.gen_binarized_state(
                approach="threshold", input_c=np.array(test_subset_c)
            )
            auc, roc_curve = self.bin_classif(binarized_c, predicted_c)
            aucs.append(auc)
            roc_curves.append(roc_curve)

            coinflip = np.random.choice(
                [0, 1], size=(len(test_subset_pt),), p=[0.5, 0.5]
            )

            n_auc, n_roc = self.bin_classif(binarized_c, coinflip)
            null_aucs.append(n_auc)
            null_roc_curves.append(n_roc)

        self.plot_classif_runs(aucs, roc_curves)
        # self.plot_classif_runs(null_aucs,null_roc_curves) # if you want a sanity check with a coinflip null

    """Here we'll do a 2-d density plot for error rates using both DR-SCC and nHDRS"""

    def density_plot(self):
        pass

    def plot_classif_runs(self, aucs, roc_curves, **kwargs):
        plt.figure()
        plt.hist(aucs)
        plt.vlines(np.mean(aucs), -1, 10, linewidth=10)
        plt.xlim((0.0, 1.0))
        plt.title(kwargs["title"])

        fig, ax = plt.subplots()
        mean_fpr = np.linspace(0, 1, 100)
        interp_tpr = []
        for aa in roc_curves:
            interp_tpr_individ = interp(mean_fpr, aa[0], aa[1])
            interp_tpr_individ[0] = 0
            interp_tpr.append(interp_tpr_individ)

        mean_tpr = np.mean(interp_tpr, axis=0)
        std_tpr = np.std(interp_tpr, axis=0)

        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)

        ax.plot(mean_fpr, mean_tpr)
        ax.fill_between(mean_fpr, tprs_lower, tprs_upper, alpha=0.2)
        ax.plot(mean_fpr, mean_fpr, linestyle="dotted")
        plt.plot([0, 1], [0, 1], linestyle="dotted")
        if "title" in kwargs:
            plt.title(kwargs["title"])

        # for aa in roc_curves:
        #    plt.plot(aa[0],aa[1],alpha=0.2)

    def roc_auc(self):
        pass


#%%
class feat_check(base_decoder):
    def __init__(self, **kwargs):
        self.YFrame = kwargs["BRFrame"]
        self.fvect = self.YFrame.data_basis["F"]
        self.do_feats = [
            "Delta",
            "Theta",
            "Alpha",
            "Beta*",
            "Gamma1",
            "SHarm",
            "THarm",
            "Stim",
            "Clock",
            "fSlope",
            "nFloor",
            "GCratio",
        ]
        self.pts = ["901", "903", "905", "906", "907", "908"]

        self.filter_recs()

    def calculate_states_in_set(self, data_set):
        state_vector = []

        for rr in data_set:
            psd_poly_done = {
                ch: dbo.poly_subtr_vect(
                    fvect=self.fvect, inp_psd=rr["Data"][ch], polyord=5
                )[0]
                for ch in rr["Data"].keys()
            }

            feat_vect = np.zeros(shape=(len(self.do_feats), self.ch_num))
            for ff, featname in enumerate(self.do_feats):
                dofunc = dbo.FEAT_DICT[featname]
                feat_calc = dofunc["fn"](psd_poly_done, self.fvect, dofunc["param"])
                feat_vect[ff, :] = np.array([feat_calc[ch] for ch in ["Left", "Right"]])

            # We need to flatten the state between channels...
            # Then we go ahead and append it to the state vector
            state_vector.append(
                np.reshape(feat_vect, -1, order="F")
            )  # we want our FEATURE index to change quickest so we go (0,0) -> (1,0) -> (2,0) -> ... (4,1)

            # pt_ph_vector.append(rr['Patient'],rr['Phase'])
        return np.array(state_vector)

    def check_stim_corr(self, band="Beta*", artifact="THarm", pt="ALL"):
        band_idx = self.do_feats.index(band)
        stim_idx = self.do_feats.index(artifact)

        if pt == "ALL":
            do_pts = self.pts
        else:
            do_pts = [pt]

        use_rec_list = [x for x in self.active_rec_list if x["Patient"] in do_pts]
        B_rec_list = [x for x in use_rec_list if x["Phase"][0] == "B"]
        C_rec_list = [x for x in use_rec_list if x["Phase"][0] == "C"]

        B_state = self.calculate_states_in_set(B_rec_list)
        C_state = self.calculate_states_in_set(C_rec_list)

        # get this vector
        colors = ["r", "g"]
        phases = ["B", "C"]
        plt.figure()
        for ss, side in enumerate(["Left", "Right"]):
            plt.subplot(1, 2, ss + 1)
            for seti, state_vect in enumerate([B_state, C_state]):
                band_vect = state_vect[:, band_idx + ss * len(self.do_feats)]
                stim_vect = state_vect[:, stim_idx + ss * len(self.do_feats)]
                corr = stats.pearsonr(band_vect, stim_vect)
                plt.scatter(
                    band_vect, stim_vect, color=colors[seti], label=phases[seti]
                )
                plt.xlabel(band)
                plt.ylabel(artifact)
                plt.title(side + " " + str(corr))
                plt.legend()

        plt.suptitle("Patient: " + pt)


class weakly_decoderCV_Lasso(weekly_decoderCV):
    def __init__(self, *args, **kwargs):
        print("Initialized the Weekly CV decoder")
        super().__init__(*args, **kwargs)

        if kwargs["algo"] == "ENR":

            self.regression_algo = linear_model.Lasso
            # self.model_args = {'alpha':np.e **-3.4,'l1_ratio':0.8} #5/28, good stats, visual is small
            self.model_args = {"alpha": np.e ** kwargs["alpha"]}
            print("Running Lasso")

            # BEFORE 5/28
            # self.regression_algo = linear_model.ElasticNetCV
            # self.model_args = {'alphas':np.linspace(0.01,0.04,20),'l1_ratio':np.linspace(0.1,0.3,10),'cv':10}

        self.pt_CV_sets(n=3)

    def _path_slope_regression(self, do_plot=False, suppress_vars=0.2):
        alpha = 0.1
        return alpha
