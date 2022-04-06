import numpy as np
from sklearn import metrics
from sklearn.metrics import (
    auc,
    average_precision_score,
    mean_absolute_error,
    mean_squared_error,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.utils import shuffle
from dbspace.utils.structures import nestdict
import matplotlib as plt
from scipy import interp
import random


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

    def controller_simulations(self):
        """
        Controller Types:
            "Readout": The main DR-SCC
            "Empirical + Readout": The nHDRS with the readout
            "Empirical + inv_readout": The nHDRS with the inverse of the readout
            "Oracle": The best case scenario along with some noise
            "Null": Pure chance
            "Empirical": The nHDRS

        """
        controllers = nestdict()
        controllers = nestdict()
        
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
            controllers["inv_readout"].append(self.pr_classif(binarized_c, 1/predicted_c))
            controllers["empirical+readout"].append(
                self.pr_classif_2pred(binarized_c, predicted_c, test_subset_c)
            )
            controllers["empirical+inv_readout"].append(
                self.pr_classif_2pred(binarized_c, 1/predicted_c, test_subset_c)
            )
            controllers["oracle"].append(self.pr_oracle(binarized_c, level=0.5))
            controllers["empirical"].append(self.pr_classif(binarized_c, test_subset_c))
            controllers["null"].append(self.pr_classif(binarized_c, coinflip))

        self.controllers = controllers

    def controller_sim_metrics(self):
        # organize results
        controllers = self.controllers
        aucs = nestdict()
        pr_curves = nestdict()

        plot_controllers = ["readout","empirical","null","oracle"]
        for kk in plot_controllers:
            for ii in range(100):
                aucs[kk].append(
                    metrics.auc(controllers[kk][ii][1], controllers[kk][ii][0])
                )
                pr_curves[kk].append((controllers[kk][ii][0], controllers[kk][ii][1]))

            self.plot_controller_runs(aucs[kk], pr_curves[kk], title=kk)

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

        self.plot_controller_runs(aucs, roc_curves)

    def plot_controller_runs(self, aucs, roc_curves, **kwargs):
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

    def plot_controller_simulations(self, plot_controller_list):
        """
        Generate a compound plot of all simulations desired
        """
        fig, ax = plt.subplots()
        if set(plot_controller_list) != set(self.controllers.keys()):
            raise ValueError("There's a mismatch in the controllers you want...")
        for controller in plot_controller_list:
            pass
