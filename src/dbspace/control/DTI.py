#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 20:06:36 2018

@author: virati
Tractography Class to preprocess and package DTI data relevant to project
"""


import json
from tkinter import W
import tqdm

import dbspace as dbo
import matplotlib.pyplot as plt
import nibabel
import nilearn
import nilearn.image as image
import numpy as np
from dbspace.utils.structures import nestdict
from dbspace.utils.functions import unity
from nilearn import image, plotting
import copy
import itertools


class engaged_tractography:
    def __init__(
        self,
        target_electrode_map,
        do_pts=["901", "903", "905", "906", "907", "908"],
        v_list=range(2, 8),
        do_condits=["OnT", "OffT"],
    ):
        self.do_pts = do_pts
        self.v_list = list(v_list)
        self.do_condits = do_condits
        self.stim_configurations = ["L", "R"]

        self.load_electrode_map(target_electrode_map)

    def load_dti(self, hide_progress=False):
        data_arr = nestdict()

        if hide_progress:
            progress_fn = unity
        else:
            progress_fn = tqdm.tqdm
        for pp, pt in enumerate(self.do_pts):
            print(pt)
            for cc, condit in enumerate(self.do_condits):
                for vv, vstim in progress_fn(enumerate(self.v_list)):
                    dti_file = {key: [] for key in self.stim_configurations}
                    for ss, side in enumerate(self.stim_configurations):
                        cntct = dbo.Etrode_map[condit][pt][ss] + 1
                        fname = (
                            "/home/virati/Dropbox/projects/Research/MDD-DBS/Data/Anatomy/DTI/MDT_DBS_2_7V_Tractography/DBS"
                            + str(pt)
                            + "."
                            + side
                            + str(cntct)
                            + "."
                            + str(vstim)
                            + "V.bin.nii.gz"
                        )
                        dti_file[side] = fname

                    dti_data = {
                        side: image.load_img(dti_file[side])
                        for side in self.stim_configurations
                    }

                    bilateral_dti_data_at_v = image.math_img(
                        "img1+img2",
                        img1=dti_data[self.stim_configurations[0]],
                        img2=dti_data[self.stim_configurations[1]],
                    )
                    if vstim != 2:
                        data_arr[pt][condit] = image.math_img(
                            "img1+img2",
                            img1=data_arr[pt][condit],
                            img2=bilateral_dti_data_at_v,
                        )
                    else:
                        data_arr[pt][condit] = copy.copy(bilateral_dti_data_at_v)

        self.dti_data = data_arr

    def load_electrode_map(self, target_map_config):
        with open(target_map_config, "r") as electrode_map:
            self.electrode_map = json.load(electrode_map)

    def plot_V_thresh(self, pt="906", condit="OnT"):
        vstim = 2

        new_img = nilearn.image.new_img_like(
            self.data[pt][condit][vstim]["L"], (self.middle_idx)
        )
        plotting.plot_glass_brain(new_img)

    """
    This method plots the DTI for a given patient x condition combination
    """

    def plot_engaged_tractography(self, condits=["OnT", "OffT"]):

        for cc, condit in enumerate(condits):

            engaged_tracto = self.get_engaged_tractography(condit=condit)

            plotting.plot_glass_brain(
                engaged_tracto,
                black_bg=True,
                title=condit + " Tractography",
                vmin=-15,
                vmax=15,
            )

    def get_engaged_tractography(self, condit):
        dti_data = self.dti_data

        avg_image = [dti_data[pt][condit] for pt in self.do_pts]

        keys = [f"img{n}" for n, pt in enumerate(self.do_pts)]
        # sum_string = "+".join(keys)
        mean_string = "np.mean(np.array([" + ",".join(keys) + "]),axis=0)"
        sum_args = {f"img{n}": dti_data[pt][condit] for n, pt in enumerate(self.do_pts)}

        return image.math_img(mean_string, **sum_args)

    def plot_preference_mask(self, condits=["OnT", "OffT"], threshold=0.05):
        if len(condits) != 2:
            raise ValueError("Preference Mask needs two conditions to compare...")

        diff_map = nestdict()
        dti_data = {key: self.get_engaged_tractography(condit=key) for key in condits}
        for ccs in itertools.product(condits, repeat=2):
            if ccs[0] == ccs[1]:
                continue
            diff_map[ccs[0]] = image.math_img(
                "img1 > img2+" + str(threshold),
                img1=dti_data[ccs[0]],
                img2=dti_data[ccs[1]],
            )

        for target in condits:
            plotting.plot_glass_brain(
                diff_map[target],
                black_bg=True,
                title=target + " Engaged Preference Mask",
                vmin=-1,
                vmax=1,
            )

    plt.show()


if __name__ == "__main__":
    do_pts = ["907"]
    V_DTI = DTI(do_pts=do_pts)
    V_DTI.load_data()
    V_DTI.plot_V_DTI(pt=do_pts[0], merged=True)
