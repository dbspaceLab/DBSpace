#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 20:06:36 2018

@author: virati
Tractography Class to preprocess and package DTI data relevant to project
"""


import json
from tkinter import W

import dbspace as dbo
import matplotlib.pyplot as plt
import nibabel
import nilearn
import nilearn.image as image
import numpy as np
from dbspace.utils.structures import nestdict
from nilearn import image, plotting


class engaged_tractography:
    def __init__(
        self,
        do_pts=["901", "903", "905", "906", "907", "908"],
        v_list=range(2,8),:W
        do_condits=["OnT", "OffT"],
        electrode_map=None,
    ):
        self.do_pts = do_pts
        self.v_list = list(v_list)
        self.do_condits = do_condits
        self.stim_configurations = ["L", "R"]

        self.load_electrode_map(electrode_map)

    def load_dti(self):
        dti_file = nestdict()
        data_arr = nestdict()

        for pp, pt in enumerate(self.do_pts):
            for cc, condit in enumerate(self.do_condits):
                for vv, vstim in enumerate(self.v_list):
                    dti_file = {key:[] for key in self.stim_configurations}
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
                        dti_file[pt][condit][vstim][side] = fname
                        dti_file[side] = fname
                    
                    dti_data = {side: image.load_img(dti_file[side]) for side in self.stim_configurations}

                    bilateral_dti_data_at_v = image.math_img("img1+img2",img1=dti_data[self.stim_configurations[0]], img2=dti_data[self.stim_configurations[1]])
                    if vstim != 2:
                        data_arr[pt][condit] = image.math_img("img1+img2", img1=data_arr[pt][condit], img2=bilateral_dti_data_at_v)
                    else:
                        data_arr[pt][condit] = copy(bilateral_dti_data_at_v)
        
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

    def plot_engaged_DTI(self, pt, condit="OnT"):
        """
        Engaged DTI is defined to be the average throughout 2-7V
        Previous called 'DTI Flow'

        """
        combined = self.combined
        stacked = image.math_img()

    def plot_V_DTI(self, pt="906", condit="OnT", v_select=2, merged=False):
        combined = self.combined
        vidx = self.v_list.index(v_select)

        condit = "OnT"

        if merged:
            stacked = image.math_img(
                "img1+img2+img3+img4+img5+img6",
                img1=combined[pt][condit][2],
                img2=combined[pt][condit][3],
                img3=combined[pt][condit][4],
                img4=combined[pt][condit][5],
                img5=combined[pt][condit][6],
                img6=combined[pt][condit][7],
            )

            plotting.plot_glass_brain(
                stacked,
                black_bg=True,
                title=condit + " Tractography",
                vmin=-15,
                vmax=15,
            )
        else:
            new_img = nilearn.image.new_img_like(
                self.data[pt][condit][v_select]["L"], (self.middle_idx)
            )
            plotting.plot_glass_brain(new_img)


if __name__ == "__main__":
    do_pts = ["907"]
    V_DTI = DTI(do_pts=do_pts)
    V_DTI.load_data()
    V_DTI.plot_V_DTI(pt=do_pts[0], merged=True)
