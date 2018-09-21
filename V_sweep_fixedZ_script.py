#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 11:35:00 2018

@author: virati
This script does a fixed voltage sweep check in the saline experiments.

"""
import sys
sys.path.append('/home/virati/Dropbox/projects/Research/MDD-DBS/Ephys/IntegratedAnalysis')
import DEPR_DBSOsc as DBSOsc
import DBS_Osc as dbo

import matplotlib
import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
import scipy.io as io
from collections import defaultdict
from DBS_Osc import nestdict

import pdb

from tkinter.filedialog import askopenfilename
import tkinter as tk

import matplotlib.pyplot as plt

from spot_check import *


chann_label = ['Left','Right']

fname = {'IF-300':'/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_16_53_36__MR_0.txt', 
         'SA-100':'/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_32_09__MR_0.txt',
         '2SA-100': '/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_18_02_40__MR_2.txt'}

do_exp = 'IF-300'
V_list = [2,4,6,8]
V_dict = {2:(40,50),4:(110,120),6:(170,180),8:(230,250)} #The times *in seconds) that the voltages are active

#%%

#import ipdb as pdb

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)
matplotlib.rcParams['svg.fonttype'] = 'none'

plt.rcParams['image.cmap'] = 'jet'

#Parameters for analysis
band_scheme = 'Adjusted'
band_compute = 'median'

if __name__ == '__main__':

    root = tk.Tk()
    root.withdraw()
    
    plt.ion()
    
    results = defaultdict(dict)

    #Generate a single SG for the entire experiment
    _ = spot_check(fname[do_exp],plot_sg=True)            

    results = spot_check(fname[do_exp],plot_sg=False)

    root.destroy()
    
    
    #%%
    bigmed = plt.figure()
    osc_feat = plt.figure()

    for vv in V_list:
        grab_median(results,bigmed,osc_feat,tlim=V_dict[vv],title='Voltage ' + str(vv),do_corr=False)
    
    plt.legend()
