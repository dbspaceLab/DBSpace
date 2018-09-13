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

flist = '/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_16_53_36__MR_0.txt'
expname = 'IF-300'
V_list = [2,4,6,8]
V_dict = {2:(0,0),4:(0,0),6:(0,0),8:(230,250)} #The times *in seconds) that the voltages are active

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
            
    for vv in V_list:
        results[expname[ff]] = spot_check(fname,plot_sg=True)

    root.destroy()
    
    
    #%%
    bigmed = plt.figure()
    osc_feat = plt.figure()
    
    for key in expname:
        print(key)
        for vv in V_list:
            grab_median(results[key],tlim=do_voltage,title=key,do_corr=False)
    
    plt.legend()
