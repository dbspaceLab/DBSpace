#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 11:34:21 2018

@author: virati
Script for VOLTAGE SPOT CHECK Only
This script checks across experiments/impedance mismatch conditions

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

#json file for major experiments
experiments = ['Targeting','Amplitude','Frequency','Resting']


flist = []
#flist = ['/home/virati/MDD_Data/BR/908/Session_2016_01_19_Tuesday/908_2016_01_13_16_27_06__MR_2.txt','/home/virati/MDD_Data/BR/908/Session_2016_04_18_Monday/DBS908_2016_04_11_21_15_45__MR_2.txt']
#expname = ['Exp0','Exp1']
#flist = ['/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_04_Wednesday/PCSTES_2018_04_04_15_37_06__MR_1.txt', '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_04_Wednesday/PCSTES_2018_04_04_15_33_36__MR_0.txt']
flist = ['/home/virati/MDD_Data/BR/907/Session_2015_12_17_Thursday/DBS907_2015_12_17_11_39_26__MR_0.txt']
expname=['test']

chann_label = ['Left','Right']


#Below is the saline testing recordings
flist = ['/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_16_53_36__MR_0.txt',
 '/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_15_20__MR_0.txt',
 '/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_32_09__MR_0.txt',
 '/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_49_21__MR_0.txt',
 '/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_56_09__MR_1.txt',
 '/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_18_02_40__MR_2.txt']
expname = ['IF-300','GE-300','SA-100','2IF-300','2GE-100','2SA-100']
#(Interface, pure gel, pure saline; second interface, second gel, second saline)
do_voltage = (230,250)

#flist = ['/home/virati/MDD_Data/Benchtop/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_16_53_36__MR_0.txt']
#expname = ['IF-300']

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

#%%



if __name__ == '__main__':
    patients = ['901']#,'903','905','906','907','908']
    #file chooser
    patient=  '905'
    
    
    root = tk.Tk()
    root.withdraw()
    
    notdone = 1
    plt.ion()
    
    #curr_dir = '/home/virati/MDD_Data/BR/'
    curr_dir = '/run/media/virati/'
    
    results = defaultdict(dict)
    
    if flist == []:
        while notdone:
            fname = askopenfilename(initialdir=curr_dir)
            
            if fname == None or fname == '':
                notdone = 0
                print('No File Selected: Goodbye!')
            else:
                print('File selected: ' + fname[-20:])
                flist.append(fname)
                #results[fname] = spot_check(fname)
                curr_dir = '/'.join(fname.split('/')[:-1])
                print(curr_dir)
            
            
    for ff,fname in enumerate(flist):
        results[expname[ff]] = spot_check(fname,plot_sg=True)
    
    
    root.destroy()
    
    
    #%%
    bigmed = plt.figure()
    osc_feat = plt.figure()
    
    for key in expname:
        print(key)
        #6v
        #grab_median(val,tlim=(223,254),title=key,do_corr=False)
        #2v
        grab_median(results[key],bigmed,osc_feat,tlim=do_voltage,title=key,do_corr=False)
    
    plt.legend()    
    
    bigmed = plt.figure()
    osc_feat = plt.figure()
    for key,val in results.items():
        print(key)
        #6v
        #grab_median(val,tlim=(223,254),title=key,do_corr=False)
        #2v
        grab_median(val,bigmed,osc_feat,tlim=do_voltage,title=key,do_corr=True)
    
    plt.legend()    
    #%%
    
    plt.show()

#grab_median(results['/home/virati/MDD_Data/BR/907/Session_2015_12_17_Thursday/DBS907_2015_12_17_11_39_26__MR_0.txt'],tlim=(70,120))
#grab_median(results['/home/virati/MDD_Data/BR/907/Session_2015_12_17_Thursday/DBS907_2015_12_17_11_39_26__MR_0.txt'],tlim=(880,930))


#plt.ion()
#data = defaultdict(dict)
#for pt in patients:
#    plt.figure()
#    data[pt] = spot_check('/home/virati/temp_pavel_6mo/' + pt + '_LFP_OnTarget')
#    plt.show()
#    _ = input('press enter')
#    
##io.savemat('/tmp/Recording_Bryan.mat',data)
# 
