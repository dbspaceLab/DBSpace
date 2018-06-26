#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:23:07 2017

@author: virati
"""
import sys
#sys.path.append('/home/virati/Dropbox/projects/Research/MDD-DBS/Ephys/IntegratedAnalysis')
import DBSOsc
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


#json file for major experiments
experiments = ['Targeting','Amplitude','Frequency','Resting']


flist = []
#flist = ['/home/virati/MDD_Data/BR/908/Session_2016_01_19_Tuesday/908_2016_01_13_16_27_06__MR_2.txt','/home/virati/MDD_Data/BR/908/Session_2016_04_18_Monday/DBS908_2016_04_11_21_15_45__MR_2.txt']
#expname = ['Exp0','Exp1']
#flist = ['/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_04_Wednesday/PCSTES_2018_04_04_15_37_06__MR_1.txt', '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_04_Wednesday/PCSTES_2018_04_04_15_33_36__MR_0.txt']
#flist = ['/home/virati/MDD_Data/BR/907/Session_2015_12_17_Thursday/DBS907_2015_12_17_11_39_26__MR_0.txt']

#Below is the saline testing recordings

flist = ['/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_16_53_36__MR_0.txt',
 '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_15_20__MR_0.txt',
 '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_32_09__MR_0.txt',
 '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_49_21__MR_0.txt',
 '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_56_09__MR_1.txt',
 '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_18_02_40__MR_2.txt']
expname = ['IF-300','GE-300','SA-100','2IF-300','2GE-100','2SA-100']

#%%

#import ipdb as pdb

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)
matplotlib.rcParams['svg.fonttype'] = 'none'

plt.rcParams['image.cmap'] = 'jet'

#plt.close('all')


def spot_check(fname):
    
    Container = DBSOsc.load_BR_feats(fname,snippet=False)
    #plt.plot(Container['TS']['T'][:-1],Container['TS']['Y'][:,0])
    
    #plt.figure()
    #plt.subplot(2,1,1)
    #plt.plot(Container['TS']['T'],Container['TS']['Y'][:,0])
    #plt.ylim((-0.001,0.005))
    
    NFFT = 2**10
    SG = defaultdict(dict)
    for cc in range(2):
        
        F,T,SG[cc] = sig.spectrogram(Container['TS']['Y'][:,cc],nperseg=NFFT,noverlap=NFFT*0.5,window=sig.get_window('blackmanharris',NFFT),fs=422)
        #plt.subplot(2,2,3+cc)
        #plt.subplot(212)
        
        plt.pcolormesh(T,F,10*np.log10(SG[cc]),rasterized=True)
        #plt.clim((-200,-100))
    
    #plt.colorbar()
    #plt.suptitle('Raw TS: ' + fname.split('/')[-1])
    
    

    return {'TS':Container['TS']['Y'],'TF':{'SG':SG,'F':F,'T':T}}

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
    results[expname[ff]] = spot_check(fname)


root.destroy()


#%%
def grab_median(TFcont,tlim=(880,900),title='',do_corr=True):
    #Plot some PSDs
    #plt.figure()
    chann_label = ['Left','Right']
    pf_lPSD = nestdict()
    
    if do_corr:
        psd_lim = (-20,50)
    else:
        psd_lim = (-220,-70)
    
    #Make the big figure that will have both channels
    plt.figure(bigmed.number)
    for cc in range(2):
        plt.subplot(2,2,cc+1)
        T = TFcont['TF']['T']
        F = TFcont['TF']['F']
        SG = TFcont['TF']['SG']
        
        t_idxs = np.where(np.logical_and(T > tlim[0], T < tlim[1]))
    
        med_psd = np.median(10*np.log10(SG[cc][:,t_idxs]).squeeze(),axis=1)
        var_psd = np.var(10*np.log10(SG[cc][:,t_idxs]).squeeze(),axis=1).reshape(-1,1)
        corr_psd = {chann_label[cc]:10**(med_psd/10)}
        
        if do_corr:
            #do polynomial subtraction
            pf_lPSD[chann_label[cc]] = dbo.poly_subtr(corr_psd,F)[chann_label[cc]].reshape(-1,1)
            #plt.ylim((0,10))
        else:
            pf_lPSD[chann_label[cc]] = 10**(med_psd/10).reshape(-1,1)
            
        plt.plot(F,10*np.log10(pf_lPSD[chann_label[cc]]),label=title)
        plt.title('Channel ' + chann_label[cc] + ' psd')
        #try: plt.fill_between(F,(10*np.log10(pf_lPSD[chann_label[cc]]))+var_psd,(10*np.log10(pf_lPSD[chann_label[cc]]))-var_psd)
        #except e: print(e);pdb.set_trace()
        
        plt.ylim(psd_lim)
        plt.legend()
        
        plt.subplot(2,2,2 + (cc+1))
        plt.plot(F,10*np.log10(var_psd),label=title)
        plt.title('Variance in PSD across time: ' + chann_label[cc])
        plt.legend()
        
    fcalced = dbo.calc_feats(pf_lPSD,F,modality='lfp')
    
    plt.figure(osc_feat.number)
    plt.subplot(1,2,1)
    plt.plot(fcalced[:,0],label=title)
    plt.title('Left')
    plt.subplot(1,2,2)
    plt.plot(fcalced[:,1],label=title)
    plt.title('Right')
    plt.suptitle('Features')

#%%
bigmed = plt.figure()
osc_feat = plt.figure()
for key in expname:
    print(key)
    #6v
    #grab_median(val,tlim=(223,254),title=key,do_corr=False)
    #2v
    grab_median(results[key],tlim=(40,50),title=key,do_corr=False)

plt.legend()    

bigmed = plt.figure()
osc_feat = plt.figure()
for key,val in results.items():
    print(key)
    #6v
    #grab_median(val,tlim=(223,254),title=key,do_corr=False)
    #2v
    grab_median(val,tlim=(40,50),title=key,do_corr=True)

plt.legend()    
#%%


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