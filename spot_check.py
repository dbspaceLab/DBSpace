#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:23:07 2017

@author: virati
"""
import sys
sys.path.append('/home/virati/Dropbox/projects/Research/MDD-DBS/Ephys/IntegratedAnalysis')
#import DBSOsc as dbo
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


#json file for major experiments
experiments = ['Targeting','Amplitude','Frequency','Resting']


flist = []
#flist = ['/home/virati/MDD_Data/BR/908/Session_2016_01_19_Tuesday/908_2016_01_13_16_27_06__MR_2.txt','/home/virati/MDD_Data/BR/908/Session_2016_04_18_Monday/DBS908_2016_04_11_21_15_45__MR_2.txt']
#expname = ['Exp0','Exp1']
#flist = ['/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_04_Wednesday/PCSTES_2018_04_04_15_37_06__MR_1.txt', '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_04_Wednesday/PCSTES_2018_04_04_15_33_36__MR_0.txt']
flist = ['/home/virati/MDD_Data/BR/907/Session_2015_12_17_Thursday/DBS907_2015_12_17_11_39_26__MR_0.txt']
expname=['test']
#Below is the saline testing recordings

#flist = ['/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_16_53_36__MR_0.txt',
# '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_15_20__MR_0.txt',
# '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_32_09__MR_0.txt',
# '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_49_21__MR_0.txt',
# '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_17_56_09__MR_1.txt',
# '/home/virati/MDD_Data/VRT_Impedance_RB/Session_2018_04_24_Tuesday/demo_2018_04_24_18_02_40__MR_2.txt']
#expname = ['IF-300','GE-300','SA-100','2IF-300','2GE-100','2SA-100']

#%%

#import ipdb as pdb

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)
matplotlib.rcParams['svg.fonttype'] = 'none'

plt.rcParams['image.cmap'] = 'jet'

#plt.close('all')

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
        chann = chann_label[cc]
        plt.subplot(2,2,cc+1)
        T = TFcont['TF']['T']
        F = TFcont['TF']['F']
        SG = TFcont['TF']['SG']
        
        t_idxs = np.where(np.logical_and(T > tlim[0], T < tlim[1]))
    
        med_psd = np.median(10*np.log10(SG[chann][:,t_idxs]).squeeze(),axis=1)
        var_psd = np.var(10*np.log10(SG[chann][:,t_idxs]).squeeze(),axis=1).reshape(-1,1)
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
        
    fcalced,bands = dbo.calc_feats(pf_lPSD,F,modality='lfp')
    
    plt.figure(osc_feat.number)
    plt.subplot(1,2,1)
    plt.plot(fcalced[:,0],label=title)
    plt.title('Left')
    plt.subplot(1,2,2)
    plt.plot(fcalced[:,1],label=title)
    plt.title('Right')
    plt.suptitle('Features')

#%%


def spot_check(fname,tlims=(0,-1),plot_sg=False,chann_labels=['Left','Right']):
    ''' Spotcheck function
    tlims - In seconds. -1 implies end of the recording
    '''
    
    Container = DBSOsc.load_BR_feats(fname,snippet=False)
    
    NFFT = 2**10
    fs = 422 #hardcoded for brLFP for now
    
    inv_try = [None] * 100
    # Try inverse tanh
    for cc,cs in enumerate(np.linspace(0.002,1,100)):
        inv_try[cc] = np.arctanh(Container['TS']['Y']/cs)
    
    #go to each and compute how much power there is in 32/64 Hz band and find the smallest
    #for cc in range(len(inv_try)):
     
    inv_x = inv_try[0]
    
    #What are our time limits?
    nlims = np.array(tlims) * fs
    
    if tlims[1] == -1:
        nlims[1] == -1
    
    
    ## Do spectrogram stuff
    SG = defaultdict(dict)
    Pxx = defaultdict(dict)
    for cc in range(2):
        #first, let's do the PWelch
        Fpsd,Pxx[chann_labels[cc]] = sig.welch(Container['TS']['Y'][nlims[0]:nlims[1],cc],fs,window='blackmanharris',nperseg=NFFT,noverlap=0,nfft=NFFT)
        
        
        F,T,SG[chann_labels[cc]] = sig.spectrogram(Container['TS']['Y'][nlims[0]:nlims[1],cc],nperseg=NFFT,noverlap=NFFT*0.5,window=sig.get_window('blackmanharris',NFFT),fs=422)    
        #Need to transpose for the segmentation approach to work, retranspose in the plotting
    
    polycorr = True
    if plot_sg:
        plt.figure()
        #if we want to do polynomial corrections
        if polycorr:
            print('Polynom Correction!')
            corr_sg = dbo.poly_SG(SG,F)
            
            pdb.set_trace()
            
        #Now we want to plot
        else:
            for cc in range(2):
                    
                    
                    
                plt.pcolormesh(T,F,10*np.log10(SG[chann_labels[cc]].T),rasterized=True)
                plt.clim((-200,-100))
            plt.suptitle('Raw TS: ' + fname.split('/')[-1])
        
        
      
    
    

    #if you want to return the raw
    #return {'TS':Container['TS']['Y'],'TF':{'SG':SG,'F':F,'T':T},'F':{'Pxx':Pxx,'F':Fpsd},'InvTS':inv_x}
    #if you want to try to work with the inverse
    print('RETURNING INVERSE!!')
    return {'TS':0.004*inv_x/(3),'TF':{'SG':SG,'F':F,'T':T},'F':{'Pxx':Pxx,'F':Fpsd},'InvTS':inv_x}



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
