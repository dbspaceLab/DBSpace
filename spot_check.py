#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:23:07 2017

@author: virati
Spot check GUI. THIS STILL USES THE OLD DBSOsc and needs to be ported to new DBS_Osc library. But everything breaks for 50 reasons, so might be best to just start from scratch
TODO This needs to be stripped to bare-minimum spot checking of an arbitrary BRadio File
"""


import DBSpace as dbo

import matplotlib
import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
import scipy.io as io
from collections import defaultdict
from DBSpace import nestdict

import pdb

from tkinter.filedialog import askopenfilename
import tkinter as tk

import matplotlib.pyplot as plt

# Flist will be all the files we want to spotcheck, with the key as the experiment/info and the fname as the file loaded in
flist = {'DBS905_VSweep':{'fname':'/home/extend/MDD_Data/BR/905/Session_2015_09_29_Tuesday/Dbs905_2015_09_29_13_19_27__MR_0.txt'}}


#%%

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)
matplotlib.rcParams['svg.fonttype'] = 'none'

plt.rcParams['image.cmap'] = 'jet'

#Parameters for analysis
#band_scheme = 'Adjusted'
#band_compute = 'median'



#%%
# General methods

def gui_file_select():
    curr_dir = '/run/media/virati/'
    notdone = True
    flist = []
    
    while notdone:
        fname = askopenfilename(initialdir=curr_dir)
        if fname == None or fname == '':
            notdone = False
        else:
            flist.append(fname)
            curr_dir = '/'.join(fname.split('/')[:-1])

    return flist

def grab_median(TFcont,bigmed,osc_feat,tlim=(880,900),title='',do_corr=True,band_compute='median',band_scheme='Adjusted'):
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
            fixed_psd, polyitself = dbo.poly_subtr(corr_psd,F)
            pf_lPSD[chann_label[cc]] = fixed_psd[chann_label[cc]].reshape(-1,1)
        else:
            correct_psd, polyitself = dbo.poly_subtr(corr_psd,F)

            pf_lPSD[chann_label[cc]] = 10**(med_psd/10).reshape(-1,1)
            plt.plot(F,polyitself,label='Polynomial Fit',color='black')
            
        plt.plot(F,10*np.log10(pf_lPSD[chann_label[cc]]),label=title)
        plt.title('Channel ' + chann_label[cc] + ' psd')
        #try: plt.fill_between(F,(10*np.log10(pf_lPSD[chann_label[cc]]))+var_psd,(10*np.log10(pf_lPSD[chann_label[cc]]))-var_psd)
        #except e: print(e);pdb.set_trace()
        
        plt.ylim(psd_lim)
        
        plt.subplot(2,2,2 + (cc+1))
        plt.plot(F,10*np.log10(var_psd),label=title)
        plt.title('Variance in PSD across time: ' + chann_label[cc])
    plt.subplot(2,2,4)
    plt.legend()
    

    if band_scheme == 'Standard':
        band_wins = ['Delta','Theta','Alpha','Beta','Gamma']
    elif band_scheme == 'Adjusted':
        band_wins = ['Delta','Theta','Alpha','Beta*','Gamma1']
    
    
    fcalced,bands = dbo.calc_feats(pf_lPSD,F,dofeats=band_wins, modality='lfp',compute_method=band_compute)
    
    plt.figure(osc_feat.number)
    plt.subplot(1,2,1)
    plt.plot(fcalced[:,0],label=title)
    
    plt.title('Left')
    plt.subplot(1,2,2)
    plt.plot(fcalced[:,1],label=title)
    plt.title('Right')
    plt.suptitle('Features ' + band_compute + ' ' + band_scheme)

#%%

def spot_SG(fname,chann_labels=['Left','Right']):
    F,T,SG[chann_labels[cc]] = sig.spectrogram(Container['TS']['Y'][nlims[0]:nlims[1],cc],nperseg=NFFT,noverlap=NFFT*0.5,window=sig.get_window('blackmanharris',NFFT),fs=422)
    

def spot_check(fname=[],tlims=(0,-1),plot_sg=False,chann_labels=['Left','Right']):
    ''' Spotcheck function
    tlims - In seconds. -1 implies end of the recording
    '''
    NFFT = 2**10
    fs = 422 #hardcoded for brLFP for now
    
    if 'flist' in globals():
        curr_exp = flist.keys()
    else:
        curr_exp = 'Generic'
        
    for key in flist.keys():
        #Container = DBSOsc.load_BR_feats(fname,snippet=False)
        Container = dbo.load_BR_dict(flist[key]['fname'],sec_offset=0)
    
    
    nlims = np.array(tlims) * fs
    
    if tlims[1] == -1:
        nlims[1] == -1
    
    
    ## Do spectrogram stuff
    SG = defaultdict(dict)
    Pxx = defaultdict(dict)
    for cc,side in enumerate(['Left','Right']):
        #first, let's do the PWelch
        Fpsd,Pxx[chann_labels[cc]] = sig.welch(Container[side][nlims[0]:nlims[1]],fs,window='blackmanharris',nperseg=NFFT,noverlap=0,nfft=NFFT)
        
        
        F,T,SG[side] = sig.spectrogram(Container[side][nlims[0]:nlims[1]],nperseg=NFFT,noverlap=0,window=sig.get_window('blackmanharris',NFFT),fs=422)    
        #Need to transpose for the segmentation approach to work, retranspose in the plotting
    #%%
    polycorr = False
    if plot_sg:
        plt.figure()
        #if we want to do polynomial corrections
        if polycorr:
            print('Polynom Correction!')
            for chann in SG.keys():
                corr_sg = dbo.poly_SG(SG[chann],F)
            
        #Now we want to plot
        else:
            for cc,side in enumerate(['Left']):
                plt.subplot(2,1,1)
                plt.plot(Container[side][:])
                plt.subplot(2,1,2)
                plt.pcolormesh(T,F,10*np.log10(SG[side]),rasterized=True)
                #plt.clim((-200,-100))
                plt.title('Channel ' + chann_label[cc])
                
            #plt.suptitle('Raw TS: ' + fname.split('/')[-1])
            
            plt.suptitle('Raw TS: ' + curr_exp)
    
    #TODO strip out the inverses from spot_check. This is focused on just displaying as a tool. Can make a separate script that wraps spot_check with preprocessing steps
    print('RETURNING INVERSE!!')
    return {'TS':Container,'TF':{'SG':SG,'F':F,'T':T},'F':{'Pxx':Pxx,'F':Fpsd}}




#%%
#Unit test for the methods above. TODO make this all more OOP
    
if __name__ == '__main__':
    patients = ['901']#,'903','905','906','907','908']
    #file chooser
    patient=  '905'
    
    
    root = tk.Tk()
    root.withdraw()
    
    plt.ion()
    
    results = defaultdict(dict)
    
    if flist == []:
        flist = gui_file_select()
            
    
    #Here, we do the actual spot_checking method
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
