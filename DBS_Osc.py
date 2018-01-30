#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 13:44:38 2018

@author: virati
DBS_Osc File for the analytical methodologies used in ALL Vineet project related analyses
"""

import numpy as np
import pandas as pd
from collections import defaultdict
import scipy.signal as sig

#Definitions for oscillatory features -> Oscillatory Vector
band_dict = {'Delta':(1,4), 'Theta':(4,8), 'Alpha':(8,14), 'Beta*':(14,20), 'Beta+':(25,30), 'Gamma1':(30,50), 'Gamma2':(50,80), 'Gamma3':(80,100), 'Stim':(128,132)}
band_order = ['Delta','Theta','Alpha','Beta*','Beta+','Gamma1','Gamma2','Gamma3','Stim']


#Method to load in brainradio file
def load_br_file(fname):
    rawdata = np.array(pd.read_csv(fname,sep=',',header=None))
    return rawdata
    
def gen_psd(inpX):
    #inp X is going to be assumed to be a dictionary with different keys for different channels
    outPSD = defaultdict(dict)
    for chann in inpX.keys():
        outPSD[chann] = F_Domain(inpX[chann])['PSD']   
    
    return outPSD

def gen_TF(ts):
    pass

#extract the oscillatory state    
def extr_Ostate(PSD,norm=None):
    #norm is the normalization method
    pass


def F_Domain(timeser):
    #pdb.set_trace()
    #assert isinstance(timeser,dbs.timeseries)
    #Window size is about 1 second (512 samples is over 1 sec)
    
    Fvect,Pxx = sig.welch(timeser,422,window='blackmanharris',nperseg=512,noverlap=128,nfft=2**10)
    FreqReturn = {'F': Fvect,'PSD': Pxx}
    
    return FreqReturn

def TF_Domain(timeser):
    #assert isinstance(timeser,dbs.timeseries)
    
    F,T,SG = sig.spectrogram(timeser,nperseg=2**10,noverlap=2**10-50,window=sig.get_window('blackmanharris',2**10),fs=422)
    
    TFreqReturn = {'T': T,'F':F,'SG': SG}
    
    return TFreqReturn

def plot_TF(TFR,ch):
    plt.figure()
    
    plt.pcolormesh(TFR['T'],TFR['F'],10*(TFR['SG'][ch,:,:]))
    plt.xlabel('Time')
    plt.ylabel('Frequency')