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

import pdb

np.seterr(divide='raise')

#Definitions for oscillatory features -> Oscillatory Vector
band_dict = {'Delta':(1,4), 'Theta':(4,8), 'Alpha':(8,14), 'Beta*':(14,20), 'Beta+':(25,30), 'Gamma1':(30,50), 'Gamma2':(50,80), 'Gamma3':(80,100), 'Stim':(128,132)}
band_order = ['Delta','Theta','Alpha','Beta*','Beta+','Gamma1','Gamma2','Gamma3','Stim']

all_pts = ['901','903','905','906','907','908']

#Method to load in brainradio file
def load_br_file(fname):
    rawdata = np.array(pd.read_csv(fname,sep=',',header=None))
    return rawdata
    
def gen_psd(inpX,Fs=422,nfft=2**10):
    #inp X is going to be assumed to be a dictionary with different keys for different channels
    outPSD = defaultdict(dict)
    for chann in inpX.keys():
        #The return here is a dictionary with two keys: F and PSD
        outPSD[chann] = F_Domain(inpX[chann],Fs=Fs,nfft=nfft)['Pxx']

    #Return here is a dictionary with Nchann keys
    return outPSD

def get_pow(Pxx,F,frange,cmode=np.median):
    #Pxx is a dictionary where the keys are the channels, the values are the [Pxx desired]
    
    #check if Pxx is NOT a dict
    if isinstance(Pxx,np.ndarray):
        Pxx = {0:Pxx}
    
    #find the power in the range of the PSD
    #Always assume PSD is a dictionary of channels, and each value is a dictionary with Pxx and F
    
    #frange should just be a tuple with the low and high bounds of the band
    out_feats = {keys:0 for keys in Pxx.keys()}
    
    Fidxs = np.where(np.logical_and(F > frange[0],F < frange[1]))
    
    for chans,psd in Pxx.items():
        #if we want the sum
        #out_feats[chans] = np.sum(psd[Fidxs])
        #if we want the MEDIAN instead
        out_feats[chans] = cmode(psd[Fidxs])
    
    #return is going to be a dictionary with same elements
    return out_feats

def F_Domain(timeser,nperseg=512,noverlap=128,nfft=2**10,Fs=422):
    #pdb.set_trace()
    #assert isinstance(timeser,dbs.timeseries)
    #Window size is about 1 second (512 samples is over 1 sec)
    
    Fvect,Pxx = sig.welch(timeser,Fs,window='blackmanharris',nperseg=nperseg,noverlap=noverlap,nfft=nfft)
    FreqReturn = {'F': Fvect,'Pxx': Pxx}
    
    return FreqReturn

def TF_Domain(timeser,fs=422,noverlap=2**10-50):
    #assert isinstance(timeser,dbs.timeseries)
    
    F,T,SG = sig.spectrogram(timeser,nperseg=2**10,noverlap=noverlap,window=sig.get_window('blackmanharris',2**10),fs=fs)
    
    TFreqReturn = {'T': T,'F':F,'SG': SG}
    
    return TFreqReturn

def plot_TF(TFR,ch):
    plt.figure()
    
    plt.pcolormesh(TFR['T'],TFR['F'],10*(TFR['SG'][ch,:,:]))
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    
## BELOW THIS ARE THE STUDY RELATED FUNCTIONS/CLASSES
all_phases = ['A0'+str(num) for num in range(4,0,-1)] + ['B0' + str(num) for num in range(1,5)] + ['C0' + str(num) for num in range(1,10)] + ['C' + str(num) for num in range(10,25)]
def Phase_List(exprs='all'):
    if exprs=='all':
        return all_phases
    elif exprs=='ephys':
        return all_phases[4:]
    elif exprs=='therapy':
        return all_phases[8:]
    elif exprs=='notherapy':
        return all_phases[0:8]
    
    
def get_slope(Pxx,F,params):
    #method to get the fitted polynomial for the range desired
    frange = params['frange']
    linorder = params['linorder']
    
    if isinstance(Pxx,np.ndarray):
        Pxx = {0:Pxx}
    
    out_feats = {keys:0 for keys in Pxx.keys()}
    
    Fidxs = np.where(np.logical_and(F > frange[0],F < frange[1]))
    
    for chans,psd in Pxx.items():
        try:
            logpsd = np.log10(psd[Fidxs])
            logF = np.log10(F[Fidxs])
        except FloatingPointError:
            pdb.set_trace()
            
        fitcoeffs = np.polyfit(logF,logpsd,linorder)
        
        out_feats[chans] = fitcoeffs[-linorder]
    
    #return is going to be a dictionary with same channel keys
    return out_feats
    

## SIMPLE FUNCTIONS
def unity(invar):
    return invar


def displog(values):
    return 10*np.log10(values)

#Variable definitions
#Generalized Features
#Need to regen this based off of the bands up there
feat_dict = {
                'Delta':{'fn':get_pow,'param':(1,4)},
                'Alpha':{'fn':get_pow,'param':(8,14)},
                'Theta':{'fn':get_pow,'param':(4,8)},
                'Beta':{'fn':get_pow,'param':(14,30)},
                'Gamma1':{'fn':get_pow,'param':(30,60)},
                'Gamma2':{'fn':get_pow,'param':(60,100)},
                'Stim':{'fn':get_pow,'param':(129,131)},
                'SHarm':{'fn':get_pow,'param':(31,33)},
                'Clock':{'fn':get_pow,'param':(104.5,106.5)},
                'fSlope':{'fn':get_slope,'param':{'frange':(1,20),'linorder':1}},
                'nFloor':{'fn':get_slope,'param':{'frange':(50,200),'linorder':0}}
            }
feat_order = ['Delta','Theta','Alpha','Beta','Gamma1','fSlope','nFloor']

#Function to go through and find all the features from the PSD structure of dbo
def calc_feats(psdIn,yvect):
    #psdIn is a VECTOR, yvect is the basis vector
    
    feat_vect = []
    for feat in feat_order:
        #dofunc = feat_dict[feat]['fn']
        feat_vect.append(feat_dict[feat]['fn'](psdIn,yvect,feat_dict[feat]['param'])[0])
        #feat_dict[feat] = dofunc['fn'](datacontainer,yvect,dofunc['param'])[0]

    feat_vect = np.array(feat_vect)
    
    return feat_vect

#Convert a feat dict that comes from a get feature function (WHERE IS IT?!)
def featDict_to_Matr(featDict):
    #structure of feat dict is featDict[FEATURE][CHANNEL] = VALUE
    ret_matr = np.array([(featDict[feat]['Left'],featDict[feat]['Right']) for feat in feat_order])
    
    #assert that the size is as expected?
    assert ret_matr.shape == (len(feat_order),2)
    
    return ret_matr
