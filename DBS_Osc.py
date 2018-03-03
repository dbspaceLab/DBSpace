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

import matplotlib.pyplot as plt
plt.rcParams['image.cmap'] = 'jet'

np.seterr(divide='raise')

#Definitions for oscillatory features -> Oscillatory Vector
band_dict = {'Delta':(1,4), 'Theta':(4,8), 'Alpha':(8,14), 'Beta*':(14,20), 'Beta+':(25,30), 'Gamma1':(30,50), 'Gamma2':(50,80), 'Gamma3':(80,100), 'Stim':(128,132)}
band_order = ['Delta','Theta','Alpha','Beta*','Beta+','Gamma1','Gamma2','Gamma3','Stim']

all_pts = ['901','903','905','906','907','908']

#Method to load in brainradio file
def load_br_file(fname):
    rawdata = np.array(pd.read_csv(fname,sep=',',header=None))
    return rawdata

def load_BR_dict(fname,sec_end=10):
    txtdata = load_br_file(fname)[:,[0,2]]
    
    X = {'Left':txtdata[-(422*sec_end):-1,0],'Right':txtdata[-(422*sec_end):-1,1]}
    
    return X

def gen_T(inpX,Fs=422,nfft=2**10):
    outT = defaultdict(dict)
    for chann in inpX.keys():
        outT[chann] = {'T':np.linspace(0,inpX[chann].shape[0]/Fs,inpX[chann].shape[0]),'V':inpX[chann]}
        
    return outT
        

#ALERT    
#gen_psd has gotten super complicated.... need to check if it still works with DSV and Regression stuff
    
def gen_psd(inpX,Fs=422,nfft=2**10,polyord=0):
    #inp X is going to be assumed to be a dictionary with different keys for different channels
    outPSD = defaultdict(dict)
    outPoly = defaultdict(dict)
    #assume input is time x seg
    for chann in inpX.keys():
        #The return here is a dictionary with two keys: F and PSD
        #check the size of the matrix now; it could be that we need to do this many times for each "segment"
        fmatr = np.zeros((inpX[chann].shape[-1],int(nfft/2)+1))
        polysub = np.zeros((inpX[chann].shape[-1],polyord+1))
        
        for seg in range(inpX[chann].shape[-1]):
            psd = F_Domain(inpX[chann][:,seg].squeeze(),Fs=Fs,nfft=nfft)['Pxx']
            
            if polyord != 0:
                intermed = 10*np.log10(psd)
                fVect = np.linspace(0,Fs/2,int(nfft/2)+1)
                
                polyCoeff = np.polyfit(fVect,psd,polyord)
                polyfunc = np.poly1d(polyCoeff)
                polyitself = polyfunc(fVect)
                
                
                fmatr[seg,:] = 10**((intermed - polyitself)/10)
                polysub[seg,:] = polyCoeff
                
            else:
                fmatr[seg,:] = psd
            
        
        outPSD[chann] = fmatr.squeeze()
        outPoly[chann] = polysub.squeeze()

    #Return here is a dictionary with Nchann keys
    return outPSD,outPoly

#This function takes a PSD and subtracts out the PSD's fourth order polynomial fit
def poly_subtr(inpPSD,fVect,order=4):
    #What's our feature/frequency vector?
    
    for chann in inpPSD.keys():
        for seg in range(inpX[chann].shape[-1]):
            polyCoeff = np.polyfit(fVect,inpPSD[chann][seg,:],order)
            polyfunc = np.poly1d(polyCoeff)
            polyitself = polyfunc(fVect)
            
            fix_psd[seg,:] = inpPSD[chann][seg,:] - polyfit
            fix_poly[seg,:] = polyCoeff
        
        
    
    return fix_psd,fix_poly
    

def gen_SG(inpX,Fs=422,nfft=2**10,plot=False):
    outSG = defaultdict(dict)
    for chann in inpX.keys():
        outSG[chann] = TF_Domain(inpX[chann])
    
    if plot:
        plot_TF(outSG,chs=inpX.keys())
    
    return outSG

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
    
    #what are the dimensions of the timeser we're dealing with?
    
    Fvect,Pxx = sig.welch(timeser,Fs,window='blackmanharris',nperseg=nperseg,noverlap=noverlap,nfft=nfft)
    FreqReturn = {'F': Fvect,'Pxx': Pxx}
    
    return FreqReturn

def TF_Domain(timeser,fs=422,nperseg=2**10,noverlap=2**10-50):
    #assert isinstance(timeser,dbs.timeseries)
    F,T,SG = sig.spectrogram(timeser,nperseg=nperseg,noverlap=noverlap,window=sig.get_window('blackmanharris',nperseg),fs=fs)
    
    TFreqReturn = {'T': T,'F':F,'SG': SG}
    
    return TFreqReturn

def plot_TF(TFR,chs=['Left','Right']):
    plt.figure()
    for cc,chann in enumerate(chs):
        plt.subplot(1,len(chs),cc+1)
        aTFR = TFR[chann]
        
        plt.pcolormesh(aTFR['T'],aTFR['F'],10*np.log10(aTFR['SG']))
        plt.xlabel('Time')
        plt.ylabel('Frequency')
        

def plot_T(Tser):
    plt.figure()
    for cc,chann in enumerate(Tser.keys()):
        plt.subplot(1,len(Tser.keys()),cc+1)
        aT = Tser[chann]
        plt.plot(aT['T'],aT['V'])
        
        plt.xlabel('Time')
        
    
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
                'SHarm':{'fn':get_pow,'param':(31,33)}, #Secondary Harmonic
                'THarm':{'fn':get_pow,'param':(63,66)}, #Tertiary Harmonic
                'Clock':{'fn':get_pow,'param':(104.5,106.5)},
                'fSlope':{'fn':get_slope,'param':{'frange':(1,20),'linorder':1}},
                'nFloor':{'fn':get_slope,'param':{'frange':(50,200),'linorder':0}}
            }
feat_order = ['Delta','Theta','Alpha','Beta','Gamma1']#,'fSlope','nFloor']

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
    #should be number of feats x number of channels!
    assert ret_matr.shape == (len(feat_order),2)
    
    return ret_matr






#%%
## SIMPLE FUNCTIONS
def unity(invar):
    return invar


def displog(values):
    return 10*np.log10(values)


def nestdict():
    return defaultdict(nestdict)