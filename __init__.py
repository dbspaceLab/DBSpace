#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 13:44:38 2018

@author: virati
DBS_Osc File for the analytical methodologies used in ALL Vineet project related analyses

Copyright (C) 2018 Vineet Ravi Tiruvadi

This is the primary file with the methods needed for the rest of the libraries. Used to be called "DBS_Osc" or "DBSOsc".
    
"""
print('Using DBSpace LATEST')
import numpy as np
import pandas as pd
from collections import defaultdict
import scipy.signal as sig
import random

# IF you want to do OR related analyses, this needs to be uncommented
#from brpylib import NsxFile

import pdb

import matplotlib.pyplot as plt
plt.rcParams['image.cmap'] = 'jet'

np.seterr(divide='raise')

all_pts = ['901','903','905','906','907','908']

# This is our map for the electrodes that each patient has for ONTarget and OFFTarget
Etrode_map = {'OnT':{'901':(2,1),'903':(2,2),'905':(2,1),'906':(2,2),'907':(1,1),'908':(2,1)},'OffT':{'901':(1,2),'903':(1,1),'905':(1,2),'906':(1,1),'907':(2,2),'908':(1,2)}}

#%%
# BlackRock Methods

def load_or_file(fname,**kwargs):
    #nsx_file = NsxFile(fname)
    
    arg_list = ['elec_ids','start_time_s','data_time_s','downsample','plot_chan']
    
    for aa in arg_list:
        print(arg_list[aa])

    # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
    cont_data = nsx_file.getdata(elec_ids, start_time_s, data_time_s, downsample)
    
    # Close the nsx file now that all data is out
    nsx_file.close()
    
    return cont_data

def load_or_dict(fname,sec_win=(0,10),channels=[0]):
    pass

#%%
# BRAIN RADIO METHODS
#Method to load in brainradio file
def load_br_file(fname):
    return np.array(pd.read_csv(fname,sep=',',header=None))

# Load BR file and return it as a dictionary
def load_BR_dict(fname,sec_offset=10,channels=['Left','Right']):
    txtdata = load_br_file(fname)[:,[0,2]]
    
    return {chann:txtdata[-(422*sec_offset):-1,cidx] for cidx,chann in enumerate(channels)}

def gen_T(inpX,Fs=422,nfft=2**10):
    outT = defaultdict(dict)
    for chann in inpX.keys():
        outT[chann] = {'T':np.linspace(0,inpX[chann].shape[0]/Fs,inpX[chann].shape[0]),'V':inpX[chann]}
        
    return outT

''' gen_psd outputs a PSD, not a LogPSD '''
''' This function WRAPS F_Domain'''
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
        
        if inpX[chann].ndim > 1:
            for seg in range(inpX[chann].shape[-1]):
                
                psd = F_Domain(inpX[chann][:,seg].squeeze(),Fs=Fs,nfft=nfft)['Pxx']
                
                    
                fmatr[seg,:] = psd
        else:
            psd = F_Domain(inpX[chann],Fs=Fs,nfft=nfft)['Pxx']
            fmatr = psd
                    
        outPSD[chann] = fmatr.squeeze()


        #do polysub here
    if polyord != 0:
        print("Polynomial Correcting Stack")
        outPSD = poly_subtr(outPSD,np.linspace(0,Fs/2,nfft/2 + 1))
        
    #Return here is a dictionary with Nchann keys
    return outPSD

def cmedian(inArray,axis=-1):
    return  np.median(np.real(inArray),axis=axis) + 1j * np.median(np.imag(inArray),axis=axis)

def l2_pow(x):
    return np.sqrt(np.sum(x**2))

#This function takes a PSD and subtracts out the PSD's fourth order polynomial fit
#I THINK this is only used in EEG cortical signatures
def poly_subtrEEG(inpPSD,fVect,order=5):
    #What's our feature/frequency vector?
    fix_psd = defaultdict()
    for chann in inpPSD.keys():
        fix_psd[chann] = []
        inpPSD[chann] = inpPSD[chann].reshape(-1,1).T
        #This should now be (513(nfft) x segments).T
        
        
        #SEGMENTS x PSD
        postpsd_matr = np.zeros((inpPSD[chann].shape[0],inpPSD[chann].shape[1]))
        
        for seg in range(inpPSD[chann].shape[0]):
            curr_psd = 10*np.log10(inpPSD[chann][seg,:])
            #print(fVect.shape)
            #print(curr_psd.shape)
            polyCoeff = np.polyfit(fVect,curr_psd,order)
            
            polyfunc = np.poly1d(polyCoeff)
            polyitself = polyfunc(fVect)
            
            postpsd_matr[seg,:] = curr_psd - polyitself
            #fix_poly[seg,:] = polyCoeff
        
        fix_psd[chann] = 10**(postpsd_matr/10).T
        
    
    return fix_psd, polyitself

def poly_subtrLFP(fvect,inp_psd,polyord=4):
    #fvect HAS to be a vector function
    # This function takes in a raw PSD, Log transforms it, poly subtracts, and then returns the unloged version.
    #log10 in_psd first

    try:
        log_psd = 10*np.log10(inp_psd)
        pfit = np.polyfit(fvect,log_psd,polyord)
        pchann = np.poly1d(pfit)
        
        bl_correction = pchann(fvect)
    except Exception as e:
        print(e)
        pdb.set_trace()
        
    return 10**((log_psd - bl_correction)/10), pfit

def poly_SG(inSG,fVect,order=4):
    out_sg = np.zeros_like(inSG)
    
    for seg in range(inSG.shape[1]):
        inpsd = 10*np.log10(inpPSD[chann][seg,:])
        polyCoeff = np.polyfit(fVect,inpsd,order)
        polyfunc = np.poly1d(polyCoeff)
        polyitself = polyfunc(fVect)
        out_sg[:,seg] = 10**((curr_psd - polyitself)/10)
        
    return out_sg
    
    
def gen_SG(inpX,Fs=422,nfft=2**10,plot=False,overlap=True):
    outSG = defaultdict(dict)
    for chann in inpX.keys():
        if overlap==True:
            outSG[chann] = TF_Domain(inpX[chann])
        else:
            outSG[chann] = TF_Domain(inpX[chann],noverlap=0,nperseg=422*2)
    
    if plot:
        plot_TF(outSG,chs=inpX.keys())
    
    return outSG

'''
Single method that takes an inputX.channels dictionary and outputs the oscillatory state, no fuss
'''
def osc_state(inpX):
    pass

''' Return to us the power in an oscillatory feature'''
def get_pow(Pxx,F,frange,cmode=np.median):
    #Pxx is a dictionary where the keys are the channels, the values are the [Pxx desired]
    #Pxx is assumed to NOT be log transformed, so "positive semi-def"

    #check if Pxx is NOT a dict
    if isinstance(Pxx,np.ndarray):
        #JUST ADDED THIS
        chann_order = range(Pxx.shape[0])
        Pxx = {ch:Pxx[ch,:] for ch in chann_order}
        
        #ThIS WAS WORKING BEFORE
        #Pxx = {0:Pxx}
    elif len(Pxx.keys()) > 2:
        chann_order = np.arange(0,257)
    else:
        chann_order = ['Left','Right']
        
    
    #find the power in the range of the PSD
    #Always assume PSD is a dictionary of channels, and each value is a dictionary with Pxx and F
    
    #frange should just be a tuple with the low and high bounds of the band
    out_feats = {keys:0 for keys in Pxx.keys()}
    
    Fidxs = np.where(np.logical_and(F > frange[0],F < frange[1]))[0]
    
    #for chans,psd in Pxx.items():
    for cc,chann in enumerate(chann_order):
        #let's make sure the Pxx we're dealing with is as expected and a true PSD
        assert (Pxx[chann] > 0).all()
    
        #if we want the sum
        #out_feats[chans] = np.sum(psd[Fidxs])
        #if we want the MEDIAN instead
        
        #log transforming this makes sense, since we find the median of the POLYNOMIAL CORRECTED Pxx, which is still ALWAYS positive
        try:
            out_feats[chann] = 10*np.log10(cmode(Pxx[chann][Fidxs]))
        except Exception as e:
            print(e)
            pdb.set_trace()
    
    #return is going to be a dictionary with same elements

    return out_feats


def get_ratio(Pxx,F,f_r_set,cmode=np.median):
    bandpow = [None] * len(f_r_set)
    #first get the power for each of the individual bands
    for bb, frange in enumerate(f_r_set):
        bandpow[bb] = get_pow(Pxx,F,frange,cmode=cmode)
        
    
    ret_ratio = {ch:bandpow[1][ch]/bandpow[0][ch] for ch in bandpow[0].keys()}
    return ret_ratio


'''Rotating into frequency-space'''
''' WE USE THIS FOR EEG classifier analysis'''
def F_Domain(timeser,nperseg=512,noverlap=128,nfft=2**10,Fs=422):
    
    #assert isinstance(timeser,dbs.timeseries)
    #Window size is about 1 second (512 samples is over 1 sec)
    
    #what are the dimensions of the timeser we're dealing with?
    
    Fvect,Pxx = sig.welch(timeser,Fs,window='blackmanharris',nperseg=nperseg,noverlap=noverlap,nfft=nfft)
    
    FreqReturn = {'F': Fvect,'Pxx': Pxx}
    
    return FreqReturn

def TF_Domain(timeser,fs=422,nperseg=2**10,noverlap=2**10-50):
    #raise Exception
    #assert isinstance(timeser,dbs.timeseries)
    F,T,SG = sig.spectrogram(timeser,nperseg=nperseg,noverlap=noverlap,window=sig.get_window('blackmanharris',nperseg),fs=fs)
    
    TFreqReturn = {'T': T,'F':F,'SG': SG}
    
    return TFreqReturn

'''
F and TF domain plotting 
'''
def plot_TF(TFR,chs=['Left','Right']):
    plt.figure()
    for cc,chann in enumerate(chs):
        plt.subplot(1,len(chs),cc+1)
        aTFR = TFR[chann]
        
        plt.pcolormesh(aTFR['T'],aTFR['F'],10*np.log10(aTFR['SG']))
        plt.xlabel('Time')
        plt.ylabel('Frequency')


#This function plots the median/mean across time of the TF representation to get the Frequency representation
#Slightly different than doing the Welch directly
def plot_F_fromTF(TFR,chs=['Left','Right']):
    plt.figure()
    for cc,chann in enumerate(chs):
        plt.subplot(1,len(chs),cc+1)
        aTFR = TFR[chann]
        
        for ss in range(aTFR['SG'].shape[1]):
            plt.plot(aTFR['F'],10*np.log10(aTFR['SG'])[:,ss],alpha=0.1)
            
        plt.plot(aTFR['F'],np.median(10*np.log10(aTFR['SG']),axis=1))
        plt.xlabel('Frequency')
        plt.ylabel('Power')
        

def plot_T(Tser):
    plt.figure()
    for cc,chann in enumerate(Tser.keys()):
        plt.subplot(1,len(Tser.keys()),cc+1)
        aT = Tser[chann]
        plt.plot(aT['T'],aT['V'])
        
        plt.xlabel('Time')
        
''' Return the phases and the labels associated with them'''
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
    
''' Function that retrieves the PSD slope'''
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

'''
FEATURE DEFINITIONS HERE
Important block, this is where we define our features
'''
#Variable definitions
#Generalized Features
#Need to regen this based off of the bands up there
feat_dict = {
                'Delta':{'fn':get_pow,'param':(1,4)},
                'Alpha':{'fn':get_pow,'param':(8,14)},
                'Theta':{'fn':get_pow,'param':(4,8)},
                'Beta*':{'fn':get_pow,'param':(14,20)},
                'Beta':{'fn':get_pow,'param':(14,30)},
                'Gamma1':{'fn':get_pow,'param':(35,60)},
                'Gamma2':{'fn':get_pow,'param':(60,100)},
                'Gamma':{'fn':get_pow,'param':(30,100)},
                'Stim':{'fn':get_pow,'param':(129,131)},
                'SHarm':{'fn':get_pow,'param':(30,34)}, #Secondary Harmonic
                'THarm':{'fn':get_pow,'param':(62,66)}, #Tertiary Harmonic
                'Clock':{'fn':get_pow,'param':(104.5,106.5)},
                'fSlope':{'fn':get_slope,'param':{'frange':(1,20),'linorder':1}},
                'nFloor':{'fn':get_slope,'param':{'frange':(50,200),'linorder':0}},
                'GCratio':{'fn':get_ratio,'param':((30,34),(62,66))}
            }


feat_order = ['Delta','Theta','Alpha','Beta*','Gamma1']#,'fSlope','nFloor']

#Function to go through and find all the features from the PSD structure of dbo
def calc_feats(psdIn,yvect,dofeats='',modality='eeg',compute_method='median'):
    #psdIn is a VECTOR, yvect is the basis vector
    if dofeats == '':
        dofeats = feat_order
    
    if modality == 'eeg':
        ch_list = np.arange(0,257)
    elif modality == 'lfp':
        ch_list = ['Left','Right']
    
    feat_vect = []
    for feat in dofeats:
        #print(feat_dict[feat]['param'])
        #dofunc = feat_dict[feat]['fn']
        if compute_method == 'median':
            computed_featinspace = feat_dict[feat]['fn'](psdIn,yvect,feat_dict[feat]['param'])
        elif compute_method == 'mean':
            computed_featinspace = feat_dict[feat]['fn'](psdIn,yvect,feat_dict[feat]['param'],cmode=np.mean)
        
        cfis_matrix = [computed_featinspace[ch] for ch in ch_list]
        feat_vect.append(cfis_matrix)
        #feat_dict[feat] = dofunc['fn'](datacontainer,yvect,dofunc['param'])[0]

    feat_vect = np.array(feat_vect).squeeze()
    
    return feat_vect, dofeats


#Convert a feat dict that comes from a get feature function (WHERE IS IT?!)
def featDict_to_Matr(featDict):
    #structure of feat dict is featDict[FEATURE][CHANNEL] = VALUE
    ret_matr = np.array([(featDict[feat]['Left'],featDict[feat]['Right']) for feat in feat_order])
    
    #assert that the size is as expected?
    #should be number of feats x number of channels!
    assert ret_matr.shape == (len(feat_order),2)
    
    return ret_matr


''' Function to plot bands since this has been annoying every time I've had to recode the thing from scratch'''
def plot_bands(bandM,bandLabels):
    plt.figure()
    for cc in bandM:
        plt.bar(cc)
        plt.xticks(range(len(cc)))
        plt.xticklabels(bandLabels)

'''
SIMPLE FUNCTIONS
Should probably be split off into a utils.py file...
'''
def unity(invar):
    return invar


def displog(values):
    return 10*np.log10(values)

''' Powerful structure that saves my ass'''
def nestdict():
    return defaultdict(nestdict)


''' Using Jack-knife for the median calculation'''
def jk_median(inX,niter=100):
    #assume first dim is SEGMENTS/OBSERVATIONS
    med_vec = []
    for ii in range(niter):
        choose_idxs = random.sample(range(0,inX.shape[0]),np.floor(inX.shape[0]/2).astype(np.int))
        med_vec.append(np.median(inX[choose_idxs,:],axis=0))
        
    return np.array(med_vec)

''' Higher order measures here'''
# Make a coherence generation function
def gen_coher(inpX,Fs=422,nfft=2**10,polyord=0):
    print('Starting a coherence run...')
    outPLV = nestdict()
    outCSD = nestdict()
    
    #What's our goddamn fvector
    fvect = np.linspace(0,Fs/2,nfft/2+1)

    for chann_i in inpX.keys():
        print(chann_i)
        for chann_j in inpX.keys():
            csd_ensemble = np.zeros((inpX[chann_i].shape[1],len(feat_order)),dtype=complex)
            plv = np.zeros((inpX[chann_i].shape[1],len(feat_order)))
                
            for seg in range(inpX[chann_i].shape[1]):
                #First we get the cross spectral density
                csd_out = sig.csd(inpX[chann_i][:,seg],inpX[chann_j][:,seg],fs=Fs,nperseg=512)[1]
                
                #normalize the entire CSD for the total power in input signals
                norm_ms_csd = np.abs(csd_out) / np.sqrt(l2_pow(inpX[chann_i][:,seg]) * l2_pow(inpX[chann_j][:,seg]))
                
                #Are we focusing on a band or doing the entire CSD?

                for bb,band in enumerate(feat_order):
                    #what are our bounds?
                    band_bounds = feat_dict[band]['param']
                    band_idxs =  np.where(np.logical_and(fvect >= band_bounds[0],fvect <= band_bounds[1]))
                    csd_ensemble[seg,bb] = cmedian(csd_out[band_idxs])
                    plv[seg,bb] = np.max(norm_ms_csd[band_idxs])
                    
                # Below brings in the entire csd, but this is dumb
                #csd_ensemble[seg] = csd_out

                #Compute the PLV
                
            # Here we find the median across segments
            #outCSD[chann_i][chann_j] = cmedian(csd_ensemble,axis=0)
            outCSD[chann_i][chann_j] = csd_ensemble
            
            #Compute the normalized coherence/PLV
            outPLV[chann_i][chann_j] = plv
            #outPLV[chann_i][chann_j] = np.median(plv,axis=0)
            #if chann_j == 176: pdb.set_trace()
            #if chann_i > 10: pdb.set_trace()
            ## PLV abs EEG -> 
            ## Coherence value
            
    return outCSD, outPLV
    
''' PCA function for merging rPCA results into final DSC 
This is borrowed from online stackexchange somewhere, not the most efficient way to do it
'''
def simple_pca(data,numComps=None):
    m,n = data.shape
    data -= data.mean(axis=0)
    R = np.cov(data,rowvar=False)
    evals,evecs = np.linalg.eigh(R)
    idx=np.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    evals = evals[idx]
    
    if numComps is not None:
        evecs = evecs[:,:numComponents]
    
    return np.dot(evecs.T,data.T).T,evals,evecs