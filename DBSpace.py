#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 17:47:02 2018

@author: virati
DBSpace Library

Copyright (C) 2018 Vineet Ravi Tiruvadi
Main library for DBSpace
Has basic loading functions and large-scale dataframes
"""

#import sys
#sys.path.append('/home/virati/Dropbox/projects/Research/MDD-DBS/Ephys/IntegratedAnalysis/')

import numpy as np
import pandas as pd

from collections import defaultdict
import scipy.signal as sig
import matplotlib.pyplot as plt


# Basic loading functions for text file based recordings (such as those used by the Activa PC+S)

def load_br(fname,from_end=10,fs=422):
    #Bring in raw data
    rawdata = np.array(pd.read_csv(fname,sep=',',header=None))
    
    #calculate what the sample number is for the from_end clipping: this stage is to avoid the first ~5 second settling artifact
    start_sample = -422*from_end
    #The PC+S has channels 1 and 2 in the first and third columns.
    txtdata = rawdata[from_end:-1,[0,2]]
    
    #make a dictionary structure for the two channels in the PC+S
    Y = {'Left':txtdata[:,0],'Right':txtdata[:,1]}
    
    return Y

# This function will parse out the filename from BrainRadio Recordings
def parse_filename(fname):
    pass

def poly_subtr(inPSD,order=4):
    #inPSD should be a dictionary with the PSD and fvect in it
    fixed_psd = {chann:[] for chann in inPSD.keys()}
    
    #Main loop to go through channels
    # KEEP IN MIND, .keys() does not guarantee order
    for chann in inPSD.keys():
        active_PSD = 10*np.log10(inPSD[chann][seg,:])


def feat_extract(inY):

    
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


#%%
def plot_ts(y_dict,x,fs,x_offset=0):
    plt.figure()
    for cc in y_dict.keys():
        plt.subplot(n,1,cc)
        plt.plot(x,y_dict[cc])

#%%
# Finally, we'll have a class that performs a fixed preprocessing pipeline for a single recording
class BR_Pipeline:
    def __init__(self):
        pass
    
#A class that wraps a BR recording
class BRec:
    def __init__(self,fname=''):
       self.filename = fname #what file is this recording associated with?
       self.Tser = {} #dictionary with keys for channels
       self.Data = {} #
       
       # Apply a partition on the recording
       self.partition(profile='chronic')
       
       if fname != '':
           self.load_recording(fname)
       
    #wrapper for external call to main BR loading function
    def load_recording(self,fname,partition=-1):
        data_dict = load_br(fname)
        
        
    # All recordings can be PARTITIONED, with the default partition being two partitions
    # 0 - settling partition; 1 - remaining
    
    @property
    def partition(self,profile='chronic'):
        self.partition
    
    @partition.setter
    def partition(self):
        pass    
        
        
    @property
    def preproc(self):
        return self.__preproc
    
    @preproc.setter
    def preproc(self):
        self.__preproc = []
      
    # Display function for the raw timeseries
    def view_ts(self):
        plt.figure()
        plt.plot(self.Tser)