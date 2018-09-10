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


# Finally, we'll have a class that performs a fixed preprocessing pipeline for a single recording
class BR_Pipeline:
    def __init__(self):
        pass