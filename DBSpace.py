#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 17:47:02 2018

@author: virati
DBSpace Library

Copyright (C) 2018 Vineet Ravi Tiruvadi
"""

#import sys
#sys.path.append('/home/virati/Dropbox/projects/Research/MDD-DBS/Ephys/IntegratedAnalysis/')

import DBS_Osc as dbo
from collections import defaultdict
import scipy.signal as sig
import matplotlib.pyplot as plt

#Uses DBSOsc for the analytical methods
class timeseries:
    def __init__(self,x,t):
        self.ts = defaultdict(dict)
        self.ts['x'] = x
        self.ts['t'] = t
        self.ts['fs'] = x.shape[0]/(t[-1] - t[0])
        
    def transforms(self,domains=['F','TF']):
        for dd,dom in enumerate(domains):
            if dom == 'F':
                #self.F['F'],self.F['PSD'] = sig.welch(self.ts['x'])
                #f,P[:,ii] = sig.welch(x_in[:,ii],fs,window='blackmanharris',nperseg=512,noverlap=128,nfft=nfft)
                self.Freq = dbo.F_Domain(self.ts)
            elif dom == 'TF':
                #self.TF['F'],self.TF['T'],self.TF['SG'] = sig.spectrogram(self.ts['x'],nperseg=512,noverlap=256,window=sig.get_window('blackmanharris',512),fs=fs)
                self.TimeFreq = dbo.TF_Domain(self.ts)
            
        
    def display(self,domains=['T','F','TF']):
        plt.figure()
        numdoms = len(domains)
        
        for dd,dom in enumerate(domains):
        
            plt.subplot(numdoms,1,dd)
            
            if dom == 'T':
                plt.plot(self.ts['t'],self.ts['x'])
            elif dom == 'TF':
                pass
        
