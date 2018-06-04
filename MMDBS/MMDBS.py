# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:05:48 2016

@author: virati

Description: Library to load in multimodal DBS data, specifically for data handled by Mayberg Lab in ~2014-2016

"""

#Numerical Data libraries
import pandas as pd
import numpy as nm

#Scipy and signals
import scipy.signal as sig
import scipy.io

#Plotting
import matplotlib.pyplot as plt
#import pyqtgraph as pg
import pdb

class signal_container():
    timeseries = []
    interval = []
    
    def __init__(self):
        self.timeseries = timeseries(422)

class timeseries():
    raw_data = []
    proc_data = []
    raw_Fs = 0
    proc_Fs = 0
    modality = []

    def __init__(self,Fs,modality='generic'):
        self.modality = modality
        self.raw_Fs = Fs
        
    def addRaw(self,data):
        #Data should be a numpy array
        self.raw_data = data
        
    def decimate(self,decim=2):
        data_size = self.raw_data.shape;
        raw_subsample_size = data_size[0];
        self.plotdata = nm.zeros((raw_subsample_size,2))
        #plotdata[:,1] = sig.decimate(self.raw_data[:,1],2)
        
        self.plotdata[:,0] = sig.decimate(self.raw_data[:,0],1)
        self.plotdata[:,1] = sig.decimate(self.raw_data[:,2],1)

        
    def pplot(self):
        print('Plotting Data...')
                
        tvect = nm.arange(nm.shape(self.plotdata)[0])

        plt.subplot(1,2,1)
        plt.plot(self.plotdata[:,0])
        plt.subplot(1,2,2)
        plt.plot(self.plotdata[:,1])
        plt.show()
        
    def SGPlot(self):
        
        plt.specgram(self.plotdata[:,0],window=sig.blackmanharris(2**10),NFFT=2**10,Fs=422,noverlap=250)
        
        plt.show()
        
    def osc_pow(self):
        band_lims = [[1,4],[4,8],[8,14],[14,20],[20,30],[30,50],[50,100]]
        bands = {'DC': [0,1], 'Delta': [1,4], 'Theta': [4,8], 'Alpha': [8,14], 'Low Beta': [14,20], 'High Beta': [20,30], 'Low Gamma': [30,50], 'High Gamma': [50,100]}

        for key,value in bands.items():
            print('Band is:' + key)
            b,a = sig.ellip(10,0.1,10,[])
            
            self.featdata[:,key]
        
        

def loadData(fname):
    #Extract extension of fname
    extension = fname[-3:]
    print('Loading ' + extension + ' type')

    if extension == 'txt':
        #Load in the brainradio file
        x = importbrLFP(fname)
        #place the data inside a timeseries object, CUSTOM
        datac = timeseries(422,modality='brLFP')
        datac.addRaw(x)
    elif extension == 'mat':
        x = importhdEEG(fname)
        datac = timeseries(1000,modality='hdEEG')
        datac.addRaw(x)
        
    print('Done Loading ' + fname)
    return datac
    #datac.pplot()

def importbrLFP(fname):
    data = pd.read_csv(fname,header=None)
    
    print('Loading Raw brLFP Data...')
    return nm.array(data)
    

def importhdEEG(fname):
    print('Loading Raw hdEEG Data...')
    #import h5py
    #f = h5py.File(fname,'r')
    #ipdb.set_trace()
    
    #data = f.get('x')
    
    #structfields = scipy.io.whosmat(fname)
    loadc = scipy.io.loadmat(fname)
    for fit in loadc.keys():
        if fit[0:3] == 'DBS':
            print(fit)
            data = loadc[fit]
    
    return nm.array(data)