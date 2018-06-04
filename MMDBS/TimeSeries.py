#Vineet Tiruvadi
# Library file for TimeSeries analysis

import numpy as np
import pandas as pd
import scipy.signal as sig
import matplotlib.pyplot as plt

from collections import defaultdict, OrderedDict

#Debug stuff
import pdb

band_dict = {'Delta':(1,4),'Theta':(4,8),'Alpha':(8,14),'Beta':(14,30),'Beta*':(14,20),'Beta+':(25,30),'Gamma':(30,50),'Gamma+':(50,70)}
do_bands = ['Delta','Theta','Alpha','Beta*']

def band_structs():
    return band_dict

class timeser:
    data = np.array([])
    n_chann = 0
    Fs = 0
    t_vect = np.array([])

    filt_SG = np.array([])
    
    segments = defaultdict(list)
    
    def __init__(self,data,Fs):
        self.data = data
        self.Fs = Fs
        self.n_chann = data.shape[1]
        self.t_vect = np.linspace(0,data.shape[0]/Fs,data.shape[0]) #this assumes data is (channels,observations)
        
    
    def view_raw_ts(self,title='TimeSeries',channs=range(1)):
        plt.figure()
        plt.plot(self.data[:,channs])
        plt.title(title)
        plt.show()
        
    def raw_ts(self):
        return [self.t_vect,self.data]
    
    def view_raw_hist(self,chann):
        plt.figure()
        for ii in chann:        
            plt.hist(self.data[:,ii],bins=50)
        plt.show()

    def define_segments(self,segs):
        self.segs = segs
    
    def seg_ts(self,seg_name):
        return self.data[self.segs[seg_name][0]:self.segs[seg_name][1],:]
            
    def stim_segment(self):
        #find the stim sites via 130Hz filtering
        
        #Segment based on a set policy off of stim sites
        
        #construct the dictionary
        for ss in range(stim_seg_n):
            self.segments[ss] = 5
    
    def view_tf(self,title='Spectrogram',channs=np.arange(1)):
        for ii in channs:
            plt.figure()
            F,T,SG,bands = self.compute_tf()
            plt.pcolormesh(T,F,10*np.log10(np.abs(SG[ii])))
            plt.title('Channel ' + str(ii))            
            plt.ylim((0,50))
        
    def PCA_Comps(self):
        for ii in channs:
            F,T,SG = self.compute_tf()
            runSG.append(SG)
        bigSG = np.array(runSG)
        
    def compute_tf(self,title='Spectrogram',channs=np.arange(2)):
        SG = []
        BANDS = []
        for ii in channs:
            #Pxx,f,b,i = plt.specgram(self.data[:,ii],NFFT=1024,Fs=self.Fs,noverlap=512,cmap=plt.cm.jet,hold=True)
            f,t,Sxx = sig.spectrogram(self.data[:,ii],nperseg=512,noverlap=0,window=sig.get_window('blackmanharris',512),fs=self.Fs)            
            SG.append(Sxx)
            
            BANDS.append(self.compute_bands(f,Sxx))
        return f,t,SG,BANDS

    def compute_bands(self,f,SG):
        BandPow = defaultdict(dict)
        for bb, band in enumerate(band_dict):
            f_idxs = np.where(np.logical_and(f > band_dict[band][0],f < band_dict[band][1]))
            BandPow[band] = np.mean(np.squeeze(SG[f_idxs,:]),0)
        return BandPow
            
    def apply_filt(self,band):
        for ii in range(self.n_chann):
            f,t,Sxx = sig.spectrogram(self.data[:,ii],self.Fs,window=sig.get_window('blackmanharris',512),nperseg=512,nfft=1024)
            filt_SG = Sxx[np.where(np.logical_and(f >= band[0],f <= band[1])),:]        
        return filt_SG
    
    def compute_PSD(self,chann):
        f,P = sig.welch(self.data[:,chann],self.Fs)
        print(str(self.Fs))
        return f,P
        
    def ret_n_chann(self):
        return self.n_chann
    
def import_BR(file,Fs=422,snip=(0,0)):
    rawdata = np.array(pd.read_csv(file,sep=',',header=None))
    data = rawdata[:,[0,2]]
    #generate a time vector based off of the above
    tvect = np.linspace(0,data.shape[0]/Fs,data.shape[0])
    
    if snip == (0,0):
        data_snip = data[0:-1,:]
    else:
        data_snip = data[snip[0]*Fs:snip[1]*Fs,:]
        
    data_struct = timeser(data_snip,Fs)    
    return data_struct
    
def Diff_PSD(T2,T1,chann):
    f,Psd2 = T2.compute_PSD(chann)
    f,Psd1 = T1.compute_PSD(chann)
    
    return f,10*np.log10(np.abs(Psd2)) - 10*np.log10(np.abs(Psd1))