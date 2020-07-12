#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 25 17:41:48 2018

@author: virati
Streaming Class
"""
import scipy.io as scio
import numpy as np
import pandas as pds
from collections import defaultdict
import scipy.signal as sig
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams['image.cmap'] = 'jet'


import DBSpace.control.neigh_mont as neigh_mont

import sys
sys.path.append('/home/virati/Dropbox/projects/Research/MDD-DBS/Ephys/DBSpace/')
import DBSpace as dbo
from DBSpace import nestdict
from DBSpace.visualizations import EEG_Viz

import pdb

import pickle

from sklearn import mixture
from sklearn.decomposition import PCA
from sklearn import svm



Targeting = defaultdict(dict)
Targeting['All'] = {
        '901':{
                'OnT':{
                        'fname':'/home/extend/MDD_Data/hdEEG/Continuous/ALLMATS/DBS901_E52_On_Target_20151030_015625.mat',
                        'lfp':'/home/virati/MDD_Data/BR/901/Session_2014_05_16_Friday/DBS901_2014_05_16_17_10_31__MR_0.txt',
                        'epochs':{'Bilat':(600,630),'PreBilat':(500,530)}},
                'OffT':{
                        'fname':'/home/extend/MDD_Data/hdEEG/Continuous/ALLMATS/DBS901_E52_Off_Target_20151030_022924.mat',
                        'lfp':'/home/virati/MDD_Data/BR/901/Session_2014_05_16_Friday/DBS901_2014_05_16_16_25_07__MR_0.txt',
                        'epochs':{'Bilat':(600,630),'PreBilat':(480,510)},
                'Volt':{}}},
        '903':{
                'OnT':{
                        'fname':'',
                        'lfp':'/home/virati/MDD_Data/BR/903/Session_2014_09_03_Wednesday/DBS903_2014_09_03_14_16_57__MR_0.txt',
                        'epochs':{'Bilat':(550,580),'PreBilat':(501,531)}},
                'OffT':{
                        'fname':'',
                        'lfp':'/home/virati/MDD_Data/BR/903/Session_2014_09_04_Thursday/DBS903_2014_09_04_12_53_09__MR_0.txt' ,
                        'epochs':{'Bilat':(550,580),'PreBilat':(501,531)},
                'Volt':{}}},
        '905':{
                'OnT':{
                        'fname':'/home/extend/MDD_Data/hdEEG/Continuous/ALLMATS/DBS905_TurnOn_Day1_onTARGET_20150928_015403.mat',
                        'lfp':'/home/virati/MDD_Data/BR/905/Session_2015_09_28_Monday/Dbs905_2015_09_28_13_51_42__MR_0.txt',
                        'epochs':{'Bilat':(610,640),'PreBilat':(561,591)}},
                'OffT':{
                        'fname':'/home/extend/MDD_Data/hdEEG/Continuous/ALLMATS/DBS905_TurnOn_OffTargetStims_20150929_123449.mat',
                        'lfp':'/home/virati/MDD_Data/BR/905/Session_2015_09_29_Tuesday/Dbs905_2015_09_29_12_32_47__MR_0.txt' ,
                        'epochs':{'Bilat':(610,640),'PreBilat':(561,591)}},
                'Volt':{}},
        '906':{
                'OnT':{
                        #'fname':'/home/virati/MDD_Data/hdEEG/Continuous/DS500/DBS906_TurnOn_Day1_Sess1_20150827_024013_tds.mat'
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/Targeting/B04/DBS906/DBS906_TurnOn_Day1_Sess1_20150827_024013_OnTarget.mat',
                        'lfp':'/home/virati/MDD_Data/BR/906/Session_2015_08_27_Thursday/DBS906_2015_08_27_15_10_44__MR_0.txt',
                        'epochs':{'Bilat':(610,640),'PreBilat':(561,591)}},
                'OffT':{
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/Targeting/B04/DBS906/DBS906_TurnOn_Day1_Sess2_20150827_041726_OffTarget.mat',
                        'lfp':'/home/virati/MDD_Data/BR/906/Session_2015_08_27_Thursday/DBS906_2015_08_27_16_20_23__MR_0.txt',
                        'epochs':{'Bilat':(610,640),'PreBilat':(561,591)},
                        },
                
                'Volt':{
                        #
                        #'fname':'/home/virati/MDD_Data/hdEEG/Continuous/ALLMATS/DBS906_TurnOn_Day2_Sess3_Sess4_20150828_043231_VoltageAndFreq.mat'
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/ALLMATS/DBS906_TurnOn_Day2_Sess2_20150828_032515_CurrentSweep.mat',
                        'lfp':''
                        }
                },
        '907':{
                'OnT':{
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/ALLMATS/DBS907_TurnOn_Day1_onTARGET_20151216_105913.mat',
                        'lfp':'/home/virati/MDD_Data/BR/907/Session_2015_12_16_Wednesday/DBS907_2015_12_16_12_09_04__MR_0.txt',
                        'epochs':{'Bilat':(640,670),'PreBilat':(590,620)},
                        },
                'OffT':{
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/ALLMATS/DBS907_TurnOn_Day2_offTARGET_20151217_094245.mat',
                        'lfp':'/home/virati/MDD_Data/BR/907/Session_2015_12_17_Thursday/DBS907_2015_12_17_10_53_08__MR_0.txt',
                        'epochs':{'Bilat':(625,655),'PreBilat':(560,590)},
                        },
                'Volt':{
                        #'fname':'/home/virati/MDD_Data/hdEEG/Continuous/ALLMATS/DBS907_TurnOn_Day2_Voltage_20151217_102952.mat'
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/ALLMATS/DBS907_TurnOn_Day3_Current_20151218_092443.mat',
                        'lfp':''
                        }
                },
        '908':{
                'OnT':{
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/Targeting/B04/DBS908/DBS908_TurnOn_Day1_onTARGET_20160210_125231.mat',
                        'lfp':'/home/virati/MDD_Data/BR/908/Session_2016_02_10_Wednesday/DBS908_2016_02_10_13_03_10__MR_0.txt',
                        'epochs':{'Bilat':(611,641),'PreBilat':(551,581)}},
                'OffT':{
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/Targeting/B04/DBS908/DBS908_TurnOn_Day2_offTARGET_20160211_123540.mat',
                        'lfp':'/home/virati/MDD_Data/BR/908/Session_2016_02_11_Thursday/DBS908_2016_02_11_12_34_21__MR_0.txt',
                        'epochs':{'Bilat':(611,641),'PreBilat':(551,581)},
                        },
                'Volt':{'fname':'',
                        'lfp':''
                        }
                },
        '910':{
                'OnT':{
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/ALLMATS/DBS910_TurnOn_OnTarget_20180530_022545.mat',
                        'lfp':''
                        
                        },
                'OffT':{
                        'fname':'/home/virati/MDD_Data/hdEEG/Continuous/ALLMATS/DBS910_TurnOn_OffTarget_TO_20180530_014051.mat',
                        'lfp':''},
                'Volt':{}
                }
            }

class EEG_check:
    def __init__(self,pt='908',condit='OnT',ds_fact=1,fs=500,spotcheck=False):
        pass

class responseLFP:
    def __init__(self):
        pass

class streamLFP:
    def __init__(self,pt='908',condit='OnT',ds_fact=1,spotcheck=False):
        self.donfft=2**10
        
        self.pt = pt
        self.condit = condit
        
        try: container = dbo.load_BR_dict(Targeting['All'][pt][condit]['lfp'],sec_offset=0)
        except: pdb.set_trace()
        fs = 422
        
        
        rec_length = container['Left'].shape
        
        self.tvect = np.linspace(0,rec_length[0]* 1/fs,rec_length[0])
        self.data_dict = container
        self.Fs = fs
        self.gen_epochs()
        
    def gen_epochs(self):
        print('Generating Epochs...')
        #Go in and generate the epochs associated with this recording
        self.epochs = Targeting['All'][self.pt][self.condit]['epochs']
        self.epochs['ALL'] = (0,-1)
        
    '''
    This should return the timeseries associated with a particular segment
    '''
    def time_series(self,epoch_name='All',full_stim=False):

        #Find the indices we need
        # Do adjustments here if you want LARGER or SMALLER epochs
        if full_stim:
            rec_idxs = np.where(np.logical_and(self.tvect < self.epochs[epoch_name][1]+140,self.tvect > self.epochs[epoch_name][0]-5))
        else:
            rec_idxs = np.where(np.logical_and(self.tvect < self.epochs[epoch_name][1],self.tvect > self.epochs[epoch_name][0]))

        return {key:self.data_dict[key][rec_idxs] for key in self.data_dict}
    
    def tf_transform(self,epoch_name):
        ts_dict = self.time_series(epoch_name)
        
        SG = dbo.gen_SG(ts_dict,Fs=422,overlap=False)
        
        return SG
    '''
    Plots and transforms below
    '''
    def osc_plot(self,epoch_name='All'):
        Osc_state = self.osc_transform(epoch_name)
        
    def osc_transform(self,epoch_name):
        #psds = dbo.gen_psd()
        firstpsds = dbo.gen_SG(self.time_series(epoch_name = epoch_name),Fs=self.Fs,nfft=self.donfft,overlap=False)
        fvect = np.linspace(0,self.Fs/2,firstpsds['Left']['F'].shape[0])
        psds = {chann:firstpsds[chann]['SG'] for chann in firstpsds.keys()}
        
        out_vect = nestdict()
        #need to manually go into segments
        for seg in range(psds['Left'].shape[1]):
            try:
                calcpsds = {chann:psds[chann][:,seg] for chann in ['Left','Right']}
            except: pdb.set_trace()
            out_vect[seg] = dbo.calc_feats(calcpsds,fvect,dofeats=['Delta','Theta','Alpha','Beta*','Gamma1'],modality='lfp')
        try:
            ret_vect = {chann: np.array([out_vect[seg][0][:,cc] for seg in range(psds[chann].shape[1])]) for cc,chann in enumerate(['Left','Right'])}
        except:
            pdb.set_trace()
        return ret_vect
    
    def plot_tf(self,epoch_name='All'):
        TF_dict = self.tf_transform(epoch_name) 
        
        dbo.plot_TF(TF_dict)
        
    def plot_f(self,epoch_name='All'):
        TF_dict = self.tf_transform(epoch_name)
        
        dbo.plot_F_fromTF(TF_dict)
        
        
        
class streamEEG:
    def __init__(self,pt='908',condit='OnT',ds_fact=1,spotcheck=False,reref_class='local',full_experiment=True):
        #self.data_dict = {ev:{condit:[] for condit in do_condits} for ev in do_pts}
        
        self.donfft = 2**10
        
        
        data_dict = defaultdict(dict)
        container = scio.loadmat(Targeting['All'][pt][condit]['fname'])
        dkey = [key for key in container.keys() if key[0:3] == 'DBS'][-1]
        fs = container['EEGSamplingRate']
               
                   
        #data_dict = np.zeros((257,6*60*fs))
        #THIS IS FINE SINCE it's like a highpass with a DCish cutoff
        #10 * 60 * fs:18*60*fs
        snippet = True
        
        #THIS CAPTURES THE ENTIRE EXPERIMENT INCLUDING UNILATERAL STIMULATION
        start_time=0
        
        if snippet:
            #906
            if pt == '906' and condit == 'OnT':
                #tint = (np.array([238,1090]) * self.fs).astype(np.int)
                start_time = 4
            elif pt == '906' and condit =='OffT':
                start_time = 3
            elif pt == '908' and condit == 'OnT':
                #tint = (np.array([1000,1800]) * self.fs).astype(np.int)
                start_time = 16
            elif pt == '908' and condit == 'OffT':
                start_time = 0
            elif pt == '907':
                start_time = 0
            elif pt == '901':
                start_time = 0
            elif pt == '910':
                start_time = 5
            
            if full_experiment:
                tlim = np.array((start_time,start_time + 16)) * 60 #in seconds
            else:
                tlim = np.array((start_time+8,start_time + 16)) * 60 #in seconds
                        
        tint = (tlim*fs).astype(np.int)[0]
        
        data_matr = sig.detrend(sig.decimate(container[dkey][:,tint[0]:tint[1]],ds_fact,zero_phase=True))
        #data_dict = data_dict - np.mean(data_dict,0)
        
        self.fs = container['EEGSamplingRate'][0] / ds_fact
        del(container)
                
        #self.data_dict[pt][condit] = data_dict
        self.data_matr = data_matr
        
        # Make a random timeseries for 256.... for some reason...?
        self.data_matr[256,:] = np.random.normal(size=self.data_matr[256,:].shape)
        
        
        
        # Do local re-referencing and set the datamatrix to the re-referenced data
        if reref_class:
            print('Doing Local Re-referencing...')
            self.data_matr = self.re_ref(scheme=reref_class)
        
        #can add LFP and interpolate HERE?!
        
        
        
        #Proceed with standard code for EEG that should not break with LFP addition
        self.tvect = np.linspace(tlim[0],tlim[1],data_matr.shape[1])
        self.fvect = np.linspace(0,self.fs/2,self.donfft/2+1)
                

        self.re_ref_data = nestdict()
        self.pt = pt
        self.condit = condit
        
    
        
    def re_ref(self,scheme='local'):
        #do a very simple lowpass filter at 1Hz
        hpf_cutoff=1/(self.fs/2)
        bc,ac = sig.butter(3,hpf_cutoff,btype='highpass',output='ba')
        
        if scheme == 'local':
            dist_matr = neigh_mont.return_cap_L(dth=3)
            
            dataref = self.data_matr
            post_ref = neigh_mont.reref_data(dataref,dist_matr)     
        elif scheme == 'avg':
            post_ref = self.data_matr - np.mean(self.data_matr,axis=0)
        elif scheme == 'none':
            post_ref = self.data_matr
            
        
        #self.re_ref_data = post_ref
        return post_ref
                
    def Osc_state(self):
        pass
    
    def make_epochs(self):
        pass
    
    def median_state(self,intv=(0,1),do_plot=False,use_maya=False,band='Alpha'):
        #Intervals NEED TO BE IN THE SEGMENT NUMBER
        
        band_i = dbo.feat_order.index(band)
        
        medians = np.median(self.osc_matr[:,intv[0]:intv[1],:],axis=1)
        
        #pdb.set_trace()
        if do_plot:
            if use_maya:
                EEG_Viz.maya_band_display(medians[:,band_i])
            else:
                EEG_Viz.plot_3d_scalp(medians[:,band_i],plt.figure(),label='Volt Mean Response ' + band + ' | ',unwrap=True,scale=100,clims=(-1,1),alpha=0.3,marker_scale=5)
    
                plt.suptitle(self.pt)
                
        
        return medians
    
    def median_response(self,intv=(0,1),do_plot=False,use_maya=False,band='Alpha'):
        #Intervals NEED TO BE IN THE SEGMENT NUMBER
        
        band_i = dbo.feat_order.index(band)
        
        medians = np.median(self.osc_matr[:,intv[0]:intv[1],:],axis=1) - self.baseline_state
        
        #pdb.set_trace()
        if do_plot:
            if use_maya:
                EEG_Viz.maya_band_display(medians[:,band_i])
            else:
                EEG_Viz.plot_3d_scalp(medians[:,band_i],plt.figure(),label='Volt Mean Response ' + band + ' | ',unwrap=True,scale=100,clims=(-5,5),alpha=0.5,marker_scale=5)
    
                plt.suptitle(self.pt)
                
        
        return medians
        
    def seg_PSDs(self):
        tvect = self.tvect
        max_idx = tvect.shape[0]
        int_len = 6*int(self.fs)
        
        
        idxs = range(0,max_idx,int_len)
        idxs = idxs[:-1]
        num_segs = len(idxs)
        
        self.psd_matr = np.zeros((257,num_segs,int(self.donfft/2)+1))
        self.osc_matr = np.zeros((257,num_segs,len(dbo.feat_order)))
        self.stim_feat = np.zeros((num_segs))
        
        for ii,idx in enumerate(idxs):
            #print('Transforming segment ' + str(ii) + ' at ' + str(idx))
            seg_dict = {ch:self.data_matr[ch,idx:idx+int_len].squeeze().reshape(-1,1) for ch in range(257)}
            #generate the psd
            psd_vect = dbo.gen_psd(seg_dict,Fs=self.fs,nfft=self.donfft)
            self.psd_matr[:,ii,:] = np.array([psd_vect[ch] for ch in range(257)])
            #self.stim_feat[ii] = dbo.calc_feats(self.psd_matr[:,ii,:],self.fvect,dofeats=['Stim'])[0]
            
            #subtract out the polynom
            #pdb.set_trace()
            postpoly = dbo.poly_subtrEEG(psd_vect,self.fvect.squeeze())[0]
            
            out_vect = dbo.calc_feats(postpoly,self.fvect,dofeats=['Delta','Theta','Alpha','Beta*','Gamma1','Stim'])[0]
            
            self.osc_matr[:,ii,:] = out_vect[0:5,:].T
            
            self.stim_feat[ii] = out_vect[-1,0]
            
            
        #seg_starts = tvect[0::2*self.fs]
    def load_classifier(self,ctype,train_type='cleaned'):
        if train_type == 'cleaned':
            self.clf = pickle.load(open('/home/virati/SVMModel_' + ctype,'rb'))
        elif train_type == 'stream':
            self.clf = pickle.load(open('/home/virati/Stream_SVMModel_' + ctype,'rb'))
   
    def calc_baseline(self,intv=(20,40)):
        #Which segments are with stim off?
        baseline_state = self.median_state(intv=intv) #large amount of time for us to average.
        
        self.baseline_state = baseline_state
        return 1

    def plot_segment_labels(self):
        plt.figure()
        elements = len(self.true_labels)
        
        plt.stem(self.true_labels)
        plt.xlabel('Segment number')
    
    def label_segments(self,baseline_calibration = True):
        #go to every stim_feat segment WITHOUT stimulation and average them together. This is like a calibration
        
        no_stim_segs = self.stim_feat < 10
        stim_segs = np.logical_not( no_stim_segs)
        
        self.label_time = np.zeros((self.osc_matr.shape[1]))
        self.stim_matr = np.zeros_like(self.osc_matr)
        self.true_labels = np.zeros((self.osc_matr.shape[1])) # setup our labels
        
        # Find the median of the epochs without stimulation along the axis of epochs
        self.no_stim_median = np.median(self.osc_matr[:,no_stim_segs,:],axis=1)
        
        #Go through each segment and subtract out the median of the stim
        for ss in range(self.osc_matr.shape[1]):
            self.label_time[ss] = ss

        
        if self.condit == 'OnT':
            label_val = 2
        elif self.condit == 'OffT':
            label_val = 1
        elif self.condit == 'Volt':
            label_val = 2
        
        # Transform from our labels to the integer labels
        self.label_val = label_val
        
        # label each segment with the condition of its time
        self.true_labels[stim_segs] = label_val
        
        
    def gen_test_matrix(self):
        #We copy the STIM matrix here, not the osc matrix.... TODO
        test_matr = np.copy(self.stim_matr)
        
        #Swap the 0 and 1's axes, I think to give us...?
        test_matr = np.swapaxes(test_matr,0,1)
        
        # Reshape to flatten our features TODO CHECK THE ORDER OF THIS AND MAKE SURE WE'RE KOSHER
        test_matr = test_matr.reshape(-1,257*5,order='F')
        
        # Directly return to us our matrix of interest
        return test_matr
    
    def classify_segs(self,ctype='l2',train_type='cleaned'):
        self.load_classifier(ctype,train_type)
        
        test_matr = self.gen_test_matrix()
        
        
        self.pred_labels = self.clf.predict(test_matr)
        
        labmap = {'OnTON':2,'OffTON':1,'OFF':0}
        
        
        
        pred_nums = np.array([labmap[label] for label in self.pred_labels])
        #pred_nums = sig.medfilt(pred_nums.astype(np.float64),5)
        
        print('Accuracy: ' + str(sum(pred_nums == self.true_labels)/len(pred_nums)))
        #print(stats.mode(pred_nums))
        
        plt.figure()
        #plt.plot(self.pred_labels)
        plt.plot(pred_nums,label='Predicted',linewidth=3)
        plt.plot(self.true_labels,label='True',linewidth=5,alpha=0.6)
        plt.legend()
        
        
        #What percentage of the time where it's "ONT" is 
        loc_true = self.true_labels == self.label_val
        print('Prob Measure OnT | OnT Stim - True Positive')
        print(sum(pred_nums[loc_true] == self.label_val) / sum(loc_true))
        print('Prob Measure OnT | NOT OnT Stim - False Positive')
        print(sum(pred_nums[np.logical_not(loc_true)] == self.label_val) / sum(np.logical_not(loc_true)))
    
        #BAYES FLIP
        loc_pos = pred_nums == self.label_val
        print('Prob OnT | predicted On Target - PPV')
        print(sum(self.true_labels[loc_pos] == self.label_val)/sum(loc_pos))
        print('Prob OnT | predicted NOT OnTarget - NPV')
        print(sum(self.true_labels[np.logical_not(loc_pos)] == self.label_val)/sum(np.logical_not(loc_pos)))
    
        return (pred_nums,self.true_labels)
    
    def plot_TF(self,chann):
        in_x = self.data_matr[chann,:]
        
        nperseg = 2**10
        noverlap = 512
        
        freq,time,specg = sig.spectrogram(in_x,nperseg=nperseg,noverlap=noverlap,window=sig.get_window('blackmanharris',nperseg),fs=self.fs)
        
        
        #self.F = F
        #self.T = T
        #self.SG = SG
        
        
        plt.figure()
        plt.subplot(211)
        plt.pcolormesh(time,freq.squeeze(),10*np.log10(specg))
        plt.colorbar()
            
    