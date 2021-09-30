#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 21:02:09 2018

@author: virati
Main Class for Processed dEEG Data

"""

import DBSpace as dbo
from DBSpace import simple_pca

from collections import defaultdict
import mne
from scipy.io import loadmat
import ipdb
import numpy as np

import scipy.signal as sig

import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from sklearn.utils import resample
plt.close('all')

import random
from DBSpace.visualizations import EEG_Viz
#import DBSpace.visualizations.EEG_Viz.EEG_Viz.plot_3d_scalp as EEG_Viz.plot_3d_scalp

#import EEG_Viz.EEG_Viz.plot_3d_scalp as EEG_Viz.plot_3d_scalp
#from EEG_Viz import EEG_Viz.plot_3d_scalp

import seaborn as sns

sns.set_context('paper')
sns.set(font_scale=4)
sns.set_style('white')


from DBSpace import nestdict

from statsmodels import robust

from sklearn import mixture
from sklearn.decomposition import PCA, FastICA
from sklearn import svm
import sklearn
from sklearn.metrics import confusion_matrix, roc_curve, auc, roc_auc_score
from sklearn.model_selection import learning_curve, StratifiedKFold

from sklearn.decomposition import FactorAnalysis, PCA
#import tensortools as tt
#from tensortools.operations import unfold as tt_unfold, khatri_rao
import tensorly as tl
from tensorly import unfold as tl_unfold
from tensorly.decomposition import parafac, tucker

# import some useful functions (they are available in utils.py)
#from utils import *

import pickle

import sys
sys.path.append('/home/virati/Dropbox/projects/libs/robust-pca/')
import r_pca
import pdb
#%%

TargetingEXP = defaultdict(dict)
#TargetingEXP['conservative'] = {'905':0,'906':0,'907':0,'908':0}
TargetingEXP['conservative'] = {
        '905':{'OnT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/conservative/DBS905_B4_OnTarget_HP_LP_seg_mff_cln_ref_con.mat','OffT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/conservative/DBS905_B4_OffTar_HP_LP_seg_mff_cln_ref_con.mat'},
        '906':{
                'OnT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/conservative/DBS906_TO_onTAR_MU_HP_LP_seg_mff_cln_ref_1.mat',
                'OffT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/conservative/DBS906_TO_offTAR_bcr_LP_HP_seg_bcr_ref.mat'
                },
        '907':{
                'OnT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/conservative/DBS907_TO_onTAR_MU_seg_mff_cln_ref.mat',
                'OffT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/conservative/DBS907_TO_offTAR_MU_seg_mff_cln_ref.mat'
                },
        '908':{
                'OnT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/conservative/DBS908_TO_onTAR_bcr_LP_seg_mff_cln_bcr_ref.mat',
                'OffT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/conservative/DBS908_TO_offTAR_bcr_MU_seg_mff_cln_ref_1.mat'
                }
        }
TargetingEXP['liberal'] = {
        '905':{'OnT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/liberal/DBS905_B4_OnTarget_HP_LP_seg_mff_cln_ref_lib.mat','OffT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/liberal/DBS905_B4_OffTar_HP_LP_seg_mff_cln_ref_lib.mat'},
        '906':{
                'OnT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/liberal/DBS906_TO_onTAR_MU_HP_LP_seg_mff_cln_ref.mat',
                'OffT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/liberal/DBS906_TO_offTAR_LP_seg_mff_cln_ref_1.mat'
                },
        '907':{
                'OnT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/liberal/DBS907_TO_onTAR_MU_seg_mff_cln_2ref.mat',
                'OffT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/liberal/DBS907_TO_offTAR_MU_seg_mff_cln_2ref.mat'
                },
        '908':{
                'OnT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/liberal/DBS908_TO_onTAR_bcr_LP_seg_mff_cln_ref.mat',
                'OffT':'/home/virati/MDD_Data/hdEEG/Segmented/Targeting_B4/liberal/DBS908_TO_offTAR_bcr_MU_seg_mff_cln_ref.mat'
                }
        }

keys_oi = {'OnT':["Off_3","BONT"],'OffT':["Off_3","BOFT"]}

class proc_dEEG:
    def __init__(self,pts,procsteps='liberal',condits=['OnT','OffT'],pretty_mode=False,polyfix=0):

        self.chann_dim = 257
        self.ch_order_list = range(self.chann_dim)
        self.procsteps = procsteps
        
        self.do_pts = pts
        self.condits = condits
        
        self.polyorder = polyfix
        self.pretty = pretty_mode
        
        self.label_map = {'OnT':1,'OffT':0}
        
        # Load in the data
        self.ts_data = self.load_data(pts)
        
        self.eeg_locs = mne.channels.read_montage('/home/virati/Dropbox/GSN-HydroCel-257.sfp')
        
        self.gen_output_variables()
    
    '''Setup all the output variables we need'''
    def gen_output_variables(self):
        # CHECK IF we're still using ANY of these
        
        #sloppy containers for the outputs of our analyses        
        self.psd_trans = {pt:{condit:{epoch:[] for epoch in keys_oi} for condit in self.condits} for pt in self.do_pts}
        self.PSD_diff = {pt:{condit:[] for condit in self.condits} for pt in self.do_pts}
        self.PSD_var = {pt:{condit:[] for condit in self.condits} for pt in self.do_pts}
        
        self.Feat_trans = {pt:{condit:{epoch:[] for epoch in keys_oi} for condit in self.condits} for pt in self.do_pts}
        self.Feat_diff = {pt:{condit:[] for condit in self.condits} for pt in self.do_pts}
        self.Feat_var = {pt:{condit:[] for condit in self.condits} for pt in self.do_pts}
    
    '''Run the standard pipeline to prepare segments'''
    def standard_pipeline(self):
        print('Doing standard init pipeline')
        self.extract_feats(polyorder=self.polyorder)
        self.pool_patients() #pool all the DBS RESPONSE vectors
        #self.median_responses = self.median_response(pt=self.do_pts)
        self.median_responses = self.median_bootstrap_response(pt='POOL',bootstrap=100)['mean']
        
    '''Load in the MAT data for preprocessed EEG recordings'''
    def load_data(self,pts):
        ts_data = defaultdict(dict)
        for pt in pts:
            ts_data[pt] = defaultdict(dict)
            for condit in self.condits:
                ts_data[pt][condit] = defaultdict(dict)
                
                temp_data = loadmat(TargetingEXP[self.procsteps][pt][condit])
                
                for epoch in keys_oi[condit]:
                    ts_data[pt][condit][epoch] = temp_data[epoch]

        self.fs = temp_data['EEGSamplingRate'][0][0]
        self.donfft = 2**11
        self.fvect = np.linspace(0,np.round(self.fs/2).astype(np.int),np.round(self.donfft/2+1).astype(np.int))
        
        return ts_data
        
    '''Extract features from the EEG datasegments'''
    def extract_feats(self,polyorder=4):
        pts = self.do_pts
        
        psd_dict = nestdict()
        osc_dict = nestdict()
        
        for pt in pts:
            #feat_dict[pt] = defaultdict(dict)
            
            for condit in self.condits:
                #feat_dict[pt][condit] = defaultdict(dict)
                for epoch in keys_oi[condit]:
                    #find the mean for all segments
                    data_matr = self.ts_data[pt][condit][epoch] #should give us a big matrix with all the crap we care about
                    data_dict = {ch:data_matr[ch,:,:].squeeze() for ch in range(data_matr.shape[0])} #transpose is done to make it seg x time
                    
                    #TODO check if this is in the right units
                    seg_psds = dbo.gen_psd(data_dict,Fs=self.fs,nfft=self.donfft,polyord=polyorder)
                    
                    #gotta flatten the DICTIONARY, so have to do it carefully
                    PSD_matr = np.array([seg_psds[ch] for ch in self.ch_order_list])
                    
                    OSC_matr = np.zeros((seg_psds[0].shape[0],257,len(dbo.feat_order)))
                    #middle_osc = {chann:seg_psd for chann,seg_psd in seg_psds.items}
                    middle_osc = np.array([seg_psds[ch] for ch in range(257)])
                    
                    #have to go to each segment due to code
                    for ss in range(seg_psds[0].shape[0]):
                        try:
                            state_return = dbo.calc_feats(middle_osc[:,ss,:],self.fvect)[0].T
                            state_return[:,4] = 0 #we know gamma is nothing for this entire analysis
                            #pdb.set_trace()
                            OSC_matr[ss,:,:] = np.array([state_return[ch] for ch in range(257)])
                        except Exception as e:
                            print('CRAP')
                            print(e)
                            ipdb.set_trace()
                    
                    #find the variance for all segments
                    psd_dict[pt][condit][epoch] = PSD_matr
                    osc_dict[pt][condit][epoch] = OSC_matr
                    
                    #need to do OSCILLATIONS here
        
        #THIS IS THE PSDs RAW, not log transformed
        self.psd_dict = psd_dict
        #Below is the oscillatory power
        #we should go through each osc_dict element and zero out the gamma
        self.osc_dict = osc_dict
        
    '''Generate Distributions across all channels for each band'''
    def band_distrs(self,pt='POOL'):
        print('Plotting Distribution for Bands')
        
        meds = nestdict()
        mads = nestdict()
        
        marker=['o','s']
        color = ['b','g']
        plt.figure()
        ax2 = plt.subplot(111)
        for cc,condit in enumerate(['OnT','OffT']):
            allsegs = np.median(self.osc_bl_norm[pt][condit][:,:],axis=0).squeeze()
            #pdb.set_trace()
            parts = ax2.violinplot(allsegs,positions=np.arange(5)+(cc-0.5)/10)
            for pc in parts['bodies']:
                pc.set_facecolor(color[cc])
                pc.set_edgecolor(color[cc])
                #pc.set_linecolor(color[cc])
                                 
            #plt.ylim((-0.5,0.5))
        
        for bb in range(5):
            pass
            #rsres = stats.ranksums(meds['OnT'][:,bb],meds['OffT'][:,bb])
            #rsres = stats.wilcoxon(meds['OnT'][:,bb],meds['OffT'][:,bb])
            #rsres = stats.ttest_ind(10**(meds['OnT'][:,bb]/10),10**(meds['OffT'][:,bb]/10))
            #print(rsres)
        
        #plt.suptitle(condit)
    
    ''' This function sets the response vectors to the targets x patient'''
    def compute_response(self,do_pts=[],condits=['OnT','OffT']):
        if do_pts == []:
            do_pts = self.do_pts
            
        BL = {pt:{condit:[] for condit in condits} for pt in do_pts}
        response = {pt:{condit:[] for condit in condits} for pt in do_pts}
        
        for pt in do_pts:
            for condit in condits:
                #first, compute the median state during baseline
                try: BL[pt][condit] = np.median(self.osc_dict[pt][condit]['Off_3'],axis=0)
                except: pdb.set_trace()               
                #Now, go to each segment during stim and subtract the BL for that
                response[pt][condit] = self.osc_dict[pt][condit][keys_oi[condit][1]] - BL[pt][condit]
                
        self.targ_response = response
        
    def DEPRtrain_SVM(self,mask=False):
        #Bring in and flatten our stack
        SVM_stack = 1
    
            
    def response_stats(self,band='Alpha',plot=False):
        band_idx = dbo.feat_order.index(band)
        response_diff_stats = {pt:[] for pt in self.do_pts}
        
        ## First, check to see if per-channel h-testing rejects the null
        for pt in self.do_pts:
            for cc in range(256):
                response_diff_stats[pt].append(stats.mannwhitneyu(self.targ_response[pt]['OnT'][:,cc,band_idx],self.targ_response[pt]['OffT'][:,cc,band_idx])[1])
    
        self.response_diff_stats = response_diff_stats
        
        ## Now check variances\
        ONT_var = {pt:[] for pt in self.do_pts}
        OFFT_var = {pt:[] for pt in self.do_pts}
        pool_ONT = []
        pool_OFFT = []
        for pt in self.do_pts:
            for cc in range(256):
                ONT_var[pt].append(np.var(self.targ_response[pt]['OnT'][:,cc,band_idx]))
                OFFT_var[pt].append(np.var(self.targ_response[pt]['OffT'][:,cc,band_idx]))
                
            pool_ONT.append(self.targ_response[pt]['OnT'][:,:,band_idx])
            pool_OFFT.append(self.targ_response[pt]['OffT'][:,:,band_idx])
                
        # Now stack across all patients
        pool_ONT_var = np.var(np.concatenate(pool_ONT,axis=0),axis=0)
        pool_OFFT_var = np.var(np.concatenate(pool_OFFT,axis=0),axis=0)
        
        ch_response_sig = {pt:np.array(response_diff_stats[pt]) for pt in self.do_pts}
        aggr_resp_sig = np.array([(resp < 0.05/256).astype(np.int) for pt,resp in ch_response_sig.items()])
        union_sig = np.sum(aggr_resp_sig,axis=0) >= 2
        
        if plot:
            for pt in self.do_pts:
                if 0:
                    pass
                
                    # Look at each patient's ONT and OFFT VARIANCE
                    bins = np.linspace(0,40,100)
                    plt.figure()
                    plt.violinplot(ONT_var[pt])#,bins=bins)
                    print(np.median(ONT_var[pt]))
                    plt.violinplot(OFFT_var[pt])#,bins=bins)
                    print(np.median(OFFT_var[pt]))
                    
                #Stats for ONT vs OFFT within each patient
                plt.figure()
                plt.plot(response_diff_stats[pt])
                plt.hlines(0.05/256,0,256)
                n_sig = np.sum((ch_response_sig[pt] < 0.05/256).astype(np.int))
                plt.suptitle(pt + ' ' + str(n_sig))
                    
                
                    
            plt.figure()
            plt.plot(pool_ONT_var)
            plt.plot(pool_OFFT_var)
            print(np.median(pool_ONT_var))
            print(np.median(pool_OFFT_var))
            plt.suptitle('Pooled stats for ONT/OFFT consistency check')
    
    '''Generate the pooled/ensemble segment response matrices -> dict'''
    def pool_patients(self):
        print('Pooling Patient Observations')
        self.osc_bl_norm = {pt:{condit:self.osc_dict[pt][condit][keys_oi[condit][1]] - np.median(self.osc_dict[pt][condit][keys_oi[condit][0]],axis=0) for condit in self.condits} for pt in self.do_pts}
        self.osc_bl_norm['POOL'] = {condit:np.concatenate([self.osc_dict[pt][condit][keys_oi[condit][1]] - np.median(self.osc_dict[pt][condit][keys_oi[condit][0]],axis=0) for pt in self.do_pts]) for condit in self.condits}
   
        self.osc_stim = nestdict()
        self.osc_stim = {pt:{condit:10**(self.osc_dict[pt][condit][keys_oi[condit][1]]/10) for condit in self.condits} for pt in self.do_pts}
        self.osc_stim['POOL'] = {condit:np.concatenate([10**(self.osc_dict[pt][condit][keys_oi[condit][1]]/10) for pt in self.do_pts]) for condit in self.condits}
     
    def plot_psd(self,pt,condit,epoch):
        plt.figure()
        #reshape
        reshaped = self.psd_dict[pt][condit][epoch].reshape(-1,1025)
        print(reshaped.shape)
        plt.plot(np.linspace(0,1000/2,1025),np.mean(20*np.log10(reshaped.T),axis=1),alpha=0.9)
        plt.plot(np.linspace(0,1000/2,1025),20*np.log10(reshaped.T)[:,np.random.randint(0,256,size=(30,))],alpha=0.2)
        plt.xlim((0,300))
        
    #Median dimensionality reduction here; for now rPCA
    def distr_response(self,pt='POOL'):
        return {condit:self.osc_bl_norm[pt][condit] for condit in self.condits}
    
    '''Compute the median response using bootstrap'''
    def median_bootstrap_response(self,pt='POOL',mfunc=np.mean,bootstrap=100):
        print('Computing Bootstrap Median Response for ' + pt)

        bs_mean = []
        bs_var = []
        for ii in range(bootstrap):
            rnd_idxs = {condit:random.sample(range(self.osc_bl_norm[pt][condit].shape[0]),100) for condit in self.condits}
            bs_mean.append({condit:mfunc(self.osc_bl_norm[pt][condit][rnd_idxs[condit],:,:],axis=0) for condit in self.condits})
            #bs_var.append({condit:np.var(self.osc_bl_norm[pt][condit][rnd_idxs[condit],:,:],axis=0) for condit in self.condits})
        
        mean_of_means = {condit:np.mean([iteration[condit] for iteration in bs_mean],axis=0) for condit in self.condits}
        var_of_means = {condit:np.var([iteration[condit] for iteration in bs_mean],axis=0) for condit in self.condits}

        return {'mean':mean_of_means, 'var':var_of_means}
    
    def OBSmedian_response(self,pt='POOL',mfunc = np.median):
        print('Computing Median Response for ' + pt)
        print('Doing ' + str(mfunc))
        return {condit:mfunc(self.osc_bl_norm[pt][condit],axis=0) for condit in self.condits}

        
    #In this function, we stack ONT_Off3 and OFFT_Off3 together to DEFINE the null distribution
    def OBScombined_bl(self):
        self.combined_BL = nestdict()
        for pt in self.do_pts:
            ONT_BL = self.osc_dict[pt]['OnT'][keys_oi['OnT'][0]]
            OFFT_BL = self.osc_dict[pt]['OffT'][keys_oi['OffT'][0]]
            
            if pt != '905':
                self.combined_BL[pt] = np.concatenate((ONT_BL,OFFT_BL),axis=0)
            else:
                self.combined_BL[pt] = ONT_BL
    
    def OBScombined_bl_distr(self,band='Alpha'):
        band_idx = dbo.feat_order.index(band)
        
        for pt in self.do_pts:
            plt.figure()
            for ch in range(256):
                plt.violinplot(self.combined_BL[pt][:,ch,band_idx])
    
    # Compare ONTarget and OFFTarget distributions
    def OBSONTvsOFFT(self,band='Alpha',stim=0):
        band_idx = dbo.feat_order.index(band)
        
        for pt in self.do_pts:
            ch_stat = np.zeros((257,))
            ch_ont = []
            ch_offt = []
            for ch in range(256):
                #distribution for pre-stimulation period
                #pdb.set_trace()
                ONT_distr = []
                OFFT_distr = []
                for ii in range(10):
                    ont_rand_idx = random.sample(range(0,self.osc_dict[pt]['OnT'][keys_oi['OnT'][stim]].shape[0]),10)
                    offt_rand_idx = random.sample(range(0,self.osc_dict[pt]['OffT'][keys_oi['OffT'][stim]].shape[0]),10)

                    ONT_distr.append(np.mean(self.osc_dict[pt]['OnT'][keys_oi['OnT'][stim]][ont_rand_idx,ch,band_idx]))
                    OFFT_distr.append(np.mean(self.osc_dict[pt]['OffT'][keys_oi['OffT'][stim]][offt_rand_idx,ch,band_idx]))
                    
                #baseline_distr = self.osc_dict[pt][condit][keys_oi[condit][0]][0:20,ch,band_idx]#should be segments x bands
                #stim_distr = self.osc_dict[pt][condit][keys_oi[condit][1]][0:20,ch,band_idx]
                diff_stat = stats.ranksums(ONT_distr,OFFT_distr)
                #diff_stat = stats.f_oneway(baseline_distr,stim_distr)
                print(str(ch) + ':' + str(diff_stat))
                ch_stat[ch] = diff_stat[1]
                
                ch_ont.append(np.mean(ONT_distr))
                ch_offt.append(np.mean(OFFT_distr))
                
            
            plt.figure()
            plt.violinplot(ch_ont)
            plt.violinplot(ch_offt)
            plt.ylim((-10,10))
            plt.suptitle(pt + ' stim: ' + str(stim))
            
            plt.figure()
            plt.plot(ch_stat)
            plt.axhline(0.05/256,0,256)
            
            plt.suptitle(pt + ' stim: '+ str(stim))
            
            
    #Do per-channel, standard stats. Compare pre-stim to stim condition
    def per_chann_stats(self,condit='OnT',band='Alpha'):
        band_idx = dbo.feat_order.index(band)
        
        for pt in self.do_pts:
            ch_stat = np.zeros((257,))
            ch_bl_mean = []
            ch_stim_mean = []
            for ch in range(256):
                #distribution for pre-stimulation period
                #pdb.set_trace()
                baseline_distr = []
                stim_distr = []
                for ii in range(100):
                    bl_rand_idx = random.sample(range(0,self.osc_dict[pt][condit][keys_oi[condit][0]].shape[0]),10)
                    stim_rand_idx = random.sample(range(0,self.osc_dict[pt][condit][keys_oi[condit][1]].shape[0]),10)

                    baseline_distr.append(np.mean(self.osc_dict[pt][condit][keys_oi[condit][0]][bl_rand_idx,ch,band_idx]))
                    stim_distr.append(np.mean(self.osc_dict[pt][condit][keys_oi[condit][1]][stim_rand_idx,ch,band_idx]))
                    
                #baseline_distr = self.osc_dict[pt][condit][keys_oi[condit][0]][0:20,ch,band_idx]#should be segments x bands
                #stim_distr = self.osc_dict[pt][condit][keys_oi[condit][1]][0:20,ch,band_idx]
                diff_stat = stats.mannwhitneyu(baseline_distr,stim_distr)
                #diff_stat = stats.f_oneway(baseline_distr,stim_distr)
                print(str(ch) + ':' + str(diff_stat))
                ch_stat[ch] = diff_stat[1]
                
            
            
                #plt.violinplot(baseline_distr)
                #plt.violinplot(stim_distr)
                ch_bl_mean.append(np.mean(baseline_distr))
                ch_stim_mean.append(np.mean(stim_distr))
            
            plt.figure()
            plt.plot(ch_stat)
            plt.axhline(0.05/256,0,256)
            
            plt.figure()
            plt.violinplot(ch_bl_mean)
            plt.violinplot(ch_stim_mean)
    
    
    def topo_median_variability(self,pt='POOL',band='Alpha',do_condits=[],use_maya=False):
        band_i = dbo.feat_order.index(band)
       
        #medians = self.median_response(pt=pt)

        for condit in do_condits:
            response_dict = np.median(self.osc_bl_norm[pt][condit][:,:,:],axis=0).squeeze()   

            
            var_dict = robust.mad(self.osc_bl_norm[pt][condit][:,:,:],axis=0).squeeze()      
            var_mask = ((var_dict) > 2.5).astype(np.int)
            #The old scatterplot approach
            bins = np.linspace(0,5,50)
            plt.figure();plt.hist(var_dict[:,band_i],bins=bins)
            plt.ylim((0,70))
            if use_maya:
                EEG_Viz.maya_band_display(var_dict[:,band_i])
            else:
                #EEG_Viz.plot_3d_scalp(response_mask[:,band_i]*var_dict[:,band_i],plt.figure(),label=condit + ' Response Var ' + band + ' | ' + pt,unwrap=True,scale=100,clims=(0,4),alpha=0.3,marker_scale=10)
                EEG_Viz.plot_3d_scalp(var_mask[:,band_i],plt.figure(),label=condit + ' Response Var ' + band + ' | ' + pt,unwrap=True,scale=100,clims=(0,4),alpha=0.3,marker_scale=10,anno_top=False,binary_mask=True)
                plt.suptitle(pt)
    
    def topo_median_response(self,pt='POOL',band='Alpha',do_condits=[],use_maya=False,scale_w_mad=False,avg_func=np.median):
        band_i = dbo.feat_order.index(band)
        #medians = self.median_response(pt=pt)

        for condit in do_condits:
            response_dict = avg_func(self.osc_bl_norm[pt][condit][:,:,:],axis=0).squeeze()      
            #The old scatterplot approach
            if use_maya:
                EEG_Viz.maya_band_display(response_dict[:,band_i])
            else:
                if scale_w_mad:
                    mad_scale = robust.mad(self.osc_bl_norm[pt][condit][:,:,:],axis=0).squeeze() 
                    
                    EEG_Viz.plot_3d_scalp(response_dict[:,band_i],plt.figure(),label=condit + ' Mean Response ' + band + ' | ' + pt,unwrap=True,scale=100,clims=(-1,1),alpha=0.3,marker_scale=10*np.tanh(mad_scale[:,band_i]-5))
                    
                else:
                    EEG_Viz.plot_3d_scalp(response_dict[:,band_i],plt.figure(),label=condit + ' Mean Response ' + band + ' | ' + pt,unwrap=True,scale=100,clims=(-1,1),alpha=0.3,marker_scale=5)
                plt.suptitle(pt)
    
    
    def OBSOBSsupport_analysis(self,pt='POOL',condit='OnT',voltage='3',band='Alpha'):
        support_struct = pickle.load(open('/tmp/'+ pt + '_' + condit + '_' + voltage,'rb'))
        distr = self.distr_response(pt=pt)
        #medians = np.median(self.targ_response[pt][condit],axis=0)
        fig = plt.figure()
        #First, we'll plot what the medians actually are
        band_i = dbo.feat_order.index(band)
        EEG_Viz.plot_3d_scalp(np.median(distr['OnT'][:,:,band_i],axis=0),fig,label='OnT Mean Response ' + band,unwrap=True,scale=10)
        plt.suptitle(pt)
        
        band_i = dbo.feat_order.index(band)
        
        full_distr = distr['OnT'][:,:,band_i]# - np.mean(medians['OnT'][:,band_i]) #this zeros the means of the distribution
        
        primary_distr = full_distr[:,support_struct['primary'] == 1]
        secondary_distr = full_distr[:,support_struct['secondary'] == 1]
        
        for cc in range(257):
            p_val[cc] = stats.ks_2samp(primary_distr[:,cc],secondary_distr[:,cc])
        
    '''
    Support analysis involves looking at forward-modeled EEG changes for Primary and Secondary nodes built from tractography
    
    '''
    def support_analysis(self,support_struct,pt='POOL',condit='OnT',voltage='3',band='Alpha'):
        #support_struct = pickle.load(open('/tmp/'+ pt + '_' + condit + '_' + voltage,'rb'))
        if band == 'rP0':    
            medians = self.dyn_L.swapaxes(0,1) #if we want to use the 0th component of the dyn_rPCA eigenvector
            band_i = 0
        else:
            medians = self.OBSmedian_response(pt=pt)[condit] #if we want to use the standard median Alpha change
            band_i = dbo.feat_order.index(band)
            
        #medians = np.median(self.targ_response[pt][condit],axis=0)
        fig = plt.figure()
        #First, we'll plot what the medians actually are
        
        EEG_Viz.plot_3d_scalp(medians[:,band_i],fig,label=condit + ' Mean Response ' + band,unwrap=True,scale=10)
        plt.suptitle(pt)
        
        full_distr = medians[:,band_i]# - np.mean(medians[:,band_i]) #this zeros the means of the distribution
        
        primary_distr = full_distr[support_struct['primary'] == 1]
        #now we'll circle where the primary nodes are
        
        print(np.sum((support_struct['primary'] == 1).astype(np.int)))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        EEG_Viz.plot_3d_scalp(support_struct['primary'],ax,scale=10,alpha=0.5,unwrap=True)
        plt.title('Primary Channels')
        
        secondary_distr = full_distr[support_struct['secondary'] == 1]
        print(np.sum((support_struct['secondary'] == 1).astype(np.int)))
        fig = plt.figure()
        EEG_Viz.plot_3d_scalp(support_struct['secondary'],fig,scale=10,alpha=0.5,unwrap=True)
        plt.title('Secondary Channels')
        
        labels = []
        def add_label(violin, label):
            color = violin["bodies"][0].get_facecolor().flatten()
            labels.append((mpatches.Patch(color=color), label))
            
        plt.figure()
        bins = np.linspace(-2,2,20)
        #plt.hist(primary_distr,bins=bins,alpha=0.5,label='Primary')
        print('Primary mean: ' + str(np.median(primary_distr)))
        add_label(plt.violinplot(primary_distr),'Primary Nodes')
        #pdb.set_trace()
        
        #plt.hist(secondary_distr,bins=bins,alpha=0.5,label='Secondary')
        print('Secondary mean: ' + str(np.median(secondary_distr)))
        add_label(plt.violinplot(secondary_distr),'Secondary Nodes')
        plt.legend(*zip(*labels),loc=2)
        
        print(stats.ks_2samp(primary_distr,secondary_distr))
        
        #plt.hist(full_distr,bins=bins,alpha=0.5,label='FULL')
        #plt.legend(['Primary','','','Secondary'])
        plt.title(pt + ' ' + condit + ' ' + band)
    
    '''I guess this is about developing a rPCA approach to *dynamic* response without oscillations?'''
    def OnT_ctrl_dyn(self,pt='POOL',condit='OnT',do_plot=False):
        source_label = 'Dyn PCA'
        
        response_stack = self.osc_bl_norm['POOL'][condit][:,:,2]
        # Focusing just on alpha
        #response_stack = np.dot(response_stack.T,response_stack)
        
        #pdb.set_trace()
        rpca = r_pca.R_pca(response_stack)
        L,S = rpca.fit()
        
        svm_pca = PCA()
        svm_pca.fit(L)
        
        #pdb.set_trace()
        svm_pca_coeffs = svm_pca.components_
        # ALL PLOTTING BELOW
        if do_plot:
            for comp in range(5):
                fig = plt.figure()
                EEG_Viz.plot_3d_scalp((L[comp,:]),fig,label='OnT Mean Response',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
                plt.title('rPCA Component ' + str(comp))
                
            
            plt.figure()
            plt.subplot(221)
            plt.plot(svm_pca.explained_variance_ratio_)
            plt.ylim((0,1))
            plt.subplot(222)
            plt.plot(np.array(svm_pca_coeffs)[0:4])
            plt.legend(['PC1','PC2','PC3','PC4','PC5'])
            plt.title('rPCA Components ' + source_label)
            
        self.dyn_pca = svm_pca
        self.dyn_L = L
    
    def OnT_ctrl_modes_segs_ICA(self,pt='POOL',do_plot=False):   

        seg_responses = self.osc_bl_norm[pt]['OnT'][:,:,0:4]
        source_label = 'Segment Responses'
        
        svm_ica_coeffs = []
        for ii in range(seg_responses.shape[0]):
            
            #pdb.set_trace()
            rpca = r_pca.R_pca(seg_responses[ii,:,:])
            L,S = rpca.fit()
            
            #L = seg_responses[ii,:,:]
            #S = 0
            svm_ica = FastICA(tol=0.1,max_iter=1000)
    
            svm_ica.fit(L)
            rotated_L = svm_ica.fit_transform(L)
            
            svm_ica_coeffs.append(svm_ica.components_)
        
        mode_model = {'L':L, 'S':S, 'Vectors':svm_ica_coeffs, 'Model':svm_ica, 'RotatedL':rotated_L}
        return mode_model
    
    
    def OnT_alpha_modes_segs(self,pt='POOL',data_source=[],do_plot=False,band='Alpha'):
        seg_responses = self.osc_bl_norm[pt]['OnT'][:,:,dbo.feat_order.index(band)].squeeze()
        source_label = 'Alpha Segmental Response'
        
        lr_pca_coeffs = []
        S_pca_coeffs = []
        rot_L = []
        rot_S = []
        
        #pdb.set_trace()
        rpca = r_pca.R_pca(seg_responses[:,:].T)
        L,S = rpca.fit()
        
        #L = seg_responses[ii,:,:]
        #S = 0
        lr_pca = PCA()

        lr_pca.fit(L)
        rotated_L = lr_pca.fit_transform(L)
        rot_L.append(rotated_L)
        
        lr_pca_coeffs.append(lr_pca.components_)
    
        S_pca = PCA()
        S_pca.fit(S)
        rotated_S = S_pca.fit_transform(S)
        rot_S.append(rotated_S)
        S_pca_coeffs.append(S_pca.components_)
        
        mode_model = {'L':L, 'S':S, 'SModel':S_pca, 'RotatedS':np.array(rot_S), 'SVectors':S_pca_coeffs,'Vectors':lr_pca_coeffs, 'Model':lr_pca, 'RotatedL':np.array(rot_L)}
        return mode_model

    def topo_OnT_alpha_ctrl(self,**kwargs):
        self.alpha_ctrl_model = self.OnT_alpha_modes_segs(**kwargs)
        #model = self.alpha_ctrl_model
        
        
    def plot_alpha_ctrl_L(self,top_comp=5):
        model = self.alpha_ctrl_model
        
        L = model['RotatedL']
        expl_var = model['Model'].explained_variance_ratio_
        coeffs = np.array(model['Vectors'])

        
        # ALL PLOTTING BELOW
        #Plot the topo for our low-rank component
        for comp in range(0,top_comp):
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(np.median(L,axis=0)[:,comp],fig,label='OnT Mean Response',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
            plt.title('rPCA Component ' + str(comp) + ' ' + str(expl_var[comp]))
        
        
        plt.figure()
        plt.subplot(221)
        plt.plot(expl_var)
        plt.ylim((0,1))
        plt.subplot(222)
        for ii in range(4): #this loops through our COMPONENTS to find the end
            plt.plot(np.mean(coeffs,axis=0)[ii,:].T,linewidth=5-ii,alpha=expl_var[ii])
        plt.ylim((-0.3,0.3))
        #plt.hlines(0,0,5)
        plt.legend(['PC0','PC1','PC2','PC3','PC4'])
        plt.title('rPCA Components')
        
    def plot_alpha_ctrl_S(self,top_comp=10):
        model = self.alpha_ctrl_model
        
        S = model['RotatedS']
        s_coeffs = np.array(model['SVectors'])
        s_expl_var = model['SModel'].explained_variance_ratio_
        
        for comp in range(0,top_comp):
            #Plot sparse component next
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(np.median(S,axis=0)[:,comp],fig,label='OnT Mean Response',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
            plt.title('Sparse Component ' + str(comp))
            
        plt.figure()
        plt.subplot(221)
        plt.plot(s_expl_var)
        plt.ylim((0,1))
        plt.subplot(222)
        
        for ii in range(4): #this loops through our COMPONENTS to find the end
            plt.plot(np.mean(s_coeffs,axis=0)[ii,:].T,linewidth=5-ii)#,alpha=expl_var[ii])
        plt.ylim((-0.3,0.3))
        plt.hlines(0,0,5)
        plt.legend(['PC1','PC2','PC3','PC4','PC5'])
        plt.title('rPCA Components')   
    
    #def OnT_ctrl_modes_segs(self,pt='POOL',data_source=[],do_plot=False):
    
    ''' We do naive rPCA on each segment and then average things together: this has been superceded by the tensor decomposition'''
    def OnT_ctrl_modes_segs(self,**kwargs):
        pt = kwargs['pt']
        do_plot = kwargs['do_plot']

        print('Using BL Norm Segments - RAW')
        seg_responses = self.osc_bl_norm[pt]['OnT'][:,:,0:4]
        source_label = 'Segment Responses'
        
        lr_pca_coeffs = []
        S_pca_coeffs = []
        rot_L = []
        rot_S = []
        for ii in range(seg_responses.shape[0]):
            
            #pdb.set_trace()
            rpca = r_pca.R_pca(seg_responses[ii,:,:])
            L,S = rpca.fit()
            
            #L = seg_responses[ii,:,:]
            #S = 0
            lr_pca = PCA()
    
            lr_pca.fit(L)
            rotated_L = lr_pca.fit_transform(L)
            rot_L.append(rotated_L)
            
            lr_pca_coeffs.append(lr_pca.components_)
        
            S_pca = PCA()
            S_pca.fit(S)
            rotated_S = S_pca.fit_transform(S)
            rot_S.append(rotated_S)
            S_pca_coeffs.append(S_pca.components_)
        
        mode_model = {'L':L, 'S':S, 'SModel':S_pca, 'RotatedS':np.array(rot_S), 'SVectors':S_pca_coeffs,'Vectors':lr_pca_coeffs, 'Model':lr_pca, 'RotatedL':np.array(rot_L)}
        return mode_model
    
    '''THIS IS OBSOLETE FOR TENSOR ANALYSES'''
    def topo_OnT_ctrl_segs(self,**kwargs):
        model = self.OnT_ctrl_modes_segs(**kwargs)
        
        L = model['RotatedL']
        expl_var = model['Model'].explained_variance_ratio_
        coeffs = np.array(model['Vectors'])
        S = model['RotatedS']
        s_coeffs = np.array(model['SVectors'])
        #s_expl_var = model['SModel'].explained_variance_ratio_
        
        # ALL PLOTTING BELOW
        #Plot the topo for our low-rank component
        for comp in range(1):
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(np.mean(L,axis=0)[:,comp],fig,label='OnT Mean Response',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
            #plt.title('rPCA Component ' + str(comp))
            plt.title('Median Low-Rank Component ' + str(comp))
        
        # Plot our Components
        plt.figure()
        plt.subplot(221)
        plt.plot(expl_var)
        plt.ylim((0,1))
        plt.subplot(222)
        for ii in range(4): #this loops through our COMPONENTS to find the end
            plt.plot(coeffs[ii,:].T,linewidth=5-ii,alpha=expl_var[ii])#,alpha=expl_var[ii])
        plt.ylim((-1,1))
        plt.hlines(0,0,3)
        plt.legend(['PC0','PC1','PC2','PC3','PC4'])
        plt.title('Low-Rank Components')
        
        plot_sparse = False
        if plot_sparse:
            for comp in range(4):
                #Plot sparse component next
                fig = plt.figure()
                EEG_Viz.plot_3d_scalp(np.median(S,axis=0)[:,comp],fig,label='OnT Mean Response',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
                plt.title('Sparse Component ' + str(comp))
                
            plt.figure()
            plt.subplot(221)
            plt.plot(s_expl_var)
            plt.ylim((0,1))
            plt.subplot(222)
            
            for ii in range(4): #this loops through our COMPONENTS to find the end
                plt.plot(np.mean(s_coeffs,axis=0)[ii,:].T,linewidth=5-ii)#,alpha=expl_var[ii])
            plt.ylim((-0.3,0.3))
            plt.hlines(0,0,3)
            plt.legend(['PC0','PC1','PC2','PC3','PC4'])
            plt.title('Sparse Components')
            
        if kwargs['plot_maya']:
            response_dict = np.median(L,axis=0)[:,comp].squeeze()      
            EEG_Viz.maya_band_display(response_dict)
    
    
    #Dimensionality reduction of ONTarget response; for now rPCA
    def topo_OnT_ctrl_tensor(self,**kwargs):
        pt = kwargs['pt']
        seg_responses = self.osc_bl_norm[pt]['OnT'][:,:,0:4]
        
        #factors = parafac(seg_responses,rank=4)
        print(seg_responses.shape)
        #core, factors = tucker(seg_responses)
        factors = parafac(seg_responses.swapaxes(0,1),rank=4) # factors gives us weight, factors
        #plt.plot(core[0])
        #print(len(factors))
        print(factors)
        for ii in range(4):
            fig = plt.figure()
            #EEG_Viz.plot_3d_scalp(seg_responses[20,:,2],fig,label='Raw Segment',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
            EEG_Viz.plot_3d_scalp(factors[1][0][:,ii],fig,label='Tensor Decomp ' + str(ii),unwrap=True,scale=100,alpha=0.3,marker_scale=5)
        print(factors[0])
        #print(core.shape)
        #print((factors))


    def OnT_ctrl_modes(self,pt='POOL',data_source=[],do_plot=False,plot_maya=True):
        
        print('Using BL Norm Segments - RAW')
        med_response = np.median(self.osc_bl_norm[pt]['OnT'],axis=0).squeeze()
        source_label = 'BL Normed Segments'
        
        svm_pca_coeffs = []
        rpca = r_pca.R_pca(med_response)
        L,S = rpca.fit()
        
        #L = med
        svm_pca = PCA()
        svm_pca.fit(L)
        rotated_L = svm_pca.fit_transform(L)
        
        svm_pca_coeffs = svm_pca.components_
        
        #plt.scatter(med_response[:,2],med_response[:,3])
  
        mode_model = {'L':L, 'S':S, 'SModel':[], 'RotatedS':[], 'SVectors':[], 'Vectors':svm_pca_coeffs, 'Model':svm_pca, 'RotatedL':rotated_L}
        return mode_model
    
    
    def topo_OnT_ctrl(self,**kwargs):
        model = self.OnT_ctrl_modes(**kwargs)
        
        L = model['RotatedL']
        expl_var = model['Model'].explained_variance_ratio_
        coeffs = np.array(model['Vectors'])
        #s_expl_var = model['SModel'].explained_variance_ratio_
        
        # ALL PLOTTING BELOW
        #Plot the topo for our low-rank component
        for comp in range(2):
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(L[:,comp],fig,label='OnT Mean Response',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
            #plt.title('rPCA Component ' + str(comp))
            plt.title('Median Low-Rank Component ' + str(comp))
        
        # Plot our Components
        plt.figure()
        plt.subplot(221)
        plt.plot(expl_var)
        plt.ylim((0,1))
        plt.subplot(222)
        for cc in range(4): #this loops through our COMPONENTS to find the end
            plt.plot(coeffs[cc,:],linewidth=5-cc,alpha=0.2)#,alpha=expl_var[ii])
        plt.ylim((-1,1))
        plt.hlines(0,0,3)
        plt.legend(['PC0','PC1','PC2','PC3','PC4'])
        plt.title('Low-Rank Components')
        
        plot_sparse = False
        if plot_sparse:
            for comp in range(4):
                #Plot sparse component next
                fig = plt.figure()
                EEG_Viz.plot_3d_scalp(np.median(S,axis=0)[:,comp],fig,label='OnT Mean Response',unwrap=True,scale=100,alpha=0.3,marker_scale=5)
                plt.title('Sparse Component ' + str(comp))
                
            plt.figure()
            plt.subplot(221)
            plt.plot(s_expl_var)
            plt.ylim((0,1))
            plt.subplot(222)
            
            for ii in range(4): #this loops through our COMPONENTS to find the end
                plt.plot(np.mean(s_coeffs,axis=0)[ii,:].T,linewidth=5-ii)#,alpha=expl_var[ii])
            plt.ylim((-0.3,0.3))
            plt.hlines(0,0,3)
            plt.legend(['PC0','PC1','PC2','PC3','PC4'])
            plt.title('Sparse Components')
            
        if kwargs['plot_maya']:
            #response_dict = np.median(L,axis=0)#[:,comp].squeeze()
            response = L[:,1].squeeze()
            EEG_Viz.maya_band_display(response)
            #EEG_Viz.plot_3d_scalp(response)
    
    def dict_all_obs(self,condits=['OnT']):
        full_stack = nestdict()
        label_map = self.label_map
        

        for condit in condits:
            full_stack[condit]['x'] = self.osc_bl_norm['POOL'][condit]
            full_stack[condit]['g'] = [label_map[condit] for seg in self.osc_bl_norm['POOL'][condit]]
            
        # List comprehend version
        #full_stack = {condit:}
        
        return full_stack
    
    def control_rotate(self,condits=['OnT','OffT']):
        #get our bases
        control_bases = self.control_bases
        control_modes = self.control_modes
        #Get our observations, ONT and OFFT alike
        obs_dict = self.dict_all_obs(condits=condits)
        
        trajectories = nestdict()
        rotated_stack = nestdict()
        for condit in condits:
            rotated_stack[condit] = np.dot(obs_dict[condit]['x'],control_bases)
            
            for ii in range(2):
                trajectories[condit][ii] = np.dot(rotated_stack[condit][:,:,ii].squeeze(),control_modes[:,ii].squeeze())
        
        # PLOTTING HERE
        plt.figure()
        #plt.subplot(121)
        for ii in [0,-1]:
            plt.scatter(trajectories['OnT'][0][ii],trajectories['OnT'][1][ii],cmap='jet',marker='s')
        plt.scatter(trajectories['OnT'][0],trajectories['OnT'][1],cmap='jet')
        #plt.plot(trajectories['OnT'][0],trajectories['OnT'][1],alpha=0.8)
        plt.xlim([-100,100])
        plt.ylim([-100,100])
        
        plt.figure()
        #plt.subplot(122)
        plt.scatter(trajectories['OnT'][0],trajectories['OnT'][1],c=np.linspace(0,1,trajectories['OnT'][0].shape[0]),cmap='jet',alpha=0.1)
        #plt.plot(trajectories['OnT'][0],trajectories['OnT'][1],alpha=0.1)

        for ii in [0,-1]:
            plt.scatter(trajectories['OffT'][0][ii],trajectories['OffT'][1][ii],cmap='jet',marker='s')
        plt.scatter(trajectories['OffT'][0],trajectories['OffT'][1],cmap='jet')
        #plt.plot(trajectories['OffT'][0],trajectories['OffT'][1],alpha=0.8)
        plt.xlim([-100,100])
        plt.ylim([-100,100])
        
    def OBStrain_simple(self):
        #Train our simple classifier that just finds the shortest distance
        self.signature = {'OnT':0,'OffT':0}
        self.signature['OnT'] = self.pop_osc_change['OnT'][dbo.feat_order.index('Alpha')] / np.linalg.norm(self.pop_osc_change['OnT'][dbo.feat_order.index('Alpha')])
        self.signature['OffT'] = self.pop_osc_change['OffT'][dbo.feat_order.index('Alpha')] / np.linalg.norm(self.pop_osc_change['OffT'][dbo.feat_order.index('Alpha')])
        
    
    def OBStest_simple(self):
        #go to our GMM stack and, for each segment, determine the distance to the two conditions
        
        #Set up our signature
        cort_sig = {'OnT':self.Seg_Med[0]['OnT'],'OffT':self.Seg_Med[0]['OffT']}
        
        #Now let's generate a stacked set of the SAME 
        
        for condit in self.condits:
            seg_stack = self.GMM_Osc_stack[condit]
            seg_num = seg_stack.shape[1]
            OnT_sim[condit] = [None] * seg_num
            OffT_sim[condit] = [None] * seg_num
            
            for seg in range(seg_num):
                net_vect = seg_stack[:,seg,dbo.feat_order.index('Alpha')].reshape(-1,1)
                net_vect = net_vect/np.linalg.norm(net_vect)
                
                OnT_sim[condit][seg] = np.arccos(np.dot(net_vect.T,self.signature['OnT'].reshape(-1,1)) / np.linalg.norm(net_vect))
                OffT_sim[condit][seg] = np.arccos(np.dot(net_vect.T,self.signature['OffT'].reshape(-1,1)) / np.linalg.norm(net_vect))
                
            OnT_sim[condit] = np.array(OnT_sim[condit])
            OffT_sim[condit] = np.array(OffT_sim[condit])
            
        return (OnT_sim, OffT_sim)
    

    
    def OBSinterval_stats(self,do_band='Alpha'):
        big_stack = self.osc_dict
        band_idx = dbo.feat_order.index(do_band)
        plag = nestdict()
        
        for pt in self.do_pts:
            for condit in self.condits:
                plag[pt][condit] = big_stack[pt][condit][keys_oi[condit][1]][:,:,band_idx] - np.median(big_stack[pt][condit]['Off_3'][:,:,band_idx],axis=0)
        
            sig_chann_list = []
        
            plt.figure()
            bins = np.linspace(-10,10,20)
            for ch in range(257):
                #Unfortunately, wilxocon won't work since we're not matched :(
                #wrst = stats.wilcoxon(plag[pt]['OnT'][:,ch],plag[pt]['OffT'][:,ch])
                #mwut = stats.mannwhitneyu(plag[pt]['OnT'][:,ch],plag[pt]['OffT'][:,ch]) #to compare OnT with OffT
                mwut = stats.wilcoxon(plag[pt]['OnT'][:,ch])
                #kstest = stats.kstest(plag[pt]['OnT'][:,ch],'norm')
                
                
                usestat = mwut
                if usestat[1] < 0.05/256:
                    sig_chann_list.append(ch)
                    print(str(ch) + ' ' + str(usestat))
                    plt.hist(plag[pt]['OnT'][:,ch],bins,alpha=1)
                    
            #now we want 3d plot of significant channels!
            sig_chann_list = np.array(sig_chann_list) 
            sig_stat_mask = np.zeros((257,))
            sig_stat_mask[sig_chann_list] = 1
            
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(sig_stat_mask,fig,animate=False,unwrap=True)
    
    #in this method, we're going to do per-channel statistics for each patient, channel, band

    def band_stats(self,do_band='Alpha'):
        self.pop_meds()
        
    def plot_band_stats(self,do_band='Alpha'):
        self.plot_meds(band=do_band,flatten=not self.pretty)
          
    

    def OBSsimple_stats(self):
        # We have a bit stack of the segments oscillatory powers
        
        ref_stack = self.big_stack_dict
        #Work with the Osc Dict data
        for condit in self.condits:
            ref_stack[condit]['Diff'] = defaultdict()#{key:[] for key in ([self.do_pts] + ['All'])}
            for epoch in ['OF','ON']:
                for pt in self.do_pts:
                    ref_stack[condit][epoch][pt + 'med'] = np.median(ref_stack[condit][epoch][pt],axis=0)
                    ref_stack[condit][epoch][pt + 'mad'] = robust.mad(ref_stack[condit][epoch][pt],axis=0)
                #stack all
                all_stack = [ref_stack[condit][epoch][pt] for pt in self.do_pts]
                
                ref_stack[condit][epoch]['MED'] = np.median(np.concatenate(all_stack,axis=0),axis=0)
                ref_stack[condit][epoch]['MAD'] = robust.mad(np.concatenate(all_stack,axis=0),axis=0)
            
            for pt in self.do_pts:
                ref_stack[condit]['Diff'][pt] = ref_stack[condit]['ON'][pt+'med'] - ref_stack[condit]['OF'][pt+'med']
            ref_stack[condit]['Diff']['All'] = ref_stack[condit]['ON']['MED'] - ref_stack[condit]['OF']['MED']
        
        all_stack = np.concatenate([np.concatenate([np.concatenate([ref_stack[condit][epoch][pt] for pt in self.do_pts],axis=0) for epoch in ['OF','ON']],axis=0) for condit in self.condits],axis=0)
        label_stack = [[[[condit+epoch for seg in ref_stack[condit][epoch][pt]] for pt in self.do_pts] for epoch in ['OF','ON']] for condit in self.condits]
        label_list = [item for sublist in label_stack for item in sublist]
        label_list = [item for sublist in label_list for item in sublist]
        label_list = [item for sublist in label_list for item in sublist]
        
        #go through each and if the last two characters are 'OF' -> 'OF' is label
        newlab = {'OffTOF':'OFF','OnTOF':'OFF','OnTON':'OnTON','OffTON':'OffTON'}
        label_list = [newlab[item] for item in label_list]
        
        self.SVM_stack = all_stack
        self.SVM_labels = np.array(label_list)
        
    
    def find_seg_covar(self):
        covar_matrix = nestdict()
        
        for condit in self.condits:
            seg_stack = sig.detrend(self.GMM_Osc_stack[condit],axis=1,type='constant')
            seg_num = seg_stack.shape[1]
            
            for bb,band in enumerate(dbo.feat_order):
                covar_matrix[condit][band] = []
                for seg in range(seg_num):
                    
                    net_vect = seg_stack[:,seg,bb].reshape(-1,1)
                    
                    cov_matr = np.dot(net_vect,net_vect.T)
                    covar_matrix[condit][band].append(cov_matr)
                    
                covar_matrix[condit][band] = np.array(covar_matrix[condit][band])
                    
        
        self.cov_segs = covar_matrix
                    
    def plot_seg_covar(self,band='Alpha'):
        for condit in self.condits:
            plt.figure()
            plt.subplot(211)
            plt.imshow(np.median(self.cov_segs[condit][band],axis=0),vmin=0,vmax=0.5)
            plt.colorbar()
            plt.title('Mean Covar')
            
            plt.subplot(212)
            plt.imshow(np.var(self.cov_segs[condit][band],axis=0),vmin=0,vmax=2)
            plt.colorbar()
            plt.title('Var Covar')
            
            plt.suptitle(condit)
    
    def OnT_v_OffT_MAD(self,band='Alpha'):
        Xdsgn = self.SVM_stack
        #THIS DOES MAD across the segments!!!!!!
        X_onT = Xdsgn[self.SVM_labels == 'OnTON',:,:].squeeze()
        X_offT = Xdsgn[self.SVM_labels == 'OffTON',:,:].squeeze()
        X_NONE = Xdsgn[self.SVM_labels == 'OFF',:,:].squeeze()
        
        OnT_MAD = robust.mad(X_onT,axis=0)
        OffT_MAD = robust.mad(X_offT,axis=0)
        NONE_MAD = robust.mad(X_NONE,axis=0)
        
        print('OnT segs: ' + str(X_onT.shape[0]))
        print('OffT segs: ' + str(X_offT.shape[0]))
        print('OFF segs: ' + str(X_NONE.shape[0]))
        
        self.Var_Meas = {'OnT':{'Med':np.median(X_onT,axis=0),'MAD':OnT_MAD},'OffT':{'Med':np.median(X_offT,axis=0),'MAD':OffT_MAD},'OFF':{'Med':np.median(X_NONE,axis=0),'MAD':NONE_MAD}}
        
    def plot_pca_decomp(self,pca_condit='OnT',approach = 'rpca'):
        self.pca_decomp(direction='channels',condit=pca_condit,bl_correct=True,pca_type=approach,plot_distr=True)
        
        plt.figure();
        plt.subplot(221)
        plt.imshow(self.PCA_d.components_,cmap=plt.cm.jet,vmax=1,vmin=-1)
        plt.colorbar()
        plt.subplot(222)
        plt.plot(self.PCA_d.components_)
        plt.ylim((-1,1))
        plt.legend(['PC0','PC1','PC2','PC3','PC4'])
        plt.xticks(np.arange(0,5),['Delta','Theta','Alpha','Beta','Gamma1'])
        plt.subplot(223)
        
        plt.plot(self.PCA_d.explained_variance_ratio_)
        plt.ylim((0,1))
        
        for cc in range(2):
            
            #plot the boring views first
            plt.figure()
            plt.subplot(211)
            plt.plot(self.PCA_x[:,cc])
            plt.subplot(212)
            plt.hist(self.PCA_x[:,cc],bins=np.linspace(-1,1,50))
            #find the top mode
            chann_high = np.where(np.abs(self.PCA_x[:,cc]) > 0.7)
            print(chann_high)
            
            #Plot the 3d scalp distribution
            fig=plt.figure()
            EEG_Viz.plot_3d_scalp(self.PCA_x[:,cc],fig,animate=False,unwrap=True,highlight=chann_high)
            plt.title('Plotting component ' + str(cc))
            plt.suptitle(approach + ' rotated results for ' + pca_condit)

            
            
    def DEPRpca_decomp(self,direction='channels',band='Alpha',condit='OnT',bl_correct=False,pca_type='pca',plot_distr=False):
        print('Doing PCA on the SVM Oscillatory Stack')
        #check to see if we have what variables we need
        Xdsgn = self.SVM_stack
        lbls = self.SVM_labels
        
        lblsdo = lbls == condit + 'ON'
        
        Xdo = Xdsgn[lblsdo,:,:]
        
        #what do we want to do with this now?
        #spatiotemporal PCA
        Xdo = np.median(Xdo,axis=0) 
        
        #BL_correct here is NOT patient specific bl correction
        if bl_correct:
            print('Correcting with baseline (OFF EEG)')
            #find the stim Off timepoints to subtract
            X_bl = np.median(Xdsgn[lbls=='OFF',:,:],axis=0)
            
            Xdo = Xdo - X_bl
            
            
        PCAdsgn = sig.detrend(Xdo,axis=0,type='constant')
        PCAdsgn = sig.detrend(PCAdsgn,axis=1,type='constant')
        
        bins = np.linspace(-3,3,50)
        
        # If we want to do PCA here
        if pca_type == 'pca':
            pca = PCA()
            pca.fit(PCAdsgn)
            
            if plot_distr:
                plt.figure()
                for bb in range(5):
                    plt.hist(PCAdsgn[:,bb],bins=bins,alpha=0.2)
                plt.suptitle('PCA inputs')
    
            self.PCA_d = pca
            self.PCA_inX = Xdo
            
            PCA_X = pca.fit_transform(PCAdsgn)
            self.PCA_x = PCA_X
        elif pca_type == 'rpca':
            # if we want to do rPCA here
            rpca = r_pca.R_pca(PCAdsgn)
            L,S = rpca.fit()
            
            if plot_distr:
                plt.figure()
                for bb in range(5):
                    plt.hist(L[:,bb],bins=bins,alpha=0.2)
                plt.suptitle('rPCA outputs')
            ##We treated rpca as a filtering step, so now we work solely with the low-rank component using the same procedure as above in the 'pca' block
            #definitely a more elegant way of merging these steps, good luck next grad student
            
            #Srcomp, Srevals, Srevecs = simple_pca(S)
            #Lrcomp, Lrevals, Lrevecs = simple_pca(L)
            pca = PCA()
            pca.fit(L)
#            
#            print('Using fit-transformed L')
            #Below shouldn't actually DO anything, since L is already from the output of rPCA and should already be aligned along its principal axes
            #But need pca() wrapper to get the coefficients, since I don't think r_pca includes it
            PCA_L = pca.fit_transform(L)
            
            self.PCA_d = pca
            self.PCA_inX = Xdo
            
            self.PCA_x = PCA_L
        
        
    
    def DEPRgen_GMM_priors(self,condit='OnT',mask_chann=False,band='Alpha'):
        #prior_change = self.pop_osc_mask[condit][dbo.feat_order.index(band)] * self.pop_osc_change[condit][dbo.feat_order.index(band)].reshape(-1,1)
        prior_change = self.pop_osc_change[condit][dbo.feat_order.index(band)].reshape(-1,1)
        if mask_chann:
            mask = self.median_mask
            prior_covar=np.dot(prior_change[mask],prior_change[mask].T) 
        else:
            prior_covar = np.dot(prior_change,prior_change.T)
        
        self.prior_covar = prior_covar
        return prior_covar
        
    
    def DEPRgen_GMM_Osc(self,GMM_stack):
        fvect = self.fvect
        
        feat_out = nestdict()
        
        for condit in self.condits:
            num_segs = GMM_stack[condit].shape[1]
            feat_out[condit] = np.zeros((257,num_segs,len(dbo.feat_order)))
            
            for ss in range(num_segs):
                
                feat_out[condit][:,ss,:] = dbo.calc_feats(10**(GMM_stack[condit][:,ss,:]/10).squeeze(),fvect)[0].T
                
                
        #THIS CURRENTLY HAS nans
        GMM_Osc_stack = feat_out
        
        return {'Stack':GMM_Osc_stack}
        
    def shape_GMM_dsgn(self,inStack_dict,band='Alpha',mask_channs=False):
        segs_feats = nestdict()
        
        inStack = inStack_dict['Stack']
        
        
        assert inStack['OnT'].shape[2] < 10        
        for condit in self.condits:
            num_segs = inStack[condit].shape[1]
            
            #I think this makes it nsegs x nchann x nfeats?
            segs_chann_feats = np.swapaxes(inStack[condit],1,0)
            
            #this stacks all the ALPHAS for all channels together
            
            
            if mask_channs:
                chann_mask = self.median_mask
            else:
                chann_mask = np.array([True] * 257)
            
            #CHOOSE ALPHA HERE
            
            if band == 'Alpha':
                band_idx = dbo.feat_order.index(band)
                segs_feats[condit] = segs_chann_feats[:,chann_mask,band_idx]
                
            elif band == 'All':
                segs_feats[condit] = segs_chann_feats[:,chann_mask,:]
                
            #segs_feats = np.reshape(segs_chann_feats,(num_segs,-1),order='F')
            
            #We want a 257 dimensional vector with num_segs observations
        
        return segs_feats

    '''Plot the population medians'''
    def pop_meds(self,response=True,pt='POOL'):
        print('Doing Population Meds/Mads on Oscillatory RESPONSES')
        
        #THIS IS THE OLD WAY: #dsgn_X = self.shape_GMM_dsgn(self.gen_GMM_Osc(self.gen_GMM_stack(stack_bl='normalize')['Stack']),band='All')
        if response:
            dsgn_X = self.osc_bl_norm[pt]
        else:
            dsgn_X = self.osc_stim[pt]
        
        X_med = nestdict()
        X_mad = nestdict()
        X_segnum = nestdict()
        #do some small simple crap here
        
        #Here we're averaging across axis zero which corresponds to 'averaging' across SEGMENTS
        for condit in self.condits:
            #Old version just does one shot median
            X_med[condit]= 10*np.median(dsgn_X[condit],axis=0)
            X_med[condit]= 10*np.mean(dsgn_X[condit],axis=0)
            
            
            # VARIANCE HERE

            X_mad[condit] = robust.mad(dsgn_X[condit],axis=0)
            #X_mad[condit] = np.var(dsgn_X[condit],axis=0)
            X_segnum[condit] = dsgn_X[condit].shape[0]
        
        self.Seg_Med = (X_med,X_mad,X_segnum)
        
        weigh_mad = 0.3
        try:
            self.median_mask = (np.abs(self.Seg_Med[0]['OnT'][:,2]) - weigh_mad*self.Seg_Med[1]['OnT'][:,2] >= 0)
        except:
            pdb.set_trace()
            
        #Do a quick zscore to zero out the problem channels
        chann_patt_zs = stats.zscore(X_med['OnT'],axis=0)
        outlier_channs = np.where(chann_patt_zs > 3)
        
    def pop_meds_jk(self,response=True):
        print('Doing Population Meds/Mads on Oscillatory RESPONSES - JackKnife Version')
        
        #THIS IS THE OLD WAY: #dsgn_X = self.shape_GMM_dsgn(self.gen_GMM_Osc(self.gen_GMM_stack(stack_bl='normalize')['Stack']),band='All')
        if response:
            dsgn_X = self.osc_bl_norm['POOL']
        else:
            dsgn_X = self.osc_stim['POOL']
        
        X_med = nestdict()
        X_mad = nestdict()
        X_segnum = nestdict()
        #do some small simple crap here
        
        #Here we're averaging across axis zero which corresponds to 'averaging' across SEGMENTS
        for condit in self.condits:
            # this version does jackknifing of the median estimate
            ensemble_med = dbo.jk_median(dsgn_X[condit])
            X_med[condit] = np.median(ensemble_med,axis=0)
            
            # VARIANCE HERE
            X_mad[condit] = np.std(ensemble_med,axis=0)
            #X_mad[condit] = robust.mad(dsgn_X[condit],axis=0)
            #X_mad[condit] = np.var(dsgn_X[condit],axis=0)
            X_segnum[condit] = dsgn_X[condit].shape[0]
        
        self.Seg_Med = (X_med,X_mad,X_segnum)
        
        weigh_mad = 0.3
        try:
            self.median_mask = (np.abs(self.Seg_Med[0]['OnT'][:,2]) - weigh_mad*self.Seg_Med[1]['OnT'][:,2] >= 0)
        except:
            pdb.set_trace()
            
        #Do a quick zscore to zero out the problem channels
        chann_patt_zs = stats.zscore(X_med['OnT'],axis=0)
        outlier_channs = np.where(chann_patt_zs > 3)

    def do_ICA_fullstack(self):
        rem_channs = False     
        print('ICA Time')
        ICA_inX = X_med['OnT']
        if rem_channs:
            ICA_inX[outlier_channs,:] = np.zeros_like(ICA_inX[outlier_channs,:])
        
        PCA_inX = np.copy(ICA_inX)
        
        #Go ahead and do PCA here since the variables are already here
        
        #PCA SECTION
        #pca = PCA()

        
        #ICA
        ica = FastICA(n_components=5)
        ica.fit(ICA_inX)
        self.ICA_d = ica
        self.ICA_inX = ICA_inX
        self.ICA_x = ica.fit_transform(ICA_inX)
        
        
    def DEPRgen_GMM_stack(self,stack_bl=''):
        state_stack = nestdict()
        state_labels = nestdict()
        
        
        for condit in self.condits:
            state_stack[condit] = []
            keyoi = keys_oi[condit][1]
            
            if stack_bl == 'add':
                raise ValueError
                pt_stacks = [val[condit][keyoi] for pt,val in self.feat_dict.items()] 
                base_len = len(pt_stacks)
                pt_labels = [keyoi] * base_len
                
                pt_stacks += [val[condit]['Off_3'] for pt,val in self.feat_dict.items()]
                pt_labels += ['Off_3'] * (len(self.feat_dict.keys()) - base_len)
                #pt_labels = [keyoi for pt,val in self.feat_dict.items()] + ['Off_3' for pt,val in self.feat_dict.items()]
                self.num_gmm_comps = 3
            elif stack_bl == 'normalize':
                #bl_stack = [val[condit]['Off_3'] for pt,val in self.feat_dict.items()]
                #bl_stack = np.concatenate(bl_stack,axis=1)
                bl_stack = {pt:val[condit]['Off_3'] for pt,val in self.feat_dict.items()}
                
                pt_stacks = [val[condit][keyoi] - np.repeat(np.median(bl_stack[pt],axis=1).reshape(257,1,1025),val[condit][keyoi].shape[1],axis=1) for pt,val in self.feat_dict.items()]
                pt_labels = [[keyoi] * gork.shape[1] for gork in pt_stacks]
                self.num_gmm_comps = 2
            else:
                pt_stacks = [val[condit][keyoi] for pt,val in self.feat_dict.items()]
                pt_labels = [keyoi] * len(pt_stacks)
                self.num_gmm_comps = 2
                
            state_stack[condit] = np.concatenate(pt_stacks,axis=1)
            state_labels[condit] = pt_labels
            
        
        GMM_stack = state_stack
        GMM_stack_labels = state_labels
        self.GMM_stack_labels = state_labels
        return {'Stack':GMM_stack,'Labels':GMM_stack_labels}
    
    
    def DEPRdo_response_PCA(self):
        pca = PCA()
        
        print("Using GMM Stack, which is baseline normalized within each patient")
        
        #This has both OnTarget and OffTarget conditions in the stack
        
        GMMStack = self.gen_GMM_stack(stack_bl='normalize')
        GMMOsc = self.gen_GMM_Osc(GMMStack['Stack'])
        
        
        
        PCA_inX = np.median(np.swapaxes(GMMOsc['Stack']['OnT'],0,1),axis=0)
        
        pca.fit(PCA_inX)
        self.PCA_d = pca
        self.PCA_inX = PCA_inX
        self.PCA_x = pca.fit_transform(PCA_inX)
        
    
    def plot_PCA_stuff(self):
        plt.figure();
        plt.subplot(221)
        plt.imshow(self.PCA_d.components_,cmap=plt.cm.jet,vmax=1,vmin=-1)
        plt.colorbar()
        
        plt.subplot(222)
        plt.plot(self.PCA_d.components_)
        plt.ylim((-1,1))
        plt.legend(['PC0','PC1','PC2','PC3','PC4'])
        plt.xticks(np.arange(0,5),['Delta','Theta','Alpha','Beta','Gamma1'])
        
        plt.subplot(223)
        plt.plot(self.PCA_d.explained_variance_ratio_)
        plt.ylim((0,1))
        
        for cc in range(2):
            fig=plt.figure()
            EEG_Viz.plot_3d_scalp(self.PCA_x[:,cc],fig,animate=False,unwrap=True)
            plt.title('Plotting component ' + str(cc))
            plt.suptitle('PCA rotated results for OnT')

    
    def plot_ICA_stuff(self):
        plt.figure()
        plt.subplot(221)
        plt.imshow(self.ICA_d.components_[:,:-1],cmap=plt.cm.jet,vmax=1,vmin=-1)
        plt.colorbar()
        plt.subplot(222)
        plt.plot(self.ICA_d.components_[:,:-1])
        plt.legend(['IC0','IC1','IC2','IC3','IC4'])
        plt.xticks(np.arange(0,5),['Delta','Theta','Alpha','Beta','Gamma1'])
        plt.subplot(223)
        
        plt.plot(self.ICA_d.mixing_)
        
        for cc in range(2):
            fig=plt.figure()
            EEG_Viz.plot_3d_scalp(self.ICA_x[:,cc],fig,animate=False,unwrap=True)
            plt.title('Plotting component ' + str(cc))
            plt.suptitle('ICA rotated results for OnT')
            
    
    '''gets us a histogram of the mads'''
    def band_mads(self):
        pass
    
    '''Plot distribution of change for individual bands'''
    def band_distr(self,do_moment='meds'):
        print('Plotting Distribution for Bands')
        
        meds = nestdict()
        mads = nestdict()
        
        marker=['o','s']
        color = ['b','g']
        plt.figure()
        ax2 = plt.subplot(111)
        for cc,condit in enumerate(['OnT','OffT']):
            #for bb in range(5):
            meds[condit] = self.Seg_Med[0][condit][:,:] #result here is 257(chann) x 5(bands)
            mads[condit] = self.Seg_Med[1][condit][:,:]
                #band_segnum[condit] = self.Seg_Med[2][condit]
            #plt.plot(np.arange(0,5)+(cc-0.5)/10,meds[condit][:,:].T,color[cc]+'.',markersize=20,alpha=0.05)
            #There's a way to do MATCHED CHANGES here!! TODO
                
                #plt.scatter((bb+(cc-0.5)/10)*np.ones_like(meds[condit][:,bb]),meds[condit][:,bb],marker=marker[cc],color=color[cc],s=100,alpha=0.2)
            #plt.boxplot(meds[condit][:,:],positions=np.arange(5)+(cc-0.5)/10,labels=dbo.feat_order)
            if do_moment == 'meds':
                parts = ax2.violinplot(meds[condit][:,:],positions=np.arange(5)+(cc-0.5)/10)
            elif do_moment == 'mads':
                parts = ax2.violinplot(mads[condit][:,:],positions=np.arange(5)+(cc - 0.5)/10)
                
            for pc in parts['bodies']:
                pc.set_facecolor(color[cc])
                pc.set_edgecolor(color[cc])
                #pc.set_linecolor(color[cc])
                                 
            #plt.ylim((-0.5,0.5))
        
        
        for bb in range(5):
            #rsres = stats.ranksums(meds['OnT'][:,bb],meds['OffT'][:,bb])
            #rsres = stats.wilcoxon(meds['OnT'][:,bb],meds['OffT'][:,bb])
            rsres = stats.ttest_ind(10**(meds['OnT'][:,bb]/10),10**(meds['OffT'][:,bb]/10))
            print(rsres)
        
        #plt.suptitle(condit)
        plt.ylim((-50,50))
        plt.hlines(0,-1,4,linestyle='dotted')
        plt.legend(['OnTarget','OffTarget'])
        
    def plot_meds(self,band='Alpha',flatten=True,condits=['OnT','OffT']):
        print('Doing Population Level Medians and MADs')
        
        band_median = {key:0 for key in self.condits}
        band_mad = {key:0 for key in self.condits}
        band_segnum = {key:0 for key in self.condits}
        
        
        if band == 'DSV':
            #lridge = [-0.00583578, -0.00279751,  0.00131825,  0.01770169,  0.01166687]
            #rridge = [-1.06586005e-02,  2.42700023e-05,  7.31445236e-03,  2.68723035e-03,-3.90440108e-06]
            doridge = np.array([-0.00583578, -0.00279751,  0.00131825,  0.01770169,  0.01166687])
            #doridge = np.array([-1.06586005e-02,  2.42700023e-05,  7.31445236e-03,  2.68723035e-03,-3.90440108e-06])
            doridge = doridge/np.linalg.norm(doridge)
            band_idx = np.array([0,1,2,3,4])
        else:
            band_idx = dbo.feat_order.index(band)
            doridge = [0,0,0,0,0]
            doridge[band_idx] = 1
            #rridge = [0,0,0,0,0]
            
        
        
        for condit in self.condits:
            band_median[condit] = np.dot(self.Seg_Med[0][condit][:,:],doridge)
            band_mad[condit] = self.Seg_Med[1][condit][:,band_idx]
            band_segnum[condit] = self.Seg_Med[2][condit]
            
            #band_mad[condit] = self.Seg_Med[1][condit][:,band_idx]
            #band_segnum[condit] = self.Seg_Med[2][condit]
        
        
        
        #Plot the Medians across channels
        plt.figure()
        plt.subplot(211)
        serr_med = {key:0 for key in self.condits}
        for condit in self.condits:
            
            plt.plot(band_median[condit],label=condit)
            serr_med[condit] = 1.48*band_mad[condit]/np.sqrt(band_segnum[condit])
            
            plt.fill_between(np.arange(257),band_median[condit] - serr_med[condit],band_median[condit] + serr_med[condit],alpha=0.4)

        plt.hlines(0,0,256)                    
        plt.title('Medians across Channels')
        plt.legend()
        plt.suptitle(band)
        
        ##
        # Do Wilcoxon signed rank test
        if 'OffT' in band_median.keys():
            WCSRtest = stats.wilcoxon(band_median['OnT'],band_median['OffT'])
            print(WCSRtest)
        
            # This is the plot of MADs
            plt.subplot(212)
            plt.plot(serr_med['OnT'],label='OnT')
            plt.plot(serr_med['OffT'],label='OffT')
            plt.title('Normed MADs across Channels')
            plt.legend()
        
        #plot EEG TOPOLOGICAL change for conditions
        for condit in self.condits:
            fig = plt.figure()
            #This is MEDS
            EEG_Viz.plot_3d_scalp(band_median[condit],fig,label=condit + '_med',animate=False,clims=(-0.2,0.2),unwrap=flatten)
            plt.suptitle('Median of Cortical Response across all ' + condit + ' segments | Band is ' + band)
            
        for condit in self.condits:
            #let's plot the exterior matrix for this
            fig = plt.figure()
            band_corr_matr = band_median[condit].reshape(-1,1) * band_median[condit].reshape(-1,1).T
            #pdb.set_trace()
            plt.imshow(band_corr_matr,vmin=-0.01,vmax=0.05)
            plt.colorbar()
        
        #plot the scalp EEG changes
        for condit in self.condits:
            fig = plt.figure()
            #this is MADs
            EEG_Viz.plot_3d_scalp(band_mad[condit],fig,label=condit + '_mad',animate=False,unwrap=flatten,clims=(0,1.0))
            plt.suptitle('MADs of Cortical Response across all ' + condit + ' segments | Band is ' + band)
        
        plt.suptitle(band)
        
        #Finally, for qualitative, let's look at the most consistent changes
        #This is the MASKED EEG channels
        if 0:
            for condit in self.condits:
                weigh_mad = 0.4
                fig = plt.figure()
                masked_median = self.Seg_Med[0][condit][:,band_idx] * (np.abs(self.Seg_Med[0][condit][:,band_idx]) - weigh_mad*self.Seg_Med[1][condit][:,band_idx] >= 0).astype(np.int)
                EEG_Viz.plot_3d_scalp(masked_median,fig,label=condit + '_masked_median',animate=False,clims=(-0.1,0.1))
                plt.suptitle('Medians with small variances (' + str(weigh_mad) + ') ' + condit + ' segments | Band is ' + band)
            
        
        olap = {key:0 for key in self.condits}
        
        ## Figure out which channels have overlap
        for condit in self.condits:
            olap[condit] = np.array((band_median[condit],band_median[condit] - band_mad[condit]/np.sqrt(band_segnum[condit]),band_median[condit] + band_mad[condit]/np.sqrt(band_segnum[condit])))
        
        for cc in range(257):
            np.hstack((olap['OnT'][1][cc],olap['OnT'][2][cc]))
        
        #find the channels that do overlap
        #allcond_vec = np.array([[self.Seg_Med[0]['OnT'][:,band_idx] + ont_semed,self.Seg_Med[0]['OnT'][:,band_idx] - ont_semed],[self.Seg_Med[0]['OffT'][:,band_idx] + offt_semed,self.Seg_Med[0]['OffT'][:,band_idx] - offt_semed]])
        
        
        # TODO
        # THID SHOULD BE MOVED TO A SEPARATE METHOD
        #do a sweep through to find the channels that don't overlap
#        sweep_range = np.linspace(0,0.11,100)
#        
#        ont_over = np.zeros_like(sweep_range)
#        offt_over = np.zeros_like(sweep_range)
#        
#        for ss,sr in enumerate(sweep_range):
#            ont_over[ss] = sum(serr_med['OnT'] > sr)
#            offt_over[ss] = sum(serr_med['OffT'] > sr)
#        plt.suptitle(band)
#        
#        plt.figure()
#        plt.plot(sweep_range,ont_over)
#        plt.plot(sweep_range,offt_over)
#        plt.suptitle(band)
        
    def train_GMM(self):
        #shape our dsgn matrix properly
        intermed_X = self.gen_GMM_Osc(self.gen_GMM_stack(stack_bl='normalize')['Stack'])
        dsgn_X = self.shape_GMM_dsgn(intermed_X,mask_channs=True)
        

        #let's stack the two together and expect 3 components?
        #this is right shape: (n_samples x n_features)
        
        #break out our dsgn matrix
        condit_dict = [val for key,val in dsgn_X.items()]
        
        full_X = np.concatenate(condit_dict,axis=0)
        
        #setup our covariance prior from OUR OTHER ANALYSES
        
        covariance_prior = self.gen_GMM_priors(mask_chann=True)
        
        #when below reg_covar was 1e-1 I got SOMETHING to work
        #gMod = mixture.BayesianGaussianMixture(n_components=self.num_gmm_comps,mean_prior=self.Seg_Med[0]['OnT'],mean_precision_prior=0.1,covariance_type='full',covariance_prior=np.dot(self.Seg_Med[1]['OnT'].reshape(-1,1),self.Seg_Med[1]['OnT'].reshape(-1,1).T),reg_covar=1,tol=1e-2)#,covariance_prior=covariance_prior_alpha)
        
        #BAYESIAN GMM version
        gMod = mixture.BayesianGaussianMixture(n_components=self.num_gmm_comps,mean_precision_prior=0.1,covariance_type='full',reg_covar=1e-6,tol=1e-2)#,covariance_prior=covariance_prior_alpha)
        condit_mean_priors = [np.median(rr) for rr in condit_dict]
        
        gMod.means_ = condit_mean_priors
        gMod.covariance_prior_ = covariance_prior
        
        #STANDARD GMM        
        #gMod = mixture.GaussianMixture(n_components=self.num_gmm_comps)
        
        
        try:
            gMod.fit(full_X)
        except Exception as e:
            print(e)
            pdb.set_trace()
            
        
        self.GMM = gMod
        
        self.predictions = gMod.predict(full_X)
        self.posteriors = gMod.predict_proba(full_X)
    
    def train_newGMM(self,mask=False):
        num_segs = self.SVM_stack.shape[0]
        
        #generate a mask
        if mask:
            sub_X = self.SVM_stack[:,self.median_mask,:]
            dsgn_X = sub_X.reshape(num_segs,-1,order='C')
        else:
            dsgn_X = self.SVM_stack.reshape(num_segs,-1,order='C')
        
        
        #doing a one class SVM
        #clf = svm.OneClassSVM(nu=0.1,kernel="rbf", gamma=0.1)
        #clf = svm.LinearSVC(penalty='l2',dual=False)
        #gMod = mixture.BayesianGaussianMixture(n_components=3,mean_precision_prior=0.1,covariance_type='full',reg_covar=1e-6,tol=1e-2)#,covariance_prior=covariance_prior_alpha)
        gMod = mixture.BayesianGaussianMixture(n_components=3)
        
        #split out into test and train
        Xtr,Xte,Ytr,Yte = sklearn.model_selection.train_test_split(dsgn_X,self.SVM_labels,test_size=0.33)
        

        
        gMod.fit(Xtr,Ytr)
        
        #predict IN training set
        predlabels = gMod.predict(Xte)
        
        plt.figure()
        plt.plot(Yte,label='test')
        plt.plot(predlabels,label='predict')
        
        print(np.sum(np.array(Yte) == np.array(predlabels))/len(Yte))
        
    def OBStrain_SVM(self,mask=False):
        num_segs = self.SVM_stack.shape[0]
        
        #generate a mask
        if mask:
            #what mask do we want?
            #self.SVM_Mask = self.median_mask
            self.SVM_Mask = np.zeros((257,)).astype(bool)
            self.SVM_Mask[[238,237]] = True
            
            sub_X = self.SVM_stack[:,self.SVM_Mask,:]
            dsgn_X = sub_X.reshape(num_segs,-1,order='C')
        else:
            dsgn_X = self.SVM_stack.reshape(num_segs,-1,order='C')
        
        #Learning curve
        print('Learning Curve')
        tsize,tscore,vscore = learning_curve(svm.LinearSVC(penalty='l2',dual=False),dsgn_X,self.SVM_labels,train_sizes=np.linspace(0.2,1.0,5),cv=5)
        plt.figure()
        plt.plot(tsize,np.mean(tscore,axis=1))
        plt.plot(tsize,np.mean(vscore,axis=1))



        #doing a one class SVM
        #clf = svm.OneClassSVM(nu=0.1,kernel="rbf", gamma=0.1)
        clf = svm.LinearSVC(penalty='l2',dual=False)
        
        #split out into test and train
        Xtr,Xte,Ytr,Yte = sklearn.model_selection.train_test_split(dsgn_X,self.SVM_labels,test_size=0.33)
        

        
        clf.fit(Xtr,Ytr)
        
        #predict IN training set
        predlabels = clf.predict(Xte)
        
        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(Yte,label='test')
        plt.plot(predlabels,label='predict')
        simple_accuracy = np.sum(np.array(Yte) == np.array(predlabels))/len(Yte)
        plt.title('SVM Results with Mask:' + str(mask) + ' ; Accuracy: ' + str(simple_accuracy))
        plt.legend()
        
        pickle.dump(clf,open('/tmp/SVMModel_l2','wb'))
        
        plt.subplot(2,1,2)
        conf_matrix = confusion_matrix(predlabels,Yte)
        plt.imshow(conf_matrix)
        plt.yticks(np.arange(0,3),['OFF','OffT','OnT'])
        plt.xticks(np.arange(0,3),['OFF','OffT','OnT'])
        plt.colorbar()
        
        
        self.SVM = clf
        self.SVM_dsgn_X = dsgn_X
        self.SVM_test_labels = predlabels
    
    def assess_dynamics(self,band='Alpha'):
        band_idx = dbo.feat_order.index(band)
        self.OnT_v_OffT_MAD()
        
        
        # Now, move on to plotting
        for stat in ['Med','MAD']:
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(self.Var_Meas['OnT'][stat][:,band_idx],fig,clims=(0,0),label='OnT '+ stat,unwrap=True)
            plt.suptitle('Non-normalized Power ' + stat + ' in ' + band + ' OnT')
            
            plt.figure()
            plt.bar(np.arange(1,258),self.Var_Meas['OnT'][stat][:,band_idx])
            
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(self.Var_Meas['OffT'][stat][:,band_idx],fig,clims=(0,0),label='OffT ' + stat,unwrap=True)
            plt.suptitle('Non-normalized Power ' + stat + ' in ' + band + ' OffT')
            fig = plt.figure()
            EEG_Viz.plot_3d_scalp(self.Var_Meas['OFF'][stat][:,band_idx],fig,clims=(0,0),label='OFF ' + stat,unwrap=True)
            plt.suptitle('Non-normalized Power ' + stat + ' in ' + band + ' OFF')
            
            plt.figure()
            plt.subplot(211)
            plt.hist(self.Var_Meas['OnT']['Med'][:,band_idx],label='OnT',bins=30)
            plt.hist(self.Var_Meas['OFF']['Med'][:,band_idx],label='OFF',bins=30)
            plt.title('Distributions of Medians')
            
            plt.subplot(212)
            plt.hist([self.Var_Meas['OnT']['MAD'][:,band_idx],self.Var_Meas['OFF']['MAD'][:,band_idx]],label=['OnT','OFF'],bins=30)
            #plt.hist(self.Var_Meas['OFF']['MAD'][:,band_idx],label='OFF',bins=30)
            plt.title('Distributions of MADs')
            plt.legend()
    
    def assess_binSVM(self,mask=False):
        num_segs = self.SVM_stack.shape[0]
        print('DOING BINARY')
        #generate a mask
        if mask:
            #what mask do we want?
            #self.SVM_Mask = self.median_mask
            self.SVM_Mask = np.zeros((257,)).astype(bool)
            self.SVM_Mask[np.arange(216,239)] = True
            
            sub_X = self.SVM_stack[:,self.SVM_Mask,:]
            dsgn_X = sub_X.reshape(num_segs,-1,order='C')
        else:
            dsgn_X = self.SVM_stack.reshape(num_segs,-1,order='C')
        
        #doing a one class SVM
        #clf = svm.OneClassSVM(nu=0.1,kernel="rbf", gamma=0.1)
    
        #get rid of ALL OFF, and only do two labels
        OFFs = np.where(self.SVM_labels == 'OFF')
        dsgn_X = np.delete(dsgn_X,OFFs,0)
        SVM_labels = np.delete(self.SVM_labels,OFFs,0)

        #Just doing a learning curve on the training data
        tsize,tscore,vscore = learning_curve(svm.LinearSVC(penalty='l2',dual=False,C=1),dsgn_X,SVM_labels,train_sizes=np.linspace(0.1,1,20),shuffle=True,cv=5,random_state=12342)
        plt.figure()
        plt.plot(tsize,np.mean(tscore,axis=1))
        plt.plot(tsize,np.mean(vscore,axis=1))
        plt.legend(['Training Score','Cross-validation Score'])
    
    ''' Learning curve for the Binary SVM'''
    def learning_binSVM(self,mask=False):
        label_map = {'OnT':1,'OffT':0}
    
        SVM_stack = np.concatenate([self.osc_bl_norm['POOL'][condit] for condit in self.condits],axis=0)
        SVM_labels = np.concatenate([[label_map[condit] for seg in self.osc_bl_norm['POOL'][condit]] for condit in self.condits],axis=0)
        num_segs = SVM_stack.shape[0]
        
        dsgn_X = SVM_stack.reshape(num_segs,-1,order='C')
        
        print('DOING BINARY - Learning Curve')
        tsize,tscore,vscore = learning_curve(svm.LinearSVC(penalty='l2',dual=False,C=1),dsgn_X,SVM_labels,train_sizes=np.linspace(0.1,1,20),shuffle=True,cv=10,random_state=12342)
        plt.figure()
        plt.plot(tsize,np.mean(tscore,axis=1))
        plt.plot(tsize,np.mean(vscore,axis=1))
        plt.legend(['Training Score','Cross-validation Score'])
    
    ''' The old method of generating a stack from our observations
    May more convoluted, vestiges of code from attempts to do GMM classification
    Need to phase this out completely
    Except it gives results in $\gamma$ that make more sense superficially
    '''
    def OBS_SVM_dsgn(self,do_plot=False):
        label_map = self.label_map
        
        # PREEMPT WITH OLD WAY HERE
        dsgn_X = self.shape_GMM_dsgn(self.gen_GMM_Osc(self.gen_GMM_stack(stack_bl='normalize')['Stack']),band='All')
        ALLT_dsgn_X = np.concatenate([dsgn_X[c] for c in ['OnT','OffT']],axis=0)
        
        flat_dsgn_Y = np.concatenate([[c for a in dsgn_X[c]] for c in ['OnT','OffT']],axis=0) # attempt to handle labels
        num_segs = ALLT_dsgn_X.shape[0]
        flat_dsgn_X = ALLT_dsgn_X.reshape(num_segs,-1,order='C')
        #dsgn_Y = np.concatenate([[label_map[condit] for seg in self.osc_bl_norm['POOL'][condit]] for condit in self.condits],axis=0)
        flat_dsgn_X[flat_dsgn_X > 1e300] = 0
        #ipdb.set_trace()
        if do_plot:
            #collapse along all segments and channels
            plot_stack = flat_dsgn_X.swapaxes(0,2).reshape(5,-1,order='C')
            plt.figure()
            try: sns.violinplot(data=plot_stack,positions=np.arange(5))
            except: ipdb.set_trace()
        
        return flat_dsgn_X, flat_dsgn_Y, num_segs
    
    ''' This function retrieves a design matrix from the pooled observations '''
    def stack_dsgn(self,do_plot=False):
        label_map = self.label_map
        
        # New method here...
        SVM_stack = np.concatenate([self.osc_bl_norm['POOL'][condit] for condit in self.condits],axis=0)
        num_segs = SVM_stack.shape[0]
        flat_dsgn_X = SVM_stack.reshape(num_segs,-1,order='C') #Here we're collapsing all the 5 Osc x 256 Chann features
         
        dsgn_Y = np.concatenate([[label_map[condit] for seg in self.osc_bl_norm['POOL'][condit]] for condit in self.condits],axis=0)
        
        
        if do_plot:
            #collapse along all segments and channels
            plot_stack = dsgn_X['OnT'].swapaxes(0,2).reshape(5,-1,order='C')
            plt.figure()
            for ii in range(5):
                sns.violinplot(data=plot_stack.T,positions=np.arange(5))
        
        self.SVM_raw_stack = SVM_stack
        return flat_dsgn_X, dsgn_Y, num_segs
                                    
    ''' WIP SVM masked classifier where channels can be toggled '''
    def masked_SVM(self):
        #generate a mask
    
        #what mask do we want?
        self.SVM_Mask = np.zeros((257,)).astype(bool)
        self.SVM_Mask[np.arange(216,239)] = True
        
        sub_X = self.SVM_stack[:,self.SVM_Mask,:]
        dsgn_X = sub_X.reshape(num_segs,-1,order='C')

    ''' Train our Binary SVM '''
    def train_binSVM(self,mask=False):
        self.bin_classif = nestdict()
        
        dsgn_X, SVM_labels, num_segs = self.stack_dsgn()
        
        #Next, we want to split out a validation set
        Xtr,self.Xva,Ytr,self.Yva = sklearn.model_selection.train_test_split(dsgn_X,SVM_labels,test_size=0.7,shuffle=True,random_state=None)
        
        #Next, we want to do CV learning on just the training set
        #Ensemble variables
        big_score = []
        coeffs = []
        models = []
        
        #Parameters for CV
        nfold = 10
        cv = StratifiedKFold(n_splits=nfold)
        for train,test in cv.split(Xtr,Ytr):
            clf = svm.LinearSVC(penalty='l1',dual=False,C=1) #l2 works well here...
            mod_score = clf.fit(Xtr[train],Ytr[train]).score(Xtr[test],Ytr[test])
            outpred = clf.predict(Xtr[test])
            coeffs.append(clf.coef_)
            big_score.append(mod_score)
            models.append(clf)
            #Maybe do ROC stuff HERE? TODO
            
        #Plot the big score for the CVs
        plt.figure()
        plt.plot(big_score)
        plt.title('Plotting the fit scores for the CV training procedure')
        
        # Find the best model
        best_model_idx = np.argmax(big_score)
        best_model = models[best_model_idx]
        self.bin_classif['Model'] = best_model
        self.bin_classif['Coeffs'] = coeffs
        self.cv_folding = nfold
        
        #if we want to keep all of our models we can do that here
        self.cv_bin_classif = models
    
    '''This function ASSESSES the best SVM classifier using a bootstrap procedure'''
    def bootstrap_binSVM(self):
        best_model = self.bin_classif
        
        #randomlt sample the validation set
        validation_accuracy = []
        rocs_auc = []
        total_segments = self.Xva.shape[0]
        for ii in range(100):
            Xva_ss,Yva_ss = resample(self.Xva,self.Yva,n_samples=int(round(total_segments * 0.6)),replace=False)
            validation_accuracy.append(best_model['Model'].score(Xva_ss,Yva_ss))
            predicted_Y = best_model['Model'].predict(Xva_ss)
            fpr,tpr,_ = roc_curve(Yva_ss,predicted_Y)
            rocs_auc.append(auc(fpr, tpr))
            
        plt.figure()
        plt.subplot(311)
        plt.hist(validation_accuracy)
        plt.subplot(312)
        plt.plot(fpr,tpr)
        plt.subplot(313)
        plt.hist(rocs_auc)
        plt.suptitle('Bootstrapped SVM Assessment')
    
    '''One shot assessment of SVM classifier'''
    def oneshot_binSVM(self):
        best_model = self.bin_classif
        #Plotting of confusion matrix and coefficients
        # Validation set assessment now
        
        validation_accuracy = best_model['Model'].score(self.Xva,self.Yva)
        Ypred = best_model['Model'].predict(self.Xva)
        print(validation_accuracy)
        
        plt.figure()
        plt.subplot(1,2,1)
        #confusion matrix here
        conf_matrix = confusion_matrix(Ypred,self.Yva)
        plt.imshow(conf_matrix)
        plt.yticks(np.arange(0,2),['OffT','OnT'])
        plt.xticks(np.arange(0,2),['OffT','OnT'])
        plt.colorbar()
        
        plt.subplot(2,2,2)
        coeffs = np.array(best_model['Coeffs']).squeeze().reshape(self.cv_folding,257,5,order='C')
        #plt.plot(coeffs,alpha=0.2)
        plt.plot(np.median(coeffs,axis=0))
        plt.title('Plotting Median Coefficients for CV-best Model performance')
        
        plt.subplot(2,2,4)
        plt.plot(np.median(np.median(coeffs,axis=0),axis=0))
        plt.suptitle('Oneshot SVM Assessment')
        
        #self.SVM_coeffs = coeffs
    
    def assess_binSVM(self):
        best_model = self.bin_classif
        
        for ii in range(100):
            Xrs,Yrs = resample(self.Xva, self.Yva,100)
            valid_accuracy = best_model.score(Xva,Yva)
            
    '''Analysis of the binary SVM coefficients should be here'''    
    def analyse_binSVM(self, feature_weigh = False):

        #What exactly is this trying to do here??
        # plt.figure()
        # for ii in range(4):
        #     #plt.hist(self.bin_classif['Model'].coef_[:,ii])
        #     plt.scatter(ii,self.bin_classif['Model'].coef_[:,ii])
        
        #BELOW (order = 'C') IS CORRECT since before, in the features, we collapse to a feature vector that is all 257 deltas, then all 257 thetas, etc...
        #So when we want to reshape that to where we are now, we have to either 'C': (5,257) where C means the last index changes fastest; or 'F': (257,5) where the first index changes fastest.
        
        if not feature_weigh:
            coeffs = stats.zscore(np.sum(np.abs(self.bin_classif['Model'].coef_.reshape(5,257,order='C')),axis=0)) #what we have here is a reshape where the FEATURE VECTOR is [257 deltas... 257 gammas]
            self.import_mask = coeffs > 0.2
            analysis_title = 'Pure Coefficients'
        else:
            avg_feat_value = np.abs(np.median(self.SVM_raw_stack,axis=0))
            #pdb.set_trace()
            coeffs = stats.zscore(np.sum(np.multiply(np.abs(self.bin_classif['Model'].coef_.reshape(5,257,order='C')),avg_feat_value.T),axis=0))
            self.import_mask = coeffs > 0.1
            analysis_title = 'Empirical Feature-weighed'
            
            
        #WTF is this shit?
        #avg_coeffs = np.mean(np.array(self.bin_classif['Coeffs']),axis=0).reshape(5,257,order='C')
        #coeffs = stats.zscore(np.sum(avg_coeffs**2,axis=0))

        plt.figure();plt.hist(coeffs,bins=50)#,range=(0,1))
        #plt.figure()

        EEG_Viz.plot_3d_scalp(coeffs,unwrap=True,alpha=0.4)
        EEG_Viz.plot_3d_scalp(self.import_mask.astype(np.int),unwrap=True,alpha=0.2)
        plt.suptitle('SVM Coefficients ' + analysis_title)
            
    '''Do a CV-level analysis here'''
    def analysis_CV_binSVM(self):
        pass
            
            
    # THE BELOW FUNCTION DOES NOT RUN, JUST HERE FOR REFERENCE AS THE SVM IS BEING RECODED ABOVE
    def OLDtrain_binSVM(self):
        #%% PLOT THE WHOLE DATA
        plt.figure()
        Yall = np.zeros(SVM_labels.shape[0]).astype(np.float)
        
        Yall[SVM_labels == 'OffTON'] = 0
        Yall[SVM_labels == 'OnTON'] = 1
        plt.imshow(Yall.reshape(1,-1),aspect='auto')
               
        
        #split out into test and train
        testing_size = 200/310
        print('Total segments: ' + str(dsgn_X.shape))
       
        Xtr,Xte,Ytr,Yte = sklearn.model_selection.train_test_split(dsgn_X,SVM_labels,test_size=testing_size,random_state=1230,shuffle=True)
        print('Training size ' + str(Xtr.shape))
        
        plt.figure()
        trYall = np.zeros(Ytr.shape[0]).astype(np.float)
        trYall[Ytr == 'OffTON'] = 0
        trYall[Ytr == 'OnTON'] = 1
        
        teYall = np.zeros(Yte.shape[0]).astype(np.float)
        teYall[Yte == 'OffTON'] = 0
        teYall[Yte == 'OnTON'] = 1
        
        plt.imshow(np.hstack((trYall,teYall)).reshape(1,-1),aspect='auto')
        
        #THIS HAS BEEN MOVED TO SEPARATE FUNCTION/METHOD IN THIS CLASS
        #Just doing a learning curve on the training data
        #tsize,tscore,vscore = learning_curve(svm.LinearSVC(penalty='l2',dual=False,C=1),Xtr,Ytr,train_sizes=np.linspace(0.4,1,10),shuffle=True,cv=5,random_state=0)
        #plt.figure()
        #plt.plot(tsize,np.mean(tscore,axis=1))
        #plt.plot(tsize,np.mean(vscore,axis=1))

        
        #classifier time itself
        clf = svm.LinearSVC(penalty='l2',dual=False,C=1)
        #Fit the actual algorithm
        
        big_score = []
        coeffs = []
        nfold = 50
        cv = StratifiedKFold(n_splits=nfold)
        
        rocs = []
        aucs = []
        plt.figure()
        for train,test in cv.split(Xtr,Ytr):
            mod_score = clf.fit(Xtr[train],Ytr[train]).score(Xtr[test],Ytr[test])
            outpred = clf.predict(Xtr[test])
            
            Ytestr = np.zeros(Ytr[test].shape[0]).astype(np.float)
            Ytestr[Ytr[test] == 'OffTON'] = 0
            Ytestr[Ytr[test] == 'OnTON'] = 1
            
            outpred[outpred == 'OffTON'] = 0
            outpred[outpred == 'OnTON'] = 1
            
            #pdb.set_trace()
            outpred = outpred.astype(np.float)
            fpr,tpr,thresholds = roc_curve(Ytestr,outpred)
            auc_perf = roc_auc_score(Ytestr,outpred)
            #rocs.append(roc_perf)
            aucs.append(auc_perf)
            plt.plot(fpr,tpr)
            
            coeffs.append(clf.coef_)
            #fpr,tpr,threshold = roc_curve(Ytr[test],probas)
            #roc_auc = auc(fpr,tpr)
            big_score.append(mod_score)
        plt.ylim((0,1))
        coeffs = np.array(coeffs).squeeze().reshape(nfold,5,-1)
        print(aucs)
        print('CV Scores: ' + str(big_score))
        plt.figure()
        #pdb.set_trace()
        #plt.plot(np.mean(coeffs,axis=0))
        for ii in range(nfold):
            plt.plot(coeffs[ii,:,:].T,alpha=0.2)
        plt.legend(['Delta','Theta','Alpha','Beta','Gamma'])
        plt.plot(np.mean(coeffs,axis=0).T,alpha=1)
        
        
        #%% Now do fit on the full training set
            
        #clf.fit(Xtr,Ytr)
        
        #predict IN training set
        predlabels = clf.predict(Xte)
        
        plt.figure()
        plt.subplot(2,1,1)
        
        Yten = np.copy(Yte)
        predlabelsn = np.copy(predlabels)
        
        #Fix labels for PLOTTING
        Yten[Yte == 'OffTON'] = 0
        Yten[Yte == 'OnTON'] = 1
        predlabelsn[predlabels == 'OffTON'] = 0
        predlabelsn[predlabels == 'OnTON'] = 1
        
        
        plt.imshow(np.vstack((Yten.astype(np.int),predlabelsn.astype(np.int))),aspect='auto')
        #plt.plot(Yte,label='test')
        #plt.plot(predlabels,label='predict')
        #simple_accuracy = np.sum(np.array(Yte) == np.array(predlabels))/len(Yte)
        score = clf.score(Xte,Yte)
        plt.title('SVM Results with Mask:' + str(mask) + ' ; Accuracy: ' + str(score))
        plt.legend()
        
        pickle.dump(clf,open('/tmp/SVMModel_l2','wb'))
        
        plt.subplot(2,2,3)
        plt.plot(clf.coef_.reshape(5,-1).T)
        plt.subplot(2,2,4)
        conf_matrix = confusion_matrix(predlabels,Yte)
        plt.imshow(conf_matrix)
        plt.yticks(np.arange(0,2),['OffT','OnT'])
        plt.xticks(np.arange(0,2),['OffT','OnT'])
        plt.colorbar()
        
        
        self.binSVM = clf
        self.binSVM_dsgn_X = dsgn_X
        self.binSVM_test_labels = predlabels
        
    
    def OBScompute_diff(self,take_mean=True):
        print('Computing Difference')
        assert len(self.condits) >= 2
        avg_psd = nestdict()
        avg_change = nestdict()
        var_psd = nestdict()
        
        for pt in self.do_pts:
            #avg_psd[pt] = defaultdict(dict)
            #avg_change[pt] = defaultdict(dict)
            for condit in self.condits:
               #average all the epochs together
                avg_psd[pt][condit] = {epoch:np.median(self.feat_dict[pt][condit][epoch],axis=1) for epoch in self.feat_dict[pt][condit].keys()}
                #if you want variance
                #var_psd[pt][condit] = {epoch:np.var(self.feat_dict[pt][condit][epoch],axis=1) for epoch in self.feat_dict[pt][condit].keys()}
                #if you want Mean Absolute Deviance
                var_psd[pt][condit] = {epoch:robust.mad(self.feat_dict[pt][condit][epoch],axis=1) for epoch in self.feat_dict[pt][condit].keys()}
                
                
                keyoi = keys_oi[condit][1]
                
                avg_change[pt][condit] = 10*(np.log10(avg_psd[pt][condit][keyoi]) - np.log10(avg_psd[pt][condit]['Off_3']))
                
                
        self.psd_change = avg_change
        self.psd_avg = avg_psd
        #This is really just a measure of how dynamic the underlying process is, not of particular interest for Aim 3.1, maybe 3.3
        self.psd_var = var_psd
        
    def NEWcompute_diff(self):
        avg_change = {pt:{condit:10*(np.log10(avg_psd[pt][condit][keys_oi[condit][1]]) - np.log10(avg_psd[pt][condit]['Off_3'])) for pt,condit in itertools.product(self.do_pts,self.condits)}}
   
    
    ''' Below are functions related to the oscillatory response characterization'''
    #This goes to the psd change average and computed average PSD across all available patients
    def pop_response(self):
        psd_change_matrix = nestdict()
        population_change = nestdict()
        #pop_psds = defaultdict(dict)
        #pop_psds_var = defaultdict(dict)
        
        #Generate the full PSDs matrix and then find the MEAN along the axes for the CHANGE
        for condit in self.condits:
            #Get us a matrix of the PSD Changes
            psd_change_matrix[condit] = np.array([rr[condit] for pt,rr in self.psd_change.items()])
            
            population_change[condit]['Mean'] = np.mean(psd_change_matrix[condit],axis=0)
            population_change[condit]['Var'] = np.var(psd_change_matrix[condit],axis=0)
            
            
        self.pop_psd_change = population_change
    
    def do_pop_stats(self,band='Alpha',threshold=0.4):
        #self.reliablePSD = nestdict()
        avgOsc = nestdict()
        varOsc = nestdict()
        varOsc_mask = nestdict()
        
        for condit in self.condits:
            
            #First, let's average across patients
            
            #band = np.where(np.logical_and(self.fvect > 14, self.fvect < 30))
            
            avgOsc[condit] = dbo.calc_feats(10**(self.pop_psd_change[condit]['Mean']/10),self.fvect)
            varOsc[condit] = dbo.calc_feats(10**(self.pop_psd_change[condit]['Var']),self.fvect)
            #varOsc_mask[condit] = np.array(np.sqrt(varOsc[condit]) < threshold).astype(np.int)
            
        self.pop_osc_change = avgOsc
        self.pop_osc_var = varOsc
        #self.pop_osc_mask = varOsc_mask

    def do_pop_mask(self,threshold):
        cMask = nestdict()
        
        #weighedPSD = np.divide(self.pop_stats['Mean'][condit].T,np.sqrt(self.pop_stats['Var'][condit].T))
        #do subtract weight
        #THIS SHOULD ONLY BE USED TO DETERMINE THE MOST USEFUL CHANNELS
        #THIS IS A SHITTY THING TO DO
        for condit in self.condits:
            pre_mask = np.abs(self.pop_change['Mean'][condit].T) - threshold*np.sqrt(self.pop_change['Var'][condit].T)
            #which channels survive?
            cMask[condit] = np.array(pre_mask > 0).astype(np.int)
            
        cMask['Threshold'] = threshold
        
        self.reliability_mask = cMask
        
            
    def plot_pop_stats(self):
        for condit in self.condits:
            plt.figure()
            plt.subplot(311)
            plt.plot(self.fvect,self.pop_stats['Mean'][condit].T)
            plt.title('Mean')
            plt.subplot(312)
            plt.plot(self.fvect,np.sqrt(self.pop_stats['Var'][condit].T)/np.sqrt(3))
            plt.title('Standard Error of the Mean; n=3')
            
            
            plt.subplot(313)
            #do divide weight
            
            
            
            reliablePSD = self.reliablePSD[condit]['cPSD']
            
            plt.plot(self.fvect,reliablePSD)
            plt.title('SubtrPSD by Pop 6*std')
            
            plt.xlabel('Frequency (Hz)')
            #plt.subplot(313)
            #plt.hist(weighedPSD,bins=50)
            
            plt.suptitle(condit + ' population level')
            
    
    def topo_wrap(self,band='Alpha',label='',condit='OnT',mask=False,animate=False):
        mainfig = plt.figure()
        
        if not mask:
            EEG_Viz.plot_3d_scalp(self.pop_osc_change[condit][dbo.feat_order.index(band)],mainfig,animate=animate,label=label,clims=(-1,1))
        else:
            try:
                EEG_Viz.plot_3d_scalp(self.pop_osc_mask[condit][dbo.feat_order.index(band)] * self.pop_osc_change[condit][dbo.feat_order.index(band)],mainfig,clims=(-1,1),animate=animate,label=label)
            except:
                pdb.set_trace()
        
            
        plt.suptitle(label)
        
        
    def topo_3d_chann_mask(self,band='Alpha',animate=True,pt='all',var_fixed=False,label=''):
        for condit in ['OnT','OffT']:
            #plot_vect = np.zeros((257))
            #plot_vect[self.reliablePSD[condit]['CMask']] = 1
            #plot_vect = self.reliablePSD[condit]['CMask']
            
            mainfig = plt.figure()
            
            band_idxs = np.where(np.logical_and(self.fvect > dbo.feat_dict[band]['param'][0], self.fvect < dbo.feat_dict[band]['param'][1]))
            if pt == 'all':
                if var_fixed:
                    plot_vect = self.reliablePSD[condit]['BandVect']['Vect']
                else:
                    
                    plot_vect = np.median(self.avgPSD[condit]['PSD'][:,band_idxs].squeeze(),axis=1)
                    #abov_thresh = np.median(self.reliablePSD[condit]['cPSD'][band_idxs,:],axis=1).T.reshape(-1)
                    
            else:
                #just get the patient's diff
                
                plot_vect = np.median(self.feat_diff[pt][condit][:,band_idxs].squeeze(),axis=1)
                
            
            EEG_Viz.plot_3d_scalp(plot_vect,mainfig,clims=(-3,3),label=condit+label,animate=animate)
            plt.title(condit + ' ' + pt + ' ' + str(var_fixed))
            plt.suptitle(label)
            
            self.write_sig(plot_vect)
            
    def write_sig(self,signature):
        for condit in ['OnT','OffT']:
            #plot_vect = np.zeros((257))
            #plot_vect[self.reliablePSD[condit]['CMask']] = 1
            #plot_vect = self.reliablePSD[condit]['CMask']
            
            #write this to a text file
            np.save('/tmp/' + condit + '_sig.npy',signature)
    
    def plot_diff(self):
        for pt in self.do_pts:
            plt.figure()
            plt.subplot(221)
            plt.plot(self.fvect,self.psd_change[pt]['OnT'].T)
            plt.title('OnT')
            
            plt.subplot(222)
            plt.plot(self.fvect,self.psd_change[pt]['OffT'].T)
            plt.title('OffT')
            
            plt.subplot(223)
            plt.plot(self.fvect,10*np.log10(self.psd_var[pt]['OnT']['BONT'].T))
            plt.title('BONT Variance')
            
            plt.subplot(224)
            plt.plot(self.fvect,10*np.log10(self.psd_var[pt]['OffT']['BOFT'].T))
            plt.title('BOFT Variance')
            
   
            plt.suptitle(pt)
    
    '''GMM Classification Modules here - Analysis is obsolete'''
    def GMM_train(self,condit='OnT'):
        #gnerate our big matrix of observations; Should be 256(chann)x4(feats)x(segxpatients)(observations)
        pass
    
        #this function will generate a big stack of all observations for a given condition across all patients
    
    def DEPRgen_OSC_stack(self,stack_type='all'):
        remap = {'Off_3':'OF','BONT':'ON','BOFT':'ON'}
        #big_stack = {key:0 for key in self.condits}
        big_stack = nestdict()
        for condit in self.condits:
            #want a stack for OnT and a stack for OffT
            big_stack[condit] = {remap[epoch]:{pt:self.osc_dict[pt][condit][epoch] for pt in self.do_pts} for epoch in keys_oi[condit]}
            
        if stack_type == 'all':
            pass
        
        self.big_stack_dict = big_stack
        return big_stack
    
    def DEPRextract_feats(self):
        donfft = self.donfft
        pts = self.do_pts
        for pt in pts:
            var_matr = defaultdict(dict)
            for condit in self.condits:
                med_matr = defaultdict(dict)
                var_matr = defaultdict(dict)
                
                med_feat_matr = defaultdict(dict)
                var_feat_matr = defaultdict(dict)
                
                #the size of this object should be
                for epoch in keys_oi[condit]:
                    data_matr = self.ts_data[pt][condit][epoch]
                    
                    #pdb.set_trace()
                    #plt.plot(data_dict[100].T)
                    #print(data_dict[100].shape)
                    # METHOD 1
                    method = 2
                    #Method one concatenates all segments then does the PSD
                    #this is not a good idea with the liberal preprocessing, since there are still artifacts
                    if method == 1:
                        #this next step gets rid of all the segments
                        data_dict = {ch:data_matr[ch,:,:].reshape(1,-1,order='F') for ch in range(data_matr.shape[0])}
                        
                        ret_Fsegs = dbo.gen_psd(data_dict,Fs=1000,nfft=donfft)
                        
                        #need to ordered, build list
                        allchanns = [ret_Fsegs[ch] for ch in range(257)]
                        
                        med_matr[epoch] = 10*np.log10(np.squeeze(np.array(allchanns)))
                        self.psd_trans[pt][condit][epoch] = ret_Fsegs
                    #in methods 2, we do the PSD of each segment individually then we find the median across all segments
                    elif method == 2:
                        
                        chann_seglist = np.zeros((257,int(donfft/2+1),data_matr.shape[2]))
                        chann_segfeats = np.zeros((257,len(dbo.feat_order),data_matr.shape[2]))
                        for ss in range(data_matr.shape[2]):
                            data_dict = {ch:data_matr[ch,:,ss] for ch in range(data_matr.shape[0])}
                            ret_f = dbo.gen_psd(data_dict,Fs=1000,nfft=donfft)
                            #Push ret_f, the dictionary, into the matrix we need for further processing
                            for cc in range(257):
                                chann_seglist[cc,:,ss] = ret_f[cc]
                                #go to each channel and find the oscillatory band feature vector
                                chann_segfeats[cc,:,ss] = dbo.calc_feats(ret_f[cc],self.fvect)
                            
                        self.psd_trans[pt][condit][epoch] = chann_seglist
                        self.Feat_trans[pt][condit][epoch] = chann_segfeats
                        
                        med_matr[epoch] = np.median(10*np.log10(chann_seglist),axis=2)
                        var_matr[epoch] = np.var(10*np.log10(chann_seglist),axis=2)
                    
                        med_feat_matr[epoch] = np.median(10*np.log10(chann_segfeats),axis=2)
                        var_feat_matr[epoch] = np.var(10*np.log10(chann_segfeats),axis=2)
                
                assert med_matr[keys_oi[condit][1]].shape == (257,donfft/2+1)
                diff_matr = med_matr[keys_oi[condit][1]] - med_matr['Off_3']
                diff_feat = med_feat_matr[keys_oi[condit][1]] - med_feat_matr['Off_3']
                
                #Put the PSD information in the class structures
                self.PSD_var[pt][condit] = var_matr
                self.PSD_diff[pt][condit] = diff_matr
                #Put the Oscillator feature vector in the class structure
                self.Feat_diff[pt][condit] = diff_feat
                self.Feat_var[pt][condit] = var_feat_matr

        plot = False
        if plot:
            plt.figure()
            plt.subplot(211)
            plt.plot(self.fvect,med_matr['Off_3'].T)
            plt.subplot(212)
            plt.plot(self.fvect,med_matr[keys_oi[condit][1]].T)        
        
    def DEPRplot_diff(self,pt='906',condit='OnT',varweigh=False):
        diff_matr = self.PSD_diff[pt][condit]
        
        if varweigh:
            var_matr = self.PSD_var[pt][condit][keys_oi[condit][1]]
            plotdiff = (diff_matr / np.sqrt(var_matr)).T
            varweightext = ' WEIGHTED WITH 1/VAR'
        else:
            plotdiff = diff_matr.T
            varweightext = ''
            
        plt.figure()
        plt.plot(self.fvect,plotdiff,alpha=0.2)
        
        plt.xlim((0,150))
        plt.title('Difference from Pre-Stim to Stim')
        plt.suptitle(pt + ' ' + condit + varweightext)
    
    def plot_ontvsofft(self,pt='906'):
        if 'OffT' not in self.condits:
            raise ValueError
            
        condit_diff = (self.PSD_diff[pt]['OnT'] - self.PSD_diff[pt]['OffT']).T
        plt.figure()
        plt.plot(self.fvect,condit_diff,alpha=0.2)
        plt.xlim((0,150))
        plt.title('Difference from OnT and OffT')
        plt.suptitle(pt)
    
    def plot_chann_var(self,pt='906',condit='OnT'):
        plt.figure()
        plt.subplot(121)
        plt.plot(self.PSD_var[pt][condit]['Off_3'].T)
        plt.xlim((0,150))
        plt.subplot(122)
        plt.plot(self.PSD_var[pt][condit][keys_oi[condit][1]].T)
        plt.xlim((0,150))
        
        plt.suptitle(pt + ' ' + condit)
        
    #This function quickly gets the power for all channels in each band
    def Full_feat_band(self):
        pass
    
    def extract_coher_feats(self,do_pts=[],do_condits=['OnT'],epochs='all'):
        if do_pts == []:
            do_pts=self.do_pts
            
        PLV_dict = nestdict()
        CSD_dict = nestdict()
        if do_condits == []:
            do_condits = self.condits

        for pt in do_pts:
            for condit in do_condits:
                if epochs == 'all':
                    do_epochs = keys_oi[condit]
                else:
                    do_epochs = epochs
                    
                for epoch in do_epochs:
                    print('Doing ' + pt + condit + epoch)
                    data_matr = self.ts_data[pt][condit][epoch]
                    data_dict = {ch:data_matr[ch,:,:].squeeze() for ch in range(self.chann_dim)}
                    CSD_dict[pt][condit][epoch], PLV_dict[pt][condit][epoch] = dbo.gen_coher(data_dict,Fs=self.fs,nfft=2**9,polyord=self.polyorder)
                    
        print('Done with coherence... I guess...')
        return CSD_dict,PLV_dict
        

    
    #This is the main wrapper function that leads us to the MNE python plot_topomap
    # The topomap is very misleading, so we try not to use it due to its strange interpolation
    def plot_topo(self,vect,vmax=2,vmin=-2,label='',):
        plt.figure()
        mne.viz.plot_topomap(vect,pos=self.eeg_locs.pos[:,[0,1]],vmax=vmax,vmin=vmin,image_interp='none')
        plt.suptitle(label)
        
    
    #Here, we'll plot the PSDs for channels of interest for the conditions of interest
    def psd_stats(self,chann_list=[]):
        self.view_PSDs(chann_list=chann_list,zoom_in=True)
        
    #This function is to just show the raw PSDs of each of the experimental conditions collected
    def view_PSDs(self,zoom_in=True,chann_list=[],plot_var=False):
        print('Showing raw PSDs')
        avg_psd = nestdict()
        var_psd = nestdict()
        
        if chann_list == []:
            chann_list = np.arange(256)
        else:
            chann_list = np.array(chann_list)
            
        f_vect = np.linspace(0,500,2**10+1)
        
        for pt in self.do_pts:
            #avg_psd[pt] = defaultdict(dict)
            #avg_change[pt] = defaultdict(dict)
            for condit in self.condits:
               #average all the epochs together
                avg_psd[pt][condit] = {epoch:np.median(self.feat_dict[pt][condit][epoch],axis=1) for epoch in self.feat_dict[pt][condit].keys()}
                if plot_var:
                    var_psd[pt][condit] = {epoch:robust.mad(self.feat_dict[pt][condit][epoch],axis=1) for epoch in self.feat_dict[pt][condit].keys()}
                
            
            psd_fig =  plt.figure()
            plt.subplot(2,2,1)
            plt.plot(f_vect,10*np.log10(avg_psd[pt]['OnT']['Off_3'][chann_list,:].T))
            plt.title('OnT-Pre')
            
            plt.subplot(2,2,2)
            plt.plot(f_vect,10*np.log10(avg_psd[pt]['OffT']['Off_3'][chann_list,:].T))
            plt.title('OffT-Pre')
            
            plt.subplot(2,2,3)
            plt.plot(f_vect,10*np.log10(avg_psd[pt]['OnT']['BONT'][chann_list,:].T))
            plt.title('BONT')
            
            plt.subplot(2,2,4)
            plt.plot(f_vect,10*np.log10(avg_psd[pt]['OffT']['BOFT'][chann_list,:].T))
            plt.title('BOFT')
            
            if plot_var:
                plt.figure(psd_fig.number)
                plt.subplot(2,2,1)
                plt.fill_between(f_vect,10*np.log10(var_psd[pt]['OnT']['Off_3'][chann_list,:].T))
            
            if zoom_in:
                for ii in range(1,5):
                    plt.subplot(2,2,ii)
                    plt.xlim(0,40)
                    plt.ylim(-20,10)
            
            plt.suptitle(pt)
            
    '''Coherence Statistics here'''
    def coher_stat(self,pt_list=[],chann_list=[]):
        return self.extract_coher_feats(do_pts=pt_list,do_condits=['OnT','OffT'])
    
    
    ''' GRAVEYARD '''
    
    def OBSBLWEIRDcompute_response(self,combine_baselines=True,plot=False):
        if combine_baselines:
            baseline = {pt:np.median(self.combined_BL[pt],axis=0) for pt in self.do_pts}
        else:
            baseline = {pt:np.median(self.osc_dict[pt][condit][keys_oi[condit][0]],axis=0) for pt in self.do_pts}
            
        self.osc_bl_norm = {pt:{condit:(self.osc_dict[pt][condit][keys_oi[condit][1]] - baseline[pt]) for condit in self.condits} for pt in self.do_pts}
        
        if plot:
            plt.figure()