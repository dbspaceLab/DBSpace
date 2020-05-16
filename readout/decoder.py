#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:47:25 2020

@author: virati
NEW classes for readout training, testing, and validation
"""
import sklearn
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve, average_precision_score, auc
from sklearn.metrics import roc_auc_score

import warnings
from collections import defaultdict
import itertools as itt
from itertools import compress

import json

import ipdb

import numpy as np
import scipy.stats as stats
import scipy.signal as sig

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import random

#import sys
#sys.path.append('/home/virati/Dropbox/projects/Research/MDD-DBS/Ephys/DBSpace/')
import DBSpace as dbo
from DBSpace import nestdict

from sklearn import linear_model

default_params = {'CrossValid':10}

import seaborn as sns
#sns.set_context("paper")

sns.set(font_scale=4)
sns.set_style("white")

import time
import copy
import pdb

#%%            
class RO:
    #Parent readout class
    # Some very fixed constants up here
    circ = 'day'
    ch_num = 2
    
    def __init__(self,BRFrame,ClinFrame,pts,clin_measure='HDRS17'):
        self.YFrame = BRFrame
        self.CFrame = ClinFrame
        self.pts = pts
        self.c_meas = clin_measure
        self.fvect = self.YFrame.data_basis['F']
    
    '''Filter out the recordings we want'''
    def filter_recs(self,rec_class='main_study'):
        if rec_class == 'main_study':
            filter_phases = dbo.Phase_List(exprs='ephys')
            self.active_rec_list = [rec for rec in self.YFrame.file_meta if rec['Phase'] in filter_phases and rec['Patient'] in self.pts and rec['Circadian'] in self.circ]
    
    def y_c_pair(self,rec_list):
        scale_lookup = self.CFrame.clin_dict
        
        #self.data = [(rec,scale_lookup[pt][phase]['nHDRS'] for rec in rec_list if rec]
        
    ''' Plot things we care about when it comes to how many recordings each patient x phase has, etc.'''
    def rec_set_size(self):
        filter_phases = dbo.Phase_List(exprs='ephys')
        accounting = np.zeros((len(self.pts),len(filter_phases)))
        detailed_dict = nestdict()

        for pp,pt in enumerate(self.pts):
            
            print(pt + ' has ' + str(len([rec for rec in self.YFrame.file_meta if rec['Phase'] in filter_phases and rec['Patient'] == pt])) + ' recordings')
            for ph,phase in enumerate(filter_phases):
                
                detailed_dict[pt][phase] = [rec for rec in self.YFrame.file_meta if rec['Phase'] == phase and rec['Patient'] == pt]
                print(pt + ' has ' + str(len([rec for rec in self.YFrame.file_meta if rec['Phase'] == phase and rec['Patient'] == pt])) + ' recordings in Phase ' + phase)
                
                accounting[pp,ph] = len(detailed_dict[pt][phase])
        
        #Plot the accounting
        plt.figure()
        plt.imshow(accounting)
        plt.figure()
        plt.plot(accounting[:,:].T)
        
    def plot_psds(self,upper_lim=10):
        plt.figure()
        for ii in range(upper_lim):
            plt.subplot(121)
            plt.plot(np.linspace(0,211,513),np.log10(self.train_set[ii]['Data']['Left']))
            plt.subplot(122)
            plt.plot(np.log10(self.train_set[ii]['Data']['Right']))
        
    def split_validation_set(self,train_ratio=0.6):
        self.train_set, self.validation_set = train_test_split(self.active_rec_list,train_size=train_ratio,shuffle=True)
    
    ''' Calculate oscillatory states for a set of recordings'''
    def calculate_oscillatory_states(self,data_set):
        state_vector = []
        for rr in data_set:
            psd_poly_done = {ch: dbo.poly_subtrLFP(fvect=self.fvect,inp_psd=rr['Data'][ch],polyord=5)[0] for ch in rr['Data'].keys()}
            
            feat_vect = np.zeros(shape=(len(dbo.feat_order),self.ch_num))
            for ff,featname in enumerate(dbo.feat_order):
                dofunc = dbo.feat_dict[featname]
                try: feat_calc = dofunc['fn'](psd_poly_done,self.fvect,dofunc['param'])
                except Exception as e: print(e);pdb.set_trace()
                feat_vect[ff,:] = np.array([feat_calc[ch] for ch in ['Left','Right']])
                
            state_vector.append(feat_vect)
            
            # now we need to get a vector of the clinical states
            
        return np.array(state_vector)
        
    def vectorize_set(self,dataset):
        #First, we're going to go through all the recordings and pull out the patients we want
        # Now, go through each recording and extract the features we want.
        self.feature_extract() #This function works directly with the list and data to populate the feature/state vectors
        self.clin_extract()
        
        #so, basically.... now we have every recording with its associated c_score....
        #The output of this should be a (Feats x Channs) x NRECS matrix
    
    def OBSfeature_extract(self):
        print('Extracting Oscillatory Features')

        big_list = self.YFrame.file_meta
        fvect = self.YFrame.data_basis['F']
        for rr in big_list:
            feat_dict = {key:[] for key in dbo.feat_dict.keys()}
            for featname,dofunc in dbo.feat_dict.items():
                #pdb.set_trace()
                #Choose the zero index of the poly_subtr return because we are not messing with the polynomial vector itself
                #rr['data'] is NOT log10 transformed, so makes no sense to do the poly subtr
                datacontainer = {ch: poly_subtr(fvect=self.fvect,inp_psd=rr['Data'][ch],polyord=5)[0] for ch in rr['Data'].keys()} #THIS RETURNS a PSD, un-log transformed
                
                feat_dict[featname] = dofunc['fn'](datacontainer,self.YFrame.data_basis['F'],dofunc['param'])
            rr.update({'FeatVect':feat_dict})

    def default_run(self):
        #This method captures the default run used to generate the figures from the paper
        self.split_validation_set()
        self.filter_recs(dataset=self.train_set)
        #Here we want to vectorize our training set
        self.vectorize_set(self.train_set)

        #Once our training set is vectorized, we want to train the model
        self.train_model(self.train_set)

        self.validate_model(self.valid_set) #Gives us accuracy measurements
        self.validate_controller(policy='BINARY')
    
    def train_model(self):
        #By this point, we assume we have an active matrix of the data from the recordings of interest
        #We have a feature vector consisting of recordings and clinicals
        #N obs x F feats + N obs of clin
        
       self.learned_model = self.regress_function(X,Y)
    
    def validate_model(self):
        
        summary_stats = self.prediction_accuracy(X_v,Y_v)
        
    def validate_controller(self):
        pass

    def regress_function(self,X,Y):
        pass
    
    
    def readout(self,X):
        print(X.shape[1])
        
class controller_analysis:
    def __init__(self,decoder):
        self.decoder_model = decoder

