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
from sklearn.metrics import precision_recall_curve, average_precision_score, auc, mean_squared_error, mean_absolute_error
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

import itertools

import time
import copy
import pdb

#%%            
class base_decoder:
    #Parent readout class
    # Some very fixed constants up here
    circ = 'day'
    ch_num = 2
    
    def __init__(self,*args,**kwargs):#BRFrame,ClinFrame,pts,clin_measure='HDRS17'):
        self.YFrame = kwargs['BRFrame']
        self.CFrame = kwargs['ClinFrame']
        self.pts = kwargs['pts']
        self.c_meas = kwargs['clin_measure']
        self.fvect = self.YFrame.data_basis['F']
        
        self.regression_algo = linear_model.LinearRegression
    
    '''Filter out the recordings we want'''
    def filter_recs(self,rec_class='main_study'):
        if rec_class == 'main_study':
            filter_phases = dbo.Phase_List(exprs='ephys')
            self.active_rec_list = [rec for rec in self.YFrame.file_meta if rec['Phase'] in filter_phases and rec['Patient'] in self.pts and rec['Circadian'] in self.circ]
    
        self.filter_phases = filter_phases
    
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
    
    '''Plot PSDs for the first N recordings, sanity check'''
    def plot_psds(self,upper_lim=10):
        plt.figure()
        for ii in range(upper_lim):
            plt.subplot(121)
            plt.plot(np.linspace(0,211,513),np.log10(self.train_set[ii]['Data']['Left']))
            plt.subplot(122)
            plt.plot(np.log10(self.train_set[ii]['Data']['Right']))
    
    '''split out our training and validation set recordings'''
    def split_train_set(self,train_ratio=0.6):
        self.train_set, self.test_set = train_test_split(self.active_rec_list,train_size=train_ratio,shuffle=True)

    '''Setup our data for training'''
    def train_setup(self):
        self.train_set_y, self.train_set_c = self.calculate_states_in_set(self.train_set)

    ''' Train our model'''
    def train_model(self,do_null=False):
        if do_null:
            shuffled_c = copy.deepcopy(self.train_set_c)
            np.random.shuffle(shuffled_c)
            #pdb.set_trace()
            self.decode_model = self.regression_algo().fit(self.train_set_y,shuffled_c)    
        else:
            self.decode_model = self.regression_algo().fit(self.train_set_y,self.train_set_c)
        
    def test_setup(self):
        self.test_set_y, self.test_set_c = self.calculate_states_in_set(self.test_set)

    def test_model(self):
        predicted_c = self.decode_model.predict(self.test_set_y)
        test_stats = self.get_test_stats(self.test_set_y,self.test_set_c,predicted_c)

        plt.figure()
        plt.plot(self.test_set_c,predicted_c,'r.')
        return predicted_c, test_stats
    
    def get_test_stats(self,test_y,true_c,predicted_c):        
        #Pearson
        test_y = test_y.squeeze()
        true_c = true_c.squeeze()
        predicted_c = predicted_c.squeeze()
        p_stats = stats.pearsonr(true_c,predicted_c)
        #Spearman
        s_stats = stats.spearmanr(true_c,predicted_c)
        
        #Linear Regression
        regr_model = linear_model.LinearRegression().fit(true_c.reshape(-1,1),predicted_c.reshape(-1,1))
        lr_slope = regr_model.coef_[0]
        
        stat_dict = {'Score':self.decode_model.score(test_y,true_c),'Pearson':p_stats,'Spearman':s_stats,'Slope':lr_slope}
        return stat_dict
    
    '''Plot the test statistics'''
    def plot_test_stats(self):
        predicted_c = self.test_model()
        
        plt.scatter(self.test_set_c,predicted_c)
        
        #except Exception as e: print(e); pdb.set_trace()
        print(corr)
        
    '''Plot the regression visualization of the test procedure'''
    def plot_test_regression(self):
        #do a final test on *all* the data for plotting purposes
        predicted_c = self.decode_model.predict(self.test_set_y)
        r2score = self.decode_model.score(self.test_set_y,self.test_set_c)
        mse = mean_squared_error(self.test_set_c,predicted_c)
        corr = stats.pearsonr(self.test_set_c.squeeze(),predicted_c.squeeze())
        
        plt.plot([0,1],[0,1],color='gray',linestyle='dotted')
        ax = sns.regplot(x=self.test_set_c,y=predicted_c)
        plt.title('R^2:' + str(r2score) + '\n' + ' MSE:' + str(mse) + '\n Corr:' + str(corr))
        plt.xlim((0,1.1))
        plt.ylim((0,1.1))
        
    ''' Calculate oscillatory states for a set of recordings'''
    def calculate_states_in_set(self,data_set):
        state_vector = []
        depr_vector = []
        
        for rr in data_set:
            psd_poly_done = {ch: dbo.poly_subtrLFP(fvect=self.fvect,inp_psd=rr['Data'][ch],polyord=5)[0] for ch in rr['Data'].keys()}
            
            feat_vect = np.zeros(shape=(len(dbo.feat_order),self.ch_num))
            for ff,featname in enumerate(dbo.feat_order):
                dofunc = dbo.feat_dict[featname]
                feat_calc = dofunc['fn'](psd_poly_done,self.fvect,dofunc['param'])
                #except Exception as e: print(e);pdb.set_trace()
                feat_vect[ff,:] = np.array([feat_calc[ch] for ch in ['Left','Right']])
                
            # We need to flatten the state between channels...
            #Then we go ahead and append it to the state vector
            state_vector.append(np.reshape(feat_vect,-1,order='F')) #we want our FEATURE index to change quickest so we go (0,0) -> (1,0) -> (2,0) -> ... (4,1)
            
            # now we need to get a vector of the clinical states
            depr_value = self.CFrame.get_depression_measure('DBS'+rr['Patient'],self.c_meas,rr['Phase'])
            depr_vector.append(depr_value)
            
        return np.array(state_vector), np.array(depr_vector)
    
    def get_coeffs(self):
        model = self.decode_model
        active_coeffs = self.decode_model.coef_
        
        return active_coeffs
    '''Plot coefficients of our model'''
    def plot_decode_coeffs(self):
        active_coeffs = np.array(self.get_coeffs).squeeze()
        #plt.subplot(1,2,side+1)
        plt.plot(active_coeffs)
        plt.hlines(0,-2,10,linestyle='dotted')
        plt.vlines(5,-1,1,linestyle='solid',color='blue')
        plt.ylim((np.max(np.abs(active_coeffs)) + 0.01,-np.max(np.abs(active_coeffs)) - 0.01))
        plt.xlim((-1,10))
            
    def plot_test_ensemble(self):
        plt.figure()
        plt.subplot(211)
        self.plot_test_regression()
        plt.subplot(212)
        self.plot_decode_coeffs()

class weekly_decoder(base_decoder):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        
        if kwargs['algo'] == 'ENR':
            self.regression_algo = linear_model.ElasticNet(alpha=0.2)
        elif kwargs['algo'] == 'Ridge':
            self.regression_algo = linear_model.RidgeCV()
        elif kwargs['algo'] == 'Lasso':
            self.regression_algo = linear_model.LassoCV()
    
    def train_model(self):
        self.decode_model = self.regression_algo.fit(self.train_set_y,self.train_set_c)
        print('Alpha: ' + str(self.decode_model.alpha_) + ' | L1r: ' + str(self.decode_model.l1_ratio_))
        
        #self.plot_decode_coeffs(self.decode_model)
        
    def aggregate_weeks(self,dataset):
        #print('Performing Training Setup for Weekly Decoder')
        #go through our training set and aggregate every recording within a given week
        #train_set_y,train_set_c = self.calculate_states_in_set(self.train_set)
        
        running_list = []
        for pt in self.pts:
            for phase in self.filter_phases:
                block_set = [rr for rr in dataset if rr['Patient'] == pt and rr['Phase'] == phase]
                if block_set != []:
                    y_set,c_set = self.calculate_states_in_set(block_set)
                    weekly_y_set = np.mean(y_set,axis=0)
                    
                    running_list.append((weekly_y_set,c_set[0],pt,phase)) #all the c_set values should be the exact same
        
        y_state = np.array([a for (a,b,c,d) in running_list]) #outputs ~168 observed weeks x 10 features
        c_state = np.array([b for (a,b,c,d) in running_list]).reshape(-1,1) #outputs ~168 observed weeks
        pt_name = np.array([c for (a,b,c,d) in running_list])
        phase_label = np.array([d for (a,b,c,d) in running_list])
        
        return y_state, c_state, pt_name, phase_label
    
    def train_setup(self):
        print('Performing Training Setup for Weekly Decoder')
    
        self.train_set_y, self.train_set_c, self.train_set_pt, self.train_set_ph  = self.aggregate_weeks(self.train_set)
    def test_setup(self):
        print('Performing TESTING Setup for Weekly Decoder')
        
        self.test_set_y, self.test_set_c, self.train_set_pt, self.train_set_ph = self.aggregate_weeks(self.test_set)

class weekly_decoderCV(weekly_decoder):
    def __init__(self,*args,**kwargs):
        print('Initialized the Weekly CV decoder')
        super().__init__(*args,**kwargs)
        
        if kwargs['algo'] == 'ENR':
            self.regression_algo = linear_model.ElasticNetCV
            self.model_args = {'alphas':np.linspace(0.01,0.04,20),'l1_ratio':np.linspace(0.1,0.3,10),'cv':10}
            
        self.pt_CV_sets(n=3)

        
    def pt_CV_sets(self,n=3):
        pt_combos = list(itertools.combinations(self.pts,n))
        
        self.CV_num_combos = len(pt_combos)
        self.CV_pt_combos = pt_combos
        
    def train_setup(self):
        print('Performing Training Setup for Weekly Decoder')
    
        self.train_set_y, self.train_set_c, self.train_set_pt  = self.aggregate_weeks(self.train_set)
    
    ''' Train our model'''
    def train_model(self):
        #Our first goal is to learn a model for each patient combination
        decode_model_combos = [None] * self.CV_num_combos
        model_performance_combos = [None] * self.CV_num_combos
        for run,pt_combo in enumerate(self.CV_pt_combos):
            print(pt_combo)
            combo_train_y = [a for (a,c) in zip(self.train_set_y,self.train_set_pt) if c in pt_combo]
            combo_train_c = [b for (b,c) in zip(self.train_set_c,self.train_set_pt) if c in pt_combo]
            
            decode_model_combos[run] = self.regression_algo(**self.model_args).fit(combo_train_y,combo_train_c)
            
            combo_test_y = [a for (a,c) in zip(self.train_set_y,self.train_set_pt) if c not in pt_combo]
            combo_test_c = [b for (b,c) in zip(self.train_set_c,self.train_set_pt) if c not in pt_combo]
            
            #model_performance_combos[run] = decode_model_combos[run].score(combo_test_y,combo_test_c)
            #pred_c = decode_model_combos[run].predict(combo_test_y)
            #model_performance_combos[run] = mean_absolute_error(combo_test_c,pred_c)
            
        
        self.decode_model_combos_ = decode_model_combos
        
        average_model_coeffs,_ = self.get_average_model(self.decode_model_combos_)
        self.decode_model = linear_model.LinearRegression()
        self.decode_model.coef_ = average_model_coeffs
        self.decode_model.intercept_ = np.mean([m.intercept_ for m in self.decode_model_combos_])
        
    def get_average_model(self,model):
        active_coeffs = []
        for ii in self.decode_model_combos_:
            active_coeffs.append([ii.coef_])
        
        active_coeffs = np.array(active_coeffs).squeeze()
        average_model = np.mean(active_coeffs,axis=0)
        #average_model = np.zeros(shape=active_coeffs.shape)
        
        #do some stats
        
        
        #return the average model with the stats for each coefficient
        return average_model, stats
    
    def test_model(self):
        ensemble_score = []
        ensemble_corr = []
        self.test_stats = [] #{'Prediction Score': [], 'Pearson Corr Score': [], 'Spearman Corr Score': []}

        for tt in range(1000):
            test_subset_y,test_subset_c = zip(*random.sample(list(zip(self.test_set_y,self.test_set_c)),np.ceil(0.8 * len(self.test_set_y)).astype(np.int)))
            test_subset_y = np.array(test_subset_y)
            test_subset_c = np.array(test_subset_c)
            
            predicted_c = self.decode_model.predict(test_subset_y)
            self.test_stats.append(self.get_test_stats(test_subset_y,test_subset_c,predicted_c))
    
    def plot_test_regression_figure(self):
        #do a final test on *all* the data for plotting purposes
        predicted_c = self.decode_model.predict(self.test_set_y)
        r2score = self.decode_model.score(self.test_set_y,self.test_set_c)
        mse = mean_squared_error(self.test_set_c,predicted_c)
        plt.figure()
        plt.plot([0,1],[0,1],color='gray',linestyle='dotted')
        ax = sns.regplot(x=self.test_set_c,y=predicted_c)
        plt.title('R2:' + str(r2score) + '\n' + ' MSE:' + str(mse))
        plt.xlim((0,1.1))
        plt.ylim((0,1.1))
        
            
    def plot_test_stats(self):
        plt.figure()
        plt.subplot(311)
        #plt.scatter(test_subset_c,predicted_c)
        plt.hist([a['Score'] for a in self.test_stats]);plt.title('R2 Score')
        plt.subplot(312)
        plt.hist([a['Pearson'][0] for a in self.test_stats]);plt.title('Pearson')
        plt.subplot(313)
        plt.hist([a['Spearman'][0] for a in self.test_stats]);plt.title('Spearman')
            
            
    '''PLOTTING--------------------------------------------------------'''

    '''Plot the decoding CV coefficients'''
    def plot_decode_CV(self):
        plt.figure()
        for side in range(2):
            active_coeffs = []
            for ii in self.decode_model_combos_:
                active_coeffs.append([ii.coef_[side*5:side*5+5]])
            
            active_coeffs = np.array(active_coeffs).squeeze()
            plt.subplot(1,2,side+1)
            #plt.plot(ii.coef_[side*5:side*5+5])
            #pdb.set_trace()
            vp_obj = sns.violinplot(data=active_coeffs,scale='width')
            plt.setp(vp_obj.collections,alpha=0.3)
        
            average_model, _ = self.get_average_model(self.decode_model_combos_)
            plt.plot(average_model[side*5:side*5+5])
            plt.hlines(0,-2,7,linestyle='dotted')
            plt.ylim((-0.3,0.3))
            plt.xlim((-1,5))
            
        
class controller_analysis:
    def __init__(self,decoder):
        self.decoder_model = decoder
        # get our binarized disease states
        
    def roc_auc(self):
        pass
