#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 11:23:06 2017

@author: virati
This library is a small quick library for 3d plotting of EEG
"""

import matplotlib.pyplot as plt
import mne
import numpy as np
import scipy.stats as stats

import pdb

import time
import pylab

from mpl_toolkits.mplot3d import Axes3D

import mayavi.mlab as mlab
from mayavi.mlab import *


def return_adj_net(dist_thresh = 3):
    egipos = mne.channels.read_montage('/tmp/GSN-HydroCel-257.sfp')
    etrodes = egipos.pos
    
    dist = np.zeros((257,257))
    for ii,ipos in enumerate(etrodes):
        #loop through all others and find the distances
        for jj,jpos in enumerate(etrodes):
            dist[ii][jj] = np.linalg.norm(ipos - jpos)
            
    mask = (dist <= dist_thresh).astype(int)
    
    return mask

def get_coords(scale,montage='dense'):
    if montage == 'dense':
        fname = '/home/virati/Dropbox/GSN-HydroCel-257.sfp'
    elif montage == 'standard':
        fname = '/home/virati/Dropbox/standard_postfixed.elc'
    
    egipos = mne.channels.read_montage(fname)
    etrodes = scale*egipos.pos
    
    etrodes[:,2] = etrodes[:,2]
    
    return etrodes
    

#This function is to plot vector data for each channel at the channel's coordinates
def plot_3d_locs(band,ax,n=1,scale=1,clims=(0,0),label='generic',animate=False,unwrap=False,sparse_labels = True,highlight=[],montage='dense'):
    #fig = plt.figure()
    
    etrodes = get_coords(scale=scale)

    #gotta normalize the color
    #band = np.tanh(band / 10) #5dB seems to be reasonable
    
    cm = plt.cm.get_cmap('jet')
    
    if clims == (0,0):
        clims = (np.min(band),np.max(band))
    
    linewidths = np.ones_like(etrodes[:,0])
    linewidths[highlight] = 5
    sc = ax.scatter(etrodes[:,0],etrodes[:,1],etrodes[:,2],color='#ffffff',s=100,linewidth=3,alpha=0.1,edgecolors='k')
 
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # Get rid of the spines                         
    ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
    ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    #ax.xlim((-10,10))
    ax.set_xticks([])     
    ax.set_yticks([])     
    ax.set_zticks([])
    
    ims = []
    plt.title(label)
    
    print('Animation: ' + str(animate))
    if animate:
        
        for angl in range(0,360,10):
            print('Animating frame ' + str(angl))
            ax.view_init(azim=angl)
            strangl = '000' + str(angl)
            plt.savefig('/tmp/'+ label + '_' + strangl[-3:] + '.png')
            time.sleep(.3)

# The goal of this is to plot the bands on the scalp
def plot_3d_scalp(band,infig=[],n=1,clims=(0,0),scale=1,label='generic',animate=False,unwrap=False,sparse_labels = True,highlight=[],montage='dense',alpha=1,marker_scale=5,anno_top=True,binary_mask=False):
    #fig = plt.figure()
    
    if montage == 'dense':
        fname = '/home/virati/Dropbox/GSN-HydroCel-257.sfp'
    elif montage == 'standard':
        fname = '/home/virati/Dropbox/standard_postfixed.elc'
    
    egipos = mne.channels.read_montage(fname)
    etrodes = scale * egipos.pos
    
    #gotta normalize the color
    #band = np.tanh(band / 10) #5dB seems to be reasonable
    
    cm = plt.cm.get_cmap('jet')
    
    if clims == (0,0):
        clims = (np.min(band),np.max(band))
    
    if unwrap:
        flat_etrodes = np.copy(etrodes)
        flat_etrodes[:,2] = flat_etrodes[:,2] - np.max(flat_etrodes[:,2]) + 0.01
    
        flat_etrodes[:,0] = flat_etrodes[:,0] * -10*(flat_etrodes[:,2] + 3*1/(flat_etrodes[:,2] - 0.6) + 0.5)
        flat_etrodes[:,1] = flat_etrodes[:,1] * -10*(flat_etrodes[:,2] + 3*1/(flat_etrodes[:,2] - 0.6) + 0.5)
        
        if infig == []:
            fig=plt.figure()
            ax = fig.add_subplot(1,1,n)
        else:
            ax = infig
        
        linewidths = marker_scale * np.ones_like(flat_etrodes[:,0])
        linewidths[highlight] = 3
        #below changes can be: linewidth to only do the highlights, or fixed at 2 or something
        sc = plt.scatter(flat_etrodes[:,0],flat_etrodes[:,1],c=band,vmin=clims[0],vmax=clims[1],s=300,cmap=cm,alpha=alpha,linewidth=linewidths,marker='o')
        #this adds x's over the highlights
        #plt.scatter(flat_etrodes[:,0],flat_etrodes[:,1],c=None,vmin=clims[0],vmax=clims[1],s=300,cmap=cm,alpha=1,linewidth=linewidths,marker='x')
        
        #Which channels are above two stds?
        zsc_band = stats.zscore(band)
        top_etrodes = np.where(np.abs(zsc_band) > 1)[0]
        
        if sparse_labels:
            annotate_list = top_etrodes
        else:
            annotate_list = range(257)
        
        for ii in annotate_list:
            plt.annotate('E'+str(ii+1),(flat_etrodes[ii,0],flat_etrodes[ii,1]),size=12)
        
        #THIS IS NEW 01/11/2021 see if it works properly
        if anno_top: sc_top = plt.scatter(flat_etrodes[annotate_list,0],flat_etrodes[annotate_list,1],c=band[annotate_list],vmin=clims[0],vmax=clims[1],s=300,cmap=cm,alpha=1.0,linewidth=linewidths,marker='o')
        if binary_mask: 
            ch_ones = np.where(band == 0)[0]
            bin_mask = plt.scatter(flat_etrodes[ch_ones,0],flat_etrodes[ch_ones,1],c='black',vmin=0,vmax=1,s=300,cmap=cm,alpha=1.0,linewidth=linewidths,marker='o')
        
        plt.axis('off')        
        
        plt.colorbar(sc)
        plt.title(label)
        
    else:
        if infig == []:
            fig = plt.figure()
            ax = fig.add_subplot(1,1,n,projection='3d')
        else:
            ax = infig #infig.add_subplot(1,1,n,projection='3d')
            
        linewidths = np.ones_like(etrodes[:,0])
        linewidths[highlight] = 5
        #REMOVED a 10* z component here, I think it was originally added to help visualization
        sc = ax.scatter(etrodes[:,0],etrodes[:,1],etrodes[:,2],c=band,vmin=clims[0],vmax=clims[1],s=300,cmap=cm,linewidth=linewidths,alpha=alpha)
    
        try:plt.colorbar(sc)
        except: pdb.set_trace()
     
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # Get rid of the spines                         
        ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
        ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
        ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        #ax.xlim((-10,10))
        ax.set_xticks([])                               
        ax.set_yticks([])                               
        ax.set_zticks([])
        
        ims = []
        plt.title(label)
        
        print('Animation: ' + str(animate))
        if animate:
            
            for angl in range(0,360,10):
                print('Animating frame ' + str(angl))
                ax.view_init(azim=angl)
                strangl = '000' + str(angl)
                plt.savefig('/tmp/'+ label + '_' + strangl[-3:] + '.png')
                time.sleep(.3)

## DO UNIT TEST HERE

'''
This function takes in COORDINATES and plots dots at those COORDINATES
'''
def plot_coords(band, active_mask=[],rad= [],color=[0.,0.,0.],plot_overlay = True,alpha=0.8):
    if active_mask == []:
        active_mask = np.ones_like(band[:,-1]).astype(np.bool)
        rad_factor = 10
    else:
        band = band[active_mask,:]
        rad_factor = 20
        
    if not rad:
        #rad = np.zeros_like(band[:,2])
        #rad[active_mask] = 20
        rad = rad_factor*np.ones_like(band[:,-1])
        #rad[active_mask] = 0
        #rad = np.random.normal(size=band[:,2].shape)

    
    #figure(bgcolor=(1,1,1))
    nodes = points3d(band[:,0],band[:,1],band[:,2], rad,color=color, scale_factor=10,opacity=alpha)
    nodes.glyph.scale_mode = 'scale_by_vector'
    
    if plot_overlay:
        points3d(band[:,0],band[:,1],band[:,2], rad2,color=(0.,0.,0.),colormap="copper", scale_factor=.4,opacity=alpha/2)

def plot_tracts(band, active_mask=[],rad= [],color=[0.,0.,0.],alpha=1):
    if not rad:
        rad = np.zeros_like(band[:,2])
        rad[active_mask] = 20
        rad2 = 20*np.ones_like(rad)
        rad2[active_mask] = 0
        #rad = np.random.normal(size=band[:,2].shape)
    
    #plot3d(band[:,0],band[:,1],band[:,2],color=color,opacity=0.8)
    points3d(band[:,0],band[:,1],band[:,2], rad2,color=color,colormap="copper", scale_factor=.4,opacity=alpha/2)
    
def maya_band_display(band,montage='dense',label=''):
    if montage == 'dense':
        fname = '/home/virati/Dropbox/GSN-HydroCel-257.sfp'
    elif montage == 'standard':
        fname = '/home/virati/Dropbox/standard_postfixed.elc'
    
    mlab.figure(bgcolor=(1.0,1.0,1.0))
    
    egipos = mne.channels.read_montage(fname)
    etrodes = egipos.pos
    
    
    # Make a single sphere for the head
    head = points3d(0,0,0,scale_factor=15)
    # Setup electrodes as spheres around head
    nodes = points3d(etrodes[:,0],etrodes[:,1],etrodes[:,2], scale_factor=2)
    nodes.glyph.scale_mode = 'scale_by_vector'
    
    #Have to bring band from (-1,1) to (0,1) for mayavi color bullshit
    band_norm = band / np.max(band)
    band_norm += 1/2
    #band_norm = 0.5 * np.tanh(band * 5) + 0.5
    

    #This sets the colors for the nodes themselves to the band  changes after normalization into [0,1]
    nodes.mlab_source.dataset.point_data.scalars = (band_norm)
    #show()
    mlab.title(label)
    
def plot_maya_scalp(band,n=1,clims=(0,0),color=(1.,0.,0.),scale=1,label='generic',animate=False,unwrap=False,sparse_labels = True,highlight=[],montage='dense',alpha=1):
    
    if montage == 'dense':
        fname = '/home/virati/Dropbox/GSN-HydroCel-257.sfp'
    elif montage == 'standard':
        fname = '/home/virati/Dropbox/standard_postfixed.elc'
    
    egipos = mne.channels.read_montage(fname)
    etrodes = scale * egipos.pos
    
    #gotta normalize the color
    #band = np.tanh(band / 10) #5dB seems to be reasonable

    
    if clims == (0,0):
        clims = (np.min(band),np.max(band))
    
    linewidths = np.ones_like(etrodes[:,0])
    linewidths[highlight] = 5
    #REMOVED a 10* z component here, I think it was originally added to help visualization
    
    points3d(etrodes[:,0],etrodes[:,1],etrodes[:,2], 20*band, color=color, scale_factor=.25,opacity=alpha)
    points3d(etrodes[:,0],etrodes[:,1],etrodes[:,2], 20*np.abs(1-band), color=(1.,1.,1.), scale_factor=.25,opacity=alpha)
