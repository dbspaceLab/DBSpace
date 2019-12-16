#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 20:21:24 2018

@author: virati
Mess with IPG waveform
"""
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

import allantools
from allantools.noise import pink as f_noise

plt.close('all')


decay = 100

startfs = 1e6
endfs = 4220

#bvasic load of a 1/10th stim snippet
tenth_sec_stim = np.load('/home/virati/tenth_sec_ipg.npy')



#gaussian filter versus
tframe = (-10,15)
tlen = tframe[1] - tframe[0]
full_stim = np.tile(tenth_sec_stim,10*21)[0:int(startfs*tlen)]
tvect = np.linspace(-10,11,full_stim.shape[0])



if 0:
    expon = sig.exponential(1000,0,500,False)
    conv_out = (np.convolve(full_stim,expon) / 50)
    conv_stim = conv_out
else:
    conv_stim = full_stim

#do some bare-minimum filtering before the sampling
bl,al = sig.butter(2,5e2/startfs,btype='lowpass')
conv_stim = sig.lfilter(bl,al,conv_stim)

bh,ah = sig.butter(2,10/startfs,btype='highpass')
conv_stim = sig.lfilter(bh,ah,conv_stim)

plt.figure()
plt.plot(tvect[0:int(startfs*20)],full_stim[0:int(startfs*20)],label='Full Stim')
plt.plot(tvect[0:int(startfs*20)],conv_stim[0:int(startfs*20)],label = 'Post Conv')

#Now we're sampling the above convolved waveform
skip_ts = 238
meas_stim = conv_stim[0::skip_ts][0:(endfs*20)]
meas_stim[0:endfs*10] = 0
ds_tvect = np.linspace(-10,10,meas_stim.shape[0])

plt.plot(ds_tvect,meas_stim,label='Sampled 4220Hz (Post Conv)')
plt.legend()

np.save('/home/virati/Dropbox/projects/Research/MDD-DBS/Data/StimEphys/stim_wform',meas_stim)

#%%


#%%

nois_meas_stim = meas_stim + 0.01 * np.array(f_noise(endfs*20))

# b,a = sig.butter(5,100/422,btype='lowpass')
# meas_stim = sig.lfilter(b,a,meas_stim)

# b,a = sig.butter(5,1/422,btype='highpass')
# meas_stim = sig.lfilter(b,a,meas_stim)

#meas_stim = np.tanh(meas_stim)


plt.figure()
plt.subplot(311)
plt.plot(ds_tvect,nois_meas_stim)

F,Pxx = sig.welch(nois_meas_stim[endfs*10:],fs=endfs)
plt.subplot(312)
plt.plot(F,10*np.log10(Pxx))

plt.subplot(313)
F,T,SG = sig.spectrogram(nois_meas_stim,nperseg=2**8,noverlap=(2**8) - 10,window=sig.get_window('blackmanharris',2**8),fs=endfs)
plt.pcolormesh(T,F,10*np.log10(SG))

#full_stim = gaussian_filter1d(full_stim,100)
#237 gets us up to around 21 seconds at 1Mhz

#finally, need to highpass filter this
#b,a = sig.butter(order,Wn,btype='highpass',analog=True)
#w,h = sig.freqz(b,a)
#plt.figure()
#plt.plot(w,20*np.log10(abs(h)))

#full_stim = sig.filtfilt(b,a,full_stim)

#%%

