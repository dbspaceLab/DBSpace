"""
object-oriented version of spot_check original

"""
from tkinter import ttk
from tkinter import *
import tkinter
import os
from tkinter import filedialog

# libraries needed to run spot_check
import sys
import DBS_Osc
import matplotlib
import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
import scipy.io as io
from collections import defaultdict

class GUI(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.grid()
        self.mainWindow = tkinter.Tk()
        self.mainWindow.title("Select LFP File")
        self.mainWindow.geometry('640x480-8-200')

        self.rightFrame = tkinter.Frame(self.mainWindow)
        self.rightFrame.grid(row=1, column=2, sticky='n')
        self.rightFrame.columnconfigure(0,weight=1)

        ttk.Button(self.rightFrame, text = "Select Directory", command=self.askdirectory).grid(row=0, column=0, padx=5, pady=5)

        # configure the columns
        self.mainWindow.columnconfigure(0, weight=1)
        self.mainWindow.columnconfigure(1, weight=1)
        self.mainWindow.grid_columnconfigure(2, weight=1)

        self.fileList = tkinter.Listbox(self.mainWindow)
        self.fileList.grid(row=1, column=0, sticky='nsew', rowspan=2)
        self.fileList.config(border=2, relief='sunken')

        self.listScroll = tkinter.Scrollbar(self.mainWindow, orient=tkinter.VERTICAL, command=self.fileList.yview)
        self.listScroll.grid(row=1, column=1, sticky='nsw', rowspan=2)
        self.fileList['yscrollcommand'] = self.listScroll.set

    # User selects directory with files to be analyzed
    def askdirectory(self):
        self.clear()
        self.directory = filedialog.askdirectory()
        print(self.directory)

        try:
            for zone in os.listdir(self.directory):
                self.fileList.insert(tkinter.END, zone)
        except FileNotFoundError:
            print('Try connecting to Server and try again')

    # Clear contents of the file list
    def clear(self):
        self.fileList.delete(0, 'end')


if __name__ == "__main__":
    guiFrame = GUI()
    guiFrame.mainloop()


#
# json file for major experiments
experiments = ['Targeting','Amplitude','Frequency','Resting']
#
#
# flist = []
#
# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 20}
#
# matplotlib.rc('font', **font)
# matplotlib.rcParams['svg.fonttype'] = 'none'
#
# plt.rcParams['image.cmap'] = 'jet'
#
#
# def spot_check(fname):
#
#     Container = DBS_Osc.load_BR_feats(fname,snippet=False)
#     plt.figure()
#     plt.subplot(2,1,1)
#     plt.plot(Container['TS']['T'],Container['TS']['Y'][:,0])
#     plt.ylim((-0.001,0.005))
#     NFFT = 2**10
#     SG = defaultdict(dict)
#     for cc in range(2):
#
#         F,T,SG[cc] = sig.spectrogram(Container['TS']['Y'][:,cc],nperseg=NFFT,noverlap=NFFT*0.5,window=sig.get_window('blackmanharris',NFFT),fs=422)
#         plt.subplot(2,2,3+cc)
#         plt.pcolormesh(T,F,10*np.log10(SG[cc]),rasterized=True)
#         plt.clim((-200,-100))
#     plt.colorbar()
#     plt.suptitle('Raw TS: ' + fname.split('/')[-1])
#     return {'TS':Container['TS']['Y'],'TF':{'SG':SG,'F':F,'T':T}}
#
#
# results = defaultdict(dict)
#
# for fname in flist:
#     results[fname] = spot_check(fname)
#
#
# def grab_median(TFcont,tlim=(880,900)):
#     # Plot some PSDs
#     plt.figure()
#
#     for cc in range(2):
#         plt.subplot(1,2,cc+1)
#         T = TFcont['TF']['T']
#         F = TFcont['TF']['F']
#         SG = TFcont['TF']['SG']
#         t_idxs = np.where(np.logical_and(T > tlim[0], T < tlim[1]))
#         med_psd = np.median(10*np.log10(SG[cc][:,t_idxs]).squeeze(),axis=1)
#         plt.plot(F,med_psd)
#
#
# for key,val in results.items():
#     print(key)
#     grab_median(val,tlim=(0,100))
