#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:54:07 2018

@author: virati
Rewrite of spot_check using only the proper DBS_Osc
"""

#import sys
#sys.path.append('/home/virati/Dropbox/projects/Research/MDD-DBS/Ephys/DBSpace')
import DBS_Osc as dbo

import matplotlib
import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
import scipy.io as io
from collections import defaultdict
from DBS_Osc import nestdict

import pdb

from tkinter.filedialog import askopenfilename
import tkinter as tk

import seaborn as sns

flist = []

def gui_file_select():
    notdone = True
    flist = []
    
    while notdone:
        fname = askopenfilename(initialdir=curr_dir)
        if fname == None or fname == '':
            notdone = False
        else:
            flist.append(fname)
            curr_dir = '/'.join(fname.split('/')[:-1])

    return flist

def spot_check(fname):
    #Here, we just load in the file
    file_contain = dbo.load_BR_dict(fname)
    
    #this is where we decide what to do with the file that constitutes the actual "spot check"
    # from 2013-2018, at least, we sat comfortably in the time-frequency-domain with all of our recording spot checking
    

if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    
    plt.ion()
    
    curr_dir = '/home/virati/MDD_Data/BR/'
    
    results = defaultdict(dict)
    
    notdone = True
    
    #this routine selects the files
    if flist == []:
        flist = gui_file_select()
    
    
    #Now we go through and actually display the files
    for ff,fname in enumerate(flist):
        print(fname)
    