#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 15:28:10 2018

@author: virati
BR Data Library Script

The PURPOSE of this library should be to just bring in the BrainRadio data in a format that the DSV can handle
For example: Determining which phase a recording belongs to will NOT be done in this script, that is under the perview of the DSV

"""
import datetime
import glob
import json
import logging
import os
import pickle
import random
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from dbspace.readout.ClinVect import Phase_List
from dbspace.signal.oscillations import FEAT_DICT, gen_psd
from dbspace.utils.io.pcs import load_br_file

from dbspace.utils.functions import nearest


class BR_Data_Tree:
    """
    This class generates the data frame (.pickle) that is used in the SCC_Readout project

    """

    sec_end = 10
    intermediate_dir = "../../assets/intermediate_data"

    def __init__(
        self,
        frame_label=None,
        premade_frame_file=None,
        clin_vector_file=None,
        do_pts=["901", "903", "905", "906", "907", "908"],
    ):
        # Fix this and don't really make it accessible; we'll stick with a single intermediate file unless you really want to change it

        self.do_pts = do_pts
        self.fs = 422

        # Load in our clinical vector object with the data from ClinVec.json
        if clin_vector_file is None:
            raise ValueError("Need to pass in a Clinical Vector (json) file...")

        CVect = json.load(open(clin_vector_file))["HAMDs"]
        self.ClinVect = {pt["pt"]: pt for pt in CVect}

        if frame_label is None:
            # construct from date
            self._frame_label = str(datetime.date)
        if Path.Path(
            self.intermediate_dir + "/ChronicFrame_" + self._frame_label
        ).is_file():
            logging.info(
                f"Loading in {frame_label} frame from intermediate file cache..."
            )
        else:
            # make a new frame
            pass

        self.data_basis = defaultdict()

        if premade_frame_file is None:
            print("Generating the dataframe...")
            self.generate_sequence()
            # Save it now
            self.Save_Frame(name_addendum=frame_label)
            # Now just dump us out so we can do whatever we need to with the file above

        else:
            # IF we're loading in a preframe, we're probably doing a bigger analysis
            self.preFrame_file = premade_frame_file
            print("Loading in PreFrame..." + self.preFrame_file)
            self.Import_Frame(self.preFrame_file)

        # how many seconds to take from the chronic recordings

    def generate_TD_sequence(self):
        self.build_phase_dict()

        # Here we go through all of our files in the dictionaries and put them into our database
        self.list_files()
        self.meta_files()

        # Load in our data (timeseries)
        self.Load_Data(domain="T")
        self.prune_meta()
        self.check_empty_phases()

    def generate_sequence(self, domain="F"):
        # Here we build the dictionary with all of our phases and files
        self.build_phase_dict()

        # Here we go through all of our files in the dictionaries and put them into our database
        self.list_files()
        self.meta_files()

        # Load in our data (timeseries)
        self.Load_Data(domain=domain)
        # now go in and remove anything with a bad flag
        self.Remove_BadFlags()

        # Go in and compute recording fidelity measures for all recordings
        self.Check_GC()

        # take out the phases that don't exist, and any other stuff, but so far that's all this does
        self.prune_meta()
        self.check_empty_phases()

        # in case the meta-data isn't properly updated from the loaded in deta
        print("Data Loaded")

    def Check_GC(self):
        # do just the key features
        # get the stim-related and gc related measures
        print("Checking for Gain Compression...")
        for rr in self.file_meta:
            gc_measures = ["Stim", "SHarm", "THarm", "fSlope", "nFloor"]
            gc_results = {key: 0 for key in gc_measures}
            for meas in gc_measures:
                dofunc = FEAT_DICT[meas]
                gc_results[meas] = dofunc["fn"](
                    rr["Data"], self.data_basis["F"], dofunc["param"]
                )

            # let's do some logic to find out if GC is happening
            isgc = (
                gc_results["nFloor"]["Left"] < -8 or gc_results["nFloor"]["Right"] < -8
            ) and (
                gc_results["SHarm"]["Left"] / gc_results["THarm"]["Left"] < 1
                or gc_results["SHarm"]["Right"] / gc_results["THarm"]["Right"] < 1
            )

            # check if stim is on
            isstim = (
                gc_results["Stim"]["Left"] > 0.0001
                or gc_results["Stim"]["Right"] > 0.0001
            )

            rr.update({"GC_Flag": {"Flag": isgc, "Raw": gc_results, "Stim": isstim}})

    def plot_GC_distribution(self):
        gc_plot = [None] * len(self.file_meta)

        for rr, rec in enumerate(self.file_meta):
            gc_plot[rr] = {
                side: rec["GC_Flag"]["Raw"]["SHarm"][side]
                / rec["GC_Flag"]["Raw"]["THarm"][side]
                for side in ["Left", "Right"]
            }

        plt.figure()
        plt.plot(gc_plot)

    def grab_recs(self, reqs):
        fullfilt_data = [
            (rr["Data"]["Left"], rr["Data"]["Right"], rr["Phase"], rr["Patient"])
            for rr in self.file_meta
            if rr["Patient"] in reqs["Patient"]
        ]
        return fullfilt_data

    """Remove files that are Flagged to be Bad"""

    def Remove_BadFlags(self):
        try:
            self.file_meta = [rr for rr in self.file_meta if rr["BadFlag"] != True]
        except:
            raise ValueError("Error in the bad flag parsing...")

    """This method parses the root data structure and populates a list of recordings"""

    def list_files(self):
        file_list = []
        for pt in self.do_pts:
            for filename in glob.iglob(
                self.data_root_dir + "/" + pt + "/**/" + "*.txt", recursive=True
            ):
                # Append the full path to a list
                # check the file's STRUCTURE HERE

                islogf = filename[-7:] == "LOG.txt"
                isrealtf = filename[-6:] == "RT.txt"
                iseepromf = filename[-9:] == "Table.txt"

                if not (islogf or isrealtf or iseepromf):
                    file_list.append(filename)
        self.file_list = file_list

    # Extract is referring to taking information from the raw BR files
    def extract_date(self, fname):
        datestr = fname.split("_")
        # example: '/home/virati/MDD' 'Data/BR/908/Session' '2016' '08' '11' 'Thursday/DBS908' '2016' '08' '10' '17' '20' '28' 'MR' '14.txt'
        # important for date is -8 (year) -> month -> day -> hour -> minute -> second -> ... -> recording number
        return datetime.datetime.strptime(
            datestr[-8] + "/" + datestr[-7] + "/" + datestr[-9], "%m/%d/%Y"
        )

    # Get for computing time based off of file information
    def get_time(self, fname):
        datestr = fname.split("_")

        day_bound = [
            datetime.datetime.strptime("10:00", "%H:%M"),
            datetime.datetime.strptime("21:00", "%H:%M"),
        ]

        # where is the time?
        in_time = datetime.datetime.strptime(datestr[-6] + ":" + datestr[-5], "%H:%M")

        if in_time < day_bound[0] or in_time > day_bound[1]:
            return "night"
        else:
            return "day"

    # This function is tasked with returning the type of recording; right now we just care whether it's CHRONIC
    def extract_rectype(self, fname):
        filesiz = os.path.getsize(fname)

        ftype = None
        if filesiz > 1e5 and filesiz < 1e7:
            ftype = "Chronic"
        elif filesiz > 1e7:
            ftype = "Dense"
            # need more logic here to check what type of experiment it actually is
        return ftype

    # Extract the recording gain for every file
    def extract_gains(self, fname):
        xml_fname = fname.split(".")[0] + ".xml"

    def extract_pt(self, fname):
        return fname.split("BR")[1][1:4]

    def build_phase_dict(self):
        # In this method we go in and map the date of our sessions to the phase that the patient was in
        # TODO a checkpoint here so we can see the mapping between weeks/dates and phases and validate

        # CONVENTION DECIDED: if a recording is taken between week C04 and C05 -> it belongs to C05 since the clinical questionaiires ask about the LAST 7 days

        phdate_dict = defaultdict(dict)

        for pt in self.do_pts:
            alv = self.ClinVect["DBS" + pt]
            phdate_dict[pt] = defaultdict(dict)

            for phph, phase in enumerate(alv["phases"]):
                phdate_dict[pt][phase] = datetime.datetime.strptime(
                    alv["dates"][phph], "%m/%d/%Y"
                )

        self.pd_dict = phdate_dict

    def get_date_phase(self, pt, datet):
        # Given a patient and a date, return the PHASE of the study
        # importantly,
        searchstruct = self.pd_dict[pt]

        # find distance between the datetime provided and ALL phases
        dist_to_ph = {key: datet - val for key, val in searchstruct.items()}
        # only keep the phases that are AFTER the datetime provided
        phases_after = {
            key: val for key, val in dist_to_ph.items() if val <= datetime.timedelta(0)
        }

        if bool(phases_after):
            closest_phase = max(phases_after, key=phases_after.get)
        else:
            closest_phase = None

        return closest_phase

    def check_empty_phases(self):
        empty_phases = [rr["Filename"] for rr in self.file_meta if rr["Phase"] == None]

        if len(empty_phases):
            print("Some Empty Phases!")

    def meta_files(self, mode="Chronic"):
        # Here we're loading in the files that are in the MODE that we want
        # So far, the primary mode here is CHRONIC which consists of ambulatory recordings using the PC+S

        file_meta = [{} for _ in range(len(self.file_list))]

        for ff, filen in enumerate(self.file_list):
            # we're going to to each and every file now and give it its metadata

            # Here, we're going to extract the DATE
            file_dateinfo = self.extract_date(filen)
            file_typeinfo = self.extract_rectype(filen)
            # file_gaininfo = self.extract_gains(filen)
            file_ptinfo = self.extract_pt(filen)
            file_phaseinfo = self.get_date_phase(file_ptinfo, file_dateinfo)
            file_dayniteinfo = self.get_time(filen)

            if file_typeinfo == mode:
                file_meta[ff].update(
                    {
                        "Filename": filen,
                        "Date": file_dateinfo,
                        "Type": file_typeinfo,
                        "Patient": file_ptinfo,
                        "Phase": file_phaseinfo,
                        "Circadian": file_dayniteinfo,
                        "BadFlag": False,
                    }
                )
            else:
                file_meta[ff] = None

        # remove all the NONEs since they have nothing to do with the current analysis mode
        file_meta = [x for x in file_meta if x is not None]

        self.file_meta = file_meta

    def check_meta(self, prob_condit=0):
        for rr in self.file_meta:
            for ch in ["Left", "Right"]:
                if rr["Data"][ch].all() == 0:
                    print("PROBLEM: " + str(rr) + " has a zero PSD in channel " + ch)

    def prune_meta(self):
        print("Pruning out recordings that have no Phase in main study...")
        # prune out the parts of file_meta that are not in the study
        new_meta = [rr for rr in self.file_meta if rr["Phase"] != None]

        self.file_meta = new_meta

    def Import_Frame(self, preBuilt):
        print("Loading data from..." + self.im_root_dir)
        self.file_meta = np.load(self.preFrame_file, allow_pickle=True)

    """Load in the data"""

    def Load_Data(self, domain="F"):
        # This is the main function that loads in our FEATURES
        # This function should ALWAYS crawl the file structure and bring in the raw data
        # Use Import_Data to bring in an intermediate file

        if domain == "F":
            self.data_basis[domain] = np.linspace(0, self.fs / 2, 2**9 + 1)
        elif domain == "T":
            self.data_basis[domain] = np.linspace(0, self.sec_end)

        for rr in self.file_meta:
            # load in the file
            print("Loading in " + rr["Filename"])
            precheck_data = self.load_file(rr["Filename"], domain=domain)

            if precheck_data["Left"].all() != 0 and precheck_data["Right"].all() != 0:
                rr.update({"Data": precheck_data})
            else:
                rr.update({"BadFlag": True})

    """Saves the frame to the intermediate directory"""

    def Save_Frame(self, name_addendum=""):
        print("Saving File Metastructure in /tmp/")
        # np.save(self.im_root_dir + '/Chronic_Frame' + name_addendum + '.npy',self.file_meta)

        # Try pickling below
        with open("/tmp/Chronic_Frame" + name_addendum + ".pickle", "wb") as file:
            pickle.dump(
                self,
                file,
            )

    """Here we load in the file *AND* do some preliminary Fourier analysis"""

    def load_file(self, fname, load_intv=(0, -1), domain="T"):
        # call the DBSOsc br_load_method
        # should be 1:1 from file_meta to ts_data

        # this returns the full recording
        txtdata = load_br_file(fname)

        # take just the last 10 seconds
        sec_end = self.sec_end

        # extract channels
        X = {
            "Left": txtdata[-(422 * sec_end) : -1, 0].reshape(-1, 1),
            "Right": txtdata[-(422 * sec_end) : -1, 2].reshape(-1, 1),
        }

        F = defaultdict(dict)

        if domain == "T":

            return X

        elif domain == "F":
            # we want to return the FFT, not the timedomain signal
            # This saves a lot of RAM but obviously has its caveats
            # for this, we want to call the DBS_Osc method for doing FFTs
            # The return from gen_psd is a dictionary eg: {'Left':{'F','PSD'},'Right':{'F','PSD'}}
            F = gen_psd(X)  # JUST CHANGED THIS TO abs 12/15/2020

            return F

    #%%
    # Plotting methods in the class
    def plot_file_PSD(self, fname=""):
        if fname != "":
            psd_interest = [
                (rr["Data"]["Left"], rr["Data"]["Right"])
                for rr in self.file_meta
                if rr["Filename"] == fname
            ]

        plt.figure()
        plt.plot(psd_interest)

    def plot_TS(
        self,
    ):
        w = 10 / 211
        tss = {"Left": [], "Right": []}
        for ch in ["Left", "Right"]:
            for rr in self.file_meta:
                if "Data" in rr.keys() and rr["Phase"] in Phase_List("ephys"):
                    tss[ch].append(rr["Data"][ch])
        # choose a random sampling
        total_recs = len(tss["Left"])
        rnd_idxs = random.sample(range(1, total_recs), 100)
        for ii in rnd_idxs:

            b, a = signal.butter(5, w, "low")
            try:
                output = signal.filtfilt(b, a, tss["Right"][ii].squeeze())
            except:
                raise Exception("Filtfilt is erroring")
            plt.plot(output, alpha=0.5)

    def plot_PSD(self, pt="901"):
        # generate out F vector
        fvect = np.linspace(0, 211, 513)

        # quick way to plot all of a patient's recording
        # therapy_phases = dbo.all_phases
        psds = {"Left": 0, "Right": 0}
        for ch in ["Left", "Right"]:
            try:
                psds[ch] = np.array(
                    [
                        np.log10(rr["Data"][ch])
                        for rr in self.file_meta
                        if rr["Patient"] == pt and rr["Phase"] in Phase_List("ephys")
                    ]
                ).T
            except:
                raise ValueError("PSDs are erroring in log10 conversion")

        # list2 = np.array([np.log10(rr['Data']['Left']) for rr in DataFrame.file_meta if rr['Patient'] == '901' and rr['Circadian'] == 'night' and rr['Phase'] in dbo.Phase_List('notherapy')]).T

        plt.figure()
        plt.subplot(121)
        plt.plot(fvect, psds["Left"], color="r", alpha=0.01)
        plt.subplot(122)
        plt.plot(fvect, psds["Right"], color="b", alpha=0.01)

        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power (dB)")
        plt.legend({"Therapy", "NoTherapy"})

        #%%
        [
            plt.plot(fvect, np.log10(rr["Data"]["Left"]), alpha=0.1)
            for rr in self.file_meta
            if rr["Patient"] == "901" and rr["Circadian"] == "night"
        ]
        plt.title("Night")

        plt.figure()
        [
            plt.plot(fvect, np.log10(rr["Data"]["Left"]), alpha=0.1)
            for rr in self.file_meta
            if rr["Patient"] == "901" and rr["Circadian"] == "day"
        ]
        plt.title("Day")


if __name__ == "__main__":
    # Unit Test
    # Generate our dataframe
    DataFrame = BR_Data_Tree(premade_frame_file="GENERATE")
    DataFrame.generate_TD_sequence()
    DataFrame.Save_Frame(name_addendum="Dec2020_T")
