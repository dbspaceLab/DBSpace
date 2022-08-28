import numpy as np
import matplotlib.pyplot as plt

import scipy.signal as sig

from dbspace.signal.PAC.PyPAC import *
from dbspace.signal.dLFP.brain_amp_model import stim_sig, brain_sig


# plt.close('all')

import seaborn as sns

sns.set_context("paper")

sns.set(font_scale=2)
sns.set_style("ticks")
sns.set_style("white")
plt.rcParams["image.cmap"] = "jet"


class diff_amp:
    def __init__(
        self,
        Z_b=1e4,
        Ad=250,
        wform="moresine3",
        clock=False,
        zero_onset=True,
        stim_v=6,
        stim_freq=130,
        fullFs=4220,
    ):
        self.fullFs = fullFs
        self.tlims = (-10, 10)
        self.analogtvect = np.linspace(
            self.tlims[0], self.tlims[1], (self.tlims[1] - self.tlims[0]) * self.fullFs
        )
        self.tvect = np.linspace(
            self.tlims[0], self.tlims[1], (self.tlims[1] - self.tlims[0]) * self.fullFs
        )

        self.Ad = Ad
        self.Z_b = Z_b

        self.c12 = 1
        self.c23 = 1

        # self.osc_params = {'x_1':[12,3e-7],'x_3':[0,0],'x_2':[0,0]} #FOR PAPER
        self.osc_params = {"x_1": [18, 7e-7], "x_3": [0, 0], "x_2": [0, 0]}
        self.set_brain()
        self.set_stim(
            wform=wform, zero_onset=zero_onset, freq=stim_freq, stim_ampl=stim_v
        )
        self.clockflag = clock

        if self.clockflag:
            self.set_clock()

    def set_brain(self, params=[]):
        if params == []:
            params = self.osc_params

        self.X = {"x_1": [], "x_2": [], "x_3": []}

        for bb, sett in params.items():
            self.X[bb] = sig.detrend(
                brain_sig(self.fullFs, sett[0], sett[1]).ts_return(), type="constant"
            )

    def set_stim(self, wform, zero_onset, freq=130, stim_ampl=6):
        decay_factor = 1e-3
        # WARNING, need to figure out why this part is necessary
        stim_scaling = 10
        self.S = (
            stim_scaling
            * decay_factor
            * stim_sig(
                fs=self.fullFs,
                stim_ampl=stim_ampl,
                wform=wform,
                zero_onset=zero_onset,
                stim_freq=freq,
            ).ts_return()
        )

    def set_clock(self, clock_V=2e-3):
        self.clock = stim_sig(
            fs=self.fullFs,
            stim_ampl=clock_V,
            wform="sine",
            stim_freq=105.5,
            zero_onset=False,
        ).ts_return()

    def Vd_stim(self, Z1, Z3):
        self.stim_component = (
            self.Ad
            * self.Z_b
            * self.S
            * ((1 / (Z1 + self.Z_b)) - (1 / (Z3 + self.Z_b)))
        )

    def Vd_x2(self, Z1, Z3):
        self.x2_component = (
            self.Ad
            * self.Z_b
            * self.X["x_2"]
            * ((1 / (self.c12 * (Z1 + self.Z_b))) - (1 / (self.c23 * (Z3 + self.Z_b))))
        )

    def Vd_brain(self, Z1, Z3):
        self.brain_component = (
            self.Ad
            * self.Z_b
            * ((self.X["x_1"] / (Z1 + self.Z_b)) - (self.X["x_3"] / (Z3 + self.Z_b)))
        )

    def V_out(self, Z1, Z3):
        self.Vd_stim(Z1, Z3)
        self.Vd_x2(Z1, Z3)
        self.Vd_brain(Z1, Z3)

        # amplitudes should be determined HERE

        Vo = (self.brain_component + self.x2_component) + (self.stim_component)

        if self.clockflag:
            Vo += self.clock

        self.outputV = Vo
        # first, filter
        b, a = sig.butter(5, 100 / self.fullFs, btype="lowpass")
        # b,a = sig.ellip(4,4,5,100/2110,btype='lowpass')
        Vo = sig.lfilter(b, a, Vo)

        return {"sim_1": Vo / 2}

    def plot_V_out(self, Z1, Z3):
        plot_sig = self.V_out(Z1, Z3)

        plt.figure()
        # dbo.plot_T(plot_sig)
        plt.subplot(211)
        plt.plot(self.tvect, plot_sig["sim_1"])
        plt.xlim(self.tlims)

        nperseg = 2**9
        noverlap = 2**9 - 50
        F, T, SG = sig.spectrogram(
            plot_sig["sim_1"],
            nperseg=nperseg,
            noverlap=noverlap,
            window=sig.get_window("blackmanharris", nperseg),
            fs=self.Fs,
        )
        plt.subplot(212)

        plt.pcolormesh(T + self.tlims[0], F, 10 * np.log10(SG), rasterized=True)
