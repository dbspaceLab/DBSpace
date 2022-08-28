import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

from dbspace.utils.functions import unity
from dbspace.signal.dLFP.amp_transfer_functions import hard_amp


class sig_amp:
    def __init__(
        self,
        diff_inst,
        family="perfect",
        vmax=1,
        inscale=1e-3,
        tscale=(-10, 10),
        noise=0,
        sig_amp_gain=1,
        pre_amp_gain=1,
        final_fs=422,
    ):
        self.diff_inst = diff_inst
        self.family = family
        self.inscale = inscale  # This is a hack to bring the inside of the amp into whatever space we want, transform, and then bring it back out. It should really be removed...
        self.tscale = tscale
        self.Fs = final_fs

        self.downsample_factor = int(round(self.diff_inst.full_Fs / final_fs))
        self.sig_amp_gain = sig_amp_gain

        self.set_T_func()

        self.noise = noise

        # Frequency domain analysis
        self.nperseg = 2**9
        self.noverlap = 2**9 - 50

        self.pre_amp_gain = pre_amp_gain

    def set_T_func(self):
        if self.family == "perfect":
            self.Tfunc = unity
        elif self.family == "pwlinear":

            self.Tfunc = hard_amp
        elif self.family == "tanh":
            self.Tfunc = np.tanh

    # THIS JUST FOCUsES ON THE ACTUAL SCALING AND AMPLIFIER PROCESS, ignore noise here
    def V_out(self, V_in):
        # put some noise inside?
        V_out = (
            self.sig_amp_gain
            * self.inscale
            * self.Tfunc(self.pre_amp_gain * V_in / self.inscale)
        )

        return V_out

    # BELOW IS THE MEASUREMENT PROCESS< NOISE
    def gen_sig_amp_output(self, Z1, Z3):
        diff_out = self.diff_out

        y_out = self.V_out(diff_out)

        # do we want to add noise?
        if self.noise:
            y_out += self.noise * np.random.normal(size=y_out.shape)

        return y_out

    def dac_sample(self, input_signal, skips):
        ds_input = input_signal[0::skips]
        new_tvect = np.linspace(self.tscale[0], self.tscale[1], ds_input.shape[0])

        return ds_input, new_tvect

    def simulate(self, Z1, Z3):
        self.diff_out = self.diff_inst.V_out(Z1, Z3)["sim_1"]

        # Here we generate our recording, after the signal amplifier component
        V_preDC = self.gen_sig_amp_output(Z1, Z3)

        V_sampled, self.tvect = self.dac_sample(V_preDC, self.downsample_factor)

        self.simulated_lfp = V_sampled

    def plot_simulated(self, use_windowing="blackmanharris"):

        nperseg = self.nperseg
        noverlap = self.noverlap

        self.F, self.T, self.SGout = sig.spectrogram(
            self.simulated_lfp,
            nperseg=nperseg,
            noverlap=noverlap,
            window=sig.get_window("blackmanharris", nperseg),
            fs=self.Fs,
        )
        _, _, self.SGdiff = sig.spectrogram(
            self.sig_amp_gain * self.diff_out,
            nperseg=nperseg,
            noverlap=noverlap,
            window=sig.get_window(use_windowing, nperseg),
            fs=self.diff_inst.full_Fs,
        )

    """
    Plotting Functions
    """

    def plot_time_dom(self):
        V_out = self.simulated_lfp
        diff_out = self.diff_out
        diff_obj = self.diff_inst

        plt.figure()
        # Plot the input and output voltages directly over time
        plt.plot(self.diff_inst.tvect, diff_out, label="Input Voltage", alpha=0.6)
        plt.plot(self.tvect, V_out, label="Output Voltage", alpha=0.5)
        plt.legend()
        plt.ylim((-1e-2, 1e-2))

    def plot_freq_dom(self):
        V_out = self.simulated_lfp
        diff_out = self.diff_out
        diff_obj = self.diff_inst

        SGdiff = self.SGdiff
        SGout = self.SGout

        nperseg = self.nperseg
        noverlap = self.noverlap

        plt.figure()

        plt.subplot(1, 2, 1)
        t_beg = self.T + diff_obj.tlims[0] < -1
        t_end = self.T + diff_obj.tlims[0] > 1
        Pbeg = np.median(10 * np.log10(SGdiff[:, t_beg]), axis=1)
        Pend = np.median(10 * np.log10(SGdiff[:, t_end]), axis=1)
        plt.plot(self.F, Pbeg, color="black")
        plt.plot(self.F, Pend, color="green")
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power (dB)")
        plt.ylim((-200, -20))
        plt.title("Perfect Amp")

        plt.subplot(1, 2, 2)
        t_beg = self.T + diff_obj.tlims[0] < -1
        t_end = self.T + diff_obj.tlims[0] > 1
        # Below we're just taking the median of the SG. Maybe do the Welch estimate on this?
        Pbeg = np.median(10 * np.log10(SGout[:, t_beg]), axis=1)
        Pend = np.median(10 * np.log10(SGout[:, t_end]), axis=1)
        plt.plot(self.F, Pbeg, color="black")
        plt.plot(self.F, Pend, color="green")
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power (dB)")
        plt.ylim((-200, -20))
        plt.title("Realistic Amp")
        # plt.suptitle('Zdiff = ' + str(np.abs(Z1 - Z3)))

    def plot_tf_dom(self):

        V_out = self.simulated_lfp
        diff_out = self.diff_out
        diff_obj = self.diff_inst

        SGout = self.SGout
        SGdiff = self.SGdiff

        nperseg = self.nperseg
        noverlap = self.noverlap

        plt.figure()
        # Plot the T-F representation of both input and output
        plt.subplot(2, 2, 1)
        # Here, we find the spectrogram of the output from the diff_amp, should not be affected at all by the gain, I guess...
        # BUT the goal of this is to output a perfect amp... so maybe this is not ideal since the perfect amp still has the gain we want.

        plt.pcolormesh(
            self.T + diff_obj.tlims[0], self.F, 10 * np.log10(SGdiff), rasterized=True
        )
        plt.clim(-120, 0)
        plt.ylim((0, 200))
        plt.title("Perfect Amp Output")
        # plt.colorbar()

        plt.subplot(2, 2, 2)
        plt.clim(-120, 0)

        plt.pcolormesh(
            self.T + diff_obj.tlims[0], self.F, 10 * np.log10(SGout), rasterized=True
        )
        plt.title("Imperfect Amp Output")
        # plt.colorbar()

    def plot_osc_power(self):
        # Now we move on to the oscillatory analyses
        plt.figure()
        plt.suptitle("Oscillatory Analyses")

        plt.subplot(3, 2, 1)
        # plot PSDs here
        bl_ts = {0: V_out[: 5 * 422].reshape(-1, 1)}
        stim_ts = {0: V_out[: -5 * 422].reshape(-1, 1)}

        bl_psd = dbo.gen_psd(bl_ts, polyord=0)
        stim_psd = dbo.gen_psd(stim_ts, polyord=0)
        Frq = np.linspace(0, 211, bl_psd[0].shape[0])

        plt.plot(Frq, np.log10(bl_psd[0]))
        plt.plot(Frq, np.log10(stim_psd[0]))
        plt.ylim((-70, 0))

        plt.subplot(3, 2, 3)
        # bl_osc = dbo.calc_feats(bl_psd[0].reshape(-1,1),Frq)
        # stim_osc = dbo.calc_feats(stim_psd[0].reshape(-1,1),Frq)
        # plt.bar([bl_osc,stim_osc])

        plt.subplot(3, 2, 2)
        bl_psd = dbo.gen_psd(bl_ts, polyord=4)
        stim_psd = dbo.gen_psd(stim_ts, polyord=4)
        plt.plot(Frq, np.log10(bl_psd[0][0]))
        plt.plot(Frq, np.log10(stim_psd[0][0]))

        # do band power now
        sg_avg = False
        if sg_avg:
            # these plot the average of the Spectrogram
            plt.subplot(3, 2, 3)
            # plt.plot(F,Pbeg)
            # plt.plot(F,Pend)

            plt.subplot(3, 2, 4)
            bl_P = {0: 10 ** (Pbeg)}
            stim_P = {0: 10 ** (Pend)}

            corr_bl = dbo.poly_subtr(bl_P, F)
            corr_stim = dbo.poly_subtr(stim_P, F)
            # plt.plot(F,np.log10(corr_bl[0]))
            # plt.plot(F,np.log10(corr_stim[0]))

        plt.subplot(3, 2, 3)
        # do corrected PSD here
        plt.suptitle("Zdiff = " + str(np.abs(Z1 - Z3)))
