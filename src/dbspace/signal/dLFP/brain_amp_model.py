from allantools.noise import pink as f_noise
import numpy as np
import scipy.signal as sig


class brain_sig:
    def __init__(self, fs, freq, ampl, phase=0):
        self.center_freq = freq
        self.amplit = ampl
        self.bg_1f_strength = 1e-7
        self.phase = phase

        self.fs = fs
        tlims = (-10, 10)
        self.tvect = np.linspace(tlims[0], tlims[1], 20 * self.fs)

    def ts_return(self):
        self.do_1f()
        self.do_osc()

        return self.bg_1f + self.brain_osc

    def do_1f(self):
        self.bg_1f = self.bg_1f_strength * np.array(f_noise(self.tvect.shape[0]))

    def do_osc(self, smear=False):
        # this one makes the core sine wave oscillations
        if smear:
            self.brain_osc = (
                10
                * self.amplit
                * np.sin(2 * np.pi * self.center_freq * self.tvect)
                * np.exp(-(self.tvect**2) / (2))
            )
        else:
            self.brain_osc = self.amplit * np.sin(
                2 * np.pi * self.center_freq * self.tvect
            )


class stim_sig:
    def __init__(self, fs, stim_ampl=6, wform="sine", stim_freq=130, zero_onset=True):
        self.center_freq = stim_freq
        self.amplit = stim_ampl
        self.phase = 0
        self.wform = wform
        self.zero_onset = zero_onset

        self.fs = fs
        tlims = (-10, 10)
        self.tvect = np.linspace(tlims[0], tlims[1], 20 * self.fs)

    def ts_return(self):
        self.do_stim(wform=self.wform)
        return self.stim_osc

    def interp_ipg_func(self, tvect):
        # this will never really get called, but it's in here in case
        ipg_infile = "ipg_data/ssipgwave_vreg_Ra1p1kOhm_1usdt.txt"
        inmatr = np.array(pd.read_csv(ipg_infile, sep=",", header=None))

        # concatenate this to massive
        concatstim = np.tile(inmatr[:, 1], 223)
        ccstim_tvect = (
            np.linspace(0, concatstim.shape[0] / 1e6, concatstim.shape[0]) - 10
        )

        # now downsample this using interp
        artif = scipy.interpolate.interp1d(ccstim_tvect, concatstim)
        # save artif into a pickled function

        orig_big_x = artif(self.tvect)

        np.save("/tmp/ipg", orig_big_x)

    def brute_ipg_func(self, decay=30, order=15, Wn=0.5):
        tenth_sec_stim = np.load("/home/vscode/data/stim_waveform/tenth_sec_ipg.npy")
        # gaussian filter versus

        full_stim = np.tile(tenth_sec_stim, 10 * 21)
        expon = sig.exponential(101, 0, decay, False)

        full_stim = np.convolve(full_stim, expon)
        # full_stim = gaussian_filter1d(full_stim,100)
        # 237 gets us up to around 21 seconds at 1Mhz

        # finally, need to highpass filter this
        # b,a = sig.butter(order,Wn,btype='highpass',analog=True)
        # w,h = sig.freqz(b,a)
        # plt.figure()
        # plt.plot(w,20*np.log10(abs(h)))

        # full_stim = sig.filtfilt(b,a,full_stim)
        stim_osc = full_stim[0::2370][0 : self.fs * 20]

        np.save(
            "/home/vscode/data/stim_waveform/stim_wform",
            stim_osc,
        )

    def do_stim(self, wform):
        if wform == "sine":
            # print('Using Simple Sine Stim Waveform')
            self.stim_osc = self.amplit * np.sin(
                2 * np.pi * self.center_freq * self.tvect
            )
        elif wform[0:8] == "realsine":
            self.stim_osc = np.zeros_like(self.tvect)
            for hh in range(1, 2):
                self.stim_osc += (
                    2 * hh * np.sin(2 * np.pi * hh * self.center_freq * self.tvect)
                )
        elif wform[0:8] == "moresine":
            nharm = int(wform[-1])
            hamp = [1, 5, 0.5, 0.3, 0.25]
            self.stim_osc = np.zeros_like(self.tvect)
            for hh in range(0, nharm + 1):

                self.stim_osc += (2 * self.amplit / hamp[hh]) * np.sin(
                    2 * np.pi * hh * self.center_freq * self.tvect
                )

        elif wform == "interp_ipg":
            self.stim_osc = 1 / 20 * self.amplit * np.load("/tmp/ipg.npy")
        elif wform == "square":
            self.stim_osc = self.amplit * (
                np.array(
                    (
                        sig.square(
                            2 * np.pi * self.center_freq * self.tvect,
                            duty=(90e-6 / (1 / self.center_freq)),
                        )
                    ),
                    dtype=float,
                )
                / 2
                + 1 / 2
            )
        elif wform == "ipg":
            print("Using Medtronic IPG Stim Waveform")
            # tenth_sec_stim = np.load('/home/virati/tenth_sec_ipg.npy')
            # full_stim = gaussian_filter1d(np.tile(tenth_sec_stim,10*21),100)
            # self.stim_osc = 10 * self.amplit * full_stim[0::237][0:84400]
            in_wform = np.load("/home/vscode/data/stim_waveform/stim_wform.npy")

            # b,a = sig.butter(10,100/self.fs,btype='highpass')
            # stim_osc = sig.lfilter(b,a,in_wform)
            stim_osc = in_wform
            self.stim_osc = 10 * self.amplit * stim_osc

        # self.stim_osc = sig.detrend(self.stim_osc,type='constant')

        # self.stim_osc = self.amplit * self.ipg_func(self.tvect)
        if self.zero_onset:
            self.stim_osc[0 : int(self.tvect.shape[0] / 2)] = 0
