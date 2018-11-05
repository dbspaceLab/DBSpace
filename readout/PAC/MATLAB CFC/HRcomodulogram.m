function MIs = HRcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,method,option)
%HRCOMODULOGRAM Generates a Kullback-Liebler-Based Comodulogram
%   USAGE: MIs = HRcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,method,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz), bw is the bandwidth of the bandpass filters typically (4.5 Hz)
%   passbandRipl is on a linear scale (not decibel): its preferred value is 0.02
%   n is the number phasebins for the Heights-Ratio Modulation Index (MI)
%   method: there are 3 ways of doing this
%       1) 'Lakatos' -- h_max/h_min
%       2) 'Tort' -- (h_max - h_min)/h_max;
%       3) 'AM Radio' --- (h_max - h_min)/(h_max + h_min)
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram
    oscAmpMod = CFCfilt(sigForAmp,freqForAmp,freqForPhase,fs,passbandRipl);
    oscForPhase = CFCfilt(sigForPhase,freqForPhase,bw,fs,passbandRipl);
    MIs = HeightsRatioCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase,n,method,option);
end