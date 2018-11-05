function MIs = KLDivMIcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,option)
%KLDIVMICOMODULOGRAM Generates a Kullback-Liebler-Based Comodulogram
%   USAGE: MIs = KLDivMIcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz), bw is the bandwidth of the bandpass filters typically (4.5 Hz)
%   passbandRipl is on a linear scale (not decibel): its preferred value is 0.02
%   n is the number phasebins for the Kullback-Liebler Modulation Index(MI)
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram
    oscAmpMod = CFCfilt(sigForAmp,freqForAmp,freqForPhase,fs,passbandRipl);
    oscForPhase = CFCfilt(sigForPhase,freqForPhase,bw,fs,passbandRipl);
    MIs = KullLeibCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase,n,option);
end

