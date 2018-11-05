function [MIs MVLs] = zScoreMVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option)
%ZSCOREMVCOMODULOGRAM Generates a Mean Vector Length-Based Comodulogram
%   USAGE: MIs = zScoreMVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz), bw is the bandwidth of the bandpass filters typically (4.5 Hz)
%   passbandRipl is on a linear scale (not decibel): its preferred value is 0.02
%   option is either 'MVL' or 'Z-Score':
%       "MVL" gives the mean vector length based on Canolty's Work
%       "Z-Score" gives z-score of mean vector based on PCA
%       "None" no comodulogram displayed

    oscAmpMod = CFCfilt(sigForAmp,freqForAmp,freqForPhase,fs,passbandRipl);
    oscForPhase = CFCfilt(sigForPhase,freqForPhase,bw,fs,passbandRipl);
    [MIs MVLs] = ZScoredMVCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase);
    
    if strcmp(option,'MVL')
        imagesc(freqForPhase,freqForAmp,MVLs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    elseif strcmp(option,'Z-Score')
        imagesc(freqForPhase,freqForAmp,MIs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    elseif strcmp(option,'None')
    else
        error('inputarg:invalid','option must be "MVL" or "Z-Score"');
    end
end