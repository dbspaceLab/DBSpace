function [MIs MVLs] = zScoreMVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option)
%ZSCOREMVCOMODULOGRAMCWT Generates a Normalized-Envelope-to-Signal-Correlation Based Comodulogram
%   USAGE: MIs = zScoreMVcomodulogramCWT(sigForAmp,sigForPhase,freqForPhase,fs,Fb,Fc,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz) 
%   Fb and Fc are the wavelet bandwidth and center frequency parameters
%   option is either 'MVL', 'Z-Score', 'None':
%       "MVL" gives the mean vector length based on Canolty's Work
%       "Z-Score" gives z-score of mean vector based on PCA
%       "None" no comodulogram displayed

    coefsForAmp = CWTfilt(sigForAmp,fs,Fb,Fc,freqForAmp);
    coefsForPhase = CWTfilt(sigForPhase,fs,Fb,Fc,freqForPhase); 
    [MIs, MVLs] = ZScoredMVCFCwav(coefsForAmp,coefsForPhase,freqForAmp,freqForPhase);
    
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