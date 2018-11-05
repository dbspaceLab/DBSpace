function MIs = PLVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option)
%PLVCOMODULOGRAMCWT Generates a Phase-Locking-Value Based Comodulogram
%   USAGE: MIs = PLVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz) 
%   Fb and Fc are the wavelet bandwidth and center frequency parameters
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram

    coefsForAmp = CWTfilt(sigForAmp,fs,Fb,Fc,freqForAmp);
    coefsForPhase = CWTfilt(sigForPhase,fs,Fb,Fc,freqForPhase); 
    PLVs = PhaseLocValCFCwav(coefsForAmp,coefsForPhase,freqForAmp,freqForPhase);
    MIs = asin(2*PLVs-1);
    if strcmp(option,'Yes')
        imagesc(freqForPhase,freqForAmp,MIs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    end
end