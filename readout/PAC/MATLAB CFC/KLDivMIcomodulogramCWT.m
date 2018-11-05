function MIs = KLDivMIcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,option)
%KLDIVMICOMODULOGRAMCWT Generates a Kullback-Liebler-Based Comodulogram
%   USAGE: MIs = KLDivMIcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz)
%   Fb and Fc are the wavelet bandwidth and center frequency parameters
%   n is the number phasebins for the Kullback-Liebler Modulation Index(MI)
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram
    coefsForAmp = CWTfilt(sigForAmp,fs,Fb,Fc,freqForAmp);
    coefsForPhase = CWTfilt(sigForPhase,fs,Fb,Fc,freqForPhase); 
    MIs = KullLeibCFCwav(coefsForAmp,coefsForPhase,freqForAmp,freqForPhase,n,option);
end

