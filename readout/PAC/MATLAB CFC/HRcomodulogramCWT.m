function MIs = HRcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,method,option)
%HRCOMODULOGRAMCWT Generates a Heights-Ratio Based Comodulogram
%   USAGE: MIs = HRcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,method,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz)
%   Fb and Fc are the wavelet bandwidth and center frequency parameters
%   n is the number phasebins for the Heights-Ratio Modulation Index (MI)
%   method: there are 3 ways of doing this
%       1) 'Lakatos' -- h_max/h_min
%       2) 'Tort' -- (h_max - h_min)/h_max;
%       3) 'AM Radio' --- (h_max - h_min)/(h_max + h_min)
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram
    coefsForAmp = CWTfilt(sigForAmp,fs,Fb,Fc,freqForAmp);
    coefsForPhase = CWTfilt(sigForPhase,fs,Fb,Fc,freqForPhase); 
    MIs = HeightsRatioCFCwav(coefsForAmp,coefsForPhase,freqForAmp,freqForPhase,n,method,option);
end

