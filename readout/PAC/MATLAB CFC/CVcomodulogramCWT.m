function MIs = CVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option)
%CVCOMODULOGRAM Generates a Generalized-Linear-Model Based Comodulogram
%   USAGE: MIs = GLMcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz) 
%   Fb and Fc are the wavelet bandwidth and center frequency parameters
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram

    CVs = zeros(length(freqForPhase),length(freqForAmp));
    coefsForAmp = CWTfilt(sigForAmp,fs,Fb,Fc,freqForAmp);
    for cc = 1:length(freqForAmp)
        ampOsc = abs(coefsForAmp(cc,:));
        [Cxy,~] = mscohere(ampOsc,sigForPhase,[],[],freqForPhase,fs); 
        CVs(:,cc) = Cxy; disp(['Completed: cc = ' num2str(cc)]);
    end
    MIs = atanh(CVs);
    
    if strcmp(option,'Yes')
        imagesc(freqForPhase,freqForAmp,MIs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    end
end