function MIs = CVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,passbandRipl,bw,option)
%CVCOMODULOGRAM Generates a Generalized-Linear-Model Based Comodulogram
%   USAGE: MIs = CVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,passbandRipl,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz)
%   bw is bandwidth for sigForAmp:
%       if bw == 0: bandwidth = max(freqForPhase)
%       else bandwidth = bw
%   passbandRipl is on a linear scale (not decibel): its preferred value is 0.02
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram

    if bw == 0
        bandwidth = max(freqForPhase);
    else
        bandwidth = bw;
    end
    
    CVs = zeros(length(freqForPhase),length(freqForAmp));
    oscAmpMod = CFCfilt(sigForAmp,freqForAmp,bandwidth,fs,passbandRipl);
    for cc = 1:length(freqForAmp)
        ampOsc = abs(oscAmpMod{1,cc});
        [Cxy,~] = mscohere(ampOsc,sigForPhase,[],[],freqForPhase,fs); 
        CVs(:,cc) = Cxy; disp(['Completed: cc = ' num2str(cc)]);
    end
    MIs = atanh(CVs);
    
    if strcmp(option,'Yes')
        imagesc(freqForPhase,freqForAmp,MIs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    end
end