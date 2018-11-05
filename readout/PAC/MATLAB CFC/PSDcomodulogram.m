function MIs = PSDcomodulogram(signal,freqForAmp,freqForPhase,fs,nw,passbandRipl,option)
%PSDCOMODULOGRAM Generates a Generates a Comodulogram from PSD of High Frequency Envelope
%   USAGE: MIs = PSDcomodulogram(signal,freqForAmp,freqForPhase,fs,nw,passbandRipl,option)
%   signal is the input LFP to be analyzed for phase-amplitude-coupling
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz); nw is the time-halfbandwidth product
%   passbandRipl is on a linear scale (not decibel): its preferred value is 0.02
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram

    MIs = zeros(length(freqForPhase),length(freqForAmp));
    oscAmpMod = CFCfilt(signal,freqForAmp,max(freqForPhase),fs,passbandRipl);
    for cc = 1:length(freqForAmp)
        ampOsc = abs(oscAmpMod{1,cc});
        [pxx,~] = pmtm(ampOsc,nw,freqForPhase,fs);
        MIs(:,cc) = pxx; disp(['Completed: cc = ' num2str(cc)]);
    end
    
    if strcmp(option,'Yes')
        imagesc(freqForPhase,freqForAmp,MIs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    end
end