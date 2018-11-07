function MIs = PSDcomodulogramCWT(signal,freqForAmp,freqForPhase,fs,nw,Fb,Fc,option)
%PSDCOMODULOGRAMCWT Generates a Comodulogram from PSD of High Frequency Envelope
%   USAGE: MIs = PSDcomodulogramCWT(signal,freqForAmp,freqForPhase,fs,nw,Fb,Fc,option)
%   signal is the input LFP to be analyzed for phase-amplitude-coupling
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz); nw is the time-halfbandwidth product
%   Fb and Fc are the wavelet bandwidth and center frequency parameters
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram
    
    MIs = zeros(length(freqForPhase),length(freqForAmp));
    coefsForAmp = CWTfilt(signal,fs,Fb,Fc,freqForAmp);
    for cc = 1:length(freqForAmp)
        ampOsc = abs(coefsForAmp(cc,:));
        [pxx,~] = pmtm(ampOsc,nw,freqForPhase,fs);
        MIs(:,cc) = pxx; disp(['Completed: cc = ' num2str(cc)]);
    end
    
    if strcmp(option,'Yes')
        imagesc(freqForPhase,freqForAmp,MIs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    end
end