function oscillations = CFCfilt(signal,freqForAmp,freqForPhase,fs,passbandRipl)
%CFCFILT Returns a matrix of bandpass filtered LFP signals
%   USAGE: oscillations = CFCfilt(signal,freqForAmp,freqForPhase,fs,passbandRipl)
%   signal is the input LFP to be bandpassed, fs is the sampling rate
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   passbandRipl is on a linear scale (not decibel): its preferred value is 0.02
%   oscillations is a cell matrix of complex-valued time-series:
%       rows correspond to frequency for phase
%       columns correspond to frequency for amplitude

    oscillations = cell(length(freqForPhase),length(freqForAmp));
    Rp = 40*log10((1+passbandRipl)/(1-passbandRipl));

    % VARIABLE BAND-PASS FILTERING
    % Columns will vary by center frequency (frequency for amplitude), 
    % but rows will vary by bandwidth (frequency for phase)
    for jj = 1:length(freqForPhase)
        for kk = 1:length(freqForAmp)
            freq = freqForAmp(kk); delf = freqForPhase(jj); 
            if freq > 1.2*delf
                [bb aa] = cheby1(3,Rp,[freq-1.2*delf freq+1.2*delf]./(fs/2),'bandpass');
            else 
                [bb aa] = cheby1(3,Rp,(freq+1.2*delf)./(fs/2),'low');
            end
            oscillation = filtfilt(bb,aa,signal);
            oscillations{jj,kk} = hilbert(oscillation);
            disp(['Completed: jj = ' num2str(jj) ', kk = ' num2str(kk)]);
        end
    end

end

