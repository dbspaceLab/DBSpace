function MIs = KullLeibCFCwav(coefsForAmp,coefsForPhase,freqForAmp,freqForPhase,n,option)
%KULLLEIBCFCWAV Calculates and displays the CFC Comulolograms based on inputs
%   USAGE: MI = KullLeibCFCwav(oscAmpMod,oscForPhase,freqForAmp,freqForPhase,n,option)
%   coefsForAmp are wavelet coefficients at freqForAmp 
%   around freqForAmp with bandwidth specified by freqForPhase.
%   coefsForPhase are wavelet coefficients at freqForPhase
%   around freqForPhase with some small bandwidth
%   n is the number phasebins for the Kullback-Liebler Modulation Index(MI)
    
    % Applying Kullback-Leibler Divergence-based CFC to Oscillation Data
    phaseBins = -pi:(2*pi/n):pi; highFreqAmplitude = zeros(1,n);
    MIs = zeros(length(freqForPhase),length(freqForAmp));
    % Phases will change each row. Amplitudes will change each column
    for cc = 1:length(freqForAmp)
        for rr = 1:length(freqForPhase)
            amplitudes = abs(coefsForAmp(cc,:));
            phases = angle(coefsForPhase(rr,:));
            for kk = 1:n
                amps = amplitudes(phases > phaseBins(kk) & phases <= phaseBins(kk+1));
                highFreqAmplitude(kk) = mean(amps);
            end
            highFreqAmpProb = highFreqAmplitude/sum(highFreqAmplitude);
            MIs(rr,cc) = ModulationIndex(highFreqAmpProb);
            disp(['Completed: rr = ' num2str(rr) ', cc = ' num2str(cc)]);
        end
    end
    if strcmp(option,'Yes')
        imagesc(freqForPhase,freqForAmp,MIs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    end

end