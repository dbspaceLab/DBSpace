function MIs = KullLeibCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase,n,option)
%KULLLEIBCFC Calculates and displays the CFC Comulolograms based on inputs
%   USAGE: MIs = KullLeibCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase,n,option)
%   oscAmpMod is a cell matrix of time-series oscillations bandpassed 
%   around freqForAmp with bandwidth specified by freqForPhase.
%   oscForPhase is a row cell vector of time-series oscillations bandpassed
%   around freqForPhase with some small bandwidth
%   n is the number phasebins for the Kullback-Liebler Modulation Index(MI)
    
    % Applying Kullback-Leibler Divergence-based CFC to Oscillation Data
    phaseBins = -pi:(2*pi/n):pi; highFreqAmplitude = zeros(1,n);
    MIs = zeros(length(freqForPhase),length(freqForAmp));
    % Phases will change each row. Amplitudes will change each column
    for cc = 1:length(freqForAmp)
        for rr = 1:length(freqForPhase)
            amplitudes = abs(oscAmpMod{rr,cc});
            phases = angle(oscForPhase{1,rr});
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

