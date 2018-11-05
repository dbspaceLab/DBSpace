function MIs = PhaseLocValCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
%PHASELOCVALCFC Calculates and displays the CFC Comulolograms based on inputs
%   USAGE: MIs = PhaseLocValCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
%   oscAmpMod is a cell matrix of time-series oscillations bandpassed 
%   around freqForAmp with bandwidth specified by freqForPhase.
%   oscForPhase is a row cell vector of time-series oscillations bandpassed
%   around freqForPhase with some small bandwidth   

    % Applying Envelope-to-Signal-Correlation based CFC to Oscillation Data
    MIs = zeros(length(freqForPhase),length(freqForAmp));
    % Phases will change each row. Amplitudes will change each column
    for cc = 1:length(freqForAmp)
        for rr = 1:length(freqForPhase)
            ampOsc = abs(oscAmpMod{rr,cc});
            phaseOsc = angle(oscForPhase{1,rr});
            ampOscPhase = angle(hilbert(ampOsc));
            PLV = abs(mean(exp(1i*(phaseOsc - ampOscPhase))));
            MIs(rr,cc) = PLV;
            disp(['Completed: rr = ' num2str(rr) ', cc = ' num2str(cc)]);
        end
    end
end