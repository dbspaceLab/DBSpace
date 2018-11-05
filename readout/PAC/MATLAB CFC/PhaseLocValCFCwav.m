function MIs = PhaseLocValCFCwav(coefsForAmp,coefsForPhase,freqForAmp,freqForPhase)
%PHASELOCVALCFCWAV Calculates and displays the CFC Comulolograms based on inputs
%   USAGE: MIs = PhaseLocValCFCwav(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
%   coefsForAmp are wavelet coefficients at freqForAmp 
%   around freqForAmp with bandwidth specified by freqForPhase.
%   coefsForPhase are wavelet coefficients at freqForPhase
%   around freqForPhase with some small bandwidth   

    % Applying Envelope-to-Signal-Correlation based CFC to Oscillation Data
    MIs = zeros(length(freqForPhase),length(freqForAmp));
    % Phases will change each row. Amplitudes will change each column
    for cc = 1:length(freqForAmp)
        for rr = 1:length(freqForPhase)
            ampOsc = abs(coefsForAmp(cc,:));
            phaseOsc = angle(coefsForPhase(rr,:));
            ampOscPhase = angle(hilbert(ampOsc));
            PLV = abs(mean(exp(1i*(phaseOsc - ampOscPhase))));
            MIs(rr,cc) = PLV;
            disp(['Completed: rr = ' num2str(rr) ', cc = ' num2str(cc)]);
        end
    end
end