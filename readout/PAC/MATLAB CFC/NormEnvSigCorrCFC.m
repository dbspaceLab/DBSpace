function MIs = NormEnvSigCorrCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
%NORMENVSIGCORRCFC Calculates and displays the CFC Comulolograms based on inputs
%   USAGE: MIs = NormEnvSigCorrCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
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
            R = corrcoef(ampOsc,cos(phaseOsc));
            MIs(rr,cc) = R(1,2);
            disp(['Completed: rr = ' num2str(rr) ', cc = ' num2str(cc)]);
        end
    end
end