function MIs = GenLinModCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
%GENLINMODCFC Calculates and displays the CFC Comulolograms based on inputs
%   USAGE: MIs = GenLinModCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
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
            X = [cos(phaseOsc)' sin(phaseOsc)' ones(length(phaseOsc),1)];
            B = (X'*X)\(X')*(ampOsc'); ampOscTrend = X*B; 
            ampOscResid = ampOsc' - ampOscTrend;
            rsq = 1-var(ampOscResid)/var(ampOsc);
            MIs(rr,cc) = sqrt(rsq);
            disp(['Completed: rr = ' num2str(rr) ', cc = ' num2str(cc)]);
        end
    end
end