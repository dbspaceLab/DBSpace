function [MIs MVLs] = ZScoredMVCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
%ZSCOREDMVCFC Calculates and displays the CFC Comulolograms based on inputs
%   USAGE: [MIs MVLs] = ZScoredMVCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase)
%   MIs is the comodulogram based on the z-scored mean vector
%   MVLs is the comodulogram based on Canolty's mean vector length (MVL)
%   oscAmpMod is a cell matrix of time-series oscillations bandpassed 
%   around freqForAmp with bandwidth specified by freqForPhase.
%   oscForPhase is a row cell vector of time-series oscillations bandpassed
%   around freqForPhase with some small bandwidth.   

    % Applying Envelope-to-Signal-Correlation based CFC to Oscillation Data
    MIs = zeros(length(freqForPhase),length(freqForAmp));
    MVLs = zeros(length(freqForPhase),length(freqForAmp));
    % Phases will change each row. Amplitudes will change each column
    for cc = 1:length(freqForAmp)
        for rr = 1:length(freqForPhase)
            ampOsc = abs(oscAmpMod{rr,cc});
            phaseOsc = angle(oscForPhase{1,rr});
            phasor = ampOsc.*exp(1i*phaseOsc);
            MVLs(rr,cc) = abs(mean(phasor));
            phasorComponents = [real(phasor);imag(phasor)];
            MIs(rr,cc) = zScoredMVL(phasorComponents);
            disp(['Completed: rr = ' num2str(rr) ', cc = ' num2str(cc)]);
        end
    end
end