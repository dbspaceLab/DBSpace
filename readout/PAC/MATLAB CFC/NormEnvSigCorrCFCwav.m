function MIs = NormEnvSigCorrCFCwav(coefsForAmp,coefsForPhase,freqForAmp,freqForPhase)
%NORMENVSIGCORRCFCWAV Calculates and displays the CFC Comulolograms based on inputs
%   USAGE: MIs = NormEnvSigCorrCFCwav(coefsForAmp,coefsForPhase,freqForAmp,freqForPhase)
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
            R = corrcoef(ampOsc,cos(phaseOsc));
            MIs(rr,cc) = R(1,2);
            disp(['Completed: rr = ' num2str(rr) ', cc = ' num2str(cc)]);
        end
    end
end