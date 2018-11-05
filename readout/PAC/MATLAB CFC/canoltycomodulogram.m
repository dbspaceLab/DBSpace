function [MIs, manyMVLs] = canoltycomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,numSurr,option)
%CANOLTYCOMODULOGRAM Generates a comodulogram based on Canolty's method
%   USAGE: MIs = zScoreMVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option)
%   sigForAmp is the input LFP to be analyzed for amplitude
%   sigForPhase is the input LFP to be analyzed for phase
%   freqForAmp is a vector of center frequencies (frequency for amplitude)
%   freqForPhase is a vector of frequency for phase controlling bandwidth
%   fs is sampling rate (Hz), bw is the bandwidth of the bandpass filters typically (4.5 Hz)
%   passbandRipl is on a linear scale (not decibel): its preferred value is 0.02
%   numSurr is number of surrogates by Canolty's mean vector method
%   option: 'Yes' show comodulogram; 'No' don't show comodulogram

    oscAmpMod = CFCfilt(sigForAmp,freqForAmp,freqForPhase,fs,passbandRipl);
    oscForPhase = CFCfilt(sigForPhase,freqForPhase,bw,fs,passbandRipl);
    [~, MVLs] = ZScoredMVCFC(oscAmpMod,oscForPhase,freqForAmp,freqForPhase);
    
    numpoints = length(sigForAmp); % sigForPhase would work too
    minskip = fs; maxskip = numpoints - fs;
    skip = ceil(numpoints.*rand(numSurr*2));
    skip(skip > maxskip & skip < minskip) = [];
    skip = skip(1:numSurr);
    [r, c] = size(MVLs);
    MIs = zeros(r,c);
    oscAmpModOld = oscAmpMod;
    oscAmpMod = zeros(r,c,numpoints);
    for row=1:r
        for col=1:c
            oscAmpMod(row,col,:) = oscAmpModOld{row,col};
        end
    end
    manyMVLs = zeros(r, c, numSurr);
    for s=1:numSurr
        surrAmpMod = cat(3, oscAmpMod(:,:,skip(s):end), oscAmpMod(:,:,1:(skip(s)-1)));
        surrOscAmpMod = cell(r,c);
        for row=1:r
            for col=1:c
                surrOscAmpMod{row,col} = reshape(surrAmpMod(row,col,:),[1, numpoints]);
            end
        end
        [~, manyMVLs(:,:,s)] = ZScoredMVCFC(surrOscAmpMod,oscForPhase,freqForAmp,freqForPhase);
        disp(['Number of Surrogates Left to go: ', num2str(numSurr-s)]);
    end
    for row=1:r
        for col=1:c
            [mean,stddev] = normfit(reshape(manyMVLs(row,col,:),1,numSurr));
            MIs(row,col) = (MVLs(row,col)-mean)/stddev;
        end
    end
    if strcmp(option,'Yes')
        imagesc(freqForPhase,freqForAmp,MIs'); set(gca,'YDir','normal');
        xlabel('Frequency for Phase'); ylabel('Frequency for Amplitude');
    end
end

