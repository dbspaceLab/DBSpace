clear
clc

% Opening the .mat files: SbaselineRaw is preDBS, SendRaw is postDBS
fs = 1024.599795; % Sampling Rate

% Synthetic Signal as per Tort (2010).
fEnv = 4; fCar = 60; n = 10000; nn = 1:n; phi = pi/4; 
A = 1.2; B = 1.1; C = 0.4; D = 0.7;
AFM = (A+B*cos(2*pi*fEnv*nn/fs)).*cos(2*pi*fCar*nn/fs)...
    + C*cos(2*pi*fEnv*nn/fs+phi) + D*randn(1,n);

% Adding FM to Tort's Signal
n = 10000; nn = 1:n; phi = pi/4; 
fosc = 0.5; FMratio = 0.01; % Change fosc from 0.001 to 1 and FMratio from 0.35 to 0.01
FMcomp = 1-FMratio*sin(2*pi*fosc*nn/fs);
fEnv0 = 4; fCar0 = 60;
fEnv = fEnv0*FMcomp; fCar = fCar0*FMcomp;
A = 1.2; B = 1.1; C = 0.4; D = 0.3;
AFMwithFM = (A+B*cos(2*pi*fEnv.*nn/fs)).*cos(2*pi*fCar.*nn/fs)...
    + C*cos(2*pi*fEnv.*nn/fs+phi) + D*randn(1,n);

% Make FM linear
n = 10000; nn = 1:n; phi = pi/4; 
fosc = fs/n; FMratio = 0.1; % Change fosc from 0.001 to 1 and FMratio from 0.35 to 0.01
FMcomp = 1-FMratio*sawtooth(2*pi*fosc*nn/fs);
fEnv0 = 4; fCar0 = 60;
fEnv = fEnv0*FMcomp; fCar = fCar0*FMcomp;
A = 1.2; B = 1.1; C = 0.4; D = 0.3;
AFMwithFMlin = (A+B*cos(2*pi*fEnv.*nn/fs)).*cos(2*pi*fCar.*nn/fs)...
    + C*cos(2*pi*fEnv.*nn/fs+phi) + D*randn(1,n);

% Adding Break to Tort's Signal
n = 10000; nn = 1:n; phi = pi/4; 
FMratio = 0.1;
FMcomp = (1-FMratio)*(nn < n/2)+(1+FMratio)*(nn > n/2);
fEnv0 = 4; fCar0 = 60;
fEnv = fEnv0*FMcomp; fCar = fCar0;
A = 1.2; B = 1.1; C = 0.4; D = 0.3;
AFMwithBreak = (A+B*cos(2*pi*fEnv.*nn/fs)).*cos(2*pi*fCar.*nn/fs)...
    + C*cos(2*pi*fEnv.*nn/fs+phi) + D*randn(1,n);

% Add Gaussian Envelope to Tort's Signal
fEnv = 4; fCar = 60; n = 10000; nn = 1:n; phi = pi/4; 
A = 1.2; B = 1.1; C = 0.4; D = 0.3; 
sigma = 900; delay = 5000;
AFMwithGaussian = ((A+B*cos(2*pi*fEnv*nn/fs)).*cos(2*pi*fCar*nn/fs)...
    + C*cos(2*pi*fEnv*nn/fs+phi)).*exp(-((nn-delay).^2)/(sigma^2)) + D*randn(1,n);

%% Plots and Spectrograms

% Normal AM Signal
% plot(nn/fs, AFM);
% xlabel('time (t)'); ylabel('AM Signal');
% 
% plot(nn/fs, AFMwithFM);
% xlabel('time (t)'); ylabel('AM Signal');
%
plot(nn/fs, AFMwithFMlin);
xlabel('time (t)'); ylabel('AM Signal');
% 
% plot(nn/fs, AFMwithBreak);
% xlabel('time (t)'); ylabel('AM Signal');
% 
% plot(nn/fs, AFMwithGaussian);
% xlabel('time (t)'); ylabel('AM Signal');

% spectrogram(AFM, 'yaxis'); % NOT WORKING ON VINEET'S COMPUTER

%% Code that generates comodulograms Canolty Style
% sigForAmp = AFM; sigForPhase = AFM; numSurr = 1000;
% passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Yes'; % 'No'
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
% [MIs, manyMVLs] = canoltycomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,numSurr,option);
% savefig('NoFMnoCFC_Z.fig');

sigForAmp = AFMwithFMlin; sigForPhase = AFMwithFMlin; numSurr = 1000;
passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Yes'; % 'No'
bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
[MIs, manyMVLs] = canoltycomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,numSurr,option);
savefig('AnswerIsFMlin_Z.fig');

sigForAmp = AFMwithFM; sigForPhase = AFMwithFM; numSurr = 1000;
passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Yes'; % 'No'
bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
[MIs, manyMVLs] = canoltycomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,numSurr,option);
savefig('AnswerIsFM_Z.fig');

% sigForAmp = AFMwithBreak; sigForPhase = AFMwithBreak; numSurr = 1000;
% passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Yes'; % 'No'
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
% [MIs, manyMVLs] = canoltycomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,numSurr,option);
% savefig('CFCwithSignalBreak_Z.fig');

% sigForAmp = AFMwithGaussian; sigForPhase = AFMwithGaussian; numSurr = 1000;
% passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Yes'; % 'No'
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
% [MIs, manyMVLs] = canoltycomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,numSurr,option);
% savefig('CFCwithGaussianTaper_Z.fig');

%% Code that generates sample KLDIV comodulograms (wavelet and without)
% sigForAmp = AFM; sigForPhase = AFM; 
% passbandRipl = 0.02; option = 'Yes'; frequencies = 1.5*(1:60); bw = 2; n = 36; 
% freqForAmp = frequencies; freqForPhase = frequencies/6;
% MIs = KLDivMIcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,option);
% % figure(); Fb = 1; Fc = 1;
% % MIsCWT = KLDivMIcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,option);
% savefig('KLDIV_AFM_Z.fig');

sigForAmp = AFM; sigForPhase = AFM; 
passbandRipl = 0.02; option = 'Yes'; frequencies = 1.5*(1:60); bw = 2; n = 36; 
freqForAmp = frequencies; freqForPhase = frequencies/3;
MIs = KLDivMIcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,option);
% figure(); Fb = 1; Fc = 1;
% MIsCWT = KLDivMIcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,option);
savefig('KLDIV_AFMwithFMlin_Z.fig');

% sigForAmp = AFMwithFM; sigForPhase = AFMwithFM; 
% passbandRipl = 0.02; option = 'Yes'; frequencies = 1.5*(1:60); bw = 2; n = 36; 
% freqForAmp = frequencies; freqForPhase = frequencies/6;
% MIs = KLDivMIcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,option);
% % figure(); Fb = 1; Fc = 1;
% % MIsCWT = KLDivMIcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,option);
% savefig('KLDIV_AFMwithFM_Z.fig');

% sigForAmp = AFMwithBreak; sigForPhase = AFMwithBreak; 
% passbandRipl = 0.02; option = 'Yes'; frequencies = 1.5*(1:60); bw = 2; n = 36; 
% freqForAmp = frequencies; freqForPhase = frequencies/6;
% MIs = KLDivMIcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,option);
% %figure(); Fb = 1; Fc = 1;
% %MIsCWT = KLDivMIcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,option);
% savefig('KLDIV_AFMwithBreak_Z.fig');

% sigForAmp = AFMwithGaussian; sigForPhase = AFMwithGaussian; 
% passbandRipl = 0.02; option = 'Yes'; frequencies = 1.5*(1:60); bw = 2; n = 36; 
% freqForAmp = frequencies; freqForPhase = frequencies/6;
% MIs = KLDivMIcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,n,option);
% %figure(); Fb = 1; Fc = 1;
% %MIsCWT = KLDivMIcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,n,option);
% savefig('KLDIV_AFMwithGaussian_Z.fig');

%% Code that generates sample GLM comodulograms (wavelet and without)
% sigForAmp = AFM; sigForPhase = AFM; passbandRipl = 0.02; frequencies = 1.5*(1:60); 
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6; option = 'Yes';
% MIsESC = GLMcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
% %figure(); Fb = 1; Fc = 1;
% %MIs = GLMcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
% savefig('GLM_AFM_Z.fig');

sigForAmp = AFMwithFMlin; sigForPhase = AFMwithFMlin; passbandRipl = 0.02; frequencies = 1.5*(1:60); 
bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6; option = 'Yes';
MIsESC = GLMcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
%figure(); Fb = 1; Fc = 1;
%MIs = GLMcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
savefig('GLM_AFMwithFMlin_Z.fig');

sigForAmp = AFMwithFM; sigForPhase = AFMwithFM; passbandRipl = 0.02; frequencies = 1.5*(1:60); 
bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6; option = 'Yes';
MIsESC = GLMcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
%figure(); Fb = 1; Fc = 1;
%MIs = GLMcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
savefig('GLM_AFMwithFM_Z.fig');

% sigForAmp = AFMwithBreak; sigForPhase = AFMwithBreak; passbandRipl = 0.02; frequencies = 1.5*(1:60); 
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6; option = 'Yes';
% MIsESC = GLMcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
% %figure(); Fb = 1; Fc = 1;
% %MIs = GLMcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
% savefig('GLM_AFMwithBreak_Z.fig');
% 
% sigForAmp = AFMwithGaussian; sigForPhase = AFMwithGaussian; passbandRipl = 0.02; frequencies = 1.5*(1:60); 
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6; option = 'Yes';
% MIsESC = GLMcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
% %figure(); Fb = 1; Fc = 1;
% %MIs = GLMcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
% savefig('GLM_AFMwithGaussian_Z.fig');

%% Code that generates sample PLV comodulograms (wavelet and without)
% sigForAmp = AFM; sigForPhase = AFM; 
% passbandRipl = 0.02; frequencies = 1.5*(1:60); 
% bw = 4.5; freqForAmp = frequencies; freqForPhase = frequencies/12; option = 'Yes';
% MIsESC = PLVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
% figure(); Fb = 1; Fc = 1;
% MIs = PLVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);

%% Code that generates sample mean vector z-scores and mean vector lengths (wavelet and without)
% sigForAmp = AFM; sigForPhase = AFM; 
% passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Z-Score'; % 'None','MVL'
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
% [MIs MVLs] = zScoreMVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
% %figure(); Fb = 1; Fc = 1;
% %[MIsCWT MVLsCWT] = zScoreMVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
% savefig('PCA_AFM_Z.fig');

sigForAmp = AFMwithFMlin; sigForPhase = AFMwithFMlin; 
passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Z-Score'; % 'None','MVL'
bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
[MIs MVLs] = zScoreMVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
%figure(); Fb = 1; Fc = 1;
%[MIsCWT MVLsCWT] = zScoreMVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
savefig('PCA_AFMwithFMlin_Z.fig');

sigForAmp = AFMwithFM; sigForPhase = AFMwithFM; 
passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Z-Score'; % 'None','MVL'
bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
[MIs MVLs] = zScoreMVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
%figure(); Fb = 1; Fc = 1;
%[MIsCWT MVLsCWT] = zScoreMVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
savefig('PCA_AFMwithFM_Z.fig');

% sigForAmp = AFMwithBreak; sigForPhase = AFMwithBreak; 
% passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Z-Score'; % 'None','MVL'
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
% [MIs MVLs] = zScoreMVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
% %figure(); Fb = 1; Fc = 1;
% %[MIsCWT MVLsCWT] = zScoreMVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
% savefig('PCA_AFMwithBreak_Z.fig');

% sigForAmp = AFMwithGaussian; sigForPhase = AFMwithGaussian; 
% passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Z-Score'; % 'None','MVL'
% bw = 2; freqForAmp = frequencies; freqForPhase = frequencies/6;
% [MIs MVLs] = zScoreMVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,bw,passbandRipl,option);
% %figure(); Fb = 1; Fc = 1;
% %[MIsCWT MVLsCWT] = zScoreMVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
% savefig('PCA_AFMwithGaussian_Z.fig');

% %% Code that generates coherence value comodulograms (wavelet and without)
% sigForAmp = AFM; sigForPhase = AFM;
% passbandRipl = 0.02; frequencies = 1.5*(1:60); option = 'Yes';
% freqForAmp = frequencies; freqForPhase = frequencies/12;
% MIs = CVcomodulogram(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,passbandRipl,option);
% figure(); Fb = 1; Fc = 1;
% MIsCWT = CVcomodulogramCWT(sigForAmp,sigForPhase,freqForAmp,freqForPhase,fs,Fb,Fc,option);
