# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 17:11:22 2015

@author: rali

Researcher: REHMAN ALI
Principal Investigator: ROBERT BUTERA
Graduate Collaborator: VINEET TIRUVADI
Neurolab at Georgia Institute of Technology
"""

import numpy as np
import time
from scipy import signal as sig
from PyPAC import *

randn = np.random.randn
rand = np.random.uniform
fs = 1024.599795; # Sampling Rate


# Synthetic Signal as per Tort (2010)
fEnv = 6; fCar = 60; n = 10000; nn = np.arange(n)+1;
phi = np.pi/4; A = 1.2; B = 1.1; C = 0.4; D = 0.7;
AFM = (A+B*np.cos(2*np.pi*fEnv*nn/fs))*np.cos(2*np.pi*fCar*nn/fs) \
    + C*np.cos(2*np.pi*fEnv*nn/fs+phi) + D*randn(n); 

'''
# Synthetic Signal as per Tort (2010) TESTING CV vs PSD
fEnv = 6; fCar = 60; n = 10000; nn = np.arange(n)+1; fphase = 8;
phi = np.pi/4; A = 1.2; B = 1.1; C = 0.4; D = 0.7;
AFM = (A+B*np.cos(2*np.pi*fEnv*nn/fs))*np.cos(2*np.pi*fCar*nn/fs) \
    + C*np.cos(2*np.pi*fphase*nn/fs+phi) + D*randn(n); 
'''

"""    
# Synthetic Signal as per Tort (2010)
fEnv1 = 12; fCar1 = 35; n1 = 30000; nn1 = np.arange(n1)+1;
phi1 = np.pi/4; A1 = 1.2; B1 = 1.1; C1 = 0.4; D1 = 0.7;
AFM1 = (A1+B1*np.cos(2*np.pi*fEnv1*nn1/fs))*np.cos(2*np.pi*fCar1*nn1/fs) \
    + C1*np.cos(2*np.pi*fEnv1*nn1/fs+phi1) + D1*randn(n1);
fEnv2 = 4; fCar2 = 70; n2 = 30000; nn2 = np.arange(n2)+1;
phi2 = np.pi/4; A2 = 1.2; B2 = 1.1; C2 = 0.4; D2 = 0.7;
AFM2 = (A2+B2*np.cos(2*np.pi*fEnv2*nn2/fs))*np.cos(2*np.pi*fCar2*nn2/fs) \
    + C2*np.cos(2*np.pi*fEnv2*nn2/fs+phi2) + D2*randn(n2);
AFM = AFM1+AFM2;
"""
    
# Sigmoidal Coupling as per Penny (2008).
fEnv = 4; fCar = 60; n = 10000; nn = np.arange(n)+1;
phi = np.pi/4; aEnv = 1; k = 2; c = 1; tc = 0.95; 
xEnv = aEnv*np.sin(2*np.pi*fEnv*nn/fs); 
aCar = k/(1+np.exp(-c*xEnv-tc)); 
xCar = aCar*np.sin(2*np.pi*fCar*nn/fs);
xEnvDelayed = aEnv*np.sin(2*np.pi*fEnv*nn/fs+phi); stdDev = 1;
xSigCoup = xEnvDelayed + xCar + stdDev*randn(n);

# Von-Mises Coupling as per Penny (2008)
fEnv = 6; fCar = 35; n = 10000; nn = np.arange(n)+1; phi = np.pi/4; aEnv = 1; 
phiEnv = 2*np.pi*fEnv*nn/fs; xEnv = aEnv*np.sin(phiEnv); c = 1; lda = 0.95; 
aCar = (c/np.exp(lda))*np.exp(lda*np.cos(phiEnv-phi)); 
xCar = aCar*np.sin(2*np.pi*fCar*nn/fs); stdDev = 1.5;
xVonMises = xEnv + xCar + stdDev*randn(n);

# Biphasic Coupling as per Penny (2008).
fEnv = 6; fCar = 35; n = 10000; nn = np.arange(n)+1; phi = np.pi/4; 
aEnv = 1; xEnv = aEnv*np.sin(2*np.pi*fEnv*nn/fs); 
k1 = 8; k2 = 4; c1 = -10; c2 = 10; tc1 = -0.95; tc2 = 0.95; aCar = 2; 
aCar1 = k1/(1+np.exp(-c1*xEnv-tc1)); 
aCar2 = k2/(1+np.exp(-c2*xEnv-tc2));
s1 = rand(size=n) > 0.5; s2 = rand(size=n) > 0.5;
aCarTotal = s1*aCar1 + s2*aCar2 + aCar;
xCar = aCarTotal*np.sin(2*np.pi*fCar*nn/fs);
xEnvDelayed = aEnv*np.sin(2*np.pi*fEnv*nn/fs+phi); stdDev = 1;
xBiphasic = xEnvDelayed + xCar + stdDev*randn(n);

    
# Band-pass filtering
starttime = time.time(); passbandRipl = 0.02;
freqForAmp = 1.5*np.arange(1,61); freqForPhase = np.arange(1,61)/4+1;
MIs, comodplt = PSDcomod(1*AFM,1*AFM,freqForAmp,freqForPhase,fs,bw=6);
#MIs, comodplt = CVcomod(AFM,AFM,freqForAmp,freqForPhase,fs,bw=6);
#MIs, MVLs, comodplt = zScoreMVcomod(AFM,AFM,freqForAmp,freqForPhase,fs,bw=2);
#MIs, MVLs, comodplt = zScoreMVcomodCWT(AFM,AFM,freqForAmp,freqForPhase,fs,sd_rel_phase=0.14,sd_rel_amp=40);
#MIs, comodplt = PLVcomod(AFM,AFM,freqForAmp,freqForPhase,fs,bw=1.5);
#MIs, comodplt = GLMcomod(AFM,AFM,freqForAmp,freqForPhase,fs,bw=1.5);
#MIs, comodplt = GLMcomodCWT(AFM,AFM,freqForAmp,freqForPhase,fs,sd_rel_phase=0.14,sd_rel_amp=40);
#MIs, comodplt = ESCcomod(AFM,AFM,freqForAmp,freqForPhase,fs,bw=2);
#MIs, comodplt = ESCcomodCWT(AFM,AFM,freqForAmp,freqForPhase,fs,sd_rel_phase=0.14,sd_rel_amp=40);
#MIs, comodplt = NESCcomod(AFM,AFM,freqForAmp,freqForPhase,fs,bw=2);
#MIs, comodplt = NESCcomodCWT(AFM,AFM,freqForAmp,freqForPhase,fs,sd_rel_phase=0.14,sd_rel_amp=40);
#MIs, comodplt = KLDivMIcomod(AFM,AFM,freqForAmp,freqForPhase,fs,bw=1.5);

#MIs, comodplt = HRcomod(AFM,AFM,freqForAmp,freqForPhase,fs,bw=1.5);
endtime = time.time();
elapsedtime = endtime - starttime;
print("Elapsed time is: "+str(elapsedtime));
#plt.title("General Linearized Model");
plt.show();
