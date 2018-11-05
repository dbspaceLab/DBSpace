function coefs = CWTfilt(signal,fs,Fb,Fc,frequencies)
%CWTFILT Determines Complex Morlet Wavelet CWT Coefficients for Signal
%   USAGE: coefs = CWTfilt(signal,Fb,Fc,frequencies)
%   signal is the input signal, fs is sampling rate (Hz)
%   Fb and Fc are the wavelet bandwidth and center frequency parameters
%   frequencies are the center frequencies for the CWT
    wavename = ['cmor' num2str(Fb) '-' num2str(Fc)];
    scales = Fc*fs./frequencies;
    coefs = cwt(signal,scales,wavename);
    % Lowpass filtering block for noisy CWT coefficients
    for kk = 1:length(frequencies)
        wn = 5*frequencies(kk)/(fs/2); n = 5;
        [bb aa] = butter(n,wn,'low');
        coefs(kk,:) = filtfilt(bb,aa,coefs(kk,:));
    end
end

