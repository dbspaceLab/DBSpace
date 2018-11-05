function filteredLFP = LFPPreProc(rawLFP,fs)
%LFPPREPROC Preprocesses Raw Multichannel EEG Data
%   USAGE: filteredLFP = LFPPreProc(rawLFP,fs)
%   Input rawLFP must have the time-series for each channels along its row
%   Output filteredLFP will be bandpassed at from 1-100 Hz and notched at
%   60 Hz. Sampling rate fs is needed in order to do this

    [bb1 aa1] = butter(5,[1 100]./(fs/2));
    [bb2 aa2] = iirnotch(60/(fs/2),0.01);

    filteredLFP = filtfilt(bb1,aa1,rawLFP'); 
    filteredLFP = filtfilt(bb2,aa2,filteredLFP); 
    filteredLFP = filteredLFP';
end

