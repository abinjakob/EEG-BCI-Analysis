function [ChanPSD, f] = computePSD(EEG, window_length)
% function function [power, f] = computePSD(EEG, window_length)
%
% This function is used to calculate PSD using welch method 
% The function loads '.set' epoched EEG file for the respective SSVEP event 
% ad calculates the PSD for each trial. 
% 
% Inputs:
%   EEG (struct)       : EEGLab EEG epoched data 
%   window_lengh (int) : length of the window for pwelch
%
% Ouput:
%   ChanPSD (2D array) : PSD values averaged across trial for each channel (PSD values x channel)
%   f (1D array)       : frequency vector
%
% Example function call:
% erd = computePSD(EEG, binsize, base_start, base_end)

% overlap of the window
overlap = window_length / 2;

% loop over channels 
for iChan = 1:size(EEG.data,1)
    % loop over trials 
    for iTrial = 1:size(EEG.data,3)
        % computing psd usign pwelch
        [pxx, f] = pwelch(EEG.data(iChan,:,iTrial), hamming(window_length), overlap, 2^nextpow2(window_length*4), EEG.srate);
        % computing psd usign pwelch
        pxx_all(:,iChan,iTrial) = pxx;
    end
end
% average across trials
ChanPSD = mean(pxx_all,3);


