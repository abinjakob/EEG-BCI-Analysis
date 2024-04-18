function erd = computeERD(EEG, binsize, base_start, base_end, epoch_start, epoch_end)
% function erd = computeERD(EEG, binsize, base_start, base_end, epoch_start, epoch_end)
%
% This function is used to calculate Even Related Desynchronisation (ERD) 
% for Motor Imagery paradigm (MI) based on Pfurtscheller(1999). 
% The function loads '.set' epoched EEG file for the respective MI event 
% calculates the ERD. 
% 
% Inputs:
%   EEG (struct)     : EEGLab EEG epoched data 
%   binsize (int)    : size of the bin
%   base_start (int) : start time for the baseline period
%   base_end (int)   : end time for the baseline period
%   epoch_start (int): start time of the epoch
%   epoch_end (int)  : end time of the epoch
%
% Ouput:
%   erd (2D array) : Erd values for each channel (channel x erd values)
%
% Example function call:
% erd = computeERD(EEG, binsize, base_start, base_end, epoch_start, epoch_end)


% window size
window = EEG.pnts/binsize; 
% calculating time points 
tbins = [];
% loop over window
for t = 1:window
    start = 1+(t-1)*binsize;
    tbins{t} = (start:(start+binsize)-1);
end 

% find the strat and end index of baseline period
epochtime = epoch_start:epoch_end;
baseidx_start = ceil(((find(epochtime == base_start)-1) * EEG.srate) / binsize);
baseidx_end = ceil(((find(epochtime == base_end)-1) * EEG.srate) / binsize);

% squaring the amplitude to obtain power (step 02)
% loop over channels
for iChan = 1:size(EEG.data,1)
    % loop over bins 
    for iBin = 1:length(tbins)    
        % loop over trials
        for iTrial = 1:size(EEG.data,3)
            chanPower_vals(iChan,iBin,iTrial,:) = EEG.data(iChan,tbins{iBin},iTrial).^2;
        end
    end
end

% averagin across trials (step 03)
chanPower_TrialsAvg = squeeze(mean(chanPower_vals,3));
% averaging across time samples (step 04) 
chanPower = mean(chanPower_TrialsAvg,3);

% calculate ERD
% loop over channels 
for iChan = 1:size(chanPower,1)
    % calculating avg power in baseline period (reference period, R)
    baseline_avg = mean(chanPower(iChan, baseidx_start:baseidx_end));
    % calculating ERD% = ((A-R)/R)*100
    erd(iChan,:) = ((chanPower(iChan,:)-baseline_avg)/baseline_avg)*100;
end 





