% SSVEP EEG Analysis Script 
% -------------------------
% The script applied GED to SAT data to improve the ERP pattern, especially
% the N1 component of the auditory response. 
%
% Author: Abin Jacob
%         Carl von Ossietzky University Oldenburg
%         abin.jacob@uni-oldenburg.de            
% Date  : 22/04/2024

%% start fresh 
clear; clc; close all;
% eeglab path 
addpath('L:\Cloud\SW\eeglab2024.0');
% change directory
cd('L:\Cloud\NeuroCFN\RESEARCH PROJECT\Research Project 02\EEG Analysis')

%% load file
% open EEGLab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% folder path
DATAPATH = 'L:\Cloud\NeuroCFN\RESEARCH PROJECT\Research Project 02\EEG Analysis\data_preprocessed';
% filename of the .set file with chanlocs
filename = 'P04_SAT_ICAcleaned.set';
% load the file in EEGLAB
EEG = pop_loadset('filename', filename, 'filepath', DATAPATH);

%% parameters for the analysis

% filtering 
% high-pass filter 
HP = 0.5;                       
HP_order = 826;                    
% low-pass filter  
LP = 30;                        
LP_order = 776;                 

% epoching 
% extracting correct events
[trialCount, corrEvent_left, corrEvent_right] = SATCorrTrials(EEG);
event_name = [corrEvent_left, corrEvent_right];
epoch_start = -0.25;               
epoch_end = 3.25;                 

% baseline correction
fs = EEG.srate;                 
baseline = [epoch_start*fs 0]; 

% reject artefactual epochs 
PRUNE = 4;

%% pre-processing

% filtering 
% low-pass filter 
EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hann', 'forder', LP_order);
% high-pass filter 
EEG = pop_firws(EEG, 'fcutoff', HP, 'ftype', 'highpass', 'wtype', 'hann', 'forder', HP_order);

% re-referencing to CAR
EEG = pop_reref(EEG, [], 'refstate',0);

% removing unnecessary event marker
event_pos = 1;      % position counter for the events other than stim onset
event_idx = [];     % array to store the index of the event other than stim onset

for idx = 1: length(EEG.event)
    if ~ strcmp(EEG.event(idx).type, event_name)
        event_idx(event_pos) = idx;
        event_pos = event_pos + 1;
    end
end 
% remove events which are not stim onset from the data
EEG = pop_editeventvals(EEG, 'delete', event_idx);
EEG = eeg_checkset(EEG);

% epoching 
EEG = pop_epoch(EEG, event_name, [epoch_start epoch_end], 'newname', 'SET_pilot_Epoched', 'epochinfo', 'yes');

% reject artefactual epochs 
% joint probability-based artifact rejection (joint prob. > PRUNE (SD))
EEG = pop_jointprob(EEG, 1, [1:EEG.nbchan], PRUNE, PRUNE, 0, 1, 0);
% save pruned and epoched dataset 
EEG = eeg_checkset(EEG);

% baseline correction
EEG = pop_rmbase(EEG, baseline);
% save pruned and epoched dataset 
EEG = eeg_checkset(EEG);

%% weighting using GED

% selecting left attended epochs 
EEG_left = pop_selectevent(EEG, 'type',corrEvent_left,'renametype','left','deleteevents','off','deleteepochs','on','invertepochs','off');
% selecting right attended epochs 
EEG_right = pop_selectevent(EEG, 'type',corrEvent_right,'renametype','right','deleteevents','off','deleteepochs','on','invertepochs','off');

% time range for computing cov matrices 
% signal is defined as time between 0ms to 3000ms
tidx = dsearchn(EEG.times', [0 3000]');  

% preparing the signal matrix (matS)
matS = EEG_left.data(:,tidx(1):tidx(2),:);
% reshaping the data 
matS = reshape(matS, EEG.nbchan, []);
% making the data mean centred
matS = bsxfun(@minus, matS, mean(matS,2));
% calculating the cov matrix for S
covmatS = (matS * matS') / size(matS,2);

% preparing the reference matrix 
matR = EEG_right.data(:,tidx(1):tidx(2),:);
% reshaping the data 
matR = reshape(matR, EEG.nbchan, []);
% making the data mean centred 
matR = bsxfun(@minus, matR, mean(matR,2));
% calculating the cov matrix for R
covmatR = (matR * matR') / size(matR,2);

% solving GED 
[evecs, evals] = eig(covmatS, covmatR);
% find eigenvalue order 
[~,eigidx] = sort(diag(evals));


% filtleft = abs(hilbert( reshape(EEG_left.data, EEG_left.nbchan, [])' * evecs(:, eigidx(end)) ));
filtleft = reshape(EEG_left.data, EEG_left.nbchan, [])' * evecs(:, eigidx(end)).^2;
filtleft = reshape(filtleft, EEG_left.pnts, EEG_left.trials);

% filtright = abs(hilbert( reshape(EEG_right.data, EEG_right.nbchan, [])' * evecs(:, eigidx(end)) ));
filtright = reshape(EEG_right.data, EEG_right.nbchan, [])' * evecs(:, eigidx(end)).^2;
filtright = reshape(filtright, EEG_right.pnts, EEG_right.trials);


%% 

figure;
% onset of tones
leftTime = 0:3000/4:3000; 
rightTime = 0:3000/5:3000; 
% channel to select for EEG before GED
chan2plot = [31 07 42 02 37 03 32 01 36];

subplot(2,2,1)
plot(EEG.times, mean(filtleft,2), 'Color', [0, 0, 1, 0.5]);
for m = leftTime
    line([m m], ylim, 'Color', [0, 0, 1, 0.5], 'LineStyle', ':', 'LineWidth', 0.5);
end
xlim([0 3000]);
ylim([-30 30]);

subplot(2,2,2)
plot(EEG.times, mean(filtright,2), 'Color', [1, 0, 0, 0.5]);
for m = rightTime
    line([m m], ylim, 'Color', [1, 0, 0, 0.5], 'LineStyle', ':', 'LineWidth', 0.5);
end
xlim([0 3000]);
ylim([-30 30]);

subplot(2,2,3)
plot(EEG_left.times,mean(mean(EEG_left.data(chan2plot,:,:),3),1), 'Color', [0, 0, 1, 0.5])
for m = leftTime
    line([m m], ylim, 'Color', [0, 0, 1, 0.5], 'LineStyle', ':', 'LineWidth', 0.5);
end
xlim([0 3000]);
ylim([-30 30]);

subplot(2,2,4)
plot(EEG_right.times,mean(mean(EEG_right.data(chan2plot,:,:),3),1), 'Color', [1, 0, 0, 0.5])
for m = rightTime
    line([m m], ylim, 'Color', [1, 0, 0, 0.5], 'LineStyle', ':', 'LineWidth', 0.5);
end
xlim([0 3000]);
ylim([-30 30]);


%% single trials 
trial2plot = 1;

figure();
subplot(2,1,1)
plot(EEG.times, filtleft(:, trial2plot), 'Color', [0, 0, 1, 0.5]);
for m = leftTime
    line([m m], ylim, 'Color', [0, 0, 1, 0.5], 'LineStyle', ':', 'LineWidth', 0.5);
end
xlim([0 3000]);
ylim([-80 80]);

subplot(2,1,2)
plot(EEG_left.times,mean(EEG_left.data(chan2plot,:,trial2plot),1), 'Color', [0, 0, 1, 0.5])
for m = leftTime
    line([m m], ylim, 'Color', [0, 0, 1, 0.5], 'LineStyle', ':', 'LineWidth', 0.5);
end
xlim([0 3000]);
ylim([-80 80]);






