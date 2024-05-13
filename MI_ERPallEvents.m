% The script is used to coompute the ERP for different events in the
% dataset. 
%
% Author: Abin Jacob
%         Carl von Ossietzky University Oldenburg
%         abin.jacob@uni-oldenburg.de            
% Date  : 06/05/2024

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
filename = 'P04_MI_ICAcleaned.set';
% load the file in EEGLAB
EEG = pop_loadset('filename', filename, 'filepath', DATAPATH);

%% parameters for the analysis 

% broad band filtering 
% high-pass filter 
HP = 1;                         % cut-off
HP_order = 826;                 % filter order    
% low-pass filter  
LP = 40;                        % cut-off
LP_order = 776;                 % filter order 

% accounting for eeg offset from LSL (100ms)  
offset = 0.1;  
for idx = 1:length(EEG.event)
    EEG.event(idx).latency = EEG.event(idx).latency + (offset*EEG.srate);
end
% checkset 
EEG = eeg_checkset(EEG, 'eventconsistency');

% low-pass filter (broad band)
EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
% high-pass filter (broad band)
EEG = pop_firws(EEG, 'fcutoff', HP, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', HP_order);
% re-referencing to CAR
EEG = pop_reref(EEG, [], 'refstate',0);


% {'prep', 'rest', 'right_cue_execution', 'left_cue_execution', 'right_cue_imagery', 'left_cue_imagery','left_execution','right_execution','left_imagery','right_imagery'};

%% epoching for each event in the data 

% -erp of the preparation duration [-2s to -1s] 
% event markers 
events = 'prep';
epoch_start = -0.15;            
epoch_end = 1;
baseline = [epoch_start*EEG.srate 0];   
% make a copy of the data 
EEGprep = EEG;
% epoching 
EEGprep = pop_epoch(EEGprep, events, [epoch_start epoch_end], 'newname', 'MI_pilot_epoched','epochinfo', 'yes');
EEGprep = eeg_checkset(EEGprep);
% baseline correction
EEGprep = pop_rmbase(EEGprep, baseline);
EEGprep = eeg_checkset(EEGprep);
% create time vector
tprep = linspace(epoch_start, epoch_end, size(EEGprep.data,2));

% -erp of the cue duration [-1s to 0s] 
% event markers 
events = {'right_cue_execution', 'left_cue_execution', 'right_cue_imagery', 'left_cue_imagery'};
epoch_start = -0.15;            
epoch_end = 1;
baseline = [epoch_start*EEG.srate 0];   
% make a copy of the data 
EEGcue = EEG;
% epoching 
EEGcue = pop_epoch(EEGcue, events, [epoch_start epoch_end], 'newname', 'MI_pilot_epoched','epochinfo', 'yes');
EEGcue = eeg_checkset(EEGcue);
% baseline correction
EEGcue = pop_rmbase(EEGcue, baseline);
EEGcue = eeg_checkset(EEGcue);
% create time vector
tcue = linspace(epoch_start, epoch_end, size(EEGcue.data,2));

% -erp of rest period [-5s to -2s]
% event markers 
events = 'rest';
epoch_start = -0.15;            
epoch_end = 3;
baseline = [epoch_start*EEG.srate 0];   
% make a copy of the data 
EEGrest = EEG;
% epoching 
EEGrest = pop_epoch(EEGrest, events, [epoch_start epoch_end], 'newname', 'MI_pilot_epoched','epochinfo', 'yes');
EEGrest = eeg_checkset(EEGrest);
% baseline correction
EEGrest = pop_rmbase(EEGrest, baseline);
EEGrest = eeg_checkset(EEGrest);
% create time vector
trest = linspace(epoch_start, epoch_end, size(EEGrest.data,2));

% - erp of the execution/imagery period [0s to 4s]
% event markers 
events = {'left_execution','right_execution','left_imagery','right_imagery'};
epoch_start = -0.15;            
epoch_end = 4;
baseline = [epoch_start*EEG.srate 0];   
% make a copy of the data 
EEGimg = EEG;
% epoching 
EEGimg = pop_epoch(EEGimg, events, [epoch_start epoch_end], 'newname', 'MI_pilot_epoched','epochinfo', 'yes');
EEGimg = eeg_checkset(EEGimg);
% baseline correction
EEGimg = pop_rmbase(EEGimg, baseline);
EEGimg = eeg_checkset(EEGimg);
% create time vector
timg = linspace(epoch_start, epoch_end, size(EEGimg.data,2));


%% plotting ERPs 

% channel to plot 
ch = 41;

% plotting for rest period [-5s to -3s]
subplot(4,1,1)
plot(trest, mean(EEGrest.data(ch,:,:),3))
title('ERP for Rest Period [-5s to -3s]');
xlim([0 0.5])
ylim([-10 20])

% plotting for preparation period [-2s to -1s]
subplot(4,1,2)
plot(tprep, mean(EEGprep.data(ch,:,:),3))
title('ERP for Preparation Period [-2s to -1s]');
xlim([0 0.5])
ylim([-10 20])

% plotting for cue period [-1s to 0s]
subplot(4,1,3)
plot(tcue, mean(EEGcue.data(ch,:,:),3))
title('ERP for Cue Period [-1s to 0s]');
xlim([0 0.5])
ylim([-10 20])

% plotting for execution/imagey period [0s to 4s]
subplot(4,1,4)
plot(timg, mean(EEGimg.data(ch,:,:),3))
title('ERP for Execution/Imagery Period [0s to 4s]');
xlim([0 0.5])
ylim([-10 20])

%% plotting topoplots for each period 

% topoplots for rest ERP
pop_topoplot(EEGrest, 1, [0 100 150 200 250 300 350 400], '', [1 8] ,0, 'electrodes', 'on', 'chaninfo', EEG.chaninfo, 'maplimits', [-18 18]);
% topoplots for prep ERP
pop_topoplot(EEGprep, 1, [0 100 150 200 250 300 350 400], '', [1 8] ,0, 'electrodes', 'on', 'chaninfo', EEG.chaninfo, 'maplimits', [-18, 18]);
% topoplots for cue ERP
pop_topoplot(EEGcue, 1, [0 100 150 200 250 300 350 400], '', [1 8] ,0, 'electrodes', 'on', 'chaninfo', EEG.chaninfo, 'maplimits', [-18, 18]);
% topoplots for img ERP
pop_topoplot(EEGimg, 1, [0 100 150 200 250 300 350 400], '', [1 8] ,0, 'electrodes', 'on', 'chaninfo', EEG.chaninfo, 'maplimits', [-18, 18]);

    