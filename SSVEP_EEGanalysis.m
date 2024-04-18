% SSVEP EEG Analysis Script 
% -------------------------
% The script is used to analyse the SSVEP EEG data. The artefact corrected 
% EEG data file is loaded and preprocessed. Further the data is epoched and 
% power spectral density is computed using pwelch method. 
%
% Author: Abin Jacob
%         Carl von Ossietzky University Oldenburg
%         abin.jacob@uni-oldenburg.de            
% Date  : 18/04/2024

%% start fresh 
% clear; clc; close all;

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
filename = 'P04_SSVEP_ICAcleaned.set';
% load the file in EEGLAB
EEG = pop_loadset('filename', filename, 'filepath', DATAPATH);

%% parameters for the analysis

% filtering 
% high-pass filter 
HP = 0.1;                       % cut-off
HP_order = 826;                 % filter order    
% low-pass filter  
LP = 45;                        % cut-off
LP_order = 776;                 % filter order 

% epoching
% event markers 
events = {'stim_L20','stim_L15','stim_R20','stim_R15'};
epoch_start = -0.2;            
epoch_end = 4;               
% baseline correction
% defining baseline for baseline correcion
baseline = [epoch_start*EEG.srate 0];   

% reject artefactual epochs 
PRUNE = 4;

%% pre-processing 

% low-pass filter
EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
% high-pass filter
EEG = pop_firws(EEG, 'fcutoff', HP, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', HP_order);
% re-referencing to CAR
EEG = pop_reref(EEG, [], 'refstate',0);

% removing unnecessary event marker
event_pos = 1;      % position counter for the events other than stim onset
event_idx = [];     % array to store the index of the event other than stim onset
% loop over events 
for idx = 1: length(EEG.event)
    if ~ strcmp(EEG.event(idx).type, events)
        event_idx(event_pos) = idx;
        event_pos = event_pos +1;
    end
end 
% remove events which are not stim onset from the data
EEG = pop_editeventvals(EEG, 'delete', event_idx);
EEG = eeg_checkset(EEG);

% epoching 
EEG = pop_epoch(EEG, events, [epoch_start epoch_end], 'newname', 'MI_pilot_epoched','epochinfo', 'yes');
% reject artefactual epochs 
% joint probability-based artifact rejection (joint prob. > PRUNE (SD))
EEG = pop_jointprob(EEG, 1, [1:EEG.nbchan], PRUNE, PRUNE, 0, 1, 0);
EEG = eeg_checkset(EEG);
% baseline correction
EEG = pop_rmbase(EEG, baseline);
EEG = eeg_checkset(EEG);

%% epcohing and calculating PSD for each condition

% loop over events 
for iEvent = 1:length(events)
    % epoching  
    EEG_temp = pop_selectevent(EEG, 'type', events{iEvent},'renametype', events{iEvent}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
    
    % window size for pwelch 
    window_length = EEG_temp.pnts;

    % calculating psd and storing
    [event(iEvent).psd, event(iEvent).f] = calc_psd(EEG_temp, window_length);
    % setting event type
    event(iEvent).eventtype = events{iEvent};
end 




