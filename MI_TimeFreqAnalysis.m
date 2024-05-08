% Motor Imagery Time-Frequency Analysis 
% -------------------------------------
% The script is does a Morlet Wavelet Transform of the MI EEG signal to
% understand the frequency change across time. 
%
% Author: Abin Jacob
%         Carl von Ossietzky University Oldenburg
%         abin.jacob@uni-oldenburg.de            
% Date  : 29/04/2024

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

% epoching
% event markers 
events = {'left_execution','right_execution','left_imagery','right_imagery'};
epoch_start = -5;            
epoch_end = 6;

% baseline correction
% defining baseline for baseline correcion
baseline = [epoch_start*EEG.srate 0];   

% reject artefactual epochs 
PRUNE = 4;

% parameters for ERD calculation
mu = [8 12];
beta = [13 30];
binsize = 50;
base_start = -2;
base_end = -1;

% % for Madina's data 
% events = {'S  1','S  3'};
% epoch_start = -5;            
% epoch_end = 7;
% base_start = -2;
% base_end = -1;
% PRUNE = 4;
% baseline = [epoch_start*EEG.srate 0];  

%% pre-processing 

% low-pass filter (broad band)
EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
% high-pass filter (broad band)
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

%% time-frequency analysis
% plots time-freq plots for each condition for the specified channels

titlestr = {'for Left ME', 'for Right ME', 'for Left MI', 'for Right MI'};
savestr = {'LeftEx', 'RightEx', 'LeftIm', 'RightIm'};
% channels to plot tf 
chan2plot = [01 03 04 06 07 38 41];

% % for Madina's data 
% titlestr = {'for Right ME', 'for Right MI'};
% savestr = {'RightEx', 'RightIm'};
% % channels to plot tf 
% chan2plot = [01 02 03 05 06 31 28];

% filepath to save the tf images 
imagepath = '/Users/abinjacob/Documents/02. NeuroCFN/Research Module/RM02/Figures/rmweek13figures/P04_TF/';

% loop over events 
for ievent = 1:length(events)    
    EEG_new = pop_selectevent(EEG, 'type', events{ievent},'renametype', events{ievent}, 'deleteevents', ...
        'off', 'deleteepochs', 'on', 'invertepochs', 'off');
    % loop over channels to plot
    for ichan = 1:length(chan2plot)  
        ch = chan2plot(ichan);
        % new figure for each channel 
        figure;
        % plotting tf and itc 
        pop_newtimef( EEG_new, 1, ch, [EEG.xmin EEG.xmax]*1000, [3 0.8] , 'topovec', ch, 'elocs', EEG.chanlocs, 'chaninfo', ...
            EEG.chaninfo, 'caption', [EEG.chanlocs(ch).labels, ' ', titlestr{ievent}], 'plotphase', 'off', 'padratio', 16, 'winsize', EEG.srate);
        % set image name 
        imagename = [savestr{ievent}, '_', EEG.chanlocs(ch).labels];
        % save image as png
        saveas(gcf, [imagepath, imagename],'png');
        close;
    end 
end 



