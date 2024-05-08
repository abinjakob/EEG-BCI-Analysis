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
% load file in EEGLAB

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
epoch_start = -3;            
epoch_end = 5;               
% baseline correction
% defining baseline for baseline correcion
baseline = [epoch_start*EEG.srate 0];  
PRUNE = 4;

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

titlestr = {'for Left 20Hz', 'for Left 15Hz', 'for Right 20Hz', 'for Right 15Hz'};
savestr = {'L20', 'L15', 'R20', 'R15'};
% channels to plot tf 
% chan2plot = [56 55 54 06 41 07 36 35 01 32 33 03 38 04 46 47 48 08 20];
chan2plot = [77 75];

% % for Madina's data 
% titlestr = {'for Right ME', 'for Right MI'};
% savestr = {'RightEx', 'RightIm'};
% % channels to plot tf 
% chan2plot = [01 02 03 05 06 31 28];

% filepath to save the tf images 
imagepath = '/Users/abinjacob/Documents/02. NeuroCFN/Research Module/RM02/Figures/rmweek13figures/SSVEP_TF/';

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
            EEG.chaninfo, 'caption', [EEG.chanlocs(ch).labels, ' ', titlestr{ievent}], 'plotphase', 'off', 'padratio', 16, 'winsize', 500);
        % set image name 
        imagename = [savestr{ievent}, '_', EEG.chanlocs(ch).labels];
        % save image as png
        saveas(gcf, [imagepath, imagename],'png');
        close;
    end 
end 



