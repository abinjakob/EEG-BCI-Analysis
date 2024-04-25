% Motor Imagery CSP Filtering  
% ---------------------------
% The script is used to apply CSP filter to the Motor Imagery EEG data. The artefact
% corrected EEG data file is loaded and preprocessed. Further the data is
% subepoched into left and trials and CSP is computed 
%
% Uses the function CSP_new.m written by 
%
% Author: Abin Jacob
%         Carl von Ossietzky University Oldenburg
%         abin.jacob@uni-oldenburg.de            
% Date  : 24/04/2024

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
HP = 8;                         % cut-off
HP_order = 826;                 % filter order    
% low-pass filter  
LP = 12;                        % cut-off
LP_order = 776;                 % filter order 

% epoching
% event markers 
events = {'left_execution','right_execution','left_imagery','right_imagery'};
epoch_start = 0.5;            
epoch_end = 4.5;

% baseline correction
% defining baseline for baseline correcion
baseline = [epoch_start*EEG.srate 0]; 
% reject artefactual epochs 
PRUNE = 4;

% no. of CSP components
ncomp = 10; 

%% preprocessing 

% low-pass filter (broad band)
EEG = pop_eegfiltnew(EEG, [],HP,[],1,[],0);
EEG = pop_eegfiltnew(EEG, [],LP,[],0,[],0);

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

EEGL = pop_selectevent(EEG, 'type', events{1},'renametype', events{1}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
EEGR = pop_selectevent(EEG, 'type', events{2},'renametype', events{2}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');

%% computing CSP
cd('L:\Cloud\NeuroCFN\RESEARCH PROJECT\Research Project 02\EEG Analysis\CSP Matlab scripts shared by Conny\Matlab scripts');
% reshaping the data for CSP (nchans, ntpts*ntrils)
dataL = reshape(EEGL.data, size(EEGL.data,1), size(EEGL.data,2) * size(EEGL.data,3));
dataR = reshape(EEGR.data, size(EEGR.data,1), size(EEGR.data,2) * size(EEGR.data,3));

% computing CSP using the csp_new script from the NeuroCFN 
CSPfilt = csp_new(dataL, dataR);
CSPfilt = double(CSPfilt);
% csp patterns (pseudo-inverse of CSPfilt matrix)
CSPpattern = pinv(CSPfilt);

% eliminate emüty CSP
CSPfilt = CSPfilt(sum(abs(CSPfilt),2)>0,:);
CSPpattern = CSPpattern(:, sum(abs(CSPpattern),2)>0);
% totoal CSP components 
tcomp = size(CSPfilt, 1);

% indices of CSP filters to consider
csp2consider = [1:ncomp, (tcomp+1-ncomp):tcomp]; 
aux_evalcsp(EEGL.data, EEGR.data, {CSPfilt, CSPpattern}, EEG.chanlocs, ... 
EEG.srate, csp2consider, 2 *(epoch_end-epoch_start), {'left', 'right'});

