% EEG PERFORMING ICA 
% ------------------
% The code performs ICA (EEGLAB) and applies the weights on the data.
% Before computing the ICA the EEG data is prepared for ICA with the
% following pre-procesing steps:
% - 1Hz high-pass filter
% - 40Hz low-pass filter
% - downsampling to 250Hz
% - epoching (1s epochs)
% - removing artefactual epochs
% - dimensionality reduction using PCA (optional) 
% 
% Note: The data has to be converted into .set file with appropriate
%       channel locations. 
%
% Author: Abin Jacob
%         Carl von Ossietzky University Oldenburg
%         abin.jacob@uni-oldenburg.de            
% Date  : 16/04/2024

%% start fresh 
% clear; clc; close all;
addpath('L:\Cloud\SW\eeglab2024.0');
cd('L:\Cloud\NeuroCFN\RESEARCH PROJECT\Research Project 02\EEG Analysis')

%% load file
% open EEGLab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% folder path
DATAPATH = 'L:\Cloud\NeuroCFN\RESEARCH PROJECT\Research Project 02\EEG Analysis\data_preprocessed';
% filename of the .set file with chanlocs
filename = 'P04_MI_rawdata.set';
% load the file in EEGLAB
EEG = pop_loadset('filename', filename, 'filepath', DATAPATH);

%% set parameters for ICA
% high-pass filter 
HP = 1;                   % cut-off
HP_order = 500;           % filter order    
% low-pass filter  
LP = 40;                  % cut-off
LP_order = 100;           % filter order 
% downsampling freq. for ICA
SRATE = 250;
% artifact rejection threshold on SD
PRUNE = 3;
% perform PCA before ICA for dimension reduction [0 : 'No', 1: 'Yes']
PCA = 1;
% PCA dimension 
PCADIMS = 50;

%% prepare data for ICA
% apply low-pass filter 
EEG = pop_firws(EEG, 'fcutoff', LP, 'ftype', 'lowpass', 'wtype', 'hann', 'forder', LP_order);
% downsample the data
EEG = pop_resample(EEG, SRATE);
% apply high-pass filter
EEG = pop_firws(EEG, 'fcutoff', HP, 'ftype', 'highpass', 'wtype', 'hann', 'forder', HP_order);
% create epochs for IC
EEG = eeg_regepochs(EEG, 'recurrence', 1, 'eventtype', '999');
EEG = eeg_checkset(EEG, 'eventconsistency');
% remove epochs with artefacts above standard deviation PRUNE 
EEG = pop_jointprob(EEG, 1, [1:size(EEG.data,1)], PRUNE, PRUNE, 0, 1, 0);

%% performing ICA
% check if PCA is opted
if PCA == 1
    % perform ICA after reducing dimension using PCA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', PCADIMS);
else
    % perform ICA with current dimensions
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
end

% store ICA weights
icawinv = EEG.icawinv;
icas = EEG.icasphere;
icaw = EEG.icaweights;

%% apply ICA weights to the data 
% load the file in EEGLAB
EEG = pop_loadset('filename', filename, 'filepath', DATAPATH);

% apply ICA weights
EEG.icawinv = icawinv;
EEG.icasphere = icas;
EEG.icaweights = icaw;
EEG = eeg_checkset(EEG);
