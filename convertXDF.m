% CONVERT XDF to EEGLAB .SET file
% -------------------------------
% The scripts loads all the .xdf files from the filepath, loads them to
% EEGLAB using loadxdf function and saves it as .set files to the folder.
% 
% author:   Abin Jacob 
%           Carl von Ossietzky Universität Oldenburg
%           abin.jacob@uni-oldenburg.de
% date  :   17/05/2024

clear; clc; close all;


% --- Files & Folders ---
% folder name
foldername = 'Pilot08_earEEG_21062024';

% path to the rawdata
rootpath = '/Users/abinjacob/Documents/01. Calypso/Calpso 1.0/EEG Data';
filepath = fullfile(rootpath,foldername);
foldsplit = strsplit(foldername, '_');
p = foldsplit{1}; dt = 'rawdata';

% read all .xdf files from DATAPATH
filelist = dir(fullfile(filepath, '*.xdf'));
% opeing EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


% loop over files 
for file = 1:length(filelist)
    % select file to load
    file2load = fullfile(filepath,filelist(file).name);
    namesplit = strsplit(filelist(file).name, '-');
    exp = strrep(namesplit{3}, '_eeg.xdf', '');
    % load LSL files (.xdf)
    EEG = pop_loadxdf(file2load, 'streamtype', 'EEG', 'exclude_markerstreams', {});
    setname = [p, '_', exp, '_', dt];
    EEG.setname = setname;
    % add comments for the EEG dataset
    EEG.comments = '';
    % save the dataset to filepath 
    EEG = pop_saveset(EEG, [EEG.setname, '.set'], filepath);
end 


