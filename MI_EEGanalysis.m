% Motor Imagery EEG Analysis Script 
% ----------------------------------
% The script is used to analyse the Motor Imagery EEG data. The artefact
% corrected EEG data file is loaded and preprocessed. Further the data is
% epoched and ERD is calculated. 
%
% Author: Abin Jacob
%         Carl von Ossietzky University Oldenburg
%         abin.jacob@uni-oldenburg.de            
% Date  : 17/04/2024

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
base_start = -4;
base_end = -2;

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

%% sub-epcohing and calculating ERD for each condition in Mu Band

% narrow-band filtering (muband)
EEG_mu = pop_firws(EEG, 'fcutoff', mu(2), 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
EEG_mu = pop_firws(EEG_mu, 'fcutoff', mu(1), 'ftype', 'highpass', 'wtype', 'hamming', 'forder', HP_order);

% loop over events 
for iEvent = 1:length(events)    
    EEG_temp = pop_selectevent(EEG_mu, 'type', events{iEvent},'renametype', events{iEvent}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
    muBand(iEvent).erd = computeERD(EEG_temp, binsize, base_start, base_end, epoch_start, epoch_end);
    muBand(iEvent).eventtype = events{iEvent};
end 

%% sub-epcohing and calculating ERD for each condition in Beta Band

% narrow-band filtering (muband)
EEG_beta = pop_firws(EEG, 'fcutoff', beta(2), 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', LP_order);
EEG_beta = pop_firws(EEG_beta, 'fcutoff', beta(1), 'ftype', 'highpass', 'wtype', 'hamming', 'forder', HP_order);

% loop over events 
for iEvent = 1:length(events)    
    EEG_temp = pop_selectevent(EEG_beta, 'type', events{iEvent},'renametype', events{iEvent}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
    betaBand(iEvent).erd = computeERD(EEG_temp, binsize, base_start, base_end, epoch_start, epoch_end);
    betaBand(iEvent).eventtype = events{iEvent};
end

%% plotting  

% creating time points for plotting 
timevec = linspace(epoch_start,epoch_end,size(muBand(1).erd,2));
chan2plot = [56 55 41 7; 46 47 38 3];

starttime = 0;
endtime = 2;

% find the strat and end index of baseline period
epochtime = epoch_start:epoch_end;
id1 = ceil(((find(epochtime == starttime)-1) * EEG.srate) / binsize);
id2 = ceil(((find(epochtime == endtime)-1) * EEG.srate) / binsize);

% plotting ERD in Mu band
% plotting ME
plot_title = {'Left hand Motor Execution', 'Right hand Motor Execution'};
j = 0;
figure;
% loop over events 
for i = 1:2 
    subplot(2,2,j+i)
    % plotting for C3
    plot(timevec, mean(muBand(i).erd(chan2plot(1,:),:),1));
    hold on 
    % plotting for C4
    plot(timevec, mean(muBand(i).erd(chan2plot(2,:),:),1));
    title(plot_title{i})
    xlabel('time (sec)');
    ylabel('ERD (%)');
    xlim([-5 6])
    ylim([-100 400]);
    legend('Left Electrodes', 'Right Electrodes')
    hold off   
    % plot barplots for the mean ERD (500ms to 3500ms)
    j = j + 1;
    subplot(2,2,j+i)
    meanC3 = mean(mean(muBand(i).erd(chan2plot(1,:), id1:id2),2),1);
    meanC4 = mean(mean(muBand(i).erd(chan2plot(2,:), id1:id2),2),1);
    % Create a bar plot
    bar([meanC3, meanC4],'k');
    % Setting the category names for the x-axis
    category = {'C3', 'C4'};
    set(gca, 'XTick', 1:2, 'XTickLabel', category);
%     ylim([-60 0]);
    % Adding labels and title for clarity
    ylabel('ERD (%)');
    title(plot_title{i})
end 

% plotting MI
plot_title = {'Left hand Motor Imagery', 'Right hand Motor Imagery'};
j = 0;
figure;
% loop over events 
for i = 1:2    
    subplot(2,2,j+i)
    % plotting for C3
    plot(timevec, mean(muBand(i+2).erd(chan2plot(1,:),:),1));
    hold on 
    % plotting for C4
    plot(timevec, mean(muBand(i+2).erd(chan2plot(2,:),:),1));
    title(plot_title{i})
    xlabel('time (sec)');
    ylabel('ERD (%)');
    xlim([-5 6])
%     ylim([-100 400]);
    legend('C3', 'C4')
    hold off  
    % plot barplots for the mean ERD (500ms to 3500ms)
    j = j + 1;
    subplot(2,2,j+i)
    meanC3 = mean(mean(muBand(i+2).erd(chan2plot(1,:), id1:id2),2),1);
    meanC4 = mean(mean(muBand(i+2).erd(chan2plot(2,:), id1:id2),2),1);
    % Create a bar plot
    bar([meanC3, meanC4],'k');
    % Setting the category names for the x-axis
    category = {'C3', 'C4'};
    set(gca, 'XTick', 1:2, 'XTickLabel', category);
%     ylim([-70 0]);
    % Adding labels and title for clarity
    ylabel('ERD (%)');
    title(plot_title{i})
end 

%% 
event2plot = 3;
EEG_temp = pop_selectevent(EEG_mu, 'type', events{event2plot},'renametype', events{event2plot}, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
pop_topoplot(EEG_temp, 1, [0 100 200 300 400 500 600 700 800 900 1000 1200 1300 1400 1500 1600], '', [4 4] ,0, 'electrodes', 'on', 'chaninfo', EEG.chaninfo);






