clearvars; close all;  clc;
addpath 'M:\\Field Trip\\fieldtrip-20230118'; ft_defaults; % <---- UPDATE
addpath 'C:\\Users\\juadams\\Documents\\Field Trip\\fieldtrip-20230118'; ft_defaults; % <---- UPDATE
addpath(genpath('M:\\SNAP_scripts')); %<--- UPDATE
ctrA = 1; ctrB = 1; ctrC = 1; %counters to assign data below
linenoiseFreq = 50; %freq of the line noise to remove

% adjust this section according to your paths
%prepro_path = 'M:\\SNAP_scripts\\Julia_Ilvana_processing_script_comparisons_allgroups_10102023.m';
%scripts?
data_path = 'M:\\SNAP_Data_Processed\\processed'; % <---- UPDATE
cd(data_path);
list = dir('*.set');

% Specify the path to the Excel file containing subject IDs and conditions
excelFilePath = 'M:\\SNAP Groups\\SNAP_Groups.xlsx'; % <---- UPDATE

% Read the Excel file
groupData = readtable(excelFilePath, 'ReadVariableNames', false);

for q = 1:length(list)
    [~, fileName, ~] = fileparts(list(q).name);
    identifier = strsplit(fileName, '_'); % Change delimiter to '_'
    subID = identifier{1};

    % Convert numerical values to strings for comparison
    if isnumeric(groupData.Var1)
        strVar1 = cellfun(@num2str, num2cell(groupData.Var1), 'UniformOutput', false);
    else
        strVar1 = groupData.Var1;
    end

    rowIndex = find(strcmp(strtrim(strVar1), strtrim(subID)));

    if ~isempty(rowIndex)
        group_identifier = groupData.Var2{rowIndex}; % Change variable name to group_identifier
    else
        group_identifier = 'Unknown';
    end

    disp(['File: ' fileName ', Subject ID: ' subID ', Group Identifier: ' group_identifier]);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %clear EEG
    % use this to import files previously processed with fieldtrip
    %load(list(q).name);

    % use this to import EEGLAB .set files
    hdr = ft_read_header(list(q).name);
    data = ft_read_data(list(q).name, 'header', hdr);
    events = ft_read_event(list(q).name, 'header', hdr);

    %Preprocess the data
    cfg = [];
    cfg.dataset = list(q).name;
    subject = ft_preprocessing(cfg);

    %Find the AVERAGE across all trials
    cfg = [];
    subject = ft_timelockanalysis(cfg,subject);

    %APPLY LINE NOISE REMOVAL (cleanline)
    signal.data = subject.avg; %the data
    signal.srate = hdr.Fs; %sampling rate

    lineNoiseIn.fPassBand = [0, hdr.Fs/2]; %default = all
    lineNoiseIn.Fs = hdr.Fs; %sampling rate again
    lineNoiseIn.fScanBandWidth = 3; % +/- freq of interest to scan for sig lines
    lineNoiseIn.lineFrequencies = 50; % remove line noise from here
    lineNoiseIn.lineNoiseChannels = size(data, 1); %ALL, size of the dataset
    lineNoiseIn.maximumIterations = 50; %max iterations for removal
    lineNoiseIn.p = 0.05; % significance level cutoff
    lineNoiseIn.pad = 1; %padding (-1 = none, 0 is padding highest pwr2)
    lineNoiseIn.tapers = [5 9]; %standard
    lineNoiseIn.taperBandWidth = 2; % 2Hz
    lineNoiseIn.taperWindowSize = 2; %2 seconds window
    lineNoiseIn.taperWindowStep = 0.5; % overlap (if this matches above then no overlap)
    lineNoiseIn.tau = 100; %default for smoothing

    [clean_signal, lineNoiseOut] = cleanLineNoise(signal, lineNoiseIn); %run linenoise removal

    % NOW ADD THE 'LINE NOISE REMOVED' DATA BACK INTO THE STRUCTURE
    subject.avg = [];
    subject.avg = clean_signal.data;

    % Run FFT    
    % create 2 s segments with 1 s (50%) overlap
    cfg1 = [];
    cfg1.overlap = .5;
    cfg1.length  = 2;
    subject.base_rpt2sec    = ft_redefinetrial(cfg1, subject);

    %Multi-taper version of the original FFT method
    cfg2 = [];
    cfg2.output  = 'pow';
    cfg2.channel = 'all';
    cfg2.method  = 'mtmfft';
    cfg2.taper   = 'dpss';
    cfg2.tapsmofrq = 2;
    cfg2.foilim     = [0 100]
    cfg2.keeptrials ='no'; %this has to be NO to obtain group mean
    subject.base_freq2sec   = ft_freqanalysis(cfg2, subject.base_rpt2sec);

    %UNCOMMENT THIS TO SEE POWER SPECTRUM OF THIS PARTICIPANT
    % % % % code to plot the fft data, if interested
    % figure;
    % hold on;
    % plot(subject.base_freq2sec.freq, squeeze(mean((subject.base_freq2sec.powspctrm))))
    % legend('2 sec window')%,'2 sec window','4 sec window')
    % xlabel('Frequency (Hz)');
    % ylabel('absolute power (uV^2)');


    % %Plot this participants scalp map (uncomment if you want to visualise)
    % cfg = [];
    % %cfg.xlim = [0.9 1.3];
    % %cfg.ylim = [15 20];
    % %cfg.zlim = [-1e-27 1e-27];
    % cfg.baseline = 'no'; % [0 0.5];
    % cfg.baselinetype = 'absolute';
    % cfg.layout = subject
    % figure; ft_topoplotTFR(cfg,subject.base_freq2sec); colorbar
    %

    % Save results into different structures depending on what group this P belongs to
    disp 'saving...'
    save([subID '_' group_identifier '_cleanLine_rmline_pad.mat'],'subject');

    %Assign subject data to structure depending on their group
    if strcmp(group_identifier,'CC')
        group_data_CC{ctrA} = subject.base_freq2sec;
        ctrA = ctrA + 1;
    elseif strcmp(group_identifier,'HC')
        group_data_HC{ctrB} = subject.base_freq2sec;
        ctrB = ctrB + 1;
    elseif strcmp(group_identifier,'UHR')
        group_data_UHR{ctrC} = subject.base_freq2sec;
        ctrC = ctrC + 1;
    end

    aux.elec = subject.elec;
    clear subject filename new_data


end

%% Now we compute the group grand mean average

% Use as
%   [grandavg] = ft_freqgrandaverage(cfg, freq1, freq2, freq3...)
%
% The input data freq1..N are obtained from either FT_FREQANALYSIS with
% keeptrials=no or from FT_FREQDESCRIPTIVES. The configuration structure
% can contain
%   cfg.keepindividual = 'yes' or 'no' (default = 'no')
%   cfg.foilim         = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%   cfg.toilim         = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see FT_CHANNELSELECTION for details
%   cfg.parameter      = string or cell-array of strings indicating which
%                        parameter(s) to average. default is set to
%                        'powspctrm', if it is present in the data.

cfg.keepindividual = 'no'; %do not keep individual, average over participants within the group
cfg.foilim         = [30 100]; %frequency window of interest (change to gamma, 30 to 100)

[grandavg_HC] = ft_freqgrandaverage(cfg,group_data_HC{1:end}) %add up to number of participants
[grandavg_CC] = ft_freqgrandaverage(cfg,group_data_CC{1:end}) %add up to number of participants
[grandavg_UHR] = ft_freqgrandaverage(cfg,group_data_UHR{1:end}) %add up to number of participants

figure;
hold on;
plot(grandavg_HC.freq, squeeze(mean((grandavg_HC.powspctrm))), 'LineWidth',3)
plot(grandavg_CC.freq, squeeze(mean((grandavg_CC.powspctrm))), 'LineWidth',3)
plot(grandavg_UHR.freq, squeeze(mean((grandavg_UHR.powspctrm))), 'LineWidth',3)
legend({'HC','CC','UHR'})
xlabel('Frequency (Hz)');
ylabel('absolute power (uV^2)');

%Setup the grand average results to contain the electrode location data
grandavg_HC.elec = aux.elec ; %assign the electorde arrangement to this dataset (for plotting)
grandavg_CC.elec = aux.elec ; %assign the electorde arrangement to this dataset (for plotting)
grandavg_UHR.elec = aux.elec ; %assign the electorde arrangement to this dataset (for plotting)


%Plot scalp maps of the grand average power for each group
cfg = [];
%cfg.xlim = [0.9 1.3];
%cfg.ylim = [15 20];
cfg.zlim = [0.01 0.03];
cfg.baseline = 'no' ;%[0 0.5];

figure(1); ft_topoplotTFR(cfg,grandavg_HC); colorbar; sgtitle('HC'); hold on;
figure(2); ft_topoplotTFR(cfg,grandavg_CC); colorbar; sgtitle('CC'); hold on;
figure(3); ft_topoplotTFR(cfg,grandavg_UHR); colorbar; sgtitle('UHR'); hold on;