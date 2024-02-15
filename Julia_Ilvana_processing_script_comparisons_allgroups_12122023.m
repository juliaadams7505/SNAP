
% Resting state preprocessing pipeline
% Created by Julia Adams 06/09/2023
% Adapted for research computer use 17/10/2023

addpath('M:\\Troublesome_SNAP_data');

clear all
eeglab 

% Input of raw data path
%raw_path = 'C:/Users/juadams/Downloads/eeglab_current (1)/SNAP DATA/raw_rest';
%raw_path = 'C:/Users/juadams/OneDrive - ORYGEN/Documents/SNAP data';
raw_path = 'M:\\SNAP_data';

% Output of processed data path
%data_path = 'I:/Research/GAMMA/6. EEG_Data/SNAP_UHR_CC_HC_Comparisons/Processed';
%data_path = 'C:/Users/juadams/OneDrive - ORYGEN/Documents/Processed';
data_path = 'M:\\SNAP_Data_Processed';

cd(raw_path);
foldername = dir;

for i = 3: length(foldername)
name = foldername(i).name;
cd (name)
list(i) = dir('*Closed*.cdt');
cd ('..')
end


file_ext = '.cdt';


% Load the Curry EEG dataset (update with the actual dataset file)
% EEG = pop_loadcurry(fullfile(raw_path, 'raw_data'), 'ignorewarnings', 'on'););
% EEG = pop_loadcurry(fullfile(raw_path,'155_RestingEyesClosed_17Dec2021.cdt'), 'ignorewarnings', 'on'); %

%% make output folders:
%make folder for post-ICA, uninterpolated files 
if ~isdir ([data_path filesep 'cleanLineNoise'])
    mkdir ([data_path filesep 'cleanLineNoise']);
end

%make folder for post-ICA, uninterpolated files 
if ~isdir ([data_path filesep 'intermediate1_ICAclean'])
    mkdir ([data_path filesep 'intermediate1_ICAclean']);
end

%make folder for final preprocessed files
if ~isdir ([data_path filesep 'processed'])
    mkdir ([data_path filesep 'processed']);
end

%make folder for final preprocessed files
if ~isdir ([data_path filesep 'testprocessed'])
    mkdir ([data_path filesep 'testprocessed']);
end


%% One issue with trying to code this: inconsistency in file naming. 
% Sometimes 'EyesClosed', sometimes 'RestingEyesClosed'. Always ends in
% fixed. Does this matter? 

%% Batch Processing

% intialize report metrics
FileNames={zeros(1,length(list))};
Number_Good_Channels_Selected=zeros(1,length(list));
Percent_Good_Channels_Selected=zeros(1,length(list));
Interpolated_Channel_IDs=[];
File_Length_In_Secs=zeros(1,length(list));
Percent_File_Length=zeros(1,length(list));
Number_ICs_Rejected=zeros(1,length(list));
Percent_ICs_Rejected=zeros(1,length(list));
Number_Segments_Post_Segment_Rejection=zeros(1,length(list));

for l=3:length(list)

 loadName = list(l).name;
    dataName = loadName(1:end);
    FileNames(l) = {dataName};    

     %Display participant number
    disp('******************')
    disp(['Processing ' num2str(dataName) '...'])
    disp('******************')
    
    %Import data - does not work because on this computer the filesep is \ instead of /
    LoadName_Filename = fullfile(raw_path, foldername(l).name, loadName);
    EEG = loadcurry(LoadName_Filename, 'KeepTriggerChannel', 'True', 'CurryLocations', 'False');

    % place processing skip here
    EEG.setname = dataName;

    
    %remove trigger channel
    EEG = pop_select(EEG, 'nochannel', 67);
      %Event select
       EEG = pop_importevent(EEG, 'append','no',...
           'event','M:\\Event_EEG_markers.txt',...
           'fields',{'number','type','latency','urevent'},'skipline', 1,'timeunit', 1E-3);
  % Is the above just importing and not loading? pop_editeventfield =
  % load/remove fields in a dataset. 
  	% [EEG,eventnumbers] = pop_editeventfield( EEG, 'key1', 'value1', ...);
       EEG = pop_rmdat( EEG, {'9'}, [0 180], 0);
  
     %Remove baseline
    EEG = pop_rmbase( EEG, [],[]);

      %Channel edit
    EEG= pop_chanedit(EEG, 'settype',{'1:64','EEG'});

     %resample to 1000Hz
    EEG = pop_resample( EEG, 1000);

    % Filter 1-100Hz
    pop_eegfiltnew(EEG,'locutoff',1,'hicutoff',100);
    EEG = eeg_checkset(EEG);

     % store original information      
    orig_EEG_chan = EEG.nbchan;
    orig_EEG_dur = EEG.pnts/EEG.srate;
    orig_EEG_chanloc_labels = {EEG.chanlocs.labels};
    orig_EEG_chanlocs = EEG.chanlocs;

    %  %clean line noise
    % signal      = struct('data', EEG.data, 'srate', EEG.srate);
    % lineNoiseIn = struct('lineNoiseMethod', 'clean', ...
    %                      'lineNoiseChannels', 1:EEG.nbchan,...
    %                      'Fs', EEG.srate, ...
    %                      'lineFrequencies', [50 100 150 200],...
    %                      'p', 0.01, ...
    %                      'fScanBandWidth', 2, ...
    %                      'taperBandWidth', 2, ...
    %                      'taperWindowSize', 4, ...
    %                      'taperWindowStep', 1, ...
    %                      'tau', 100, ...
    %                      'pad', 2, ...
    %                      'fPassBand', [0 EEG.srate/2], ...
    %                      'maximumIterations', 10);
    % [clnOutput, lineNoiseOut] = cleanLineNoise(signal, lineNoiseIn);
    % EEG.data = clnOutput.data;
    % 
    % %save intermediate
    % EEG = pop_saveset(EEG, 'filename',strrep(dataName, file_ext,'_cleanLineNoise.set'),'filepath',[data_path filesep 'cleanLineNoise']);
    % 
%%% from EEGLab tutorials
% % Identify bad channels
% bad_channel_indices = [3, 7, 10];
% 
% % Remove bad channels
% EEG = pop_select(EEG, 'nochannel', bad_channel_indices);

 % Artifact rejection HH version
    raw_EEG = EEG;
    EEG = clean_artifacts(EEG,'burst_crit', 25, 'line_crit', 4, 'chancorr_crit',.8);
    vis_artifacts(EEG,raw_EEG)
    
    % Artifact rejection Melbourne version
    %raw_EEG = EEG;
    %EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
    %vis_artifacts(EEG,raw_EEG)
    %eeglab redraw

    % read Artifact rejection results
    Number_Good_Channels_Selected(l) = EEG.nbchan;
    Percent_Good_Channels_Selected(l)= EEG.nbchan/orig_EEG_chan*100;
    bad_channels_removed = setdiff(orig_EEG_chanloc_labels, {EEG.chanlocs.labels});
    if isempty(bad_channels_removed)
        Interpolated_Channel_IDs{l} = 'none';
    else
        Interpolated_Channel_IDs{l}=[sprintf('%s ',bad_channels_removed{1:end-1}),bad_channels_removed{end}];
    end
    File_Length_In_Secs(l) = EEG.pnts/EEG.srate;
    Percent_File_Length(l) = File_Length_In_Secs(l)/orig_EEG_dur*100;

 % Apply average reference after adding initial reference
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, [], 'exclude', [60 64 65 66]);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
    EEG = eeg_checkset(EEG);

% Data matric precision (needed?) 
% EEG.data = double(EEG.data);


% Run Independent Component Analysis (ICA)
EEG_forICA = pop_resample(EEG, 100);
    EEG_forICA = pop_runica(EEG_forICA, 'extended',1,'interupt','off');
    EEG.icaweights = EEG_forICA.icaweights;
    EEG.icasphere  = EEG_forICA.icasphere;
    EEG = eeg_checkset(EEG, 'ica');

     % Perform IC rejection using ICLabel scores and r.v. from dipole fitting.
    EEG = pop_iclabel(EEG, 'default');
    
    %save intermediate
    EEG = pop_saveset(EEG, 'filename',strrep(dataName, file_ext,'_ICA.set'),'filepath',[data_path filesep 'intermediate1_ICAclean']);
    
    % Obtain the most dominant class label and its label probability.
    [~, mostDominantClassLabelVector] = max(EEG.etc.ic_classification.ICLabel.classifications, [], 2);
    mostDominantClassLabelProbVector = zeros(length(mostDominantClassLabelVector),1);
    for icIdx = 1:length(mostDominantClassLabelVector)
             mostDominantClassLabelProbVector(icIdx)  = EEG.etc.ic_classification.ICLabel.classifications(icIdx, mostDominantClassLabelVector(icIdx));
    end 
    % Identify brain ICs. The order of the classes are {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}.
    brainLabelProbThresh  = 0.5; % [0-1]
    brainIdx = find((mostDominantClassLabelVector==1 & mostDominantClassLabelProbVector>=brainLabelProbThresh) |(mostDominantClassLabelVector==7 & mostDominantClassLabelProbVector>=brainLabelProbThresh));
    % Perform IC rejection using residual variance of the IC scalp maps.
    %rvList    = [EEG.dipfit.model.rv];
    %goodRvIdx = find(rvList < 0.15)'; % < 15% residual variance == good ICs.
 
    % Perform IC rejection using inside brain criterion.
    %load(EEG.dipfit.hdmfile); % This returns 'vol'.
    %dipoleXyz = zeros(length(EEG.dipfit.model),3);
    %for icIdx = 1:length(EEG.dipfit.model)
        %dipoleXyz(icIdx,:) = EEG.dipfit.model(icIdx).posxyz(1,:);
    %end
    %depth = ft_sourcedepth(dipoleXyz, vol);
    %depthThreshold = 1;
    %insideBrainIdx = find(depth<=depthThreshold);
 
    % Take AND across the three criteria.
    goodIcIdx = brainIdx;
    %goodIcIdx = intersect(brainIdx, goodRvIdx);
    %goodIcIdx = intersect(goodIcIdx, insideBrainIdx);
    Number_ICs_Rejected(l)=length(EEG.icaweights)-length(goodIcIdx);
    Percent_ICs_Rejected(l)=(length(EEG.icaweights)-length(goodIcIdx))/length(EEG.icaweights);
%    pv_afterICA(l)=eeg_pvaf(EEG,brainIdx'); %ilvana removed due to error
    pvaf=eeg_pvaf(EEG);
%    pv_ICA(l)=pvaf(end); %ilvana removed due to error
     
    raw_EEG = EEG;
    % Perform IC rejection.
    EEG = pop_subcomp(EEG, goodIcIdx, 0, 1);
    vis_artifacts(EEG,raw_EEG, 'OldColor',[0 1 0] )

     EEG.icaact = [];
    EEG = eeg_checkset(EEG, 'ica');
    
    % Interpolate channels.
    EEG = pop_interp(EEG, orig_EEG_chanlocs, 'spherical');
    EEG = eeg_checkset(EEG);

    %Remove channels
    EEG = pop_select( EEG, 'nochannel',{'VEO','HEO','CB1','CB2'});

     %Epoch data into 2 sec epochs
    EEG = eeg_regepochs (EEG, 'recurrence', 2, 'limits', [0 2], 'eventtype', '10');


    % Save cleaned data
    EEG = pop_saveset(EEG, 'filename',strrep(dataName, file_ext,'_processed.set'),'filepath',[data_path filesep 'processed']);
    %else 
    %    disp('nothing to do here');
    %end
end

%%% plots from EEGlab tute
% % Plot the spectra of the independent components
% pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freq', [2 29], 'freqrange', [1 50], 'electrodes', 'off');
% pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freq', [30 45], 'freqrange', [1 100], 'electrodes', 'off'); % what wider freqrange should be considered?
% pop_spectopo(EEG, 1, [], 'EEG', 'percent', 100, 'freq', [55 100], 'freqrange', [1 100], 'electrodes', 'off'); % default was 1 50 
% % avoiding 50 Hz thing with weird electricity 

% Save the cleaned EEG data
%pop_saveset(EEG, 'SUBJECTID_uhr_SNAP_final_processed.set');
% pop_saveset(EEG, '155_RestingEyesClosed_17Dec2021_processed_test.set');

% Things I am still playing with - how to get the code to follow the
% filepath pattern 

outputtable=table(FileNames',File_Length_In_Secs',Number_Segments_Post_Segment_Rejection',...
    Number_Good_Channels_Selected', Percent_Good_Channels_Selected', Interpolated_Channel_IDs',Number_ICs_Rejected',...
    Percent_ICs_Rejected');
outputtable.Properties.VariableNames ={'FileNames','File_Length_In_Secs','Number_Segments_Post_Segment_Rejection',...
    'Number_Good_Channels_Selected', 'Percent_Good_Channels_Selected', 'Interpolated_Channel_IDs','Number_ICs_Rejected',...
    'Percent_ICs_Rejected'};

cd([data_path filesep 'processed']);
writetable(outputtable, ['All_subs_output_table_ICA ',datestr(now,'dd-mm-yyyy'),'.csv']);





