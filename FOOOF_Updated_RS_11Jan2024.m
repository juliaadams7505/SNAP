clearvars
%add_paths;
%eeglab;

addpath(genpath('C:\Users\juadams\Downloads\fooof_mat-main\fooof_mat-main'));%make sure subfolders are included

env = 'fooof1'; % Replace with your actual environment name
conda_path = ['C:\Users\juadams\Documents\Anaconda3.10\envs\' env '\python.exe'];
pyenv("Version",conda_path)

%Import necessary modules
py.importlib.import_module('numpy') %first time from anaconda fooof1 terminal: conda install numpy
py.importlib.import_module('fooof') %first time from anaconda fooof1 terminal: conda install -c conda-forge fooof

% Data directories
eegpath = 'M:\SNAP_Data_Processed\processed';

fooof_dest_path = 'M:\SNAP_Data_Processed\foof_files_default';
chanlocs_dest_path = 'M:\SNAP_Data_Processed\chanlocs_default';

%% make output folders:
if ~isfolder (fooof_dest_path)
    mkdir (fooof_dest_path);
end

if ~isfolder (chanlocs_dest_path)
    mkdir (chanlocs_dest_path);
end
conditions_overall={'CC','HC','UHR','UNK'};
%upload subject data
% Read the Excel file
% filename = 'SNAP_groups.xlsx';  % Update with your actual filename
% [num, txt, raw] = xlsread(filename);
%
% % Assuming Participant ID is in the first column and Condition is in the second column
% participant_id_column = 1;
% condition_column = 2;
%
% % Extracting columns
% participant_ids = raw(:, participant_id_column);
% conditions = raw(:, condition_column);



%% Variable initialization
conditions={'hc','cc','uhr'};
%This is the output struct
allfooof={};
rlats={};


%FOOF settings, these are basically the defaults

settings = struct(...
    'peak_width_limits', [1, 6], ...,
    'peak_threshold', 2.0, ...
    'background_mode', 'fixed', ...
    'verbose', false);

freqrange=[1 45]; %Data above 40 Hz less studied with FOOOF

bands = {[1 3], [4 7], [8 13], [14 30], [31 45]}; %EEG frequency bands to analyze %45 instead of 50
bandNames = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'}; %"Adjusted" beta band to correct for noise in Beta for 15 Hz scanner noise

%Babiloni 2020 guidelines: Delta: 0.1-<4, Theta: 4-<8, Alpha: 8-13, Beta: 14-30, Gamma: 30-<65,
%High Gamma 65-<90
% But FOOOF is traditionally capped at 50 Hz, so include a Gamma 30-50

subjects = dir(eegpath);
%subjects = dir([eegpath '*.mat']);
subjects = {subjects.name}';
subjects = subjects(endsWith(subjects, '.set'));
num_subjects = length(subjects);

fileID = fopen([fooof_dest_path, 'metadata.txt'], 'w+');
format shortg
c = clock;
fprintf(fileID, ['\r', 'Clock: ']);
fprintf(fileID, ' %d', c);

fprintf(fileID, ['\r\n', 'Number of Subjects: ', num2str(num_subjects)]);

fprintf(fileID, ['\r\n', 'Bands: ']);
band_str = cellstr(bandNames.');
fprintf(fileID, ' %s', band_str{:});

fprintf(fileID, ['\r\n', 'Bands Frequencies: ']);
band_str = bands';
fprintf(fileID, ' [%d %d]', band_str{:});

%fprintf(fileID, ['\r\n', 'Peak Width Limits: ']);
%peak_width_limits = settings.peak_width_limits;
%fprintf(fileID, ' [%d %d]', peak_width_limits);

%fprintf(fileID, ['\r\n', 'Max n Peaks: ']);
%fprintf(fileID, ' %d', settings.max_n_peaks);

%fprintf(fileID, ['\r\n', 'Min Peak Amplitude: ']);
%fprintf(fileID, ' %0.1f', settings.min_peak_amplitude);

fprintf(fileID, ['\r\n', 'Peak Threshold: ']);
fprintf(fileID, ' %0.1f', settings.peak_threshold);

fprintf(fileID, ['\r\n', 'Background Mode: ']);
fprintf(fileID, ' %s', settings.background_mode);

cd(eegpath)
% pat = '(\d+)-(\w).*';
% l=1;
% m=1;
% n=1;
% for q = 1:length(subjects)
%     t = regexp(subjects{q},pat,'tokens');
%     load(subjects{q});
%
%     % what is the all_eeg (now raw) function calling on?
%
%     if strcmp(t{1}{2},"CC")
%         all_eeg.hc{l} = subject.base_freq2sec;
%         all_eeg.hc{l}.subject_name = t{1}{1};
%         l = l + 1;
%     elseif strcmp(t{1}{2},"HC")
%         all_eeg.scz{m} = subject.base_freq2sec;
%         all_eeg.scz{m}.subject_name = t{1}{1};
%         m = m + 1;
%     elseif strcmp(t{1}{2},"UHR")
%         all_eeg.uhr{n} = subject.base_freq2sec;
%         all_eeg.uhr{n}.subject_name = t{1}{1};
%         n = n + 1;
%     end
%     clear subject
% end

% Pulling conditions from excel

% Read the Excel file
% check that you are in the SNAP_Groups folder
filename = 'M:\\SNAP Groups\\SNAP_Groups.xlsx';
%range = 'A2:B128'; % Update with the actual range of your data
%dataTable = readtable(filename, 'Range', range);% Update with your actual filename
dataTable = readtable(filename);

% Extract participant IDs from column 1
participant_ids = dataTable{:, 1};

% Extract conditions from column named 'Groups'
conditions = dataTable{:, 2};

% Display the extracted data (optional)
disp('Participant IDs:');
disp(participant_ids);

disp('Conditions:');
disp(conditions);

%New conditions logic
% Assuming dataTable is a table containing 'ParticipantID' and 'Groups' columns
for row = 1:height(dataTable)
    participant_id = dataTable.ParticipantID(row);
    current_condition = dataTable.Group{row};

    % Convert to cell arrays if needed
    conditions = cellstr(conditions);

    % Initialize an empty map
    participant_condition_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
end
% Populate the map
for i = 1:length(participant_ids)
    participant_ids_array = participant_ids(i);
    condition = conditions{i};
    participant_condition_map(participant_ids_array) = condition;
end

% Display the resulting map (optional, for verification)
disp('Participant Condition Map:');
disp(participant_condition_map);

% Initialize the structure to store EEG data
all_eeg = struct;

m = 1;
n=1;
o=1;
p=1;
% Load power spectral density data and channel label information for each participant
for i = 1:length(participant_ids)
    participant_id = participant_ids(i);
    condition = participant_condition_map(participant_id);

    % Assuming the file name is constructed based on the participant ID
    file_name = sprintf('%d_*.mat', participant_id);
    matching_files = dir(file_name);

    % Display some debugging information
    fprintf('Processing participant ID %d with condition %s\n', participant_id, condition);


    % Load the EEG data from the file and store it in the all_eeg structure
    %loaded_data = load(file_name);

    % Check if any files were found
    if ~isempty(matching_files)
        % Assuming there is only one matching file, use the first one
        %file_name = matching_files(1).name;
        % Display information about matching files
        for k = 1:length(matching_files)
            fprintf('Matching file found: %s\n', matching_files(k).name);
        end
    end
    % Loop through all matching files
    for file_index = 1:length(matching_files)
        file_name = matching_files(file_index).name;

        % Display debugging information about the file being processed
        fprintf('Processing file: %s\n', file_name);

        % Load the EEG data from the file and store it in the all_eeg structure
        load(file_name);
        if strcmp(condition,'CC')
            index = m;
            m=m+1;
        elseif strcmp(condition,'HC')
            index = n;
            n=n+1;
        elseif strcmp(condition,'UHR')
            index = o;
            o=o+1;
        elseif strcmp(condition,'UNK')
            index = p;
            p=p+1;
        end
        % Continue with storing the data in the all_eeg structure
        all_eeg.(condition){index}.label = subject.base_freq2sec.label;%struct('powspctrm', EEG.data, 'label', {EEG.chanlocs.labels});
        all_eeg.(condition){index}.powspctrm = subject.base_freq2sec.powspctrm;
        all_eeg.(condition){index}.subject_name = participant_id;
    end
    % else
    % fprintf('No matching file found for participant ID %d\n', participant_id);
end
% Assuming the loaded data contains 'powspctrm' and 'label' fields
% all_eeg.(condition){i} = struct('powspctrm', loaded_data.powspctrm, 'label', loaded_data.label);
%end

% % Check if the participant ID is in the mapping
% if isKey(participant_condition_map, participant_id)
%     % Assign the condition from the mapping
%     condition_assigned = participant_condition_map(participant_id);
% 
%     % Now you can use 'condition_assigned' in your code
%     subject_name = sprintf('Participant %d', participant_id);
% 
%     % Display for verification (you can remove this line in the final code)
%     disp(['Assigned condition ', condition_assigned, ' to subject ', subject_name]);
% 
%     % ... (rest of your code)
% else
%     % Handle the case when the participant ID is not found in the mapping
%     disp(['Participant ID ', num2str(participant_id), ' not found in the mapping']);
% 
% end

tic;  


for cond=1:size(conditions_overall,2)%3 conditions CC, HC, UHR, and UNK
    for s_i=1:length(all_eeg.(conditions_overall{1,cond}))
        subject_name = all_eeg.(conditions_overall{1,cond}){1,s_i}.subject_name;
        %try
        tic
        %         if cond==1
        %             condname='_EC.set';
        %             condition = 'EC';
        %         elseif cond==2
        %             condname='_EO.set';
        %             condition = 'EO';
        %         end
        %
        %         subject_file_name = subjects{s_i};
        %         subject_name = subject_file_name(1:10);
        %         full_name_with_codname = [subject_name condname];
        %
        %         EEG = pop_loadset(full_name_with_codname, eegpath);
        %
        %         if EEG.srate ~= 250
        %             fileID = fopen('ommited_lemon_subjects.txt','a');
        %             fprintf(fileID, ['\r\n', 'Subject ', num2str(s_i), ' ', subject_name, '\r\n'], []);
        %             fclose(fileID);
        %             continue;
        %         end

        %chanlocs = EEG.chanlocs;
        %save(['Z:\MATLAB\LEMONEEG_Preprocessed\Mat-files\chanlocs\' subject_name '.mat'],'chanlocs')

        %%{
        % Loop through each EEG channel, subjects and TR interval
        All=struct; %Aind=1;
        EEG_data = all_eeg.(conditions_overall{1,cond}){1,s_i}.powspctrm;
        EEG_label = all_eeg.(conditions_overall{1,cond}){1,s_i}.label;
        %Setup EEG data to get FOOOF'd

        %         %Create 2 second epochs, simimlar to TR interval
        %         EEG2=EEG;
        %         EEG2.event=[];
        %         for i = 1:180
        %             EEG2.event(1,i).type = '2sec';
        %             EEG2.event(1,i).latency = (i*EEG.srate*2)-(EEG.srate*1)+1;
        %             EEG2.event(1,i).urevent = 2;
        %         end
        %         EEG2 = pop_epoch( EEG2, {  '2sec'  }, [0  2]);
        %         EEG2_data = EEG2.data;
        %         size_EEG2_data = size(EEG2_data);


        %
        %         chans = 1:size_EEG2_data(1);
        num_chans = size(EEG_data,2);

        clear master_fooof_results;

        for e=1:num_chans
            try
                %             chan = chans(e);
                %             ChanData = squeeze(EEG2_data(chan, :, :)); %Pz=31
                %             % this?
                %             ChanDat=reshape(ChanData,1,[]);
                %             %
                %             epochSize = size(EEG2_data,2); %size(icaEEG.data,2);
                %             bins = 0.05:0.5:EEG.srate/2; %Frequency bins for 250 Hz sampling rate
                ChanDat = squeeze(EEG_data(:,e,:));

                %%Using params from Donoghue et al. 2020

                %Power spectra were calculated for all channels, using
                %Welch’s method51 (2-s windows, 50% overlap), for a 2-min segment of extracted
                % fit  using the algorithm, using the settings {peak_width_limits =  [1,6],  max_n_peaks = 6,  min_peak_height =  0.05,  peak_threshold = 1.5,  aperiodic_mode =  ‘fixed’}.

                Aind=1; %initialize this index
                %Loop through each 2 sec TR
                %             winstart=1;
                %             winsize=500;
                %             winjump=100;
                %             bins = 0.5:0.5:125;
                bins = 1:0.5:100;
                %             freqers=1:.5:45;%fooof_results(1,e,1).freqs;
                freqers=1:0.5:45;

                %Slide through each 2 sec TR change 503 to 1
                for A = 1:size(ChanDat,1)%(length(ChanDat))%2:size(ChanData,2)%size(icaEEG.data,3) %noVols; %right now, A is 105 not 182 PROBLEM!! ICA size!!
                    %   datawin(A-2
                    %get the power spectra from FFT
                    %
                    %                 ChanData=ChanDat(winstart:winstart+winsize-1)';
                    if ~any(isnan(ChanDat(A))) %several TR intervals will have NaNs due to fMRI
                        %                     X = abs(fft(ChanData.*hanning(winsize))).^2;
                        Xdat=ChanDat(A,1:199)'; %exclude high frequency data

                        %             if s==11 & e==14 & Aind==133 %Loop to deal with the occasional problem chan/subj/trial
                        %                 %fooof_results(s,e,Aind).background_params=([NaN NaN]);
                        %                 % THERE IS NO BACKGROUND PARAMS OUTPUT FROM FOOOF! New
                        %                 % version calls it 'aperiodic_params'!
                        %                 fooof_results(s,e,Aind).aperiodic_params=([NaN NaN]);
                        %                 for b=1:length(bandNames)
                        %                     fooofers.(bandNames{b})(s,e,Aind)=NaN;
                        %                 end
                        %             else

                        %Core function
                        fooof_results = fooof(bins,Xdat,freqrange,settings,1);
                        master_fooof_results(e, A) = fooof_results;

                        %Data organization for each EEG parameter
                        for b=1:length(bandNames)
                            subtract_spec=fooof_results.power_spectrum-fooof_results.ap_fit; %bg_fit; % there is no bg_fit
                            fooofers.(bandNames{b})(e,A)=nanmean(subtract_spec(find(freqers>bands{b}(1) & freqers<bands{b}(2))));
                            if isnan(nanmean(subtract_spec(find(freqers>bands{b}(1) & freqers<bands{b}(2)))));
                                stopthis=1;
                            end
                        end
                        fooofers.PSDexp(e,A)=fooof_results.aperiodic_params(2);
                        fooofers.PSDoff(e,A)=fooof_results.aperiodic_params(1);
                        fooofers.subjects{s_i}=subject_name;

                    else
                        fooof_results.aperiodic_params=([NaN NaN]);
                        fooofers.PSDexp(e,A)=fooof_results.aperiodic_params(2);
                        fooofers.PSDoff(e,A)=fooof_results.aperiodic_params(1);
                        for b=1:length(bandNames)
                            fooofers.(bandNames{b})(e,A)=NaN;
                        end
                    end
                    %                 Aind=Aind+1;
                    %                 winstart=winstart+winjump;
                    %                 condcomplete=cond;
                end
                disp(['Condition ', conditions_overall{cond},', Subject ', num2str(s_i), ' ', num2str(subject_name), ', ' 'Chan ', EEG_label{e}, ...
                    ' number ',num2str(e),' out of ', num2str(num_chans)])
            catch
                warning(['Condition ', conditions_overall{cond},', Subject ', num2str(s_i), ' ', num2str(subject_name), ', ' 'Chan ', EEG_label{e}, ...
                    ' number ',num2str(e),' out of ', num2str(num_chans)])
            end
        end
        toc
        save([fooof_dest_path filesep num2str(subject_name) '_' conditions_overall{cond} '.mat'],'fooofers')
        save([fooof_dest_path filesep num2str(subject_name) '_' conditions_overall{cond} '_results.mat'], 'master_fooof_results')
        fooofers=struct;
        %}
        chanlocs = EEG_label;
        save([chanlocs_dest_path subject_name '_' conditions_overall{cond} '.mat'],'chanlocs')
    end
    %{
    catch
        fileID = fopen('ommited_lemon_subjects.txt','a');
        fprintf(fileID, ['\r\n', 'Subject ', num2str(s_i), ' ', subject_name, '\r\n'], []);
        fclose(fileID);
        continue
    end
    %}
end
% %%try
% %tic
% %         if cond==1
% %             condname='_EC.set';
% %             condition = 'EC';
% %         elseif cond==2
% %             condname='_EO.set';
% %             condition = 'EO';
% %         end
% %
% %         subject_file_name = subjects{s_i};
% %         subject_name = subject_file_name(1:10);
% %         full_name_with_codname = [subject_name condname];
% %
% %         EEG = pop_loadset(full_name_with_codname, eegpath);
% %
% %         if EEG.srate ~= 250
% %             fileID = fopen('ommited_lemon_subjects.txt','a');
% %             fprintf(fileID, ['\r\n', 'Subject ', num2str(s_i), ' ', subject_name, '\r\n'], []);
% %             fclose(fileID);
% %             continue;
% %         end
% 
% %chanlocs = EEG.chanlocs;
% %save(['Z:\MATLAB\LEMONEEG_Preprocessed\Mat-files\chanlocs\' subject_name '.mat'],'chanlocs')
% 
% %%{
% % Loop through each EEG channel, subjects and TR interval
% All=struct; %Aind=1;
% %EEG_data = all_eeg.(conditions{1,cond}){1,s_i}.powspctrm;
% %EEG_label = all_eeg.(conditions{1,cond}){1,s_i}.label;
% EEG_data = all_eeg.(conditions{1}).(participant_id).powspctrm;
% EEG_label = all_eeg.(conditions{1}).(participant_id).label;
% 
% %Setup EEG data to get FOOOF'd
% 
% % Display the resulting structure (optional, for verification)
% disp('Updated all_eeg structure:');
% disp(all_eeg);
% 
% %         %Create 2 second epochs, simimlar to TR interval
% %         EEG2=EEG;
% %         EEG2.event=[];
% %         for i = 1:180
% %             EEG2.event(1,i).type = '2sec';
% %             EEG2.event(1,i).latency = (i*EEG.srate*2)-(EEG.srate*1)+1;
% %             EEG2.event(1,i).urevent = 2;
% %         end
% %         EEG2 = pop_epoch( EEG2, {  '2sec'  }, [0  2]);
% %         EEG2_data = EEG2.data;
% %         size_EEG2_data = size(EEG2_data);
% 
% 
% %
% %         chans = 1:size_EEG2_data(1);
% num_chans = size(EEG_data,2);
% 
% clear master_fooof_results;
% 
% for e=1:num_chans
%     %             chan = chans(e);
%     %             ChanData = squeeze(EEG2_data(chan, :, :)); %Pz=31
%     %             % this?
%     %             ChanDat=reshape(ChanData,1,[]);
%     %             %
%     %             epochSize = size(EEG2_data,2); %size(icaEEG.data,2);
%     %             bins = 0.05:0.5:EEG.srate/2; %Frequency bins for 250 Hz sampling rate
%     ChanDat = squeeze(EEG_data(:,e,:));
% 
%     %%Using params from Donoghue et al. 2020
% 
%     %Power spectra were calculated for all channels, using
%     %Welch’s method51 (2-s windows, 50% overlap), for a 2-min segment of extracted
%     % fit  using the algorithm, using the settings {peak_width_limits =  [1,6],  max_n_peaks = 6,  min_peak_height =  0.05,  peak_threshold = 1.5,  aperiodic_mode =  ‘fixed’}.
% 
%     Aind=1; %initialize this index
%     %Loop through each 2 sec TR
%     Aind=1;
%     %             winstart=1;
%     %             winsize=500;
%     %             winjump=100;
%     %             bins = 0.5:0.5:125;
%     bins = 1:0.5:100;
%     %             freqers=1:.5:45;%fooof_results(1,e,1).freqs;
%     freqers=1:0.5:45;
% 
%     %Slide through each 2 sec TR change 503 to 1
%     for A = 1:size(ChanDat,1)%(length(ChanDat))%2:size(ChanData,2)%size(icaEEG.data,3) %noVols; %right now, A is 105 not 182 PROBLEM!! ICA size!!
%         %   datawin(A-2
%         %get the power spectra from FFT
%         %
%         %                 ChanData=ChanDat(winstart:winstart+winsize-1)';
%         if ~any(isnan(ChanDat(A))) %several TR intervals will have NaNs due to fMRI
%             %                     X = abs(fft(ChanData.*hanning(winsize))).^2;
%             Xdat=ChanDat(A,1:199)'; %exclude high frequency data
% 
%             %             if s==11 & e==14 & Aind==133 %Loop to deal with the occasional problem chan/subj/trial
%             %                 %fooof_results(s,e,Aind).background_params=([NaN NaN]);
%             %                 % THERE IS NO BACKGROUND PARAMS OUTPUT FROM FOOOF! New
%             %                 % version calls it 'aperiodic_params'!
%             %                 fooof_results(s,e,Aind).aperiodic_params=([NaN NaN]);
%             %                 for b=1:length(bandNames)
%             %                     fooofers.(bandNames{b})(s,e,Aind)=NaN;
%             %                 end
%             %             else
% 
%             %Core function
%             fooof_results = fooof(bins,Xdat,freqrange,settings,1);
%             master_fooof_results(e, A) = fooof_results;
% 
%             %Data organization for each EEG parameter
%             for b=1:length(bandNames)
%                 subtract_spec=fooof_results.power_spectrum-fooof_results.ap_fit; %bg_fit; % there is no bg_fit
%                 fooofers.(bandNames{b})(e,A)=nanmean(subtract_spec(find(freqers>bands{b}(1) & freqers<bands{b}(2))));
%                 if isnan(nanmean(subtract_spec(find(freqers>bands{b}(1) & freqers<bands{b}(2)))));
%                     stopthis=1;
%                 end
%             end
%             fooofers.PSDexp(e,A)=fooof_results.aperiodic_params(2);
%             fooofers.PSDoff(e,A)=fooof_results.aperiodic_params(1);
%             %  fooofers.subjects{s_i}=subject_name;
%             fooofers.subjects{participant_id}=subject_name;
% 
%         else
%             fooof_results.aperiodic_params=([NaN NaN]);
%             fooofers.PSDexp(e,A)=fooof_results.aperiodic_params(2);
%             fooofers.PSDoff(e,A)=fooof_results.aperiodic_params(1);
%             for b=1:length(bandNames)
%                 fooofers.(bandNames{b})(e,A)=NaN;
%             end
%         end
%         %                 Aind=Aind+1;
%         %                 winstart=winstart+winjump;
%         %                 condcomplete=cond;
%     end
%     try
%         %  disp(['Condition ', conditions{cond},', Subject ', num2str(s_i), ' ', subject_name, ', ' 'Chan ', EEG_label{e}, ...
%         disp(['Condition ', conditions{cond},', Subject ', num2str(participant_id), ' ', subject_name, ', ' 'Chan ', EEG_label{e}, ...
%             ' number ',num2str(e),' out of ', num2str(num_chans)])
%     catch
%         %warning(['Condition ', conditions{cond},', Subject ', num2str(s_i), ' ', subject_name, ', ' 'Chan ', EEG_label{e}, ...
%         warning(['Condition ', conditions{cond},', Subject ', num2str(participant_id), ' ', subject_name, ', ' 'Chan ', EEG_label{e}, ...
%             ' number ',num2str(e),' out of ', num2str(num_chans)])
%     end
% end
% toc
% save([fooof_dest_path subject_name '_' conditions{cond} '.mat'],'fooofers')
% save([fooof_dest_path subject_name '_' conditions{cond} '_results.mat'], 'master_fooof_results')
% fooofers=struct;
% %}
% chanlocs = EEG_label;
% save([chanlocs_dest_path subject_name '_' conditions{cond} '.mat'],'chanlocs')
% 
% %{
%     catch
%         fileID = fopen('ommited_lemon_subjects.txt','a');
%         fprintf(fileID, ['\r\n', 'Subject ', num2str(s_i), ' ', subject_name, '\r\n'], []);
%         fclose(fileID);
%         continue
%     end
% %}
% 
% 
