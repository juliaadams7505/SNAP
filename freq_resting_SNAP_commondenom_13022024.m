clearvars; close all;  clc;
% load averages for each individual subject, for each condition
% load freq_resting_stats.mat
%addpath '/Users/elrowe/Desktop/Julia_Data/fieldtrip'; ft_defaults; % <---- UPDATE
addpath 'C:\\Users\\juadams\\Documents\\Field Trip\\fieldtrip-20230118'; ft_defaults; % <---- UPDATE
% Defining regions of interest

FL_ROI = {'F7','F5','F3'}; 
FM_ROI = {'FP1', 'AF3', 'F1', 'FPZ', 'FP2', 'AF4', 'F2', 'FZ'};
FR_ROI = {'F4','F6','F8'};
PL_ROI = {'FT7','FC5','FC3', 'T7', 'C5', 'C3', 'TP7', 'CP5', 'CP3'};
PM_ROI = {'FC1', 'FCZ', 'FC2', 'C1', 'CZ', 'C2','CP1', 'CPZ', 'CP2'};
PR_ROI = {'FC4','FC6','FT8','C4','C6','T8','CP4', 'CP6', 'TP8'};
OL_ROI = {'P7','P5','P3','PO7','PO5','O1', 'CB1'};
OM_ROI = {'P1', 'PZ', 'P2', 'PO3','POZ', 'PO4', 'OZ'};
OR_ROI = {'P4','P6','P8','PO6','PO8','02','CB2'};
%ONE_FRONT_ROI = {'C5'};
%ONE_BACK_ROI = {'FT8'};

% adjust this section according to your paths
%data_path = '/Users/elrowe/Desktop/Julia_Data/processed'
data_path = 'M:\\SNAP_Data_Processed\\processed';
cd(data_path);
%list = dir('*_pad_bandwidth1.mat'); %change this depending on which bandwidth you are examining  <---- UPDATE
list = dir('*_cleanLine_rmline_pad.mat'); %change this depending on which bandwidth you are examining  <---- UPDATE


for pp = 1:length(list)

    load(list(pp).name); %load this participant

    %Find PNo (changes how we look for the group identifier)
    pNo = sscanf(list(pp).name,'%d_');

    %Determine what the group identifier is for this participant
    if pNo < 1000
         pGroupName = list(pp).name(5:7);
         group_identifier = strtok(pGroupName,'.');
    elseif pNo >=1000
        pGroupName = list(pp).name(6:8);
        group_identifier = strtok(pGroupName,'.');
    end

    %Assign subject data to structure depending on their group
    if strcmp(group_identifier,'CC')
        group_data_CC{ctrA,1} = subject.base_freq2sec;
        group_data_CC_all{ctrA,1} = subject.base_freq2sec;
        ctrA = ctrA + 1;
    elseif strcmp(group_identifier,'HC')
        group_data_HC{ctrB,1} = subject.base_freq2sec;
        group_data_HC_all{ctrB,1} = subject.base_freq2sec;
        ctrB = ctrB + 1;
    elseif strcmp(group_identifier,'UHR')
        group_data_UHR{ctrC,1} = subject.base_freq2sec;
        group_data_UHR_all{ctrC,1} = subject.base_freq2sec;
        ctrC = ctrC + 1;
    else
        group_data_none{ctrD,1} = subject.base_freq2sec; %record people with NO group
        group_data_none_all{ctrD,1} = subject.base_freq2sec;
        ctrD = ctrD + 1;
    end

    % load the electrode positions

if pp == 1
   elec = subject.elec;
end


    clear subject filename new_data clean_signal signal
end

% %% Now we compute the group grand mean average
% 
% % Use as
% %   [grandavg] = ft_freqgrandaverage(cfg, freq1, freq2, freq3...)
% %
% % The input data freq1..N are obtained from either FT_FREQANALYSIS with
% % keeptrials=no or from FT_FREQDESCRIPTIVES. The configuration structure
% % can contain
% %   cfg.keepindividual = 'yes' or 'no' (default = 'no')
% %   cfg.foilim         = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
% %   cfg.toilim         = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
% %   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
% %                        see FT_CHANNELSELECTION for details
% %   cfg.parameter      = string or cell-array of strings indicating which
% %                        parameter(s) to average. default is set to
% %                        'powspctrm', if it is present in the data.
% 
% cfg.keepindividual = 'no'; %do not keep individual, average over participants within the group
% cfg.foilim         = [30 100]; %frequency window of interest (change to gamma, 30 to 100)
% 
% [grandavg_HC] = ft_freqgrandaverage(cfg,group_data_HC{1:end}) %add up to number of participants
% [grandavg_CC] = ft_freqgrandaverage(cfg,group_data_CC{1:end}) %add up to number of participants
% [grandavg_UHR] = ft_freqgrandaverage(cfg,group_data_UHR{1:end}) %add up to number of participants
% 
% % cfg = [];
% % cfg.keepindividual = 'yes';
% HC = ft_freqgrandaverage(cfg, HC{:});
% %subject.base_freq2sec
% %subject.elec
% CC = ft_freqgrandaverage(cfg, CC{:});
% UHR = ft_freqgrandaverage(cfg, UHR{:});

cfg = [];
cfg.keepindividual = 'yes';
group_data_HC = ft_freqgrandaverage(cfg, group_data_HC{:});
group_data_CC = ft_freqgrandaverage(cfg, group_data_CC{:});
group_data_UHR = ft_freqgrandaverage(cfg, group_data_UHR{:});

% Get numerical indices of ROI to compute averages
sel_FL_ROI = match_str(elec.label, FL_ROI);
sel_FM_ROI = match_str(elec.label, FM_ROI);
sel_FR_ROI = match_str(elec.label, FR_ROI);
sel_PL_ROI = match_str(elec.label, PL_ROI);
sel_PM_ROI = match_str(elec.label, PM_ROI);
sel_PR_ROI = match_str(elec.label, PR_ROI);
sel_OL_ROI = match_str(elec.label, OL_ROI);
sel_OM_ROI = match_str(elec.label, OM_ROI);
sel_OR_ROI = match_str(elec.label, OR_ROI);
%sel_ONE_FRONT_ROI = match_str(elec.label, ONE_FRONT_ROI);
%sel_ONE_BACK_ROI = match_str(elec.label, ONE_BACK_ROI);

% Normalize power spectra using mean power over freq range
freq_oi   = [30 100];   % frequency range of interest
freq_norm = [1 100]; % frequency range used to normalize the spectrum
foi_norm  = nearest(group_data_HC.freq, freq_norm);

common_denominator_HC = mean(group_data_HC.powspctrm(:,:,foi_norm(1):foi_norm(2)),3);
common_denominator_CC = mean(group_data_CC.powspctrm(:,:,foi_norm(1):foi_norm(2)),3);
common_denominator_UHR = mean(group_data_UHR.powspctrm(:,:,foi_norm(1):foi_norm(2)),3);

group_data_HC.powspctrm_b = group_data_HC.powspctrm./common_denominator_HC;
group_data_CC.powspctrm_b = group_data_CC.powspctrm./common_denominator_CC;
group_data_UHR.powspctrm_b = group_data_UHR.powspctrm./common_denominator_UHR;

% plot Neighbours

cfg = [];
cfg.method = 'distance'
cfg.elec = elec
neighbours = ft_prepare_neighbours(cfg)

% Collapse data using ROI and frequency ranges defined in relevant paper(s)
cfg = [];
cfg.channel     = FL_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_FL_ROI = ft_selectdata(cfg, group_data_HC_grid_FL);
group_data_CC_FL_ROI = ft_selectdata(cfg, group_data_CC_grid_FL);
group_data_UHR_FL_ROI = ft_selectdata(cfg, group_data_UHR_grid_FL);

cfg = [];
cfg.channel     = FM_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_FM_ROI = ft_selectdata(cfg, group_data_HC_grid_FM);
group_data_CC_FM_ROI = ft_selectdata(cfg, group_data_CC_grid_FM);
group_data_UHR_FM_ROI = ft_selectdata(cfg, group_data_UHR_grid_FM);

cfg = [];
cfg.channel     = FR_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_FR_ROI = ft_selectdata(cfg, group_data_HC_grid_FR);
group_data_CC_FR_ROI = ft_selectdata(cfg, group_data_CC_grid_FR);
group_data_UHR_FR_ROI = ft_selectdata(cfg, group_data_UHR_grid_FR);

cfg = [];
cfg.channel     = PL_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_PL_ROI = ft_selectdata(cfg, group_data_HC_grid_PL);
group_data_CC_PL_ROI = ft_selectdata(cfg, group_data_CC_grid_PL);
group_data_UHR_PL_ROI = ft_selectdata(cfg, group_data_UHR_grid_PL);

cfg = [];
cfg.channel     = PM_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_PM_ROI = ft_selectdata(cfg, group_data_HC_grid_PM);
group_data_CC_PM_ROI = ft_selectdata(cfg, group_data_CC_grid_PM);
group_data_UHR_PM_ROI = ft_selectdata(cfg, group_data_UHR_grid_PM);

cfg = [];
cfg.channel     = PR_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_PR_ROI = ft_selectdata(cfg, group_data_HC_grid_PR);
group_data_CC_PR_ROI = ft_selectdata(cfg, group_data_CC_grid_PR);
group_data_UHR_PR_ROI = ft_selectdata(cfg, group_data_UHR_grid_PR);

cfg = [];
cfg.channel     = OL_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_OL_ROI = ft_selectdata(cfg, group_data_HC_grid_OL);
group_data_CC_OL_ROI = ft_selectdata(cfg, group_data_CC_grid_OL);
group_data_UHR_OL_ROI = ft_selectdata(cfg, group_data_UHR_grid_OL);


cfg = [];
cfg.channel     = OM_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_OM_ROI = ft_selectdata(cfg, group_data_HC_grid_OM);
group_data_CC_OM_ROI = ft_selectdata(cfg, group_data_CC_grid_OM);
group_data_UHR_OM_ROI = ft_selectdata(cfg, group_data_UHR_grid_OM);

cfg = [];
cfg.channel     = OM_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
group_data_HC_OR_ROI = ft_selectdata(cfg, group_data_HC_grid_OR);
group_data_CC_OR_ROI = ft_selectdata(cfg, group_data_CC_grid_OR);
group_data_UHR_OR_ROI = ft_selectdata(cfg, group_data_UHR_grid_OR);

% cfg.channel     = occipital_ROI;
% group_data_HC_oROI = ft_selectdata(cfg, group_data_HC);
% group_data_CC_oROI = ft_selectdata(cfg, group_data_CC);
% group_data_UHR_oROI = ft_selectdata(cfg, group_data_UHR);

% Plotspread function - may need to download from download server linked in
% fieldtrip tutorial page 
 
data_raw_FL_ROI    = {group_data_HC_FL_ROI.powspctrm...
  group_data_CC_FL_ROI.powspctrm...
  group_data_UHR_FL_ROI.powspctrm};

data_raw_FM_ROI    = {group_data_HC_FM_ROI.powspctrm...
  group_data_CC_FM_ROI.powspctrm...
  group_data_UHR_FM_ROI.powspctrm};
  
data_between_fROI ={group_data_HC_fROI.powspctrm_b...
  group_data_CC_fROI.powspctrm_b...
  group_data_UHR_fROI.powspctrm_b};
  
data_between_oROI ={group_data_HC_oROI.powspctrm_b...
  group_data_CC_oROI.powspctrm_b...
  group_data_UHR_oROI.powspctrm_b};

% Make dot plots for ROIs 1-4 in different conditions as (x = condition, y =
% absolute power)

figure

subplot(2,3,1);
plotSpread(data_raw_FL_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk',});
ylabel('abs. power (V^2)');
title('raw PSD Front');

subplot(2,3,2);
h3 = plotSpread(data_between_FL_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Front');
set(h3{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,3);
plotSpread(data_raw_FM_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('abs. power (V^2)');
title('raw PSD Occip');

subplot(2,3,4);
h6 = plotSpread(data_between_FM_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,5);
h6 = plotSpread(data_between_FR_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,6);
h6 = plotSpread(data_between_FR_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

figure2

subplot(2,3,1);
plotSpread(data_raw_PL_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk',});
ylabel('abs. power (V^2)');
title('raw PSD Front');

subplot(2,3,2);
h3 = plotSpread(data_between_PL_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Front');
set(h3{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,3);
plotSpread(data_raw_PM_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('abs. power (V^2)');
title('raw PSD Occip');

subplot(2,3,4);
h6 = plotSpread(data_between_PM_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,5);
h6 = plotSpread(data_between_PR_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,6);
h6 = plotSpread(data_between_PR_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

figure3

subplot(2,3,1);
plotSpread(data_raw_OL_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk',});
ylabel('abs. power (V^2)');
title('raw PSD Front');

subplot(2,3,2);
h3 = plotSpread(data_between_OL_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Front');
set(h3{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,3);
plotSpread(data_raw_OM_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('abs. power (V^2)');
title('raw PSD Occip');

subplot(2,3,4);
h6 = plotSpread(data_between_OM_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,5);
h6 = plotSpread(data_between_OR_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,3,6);
h6 = plotSpread(data_between_OR_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

% Topoplots and power spectra for reach ROI in each sedative condition
cfg = [];
cfg.elec             = elec;
cfg.parameter        = 'powspctrm_b'; % you can plot either powspctrm (default) or powspctrm_b
cfg.xlim             = [55 100]; % frequency range to make the topoplot
cfg.highlight        = 'on';

% Options for improving appearance of plots (next 60 lines)

cfg.highlightchannel = {frontal_ROI occipital_ROI};
cfg.highlightsymbol  = {'o','*'};
cfg.highlightcolor   = [0 0 0];
cfg.highlightsize    = 6;
cfg.markersymbol     = '.';
cfg.comment          = 'no';
cfg.colormap         = 'jet';

% NOTE - maybe need to fix plot numbering since one less group

figure('position',[680 240 1039 420]);
subplot(2,4,1); ft_topoplotER(cfg, group_data_HC); colorbar; title('Healthy Control');
subplot(2,4,2); ft_topoplotER(cfg, group_data_CC); colorbar; title('Clinical Control');
subplot(2,4,3); ft_topoplotER(cfg, group_data_UHR); colorbar; title('Ultra-High Risk');

subplot(2,4,5);loglog(group_data_HC.freq,...
  [squeeze(mean(mean(group_data_HC.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(group_data_HC.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
legend('Front ROI','Occip ROI','Location','southwest');
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('Healthy Control');

subplot(2,4,6);loglog(group_data_CC.freq,...
  [squeeze(mean(mean(group_data_CC.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(group_data_CC.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('Clinical Control');

subplot(2,4,7);loglog(group_data_UHR.freq,...
  [squeeze(mean(mean(group_data_UHR.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(group_data_UHR.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('Ultra High Risk');

%%%% BETWEEN-PARTICIPANT CONTRASTS %%%%%%

% Statistically evaluate power difference between groups 
% NOTE - check drowsy & response and how to change for our dataset

% [HC_respon HC_drowsy] = deal(HC);
% [CC_respon CC_drowsy] = deal(CC);
% [UHR_respon UHR_drowsy] = deal(UHR);

cfg = [];
cfg.parameter   = {'powspctrm','powspctrm_b'};

% cfg.trials      = respon_group; % cfg.trials will select the 'subj' dimension
% HC_respon = ft_selectdata(cfg, HC_respon);
% CC_respon = ft_selectdata(cfg, CC_respon);
% UHR_respon = ft_selectdata(cfg, UHR_respon);

% cfg.trials      = drowsy_group; % cfg.trials will select the 'subj' dimension
% HC_drowsy = ft_selectdata(cfg, HC_drowsy);
% CC_drowsy = ft_selectdata(cfg, CC_drowsy);
% UHR_drowsy = ft_selectdata(cfg, UHR_drowsy);

% Then, the permutation test
% Differences:
% 1) Different sample level effect statistic
% 2) Different design matrix
% 3) Independent variable = group assignment
% Using independent samples T-statistic (BETWEEN subjects)
% Design matrix includes group assignment 

cfg                 = [];
cfg.channel          = 'CZ';
cfg.avgovergchan     = 'no';
cfg.frequency        = [55 100];
cfg.avgovergfreq     = 'yes';
cfg.parameter        = 'powspctrm_b';
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesF';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 2;
cfg.tail             = 1;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.computeprob      = 'yes';
cfg.numrandomization = 1000;
cfg.neighbours       = neighbours;

design = ones(1,size(group_data_HC_all,1) + size(group_data_CC_all,1) + size(group_data_UHR_all,1))*3;
design(1,1:size(group_data_HC_all,1)) = 1;
design(1,(size(group_data_HC_all,1)+1):(size(group_data_HC_all,1)+size(group_data_CC_all,1))) = 2;

cfg.design = design;
cfg.ivar   = 1;

% Use configuration to perform analysis 

stat2 = ft_freqstatistics(cfg, group_data_HC, group_data_CC, group_data_UHR);
% stat2 = ft_freqstatistics(cfg, group_data_HC, group_data_CC);
% stat3 = ft_freqstatistics(cfg, group_data_CC, group_data_UHR);
% stat4 = ft_freqstatistics(cfg, group_data_UHR, group_data_HC);

% Extract the data manually from ONE (or multiple) ELECTRODES and run stats
% on the power spectrum within the defined frequency bands
foi = [55 100]; %frequency bands of interest
foi_idx = nearest(group_data_CC.freq(:),foi)
electrodeName = 'CZ'; %this electrode
elecIdx = find(strcmp(elec(:).label,electrodeName)==1)
%Loop through all participants and extract the data
% HEALTHY CONTROLS FIRST
pwr_data_group_HC = squeeze(group_data_HC.powspctrm_b(:,elecIdx,foi_idx(1):foi_idx(2)))
m_pwr_data_group_HC = mean(pwr_data_group_HC,2)
pwr_data_group_CC = squeeze(group_data_CC.powspctrm_b(:,elecIdx,foi_idx(1):foi_idx(2)))
m_pwr_data_group_CC = mean(pwr_data_group_CC,2)
pwr_data_group_UHR = squeeze(group_data_UHR.powspctrm_b(:,elecIdx,foi_idx(1):foi_idx(2)))
m_pwr_data_group_UHR = mean(pwr_data_group_UHR,2)

%Run ANOVA1
group = [ones(size(pwr_data_group_HC,1),1);...
    (ones(size(pwr_data_group_CC,1),1)*2);ones(size(pwr_data_group_UHR,1),1)*3]'
[p,anovatab,stats] = anova1([m_pwr_data_group_HC;m_pwr_data_group_CC;m_pwr_data_group_UHR],group)
%Now run stats on these data
[h,p,stats,CI] = ttest2(m_pwr_data_group_HC,m_pwr_data_group_UHR)

% Plot the results!!!!

cfg = [];
cfg.alpha     = stat2.cfg.alpha;
cfg.parameter = 'stat';
cfg.zlim      = [-3 3];
cfg.elec      = elec;
ft_clusterplot(cfg, stat2);

cfg = [];
cfg.elec         = elec;
cfg.zlim         = [0 0.3];
cfg.xlim         = [55 100];
cfg.parameter    = 'powspctrm_b';
cfg.markersymbol = '.';
cfg.comment      = 'no';
cfg.colormap     = 'jet';
cfg.colorbar     = 'no';

figure('position',[680 240 1039 420]);
figure; ft_topoplotER(cfg, group_data_HC); colorbar; title('Healthy Control [55 100]');
figure; ft_topoplotER(cfg, group_data_CC); colorbar; title('Clinical Control [55 100]');
figure; ft_topoplotER(cfg, group_data_UHR); colorbar; title('Ultra-High Risk [55 100]');


% Plot power function for each group as a function of sedative state 
% First, for frontal ROI 

figure;

subplot(2,4,1); loglog(HC_respon.freq,...
  [squeeze(mean(mean(HC_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(HC_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
legend('Front ROI resp','Front ROI drow','Location','southwest');
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('baseline');

subplot(2,4,2); loglog(CC_respon.freq,...
  [squeeze(mean(mean(CC_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(CC_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('mild');

subplot(2,4,3); loglog(UHR_respon.freq,...
  [squeeze(mean(mean(UHR_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(UHR_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('moderate');

subplot(2,4,4); loglog(reco_sedation_respon.freq,...
  [squeeze(mean(mean(reco_sedation_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(reco_sedation_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('recover');

% Then, the occipital ROI: 

subplot(2,4,5); loglog(HC_respon.freq,...
  [squeeze(mean(mean(HC_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
  squeeze(mean(mean(HC_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
legend('Occip ROI resp','Occip ROI drow','Location','southwest');

subplot(2,4,6); loglog(CC_respon.freq,...
  [squeeze(mean(mean(CC_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
  squeeze(mean(mean(CC_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);

subplot(2,4,7); loglog(UHR_respon.freq,...
  [squeeze(mean(mean(UHR_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
  squeeze(mean(mean(UHR_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);

subplot(2,4,8); loglog(reco_sedation_respon.freq,...
  [squeeze(mean(mean(reco_sedation_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
  squeeze(mean(mean(reco_sedation_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);

% Or, ANOVA to copmare all three groups 

cfg = [];
cfg.channel          = 'all';
cfg.frequency        = foi_contrast;
cfg.avgovergfreq     = 'no';
cfg.parameter        = 'powspctrm_b';
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesFmultivariate';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum'; %'maxsum', 'maxsize', 'wcm'
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 2;
cfg.tail             = 1; % For a F-statistic, it only make sense to calculate the right tail
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.computeprob      = 'yes';
cfg.numrandomization = 500;
cfg.neighbours       = cfg_neigh.neighbours;

nsubj = size(covariates,1);
design = zeros(2,4*nsubj);
design(1,1:nsubj)           = 1;
design(1,nsubj+1:2*nsubj)   = 2;
design(1,nsubj*2+1:3*nsubj) = 3;
design(1,nsubj*3+1:4*nsubj) = 4;
design(2,:) = repmat(1:nsubj,1,4);

cfg.design = design;
cfg.ivar   = 1; % sedation level
cfg.uvar   = 2; % subject number

% Pass the data from all four conditions as input variables 

stat3 = ft_freqstatistics(cfg, HC, CC, UHR, reco_sedation);

% Plot them results
cfg            = [];
cfg.frequency  = foi_contrast;
cfg.avgoverrpt = 'yes';
cfg.parameter  = {'powspctrm','powspctrm_b'};
HC_avg = ft_selectdata(cfg, HC);
CC_avg = ft_selectdata(cfg, CC);
UHR_avg = ft_selectdata(cfg, UHR);
reco_sedation_avg = ft_selectdata(cfg, reco_sedation);

% copy the mask field to each variable
HC_avg.mask = stat3.mask;
CC_avg.mask = stat3.mask;
UHR_avg.mask = stat3.mask;
reco_sedation_avg.mask = stat3.mask;

cfg = [];
cfg.zlim          = [0 90];
cfg.elec          = elec;
cfg.colorbar      = 'no';
cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle     = 'box';
cfg.parameter     = 'powspctrm_b';
cfg.maskfacealpha = 0.1;
figure; ft_multiplotER(cfg, HC_avg, CC_avg, UHR_avg, reco_sedation_avg);

