clearvars; close all;  clc;
% load averages for each individual subject, for each condition
% load freq_resting_stats.mat
%addpath '/Users/elrowe/Desktop/Julia_Data/fieldtrip'; ft_defaults; % <---- UPDATE
addpath 'C:\\Users\\juadams\\Documents\\Field Trip\\fieldtrip-20230118'; ft_defaults; % <---- UPDATE
% Defining regions of interest

F_ROI = {'F7','F5','F3', 'FP1', 'AF3', 'F1', 'FPZ', 'FP2', 'AF4', 'F2', 'FZ','F4','F6','F8'};
P_ROI = {'FT7','FC5','FC3', 'T7', 'C5', 'C3', 'TP7', 'CP5', 'CP3','FC1', 'FCZ', 'FC2', 'C1', 'CZ', 'C2','CP1', 'CPZ', 'CP2','FC4','FC6','FT8','C4','C6','T8','CP4', 'CP6', 'TP8'};
O_ROI = {'P7','P5','P3','PO7','PO5','O1', 'CB1','P1', 'PZ', 'P2', 'PO3','POZ', 'PO4', 'OZ','P4','P6','P8','PO6','PO8','02','CB2'};

%ONE_FRONT_ROI = {'C5'};
%ONE_BACK_ROI = {'FT8'};

% adjust this section according to your paths
%data_path = '/Users/elrowe/Desktop/Julia_Data/processed'
data_path = 'M:\\SNAP_Data_Processed\\processed\\FINAL_FOOOFED_RMLINE';
cd(data_path);
%list = dir('*_pad_bandwidth1.mat'); %change this depending on which bandwidth you are examining  <---- UPDATE
list = dir('*_results.mat'); %change this depending on which bandwidth you are examining  <---- UPDATE

ctrA = 1;
ctrB = 1;
ctrC = 1;
ctrD = 1; 

for pp = 122:length(list)

    load(list(pp).name); %load this participant

    %Find PNo (changes how we look for the group identifier)
    pNo = sscanf(list(pp).name,'%d_');

    %Determine what the group identifier is for this participant
    if pNo < 1000
         pGroupName = list(pp).name(5:7);
         group_identifier = strtok(pGroupName,'_*');
    elseif pNo >=1000
        pGroupName = list(pp).name(6:8);
        group_identifier = strtok(pGroupName,'_*');
    end

       % load the electrode positions
load('M:\SNAP_Data_Processed\processed\137_HC.mat');
   elec = subject.elec;
   template = subject.base_freq2sec
   template.freq = template.freq(1:89)

   for ch = 1:size(master_fooof_results,1);
       for tr = 1:size(master_fooof_results,2);
           tempA(ch,tr,:) = master_fooof_results(ch,tr).fooofed_spectrum;

       end
   end

      meantempA = squeeze(mean(tempA,2));
      template.powspctrm = meantempA;

    %Assign subject data to structure depending on their group
    if strcmp(group_identifier,'CC')
        group_data_CC{ctrA,1} = template;
        group_data_CC_all{ctrA,1} = template;
        ctrA = ctrA + 1;
    elseif strcmp(group_identifier,'HC')
        group_data_HC{ctrB,1} = template;
        group_data_HC_all{ctrB,1} = template;
        ctrB = ctrB + 1;
    elseif strcmp(group_identifier,'UHR')
        group_data_UHR{ctrC,1} = template;
        group_data_UHR_all{ctrC,1} = template;
        ctrC = ctrC + 1;
    else
        group_data_none{ctrD,1} = template; %record people with NO group
        group_data_none_all{ctrD,1} = template;
        ctrD = ctrD + 1;
    end


if pp == 1
  
end


    clear subject filename new_data clean_signal signal
end

  load('M:\SNAP_Data_Processed\processed\137_HC.mat');
   elec = subject.elec;

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
cfg.keepindividual = 'no'; %do not keep individual, average over participants within the group
cfg.foilim         = [30 100]; %frequency window of interest (change to gamma, 30 to 100)

[grandavg_HC] = ft_freqgrandaverage(cfg,group_data_HC{1:end}) %add up to number of participants
[grandavg_CC] = ft_freqgrandaverage(cfg,group_data_CC{1:end}) %add up to number of participants
[grandavg_UHR] = ft_freqgrandaverage(cfg,group_data_UHR{1:end}) %add up to number of participants

% % cfg = [];
% % cfg.keepindividual = 'yes';
% HC = ft_freqgrandaverage(cfg, HC{:});
% %subject.base_freq2sec
% %subject.elec
% CC = ft_freqgrandaverage(cfg, CC{:});
% UHR = ft_freqgrandaverage(cfg, UHR{:});

figure;
hold on;
plot(grandavg_HC.freq, squeeze(mean((grandavg_HC.powspctrm))), 'LineWidth',3)
plot(grandavg_CC.freq, squeeze(mean((grandavg_CC.powspctrm))), 'LineWidth',3)
plot(grandavg_UHR.freq, squeeze(mean((grandavg_UHR.powspctrm))), 'LineWidth',3)
title('Power Spectrum Model')
set(gca,'FontSize', 20); xlabel('Frequency'); ylabel('log(Power)')
legend('Healthy Control', 'Clinical Control', 'Ultra-High-Risk');  
% 
cfg = [];
cfg.keepindividual = 'yes';
group_data_HC = ft_freqgrandaverage(cfg, group_data_HC{:});
group_data_CC = ft_freqgrandaverage(cfg, group_data_CC{:});
group_data_UHR = ft_freqgrandaverage(cfg, group_data_UHR{:});

%Plot scalp maps of the grand average power for each group
cfg = [];
cfg.elec = subject.elec;
cfg.xlim = [30 45];
%cfg.ylim = [15 20];
cfg.zlim = [-1.9 -1.2];
cfg.baseline = 'no' ;%[0 0.5];
figure(41); ft_topoplotTFR(cfg,grandavg_UHR); colorbar; sgtitle('Ultra-High-Risk'); hold on;
figure(42); ft_topoplotTFR(cfg,grandavg_CC); colorbar; sgtitle('Clinical Control'); hold on;
figure(43); ft_topoplotTFR(cfg,grandavg_HC); colorbar; sgtitle('Healthy Control'); hold on;

% Get numerical indices of ROI to compute averages
sel_F_ROI = match_str(elec.label, F_ROI);
sel_P_ROI = match_str(elec.label, P_ROI);
sel_O_ROI = match_str(elec.label, O_ROI);
%sel_ONE_FRONT_ROI = match_str(elec.label, ONE_FRONT_ROI);
%sel_ONE_BACK_ROI = match_str(elec.label, ONE_BACK_ROI);

% % Normalize power spectra using mean power over freq range
freq_oi   = [30 45];   % frequency range of interest
% freq_norm = [30 45]; % frequency range used to normalize the spectrum
% foi_norm  = nearest(group_data_HC.freq, freq_norm);
% 
% common_denominator_HC = mean(group_data_HC.powspctrm(:,:,foi_norm(1):foi_norm(2)),3);
% common_denominator_CC = mean(group_data_CC.powspctrm(:,:,foi_norm(1):foi_norm(2)),3);
% common_denominator_UHR = mean(group_data_UHR.powspctrm(:,:,foi_norm(1):foi_norm(2)),3);
% 
% group_data_HC.powspctrm_b = group_data_HC.powspctrm./common_denominator_HC;
% group_data_CC.powspctrm_b = group_data_CC.powspctrm./common_denominator_CC;
% group_data_UHR.powspctrm_b = group_data_UHR.powspctrm./common_denominator_UHR;

% plot Neighbours

cfg = [];
cfg.method = 'distance'
cfg.elec = elec
neighbours = ft_prepare_neighbours(cfg)

% Collapse data using ROI and frequency ranges defined in relevant paper(s)
cfg = [];
cfg.channel     = F_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm'};
group_data_HC_F_ROI = ft_selectdata(cfg, group_data_HC);
group_data_CC_F_ROI = ft_selectdata(cfg, group_data_CC);
group_data_UHR_F_ROI = ft_selectdata(cfg, group_data_UHR);

cfg = [];
cfg.channel     = P_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm'};
group_data_HC_P_ROI = ft_selectdata(cfg, group_data_HC);
group_data_CC_P_ROI = ft_selectdata(cfg, group_data_CC);
group_data_UHR_P_ROI = ft_selectdata(cfg, group_data_UHR);

cfg = [];
cfg.channel     = O_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm'};
group_data_HC_O_ROI = ft_selectdata(cfg, group_data_HC);
group_data_CC_O_ROI = ft_selectdata(cfg, group_data_CC);
group_data_UHR_O_ROI = ft_selectdata(cfg, group_data_UHR);


% cfg.channel     = occipital_ROI;
% group_data_HC_oROI = ft_selectdata(cfg, group_data_HC);
% group_data_CC_oROI = ft_selectdata(cfg, group_data_CC);
% group_data_UHR_oROI = ft_selectdata(cfg, group_data_UHR);

% Plotspread function - may need to download from download server linked in
% % fieldtrip tutorial page 

data_raw_F_ROI = {group_data_HC_F_ROI.powspctrm,...
  group_data_CC_F_ROI.powspctrm,...
  group_data_UHR_F_ROI.powspctrm};

% 
% data_between_FL_ROI ={group_data_HC_FL_ROI.powspctrm_b...
%   group_data_CC_FL_ROI.powspctrm_b...
%   group_data_UHR_FL_ROI.powspctrm_b};
% 
% data_between_FM_ROI ={group_data_HC_FM_ROI.powspctrm_b...
%   group_data_CC_FM_ROI.powspctrm_b...
%   group_data_UHR_FM_ROI.powspctrm_b};
% 

% data_between_FR_ROI ={group_data_HC_FR_ROI.powspctrm_b...
%   group_data_CC_FR_ROI.powspctrm_b...
%   group_data_UHR_FR_ROI.powspctrm_b};

data_raw_P_ROI = {group_data_HC_P_ROI.powspctrm,...
  group_data_CC_P_ROI.powspctrm,...
  group_data_UHR_P_ROI.powspctrm};
% 
% data_between_PL_ROI ={group_data_HC_PL_ROI.powspctrm_b...
%   group_data_CC_PL_ROI.powspctrm_b...
%   group_data_UHR_PL_ROI.powspctrm_b};

% data_between_PM_ROI ={group_data_HC_PM_ROI.powspctrm_b...
%   group_data_CC_PM_ROI.powspctrm_b...
%   group_data_UHR_PM_ROI.powspctrm_b};
% 
% 
% data_between_PR_ROI ={group_data_HC_PR_ROI.powspctrm_b...
%   group_data_CC_PR_ROI.powspctrm_b...
%   group_data_UHR_PR_ROI.powspctrm_b};
% 
data_raw_O_ROI = {group_data_HC_O_ROI.powspctrm,...
  group_data_CC_O_ROI.powspctrm,...
  group_data_UHR_O_ROI.powspctrm};
% 
% data_between_OL_ROI ={group_data_HC_OL_ROI.powspctrm_b...
%   group_data_CC_OL_ROI.powspctrm_b...
%   group_data_UHR_OL_ROI.powspctrm_b};
% 

% data_between_OM_ROI ={group_data_HC_OM_ROI.powspctrm_b...
%   group_data_CC_OM_ROI.powspctrm_b...
%   group_data_UHR_OM_ROI.powspctrm_b};
% 
% data_between_OR_ROI ={group_data_HC_OR_ROI.powspctrm_b...
%   group_data_CC_OR_ROI.powspctrm_b...
%   group_data_UHR_OR_ROI.powspctrm_b};
% 
% data_between_OM_ROI ={group_data_HC_OR_ROI.powspctrm_b...
%   group_data_CC_OR_ROI.powspctrm_b...
%   group_data_UHR_OR_ROI.powspctrm_b};

% Make dot plots for ROIs 1-4 in different conditions as (x = condition, y =
% absolute power)

figure(1)

subplot(1,3,1);
plotSpread(data_raw_F_ROI,[],[],{'HC','CC','UHR'});
ylabel('abs. power (V^2)');
title('raw Frontal');
ylim([-2 0]);

subplot(1,3,2);
plotSpread(data_raw_P_ROI,[],[],{'HC','CC','UHR'});
ylabel('abs. power (V^2)');
title('raw Parietal');
ylim([-2 0]);

subplot(1,3,3);
plotSpread(data_raw_O_ROI,[],[],{'HC','CC','UHR'});
ylabel('abs. power (V^2)');
title('raw Occipital');
ylim([-2 0]);

% 
% subplot(2,3,2);
% h3 = plotSpread(data_between_OL_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
% ylabel('rel. power');
% title('between Occipital Left');
% set(h3{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

% subplot(2,3,4);
% h6 = plotSpread(data_between_OM_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
% ylabel('rel. power');
% title('between Occipital Middle');
% set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

% subplot(2,3,5);
% h6 = plotSpread(data_between_OR_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
% ylabel('rel. power');
% title('between Occipital Right');
% set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')
% 
% subplot(2,3,6);
% h6 = plotSpread(data_between_OR_ROI,[],[],{'Healthy Control','Clinical Control','Ultra-High Risk'});
% ylabel('rel. power');
% title('between Occipital Right');
% set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')


% Topoplots and power spectra for reach ROI 
cfg = [];
cfg.elec             = elec;
cfg.parameter        = 'powspctrm'; % you can plot either powspctrm (default) or powspctrm_b
cfg.xlim             = [0 45]; % frequency range to make the topoplot
cfg.highlight        = 'on';

% Options for improving appearance of plots (next 60 lines)

cfg.highlightchannel = {F_ROI P_ROI O_ROI};
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

subplot(2,3,1);loglog(group_data_HC.freq,...
  [squeeze(mean(mean(group_data_HC.(cfg.parameter)(:,sel_F_ROI,:),2),1))...
  squeeze(mean(mean(group_data_HC.(cfg.parameter)(:,sel_P_ROI,:),2),1))...
  squeeze(mean(mean(group_data_HC.(cfg.parameter)(:,sel_O_ROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
legend('F ROI', 'P_ROI', 'O_ROI', 'Location', 'southwest');

xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('Healthy Control');

subplot(2,3,2);loglog(group_data_CC.freq,...
   [squeeze(mean(mean(group_data_CC.(cfg.parameter)(:,sel_F_ROI,:),2),1))...
  squeeze(mean(mean(group_data_CC.(cfg.parameter)(:,sel_P_ROI,:),2),1))...
  squeeze(mean(mean(group_data_CC.(cfg.parameter)(:,sel_O_ROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('Clinical Control');

subplot(2,3,3);loglog(group_data_UHR.freq,...
   [squeeze(mean(mean(group_data_UHR.(cfg.parameter)(:,sel_F_ROI,:),2),1))...
  squeeze(mean(mean(group_data_UHR.(cfg.parameter)(:,sel_P_ROI,:),2),1))...
  squeeze(mean(mean(group_data_UHR.(cfg.parameter)(:,sel_O_ROI,:),2),1))]);
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

% cfg = [];
% cfg.parameter = {'powspctrm'};
% 
% % cfg.trials      = respon_group; % cfg.trials will select the 'subj' dimension
% % HC_respon = ft_selectdata(cfg, HC_respon);
% % CC_respon = ft_selectdata(cfg, CC_respon);
% % UHR_respon = ft_selectdata(cfg, UHR_respon);
% 
% % cfg.trials      = drowsy_group; % cfg.trials will select the 'subj' dimension
% % HC_drowsy = ft_selectdata(cfg, HC_drowsy);
% % CC_drowsy = ft_selectdata(cfg, CC_drowsy);
% % UHR_drowsy = ft_selectdata(cfg, UHR_drowsy);
% 
% % Then, the permutation test
% % Differences:
% % 1) Different sample level effect statistic
% % 2) Different design matrix
% % 3) Independent variable = group assignment
% % Using independent samples T-statistic (BETWEEN subjects)
% % Design matrix includes group assignment 

cfg = [];
cfg.channel          = 'CZ';
cfg.avgovergchan     = 'no';
cfg.frequency        = [30 45];
cfg.avgovergfreq     = 'yes';
cfg.parameter        = 'powspctrm';
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
% 
 design = ones(1,size(group_data_HC_all,1) + size(group_data_CC_all,1) + size(group_data_UHR_all,1))*3;
 design(1,1:size(group_data_HC_all,1)) = 1;
 design(1,(size(group_data_HC_all,1)+1):(size(group_data_HC_all,1)+size(group_data_CC_all,1))) = 2;
% 
 cfg.design = design;
 cfg.ivar   = 1;
% 
% % Use configuration to perform analysis 
% 
 stat2 = ft_freqstatistics(cfg, group_data_HC, group_data_CC, group_data_UHR);
% % stat2 = ft_freqstatistics(cfg, group_data_HC, group_data_CC);
% % stat3 = ft_freqstatistics(cfg, group_data_CC, group_data_UHR);
% % stat4 = ft_freqstatistics(cfg, group_data_UHR, group_data_HC);
% 
% % Extract the data manually from ONE (or multiple) ELECTRODES and run
% stats %%% 30042024 MAYBE RECOMMENT (To run anova)
% % on the power spectrum within the defined frequency bands
 foi = [30 45]; %frequency bands of interest
 foi_idx = nearest(group_data_CC.freq(:),foi)
 
    electrodeName = {'F7','F5','F3', 'FP1', 'AF3', 'F1', 'FPZ', 'FP2', 'AF4', 'F2', 'FZ','F4','F6','F8'}; %this region (here, frontal)
     for i = 1:length(electrodeName);
         thiselec = electrodeName(i);
        
 elecIdx(i) = find(strcmp(elec(:).label,thiselec)==1)
     end

% %Loop through all participants and extract the data
% % HEALTHY CONTROLS FIRST
 pwr_data_group_HC = squeeze(mean(group_data_HC.powspctrm(:,elecIdx,foi_idx(1):foi_idx(2)),2))
 m_pwr_data_group_HC = mean(pwr_data_group_HC,2)
 pwr_data_group_CC = squeeze(mean(group_data_CC.powspctrm(:,elecIdx,foi_idx(1):foi_idx(2)),2))
 m_pwr_data_group_CC = mean(pwr_data_group_CC,2)
 pwr_data_group_UHR = squeeze(mean(group_data_UHR.powspctrm(:,elecIdx,foi_idx(1):foi_idx(2)),2))
 m_pwr_data_group_UHR = mean(pwr_data_group_UHR,2)
% 
%Run ANOVA1
group = [ones(size(pwr_data_group_HC,1),1);...
    (ones(size(pwr_data_group_CC,1),1)*2);ones(size(pwr_data_group_UHR,1),1)*3]
[p,anovatab,anovastats] = anova1([m_pwr_data_group_HC;m_pwr_data_group_CC;m_pwr_data_group_UHR],group)
%Now run stats on these data
[h,p,tteststats,CI] = ttest2(m_pwr_data_group_HC,m_pwr_data_group_UHR)

% Plot the results!!!!

cfg = [];
%cfg.alpha    = stat2.cfg.alpha;
cfg.alpha = 0.05;
% JA 30042024 edit = make cfg.alpha 0.05 - is this ok? 
%cfg.parameter = 'stat';
cfg.parameter = 'stat';
cfg.zlim      = [-3 3];
cfg.elec      = elec;
%ft_clusterplot(cfg, stat2);
ft_clusterplot(cfg, stat2);

cfg = [];
cfg.elec         = elec;
%cfg.zlim         = [-1 -0.8];
cfg.xlim         = [30 45];
cfg.parameter    = 'powspctrm';
cfg.markersymbol = '.';
cfg.comment      = 'no';
cfg.colormap     = 'jet';
cfg.colorbar     = 'no';

figure('position',[680 240 1039 420]);
figure; ft_topoplotER(cfg, group_data_HC); colorbar; title('Healthy Control [30 45]');
figure; ft_topoplotER(cfg, group_data_CC); colorbar; title('Clinical Control [30 45]');
figure; ft_topoplotER(cfg, group_data_UHR); colorbar; title('Ultra-High Risk [30 45]');

% 
% % Plot power function for each group as a function of sedative state 
% % % First, for frontal ROI 
% % 
% % figure;
% % 
% % subplot(2,4,1); loglog(HC_respon.freq,...
% %   [squeeze(mean(mean(HC_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
% %   squeeze(mean(mean(HC_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
% % xlim([0.5 45]);
% % grid on; hold on;
% % plot([10,10],[10^-3 10^2],'--k')
% % ylim([10^-3 10^2]);
% % legend('Front ROI resp','Front ROI drow','Location','southwest');
% % xlabel('Frequency (Hz)');
% % ylabel(cfg.parameter);
% % title('baseline');
% % 
% % subplot(2,4,2); loglog(CC_respon.freq,...
% %   [squeeze(mean(mean(CC_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
% %   squeeze(mean(mean(CC_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
% % xlim([0.5 45]);
% % grid on; hold on;
% % plot([10,10],[10^-3 10^2],'--k')
% % ylim([10^-3 10^2]);
% % xlabel('Frequency (Hz)');
% % ylabel(cfg.parameter);
% % title('mild');
% % 
% % subplot(2,4,3); loglog(UHR_respon.freq,...
% %   [squeeze(mean(mean(UHR_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
% %   squeeze(mean(mean(UHR_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
% % xlim([0.5 45]);
% % grid on; hold on;
% % plot([10,10],[10^-3 10^2],'--k')
% % ylim([10^-3 10^2]);
% % xlabel('Frequency (Hz)');
% % ylabel(cfg.parameter);
% % title('moderate');
% % 
% % subplot(2,4,4); loglog(reco_sedation_respon.freq,...
% %   [squeeze(mean(mean(reco_sedation_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
% %   squeeze(mean(mean(reco_sedation_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
% % xlim([0.5 45]);
% % grid on; hold on;
% % plot([10,10],[10^-3 10^2],'--k')
% % ylim([10^-3 10^2]);
% % xlabel('Frequency (Hz)');
% % ylabel(cfg.parameter);
% % title('recover');
% % 
% % % Then, the occipital ROI: 
% % 
% % subplot(2,4,5); loglog(HC_respon.freq,...
% %   [squeeze(mean(mean(HC_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
% %   squeeze(mean(mean(HC_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
% % xlim([0.5 45]);
% % grid on; hold on;
% % plot([10,10],[10^-3 10^2],'--k')
% % ylim([10^-3 10^2]);
% % xlabel('Frequency (Hz)');
% % ylabel(cfg.parameter);
% % legend('Occip ROI resp','Occip ROI drow','Location','southwest');
% % 
% % subplot(2,4,6); loglog(CC_respon.freq,...
% %   [squeeze(mean(mean(CC_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
% %   squeeze(mean(mean(CC_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
% % xlim([0.5 45]);
% % grid on; hold on;
% % plot([10,10],[10^-3 10^2],'--k')
% % ylim([10^-3 10^2]);
% % xlabel('Frequency (Hz)');
% % ylabel(cfg.parameter);
% % 
% % subplot(2,4,7); loglog(UHR_respon.freq,...
% %   [squeeze(mean(mean(UHR_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
% %   squeeze(mean(mean(UHR_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
% % xlim([0.5 45]);
% % grid on; hold on;
% % plot([10,10],[10^-3 10^2],'--k')
% % ylim([10^-3 10^2]);
% % xlabel('Frequency (Hz)');
% % ylabel(cfg.parameter);
% % 
% % subplot(2,4,8); loglog(reco_sedation_respon.freq,...
% %   [squeeze(mean(mean(reco_sedation_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
% %   squeeze(mean(mean(reco_sedation_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
% % xlim([0.5 45]);
% % grid on; hold on;
% % plot([10,10],[10^-3 10^2],'--k')
% % ylim([10^-3 10^2]);
% % xlabel('Frequency (Hz)');
% % ylabel(cfg.parameter);
% 
% % Or, ANOVA to copmare all three groups 
% 
% cfg = [];
% cfg.channel          = 'all';
% cfg.frequency        = foi_contrast;
% cfg.avgovergfreq     = 'no';
% cfg.parameter        = 'powspctrm';
% cfg.method           = 'ft_statistics_montecarlo';
% cfg.statistic        = 'ft_statfun_depsamplesFmultivariate';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum'; %'maxsum', 'maxsize', 'wcm'
% cfg.clusterthreshold = 'nonparametric_common';
% cfg.minnbchan        = 2;
% cfg.tail             = 1; % For a F-statistic, it only make sense to calculate the right tail
% cfg.clustertail      = cfg.tail;
% cfg.alpha            = 0.05;
% cfg.computeprob      = 'yes';
% cfg.numrandomization = 500;
% cfg.neighbours       = cfg_neigh.neighbours;
% 
% nsubj = size(covariates,1);
% design = zeros(2,4*nsubj);
% design(1,1:nsubj)           = 1;
% design(1,nsubj+1:2*nsubj)   = 2;
% design(1,nsubj*2+1:3*nsubj) = 3;
% design(1,nsubj*3+1:4*nsubj) = 4;
% design(2,:) = repmat(1:nsubj,1,4);
% 
% cfg.design = design;
% cfg.ivar   = 1; % sedation level
% cfg.uvar   = 2; % subject number
% 
% % Pass the data from all four conditions as input variables 
% 
% stat3 = ft_freqstatistics(cfg, HC, CC, UHR, reco_sedation);
% 
% % Plot them results
% cfg            = [];
% cfg.frequency  = foi_contrast;
% cfg.avgoverrpt = 'yes';
% cfg.parameter  = {'powspctrm'};
% HC_avg = ft_selectdata(cfg, HC);
% CC_avg = ft_selectdata(cfg, CC);
% UHR_avg = ft_selectdata(cfg, UHR);
% %reco_sedation_avg = ft_selectdata(cfg, reco_sedation);
% 
% % copy the mask field to each variable
% HC_avg.mask = stat3.mask;
% CC_avg.mask = stat3.mask;
% UHR_avg.mask = stat3.mask;
% %reco_sedation_avg.mask = stat3.mask;
% 
% cfg = [];
% cfg.zlim          = [0 90];
% cfg.elec          = elec;
% cfg.colorbar      = 'no';
% cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
% cfg.maskstyle     = 'box';
% cfg.parameter     = 'powspctrm';
% cfg.maskfacealpha = 0.1;
% figure; ft_multiplotER(cfg, HC_avg, CC_avg, UHR_avg, reco_sedation_avg);