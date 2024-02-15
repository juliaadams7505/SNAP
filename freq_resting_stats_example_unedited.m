% load averages for each individual subject, for each condition
load freq_resting.mat

% Defining regions of interest
occipital_ROI = {'E50','T5','E59','E60','Pz','E65','E66','E67','O1','E71','E72','Oz','E76','E77','O2','E84','E85','E90','E91','T6','E101','E51','E97'};
frontal_ROI   = {'E3','E4','E5','E6','E7','Fp2','E10','Fz','E12','E13','E15','E16','E18','E19','E20','Fp1','E23','F3','E27','E28','E29','E30','E105','E106','E111','E112','E117','E118','E123','F4'};

cfg = [];
cfg.keepindividual = 'yes';
base_sedation = ft_freqgrandaverage(cfg, base_sedation{:});
mild_sedation = ft_freqgrandaverage(cfg, mild_sedation{:});
mode_sedation = ft_freqgrandaverage(cfg, mode_sedation{:});
reco_sedation = ft_freqgrandaverage(cfg, reco_sedation{:});

% Get numerical indices of ROI to compute averages
sel_oROI = match_str(base_sedation.label, occipital_ROI);
sel_fROI = match_str(base_sedation.label, frontal_ROI);

% And load the electrode positions
elec = prepare_elec_chennu2016(base_sedation.label);

% Normalize power spectra using mean power over freq range
freq_oi   = [8 15];   % frequency range of interest
freq_norm = [0.7 40]; % frequency range used to normalize the spectrum
foi_norm  = nearest(base_sedation.freq, freq_norm);

common_denominator = mean(base_sedation.powspctrm(:,:,foi_norm(1):foi_norm(2)),3);
base_sedation.powspctrm_b = bsxfun(@rdivide, base_sedation.powspctrm, common_denominator);    
mild_sedation.powspctrm_b = bsxfun(@rdivide, mild_sedation.powspctrm, common_denominator);
mode_sedation.powspctrm_b = bsxfun(@rdivide, mode_sedation.powspctrm, common_denominator);
reco_sedation.powspctrm_b = bsxfun(@rdivide, reco_sedation.powspctrm, common_denominator);

% Collapse data using ROI and frequency ranges defined in relevant paper(s)
cfg = [];
cfg.channel     = frontal_ROI;
cfg.avgoverchan = 'yes';
cfg.frequency   = freq_oi;
cfg.avgoverfreq = 'yes';
cfg.parameter   = {'powspctrm','powspctrm_b'};
base_sedation_fROI = ft_selectdata(cfg, base_sedation);
mild_sedation_fROI = ft_selectdata(cfg, mild_sedation);
mode_sedation_fROI = ft_selectdata(cfg, mode_sedation);
reco_sedation_fROI = ft_selectdata(cfg, reco_sedation);

cfg.channel     = occipital_ROI;
base_sedation_oROI = ft_selectdata(cfg, base_sedation);
mild_sedation_oROI = ft_selectdata(cfg, mild_sedation);
mode_sedation_oROI = ft_selectdata(cfg, mode_sedation);
reco_sedation_oROI = ft_selectdata(cfg, reco_sedation);

% Plotspread function - may need to download from download server linked in
% fieldtrip tutorial page 
 
data_raw_fROI    = {base_sedation_fROI.powspctrm...
  mild_sedation_fROI.powspctrm...
  mode_sedation_fROI.powspctrm...
  reco_sedation_fROI.powspctrm};
data_raw_oROI    = {base_sedation_oROI.powspctrm...
  mild_sedation_oROI.powspctrm...
  mode_sedation_oROI.powspctrm...
  reco_sedation_oROI.powspctrm};

data_between_fROI ={base_sedation_fROI.powspctrm_b...
  mild_sedation_fROI.powspctrm_b...
  mode_sedation_fROI.powspctrm_b...
  reco_sedation_fROI.powspctrm_b};
data_between_oROI ={base_sedation_oROI.powspctrm_b...
  mild_sedation_oROI.powspctrm_b...
  mode_sedation_oROI.powspctrm_b...
  reco_sedation_oROI.powspctrm_b};

% Make dot plots for ROIs 1-4 in different conditions as (x = condition, y =
% absolute power)

figure

subplot(2,2,1);
plotSpread(data_raw_fROI,[],[],{'baseline','mild','moderate','recovery'});
ylabel('abs. power (V^2)');
title('raw PSD Front');

subplot(2,2,2);
h3 = plotSpread(data_between_fROI,[],[],{'baseline','mild','moderate','recovery'});
ylabel('rel. power');
title('between PSD Front');
set(h3{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

subplot(2,2,3);
plotSpread(data_raw_oROI,[],[],{'baseline','mild','moderate','recovery'});
ylabel('abs. power (V^2)');
title('raw PSD Occip');

subplot(2,2,4);
h6 = plotSpread(data_between_oROI,[],[],{'baseline','mild','moderate','recovery'});
ylabel('rel. power');
title('between PSD Occip');
set(h6{1},'LineWidth',1,'Marker', '.','Color','k','MarkerFaceColor','k')

% Topoplots and power spectra for reach ROI in each sedative condition
cfg = [];
cfg.elec             = elec;
cfg.parameter        = 'powspctrm_b'; % you can plot either powspctrm (default) or powspctrm_b
cfg.xlim             = [8 15]; % frequency range to make the topoplot
cfg.highlight        = 'on';

% Options for improving appearance of plots (next 60 lines)

cfg.highlightchannel = {frontal_ROI occipital_ROI};
cfg.highlightsymbol  = {'o','*'};
cfg.highlightcolor   = [0 0 0];
cfg.highlightsize    = 6;
cfg.markersymbol     = '.';
cfg.comment          = 'no';
cfg.colormap         = 'jet';

figure('position',[680 240 1039 420]);
subplot(2,4,1); ft_topoplotER(cfg, base_sedation); colorbar; title('baseline');
subplot(2,4,2); ft_topoplotER(cfg, mild_sedation); colorbar; title('mild');
subplot(2,4,3); ft_topoplotER(cfg, mode_sedation); colorbar; title('moderate');
subplot(2,4,4); ft_topoplotER(cfg, reco_sedation); colorbar; title('recovery');

subplot(2,4,5);loglog(base_sedation.freq,...
  [squeeze(mean(mean(base_sedation.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(base_sedation.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
legend('Front ROI','Occip ROI','Location','southwest');
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('baseline');

subplot(2,4,6);loglog(mild_sedation.freq,...
  [squeeze(mean(mean(mild_sedation.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(mild_sedation.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('mild');

subplot(2,4,7);loglog(mode_sedation.freq,...
  [squeeze(mean(mean(mode_sedation.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(mode_sedation.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('moderate');

subplot(2,4,8);loglog(reco_sedation.freq,...
  [squeeze(mean(mean(reco_sedation.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(reco_sedation.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('recovery');

%%%% BETWEEN-PARTICIPANT CONTRASTS %%%%%%

% Statistically evaluate power difference between groups (they have

[base_sedation_respon base_sedation_drowsy] = deal(base_sedation);
[mild_sedation_respon mild_sedation_drowsy] = deal(mild_sedation);
[mode_sedation_respon mode_sedation_drowsy] = deal(mode_sedation);
[reco_sedation_respon reco_sedation_drowsy] = deal(reco_sedation);

cfg = [];
cfg.parameter   = {'powspctrm','powspctrm_b'};

cfg.trials      = respon_group; % cfg.trials will select the 'subj' dimension
base_sedation_respon = ft_selectdata(cfg, base_sedation_respon);
mild_sedation_respon = ft_selectdata(cfg, mild_sedation_respon);
mode_sedation_respon = ft_selectdata(cfg, mode_sedation_respon);
reco_sedation_respon = ft_selectdata(cfg, reco_sedation_respon);

cfg.trials      = drowsy_group; % cfg.trials will select the 'subj' dimension
base_sedation_drowsy = ft_selectdata(cfg, base_sedation_drowsy);
mild_sedation_drowsy = ft_selectdata(cfg, mild_sedation_drowsy);
mode_sedation_drowsy = ft_selectdata(cfg, mode_sedation_drowsy);
reco_sedation_drowsy = ft_selectdata(cfg, reco_sedation_drowsy);

% Then, the permutation test
% Differences:
% 1) Different sample level effect statistic
% 2) Different design matrix
% 3) Independent variable = group assignment
% Using independent samples T-statistic (BETWEEN subjects)
% Design matrix includes group assignment 

cfg                 = [];
cfg.channel          = 'all';
cfg.avgovergchan     = 'no';
cfg.frequency        = foi_contrast;
cfg.avgovergfreq     = 'yes';
cfg.parameter        = 'powspctrm_b';
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.computeprob      = 'yes';
cfg.numrandomization = 1000;
cfg.neighbours       = cfg_neigh.neighbours;

design = zeros(1,size(respon_group,1) + size(drowsy_group,1));
design(1,1:size(respon_group,1)) = 1;
design(1,(size(respon_group,1)+1):(size(respon_group,1)+size(drowsy_group,1))) = 2;

cfg.design = design;
cfg.ivar   = 1;

% Use configuration to perform analysis 

stat2 = ft_freqstatistics(cfg, reco_sedation_respon, reco_sedation_drowsy);

% Plot them results!!!!

cfg = [];
cfg.alpha     = stat2.cfg.alpha;
cfg.parameter = 'stat';
cfg.zlim      = [-3 3];
cfg.elec      = elec;
ft_clusterplot(cfg, stat2);

cfg = [];
cfg.elec         = elec;
cfg.zlim         = [1.5 3];
cfg.xlim         = [8 15];
cfg.parameter    = 'powspctrm_b';
cfg.markersymbol = '.';
cfg.comment      = 'no';
cfg.colormap     = 'jet';
cfg.colorbar     = 'no';

figure('position',[680 240 1039 420]);
subplot(2,4,1); ft_topoplotER(cfg, base_sedation_respon); colorbar; title('base Responsive');
subplot(2,4,2); ft_topoplotER(cfg, mild_sedation_respon); colorbar; title('mild Responsive');
subplot(2,4,3); ft_topoplotER(cfg, mode_sedation_respon); colorbar; title('mode Responsive');
subplot(2,4,4); ft_topoplotER(cfg, reco_sedation_respon); colorbar; title('reco Responsive');

subplot(2,4,5); ft_topoplotER(cfg, base_sedation_drowsy); colorbar; title('base Drowsy');
subplot(2,4,6); ft_topoplotER(cfg, mild_sedation_drowsy); colorbar; title('mild Drowsy');
subplot(2,4,7); ft_topoplotER(cfg, mode_sedation_drowsy); colorbar; title('mode Drowsy');
subplot(2,4,8); ft_topoplotER(cfg, reco_sedation_drowsy); colorbar; title('reco Drowsy');

% Plot power function for each group as a function of sedative state 
% First, for frontal ROI 

figure;

subplot(2,4,1); loglog(base_sedation_respon.freq,...
  [squeeze(mean(mean(base_sedation_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(base_sedation_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
legend('Front ROI resp','Front ROI drow','Location','southwest');
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('baseline');

subplot(2,4,2); loglog(mild_sedation_respon.freq,...
  [squeeze(mean(mean(mild_sedation_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(mild_sedation_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
title('mild');

subplot(2,4,3); loglog(mode_sedation_respon.freq,...
  [squeeze(mean(mean(mode_sedation_respon.(cfg.parameter)(:,sel_fROI,:),2),1))...
  squeeze(mean(mean(mode_sedation_drowsy.(cfg.parameter)(:,sel_fROI,:),2),1))]);
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

subplot(2,4,5); loglog(base_sedation_respon.freq,...
  [squeeze(mean(mean(base_sedation_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
  squeeze(mean(mean(base_sedation_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);
legend('Occip ROI resp','Occip ROI drow','Location','southwest');

subplot(2,4,6); loglog(mild_sedation_respon.freq,...
  [squeeze(mean(mean(mild_sedation_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
  squeeze(mean(mean(mild_sedation_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
xlim([0.5 45]);
grid on; hold on;
plot([10,10],[10^-3 10^2],'--k')
ylim([10^-3 10^2]);
xlabel('Frequency (Hz)');
ylabel(cfg.parameter);

subplot(2,4,7); loglog(mode_sedation_respon.freq,...
  [squeeze(mean(mean(mode_sedation_respon.(cfg.parameter)(:,sel_oROI,:),2),1))...
  squeeze(mean(mean(mode_sedation_drowsy.(cfg.parameter)(:,sel_oROI,:),2),1))]);
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

stat3 = ft_freqstatistics(cfg, base_sedation, mild_sedation, mode_sedation, reco_sedation);

% Plot them results
cfg            = [];
cfg.frequency  = foi_contrast;
cfg.avgoverrpt = 'yes';
cfg.parameter  = {'powspctrm','powspctrm_b'};
base_sedation_avg = ft_selectdata(cfg, base_sedation);
mild_sedation_avg = ft_selectdata(cfg, mild_sedation);
mode_sedation_avg = ft_selectdata(cfg, mode_sedation);
reco_sedation_avg = ft_selectdata(cfg, reco_sedation);

% copy the mask field to each variable
base_sedation_avg.mask = stat3.mask;
mild_sedation_avg.mask = stat3.mask;
mode_sedation_avg.mask = stat3.mask;
reco_sedation_avg.mask = stat3.mask;

cfg = [];
cfg.zlim          = [0 90];
cfg.elec          = elec;
cfg.colorbar      = 'no';
cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle     = 'box';
cfg.parameter     = 'powspctrm_b';
cfg.maskfacealpha = 0.1;
figure; ft_multiplotER(cfg, base_sedation_avg, mild_sedation_avg, mode_sedation_avg, reco_sedation_avg);

