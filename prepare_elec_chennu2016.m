function [elec cfg] = prepare_elec_chennu2016(label_mixed,neigh)
% this function links Chennu et al 2016 electrode labels with Fieldtrip
% GSN-HydroCel-129.sfp electrode template
% see this thread to know some position equivalences between EGI and 10-10
% https://www.researchgate.net/publication/266609828_Determination_of_the_Geodesic_Sensor_Nets'_Average_Electrode_Positions_and_Their_10_-_10_International_Equivalents
if nargin < 2;
  neigh = 0;
  cfg = [];
end

% this is the look-up table of elecctrode position equivalence
GSNto1010 = {...
  'Fp2' 'E9'; 'Fz' 'E11'; 'Fp1' 'E22'; 'F3' 'E24'; 'F7' 'E33'; 'C3' 'E36'; ...
  'T3' 'E45'; 'P3' 'E52'; 'T5'  'E58'; 'Pz' 'E62'; 'O1' 'E70'; 'Oz' 'E75';...
  'O2' 'E83'; 'P4' 'E92'; 'T6'  'E96'; 'C4' 'E104';'T4' 'E108';'F8' 'E122';...
  'F4' 'E124'; 'POz' 'E72'};

% 10-10 electrodes inside labels
e1010 = ft_channelselection('EEG',label_mixed);

% index look-up table GSNto1010 as a function of input label
[r1,r2] = match_str(e1010,GSNto1010(:,1));

% access to the Fieldtrip GSN template
ft_path = fileparts(which('ft_defaults'));
elec = ft_read_sens(fullfile(ft_path,'template','electrode','GSN-HydroCel-129.sfp'));

% substitute the 10-10 system labels by its GSN equivalent to get the GSN
% coordinates
gsn_r2 = GSNto1010(r2,:);
[s1,s2] = match_str(label_mixed,gsn_r2(:,1));
label = label_mixed;
label(s1) = gsn_r2(s2,2);

% find the GSN system-specific sensors
[sel1,sel2] = match_str(label,elec.label);

chanpos  = elec.chanpos(sel2,:);
elecpos  = elec.elecpos(sel2,:);
chantype = elec.chantype(sel2,:);
chanunit = elec.chanunit(sel2,:);
elec.chantype = chantype;
elec.chanunit = chanunit;
elec.chanpos  = chanpos;
elec.elecpos  = elecpos;
elec.label    = label_mixed;

if neigh
  cfg = [];
    cfg.method = 'distance';
    cfg.neighbourdist = 4;
    cfg.elec = elec;
    cfg.feedback = 'no';
    neigh = ft_prepare_neighbours(cfg);
    
    cfg = [];
    cfg.neighbours = neigh;
    cfg.elec = elec;
    cfg.enableedit = 'yes';
    cfg = ft_neighbourplot(cfg);
end