% function BIDS_BadChanRejectCAR(sub, ses, run, project_name, main_dir, ieegsource_dir, event_dir, diagnostics_dir, manualVar, pathstr)

function BIDS_BadChanRejectCAR(sub, ses, run, project_name, main_dir, ...
    ieegsource_dir, event_dir, diagnostics_dir, manualVar, pathstr, ...
    use_interactive, badchan_file)
% BIDS_BadChanRejectCAR(sub, ses, run, project_name, main_dir, manualVar, pathstr)
% 1. notch filtering: Removes line noise at 60 Hz and its harmonics. `noiseFiltData`
% 2. Spike Detection and Removal: Identifies channels with excessive spikes. nr_jumps
% 3. Power Spectrum-Based Filtering: Purpose: Identifies channels with arunormal power spectra. `pwelch`
% 4. Baseline Adjustment: CAR = mean(data_good, 1);
% 5. Epileptic HFO Detection: Identifies pathological channels with high-frequency oscillations (HFOs). `find_paChan`

use_interactive = false;
% badchan_file = '/data/manual_badchans/sub-01_ses-01_run-01.tsv';
% make dirs
if ~isfolder(fullfile(diagnostics_dir, sub))
    mkdir(fullfile(diagnostics_dir, sub))
end

% load global metadata file
global_dir = fullfile(ieegsource_dir, sub);
global_fpath = fullfile(global_dir, sprintf('%s_%s_task-%s_%s_global.mat', sub, ses, project_name, run));

if exist(global_fpath, 'file')
    tmp = load(global_fpath);
    if isfield(tmp, 'globalVar')
        globalVar = tmp.globalVar;
    else
        error('globalVar not found in %s', global_fpath);
    end
else
    error('File does not exist: %s', global_fpath);
end

%% Keep track of manually bad channels
%% DEP: I have new code for this: 
% code/LBCNpipeline_BIDS/BIDS_step00_badchannels_01template.py
% code/LBCNpipeline_BIDS/BIDS_step00_badchannels_02populate.py
% manualVar = ManuallyAddBadChannel(sub, project_name, run);

% TODO: 
if isfield(globalVar, 'refChan')
    % Use what's already in globalVar
    disp('Using existing globalVar.refChan and related fields.');
    globalVar.oobChan = [];
    globalVar.trigNoiseChan = [];
else
    % set these fields as empty if doesn't pre-exist in globalVar
    globalVar.refChan = [];
    globalVar.epiChan = [];
    globalVar.emptyChan = [];
    globalVar.oobChan = [];
    globalVar.trigNoiseChan = [];
end

%skip bad chan reject CAR if already computed
if isfield(globalVar,'badChan')
    fprintf('Common average already computed for block: %s...moving on\n',run)
    % continue
end

%Chan info PROMPT
globalVar.badChan = [globalVar.refChan globalVar.epiChan globalVar.emptyChan globalVar.oobChan globalVar.trigNoiseChan];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Bad channel detection based on raw power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This part saves data.mat

% Format the data
% data example: ./emotionalappraisal/results/ieeg-source/sub-01/ses-01/run-01/ieeg/sub-01_ses-01_task-emotionalappraisal_run-01_electrode-002_ieeg.mat
% sub_info.data.name = fullfile(globalVar.originalData,strcat('iEEG',run,'_'));
% sub_info.Pdiode.name = fullfile(globalVar.originalData, strcat('Pdio',run,'_02'));
sub_info.data.name = fullfile(globalVar.originalData,ses, run, 'ieeg', sprintf('%s_%s_task-%s_%s_electrode-',sub, ses, project_name, run));
sub_info.Pdiode.name = fullfile(globalVar.originalData, ses, run, 'pdio', sprintf('%s_%s_task-%s_%s_electrode-001_pdio', sub, ses, project_name, run));

data_nr=1;
load([sub_info.data(data_nr).name '001_ieeg.mat']) % HJ TODO: What does this do?
ttot=size(wave,2); % ttot = total time in units of samples
clear wave


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Load data into matrix convert to fieldtrip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
% data_mat = LoadPreprocDataMatrix_HJ(sub,project_name,run,'originalData',0,[],dirs, pathstr);
ieegsource_dir = fileparts(globalVar.originalData);
data_mat = BIDS_LoadPreprocDataMatrix_HJ(sub, ses, project_name, run, 'originalData',...
                                        0,[],main_dir, ieegsource_dir, event_dir, pathstr);
%convert to fieldtrip
data.label = data_mat.chans;
data.trial = {data_mat.wave};
data.sfreq = data_mat.fsample;
data.time = {data_mat.time};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Fieldtrip channel rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%split up segments
%{
    cfg = [];
    cfg.length = 1;
    data_segments = ft_redefinetrial(cfg, data);
    %remove noisy channels
    warning('Choosing Channels to reject...')
    cfg = [];
    cfg.method = 'summary';
    data_out = ft_rejectvisual(cfg,data_segments);
    
    %identify indexes
    bad_chan = setdiff(data_segments.label,data_out.label);
    bad_chan_ind = find(ismember(data_segments.label,bad_chan));
    good_chan_ind = find(~ismember(data_segments.label,bad_chan));
    %clean up workspace
    clear data_out data_segments
%}
data_all = data.trial{1}';

%add the peak array to the global variable
globalVar.peak_array = [];
data = data_all;

%filter out line noise
for n=1:globalVar.nchan
    disp(['notch filtering electrode ' num2str(n) ' out of ' num2str(globalVar.nchan)])
    %Filtering 60 Hz line noise & harmonics
    wave = noiseFiltData(globalVar, data_all(:,n)'); % LBCN_PREPROC
    wave = wave.'; %because has incorrect ordering for efficiency, should be time x 1
    data(:,n)=-wave; % invert: data are recorded inverted
    clear wave nl
end

data_all = data;
data=single(data);

% Run algorithm


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                STEP 1: based on deviations on the raw signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bad_chans=[globalVar.refChan globalVar.badChan globalVar.epiChan globalVar.emptyChan];
disp(['Bad electrodes inputed by user:' int2str(bad_chans)]);

a=var(data);
b=find(a>(5*median(a))); % 5 * higher than median.
c=find(a<(median(a)/5)); % 5 * lower than median.
bad_cha_tmp = [b c];

if ~isempty([b c])
    disp(['Bad electrodes based on raw power: ' int2str(bad_cha_tmp)]);
end

% %Plot bad channels
figureDim = [0 0 .5 .5];
f1 = figure('units', 'normalized', 'outerposition', figureDim);
subplot(2,1,1)
for ii = bad_cha_tmp
    hold on
    plot(zscore(data(:,ii)) + double(ii));
end
title('Bad electrodes based on raw power')
xlabel('Time (s)')
ylabel('Electrode number')
hold off

% Update the globalVar.badChan
globalVar.badChan = unique([bad_cha_tmp globalVar.badChan]);
globalVar.badChan_specs.deviation_raw = bad_cha_tmp;

% remove bad channels
chan_lbls=1:size(data,2);
data(:,globalVar.badChan)=[];
chan_lbls(globalVar.badChan)=[];



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            STEP 2: based on the number of spikes on the raw signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr_jumps=zeros(1,size(data,2));
for k=1:size(data,2)
    nr_jumps(k)=length(find(diff(data(:,k))>80)); % find jumps>80uV
end

subplot(2,1,2)
plot(nr_jumps);
ylabel('Number of spikes')
xlabel('Electrode number')
figure(f1)

% SAVE FIGURE
saveas(f1, fullfile(diagnostics_dir, sub, sprintf('%s_%s_task-%s_%s_spikes.png', sub, ses, project_name, run)));

% 1 jump per second in average
ej= floor(globalVar.chan_length/globalVar.iEEG_rate);% 1 jump per second in average
jm=find(nr_jumps>ej);% Find channels with more than ... jumps
clcl=chan_lbls(jm);% Real channel numbers of outliers
disp(['spiky electrodes: ' int2str(clcl)]);

% Update globalVar.badChan
globalVar.badChan = unique([clcl globalVar.badChan]);
globalVar.badChan_specs.number_of_spikes = clcl;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              STEP 3: based on deviations on the power spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_ov = 0; 
f = 0:250; 
data_pxx=zeros(length(f),size(data,2));

for k = 1:size(data,2)
    t_len = min(floor(100*globalVar.iEEG_rate),size(data,1));
    [Pxx,f] = pwelch(data(1:t_len,k),floor(globalVar.iEEG_rate),set_ov,f,floor(globalVar.iEEG_rate));
    data_pxx(:,k)=Pxx;
end

% Plot chnas in different colors:
plotthis=log(data_pxx);
figureDim = [0 0 .5 1];
f2 = figure('units', 'normalized', 'outerposition', figureDim);
for fi = 1:size(plotthis,2)
    hold on
    plot(f,plotthis(:,fi))
    text(0:20:size(plotthis,1), plotthis(1:20:size(plotthis,1),fi), num2str(chan_lbls(fi)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end
xlim([1 size(plotthis,1)])
ylim([min(plotthis(:)) max(plotthis(:))])
xlabel('Frequency')
ylabel('Power')
hold off
figure(f2);

datacursormode on;
saveas(f2, fullfile(diagnostics_dir, sub, sprintf('%s_%s_task-%s_%s_powerspectrum.png', sub, ses, project_name, run)));
saveas(f2, fullfile(diagnostics_dir, sub, sprintf('%s_%s_task-%s_%s_powerspectrum.fig', sub, ses, project_name, run)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Prompt for bad channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~strcmp(project_name, 'Rest')
%     bad_chan_spec = promptBadChanSpec; % of the remaining channels
% else
%     bad_chan_spec = [];
% end
%%%%%% input interactive reports. 
% if exist('use_interactive', 'var') && use_interactive
%     bad_chan_spec = promptBadChanSpec;
% elseif exist('badchan_file', 'var') && isfile(badchan_file)
%     T = readtable(badchan_file, 'FileType', 'text', 'Delimiter', '\t');
%     bad_chan_spec = T.badChanIndex;
% else
%     bad_chan_spec = [];
% end
if exist('use_interactive', 'var') && use_interactive
    disp('Interactive mode enabled — prompting for bad channel selection...');
    bad_chan_spec = promptBadChanSpec;
elseif exist('badchan_file', 'var') && isfile(badchan_file)
    disp(['Loading bad channels from file: ' badchan_file]);
    T = readtable(badchan_file, 'FileType', 'text', 'Delimiter', '\t');
    if ismember('badChanIndex', T.Properties.VariableNames)
        bad_chan_spec = double(T.badChanIndex)';
    else
        error('Column "badChanIndex" not found in badchan_file.');
    end
else
    warning('No interactive input or bad channel file provided — skipping Step 3 manual bad channel rejection.');
    bad_chan_spec = [];
end

% Update globalVar.badChan
globalVar.badChan = unique([bad_chan_spec globalVar.badChan]);
globalVar.badChan_specs.power_spectrum = bad_chan_spec;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      STEP 4: based on HFOs, Su's code 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_paChan LBCN_PREPROC
[pathological_chan_id,pathological_event] = find_paChan(data_all,globalVar.channame,globalVar.iEEG_rate, 1.5);
% pathological_event are in bipolar montage

% Inspect all bad channels
% blue channels - detected in steps 1,2,3
% red channels - detected in step 4
% green channels - detected in both steps (1,2,3) and 4
figureDim = [0 0 .5 .5];
f3 = figure('units', 'normalized', 'outerposition', figureDim);
all_bad_chans = unique([globalVar.badChan(:); pathological_chan_id(:)]);

for ii = all_bad_chans'
    hold on
% for ii = unique([globalVar.badChan pathological_chan_id])
%    hold on
    if ~isempty(find(pathological_chan_id == ii)) && ~isempty(find(globalVar.badChan == ii))
        color_plot = [0 1 0];
    elseif ~isempty(find(pathological_chan_id == ii)) && isempty(find(globalVar.badChan == ii))
        color_plot = [1 0 0];
    elseif isempty(find(pathological_chan_id == ii)) && ~isempty(find(globalVar.badChan == ii))
        color_plot = [0 0 1];
    end
    plot(zscore(data_all(:,ii)) + double(ii), 'Color', color_plot);
end
hold off
title('All bad electrodes')
xlabel('Time (s)')
ylabel('Electrode number')
figure(f3)
saveas(f3, fullfile(diagnostics_dir, sub, sprintf('%s_%s_task-%s_%s_badchannel.png', sub, ses, project_name, run)));

% Update globalVar.badChan
globalVar.badChan = unique([pathological_chan_id' globalVar.badChan]);
globalVar.pathological_event_bipolar_montage = pathological_event;
globalVar.badChan_specs.epileptic_hfo_spike = pathological_chan_id;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Eyeball the rereferenced data after removing the bad channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should be interactive - ask Su to help creating a gui.
%     data_all(:,globalVar.badChan)=[];
%     data_down_car = car(data);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Plot CAR data for eyeballing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iii = 1:size(data_down_car,2)
%         hold on
%         plot(zscore(data_down_car(1:round(globalVar.iEEG_rate*20),iii))+iii);
%     end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Re-referencing data to the common average reference CAR - and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear data
data_all = data_all';
data_good = data_all;
data_good(globalVar.badChan,:) = [];
% Demean good electrodes
for ci = 1:size(data_good,1)
    data_good(ci,:) = data_good(ci,:) - mean(data_good(ci,:));
end
CAR = mean(data_good,1); % common average reference
% Subtract CAR and save
for cii = 1:size(data_all,1)
    
    data.wave = -single(data_all(cii,:) - CAR); % subtract CAR
    % data.wave = -single(data_all(cii,:)); % subtract CAR
    %%% WHY DO I HAVE TO SWITCH SIGN TO MATCH PREVIOUS? (with separate FiltData)!!! %%%
    data.fsample = globalVar.iEEG_rate;
    
    if ~isfolder(globalVar.CARData)
        mkdir(globalVar.CARData)
    end
    write_fn = fullfile(globalVar.CARData, ...
    sprintf('%s_%s_task-%s_%s_electrode-%03i_ieeg-CAR.mat', sub, ses, project_name, run, cii));
    disp(['Writing: ' write_fn])
    save(write_fn,'data');
    clear data
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Save globalVar (For the last time; don't re-write after this point!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(global_fpath,'globalVar');
close all
end




