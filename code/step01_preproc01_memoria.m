%% Define the top directory
function step01_preproc01(SLURM_ID, main_dir, bidsdata_dir, lbcn_preproc, project_name)
    i = str2num(SLURM_ID);
    cd(main_dir);
    git_root = fileparts(main_dir);
    addpath(genpath(lbcn_preproc)); % preprocessing pipeline
    addpath(genpath(fullfile(lbcn_preproc, 'preproc_WH'))); % additional code from Weichen Huang
    addpath(genpath(fullfile(git_root, 'code', 'code_from_weichen')));
    addpath(genpath(fullfile(git_root, 'code')));

    %% Parameters
    analysis_dir = fullfile(main_dir, 'analysis');
    ieegsource_dir = fullfile(main_dir, 'analysis', 'ieeg-source');% sub/ses/run
    diagnostics_dir = fullfile(main_dir, 'analysis', 'diagnostics', 'badchannel');
    center = 'Stanford';
    data_format = 'edf';
    
    metadata_fname = fullfile(main_dir, 'participants_source2BIDSmapping.tsv');
    T = readtable(metadata_fname, 'Delimiter', '\t', 'FileType', 'text');

    sub = T.BIDS_sub{i}; 
    ses = T.BIDS_ses{i};
    run = T.BIDS_run{i};
    fprintf('---- Processing %s %s %s ----', sub, ses, run)

    % Define empty channels
    refChan = [];
    epiChan = [];
    emptyChan = [];

    %% 01 create globalVar that has metadata and will house bad channel info
    BIDS_step01_createglobalVar(sub, ses, run, project_name, center, bidsdata_dir, ieegsource_dir);
    %% 1-1 extract electrode info from EDF file
    if strcmp(data_format, 'edf')
        BIDS_SaveDataNihonKohden(bidsdata_dir, ieegsource_dir, sub, project_name, ses, run, center, refChan, epiChan, emptyChan)
    else
        error('Data format has to be either edf or TDT format');
    end
    %% 02, 03 organize events file and cross-compare with EDF data
    BIDS_step02_process_events(bidsdata_dir, sub, ses, run ); %step02_processBeh(beh_dir, sbj_name, block_name, block_num, project_name);
    BIDS_step03_organizetrialinfoHJ(sub, project_name, ses, run, main_dir, bidsdata_dir, ieegsource_dir);
    %% 04 Reject bad channels
    use_interactive = 'false'; badchan_file = 'emptyfile';
    BIDS_step04_badchannelreject_generic(project_name, sub, ses, run, main_dir, ieegsource_dir, bidsdata_dir, diagnostics_dir,[],use_interactive, badchan_file);

end
