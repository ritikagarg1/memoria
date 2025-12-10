% # sub="sub-01"
% # ses="ses-01"
% # run="run-02"
% # PROJECT_DIR="/Users/heejungj/Documents/projects_local/emotionalappraisal"
% # fig_path="${PROJECT_DIR}/analysis/diagnostics/badchannel/${sub}/${sub}_${ses}_task-emotionalappraisal_${run}_powerspectrum.png"
% # globalvar_path="${PROJECT_DIR}/analysis/ieeg-source/${sub}/${sub}_${ses}_task-emotionalappraisal_${run}_global.mat"
% # output_tsv_path="${PROJECT_DIR}/analysis/ieeg-source/${sub}"
% # step01_preproc_part02_manualbadchannel(fig_path, globalvar_path, output_tsv_path)

% Define variables
sub = "sub-259";
ses = "ses-01";
run = "run-01";
project_name = "memoria";
PROJECT_DIR = "/Users/ritikag/Documents/projects_local/memoria";

% Construct paths
fig_path = fullfile(PROJECT_DIR, "analysis", "diagnostics", "badchannel", sub, ...
    sub + "_" + ses + "_task-" + project_name + "_" + run + "_powerspectrum.fig");

globalvar_path = fullfile(PROJECT_DIR, "analysis", "ieeg-source", sub, ...
    sub + "_" + ses + "_task-" + project_name + "_" + run + "_global.mat");

output_tsv_path = fullfile(PROJECT_DIR, "analysis", "ieeg-source", sub, ...
    sub + "_" + ses + "_" + run + "_manualBadChans.tsv");

% Call the function
BIDS_step04_manual_badchannelappend(fig_path, globalvar_path, output_tsv_path);
