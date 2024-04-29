% Script to add stimulus onset regressor(also q2 onset regressors in transform trials) into the nuisance regressor txt file in preparation for the following GLM
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 22nd, April, 2024

%% General settings
clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

% define de correct path of all the folders
event_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/FIR-02/events';
FLM_rootdir = '/Volumes/extdrive/Task_Transform_GLM/GLM-02/results'; % should be the GLM-02 folder
conf_regs_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-01/events'; % which is in the GLM-01 event folder

% Define all the participants that you want to run
subjs_all = [2:5,7:11,13:17,20:44,46:49]; % should be 43 in total

subjs = subjs_all; % subjects that will be run
n_run = 8;         % number of runs in total !!!

% the simulus onset suffix in the SPM predictors' names
stim_preds_suffix = "-s*";
q2_pred_suffix = "-q2*";

%% run the FLM for all the subjects
for subj = subjs
    bids_sub = ['sub-00', num2str(subj)];

    % make subject level event folder if the subject specific folder does not exist yet
    if ~isfolder(fullfile(event_rootdir, bids_sub))
        mkdir(fullfile(event_rootdir, bids_sub))
    end

    % load SPM mat from the GLM-02 model
    clear SPM % clear the SPM file from the previous participant    
    FLM_subj_dir = fullfile(FLM_rootdir, bids_sub); 
    load(fullfile(FLM_subj_dir, 'SPM.mat'));

    nscans = SPM.nscan;
    stim_preds_bol = contains(SPM.xX.name, [stim_preds_suffix, q2_pred_suffix]);
    stim_preds_names = SPM.xX.name(stim_preds_bol);
    stim_preds_runs = SPM.xX.X(:,stim_preds_bol);

    % the folder of nuisance regressor files
    conf_regs_subj_dir = fullfile(conf_regs_rootdir, bids_sub);
       
    row_accum = 0;
    for run = 1:n_run
        
        % number of scan of the current run
        n_scan = nscans(run);

        % retrieve the stimulus onset preds for the current run
        % col_start = (run - 1) * 3 + 1;
        % col_end = col_start + 2;  % 3 columns corresponding to the stimulus onset of short, long, and error trials
        % pred_run_name_pattern = sprintf('Sn(%d)*-s*',run);
        % pred_run_name_bol = contains(a,pred_run_name_pattern);
        pred_run_name_pattern = sprintf('Sn(%d)',run);
        pred_col_bol = contains(stim_preds_names, pred_run_name_pattern);

        row_start = row_accum + 1;
        row_end = row_accum + n_scan;
        row_accum = sum(nscans(1:run)); % update the accumulated scans up until the current run
        stim_preds_run = stim_preds_runs(row_start:row_end,pred_col_bol); % all the stimulus onset predictors of the current run

        % specify all the nuisance regressors
        conf_regs_name = [bids_sub '_' 'run-' num2str(run) '_conf_regs.txt'];
        conf_regs_mat = readmatrix(fullfile(conf_regs_subj_dir, conf_regs_name));

        % create and save the new txt file combining the nuisance regressor and the stimulus onset regressors
        conf_regs_new = [conf_regs_mat, stim_preds_run];
        conf_regs_new_name = [bids_sub '_' 'run-' num2str(run) '_conf_regs_new.txt'];
        writematrix(conf_regs_new, fullfile(event_rootdir, bids_sub, conf_regs_new_name),'Delimiter','tab');
    end
end