% Script to Concatenate the 3d predictor nifti images into 4d for the multivariate analysis
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 15th, Feb, 2024
% images to be concatenated can be either beta images or t contrast images
%
% Some pre-defined parameters:
%  prefix_preds can be  {'RG-short-task\d-c1', '-short-task\d-c1', 'RG-\w*-task\d-c1', 'TF-short-task\d-c1',  'TF-\w*-task\d-c1'};
% 'RG-short-task\d-c1' and 'TF-short-task\d-c1': the CTI onset of short trials in either RG or TF blocks(should include 4 runs with each run containing 1 sample per task)
% 'RG-\w*-task\d-c1' and 'TF-\w*-task\d-c1' :    the CTI onset of both short and long trials in either RG or TF blocks (should include 4 runs with each run containing 2 samples per task)
% '-short-task\d-c1':  the CTI onset of short trials in BOTH RG or TF blocks (8 runs)
%
% prefix_4dfile can be {'RG-short','ALL-short','RG-all-c1', 'TF-all-c1};
% prefix_img can be {'beta', 'spmT'};

clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts'));  % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

%% general settings

% Define which model results to look at
FLM = 'GLM-02M'; % 'GLM-02M' for unsmoothed and 'GLM-02M-B' for smth8
FLM_root_dir = fullfile('/Volumes/extdrive/Task_Transform_GLM', FLM, 'results');

% Define all the participants that you want to run
subjs_all = [2:5,7:11,13:17,20:44,46:49];  % should be 43 in total

% Subjects that will be run
subjs_run = [0];              % subjects that are already run
subjs = setdiff(subjs_all, subjs_run);    % subjects that will be run

% Predictors of interests for concatenation
prefix_preds = 'RG-\w*-task\d-c1';     % or '-short-task\d-c1' if include both RG and TF blocks, or 'RG-\w*-task\d-c1'
prefix_4dfile = 'RG-all-c1';           % or 'ALL-short' if include both RG and TF blocks, or 'RG-all-c1'
prefix_smth = 'smthN';                 % if the images were smoothed before GLM

% which measure to concatenate: beta or T-contrast
prefix_img = 'spmT';  % or 'beta'

%% Concatenate files for each participant

for subj = subjs

    % get bids_sub 
    bids_subj = ['sub-00', num2str(subj)];
    
    % clear all the variables and SPM file from the previous participant first
    clear SPM 
    clearvars images tasks runs

    % prepare spm jobman
    spm_jobman('initcfg');
    clear matlabbatch;  % to be sure to start with a fresh job description for each participant      
    
    % load SPM mat
    subj_dir = fullfile(FLM_root_dir, bids_subj); 
    load(fullfile(subj_dir, 'SPM.mat'));

    % create the table of predictors of interest and their corresponding number in the SPM results
    regressors_spm = string(SPM.xX.name);            % the names of all the predictors, which is a string array
    regressors_spm_num = [1:length(regressors_spm)]; % the corresponding number in the spm output
    idx = ~cellfun(@isempty, regexp(regressors_spm, prefix_preds));   % a Boolean vector that give 1 to the regressors of interest;

    preds = regressors_spm(idx)';        % the names of predictors of interest
    preds_nrs = regressors_spm_num(idx); % the corresponding number in the spm output

    % concatenating images
    images = cell(length(preds_nrs),1); % initiate the cell array
    tasks = zeros(size(images));        % initiate the task labels
    runs = zeros(size(images));         % initiate the run labels

    for i = 1:length(preds_nrs)         % loop over each predictor of interest
        pred = preds(i);
        
        task = extractBetween(pred, "task", "-"); % the task id 1-9
        tasks(i) = str2double(task);
        
        run = extractBetween(pred, 4, 4); % the run id 1-8
        runs(i) = str2double(run);
        
        pred_nr = preds_nrs(i);
        if pred_nr < 10
            prefix_pred_nr = ['000', num2str(pred_nr)];
        elseif pred_nr < 100
            prefix_pred_nr = ['00', num2str(pred_nr)];
        else
            prefix_pred_nr = ['0', num2str(pred_nr)];
        end

        img_name = [prefix_img, '_', prefix_pred_nr, '.nii,1'];
        img_path = fullfile(subj_dir, img_name);
        images{i,1} = img_path;
    end

    % first of all, save tasks and runs in a .txt file for the later decoding or RSA
    label_table = table(runs, tasks);
    writetable(label_table, fullfile(subj_dir, [bids_subj, '_', prefix_smth, '_', prefix_img, '_', prefix_4dfile, '_labels.txt']));

    % start the spm batch job
    matlabbatch{1}.spm.util.cat.vols = images;
    matlabbatch{1}.spm.util.cat.name = [bids_subj, '_', prefix_smth, '_', prefix_img, '_', prefix_4dfile, '_4D.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 0;   % which is the same number of digits as the 3d images
    % matlabbatch{1}.spm.util.cat.RT = NaN;
    spm_jobman('run',matlabbatch);

    % some testing codes to check the saved 4d images, comment out the following section if not debugging
    % Fd_img = niftiread(fullfile(subj_dir, [bids_subj, '_', prefix_img, '_', prefix_4dfile, '_4D.nii']));
    % Fd_img_info = niftiinfo(fullfile(subj_dir, [bids_subj, '_', prefix_img, '_', prefix_4dfile, '_4D.nii']));
    % Fd_img_info.ImageSize(4) == length(preds_nrs) % should be 1 or True!

end
