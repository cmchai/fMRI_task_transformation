% Script to run Second level analysis(SLM_02) based on the result of GLM-02 for the Task Transformation Paradigm
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 6th, Dec, 2023

clear all; close all;
clear spm;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

%% Global setttings:

% subjects to run
subjs_all = [2:5,7:11,13:17,20:44,46:49]; % should be 43 in total
subjs = subjs_all;

% IMPORTANT: Specify which GLM model you want to run
model_B = true; 

% specify all the input/output directories
jobfile_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts/Matlab_GLM/univariate/GLM-02';

if model_B == true
    FLM_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-02B/results';
    output_dir  = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/second_level/univariate/GLM-02B/results';
else
    FLM_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-02/results';
    output_dir  = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/second_level/univariate/GLM-02/results';
end

% Specify welke FLM contrasts that we want to do second level analysis
cons = [7:10];
% -----------------------------------------

spm fmri

% write matrix with all directories of results of FLM for all participants: 
fullSBdirs = cell(length(subjs),1);

for i = 1:length(subjs)
    subj = subjs(i);
    bids_sub = ['sub-00', num2str(subj)];
    fullSBdirs{i} = {fullfile(FLM_rootdir, bids_sub)};    
end

% -----------------------------------------
% set SPM
spm('Defaults','fmri');
spm_jobman('initcfg');

% Specify which jobfile should be run
jobfile = {fullfile(jobfile_dir, 'spm_batchjob_SLM_02.m')};  % Attention: the format should be a cell array

jobs = repmat(jobfile, 1, length(cons)); % one job per contrast
inputs = cell(2, length(cons));

for i = 1:length(cons) % loop over each contrast
    con = cons(i); % the number indicator of the contrast
    
    % undefined batchjob input 1: output_dir
    contrast_folder = ['Contrast' num2str(con)];
    mkdir(fullfile(output_dir, contrast_folder));
    inputs{1, i} = {fullfile(output_dir, contrast_folder)};
    
    % undefined batchjob input 2: scans (or actually cons) of all PPs
    for ii = 1:length(subjs)
        subj = subjs(ii);
        if con < 10
            con_inputs(ii,1) = cellstr(cfg_getfile('FPList', fullSBdirs{ii}, ['con_000' num2str(con) '.nii'])); 
        elseif con >=10
            con_inputs(ii,1) = cellstr(cfg_getfile('FPList', fullSBdirs{ii}, ['con_00' num2str(con) '.nii']));
        end
    end 
    
    inputs{2, i} = con_inputs;
    
end

spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:}); %spm_jobman('serial', jobs, '', inputs{:});

