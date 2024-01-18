% Script to run Second level analysis(SLM_01) based on the result of GLM-01 for the Task Transformation Paradigm
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 20th, Nov, 2023
%% 
clear all; close all;
clear spm;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

%% Global setttings:
subjs_all = [2:5,7:11,13:17,20:44,46:49]; % should be 43 in total
subjs = subjs_all;
nrun = 6; % number of contrasts

% IMPORTANT: Specify which GLM model you want to run
model_B = true; 

% Study specific directories
jobfile_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts/Matlab_GLM/univariate/GLM-01';

if model_B == true
    FLM_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-01B/results';
    output_dir  = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/second_level/univariate/GLM-01B/results';
else
    FLM_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-01/results';
    output_dir  = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/second_level/univariate/GLM-01/results';
end
  
% -----------------------------------------

spm fmri

% write matrix with all directories for each subject: 
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
jobfile = {fullfile(jobfile_dir, 'spm_batchjob_SLM_01.m')};  % Attention: the format should be a cell array

jobs = repmat(jobfile, 1, nrun); % one job per contrast
inputs = cell(2, nrun); % cell(numsub, nrun); 

for crun = 1:nrun
    
    % undefined batchjob input 1: output_dir
    contrast_folder = ['Contrast' num2str(crun)];
    mkdir(fullfile(output_dir, contrast_folder));
    inputs{1, crun} = {fullfile(output_dir, contrast_folder)};
    
    % undefined batchjob input 2: scans (or actually cons) of all PPs
    for i = 1:length(subjs)
        subj = subjs(i);
        if crun < 10
            crun_inputs(i,1) = cellstr(cfg_getfile('FPList', fullSBdirs{i}, ['con_000' num2str(crun) '.nii'])); 
        elseif crun >=10
            crun_inputs(i,1) = cellstr(cfg_getfile('FPList', fullSBdirs{i}, ['con_00' num2str(crun) '.nii']));
        end
    end 
    
    inputs{2,crun} = crun_inputs;

end

spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:}); %spm_jobman('serial', jobs, '', inputs{:});
