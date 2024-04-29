% Script to form T contrasts for all the regressors from the result of GLM-02M 
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 26th, Jan, 2024
% Form the T contrast for the regressors, so the number of contrasts should be consistent with the number of Beta image.
% The relevant contrasts will be feed into the multivariate analysis later

clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts'));  % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

% want to check some contrast or not? the corrreponding code placed at the end of script
check_con = false;

% Define all the participants that you want to run
subjs_all = [2:5,7:11,13:17,20:44,46:49];  % should be 43 in total

% subjects that will be run
subjs_run = [];                            % subjects that are already run
subjs = setdiff(subjs_all, subjs_run);     % subjects that will be run

% the root directory of first level models(FLM):
model_B = true;    % model B used smoothed data(8mm)

if model_B == true
   FLM_root_dir = '/Volumes/extdrive/Task_Transform_GLM/GLM-02M-B/results';  
else
   FLM_root_dir = '/Volumes/extdrive/Task_Transform_GLM/GLM-02M/results';
end

% initiate a list of contrast names. 
cnames = {};

%% Running these contrasts for all participants
spm fmri
spm('defaults','FMRI')
global defaults
global UFp; UFp = 0.001;

for subj = subjs
    disp(' ')
    disp('=============================================')
    disp(['Subject: ' num2str(subj)])
    
    % -----------------------------------------
    % DEFINE PATHS
    % get bids_sub: 
    bids_subj = ['sub-00', num2str(subj)];
    
    % -----------------------------------------
    % load SPM mat
    clear SPM  % clear the SPM file from the previous participant    
    clear contrasts % clear contrast structure of the previous participant

    subj_dir = fullfile(FLM_root_dir, bids_subj);
    load(fullfile(subj_dir, 'SPM.mat'));

    % define the names of the contrasts, these names match the regressors' names in the SPM.mat file
    cnames = SPM.xX.name;
    
    % creating the contrast matrix
    C = eye(length(cnames));  % the contrast matrix (length of predictors * length of predictors)
    % figure; imagesc(C); colorbar
    
    for k = 1:length(cnames)
        contrasts(k)  = spm_FcUtil('Set',cnames{k},'T','c',C(k,:)',SPM.xX.xKXs);
    end % end contrast loop
    
    % clear out the contrast list in the SPM.mat first and add these contrasts in 
    SPM.xCon = [];
    SPM.xCon = [contrasts];
    
    % make contrasts
    spm_contrasts(SPM);

end

%% check some contrasts by hand to see if the contrast coding makes sense

if check_con == true
    if exist("SPM", "var")
        clear SPM
    end

    subj_to_test = 3;
    contrast_to_test = 4;
    
    bids_subj_to_test = ['sub-00', num2str(subj_to_test)];
    subj_to_test_dir = fullfile(FLM_root_dir, bids_subj_to_test);
    load(fullfile(subj_to_test_dir, 'SPM.mat'));
    predictors = (SPM.xX.name)';
    
    contrast_entries = num2cell(SPM.xCon(contrast_to_test).c);
    predictors(:,2) = contrast_entries;
end