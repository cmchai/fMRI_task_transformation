% Script to form contrasts for each condition of interest(the beta values) from the result of GLM-02 or GLM-02B 
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 12th, Dec, 2023
% 4 regressors of interest, including RG-long-c1, RG-long-c2, TF-long-c1, TF-long-c2
% The contrast image for each condition will naturally average the beta across relevant runs

clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

% want to check some contrast or not? the corrreponding code placed at the end of script
check_con = false;

% Define all the participants that you want to run
subjs_all = [2:5,7:11,13:17,20:44,46:49]; % should be 43 in total

% subjects that will be run
subjs = [4:5,7:11,13:17,20:44,46:49];

% define which sub-GLM model to run
model_B = true;

% the root directory of first level models(FLM):
if model_B == true
    FLM_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-02B/results';
else
    FLM_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-02/results';
end

%% Defining all the contrasts of interest

% the condition of interests
conds_interest = ["RG-long-c1", "RG-long-c2", "TF-long-c1", "TF-long-c2"];
sum_contrasts = length(conds_interest);

% create a list of contrast names. These names are not necessarily required to match up to regressor names in .hdr files.
cnames = {}; 

% Specify T-contrasts: These two cell arrays must have the same length, and the strings must equal to or being the substring of the regressors' names in SPM.mat
connames_pos = {};        

% loop over each condition of interest and define the contrast name and the corresponding label for the later search among predictors

for i = 1:length(conds_interest)
    cond_str = conds_interest(i);
    cond_char = convertStringsToChars(cond_str);

    cnames{i} = cond_char;
    connames_pos{i} = cond_str;
end

% setting the contrast weights
conweight = 1;
negweight = -1;

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
    clear SPM % clear the SPM file from the previous participant
    
    subj_dir = fullfile(FLM_root_dir, bids_subj); 
    load(fullfile(subj_dir, 'SPM.mat'));
    
    % creating the contrast matrix
    C = zeros(sum_contrasts, length(SPM.xX.name)); % the contrast matrix (number of contrasts * length of predictors
    regressors_spm_cellarray = SPM.xX.name; % which is a cell array
    regressors_spm = string(SPM.xX.name);   % which is a string array
    
    %% loop over each contrast and define them by specifying 0, -1, 1

    for i = 1:sum_contrasts % loop through all the conditions
        pos_regressor = connames_pos{i};   
        c_aux_pos = conweight * contains(regressors_spm, pos_regressor);
        %  unit testing: regressors_spm(contains(regressors_spm, pos_regressors)) --> result should be all the relevant positive regressors

        C(i,:) = c_aux_pos;
    end

    % figure; imagesc(C); colorbar
    
    for k = 1:sum_contrasts
        contrasts(k)  = spm_FcUtil('Set',cnames{k},'T','c',C(k,:)',SPM.xX.xKXs);
    end % end contrast loop
    
    % append these new contrasts(the betas for each condition of interest) to the already existing SPM.mat file: 
    SPM.xCon = [SPM.xCon, contrasts];
    
    % make contrasts
    spm_contrasts(SPM);
    
end


%% check some contrasts by hand to see if the contrast coding makes sense

if check_con == true
    if exist("SPM", "var")
        clear SPM
    end

    subj_to_test = 14;
    contrast_to_test = 4;
    
    bids_subj_to_test = ['sub-00', num2str(subj_to_test)];
    subj_to_test_dir = fullfile(FLM_root_dir, bids_subj_to_test);
    load(fullfile(subj_to_test_dir, 'SPM.mat'));
    predictors = (SPM.xX.name)';
    
    contrast_entries = num2cell(SPM.xCon(contrast_to_test).c);
    predictors(:,2) = contrast_entries;
end

