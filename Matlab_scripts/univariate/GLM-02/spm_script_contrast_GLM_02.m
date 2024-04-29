% Script to run contrasts from the result of GLM-02 or GLM-02B
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 6th, Dec, 2023
% 4 regressors of interest, including RG-long-c1, RG-long-c2, TF-long-c1, TF-long-c2
% Main contrasts:
% (1) the interaction effect: (RG-long-c2 - RG-long-c1)-(TF-long-c2 - TF-long-c1) and viceversa(4)
% (2) the long CTI window contrast: RG-long-c2 - TF-long-c2 and vice versa (5)
% (3) the main effect of block type: (RG-long-c1 + RG-long-c2 + RG-short-c1)-(TF-long-c1 + TF-long-c2 + TF-short-c1) and vice versa(6)

clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

% want to check some contrast or not? the corrreponding code placed at the end of script
check_con = true;

% Define all the participants that you want to run
subjs_all = [2:5,7:11,13:17,20:44,46:49]; % should be 43 in total

% subjects that will be run
subjs = [7:11,13:17,20:44,46:49];

% define which sub-GLM model to run
model_B = false;

% the root directory of first level models(FLM):
if model_B == true
    FLM_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-02B/results';
else
    FLM_root_dir = '/Volumes/extdrive/Task_Transform_GLM/GLM-02/results';  % '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-02/results';
end

%% Defining all the contrasts of interest

% create a list of contrast names. These names are not necessarily required to match up to regressor names in .hdr files.
cnames = {}; 

% Specify T-contrasts: These two cell arrays must have the same length, and the strings must equal to or being the substring of the regressors' names in SPM.mat
connames_pos = {};        
connames_neg = {}; 

% the first contrast: Interaction --> (RG-long-c2 - RG-long-c1)-(TF-long-c2 - TF-long-c1)
cnames{1} = 'interaction-RG>TF';
connames_pos{1} = ["RG-long-c2", "TF-long-c1"];
connames_neg{1} = ["RG-long-c1", "TF-long-c2"];

% the second contrast: (RG-long-c2 - TF-long-c2)
cnames{2} = 'long-RG>TF';
connames_pos{2} = ["RG-long-c2"];
connames_neg{2} = ["TF-long-c2"];

% the third contrast ---> (RG-long-c1 + RG-long-c2 + RG-short-c1)-(TF-long-c1 + TF-long-c2 + TF-short-c1)
cnames{3} = 'RG>TF';
connames_pos{3} = ["RG-long-c1", "RG-long-c2", "RG-short-c1"];
connames_neg{3} = ["TF-long-c1", "TF-long-c2", "TF-short-c1"];

sum_contrast_pairs = size(connames_pos, 2);  % the overall number of contrast   
sum_contrasts = sum_contrast_pairs*2;        % multiple by 2 since we always will look at the opposite contrast

% the names of opposite contrasts
cnames{1+sum_contrast_pairs} = 'interaction-RG<TF';
cnames{2+sum_contrast_pairs} = 'long-RG<TF';
cnames{3+sum_contrast_pairs} = 'TF>RG';

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

    for i = 1:sum_contrast_pairs % only loop through half the contrast matrix since the other half is the opposite of the first half
        pos_regressors = connames_pos{i};
        neg_regressors = connames_neg{i};
        
        c_aux_pos = conweight * contains(regressors_spm, pos_regressors);
        %  unit testing: regressors_spm(contains(regressors_spm, pos_regressors)) --> result should be all the relevant positive regressors
        c_aux_neg = negweight * contains(regressors_spm, neg_regressors);
        %  unit testing: regressors_spm(contains(regressors_spm, neg_regressors)) --> result should be all the relevant negative regressors
        c_aux = c_aux_pos + c_aux_neg;
        %  unit testing: sum(c_aux) --> should be zero

        C(i,:) = c_aux;
        C(i + sum_contrast_pairs, :) = (-1) * c_aux;
    end

    % add defined contrasts to SPM.mat file: 
    SPM.xCon = [] ;% this clears out the SPM.xCon in case some exist. Comment out to add on to existing contrasts
    
    for k = 1:length(cnames)
        contrasts(k)  = spm_FcUtil('Set',cnames{k},'T','c',C(k,:)',SPM.xX.xKXs);
    end % end contrast loop
    
    SPM.xCon = [contrasts];
    
    % make contrasts
    spm_contrasts(SPM);
    
end


%% check some contrasts by hand to see if the contrast coding makes sense

if check_con == true
    if exist("SPM", "var")
        clear SPM
    end

    subj_to_test = 4;
    contrast_to_test = 5;
    
    bids_subj_to_test = ['sub-00', num2str(subj_to_test)];
    subj_to_test_dir = fullfile(FLM_root_dir, bids_subj_to_test);
    load(fullfile(subj_to_test_dir, 'SPM.mat'));
    predictors = (SPM.xX.name)';
    
    contrast_entries = num2cell(SPM.xCon(contrast_to_test).c);
    predictors(:,2) = contrast_entries;
end

