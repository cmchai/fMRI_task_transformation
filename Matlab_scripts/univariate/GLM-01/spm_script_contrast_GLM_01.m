% Script to run contrasts from the result of GLM-01
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 14th, Nov, 2023
% 4 regressors of interest, including RG-long-qc, RG-short-qc, TF-long-qc, TF-short-qc
% Main contrasts:
% (1) (RG-long-qc - RG-short-qc)-(TF-long-qc - TF-short-qc) and vice versa
% (2) RG-long-qc - TF-long-qc and vice versa

clear all; % close all;

% Add spm12 + more general script-directory and load spm:  
% addpath(genpath('S:/pp02_labholroyd/Share/spm12')); % to add spm12 
% addpath(genpath('S:/pp02_labholroyd/Iris Ikink/replicating_PNASpaper/scripts/')); % to add other scripts 

spm fmri
spm('defaults','FMRI')
global defaults
global UFp; UFp = 0.001;

% the root directory of first level models(FLM): 
FLM_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-01/results';

%% Defining all the contrasts of interest

% create a list of contrast names. These names are not necessarily required to match up to regressor names in .hdr files.
cnames = {}; 

% Specify T-contrasts: These two cell arrays must have the same length, and the strings must equal to or being the substring of the regressors' names in SPM.mat
connames_pos = {};        
connames_neg = {}; 

% the first contrast: Interaction --> (RG-long-qc - RG-short-qc)-(TF-long-qc - TF-short-qc)
cnames{1} = 'interaction-RG>TF';
connames_pos{1} = ["RG-long-qc", "TF-short-qc"];
connames_neg{1} = ["RG-short-qc", "TF-long-qc"];

% the second contrast: (RG-long-qc - TF-long-qc)
cnames{2} = 'long-RG>TF';
connames_pos{2} = ["RG-long-qc"];
connames_neg{2} = ["TF-long-qc"];

sum_contrast_pairs = size(connames_pos, 2);      
sum_contrasts = sum_contrast_pairs*2;       % multiple by 2 since we always will look at the opposite contrast

% the names of opposite contrasts
cnames{1+sum_contrast_pairs} = 'interaction-RG<TF';
cnames{2+sum_contrast_pairs} = 'long-RG<TF';

% setting the contrast weights
conweight = 1;
negweight = -1;

%% Running these contrasts for all participants
% Subjects to run
subjs = [3];

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
        % unit testing: sum(c_aux) --> should be zero

        C(i,:) = c_aux;
        C(i + sum_contrast_pairs, :) = (-1) * c_aux;
    end

    % add defined contrasts to SPM.mat file: 
    SPM.xCon = [] ;% this clears out the SPM.xCon in case some exist. Comment out to add on to existing contrasts
    
    for k = 1:length(cnames)
        contrasts(k)  = spm_FcUtil('Set',cnames{k},'T','c',C(k,:)',SPM.xX.xKXs);
    end % end contrast loop
    
    SPM.xCon = [contrasts];
    
    %make contrasts
    spm_contrasts(SPM);
    
end