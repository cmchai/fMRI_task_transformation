% Script to do smoothing after fmriprep preprocessing
% made by Mengqiao Chai, 7 November 2023

%% in prepraration for Univariate first level GLM

clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12

% define all the participants that you want to run
subjs = [3]; % subject IDs to run, tussen 2 tot 49
smooth.fwhm = [8 8 8]; % define the FWHM smoothing kernel in mm
preproc_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/preprocess/results_fmriprep/new_results/prep_23.1.0';      

%% Smooth the images for each run/session/participant

smooth_prefix = ['smth' num2str(smooth.fwhm(1)) '_'];
n_subj = length(subjs);
nses = 2; % 2 scanning sessions
nruns = 4; % 4 runs per session

for subj = subjs % loop over participants

    % get bids_sub: 
    bids_sub = ['sub-00', num2str(subj)];
    
    % DEFINE PATHS
    sub_dir =  fullfile(preproc_rootdir, bids_sub); 
    
    disp('============================================='); 
    disp(['Smoothing data of ' bids_sub]);

    for ses = 1:nses % loop over scanning sessions
        bids_ses = ['ses-mri0', num2str(ses)];
        ses_dir = fullfile(sub_dir, bids_ses);
        func_dir = fullfile(ses_dir, 'func'); % functional data folder
        
        cd(func_dir) % move into the func directory

        for run = 1:nruns % loop over runs

            preproc_img_4d = [bids_sub '_' bids_ses '_task-transform_acq-ep2d_dir-COL_run-' num2str(run) '_echo-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'];
            preproc_img_4d_zipped = [preproc_img_4d '.gz'];

            %% only run the smoothing if the file does not already exists:
            if ~isfile([smooth_prefix preproc_img_4d])
                
                %% Check if files need to be unzippsed first: 
                if ~isfile(preproc_img_4d)
                    preproc_files = gunzip(preproc_img_4d_zipped); % this command will create a unzipped file in the current directory
                else
                    preproc_files = preproc_img_4d;
                end
    
                %% RUN THE SMOOTHING          
                clear matlabbatch;
                matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('expand', fullfile(preproc_files)));
                matlabbatch{1}.spm.spatial.smooth.fwhm = smooth.fwhm;
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix = smooth_prefix;
                spm_jobman('run', matlabbatch);
                clear matlabbatch;

                disp(['====== Finish smoothing data of ', bids_sub, bids_ses, "run", num2str(run)]);
            end 
        end
    end
end