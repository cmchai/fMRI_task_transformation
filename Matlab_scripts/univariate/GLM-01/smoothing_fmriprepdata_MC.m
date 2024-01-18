% Script to do smoothing after fmriprep preprocessing
% made by Mengqiao Chai, 7 November 2023

%% in prepraration for Univariate first level GLM

clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

% define all the participants that you want to run
subjs_all = [2:5,7:11,13:17,20:44,46:49]; % should be 43 in total

subjs = [49]; % subject IDs to run, tussen 2 tot 49
smooth.fwhm = [8 8 8]; % define the FWHM smoothing kernel in mm
preproc_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/preprocess/results_fmriprep/new_results/prep_23.1.0';      

% to choose if we want to unzip anatomical image
unzip_anat = false;

% to choose if we want to delete the functional images in T1 NATIVE space or not
del_img_t1 = true;

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
    
    %--- Unzip the anatomical image if needed ---%
    if unzip_anat == true
        if isfolder(fullfile(sub_dir,'anat')) % if 2 anatomical scans were acquired
            anat_dir = fullfile(sub_dir,'anat');
            anat_img_unzipped_name = [bids_sub '_acq-GIfMIT1MPRAGE_run-1_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii'];
        else
            anat_dir = fullfile(sub_dir, 'ses-mri01', 'anat'); % if only 1 anatomical scan was acquired
            anat_img_unzipped_name = [bids_sub '_ses-mri01_acq-GIfMIT1MPRAGE_run-1_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii'];
        end
        anat_img_zipped_name = [anat_img_unzipped_name '.gz'];
        if ~isfile(fullfile(anat_dir, anat_img_unzipped_name))
            anat_img_unzipped_file = gunzip(fullfile(anat_dir, anat_img_zipped_name));
            disp(['Finished unzipping the anatomical image of ' bids_sub]);
        end
    end

    disp('============================================='); 
    disp(['Starting smoothing data of ' bids_sub]);

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

            %% Delete the zipped file if they still exist in the folder since they occupy too much space
            if isfile(preproc_img_4d_zipped)
                delete(preproc_img_4d_zipped)
                disp(['====== Finish del zipped func images of', bids_sub, bids_ses, "run", num2str(run)])
            end

            %% Delete the functional images in the T1w space(which is default output in fmriprep) just to save space on the pc
            if del_img_t1 == true
                preproc_img_4d_t1_name = [bids_sub '_' bids_ses '_task-transform_acq-ep2d_dir-COL_run-' num2str(run) '_echo-1_space-T1w_desc-preproc_bold.nii.gz'];
                if isfile(preproc_img_4d_t1_name)
                    delete(preproc_img_4d_t1_name)
                    disp(['======= finish del func T1w', bids_sub, bids_ses, "run", num2str(run)]);
                end
            end
        end
    end
end