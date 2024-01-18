% Script to run FIR-00 for the Task Transformation Paradigm
% images are preprocessed using fMRIprep with SYN as the sdc approach
% For each fucntional run, there are 5 dummy scans that need to exclude before GLM
% Images should be smoothed before the GLM

%% General settings
clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));


% define de correct path of all the folders
event_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/FIR-00/events';
output_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/FIR-00/results';

preproc_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/preprocess/results_fmriprep/new_results/prep_23.1.0';      
conf_regs_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-01/events'; % which is in the GLM-01 event folder

% the directory for the functional images in bids format
bids_func = 'func';

% Define all the participants that you want to run
subjs_all = [2:5,7:11,13:17,20:44,46:49]; % should be 43 in total

subjs = [2];     % subjects that will be run
n_ses = 2;       % two scanning sessions, which is different from "sess" in SPM
n_run = 4;       % number of run per session
n_dummy = 5;     % number of dummy scans
smth = 8;        % in mm
smth_prefix = ['smth', num2str(smth)];

%% run the FLM for all the subjects
for subj = subjs
    bids_sub = ['sub-00', num2str(subj)];

    % make subject level output folder if the subject specific folder does not exist yet
    if ~isfolder(fullfile(output_rootdir, bids_sub))
        mkdir(fullfile(output_rootdir, bids_sub))
    end

    % prepare spm jobman
    spm_jobman('initcfg');

    % timing
    tic

    clear matlabbatch; % to be sure to start with a fresh job description for each participant    
    % run([dir_batchjob 'GLM_01_job.m']);

    % ******** general SPM batchjob info for the given participants *****
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(fullfile(output_rootdir, bids_sub)); % the FLM result directory
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.78;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 54;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 27;
    
    % the folder of preprocessed images of this subject
    preproc_subj_dir = fullfile(preproc_rootdir, bids_sub);

    % the folder of the event timing files in txt
    event_subj_dir = fullfile(event_rootdir, bids_sub);

    % the folder of nuisance regressor files
    conf_regs_subj_dir = fullfile(conf_regs_rootdir, bids_sub);

    for ses = 1:n_ses
        bids_ses = ['ses-mri0', num2str(ses)];
        
        for run = 1:n_run
            sess_spm = run + 4 * (ses - 1);
            % disp(sess_spm)

            % ******** general SPM batchjob info for the Given run *****
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).multi = {''};
            % matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).hpf = 128;
            
            % enter the func folder
            preproc_func_dir = fullfile(preproc_subj_dir, bids_ses, bids_func);
            
            % select the corresponding image after preprocessing and smoothing
            smth_img_4d = [preproc_func_dir '/' smth_prefix '_' bids_sub '_' bids_ses '_task-transform_acq-ep2d_dir-COL_run-' num2str(run) '_echo-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'];
            % disp(img_4d)
            smth_imgs_3d = cellstr(spm_select('expand', fullfile(smth_img_4d)));
            smth_imgs_3d = cellstr(smth_imgs_3d(n_dummy + 1 :end,:)); % IMPORTANT: get rid of dummy scans, otherwise the timing won't be correct !!!!!
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).scans = smth_imgs_3d;

            % specify the event onsets(and the corresponding durations if applied)
            events_run = dir(fullfile(event_subj_dir, sprintf('*run-%d*condition*.txt',sess_spm))); % which is a structure
            n_cond = size(events_run, 1); % the number of *.txt files corresponds to how many predictors(conditions)
            
            for cond = 1:n_cond
                event_txt = fullfile(event_subj_dir, events_run(cond).name);
                event_table_cond = readtable(event_txt);
                name_cond = char(event_table_cond.condition(1));
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).cond(cond).name = name_cond;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).cond(cond).onset = table2array(event_table_cond(:,end-2));
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).cond(cond).duration = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).cond(cond).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).cond(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).cond(cond).orth = 0;
            end

            % specify all the nuisance regressors
            % conf_regs_name = [bids_sub '_' 'run-' num2str(sess_spm) '_conf_regs.txt'];
            % conf_regs_txt = cellstr(fullfile(conf_regs_subj_dir, conf_regs_name));
            % matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).multi_reg = cellstr(conf_regs_txt);
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm).multi_reg = {''};
        end
    end
    % ********* the rest of batchjob specifications ***********
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = 20;
    matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order = 10;
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.95;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST';

    %********** step 2: model estimation ***********  
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0; % currently we do not save the residuals (saved in spm.mat file and used from there later on) 
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %- run SPM job %-----------------
    fprintf('Running first level model for %s \n',bids_sub);
    spm_jobman('serial', matlabbatch);
    fprintf('First level model - done for %s \n', bids_sub);
    toc
    disp('=============================================')
    disp(' ')
end
