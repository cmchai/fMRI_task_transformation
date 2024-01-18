% Script to run GLM to form the contrast between left and right button press
% images are preprocessed using bidspm, the 5 dummy scans were ALREADY DELETED from the data
% condition 1: task cue onset
% condition 2: stimulus onset
% condition 3: left button press
% condition 4: right button press
% Note: after using bidspm, sub-001 and sub-002 already get rid of dummy
% scans; from sub-003 onwards, dummy scans are NOT deleted
clear all; close all

subjs = 26;
n_ses = 2;
n_run = 4; % number of run per session
n_dummy = 5; % number of dummy scans
smth = 8; % in mm
smth_prefix = ['smth', num2str(smth)];

n_cond = 4; % number of conditions(or regressers of interest) in the GLM
names_conds = {'cueOnset', 'stimOnset', 'keypress-1', 'keypress-2'};

ses_exclu = 2;
run_exclu = 1;

% preproc_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/bids/derivatives/bidspm-preproc/';
preproc_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/preprocess/results_fmriprep/new_results/prep_23.1.0';
event_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/sanity_check/events';
output_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/sanity_check/results03'; % result02: with nuisance regressor, result03: with nuisance regressor and exlude problematic run
batchjob_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts/GLM/sanity_check';
conf_regs_rootdir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-01/events'; % which is in the GLM-01 event folder


bids_func = 'func';

for subj = subjs
    bids_sub = ['sub-00', num2str(subj)];
    % disp(subID)
    
    % the folder of preprocessed images of this subject
    preproc_subj_dir = fullfile(preproc_rootdir, bids_sub);

    % the folder of the event timing files in txt
    event_subj_dir = fullfile(event_rootdir, bids_sub);

    % the FLM output directory of this participant
    output_subj_dir = fullfile(output_rootdir, bids_sub);

    % the folder of nuisance regressor files
    conf_regs_subj_dir = fullfile(conf_regs_rootdir, bids_sub);

    if ~isfolder(output_subj_dir)
        mkdir(output_subj_dir)
    end

    % prepare spm jobman
    spm_jobman('initcfg');

    % record timing
    tic

    clear matlabbatch; % to be sure to start with a fresh job description    
    run('spm_batchjob_sanity_partial.m');     

    sess_spm_real = 0;
    for ses = 1:n_ses
        bids_ses = ['ses-mri0', num2str(ses)];
        
        for run_nr = 1:n_run
            if ses == ses_exclu && run_nr == run_exclu
                continue
            end

            sess_spm_real = sess_spm_real + 1;
            sess_spm = run_nr + 4 * (ses - 1);
            disp(sess_spm_real)

            % go the functional bids folder
            preproc_func_dir = fullfile(preproc_subj_dir, bids_ses, bids_func);

            % selecting the images
            smth_img_4d = [preproc_func_dir '/' smth_prefix '_' bids_sub '_' bids_ses '_task-transform_acq-ep2d_dir-COL_run-' num2str(run_nr) '_echo-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'];
            % disp(img_4d)
            smth_imgs_3d = cellstr(spm_select('expand', fullfile(smth_img_4d)));
            smth_imgs_3d = cellstr(smth_imgs_3d(n_dummy + 1 :end,:));        % IMPORTANT: get rid of dummy scans, otherwise the timing won't be correct !!!!!
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm_real).scans = smth_imgs_3d;

            % input the event timing for each condition
            for cond = 1:n_cond
                event_file_cond = [event_subj_dir '/' bids_sub '_run_' num2str(sess_spm) '_' names_conds{cond} '.txt'];
                event_table_cond = readtable(event_file_cond, 'HeaderLines', 1);
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm_real).cond(cond).onset = table2array(event_table_cond(:,end));
            end

            % specify all the nuisance regressors
            conf_regs_name = [bids_sub '_' 'run-' num2str(sess_spm) '_conf_regs.txt'];
            conf_regs_txt = cellstr(fullfile(conf_regs_subj_dir, conf_regs_name));
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess_spm_real).multi_reg = cellstr(conf_regs_txt);
        end
    end

    %********** step 2: model estimation ***********  
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0; % currently we do not save the residuals (saved in spm.mat file and used from there later on) 
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %- run SPM job %-----------------
    fprintf('Running first level model ... \n');
    spm_jobman('serial', matlabbatch);
    fprintf('First level model - done\n');
    toc
    disp('=============================================')
    disp(' ')
end





