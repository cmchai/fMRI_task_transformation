#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:53:10 2024

Master Script to run decoding using Searchlight approach

@author: mengqiao
"""

import os
import sys
# import glob
import time

cwd = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts/Python_multivariate'
sys.path.append(cwd)

import funcs
import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 

import numpy as np
import pandas as pd
import sklearn.svm as svm
from sklearn.model_selection import LeaveOneGroupOut

from nilearn import image, plotting
from nilearn.decoding import SearchLight
import nibabel as nib

#%% Define the FLM results(with corresponding suffix) to be decoded(should be in 4D format)

GLM_id = 'GLM-02M-B'  # GLM-02M used unsmoothed, GLM-02M-B used smoothed 8mm data for GLM
FLM_root_dir = os.path.join('/Volumes/extdrive/Task_Transform_GLM', GLM_id,'results')

if GLM_id == 'GLM-02M':
    smth_suffix = 'smthN'
elif GLM_id == 'GLM-02M-B':
    smth_suffix = 'smth8'
  
output_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/searchlight'
output_dir = os.path.join(output_root_dir, smth_suffix)

d4_files_suffix = 'RG-all-c1_4D.nii'

#%% settings regarding decoding and cross-validation

decoders = 'SVM'
need_regroup = True  # specify if reassigning runs into folds are necessary(from 4 runs into 8 folds)

#%% settings regarding the brain and grey matter masks
spm_mask_file = 'mask.nii' # using the spm brain mask as the brain mask in searchlight analysis

preproc_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/preprocess/results_fmriprep/new_results/prep_23.1.0'
GM_prob_suffix = '_acq-GIfMIT1MPRAGE_run-1_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz'

#%% Define all the subjects and the GLM results to run

subjs_all = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50)))
subjs = subjs_all[30:]
# subjs = np.setdiff1d(subjs_all, subjs_run)

#%% Create the numpy array to save all the searchlight results for all participants

sl_acc_all = np.zeros((78, 93, 78, subjs_all.size)) # (78, 93, 78) is the total number of voxels

#%% perform the Searchlight Analysis

for subj_idx, subj in enumerate(subjs, 30): # loop over subjects

    t_start = time.time()
    
    bids_subj = 'sub-00' + str(subj)
    subj_dir = os.path.join(FLM_root_dir, bids_subj)
    
    # load the 4D image
    d4_file = funcs.search_files(subj_dir, bids_subj, d4_files_suffix)[0] # the name of the 4d data
    d4_path = os.path.join(subj_dir, d4_file)
    d4_image = image.load_img(d4_path)
    
    # load the label file, define the labels and groups
    labels_file = d4_file.replace('4D.nii', 'labels.txt')
    labels_path = os.path.join(subj_dir, labels_file)
    labels_df = pd.read_csv(labels_path, sep=",")

    labels_all = labels_df['tasks'].to_numpy().copy()
    groups = labels_df['runs'].to_numpy().copy()
    if need_regroup:
        groups = funcs.from_4runs_to_8groups(labels_all)
        
    # load the brain mask from SPM results
    spm_mask_path = os.path.join(subj_dir, spm_mask_file)
    spm_mask = image.load_img(spm_mask_path)
    
    # create the grey matter mask
    GM_mask_dir = os.path.join(preproc_root_dir, bids_subj, 'anat')
    
    if os.path.exists(GM_mask_dir):
        pass
    else:
        GM_mask_dir = os.path.join(preproc_root_dir, bids_subj, 'ses-mri01' ,'anat')
        
    GM_mask_file = bids_subj + '_GM-mask-0.2-crop.nii'
    GM_mask_path = os.path.join(GM_mask_dir, GM_mask_file)
    
    if os.path.isfile(GM_mask_path): # if the cropped mask file already exist, then load that file 
        GM_mask_crop = image.load_img(GM_mask_path)
    else:  # create and save the cropped mask if the mask does not exist yet
        GM_prob_file = bids_subj + GM_prob_suffix
        GM_prob_path = os.path.join(GM_mask_dir, GM_prob_file)
        
        if os.path.isfile(GM_prob_path):
            pass 
        else:
            GM_prob_file = bids_subj + '_ses-mri01' + GM_prob_suffix
            GM_prob_path = os.path.join(GM_mask_dir, GM_prob_file)
        
        GM_prob_img = image.load_img(GM_prob_path)
        GM_prob_img_resample = image.resample_to_img(GM_prob_img, d4_image, interpolation='continuous')

        GM_mask = image.math_img('(img > 0.2).astype(np.int32)', img=GM_prob_img_resample)
        GM_mask_crop = image.math_img('(np.multiply(img1, img2)).astype(np.int32)', img1=GM_mask, img2=spm_mask)
        nib.nifti1.save(GM_mask_crop, GM_mask_path)
  
    # define decoder and cross-validation scheme
    clf = svm.SVC(kernel='linear')
    logo = LeaveOneGroupOut()
    
    # Run the SearchLight Analysis
    sl = SearchLight(
        mask_img=spm_mask,
        process_mask_img=GM_mask_crop,
        radius=5,  # 5 mm radius(2 voxels)
        estimator=clf,
        n_jobs=3,  # use only 1 core (for your own analyses, you might want to increase this!)
        scoring='accuracy',
        cv=logo,
        verbose=False  # print a progress bar while fitting
    )

    sl.fit(d4_image, labels_all, groups=groups)
    
    # saving the decoding result in a numpy array
    sl_acc = sl.scores_ # the accuray of decoding for each voxel for this participant in np array
    sl_acc_all[:,:,:,subj_idx] = sl_acc
    
    # saving the decoding result as a nifti image for each participant
    sl_acc_img = image.new_img_like(spm_mask, sl_acc)
    sl_acc_file = bids_subj + '_' + smth_suffix + '_spmT_' + d4_files_suffix.replace('4D.nii', 'SL_acc.nii')
    sl_acc_path = os.path.join(output_dir, sl_acc_file)
    nib.nifti1.save(sl_acc_img, sl_acc_path)
    
    t_elapsed = round((time.time() - t_start)/60, 2) # scale is minute
    print(f"finished decoding for the {subj_idx}th subject, which is {bids_subj}, takes {t_elapsed} mins")
    
#%% Saving the results
result_dir = "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/searchlight/smth8"
sl_acc_np_name = 'searchlight_acc_RG-all-c1_smth8_full'
np.save(os.path.join(result_dir, sl_acc_np_name), sl_acc_all)  
 
#%% load the result np array
sl_acc_np_file = 'searchlight_acc_RG-all-c1_smth8.npy'
sl_acc_all = np.load(os.path.join(result_dir, sl_acc_np_file))

#%% create the average decoding image

sl_acc_all = np.load(os.path.join(result_dir, 'searchlight_acc_RG-all-c1_smth8_full.npy'))
# a = sl_acc_all[33,20,34,:]
sl_acc_mean = np.mean(sl_acc_all, axis=3)
model_img = image.load_img(os.path.join(result_dir, subj_decode_acc_files[3]))
sl_acc_mean_img = nib.Nifti1Image(sl_acc_mean, model_img.affine, model_img.header)
nib.nifti1.save(sl_acc_mean_img, os.path.join(result_dir, 'sl_acc_mean.nii'))

#%% smoothing the average decoding image(fwhm = 8)i and save it
from scipy.ndimage import gaussian_filter

fwhm = 8
voxelsize = 2.5

sigma = fwhm / (np.sqrt(8 * np.log(2)) * voxelsize)

# sl_acc_all_smth8 = gaussian_filter(sl_acc_all, sigma=sigma, axes=(0,1,2))
acc_files = glob.glob(os.path.join(result_dir, '*_acc.nii'))
sl_acc_all_smth8_imgs = image.smooth_img(acc_files, fwhm=8)
sl_acc_all_smth8_4d = image.concat_imgs(sl_acc_all_smth8_imgs)
sl_acc_smth8_mean_img = image.math_img('np.mean(img, axis=3)', img=sl_acc_all_smth8_4d)
sl_acc_smth8_mean = sl_acc_smth8_mean_img.get_fdata()
nib.nifti1.save(sl_acc_smth8_mean_img, os.path.join(result_dir, 'sl_acc_smth8_mean.nii'))

sl_con_all_smth8_4d = image.math_img('img - 0.111111', img=sl_acc_all_smth8_4d)
nib.nifti1.save(sl_con_all_smth8_4d, os.path.join(contrast_ouput_dir, 'sl_smth8', 'sl_con_smth8_4d.nii'))

#%% create contrast image for each participant first (decoding accuray - chance level)
chance_level = 0.111111
contrast_ouput_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/searchlight/smth8/second_level'
subj_decode_acc_files = funcs.search_files(result_dir, 'sub-00', 'acc.nii')

for subj_acc_file in subj_decode_acc_files:
    subj_acc_img = image.load_img(os.path.join(result_dir, subj_acc_file))
    # plotting.plot_glass_brain(subj_acc_img, colorbar=True)
    subj_con_img = image.math_img('img - 0.111111', img=subj_acc_img)
    # plotting.plot_glass_brain(subj_con_img, colorbar=True, threshold=0)
    subj_con_file = subj_acc_file.replace('acc', 'con')
    nib.nifti1.save(subj_con_img, os.path.join(contrast_ouput_dir, subj_con_file))

#%% conduct second-level analysis

# generate a list of contrast file paths

import glob

con_files = glob.glob(os.path.join(contrast_ouput_dir, '*_con.nii'))

# create design matrix

design_matrix = pd.DataFrame([1] * len(con_files), columns=["intercept"])

# fit the second level model

from nilearn.glm.second_level import SecondLevelModel

second_level_model = SecondLevelModel(n_jobs=2).fit(
    con_files, design_matrix=design_matrix
)

z_map = second_level_model.compute_contrast(output_type="z_score")

display = plotting.plot_stat_map(z_map, title="Raw z map", threshold=0)

# threshold the image

from nilearn.glm import threshold_stats_img

thresholded_map1, threshold1 = threshold_stats_img(
    z_map,
    alpha=0.05,
    cluster_threshold=0,
    two_sided=False,
)

plotting.plot_stat_map(
    thresholded_map1,
    threshold=threshold1,
    cmap="binary",
    title="Thresholded z map, fpr <.001, clusters > 10 voxels",
)




