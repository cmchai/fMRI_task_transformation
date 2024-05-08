#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 14:53:10 2024

Master Script to run decoding using Searchlight approach on HPC

@author: mengqiao
"""

import os
import sys
import time

scripts_path = '/kyukon/scratch/gent/vo/001/gvo00170/vsc43896/task_transform_fmri/searchlight_task_decode/scripts'
sys.path.append(scripts_path)
import funcs_hpc

import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 

import numpy as np
import pandas as pd
import sklearn.svm as svm
from sklearn.model_selection import LeaveOneGroupOut

from nilearn import image
from nilearn.decoding import SearchLight
import nibabel as nib

#%% Define the FLM results(with corresponding suffix) to be decoded(should be in 4D format)

GLM_id = 'GLM-02M'  # GLM-02M used unsmoothed, GLM-02M-B used smoothed 8mm data for GLM
FLM_root_dir = os.path.join('/kyukon/scratch/gent/vo/001/gvo00170/vsc43896/task_transform_fmri/FLM_4d_data', GLM_id)

if GLM_id == 'GLM-02M':
    smth_suffix = 'smthN'
elif GLM_id == 'GLM-02M-B':
    smth_suffix = 'smth8'
  
output_root_dir = '/kyukon/scratch/gent/vo/001/gvo00170/vsc43896/task_transform_fmri/searchlight_task_decode/results'
output_dir = os.path.join(output_root_dir, GLM_id)

d4_files_suffix = 'ALL-short_4D.nii' # either 'ALL-short_4D.nii' or 'RG-all-c1_4D'

#%% settings regarding decoding and cross-validation

decoders = 'SVM'
need_regroup = False  # specify if reassigning runs into folds are necessary(from 4 runs into 8 folds)

#%% settings regarding the brain and grey matter masks
spm_mask_file = 'mask.nii' # using the spm brain mask as the brain mask in searchlight analysis
preproc_root_dir = '/kyukon/scratch/gent/vo/001/gvo00170/vsc43896/task_transform_fmri/results_fmriprep/preproc_sdc-syn'
GM_prob_suffix = '_acq-GIfMIT1MPRAGE_run-1_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz'

#%% Define subj from the subjs.csv file

# get input from command line and check its length
cl = sys.argv[1:]
assert len(cl) == 1

# assign the input as subject number
subj = int(cl[0])

#%% perform the Searchlight Analysis for each participant

t_start = time.time()
    
bids_subj = 'sub-00' + str(subj)
subj_dir = os.path.join(FLM_root_dir, bids_subj)

# load the 4D image
d4_file = funcs_hpc.search_files(subj_dir, bids_subj, d4_files_suffix)[0] # the name of the 4d data
d4_path = os.path.join(subj_dir, d4_file)
d4_image = image.load_img(d4_path)

# load the label file, define the labels and groups
labels_file = d4_file.replace('4D.nii', 'labels.txt')
labels_path = os.path.join(subj_dir, labels_file)
labels_df = pd.read_csv(labels_path, sep=",")

labels_all = labels_df['tasks'].to_numpy().copy()
groups = labels_df['runs'].to_numpy().copy()
if need_regroup:
    groups = funcs_hpc.from_4runs_to_8groups(labels_all)
    
# load the brain mask from SPM results
spm_mask_path = os.path.join(subj_dir, spm_mask_file)
spm_mask = image.load_img(spm_mask_path)

# create the grey matter mask
GM_prob_dir = os.path.join(preproc_root_dir, bids_subj, 'anat')

if os.path.exists(GM_prob_dir):
    pass
else:
    GM_prob_dir = os.path.join(preproc_root_dir, bids_subj, 'ses-mri01' ,'anat')
    
GM_crop_mask_file = bids_subj + '_GM-mask-0.2-crop.nii'  # cropped mask
GM_crop_mask_path = os.path.join(subj_dir, GM_crop_mask_file) # in the 4d file folder instead of in the preprocessing folder

if os.path.isfile(GM_crop_mask_path): # if the cropped mask file already exist, then load that file 
    GM_mask_crop = image.load_img(GM_crop_mask_path)
else:  # create and save the cropped mask if it does not exist yet
    GM_prob_file = bids_subj + GM_prob_suffix
    GM_prob_path = os.path.join(GM_prob_dir, GM_prob_file)
    
    if os.path.isfile(GM_prob_path):
        pass 
    else:
        GM_prob_file = bids_subj + '_ses-mri01' + GM_prob_suffix
        GM_prob_path = os.path.join(GM_prob_dir, GM_prob_file)
    
    GM_prob_img = image.load_img(GM_prob_path)
    GM_prob_img_resample = image.resample_to_img(GM_prob_img, d4_image, interpolation='continuous')

    GM_mask = image.math_img('(img > 0.2).astype(np.int32)', img=GM_prob_img_resample)
    GM_mask_crop = image.math_img('(np.multiply(img1, img2)).astype(np.int32)', img1=GM_mask, img2=spm_mask)
    nib.nifti1.save(GM_mask_crop, GM_crop_mask_path)
  
# define decoder and cross-validation scheme
clf = svm.SVC(kernel='linear')
logo = LeaveOneGroupOut()

# Run the SearchLight Analysis
sl = SearchLight(
    mask_img=spm_mask,
    process_mask_img=GM_mask_crop,
    radius=5,  # 5 mm radius(2 voxels)
    estimator=clf,
    scoring='accuracy',
    cv=logo,
    verbose=False  # print a progress bar while fitting
)

sl.fit(d4_image, labels_all, groups=groups)

# saving the decoding result in a numpy array
sl_acc = sl.scores_ # the accuray of decoding for each voxel for this participant in np array
# sl_acc_all[:,:,:,subj_idx] = sl_acc

# saving the decoding result as a nifti image for each participant
sl_acc_img = image.new_img_like(spm_mask, sl_acc)
sl_acc_file = bids_subj + '_' + smth_suffix + '_spmT_' + d4_files_suffix.replace('4D.nii', 'SL_acc.nii')
sl_acc_path = os.path.join(output_dir, sl_acc_file)
nib.nifti1.save(sl_acc_img, sl_acc_path)

t_elapsed = round((time.time() - t_start)/60, 2) # scale is minute
print(f"finished decoding for {bids_subj}, which took {t_elapsed} mins")
