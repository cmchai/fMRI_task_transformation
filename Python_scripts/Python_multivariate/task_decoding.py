#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 13:49:53 2024
ROI decoding on the task transformation paradigm
Script for trying out
@author: mengqiao
"""

import os

import numpy as np
import pandas as pd

from nilearn import image
from nilearn import plotting
from nilearn import masking
import nibabel as nib

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import sklearn.svm as svm
from sklearn.metrics import accuracy_score
from sklearn.model_selection import LeaveOneGroupOut

#%% general settings

FLM_root_dir = '/Volumes/extdrive/Task_Transform_GLM/GLM-02M-B/results'
ROIs_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI'

ROIs_schemes = ['GLM-02_con-5']
ROIs_scheme = ROIs_schemes[0]
ROIs_dir = os.path.join(ROIs_root_dir, ROIs_scheme)
ROIs = [file for file in os.listdir(ROIs_dir) if file.endswith('.nii')]
ROI = ROIs[1]

subjs_all = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50)))
subjs_torun = 3
subj = 3

bids_subj = 'sub-00' + str(subj)

subj_dir = os.path.join(FLM_root_dir, bids_subj)

#%% import and prepare data and labels for the subsequent decoding

# prepare the image
d4_file = 'sub-003_smth8_spmT_ALL-short_4D.nii' 
d4_path = os.path.join(subj_dir, d4_file)
big_4D = image.load_img(d4_path)

# load in ROI image
ROI_path = os.path.join(ROIs_dir, ROI)
ROI_3d = image.load_img(ROI_path)
plotting.plot_roi(ROI_3d,black_bg=True)
plotting.plot_glass_brain(ROI_3d)

# resample the ROI to match the 4D feature file
ROI_3d_resamp = image.resample_to_img(ROI_3d, big_4D, interpolation='nearest')
plotting.plot_roi(ROI_3d_resamp)

ROI_3d_resamp_data = image.get_data(ROI_3d_resamp)
np.unique(ROI_3d_resamp_data) # check the unique values, should be 0 and 1
np.sum(ROI_3d_resamp_data) # how many voxels in the ROI

R_lpfc = masking.apply_mask(big_4D, ROI_3d_resamp)
# aa1 = masking.apply_mask(big_4D, roi_resample)
# aat = np.transpose(aa)
# aat0 = np.nan_to_num(aat, nan=0)

#%% import the labels and run infos

labels_file = 'sub-003_smth8_spmT_ALL-short_labels.txt'
labels_path = os.path.join(subj_dir, labels_file)

labels = pd.read_csv(labels_path, sep=",")

S_all = labels['tasks'].to_numpy().copy()
groups = labels['runs'].to_numpy().copy()

#%% start the decoding with leave-one-run out approach

# initializing cross validation 
logo = LeaveOneGroupOut()
folds = logo.split(R_lpfc, S_all, groups=groups)
n_groups = np.unique(groups).size
acc_cv = np.zeros(n_groups)

# initializing decoder
clf = LogisticRegression(solver='liblinear', multi_class='auto')
clf = LinearDiscriminantAnalysis(solver='lsqr', shrinkage='auto')
clf = svm.SVC(kernel='linear')


for i, fold in enumerate(folds):
    train_idx, test_idx = fold
    # print(train_idx, test_idx)
    
    R_train = R_lpfc[train_idx]
    S_train = S_all[train_idx]
    
    R_test = R_lpfc[test_idx]
    S_test = S_all[test_idx]
    print(S_test)
    
    clf.fit(R_train, S_train)
    preds = clf.predict(R_test)
    print(preds)
    
    acc_fold = accuracy_score(S_test, preds)
    acc_cv[i] = acc_fold
    
acc_cv_average = np.mean(acc_cv)
    

#%% try an atlas

glasser = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas/HCP-MMP1/HCP-MMP1_2mm.nii.gz'
glasser1 = image.load_img(glasser)
plotting.plot_roi(glasser1)

glasser_data = glasser1.get_fdata()

glasser_data_flat = glasser_data.flatten()

rois = np.arange(1,361)

roi = glasser_data == 5
roi_img = nib.Nifti1Image(roi, glasser1.affine, glasser1.header)
roi_idx = np.where(glasser_data==5)[0]

affine = glasser1.affine
affine_4d = big_4D.affine

# resample the glasser atlas to match the 4D feature file
roi_resample = image.resample_to_img(roi_img, big_4D, interpolation='nearest')
plotting.plot_roi(roi_resample)
roi_resample_data = roi_resample.get_fdata()

big_4D_data = big_4D.get_fdata()
aa = big_4D_data[roi_resample_data]

#%% try another way to extract data from ROI

glasser = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas/HCP-MMP1/HCP-MMP1_25mm.nii'
glasser1 = image.load_img(glasser)
plotting.plot_roi(glasser1)
glasser_map = glasser1.get_fdata()
glasser_roi5_mask = glasser_map == 5

ROI_data = big_4D_data[glasser_roi5_mask]
ROI_data_tr = np.transpose(ROI_data)


os.chdir('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas/schaefer_2018/schaefer_2018')
os.getcwd()
nib.nifti1.save(brain_atlas_sch_img_resamp, 'Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2.5mm.nii.gz')

#%% try searchlight

from nilearn.decoding import SearchLight

# create the brain mask
brain_mask = image.math_img('(img.sum(axis=3) != 0).astype(np.int32)', img=big_4D)
plotting.plot_roi(brain_mask)

# the SPM derived mask
spm_mask_file = "mask.nii"
spm_mask_path = os.path.join(subj_dir, spm_mask_file)
spm_mask = image.load_img(spm_mask_path)
plotting.plot_roi(spm_mask)

# create the grey matter mask
GM_prob_dir = "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/preprocess/results_fmriprep/new_results/prep_23.1.0/sub-003/anat"
GM_prob_file = "sub-003_acq-GIfMIT1MPRAGE_run-1_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz"
GM_prob_path = os.path.join(GM_prob_dir, GM_prob_file)
GM_prob_img = image.load_img(GM_prob_path)
GM_prob_img_resample = image.resample_to_img(GM_prob_img, big_4D, interpolation='continuous')

plotting.plot_roi(GM_prob_img, threshold = None, colorbar=True, cmap='gray')
plotting.plot_roi(GM_prob_img_resample,threshold = None, colorbar=True, cmap='gray')

GM_mask = image.math_img('(img > 0.2).astype(np.int32)', img=GM_prob_img_resample)
plotting.plot_roi(GM_mask)
nib.nifti1.save(GM_mask, "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/preprocess/results_fmriprep/new_results/prep_23.1.0/sub-003/anat/GM_mask_0.2.nii.gz")

GM_mask_crop = image.math_img('(np.multiply(img1, img2)).astype(np.int32)', img1=GM_mask, img2=spm_mask)
plotting.plot_roi(GM_mask_crop)

# the grey matter mask from nilearn
GM_mask2 = masking.compute_brain_mask(big_4D, threshold=0.2, mask_type="gm")
plotting.plot_roi(GM_mask2)

# start the searchlight decoding
sl = SearchLight(
    mask_img=spm_mask,
    process_mask_img=GM_mask_crop,
    radius=5,  # 5 mm radius(2 voxels)
    estimator=clf,
    n_jobs=1,  # use only 1 core (for your own analyses, you might want to increase this!)
    scoring='accuracy',
    cv=logo,
    verbose=True  # print a progressbar while fitting
)

sl.fit(big_4D, S_all, groups=groups)

# look at the result
sl_score_img = image.new_img_like(spm_mask, sl.scores_)
plotting.plot_stat_map(sl_score_img, threshold=0.11)
nib.nifti1.save(sl_score_img, "/Users/mengqiao/Documents/fMRI_task_transform/Results/sub-003_smth8_spmT_ALL-short_SL_acc.nii")

