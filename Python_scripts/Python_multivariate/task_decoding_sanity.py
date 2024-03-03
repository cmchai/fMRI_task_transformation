#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 13:49:53 2024
ROI decoding on the task transformation paradigm --- sanity check about button press(left vs right)
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

FLM_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/sanity_check/results'
ROIs_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/sanity_check/results/sub-003'
ROI = 'b_motor.nii'

subjs_all = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50)))
subjs_torun = 3
subj = 3

bids_subj = 'sub-00' + str(subj)

subj_dir = os.path.join(FLM_root_dir, bids_subj)

#%% import and prepare data and labels for the subsequent decoding

# prepare the image
d4_file = 'sub-003_beta_keypress_4D.nii' 
d4_path = os.path.join(subj_dir, d4_file)
big_4D = image.load_img(d4_path)

# load in ROI image
ROI_path = os.path.join(ROIs_dir, ROI)
ROI_3d = image.load_img(ROI_path)
plotting.plot_roi(ROI_3d)

# resample the ROI to match the 4D feature file
ROI_3d_resamp = image.resample_to_img(ROI_3d, big_4D, interpolation='nearest')
plotting.plot_roi(ROI_3d_resamp)

ROI_3d_resamp_data = image.get_data(ROI_3d_resamp)
np.unique(ROI_3d_resamp_data) # check the unique values, should be 0 and 1
np.sum(ROI_3d_resamp_data) # how many voxels in the ROI

R_motor = masking.apply_mask(big_4D, ROI_3d_resamp)
R_motor = R_motor[:,np.sum(R_motor, axis=0) != 0] # get rid of the empty voxels

#%% import the labels and run infos

labels_file = 'sub-003_beta_keypress_labels.txt'
labels_path = os.path.join(subj_dir, labels_file)

labels = pd.read_csv(labels_path, sep=",")

S_all = labels['keys'].to_numpy().copy() # if not using copy(), when you change S_all, you also change the data frame:labels !!!!!
groups = labels['runs'].to_numpy().copy() 

#%% start the decoding with leave-one-run out approach

# initializing cross validation 
logo = LeaveOneGroupOut()
folds = logo.split(R_motor, S_all, groups=groups)
n_groups = np.unique(groups).size
acc_cv = np.zeros(n_groups)

# initializing decoder
# clf = LogisticRegression(solver='liblinear')
# clf = LinearDiscriminantAnalysis(solver='lsqr', shrinkage='auto')
clf = svm.SVC(kernel='linear')


for i, fold in enumerate(folds):
    train_idx, test_idx = fold
    print(train_idx, test_idx)
    
    R_train = R_motor[train_idx]
    S_train = S_all[train_idx]
    
    R_test = R_motor[test_idx]
    S_test = S_all[test_idx]
    print(S_test)
    
    clf.fit(R_train, S_train)
    preds = clf.predict(R_test)
    print(preds)
    
    acc_fold = accuracy_score(S_test, preds)
    acc_cv[i] = acc_fold
    
acc_cv_average = np.mean(acc_cv)
    

#%% try out all the self made functions

data_roi = extract_data_roi(big_4D, ROI_3d)
accuracies = logo_decoding(data_roi, S_all, groups, clf, confusion=False)
