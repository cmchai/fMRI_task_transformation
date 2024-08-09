#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 12:24:53 2024

RSA on fMRI task transformation paradigm

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
import matplotlib.pyplot as plt

from nilearn import image, plotting
import nibabel as nib

from sklearn.metrics import pairwise_distances
from sklearn.linear_model import LinearRegression
from scipy import stats

#%% import 4d-data and labels, groups files

sub_dir = '/Volumes/extdrive/Task_Transform_GLM/GLM-02M/results/sub-002'
d4_name = 'sub-002_smthN_spmT_RG-long-c1_4D.nii'
d4_path = os.path.join(sub_dir, d4_name)
big_4D = image.load_img(d4_path)

labels_file = d4_name.replace('4D.nii', 'labels.txt')
labels_path = os.path.join(sub_dir, labels_file)
labels_df = pd.read_csv(labels_path, sep=",")

# retrieve the task identities
labels = labels_df['tasks'].to_numpy().copy()

# retrieve the group identities
groups = labels_df['runs'].to_numpy().copy()

#%% import ROI and get data from this ROI

ROI_path = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI/GLM-02_con-5/Clus_L_IFG.nii'
ROI_img = image.load_img(ROI_path)
data_roi = funcs.extract_data_roi(big_4D, ROI_img)

#%% deal with all the missing labels and corresponding neural patterns
len_rdm = 36

full_task_labs = np.arange(1,10)

groups_idx = funcs.unique_preserve_order(groups)
labels_bygroup = [labels[groups == group] for group in groups_idx]    # return a list of arrays, each array represent the task labels that exist in this group from the data
groups_bool = np.concatenate([np.in1d(full_task_labs, labels_1group) for labels_1group in labels_bygroup]) # HERE assume the labels per group is ASCENDING in the data !!!

data_roi_rsa = data_roi.copy()

for cell_idx, group_bool in enumerate(groups_bool):
    if not group_bool:
        data_roi_rsa = np.insert(data_roi_rsa, cell_idx, np.nan, axis=0)
    
#%% creating the RDM
# distance measure can be: 'nan_euclidean','correlation'

rdm = pairwise_distances(data_roi_rsa, metric='correlation', force_all_finite = 'allow-nan')

spr_corr, p_spr_corr = stats.spearmanr(data_roi_rsa, axis=1, nan_policy='propagate')
rdm_spr = 1 - spr_corr

plt.imshow(rdm_spr, cmap='viridis', interpolation='none')

# Add a colorbar to show the scale
plt.colorbar()

# Add titles and labels if needed
plt.title('distance')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')

#%% creating model RDM

task_vec = np.tile(full_task_labs, 4).astype(str)

n = task_vec.shape[0]

# Create a coordinate grid
row_indices, col_indices = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')

# Create the matrix by concatenating strings based on coordinates
mat_1 = np.core.defchararray.add(task_vec[row_indices], task_vec[col_indices])

print(mat_1)


# Function to create conjunction model RDM (by looking at if the first digit equals the last digit)
def same_conjunc(tasks):
    return tasks[0] == tasks[1]

same_conjunc_vec = np.vectorize(same_conjunc)

# Function to create stimulus type model RDM (by looking at if both first and last digit belongs to [])
def  same_stim(tasks):
    tasks_vec = np.array([tasks[0], tasks[1]], dtype=int)
    both_animal = np.in1d(tasks_vec, np.array([1,2,3], dtype=int)).all()
    both_place = np.in1d(tasks_vec, np.array([4,5,6], dtype=int)).all()
    both_vehicle = np.in1d(tasks_vec, np.array([7,8,9], dtype=int)).all()
    return both_animal or both_place or both_vehicle

same_stim_vec = np.vectorize(same_stim)

# Function to create task rule model RDM (by looking at if both first and last digit belongs to [])
def  same_rule(tasks):
    tasks_vec = np.array([tasks[0], tasks[1]], dtype=int)
    both_age = np.in1d(tasks_vec, np.array([1,4,7], dtype=int)).all()
    both_size = np.in1d(tasks_vec, np.array([2,5,8], dtype=int)).all()
    both_location = np.in1d(tasks_vec, np.array([3,6,9], dtype=int)).all()
    return both_age or both_size or both_location

same_rule_vec = np.vectorize(same_rule)

# Apply the function to the string array
conjunc_model_rdm = 1 * (~ same_conjunc_vec(mat_1))
stim_model_rdm = 1 * (~ same_stim_vec(mat_1))
rule_model_rdm = 1 * (~ same_rule_vec(mat_1))

plt.imshow(conjunc_model_rdm, cmap='viridis', interpolation='none')
plt.imshow(stim_model_rdm, cmap='viridis', interpolation='none')
plt.imshow(rule_model_rdm, cmap='viridis', interpolation='none')

# create block model RDM (patterns within a block is more similar to patterns between blocks)
block_model_rdm = np.zeros(rdm.shape, dtype=int)

for index, value in np.ndenumerate(block_model_rdm):
    block_model_rdm[index] = (index[0]//9) != (index[1]//9)

plt.imshow(block_model_rdm, cmap='viridis', interpolation='none')    

#%% take up-triangular part of neural and model RDMs

up_bool = np.triu(rdm, k=0)
plt.imshow(up_bool, cmap='viridis', interpolation='none')

triu_idx = np.triu_indices(rdm.shape[0], k=1)
n_triu_elements = 36*35/2 # without diagonal

rdm_triu = rdm[triu_idx] # may have nan value btw.
conjunc_model_triu = conjunc_model_rdm[triu_idx]
stim_model_triu = stim_model_rdm[triu_idx]
rule_model_triu = rule_model_rdm[triu_idx]
block_model_triu = block_model_rdm[triu_idx]

#%% Using linear regression to regress the neural RDM on model RDMs

# define the dependent and independent variables
X = np.stack([conjunc_model_triu, stim_model_triu, rule_model_triu, block_model_triu], axis=1)
y = rdm_triu

# get rid of nan values
y_notnan_bool = ~np.isnan(y)
X_excl_nan = X[y_notnan_bool, :]
y_excl_nan = y[y_notnan_bool]

# standardize both X and y (z score)
X_z = stats.zscore(X_excl_nan, axis=0)
y_z = stats.zscore(y_excl_nan, axis=0)

# check collinearity between model predictors and y
y_X_excl_nan = np.concatenate((y_excl_nan.reshape((630, 1)), X_excl_nan), axis=1)
y_X_z = np.concatenate((y_z.reshape((630, 1)), X_z), axis=1)
colin = np.corrcoef(np.transpose(y_X_excl_nan))

# pre-whiten the data (regress out the block )
pre_wh_reg = LinearRegression(fit_intercept=False, positive=False).fit(block_model_triu.reshape((630, 1)), y_excl_nan)
prediction = pre_wh_reg.predict(block_model_triu.reshape((630, 1)))
residual = y_excl_nan - prediction

# define and fit the model
reg = LinearRegression(fit_intercept=False, positive=False).fit(X_excl_nan, y_excl_nan) # without pre-whitening
reg = LinearRegression(fit_intercept=False, positive=False).fit(X_excl_nan[:,:-1], residual) # with pre-whitening

# check coefficients and intercept
reg.coef_
reg.intercept_
print(reg.summary())


