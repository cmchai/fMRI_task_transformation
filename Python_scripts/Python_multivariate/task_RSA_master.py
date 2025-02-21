#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 10:06:23 2024

conjunctive task RSA script

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
from sklearn.linear_model import Ridge
from scipy import stats

#%% define ROIs from FPN using the Glasser atlas (Assem et al., 2020)
ROIs_FPN_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI/HCP-MMP1_resamp'
ROIs_FPN = [file for file in sorted(os.listdir(ROIs_FPN_dir)) if file.endswith('.nii')]
ROIs_FPN_paths = [os.path.join(ROIs_FPN_dir, ROI_FPN) for ROI_FPN in ROIs_FPN]
ROIs_FPN_imgs = [image.load_img(ROI_FPN_path) for ROI_FPN_path in ROIs_FPN_paths]

#%% or define ROIs based on the Schaefer atlas

#--- Load the ROI files from Schaefer et al., 2018
atlas_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas'
atlas_subdir = 'schaefer_2018/schaefer_2018'
atlas_name = 'Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2.5mm.nii.gz'
atlas_path = os.path.join(atlas_root_dir, atlas_subdir, atlas_name)
atlas_image = image.load_img(atlas_path)
atlas_image_map = atlas_image.get_fdata()

atlas_labels_name = 'Schaefer2018_400Parcels_17Networks_order.txt'
atlas_labels_path = os.path.join(atlas_root_dir, atlas_subdir, atlas_labels_name)
atlas_labels_df = pd.read_csv(atlas_labels_path, delimiter = "\t", header = None)
atlas_labels = np.array(atlas_labels_df[1].tolist().copy())
atlas_labels = np.insert(atlas_labels, 0, "Background")

#--- Select ROIs from specific networks from this atlas
atlas_labels_df_contAB = atlas_labels_df[atlas_labels_df[1].str.contains("ContA") | atlas_labels_df[1].str.contains("ContB")]
ROIs_contAB = atlas_labels_df_contAB[0].tolist()

atlas_labels_df_defB_PFCv = atlas_labels_df[atlas_labels_df[1].str.contains("DefaultB_PFCv")]
ROIs_defB_PFCv = atlas_labels_df_defB_PFCv[0].tolist()

ROIs_Schaefer_FPN = ROIs_contAB + ROIs_defB_PFCv

#%% define d4 files that need to be extracted (first dimension of the matrix)
full_task_labs = np.arange(1,10)
d4_files_suffixes = ['RG-long-c1_4D.nii', 'RG-long-c2_4D.nii', 'TF-long-c1_4D.nii', 'TF-long-c2_4D.nii']
FLM_root_dir = os.path.join('/Volumes/extdrive/Task_Transform_GLM','GLM-02M','results') # GLM-02M:unsmoothed, GLM-02M-B:8mm smoothed

#%% Define all the ROIs to be decoded (second dimension of the matrix)
ROIs = ROIs_FPN              # should be a list of numbers(Schaefer atlas) or a list of images(Glasser atlas)
ROI_scheme = "Glasser_fpn"   # "Schaefer", "Glasser_fpn" and etc.
ROIs_imgs = ROIs_FPN_imgs    # only applicable when the schaefer atalas is not used, otherwise, should be defined as a empty list

#%% Define all the subjects to run (third dimension of the matrix)
subjs_all = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50)))
subjs = subjs_all.copy()
# subjs = np.setdiff1d(subjs_all, subjs_run)

#%% initiate the big RDM matrix (condition * ROI * subject * 36 * 36)
big_RDM = np.zeros((len(d4_files_suffixes), len(ROIs), subjs.shape[0], 36, 36), dtype=float)

#%% fill in the RDM

for d4_idx, d4_suffix in enumerate(d4_files_suffixes):
    
    for subj_idx, subj in enumerate(subjs):

        # define subject folder
        bids_subj = 'sub-00' + str(subj)
        subj_dir = os.path.join(FLM_root_dir, bids_subj)        

        # load d4 image
        d4_name = funcs.search_files(subj_dir, bids_subj, d4_suffix)[0] # should only return one element from a size of 1 list
        d4_path = os.path.join(subj_dir, d4_name)
        d4_image = image.load_img(d4_path)
        d4_data = d4_image.get_fdata()                  
        
        # load label files, get group and task ID information
        labels_file = d4_name.replace('4D.nii', 'labels.txt')
        labels_path = os.path.join(subj_dir, labels_file)
        labels_df = pd.read_csv(labels_path, sep=",")
        labels = labels_df['tasks'].to_numpy().copy()
        groups = labels_df['runs'].to_numpy().copy()
        
        print(f"RDM building for {d4_suffix} for the {subj_idx}th subject")
        
        for ROI_idx, ROI in enumerate(ROIs):
            
            if type(ROI) == int: # if ROI belongs to Schaefer 2018 atlas          
                # get pattern from a ROI, result is an array(samples * features)
                data_ROI = funcs.extract_data_roi_byatlas(d4_data, atlas_image_map, ROI)
            
            elif type(ROI) == str: # if ROI belongs to Glasser atlas
                
                if not ROIs_imgs: # if ROIs_img is an empty list
                    print("ROIs_img is an empty list")
                    sys.exit("errors about ROI imgs")      
                else:
                    ROI_image = ROIs_imgs[ROI_idx]                
                
                # get pattern from a ROI, result is an array(samples * features)
                data_ROI = funcs.extract_data_roi(d4_image, ROI_image)
            else:
                print('error:ROI is neither a number nor a string!')
                sys.exit("errors about ROIs")
            
            # fill nan in the missing patterns
            data_rsa = funcs.fill_nan_in(data_ROI, groups, labels, full_task_labs)
            
            # construct RDM using Spearman correlation
            spr_corr, p_spr_corr = stats.spearmanr(data_rsa, axis=1, nan_policy='propagate')
            rdm_spr = 1 - spr_corr
            
            # fill in the big RDM
            big_RDM[d4_idx, ROI_idx, subj_idx, :, :] = rdm_spr

#%% save big RDM

if ROI_scheme == "Schaefer":
    results_np_name = 'big_RDM_schaefer_sorted'
    ROIs_np_name = 'ROIs_schaefer_sorted'
elif ROI_scheme == "Glasser_fpn":
    results_np_name = 'big_RDM_sorted'
    ROIs_np_name = 'ROIs_sorted'
else:
    print('ROI scheme cannot be found when trying to save the numpy array data')
    sys.error('Data naming error, data not saved')

data_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA'
np.save(os.path.join(data_dir, results_np_name), big_RDM)
np.save(os.path.join(data_dir, ROIs_np_name), ROIs)

#%% Contructing relevant model RDM

# the archetype rdm
archetype_rdm = funcs.gen_archetype_rdm(full_task_labs, 4)

## within task bool and and their uptriangular bool ##
within_task_bool = funcs.same_conjunc_vec(archetype_rdm)
within_task_bool_up = np.triu(within_task_bool, k=1)
# plt.imshow(within_task_bool_up, cmap='viridis', interpolation='none')

## same stim bool and and their uptriangular bool (! without conjunc) ##
same_stim_bool = funcs.same_stim_vec(archetype_rdm)
same_stim_no_conjunc_bool = same_stim_bool & ~within_task_bool
same_stim_no_conjunc_bool_up = np.triu(same_stim_no_conjunc_bool, k=1)
# plt.imshow(same_stim_no_conjunc_bool_up, cmap='viridis', interpolation='none')

## same rule bool and and their uptriangular bool (! without conjunc) ##
same_rule_bool = funcs.same_rule_vec(archetype_rdm)
same_rule_no_conjunc_bool = same_rule_bool & ~within_task_bool
same_rule_no_conjunc_bool_up = np.triu(same_rule_no_conjunc_bool, k=1)
plt.imshow(same_rule_no_conjunc_bool_up, cmap='viridis', interpolation='none')

## no task overlap bool and and their uptriangular bool ##
no_same_bool = ~within_task_bool & ~same_stim_bool &~same_rule_bool
no_same_bool_up = np.triu(no_same_bool, k=1)
plt.imshow(no_same_bool_up, cmap='viridis', interpolation='none')

## check
full_up = within_task_bool_up | same_stim_no_conjunc_bool_up | same_rule_no_conjunc_bool_up | no_same_bool_up
plt.imshow(full_up, cmap='viridis', interpolation='none')


###############################################################################
######## -----run the following analysis based on the big RDM ------- #########
###############################################################################

#%% import big RDM is not imported yet

data_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA'
ROI_scheme = "Glasser_fpn" # "Schaefer" or "Glasser_fpn"

if ROI_scheme == "Schaefer":
    results_np_file = 'big_RDM_schaefer_sorted.npy'
    ROIs_np_file = 'ROIs_schaefer_sorted.npy'
elif ROI_scheme == "Glasser_fpn":
    results_np_file = 'big_RDM_sorted.npy'
    ROIs_np_file = 'ROIs_sorted.npy'
else:
    print('ROI scheme cannot be identified')
    sys.error('Data naming error, data not saved')

big_RDM = np.load(os.path.join(data_dir, results_np_file))
ROIs = np.load(os.path.join(data_dir, ROIs_np_file))

#%% Contructing relevant model RDM

# the archetype rdm
archetype_rdm = funcs.gen_archetype_rdm(full_task_labs, 4)

## within task bool and and their uptriangular bool ##
within_task_bool = funcs.same_conjunc_vec(archetype_rdm)
within_task_bool_up = np.triu(within_task_bool, k=1)
# plt.imshow(within_task_bool_up, cmap='viridis', interpolation='none')
# plt.colorbar()
# aa = archetype_rdm[within_task_bool]

## same stim bool and and their uptriangular bool (! without conjunc) ##
same_stim_bool = funcs.same_stim_vec(archetype_rdm)
same_stim_no_conjunc_bool = same_stim_bool & ~within_task_bool
same_stim_no_conjunc_bool_up = np.triu(same_stim_no_conjunc_bool, k=1)
# plt.imshow(same_stim_no_conjunc_bool_up, cmap='viridis', interpolation='none')
# plt.colorbar()
# bb = archetype_rdm[same_stim_no_conjunc_bool]

## same rule bool and and their uptriangular bool (! without conjunc) ##
same_rule_bool = funcs.same_rule_vec(archetype_rdm)
same_rule_no_conjunc_bool = same_rule_bool & ~within_task_bool
same_rule_no_conjunc_bool_up = np.triu(same_rule_no_conjunc_bool, k=1)
# plt.imshow(same_rule_no_conjunc_bool_up, cmap='viridis', interpolation='none')
# plt.colorbar()
# cc = archetype_rdm[same_rule_no_conjunc_bool]

## same (rule or stim) bool and and their uptriangular bool (! without conjunc) ##
same_rule_or_stim_bool = same_rule_no_conjunc_bool | same_stim_no_conjunc_bool
same_rule_or_stim_bool_up = same_stim_no_conjunc_bool_up | same_rule_no_conjunc_bool_up
# plt.imshow(same_rule_or_stim_bool, cmap='viridis', interpolation='none')

## no task overlap bool and and their uptriangular bool ##
no_same_bool = ~within_task_bool & ~same_stim_bool &~same_rule_bool
no_same_bool_up = np.triu(no_same_bool, k=1)
# plt.imshow(no_same_bool_up, cmap='viridis', interpolation='none')

## check
except_conjunc = same_rule_or_stim_bool_up | no_same_bool_up
full_up = within_task_bool_up | same_stim_no_conjunc_bool_up | same_rule_no_conjunc_bool_up | no_same_bool_up
# plt.imshow(except_conjunc, cmap='viridis', interpolation='none')

## same run bool and their uptriangular bool
same_run_bool = funcs.same_run(archetype_rdm)
same_run_bool_up = np.triu(same_run_bool, k=1)
# plt.imshow(same_run_bool_up, cmap='viridis', interpolation='none')
# plt.colorbar()

#%% extracting within task similarity and between task smilarity values

# create the lists to store the data
condition_df = []
ROI_df = []
subj_df = []
task_relation_df = []
dist_df = []

for d4_idx, d4_suffix in enumerate(d4_files_suffixes):
    
    for subj_idx, subj in enumerate(subjs):
        
        for ROI_idx, ROI in enumerate(ROIs):
            
            RDM = np.squeeze(big_RDM[d4_idx, ROI_idx, subj_idx, : ,:])
            
            within_elements = np.round(RDM[within_task_bool_up],5)
            condition_df.extend([d4_suffix] * within_elements.shape[0])
            ROI_df.extend([ROI] * within_elements.shape[0])
            subj_df.extend([subj] * within_elements.shape[0])
            task_relation_df.extend(["same_conjunc"] * within_elements.shape[0])
            dist_df.extend(within_elements)
            
            same_stim_elements = np.round(RDM[same_stim_no_conjunc_bool_up],5)
            condition_df.extend([d4_suffix] * same_stim_elements.shape[0])
            ROI_df.extend([ROI] * same_stim_elements.shape[0])
            subj_df.extend([subj] * same_stim_elements.shape[0])
            task_relation_df.extend(["same_stim"] * same_stim_elements.shape[0])
            dist_df.extend(same_stim_elements)
            
            same_rule_elements = np.round(RDM[same_rule_no_conjunc_bool_up],5)
            condition_df.extend([d4_suffix] * same_rule_elements.shape[0])
            ROI_df.extend([ROI] * same_rule_elements.shape[0])
            subj_df.extend([subj] * same_rule_elements.shape[0])
            task_relation_df.extend(["same_rule"] * same_rule_elements.shape[0])
            dist_df.extend(same_rule_elements)
            
            same_no_elements = np.round(RDM[no_same_bool_up],5)
            condition_df.extend([d4_suffix] * same_no_elements.shape[0])
            ROI_df.extend([ROI] * same_no_elements.shape[0])
            subj_df.extend([subj] * same_no_elements.shape[0])
            task_relation_df.extend(["same_nothing"] * same_no_elements.shape[0])
            dist_df.extend(same_no_elements)


# create and save the date frame
columns = {"condition":condition_df, "ROI":ROI_df, "subject":subj_df, "task_relation":task_relation_df, "distance":dist_df}
distance_df = pd.DataFrame(columns)
distance_df.head(10)

if ROI_scheme == "Schaefer":
    distance_df_name = 'results_RDM_schaefer_sorted.csv'
elif ROI_scheme == "Glasser_fpn":
    distance_df_name = 'results_RDM_glasser_sorted.csv'
else:
    print('ROI scheme cannot be found when trying to save the distance data')
    sys.error('Data naming error, data not saved')

distance_df.to_csv(os.path.join(data_dir, distance_df_name))

#%% Correlating the BIG(36*36 matrix) model and task RDM

# contruct task-by-task(36*36) model RDM
big_model_RSM = within_task_bool
big_model_RDM = (~big_model_RSM)*1
plt.imshow(big_model_RDM, cmap='viridis', interpolation='none')
plt.colorbar()

# only retain the up-triangular part of model RDM in a 1-d numpy array
idx_uptri_big = np.triu_indices(big_model_RDM.shape[0], k=1)
big_model_RDM_uptri = big_model_RDM[idx_uptri_big] # returns a 1-d numpy array

# compute the correlation(spearman's rho) between BIG model and empirical RDM
condition_df = []
ROI_df = []
subj_df = []
mdl_corr_df = []
mdl_corr_p_df = []


for d4_idx, d4_suffix in enumerate(d4_files_suffixes):
    
   for ROI_idx, ROI in enumerate(ROIs):
    
       for subj_idx, subj in enumerate(subjs):
           
           task_RDM = np.squeeze(big_RDM[d4_idx, ROI_idx, subj_idx, :, :])
           task_RDM_uptri = task_RDM[idx_uptri_big] # returns a 1-d numpy array
           
           mdl_corr, mdl_corr_p = stats.spearmanr(task_RDM_uptri, big_model_RDM_uptri, nan_policy='omit', alternative='two-sided')
           
           # concatenate the data frame
           condition_df.extend([d4_suffix])
           ROI_df.extend([ROI])
           subj_df.extend([subj])
           mdl_corr_df.extend([mdl_corr])
           mdl_corr_p_df.extend([mdl_corr_p])

# create and save the date frame for the big RDM correlation
columns = {"condition":condition_df, "ROI":ROI_df, "subject":subj_df, "spear_rho_empirical_model":mdl_corr_df, "spear_p_empirical_model":mdl_corr_p_df}
spear_empirical_model_df = pd.DataFrame(columns)

if ROI_scheme == "Schaefer":
    spear_empirical_model_df_name = 'spear_empirical_model_conjunc_schaefer_sorted.csv'
elif ROI_scheme == "Glasser_fpn":
    spear_empirical_model_df_name = 'spear_empirical_model_conjunc_glasser_sorted.csv'
else:
    print('ROI scheme cannot be identified')
    sys.error('Data naming error, data not saved')

spear_empirical_model_df.to_csv(os.path.join(data_dir, spear_empirical_model_df_name))

#%% Correlating the SMALL(9*9 matrix) model and task RDM (or RSM) to quantify the compositional coding

# construct task-by-task(9*9) empirical RDM or RSM  by averaging within-run distance or similarity across runs
small_RDM = funcs.big_rdm_to_small_between(big_RDM)
one_slice = np.squeeze(small_RDM[2,42,3,:,:])
plt.imshow(np.squeeze(small_RDM[2,42,3,:,:]), cmap='viridis', interpolation='none')

# contruct task-by-task(9*9) model RSM or RDM
big_model_RSM = same_rule_or_stim_bool
big_model_RDM = ~big_model_RSM
plt.imshow(big_model_RDM, cmap='viridis', interpolation='none')

small_model_RDM = funcs.big_rdm_to_small_within(big_model_RDM*1)
plt.imshow(small_model_RDM, cmap='viridis', interpolation='none')
plt.colorbar()

idx_uptri_small = np.triu_indices(small_model_RDM.shape[0], k=1)
small_model_RDM_uptri = small_model_RDM[idx_uptri_small] # returns a 1-d numpy array

# compute the correlation(spearman's rho) between model and empirical RDM or RSM
condition_df = []
ROI_df = []
subj_df = []
mdl_corr_df = []
mdl_corr_p_df = []


for d4_idx, d4_suffix in enumerate(d4_files_suffixes):
    
   for ROI_idx, ROI in enumerate(ROIs):
    
       for subj_idx, subj in enumerate(subjs):
           
           task_RDM = np.squeeze(small_RDM[d4_idx, ROI_idx, subj_idx, :, :])
           task_RDM_uptri = task_RDM[idx_uptri_small] # returns a 1-d numpy array
           
           mdl_corr, mdl_corr_p = stats.spearmanr(task_RDM_uptri, small_model_RDM_uptri, nan_policy='omit', alternative='two-sided')
           
           # concatenate the data frame
           condition_df.extend([d4_suffix])
           ROI_df.extend([ROI])
           subj_df.extend([subj])
           mdl_corr_df.extend([mdl_corr])
           mdl_corr_p_df.extend([mdl_corr_p])


# create and save the date frame for the RDM correlation
columns = {"condition":condition_df, "ROI":ROI_df, "subject":subj_df, "spear_rho_empirical_model":mdl_corr_df, "spear_p_empirical_model":mdl_corr_p_df}
spear_empirical_model_df = pd.DataFrame(columns)

if ROI_scheme == "Schaefer":
    spear_empirical_model_df_name = 'spear_empirical_model_bet_run_schaefer_sorted.csv'
elif ROI_scheme == "Glasser_fpn":
    spear_empirical_model_df_name = 'spear_empirical_model_bet_run_glasser_sorted.csv'
else:
    print('ROI scheme cannot be found when trying to save the distance data')
    sys.error('Data naming error, data not saved')

spear_empirical_model_df.to_csv(os.path.join(data_dir, spear_empirical_model_df_name))

#%% Correlating the SMALL(9*9 matrix) conjunctive model and data RDM (or RSM) to quantify the conjunctive coding

# construct task-by-task(9*9) empirical RDM or RSM  by averaging between-run distance or similarity across run pairs
small_RDM_bet = funcs.big_rdm_to_small_between(big_RDM)
# plt.imshow(np.squeeze(small_RDM_bet[2,42,3,:,:]), cmap='viridis', interpolation='none')
# plt.colorbar()

# contruct task-by-task(9*9) conjunctive model RSM or RDM
big_conjunc_model_RSM = within_task_bool
big_conjunc_model_RDM = ~big_conjunc_model_RSM
# plt.imshow(big_conjunc_model_RDM, cmap='viridis', interpolation='none')
# plt.colorbar()

small_conjunc_model_RDM = funcs.big_rdm_to_small_between(big_conjunc_model_RDM*1)
# plt.imshow(small_conjunc_model_RDM, cmap='viridis', interpolation='none')
# plt.colorbar()

# since we only use between-run distance, we don't have to only retrieve the uptriangular part of it, but the whole RDM
small_conjunc_model_RDM_vec = small_conjunc_model_RDM.flatten(order='C') 

# compute the correlation(spearman's rho) between model and empirical RDM or RSM
condition_df = []
ROI_df = []
subj_df = []
mdl_corr_df = []
mdl_corr_p_df = []


for d4_idx, d4_suffix in enumerate(d4_files_suffixes):
    
   for ROI_idx, ROI in enumerate(ROIs):
    
       for subj_idx, subj in enumerate(subjs):
           
           task_RDM = np.squeeze(small_RDM_bet[d4_idx, ROI_idx, subj_idx, :, :])
           task_RDM_vec = task_RDM.flatten(order='C')
           
           mdl_corr, mdl_corr_p = stats.spearmanr(task_RDM_vec, small_conjunc_model_RDM_vec, nan_policy='omit', alternative='two-sided')
           
           # concatenate the data frame
           condition_df.extend([d4_suffix])
           ROI_df.extend([ROI])
           subj_df.extend([subj])
           mdl_corr_df.extend([mdl_corr])
           mdl_corr_p_df.extend([mdl_corr_p])

# create and save the date frame for the RDM correlation
columns = {"condition":condition_df, "ROI":ROI_df, "subject":subj_df, "spear_rho_empirical_model":mdl_corr_df, "spear_p_empirical_model":mdl_corr_p_df}
spear_empirical_conjunc_model_df = pd.DataFrame(columns)

if ROI_scheme == "Schaefer":
    spear_empirical_conjunc_model_df_name = 'spear_empirical_conjunc_model_schaefer_sorted.csv'
elif ROI_scheme == "Glasser_fpn":
    spear_empirical_conjunc_model_df_name = 'spear_empirical_conjunc_model_glasser_sorted.csv'
else:
    print('ROI scheme cannot be found when trying to save the distance data')
    sys.error('Data naming error, data not saved')

spear_empirical_conjunc_model_df.to_csv(os.path.join(data_dir, spear_empirical_conjunc_model_df_name))

#%% Try the regression apporach (empirical RDM = intercept + conjunc RDM + stim RDM + rule RDM + run RDM)

# To determine whether standardize both y and X or not
standardize_xy = True

# Create all the relevant model RDMs and put them into one design matrix
idx_uptri_big = np.triu_indices(archetype_rdm.shape[0], k=1)

conjunc_rdm_up = (1 * (~ within_task_bool))[idx_uptri_big]
compo_stim_rdm_up = 1 * (~ same_stim_no_conjunc_bool)[idx_uptri_big]
compo_rule_rdm_up = 1 * (~ same_rule_no_conjunc_bool)[idx_uptri_big]
run_rdm_up = 1 * (~ same_rule_bool)[idx_uptri_big]

X = np.stack([conjunc_rdm_up, compo_stim_rdm_up, compo_rule_rdm_up, run_rdm_up], axis=1)

# Fitting the regression model and save all the coefficients into a dataframe
condition_df = []
ROI_df = []
subj_df = []
r_2 = []
beta_conjunc = []
beta_stim = []
beta_rule = []
beta_run = []

for d4_idx, d4_suffix in enumerate(d4_files_suffixes):
    
   for ROI_idx, ROI in enumerate(ROIs):
    
       for subj_idx, subj in enumerate(subjs):
           
           task_RDM = np.squeeze(big_RDM[d4_idx, ROI_idx, subj_idx, :, :])
           y = task_RDM[idx_uptri_big] # returns a 1-d numpy array
           
           # get rid of nan values
           y_notnan_bool = ~np.isnan(y)
           X_excl_nan = X[y_notnan_bool, :]
           y_excl_nan = y[y_notnan_bool]
           
           # skip the regression if usable data is too little
           if y_excl_nan.shape[0] < 100:      
               continue           
           
           # fit the regression model
           reg = Ridge(alpha=1.0, fit_intercept=True, positive=False)
           reg.fit(X_excl_nan, y_excl_nan)         
           r2_score = reg.score(X_excl_nan, y_excl_nan)
           
           # concatenate the data frame
           r_2.extend([np.round(r2_score, 5)])
           beta_conjunc.extend([np.round(reg.coef_[0], 5)])
           beta_stim.extend([np.round(reg.coef_[1], 5)])
           beta_rule.extend([np.round(reg.coef_[2], 5)])
           beta_run.extend([np.round(reg.coef_[3], 5)])
           
           condition_df.extend([d4_suffix])
           ROI_df.extend([ROI])
           subj_df.extend([subj])

#%% tryouts
within_task_rsm = within_task_bool*1

plt.imshow(within_task_rsm, cmap='viridis', interpolation='none')
plt.colorbar()

np.fill_diagonal(within_task_rsm, 2)

small_within = funcs.big_rdm_to_small_within(within_task_rsm)
plt.imshow(small_within, cmap='viridis', interpolation='none')
plt.colorbar()

small_between = funcs.big_rdm_to_small_between(within_task_rsm)
plt.imshow(small_between, cmap='viridis', interpolation='none')
plt.colorbar()

#%% for visualization
typical_big_RDM = np.nanmean(big_RDM, axis=(0,1,2)) - 0.2
np.fill_diagonal(typical_big_RDM, 0)
plt.imshow(typical_big_RDM, cmap='viridis', interpolation='none')
plt.colorbar()

noise = np.random.normal(loc=0, scale=0.05, size=typical_big_RDM.shape)
np.fill_diagonal(noise, 0)
noisy_big_RDM = typical_big_RDM + noise


fig, ax = plt.subplots()
im = ax.imshow(noisy_big_RDM, cmap='viridis', interpolation='none')
ax.set_xticks([9,18,27,36])
ax.set_yticks([9,18,27,36])
fig.colorbar(im)


# compositional model RDM
fig, ax = plt.subplots()
im = ax.imshow(small_model_RDM, cmap='viridis', interpolation='none')
ax.set_xticks([0,1,2,3,4,5,6,7,8])
ax.set_yticks([0,1,2,3,4,5,6,7,8])
fig.colorbar(im)

# conjunctive model RDM
fig, ax = plt.subplots()
im = ax.imshow(small_conjunc_model_RDM, cmap='viridis', interpolation='none')
ax.set_xticks([0,1,2,3,4,5,6,7,8])
ax.set_yticks([0,1,2,3,4,5,6,7,8])
fig.colorbar(im)

