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
ROIs = ROIs_Schaefer_FPN     # should be a list of numbers(Schaefer atlas) or a list of images(Glasser atlas)
ROI_scheme = "Schaefer"      # "Schaefer", "Glasser_fpn" and etc.
ROIs_imgs = []               # only application when the schaefer atalas is not used, otherwise, should be defined as a empty list

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
ROI_scheme == "Schaefer" # "Schaefer" or "Glasser_fpn"

if ROI_scheme == "Schaefer":
    results_np_file = 'big_RDM_schaefer_sorted.npy'
    ROIs_np_file = 'ROIs_schaefer_sorted.npy'
elif ROI_scheme == "Glasser_fpn":
    results_np_file = 'big_RDM_sorted.npy'
    ROIs_np_file = 'ROIs_sorted.npy'
else:
    print('ROI scheme cannot be found when trying to save the numpy array data')
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

## same stim bool and and their uptriangular bool (! without conjunc) ##
same_stim_bool = funcs.same_stim_vec(archetype_rdm)
same_stim_no_conjunc_bool = same_stim_bool & ~within_task_bool
same_stim_no_conjunc_bool_up = np.triu(same_stim_no_conjunc_bool, k=1)
# plt.imshow(same_stim_no_conjunc_bool_up, cmap='viridis', interpolation='none')

## same rule bool and and their uptriangular bool (! without conjunc) ##
same_rule_bool = funcs.same_rule_vec(archetype_rdm)
same_rule_no_conjunc_bool = same_rule_bool & ~within_task_bool
same_rule_no_conjunc_bool_up = np.triu(same_rule_no_conjunc_bool, k=1)
# plt.imshow(same_rule_no_conjunc_bool_up, cmap='viridis', interpolation='none')

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


#%% extracting within task similarity and between task smilarity values

## create the lists to store the data
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


#%% create and save the date frame
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

#%% create task-by-task RDM(or RSM) 9*9 matrix by averaging distance or similarity across runs

# construct task-by-task(9*9) empirical RDM or RSM
small_RDM = funcs.big_rdm_to_small(big_RDM)
one_slice = np.squeeze(small_RDM[2,42,3,:,:])
plt.imshow(np.squeeze(small_RDM[2,42,3,:,:]), cmap='viridis', interpolation='none')

# contruct task-by-task(9*9) model RSM or RDM
big_model_RSM = same_rule_or_stim_bool
big_model_RDM = ~big_model_RSM
plt.imshow(big_model_RDM, cmap='viridis', interpolation='none')

small_model_RDM = funcs.big_rdm_to_small(big_model_RDM*1)
plt.imshow(small_model_RDM, cmap='viridis', interpolation='none')

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


#%% create and save the date frame for the RDM correlation
columns = {"condition":condition_df, "ROI":ROI_df, "subject":subj_df, "spear_rho_empirical_model":mdl_corr_df, "spear_p_empirical_model":mdl_corr_p_df}
spear_empirical_model_df = pd.DataFrame(columns)

if ROI_scheme == "Schaefer":
    spear_empirical_model_df_name = 'spear_empirical_model_schaefer_sorted.csv'
elif ROI_scheme == "Glasser_fpn":
    spear_empirical_model_df_name = 'spear_empirical_model_glasser_sorted.csv'
else:
    print('ROI scheme cannot be found when trying to save the distance data')
    sys.error('Data naming error, data not saved')

spear_empirical_model_df.to_csv(os.path.join(data_dir, spear_empirical_model_df_name))
