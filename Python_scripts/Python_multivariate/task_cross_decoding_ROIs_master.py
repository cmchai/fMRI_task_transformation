#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 11:33:03 2024

Script to run Cross-decoding(!!!!!) on a set of ROIs

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
# import numpy_indexed as npi
import pandas as pd

from nilearn import image, plotting
# from nilearn import masking
# from nilearn import datasets
import nibabel as nib

# from sklearn.preprocessing import StandardScaler
import sklearn.svm as svm
from sklearn.feature_selection import SelectKBest, f_classif
# from sklearn.metrics import accuracy_score
# from sklearn.metrics import confusion_matrix
# from sklearn.model_selection import LeaveOneGroupOut

from scipy import stats
import statsmodels.stats.multitest as multi

#%% Load the ROI files from Schaefer et al., 2018

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

#%% select ROIs from specific networks from this atlas
atlas_labels_df_contAB = atlas_labels_df[atlas_labels_df[1].str.contains("ContA") | atlas_labels_df[1].str.contains("ContB")]
ROIs_contAB = atlas_labels_df_contAB[0].tolist()

atlas_labels_df_defB_PFCv = atlas_labels_df[atlas_labels_df[1].str.contains("DefaultB_PFCv")]
ROIs_defB_PFCv = atlas_labels_df_defB_PFCv[0].tolist()

#%% define the index of the parcels that we will perform the decoding on
ROIs_Sch = [128, 129, 130, 131, 132, 330, 331, 332, 333, 334]    # 128-132 and 330-334 are DLPFC in control A network
ROIs_Sch_comb = [np.array([128, 129, 130, 131, 132]),np.array([330, 331, 332, 333, 334])]

#%% Define ROIs from the short trials decoding result
ROIs_Sch_short = [10,55,58,62,63,65,121,125,128,180,186,195,197,199,200,207,258,341]
ROIs_Sch_short_group = [10,np.array([55,58]),np.array([62,63,65]),np.array([121,125,128]),np.array([180,186]),np.array([195,197,199,200]),207,258,341]

#%% Define ROIs from AAL anatomical atlas
ROIs_AAL_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI/AAL3'
ROIs_AAL = [file for file in os.listdir(ROIs_AAL_dir) if file.endswith('.nii')]
ROIs_AAL_paths = [os.path.join(ROIs_AAL_dir, ROI_AAL) for ROI_AAL in ROIs_AAL]
ROIs_AAL_imgs = [image.load_img(ROI_AAL_path) for ROI_AAL_path in ROIs_AAL_paths]

#%% define ROIs from the first level contrast
ROIs_FLM_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI/GLM-02_con-5'
ROIs_FLM = [file for file in os.listdir(ROIs_FLM_dir) if file.endswith('.nii')]
ROIs_FLM_paths = [os.path.join(ROIs_FLM_dir, ROI_FLM) for ROI_FLM in ROIs_FLM]
ROIs_FLM_imgs = [image.load_img(ROI_FLM_path) for ROI_FLM_path in ROIs_FLM_paths]

#%% define ROIs from FPN using the Glasser atlas (Assem et al., 2020)
ROIs_FPN_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI/HCP-MMP1_resamp'
ROIs_FPN = [file for file in os.listdir(ROIs_FPN_dir) if file.endswith('.nii')]
ROIs_FPN_paths = [os.path.join(ROIs_FPN_dir, ROI_FPN) for ROI_FPN in ROIs_FPN]
ROIs_FPN_imgs = [image.load_img(ROI_FPN_path) for ROI_FPN_path in ROIs_FPN_paths]

# visualize the masking image
plotting.plot_glass_brain(ROIs_FPN_imgs[3], title='Glasser Parcel:' + ROIs_FPN[3])

#####################################
#%% Define all the ROIs to be decoded

ROIs = ROIs_contAB + ROIs_defB_PFCv
ROI_scheme = "Schaefer"   # "Schaefer", "Glasser_fpn" and etc.
# ROIs_imgs = ROIs_FPN_imgs # only application when the schaefer atalas is not used

#%% choosing whether using feature selection
feat_select = False
feat_num = 1500

#%% plotting ROIs

plotting_ROIs = False

if plotting_ROIs:
    # ploting parcels from Schaefer 2018
    for parcel_idx in ROIs_Sch:
        img_data = np.zeros(atlas_image_map.shape)
        img_data[np.isin(atlas_image_map, parcel_idx)] = 1
        img = nib.Nifti1Image(img_data, atlas_image.affine, atlas_image.header)
        plotting.plot_glass_brain(img, title='Schaefer Atlas Parcel No.' + str(parcel_idx))
    
    for region in ROIs_AAL:
        img_path = os.path.join(ROIs_AAL_dir, region)
        plotting.plot_glass_brain(img_path, title = region)

# the fpn network from Schaefer 2018
ROIs_to_plot = np.array([122,123,124,125,126, 128, 129, 130, 131, 132, 136, 137, 138, 140, 141, 142, 325, 326, 327, 328, 330, 331, 332, 333, 334, 338, 339, 340, 341, 346, 347, 348, 349, 183, 184, 185, 186, 187, 382, 383, 384])

img_data = np.zeros(atlas_image_map.shape)
img_data[np.isin(atlas_image_map, ROIs_to_plot)] = 1
img = nib.Nifti1Image(img_data, atlas_image.affine, atlas_image.header)
plotting.plot_glass_brain(img, title='Schaefer Atlas FPN network')
nib.nifti1.save(img, '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/independent_roi/conA_conB_defB/Schaefer_FPN_mask.nii')


#%% define what data suffix to be decoded(should be in 4D format)
compositional = False

if compositional:
    comp_suffix = '-Comp'
    d4_files_suffixes_1 = []
    d4_files_suffixes_2 = []
else:
    comp_suffix = ''
    d4_files_suffixes_1 = ['RG-short_4D.nii', 'RG-short_4D.nii', 'TF-short_4D.nii','TF-short_4D.nii'] # training
    d4_files_suffixes_2 = ['RG-long-c1_4D.nii', 'RG-long-c2_4D.nii', 'TF-long-c1_4D.nii', 'TF-long-c2_4D.nii'] # testing
    
n_d4_files = len(d4_files_suffixes_1)
glm_folder = 'GLM-02M' + comp_suffix
FLM_root_dir = os.path.join('/Volumes/extdrive/Task_Transform_GLM', glm_folder,'results') # GLM-02M:unsmoothed, GLM-02M-B:8mm smoothed

#%% Define all the subjects to run

subjs_all = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50)))
subjs = subjs_all[3:43]
# subjs = np.setdiff1d(subjs_all, subjs_run)

#%% Define tasks to be decoded
# can be one of the following:
# 1. np.arange(1,10), 9 conjunctive tasks
# 2. np.array(['animal', 'place', 'vehicle']), 3 compositional tasks in the stimulus type dimension
# 3. np.array(['age', 'location', 'size']), 3 compositional tasks in the task rule dimension

conjunc_labs = np.arange(1,10) # 9 conjunctive tasks
stim_labs = np.array(['animal', 'place', 'vehicle']) # 3 compositional tasks in the stimulus type dimension
rule_labs = np.array(['age', 'location', 'size'])    # 3 compositional tasks in the task rule dimension

#%% Define decoders and if the cross decoding is one sided(only train on data 1 and test on data 2 and not vice versa)
decoder = 'SVM'

oneSide = True

#%% Creating the dataframe to store the decoding results

columns = ['subject', 'decoder', 'smoothing', 'measure', 'training', 'testing', 'atlas',
           'roi', 'roi_size', 'n_folds', 'mean_accuracy']

results = pd.DataFrame(columns = columns)

#%% Performing the decoding across ROIs
t_start = time.time()
    
for subj_idx, subj in enumerate(subjs): # loop over subjects

    bids_subj = 'sub-00' + str(subj)
    subj_dir = os.path.join(FLM_root_dir, bids_subj)
    
    d4_files_1 = []
    d4_files_2 = []
    
    for i in range(len(d4_files_suffixes_1)): # loop over each pair of cross-decoding
        
        d4_files_suffix_1 = d4_files_suffixes_1[i]
        d4_files_suffix_2 = d4_files_suffixes_2[i]
        
        d4_file_name_1 = funcs.search_files(subj_dir, bids_subj, d4_files_suffix_1) # which is a list
        d4_files_1 = d4_files_1 + d4_file_name_1
        
        d4_file_name_2 = funcs.search_files(subj_dir, bids_subj, d4_files_suffix_2) # which is a list
        d4_files_2 = d4_files_2 + d4_file_name_2
    
    print(f"start cross-decoding for {bids_subj}")

    for ii in range(len(d4_files_1)): # loop over different kind of 4d data
        
        d4_file_1 = d4_files_1[ii]
        d4_file_2 = d4_files_2[ii]
        
        file_substrings_1 = d4_file_1.split('_')      # a list of substrings indicating all the information regarding the data
        file_substrings_2 = d4_file_2.split('_')      # a list of substrings indicating all the information regarding the data
        
        d4_path_1 = os.path.join(subj_dir, d4_file_1)
        d4_image_1 = image.load_img(d4_path_1)
        d4_data_1 = d4_image_1.get_fdata()
        
        d4_path_2 = os.path.join(subj_dir, d4_file_2)
        d4_image_2 = image.load_img(d4_path_2)
        d4_data_2 = d4_image_2.get_fdata()       
        
        labels_file_1 = d4_file_1.replace('4D.nii', 'labels.txt')
        labels_path_1 = os.path.join(subj_dir, labels_file_1)
        labels_df_1 = pd.read_csv(labels_path_1, sep=",")
        
        labels_file_2 = d4_file_2.replace('4D.nii', 'labels.txt')
        labels_path_2 = os.path.join(subj_dir, labels_file_2)
        labels_df_2 = pd.read_csv(labels_path_2, sep=",")

        # retrieve the task identities
        labels_all_1 = labels_df_1['tasks'].to_numpy().copy()
        labels_all_2 = labels_df_2['tasks'].to_numpy().copy()
        
        full_task_labs = conjunc_labs
        if not funcs.is_numerical(labels_all_1): # if the label are not numerical(e.g.1-9), but strings('animal', 'size'...)
            labels_all_1 = labels_all_1.astype(str)
            labels_all_2 = labels_all_2.astype(str)
            if 'animal' in labels_all_1:
                full_task_labs = stim_labs
            else:
                full_task_labs = rule_labs
        
        # retrieve the group identities
        groups_1 = labels_df_1['runs'].to_numpy().copy()
        groups_2 = labels_df_2['runs'].to_numpy().copy()

        for ROI_idx, ROI in enumerate(ROIs):
            
            if type(ROI) == int: # if ROI belongs to Schaefer 2018 atlas
                atlas = 'Schaefer'
                parcel_label = atlas_labels[ROI]              
                data_ROI_1 = funcs.extract_data_roi_byatlas(d4_data_1, atlas_image_map, ROI)
                data_ROI_2 = funcs.extract_data_roi_byatlas(d4_data_2, atlas_image_map, ROI)
            elif type(ROI) == np.ndarray: # an array for merging ROIs
                atlas = 'Schaefer'
                data_ROI_1 = funcs.extract_data_combined_roi_byatlas(d4_data_1, atlas_image_map, ROI)
                data_ROI_2 = funcs.extract_data_combined_roi_byatlas(d4_data_2, atlas_image_map, ROI)
            else: # string data
                atlas = ROI_scheme
                ROI_image = ROIs_imgs[ROI_idx]
                data_ROI_1 = funcs.extract_data_roi(d4_image_1, ROI_image)
                data_ROI_2 = funcs.extract_data_roi(d4_image_2, ROI_image)
            
            # define classifier
            clf = svm.SVC(kernel='linear')
                
            # get rid of runs that do not have all the task labels
            data_decode_1, labels_decode_1, groups_decode_1, data_decode_2, labels_decode_2, groups_decode_2 = funcs.balance_cross_decoding_data(data_ROI_1, labels_all_1, groups_1, data_ROI_2, labels_all_2, groups_2, full_task_labs)
            
            if feat_select and min(data_ROI_1.shape[1], data_ROI_2.shape[1]) > feat_num: # if using feature selection and the number of features is above the pre-defined number of retained features
                fs = SelectKBest(score_func=f_classif, k=feat_num)
                data_decode_1 = fs.fit_transform(data_decode_1, labels_decode_1)
                data_decode_2 = fs.fit_transform(data_decode_2, labels_decode_2)
                # feat_mask = fs.get_support() # returning a boolean mask showing which feature(1 in this mask) was retained
            
            # if the number of runs below 3, or there is no data in this ROI, then skip the decoding since cross-validation is impossible
            if np.unique(groups_decode_1).size < 3 or np.unique(groups_decode_2).size < 3 or data_ROI_1.shape[1] == 0 or data_ROI_2.shape[1] == 0:
                
                print('problem!')
                results = results.append({
                    'subject' : subj, 
                    'decoder': decoder, 
                    'smoothing': file_substrings_1[1], 
                    'measure': file_substrings_1[2], 
                    'training': file_substrings_1[3],
                    'testing': file_substrings_2[3],
                    'atlas': atlas,
                    'roi': ROI,
                    'roi_size' : data_decode_1.shape[1], 
                    'n_folds': np.unique(groups_decode_1).size, 
                    'mean_accuracy': float("nan")}, ignore_index=True)
                
                if not oneSide:
                    results = results.append({
                        'subject' : subj, 
                        'decoder': decoder, 
                        'smoothing': file_substrings_1[1], 
                        'measure': file_substrings_1[2], 
                        'training': file_substrings_2[3],
                        'testing': file_substrings_1[3],
                        'atlas': atlas,
                        'roi': ROI,
                        'roi_size' : data_decode_1.shape[1], 
                        'n_folds': np.unique(groups_decode_1).size, 
                        'mean_accuracy': float("nan")}, ignore_index=True)
                
                continue
            
            # running leave-one-run out cross-decoding
            accuracies = funcs.logo_cross_decoding(data_decode_1, labels_decode_1, groups_decode_1, data_decode_2, labels_decode_2, groups_decode_2, clf, one_side = oneSide)            
            acc_average = np.mean(accuracies, axis=1)            
            results = results.append({
                'subject' : subj, 
                'decoder': decoder, 
                'smoothing': file_substrings_1[1], 
                'measure': file_substrings_1[2], 
                'training': file_substrings_1[3],
                'testing': file_substrings_2[3],
                'atlas': atlas,
                'roi': ROI,
                'roi_size' : data_decode_1.shape[1], 
                'n_folds': np.unique(groups_decode_1).size, 
                'mean_accuracy': acc_average[0]}, ignore_index=True)
            
            if not oneSide:
                results = results.append({
                    'subject' : subj, 
                    'decoder': decoder, 
                    'smoothing': file_substrings_1[1], 
                    'measure': file_substrings_1[2], 
                    'training': file_substrings_2[3],
                    'testing': file_substrings_1[3],
                    'atlas': atlas,
                    'roi': ROI,
                    'roi_size' : data_decode_1.shape[1], 
                    'n_folds': np.unique(groups_decode_1).size, 
                    'mean_accuracy': acc_average[1]}, ignore_index=True)
                
    print(f"finish cross-decoding for {bids_subj}")

t_elapsed = time.time() - t_start

#%% save the data frame, the numpy array, and the confusion matrix if applicable

data_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/independent_roi/conA_conB_defB/cross_decode'
result_df_name = 'cross_cti_decodeAcc_smthN_spmT_rois_FPN_Schaefer.csv'
results.to_csv(os.path.join(data_dir, result_df_name))

#%% running one sample t test on the decoding accuracy across participants for each parcel
chance_level = 0.3333 # 0.1111 if 9 tasks, 0.3333 if 3 tasks
# accuracy_results = np.load(os.path.join(data_dir, 'decodeAcc_smthN_spmT_All_short_stim_SchPar.npy'))
accuracies_group = np.nanmean(accuracy_results,axis=1)
t_statistics, p_values = stats.ttest_1samp(a=accuracy_results, popmean=chance_level, axis=1, nan_policy='omit', alternative='greater')

#%% using FDR to control for multiple comparisons(all the parcels)
sig_level = 0.05

p_values_RG_fdr = multi.fdrcorrection(p_values[0,:])[1]
p_values_TF_fdr = multi.fdrcorrection(p_values[1,:])[1]

p_values_fdr = np.stack([p_values_RG_fdr, p_values_TF_fdr], axis=0)

#%% create Nifti images

accuracies_RG_group_img_data = np.zeros(atlas_image_map.shape)
accuracies_TF_group_img_data = np.zeros(atlas_image_map.shape)

for parcel_list_idx, parcel_idx in enumerate(atlas_idx_decode):
    accuracies_RG_group_img_data[atlas_image_map == parcel_idx] =  accuracies_group[0,parcel_list_idx]
    accuracies_TF_group_img_data[atlas_image_map == parcel_idx] =  accuracies_group[1,parcel_list_idx]

accuracies_RG_group_img= nib.Nifti1Image(accuracies_RG_group_img_data, atlas_image.affine, atlas_image.header)
plotting.plot_glass_brain(accuracies_RG_group_img,threshold=0.40)
nib.save(accuracies_RG_group_img, os.path.join(data_dir, 'decodeAcc_smthN_spmT_All_short_stim_SchPar.nii'))

accuracies_TF_group_img= nib.Nifti1Image(accuracies_TF_group_img_data, atlas_image.affine, atlas_image.header)
plotting.plot_glass_brain(accuracies_TF_group_img, threshold=0.13)
nib.save(accuracies_TF_group_img, os.path.join(data_dir, 'decodeAcc_smthN_spmT_TF_all_c1_SchPar.nii'))

# only for significant regions
sig_parcels_RG = atlas_idx_decode[p_values_RG_fdr < sig_level] # the significant parcel indexs
sig_parcels_TF = atlas_idx_decode[p_values_TF_fdr < sig_level]

sig_parcels_RG_names = atlas_labels[1:][p_values_RG_fdr < sig_level] # the significant parcel names
sig_parcels_TF_names = atlas_labels[1:][p_values_TF_fdr < sig_level]

sig_acc_RG = accuracies_group[0,:][p_values_RG_fdr < sig_level] # the decoding accuracy values of these significant parcels
sig_acc_TF = accuracies_group[1,:][p_values_TF_fdr < sig_level]

sig_p_RG = p_values_RG_fdr[p_values_RG_fdr < sig_level] # the corrected p value of these significant parcels
sig_p_TF = p_values_TF_fdr[p_values_TF_fdr < sig_level]

# create data frames that save the information of significant ROIs
regular_sig_ROIs_data ={'index':sig_parcels_RG,
                      'name':sig_parcels_RG_names,
                      'decode_acc':sig_acc_RG,
                      'p_fdr':sig_p_RG}

regular_sig_ROIs_df = pd.DataFrame(regular_sig_ROIs_data)
regular_sig_ROIs_name = 'sigROIs_RG.csv'
regular_sig_ROIs_df.to_csv(os.path.join(data_dir, regular_sig_ROIs_name))


transform_sig_ROIs_data ={'index':sig_parcels_TF,
                      'name':sig_parcels_TF_names,
                      'decode_acc':sig_acc_TF,
                      'p_fdr':sig_p_TF}

transform_sig_ROIs_df = pd.DataFrame(transform_sig_ROIs_data)
transform_sig_ROIs_name = 'sigROIs_TF.csv'
transform_sig_ROIs_df.to_csv(os.path.join(data_dir, transform_sig_ROIs_name))

## creating images that show the decoding acc for ONLY significant ROIs for regular condition

sig_accuracies_RG_group_img_data = np.zeros(atlas_image_map.shape)

for list_idx, parcel_idx_RG in enumerate(sig_parcels_RG):
    sig_accuracies_RG_group_img_data[atlas_image_map == parcel_idx_RG] = sig_acc_RG[list_idx] # accuracies_group[0,parcel_idx_RG-1]
    
sig_accuracies_RG_group_img= nib.Nifti1Image(sig_accuracies_RG_group_img_data, atlas_image.affine, atlas_image.header)
plotting.plot_glass_brain(sig_accuracies_RG_group_img)
nib.save(sig_accuracies_RG_group_img, os.path.join(data_dir, 'decodeAcc_smthN_spmT_All_short_stim_SchPar_Sig.nii'))

## creating images that show the decoding acc for ONLY significant ROIs for transform condition

sig_accuracies_TF_group_img_data = np.zeros(atlas_image_map.shape)

for list_idx, parcel_idx_TF in enumerate(sig_parcels_TF):
    sig_accuracies_TF_group_img_data[atlas_image_map == parcel_idx_TF] = sig_acc_TF[list_idx]
    
sig_accuracies_TF_group_img= nib.Nifti1Image(sig_accuracies_TF_group_img_data, atlas_image.affine, atlas_image.header)
plotting.plot_glass_brain(sig_accuracies_TF_group_img)
nib.save(sig_accuracies_TF_group_img, os.path.join(data_dir, 'decodeAcc_smthN_spmT_TF_all_c1_SchPar_Sig.nii'))

#%% check data
root_folder = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/atlas/Schaefer_2018'
smoth_folder = 'smth_8'
condition_folder = 'Block_all_c1'
file_name = 'decodeAcc_smth8_spmT_RG_TF_all_c1_SchPar.npy'

accuracy_results = np.load(os.path.join(root_folder, smoth_folder, condition_folder, file_name))

#%% check out the confusion matrix

# load confusion matrix
big_confusion_mat_file = confusion_mats_name + '.npy'
big_confusion_mat = np.load(os.path.join(data_dir, big_confusion_mat_file))

# select ONLY significant regions
sig_ROIs_df_file = ''
sig_ROIs_df_path = ''
sig_ROIs_df = pd.read_csv(sig_ROIs_df_path)
sig_ROIs_idex = np.array(sig_ROIs_df.loc[:,"index"]) - 1 # the index used to extract the data
sig_confusion_mat = big_confusion_mat[:,:,:,sig_ROIs_idex,:]

# averaging across subjects
sig_confusion_mat_sq = np.squeeze(sig_confusion_mat)
sig_confusion_mean = np.nanmean(sig_confusion_mat_sq, axis=2) # average across subjects
acc_bytask = np.diagonal(sig_confusion_mean, axis1=0, axis2=1) 

mean_acc_animal = acc_bytask[:, [0,1,2]].mean()
mean_acc_place = acc_bytask[:, [3,4,5]].mean()
mean_acc_vehicle = acc_bytask[:, [6,7,8]].mean()

mean_acc_age = acc_bytask[:, [0,3,6]].mean()
mean_acc_size = acc_bytask[:, [1,4,7]].mean()
mean_acc_location = acc_bytask[:, [2,5,8]].mean()

# plotting the result
import seaborn as sns

sns.heatmap(acc_bytask, vmin = 0.1111, vmax=0.2, annot=True,
            xticklabels = ['a|a', 'a|s', 'a|l', 'p|a', 'p|s','p|l','v|a','v|s','v|l'])


sns.heatmap(big_confusion_mean[:,:,2], vmin = 0.1111, vmax=0.2, annot=True,
            xticklabels=['a|a', 'a|s', 'a|l', 'p|a', 'p|s','p|l','v|a','v|s','v|l'],
            yticklabels=['a|a', 'a|s', 'a|l', 'p|a', 'p|s','p|l','v|a','v|s','v|l'])

#%% compare p values
p_value_smthN_RG_short_conj = np.squeeze(p_values) # 40 sig parcels w/o correction
p_value_smth8_RG_short_conj = np.squeeze(p_values) # 28 sig parcels w/o correction

np.sum(p_values < sig_level)

#%% plotting

ROI_image = image.load_img('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI/GLM-02_con-5/Clus_L_IFG.nii')
plotting.plot_glass_brain(ROI_image)

