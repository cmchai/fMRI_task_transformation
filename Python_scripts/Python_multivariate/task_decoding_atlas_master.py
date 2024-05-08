#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:11:31 2024

Master Script to run decoding on ROIs defined by an atlas

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
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import sklearn.svm as svm
# from sklearn.metrics import accuracy_score
# from sklearn.metrics import confusion_matrix
# from sklearn.model_selection import LeaveOneGroupOut

from scipy import stats
import statsmodels.stats.multitest as multi

#%% Define Atlas of interest

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

# define the index of the parcels that we will perform the decoding on
# - smth8 RG_all_c1 significant parcels [58, 143, 149, 187, 188]

atlas_idx_decode = np.arange(1, len(atlas_labels)) # all parcels
# atlas_idx_decode = np.array([58,143,149,187,188], dtype=int)

#%% define what data suffix to be decoded(should be in 4D format)
# one of following:
# 1. 'ALL-short_4D.nii', The short CTI trial regressors of all the block types(RG and TF)
# 2. 'RG-all-c1_4D.nii', The first CTI regressor of long and short trials in ONLY regular blocks, need regroup in this case
# stim dimension : ['RG-short-stim_4D.nii', 'TF-short-stim_4D.nii', 'All-short-stim_4D.nii', 'RG-all-stim-c1_4D.nii']
# rule dimension : ['RG-short-rule_4D.nii', 'TF-short-rule_4D.nii', 'All-short-rule_4D.nii', 'RG-all-rule-c1_4D.nii']

d4_files_suffix = 'RG-all-rule-c1_4D.nii'
n_d4_files = 1
need_regroup = True   # specify if reassigning runs into folds are necessary(from 4 runs into 8 folds)

smooth = True
compositional = True

if smooth:
    smth_suffix = '-B'
else:
    smth_suffix = ''
    
if compositional:
    comp_suffix = '-Comp'
else:
    comp_suffix = ''

glm_folder = 'GLM-02M' + comp_suffix + smth_suffix

#%% Define all the subjects and the GLM results to run

subjs_all = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50)))
subjs = subjs_all.copy()
# subjs = np.setdiff1d(subjs_all, subjs_run)

FLM_root_dir = os.path.join('/Volumes/extdrive/Task_Transform_GLM', glm_folder,'results') # GLM-02M:unsmoothed, GLM-02M-B:8mm smoothed

#%% Define tasks to be decoded
# can be one of the following:
# 1. np.arange(1,10), 9 conjunctive tasks
# 2. np.array(['animal', 'place', 'vehicle']), 3 compositional tasks in the stimulus type dimension
# 3. np.array(['age', 'location', 'size']), 3 compositional tasks in the task rule dimension

full_task_labs = np.array(['age', 'location', 'size'])

#%% Define decoders and if confusion matrix needed

decoders = ['SVM']
confusion = True

if confusion:
    big_confusion_mat = np.zeros((full_task_labs.size, full_task_labs.size, subjs.size, atlas_idx_decode.size, n_d4_files))    

#%% Creating the dataframe to store the decoding results

columns = ['subject', 'decoder', 'smoothing', 'measure', 'condition', 'atlas', 'parcel_index',
           'parcel', 'parcel_size', 'n_folds', 'mean_accuracy']

results = pd.DataFrame(columns = columns)

#%% Creating the np array to store the decoding results

accuracy_results = np.zeros((n_d4_files, subjs.size, atlas_idx_decode.size)) # here assume no multiple decoders are used, otherwise need to add another dimension

#%% Performing the decoding across ROIs
t_start = time.time()

    
for subj_idx, subj in enumerate(subjs): # loop over subjects

    bids_subj = 'sub-00' + str(subj)
    subj_dir = os.path.join(FLM_root_dir, bids_subj)
    d4_files = funcs.search_files(subj_dir, bids_subj, d4_files_suffix) # all the relevant 4d data
    
    print(f"start decoding for {bids_subj}")

    for d4_index, d4_file in enumerate(d4_files): # loop over different kind of 4d data
        file_substrings = d4_file.split('_')      # a list of substrings indicating all the information regarding the data
        
        d4_path = os.path.join(subj_dir, d4_file)
        d4_image = image.load_img(d4_path)
        d4_data = d4_image.get_fdata()
        
        labels_file = d4_file.replace('4D.nii', 'labels.txt')
        labels_path = os.path.join(subj_dir, labels_file)
        labels_df = pd.read_csv(labels_path, sep=",")

        # retrieve the task identities
        labels_all = labels_df['tasks'].to_numpy().copy()
        if not funcs.is_numerical(labels_all): # if the label are not numerical(e.g.1-9), but strings('animal', 'size'...)
            labels_all = labels_all.astype(str)        
        
        # retrieve the group identities
        groups = labels_df['runs'].to_numpy().copy()
        if need_regroup:
            groups = funcs.from_4runs_to_8groups(labels_all)
        
        for parcel_list_idx, parcel_idx in enumerate(atlas_idx_decode):
            
            parcel_label = atlas_labels[parcel_idx]              
            data_parcel = funcs.extract_data_roi_byatlas(d4_data, atlas_image_map, parcel_idx)
            
            for decoder_idx, decoder in enumerate(decoders): # loop over different kind of decoders
                
                if decoder == 'LDA':
                    clf = LinearDiscriminantAnalysis(solver='lsqr', shrinkage='auto')
                elif decoder == 'SVM':
                    clf = svm.SVC(kernel='linear')
                else:
                    sys.exit('decoder not properly defined!')
                    
                # get rid of runs that do not have all the task labels
                data_decode, labels_decode, groups_decode = funcs.del_group_misslabel(data_parcel, labels_all, groups, full_task_labs)
                
                # if the number of runs below 3, then skip the decoding since cross-validation is impossible
                if np.unique(groups_decode).size < 3 or data_parcel.shape[1] == 0:
                    
                    results = results.append({
                        'subject' : subj, 
                        'decoder': decoder, 
                        'smoothing': file_substrings[1], 
                        'measure': file_substrings[2], 
                        'condition': file_substrings[3],
                        'atlas': atlas_name,
                        'parcel_index': parcel_idx,
                        'parcel' : parcel_label, 
                        'parcel_size' : data_parcel.shape[1], 
                        'n_folds': np.unique(groups_decode).size, 
                        'mean_accuracy': float("nan")}, ignore_index=True)
                    
                    accuracy_results[d4_index, subj_idx, parcel_list_idx] = float("nan")
                    
                    continue
                
                # running leave-one-run out decoding
                if confusion:
                    accuracies, conf_mats = funcs.logo_decoding(data_decode, labels_decode, groups_decode, clf, confusion=True)
                    
                    # normalized the confusion matrix
                    conf_mat_sum = np.sum(conf_mats, axis=2) # sum over folds
                    n_test = np.sum(conf_mat_sum, axis=1) # the number of true categoris for each category
                    conf_mat_norm = conf_mat_sum/(n_test[:, np.newaxis]) # normalized confusion matrix
                    
                    # append the confusion matrix
                    big_confusion_mat[:,:,subj_idx, parcel_list_idx, d4_index] = conf_mat_norm               
                else:
                    accuracies = funcs.logo_decoding(data_decode, labels_decode, groups_decode, clf, confusion=False)
                
                acc_average = np.mean(accuracies)
                
                # saving the accuracy results(assume no separate decoders were used)
                accuracy_results[d4_index, subj_idx, parcel_list_idx] = acc_average
                
                results = results.append({
                    'subject' : subj, 
                    'decoder': decoder, 
                    'smoothing': file_substrings[1], 
                    'measure': file_substrings[2],
                    'condition': file_substrings[3],
                    'atlas': atlas_name,
                    'parcel_index': parcel_idx,
                    'parcel' : parcel_label, 
                    'parcel_size' : data_parcel.shape[1], 
                    'n_folds': accuracies.size, 
                    'mean_accuracy': acc_average}, ignore_index=True)
                
    print(f"finish decoding for {bids_subj}")

t_elapsed = time.time() - t_start

#%% save intermediate results
# accuracy_results_234 = accuracy_results.copy()
# accuracy_results_part1 = np.concatenate((accuracy_results_234, accuracy_results), axis=1)

#%% save the data frame, the numpy array, and the confusion matrix if applicable

data_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/compos/atlas/Schaefer_2018/smth8/task_rule/RG_all_c1_rule'

result_df_name = 'decodeAcc_smth8_spmT_RG_all_c1_rule_SchPar.csv'
results.to_csv(os.path.join(data_dir, result_df_name))

results_np_name = 'decodeAcc_smth8_spmT_RG_all_c1_rule_SchPar'
np.save(os.path.join(data_dir, results_np_name), accuracy_results)

if confusion:
    confusion_mats_name = 'conf_mat_decodeAcc_smth8_spmT_RG_all_c1_rule_SchPar'
    np.save(os.path.join(data_dir, confusion_mats_name), big_confusion_mat)

#%% running one sample t test on the decoding accuracy across participants for each parcel
chance_level = 0.3333 # 0.1111 if 9 tasks, 0.3333 if 3 tasks
accuracy_results = np.load(os.path.join(data_dir, 'decodeAcc_smthN_spmT_All_short_stim_SchPar.npy'))
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




