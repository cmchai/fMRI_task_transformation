#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:11:31 2024

Master Script to run decoding on ROIs defined by an atlas

@author: mengqiao
"""

import os
import sys
import glob

cwd = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts/Python_multivariate'
sys.path.append(cwd)

import funcs

import numpy as np
# import numpy_indexed as npi
import pandas as pd

from nilearn import image
# from nilearn import masking
# from nilearn import datasets
# import nibabel as nib

# from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import sklearn.svm as svm
# from sklearn.metrics import accuracy_score
# from sklearn.metrics import confusion_matrix

# from sklearn.model_selection import LeaveOneGroupOut

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
atlas_labels = atlas_labels_df[1].tolist().copy()
atlas_labels.insert(0, "Background")

# define the index of the parcels that we will perform the decoding on 
atlas_idx_decode = np.arange(1, len(atlas_labels))
atlas_idx_decode = atlas_idx_decode[[2,5,8]]

#%% define what data suffix to be decoded(should be in 4D format)

d4_files_suffix = '-all-c1_4D.nii'
n_d4_files = 2
need_regroup = True  # specify if reassigning runs into folds are necessary(from 4 runs into 8 folds)

#%% Define all the subjects and the GLM results to run

subjs_all = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50)))
subjs = subjs_all[[2,5,33]]
FLM_root_dir = '/Volumes/extdrive/Task_Transform_GLM/GLM-02M-B/results'

#%% Define decoders

decoders = ['SVM']

#%% Define tasks

full_task_labs = np.arange(1,10)

#%% Creating the dataframe to store the decoding results

columns = ['subject', 'decoder', 'smoothing', 'measure', 'condition', 'atlas', 'parcel_index',
           'parcel', 'parcel_size', 'n_folds', 'mean_accuracy']

results = pd.DataFrame(columns = columns)
results_name = ''

#%% Creating the dataframe to store the decoding results

accuracy_results = np.zeros((n_d4_files, subjs.size, atlas_idx_decode.size)) # here assum no multiple decoders are used, otherwise need to add another dimension

#%% Performing the decoding across ROIs 

for parcel_list_idx, parcel_idx in enumerate(atlas_idx_decode):
    
    parcel_label = atlas_labels[parcel_idx]
    
    for subj_idx, subj in enumerate(subjs): # loop over subjects
    
        bids_subj = 'sub-00' + str(subj)
        subj_dir = os.path.join(FLM_root_dir, bids_subj)
        d4_files = funcs.search_files(subj_dir, bids_subj, d4_files_suffix) # all the relevant 4d data
    
        for d4_index, d4_file in enumerate(d4_files): # loop over different kind of 4d data
            file_substrings = d4_file.split('_')      # a list of substrings indicating all the information regarding the data
            
            d4_path = os.path.join(subj_dir, d4_file)
            d4_image = image.load_img(d4_path)
            data_parcel = funcs.extract_data_roi_byatlas(d4_image, atlas_image, parcel_idx)
            
            labels_file = d4_file.replace('4D.nii', 'labels.txt')
            labels_path = os.path.join(subj_dir, labels_file)
            labels_df = pd.read_csv(labels_path, sep=",")

            labels_all = labels_df['tasks'].to_numpy().copy()
            groups = labels_df['runs'].to_numpy().copy()
            if need_regroup:
                groups = funcs.from_4runs_to_8groups(labels_all)
            
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
                if np.unique(groups_decode).size < 3:
                    
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
                    
                    continue
                
                # running leave-one-run out decoding
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

#%% save the data frame as well as the numpy array

data_dir = '/Users/mengqiao/Documents/fMRI_task_transform/Results'

result_df_name = ''
results.to_csv(os.path.join(data_dir, result_df_name))

results_np_name = ''
np.save(os.path.join(data_dir, results_np_name), accuracy_results)


#%% running one sample t test on the decoding accuracy across participants for each parcel 




#%% using FDR to control for multiple comparisons(all the parcels)