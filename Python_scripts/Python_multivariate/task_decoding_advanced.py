#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 12:14:51 2024

Task Decoding of Task Transformation paradigm

Advanced script

@author: mengqiao
"""

import os
import glob
import sys
import pprint

import numpy as np
import numpy_indexed as npi
import pandas as pd

from nilearn import image
from nilearn import plotting, surface
from nilearn import masking
from nilearn import datasets
import nibabel as nib

from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import sklearn.svm as svm
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix

from sklearn.model_selection import LeaveOneGroupOut

#%% Functions

def search_files(folder_path, start_substring, end_substring):
    """
    Function to return files' names that start and end with particular substrings
    Parameters
    ----------
    folder_path : the full path of the folder of interest
    start_substring : the particular substring the file name should start with
    end_substring : the particular substring the file name should end with.

    Returns
    -------
    file_names : return a list of ONLY the names of all the relevant files
    """
    
    # Construct the pattern for file names
    pattern = os.path.join(folder_path, f"{start_substring}*{end_substring}")
    
    # Use glob to find files matching the pattern
    matching_files = glob.glob(pattern)
    
    # Extract file names from file paths
    file_names = [os.path.basename(file) for file in matching_files]
    
    return file_names


def decoding(train_set, test_set, train_labels, test_labels, clf, n_labels, confusion=False):
    """
    Function to train a decoder and test it using separate datasets
    ----------
    Parameters
    ----------
    train_set : The training data, which is a 2-d np array [samples * features]
    test_set : The testing data
    train_labels : The correct labels of the training data, which is a 1-d np array 
    test_labels : The correct labels of the testing data
    clf : the decoder object
    n_labels : number of unique labels
    confusion : return the confusion matrix or not. The default is False.

    Returns
    -------
    the decoding accuracy on the testing data, which is between 0 and 1
    the confusion matrix, which is a 2-d array.
    """
    
    uniq_labels = np.unique(train_labels)
    if uniq_labels.shape[0] != n_labels:
        sys.exit('the training set does not have all the labels!')
    
    clf.fit(train_set, train_labels)
    pred = clf.predict(test_set)
    accuracy = accuracy_score(test_labels, pred)
    
    if confusion:
        conf_mat = confusion_matrix(test_labels, pred, labels=uniq_labels)
        return accuracy, conf_mat
    else:
        return accuracy
    
    
def logo_decoding(data_all, labels_all, groups, clf, confusion=False):
    """
    Function to run Leave-One-Group-Out(logo) cross-validation decoding
    ----------
    Parameters
    ----------
    data_all : all the data, which is a 2-d np array [samples * features]
    labels_all : The correct labels of the data, which is a 1-d np array 
    groups : the group labels of all the samples, which is a 1-d np array 
    clf : the decoder object
    confusion : return the confusion matrix or not. The default is False.

    Returns
    -------
    the decoding accuracies of all the testing folds, which is a 1-d np array
    the confusion matrics, which is a list of appended 2-d array from each fold.

    """
    # function to perform leave-one-group-out(logo) decoding
    logo = LeaveOneGroupOut()
    folds = logo.split(data_all, labels_all, groups=groups)
    n_groups = np.unique(groups).shape[0]
    n_labels = np.unique(labels_all).shape[0]
    
    accuracies = np.zeros(n_groups)
    conf_mats = []
    for idx, fold in enumerate(folds):
        train_idx, test_idx = fold
        # print(train_idx, test_idx)
        
        data_train = data_all[train_idx]
        labels_train = labels_all[train_idx]
        
        data_test = data_all[test_idx]
        labels_test = labels_all[test_idx]
        # print(S_test)
        
        if confusion:
            accuracy, conf_mat = decoding(data_train, data_test, labels_train, labels_test, clf, n_labels, confusion=True)
            accuracies[idx] = accuracy
            conf_mats.append(conf_mat)
        else:
            accuracy = decoding(data_train, data_test, labels_train, labels_test, clf, n_labels, confusion=False)
            accuracies[idx] = accuracy
            
    if confusion:
        return accuracies, conf_mats
    else:
        return accuracies
    

def extract_data_roi(d4_image, roi_image):
    """
    Function to extract data ready for decoding from a certain ROI

    Parameters
    ----------
    d4_image : 4d nifti image, the 4th dimension are the samples
    roi_image : ROI mask nifti image
    
    Returns
    -------
    data_roi : the ROI data ready for decoding, which is a 2-d np array [samples * features within the ROI]
    """
    
    if d4_image.shape[0:3] != roi_image.shape: # resample the ROI image to match the 4d image
        roi_image = image.resample_to_img(roi_image, d4_image, interpolation='nearest')
           
    data_roi = masking.apply_mask(d4_image, roi_image)
    
    # get rid of voxels that only contains zero values across samples
    non_zero_columns = np.any(data_roi != 0, axis=0)
    data_roi = data_roi[:,non_zero_columns]
    
    return data_roi

def extract_data_roi_byatlas(d4_image, atlas_image, ROI_idx):
    """
    Function to extract data ready for decoding from a certain ROI derived from an atlas

    Parameters
    ----------
    d4_image : 4d nifti image, the 4th dimension are the samples
    atlas_image : atlas map image, with different integers denoting different ROIs
    ROI_idx: the index of ROI that you want to extract data from
    
    Returns
    -------
    data_roi : the ROI data ready for decoding, which is a 2-d np array [samples * features within the ROI]
    """
    
    if d4_image.shape[0:3] != atlas_image.shape: # resample the atlas image to match the 4d image
        atlas_image = image.resample_to_img(atlas_image, d4_image, interpolation='nearest')
        
    atlas_array = atlas_image.get_fdata()
    ROI_mask_array = atlas_array == ROI_idx # return a boolean 3-d array
    
    d4_data = d4_image.get_fdata()
    data_roi = np.transpose(d4_data[ROI_mask_array])
    
    # get rid of voxels that only contains zero values across samples
    non_zero_columns = np.any(data_roi != 0, axis=0)
    data_roi = data_roi[:,non_zero_columns]
    
    return data_roi


def del_group_misslabel(data, labels, groups, labels_complete):
    """
    Function to delete data(and corresponding labels) from group(s) that does not have the data for each label

    Parameters
    ----------
    data : 2-d np array [samples * features]
    labels : 1-d np array containing the labels of all samples
    groups : 1-d np array containing group identity of all samples
    labels_complete : 1-d np array of complete label set

    Returns
    -------
    data, labels, and groups after deletion

    """
    
    groups_idx = np.unique(groups)
    labels_bygroup = npi.group_by(groups).split(labels) # return a list of arrays
    groups_bool = np.array([np.in1d(labels_complete, labels_1group).all() for labels_1group in labels_bygroup])
    
    if groups_bool.all():
        return data, labels, groups
    else:
        groups_include = groups_idx[groups_bool]
        samples_bool = np.in1d(groups, groups_include) # boolean array of which samples to include for decoding
        
        data_new = data.copy()
        labels_new = labels.copy()
        groups_new = groups.copy()
        
        data_new = data_new[samples_bool,:]
        labels_new = labels_new[samples_bool]
        groups_new = groups_new[samples_bool]
        
        return data_new, labels_new, groups_new    

def from_4runs_to_8groups(labels):
    """
    the function of reassigning data from 4 runs into 8 folds for cross-validation

    Parameters
    ----------
    labels : a numpy array with task ascendingly sorted for each fold
    
    Returns
    -------
    groups : a numpy array with corresponding group id(from 1 to 8)

    """
    groups = np.zeros(labels.size, dtype=int)
    labels_diffs = np.concatenate(([1],np.diff(labels))) # difference between neigboring elements of the labels array
    
    group = 1 # starting value of assinging group
    for label_idx, labels_diff in enumerate(labels_diffs):
        if labels_diff > 0:
            groups[label_idx] = group
        else:
            group += 1
            groups[label_idx] = group
            
    return groups        
    
#%% Creating the dataframe to store the decoding results

columns = ['subject', 'decoder', 'smoothing', 'measure', 'condition',
           'ROI', 'ROI_size', 'n_folds', 'mean_accuracy']

results = pd.DataFrame(columns = columns)

#%% Define ROIs of interest

ROIs_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI'
ROIs_scheme = 'GLM-02_con-5'
ROIs_dir = os.path.join(ROIs_root_dir, ROIs_scheme)
ROIs = [file for file in os.listdir(ROIs_dir) if file.endswith('.nii')]
plot_roi = False # plot ROIs or not

#%% Define all the subjects and the GLM results to run

subjs_all = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50)))
subjs_run = np.array([])
subjs = np.setdiff1d(subjs_all, subjs_run)
FLM_root_dir = '/Volumes/extdrive/Task_Transform_GLM/GLM-02M-B/results'

#%% Define decoders

decoders = ['LDA', 'SVM']

#%% Define full task labels

full_task_labs = np.arange(1,10)

#%% Whether regrouping is needed, true when split 4 runs into 8 groups for decoding

regrouping = True

#%% Running the decoding

for ROI_idx, ROI in enumerate(ROIs): # loop over ROIs
    
    ROI_path = os.path.join(ROIs_dir, ROI)
    ROI_image = image.load_img(ROI_path)
    
    if plot_roi:
        plotting.plot_roi(ROI_image,black_bg=True, title = ROI)
        # plotting.plot_glass_brain(ROI_image)
        
    
    for subj_idx, subj in enumerate(subjs): # loop over subjects       
        bids_subj = 'sub-00' + str(subj)
        subj_dir = os.path.join(FLM_root_dir, bids_subj)
        d4_files = search_files(subj_dir, bids_subj, 'TF-all-c1_4D.nii') # all the relevant 4d data
        
        for d4_index, d4_file in enumerate(d4_files): # loop over different kind of 4d data
            file_substrings = d4_file.split('_')      # a list of substrings indicating all the information regarding the data
            
            d4_path = os.path.join(subj_dir, d4_file)
            d4_image = image.load_img(d4_path)
            data_ROI = extract_data_roi(d4_image, ROI_image)
            
            labels_file = d4_file.replace('4D.nii', 'labels.txt')
            labels_path = os.path.join(subj_dir, labels_file)
            labels_df = pd.read_csv(labels_path, sep=",")

            labels_all = labels_df['tasks'].to_numpy().copy()
            groups = labels_df['runs'].to_numpy().copy()
            if regrouping:
                groups = from_4runs_to_8groups(labels_all)
            
            for decoder_idx, decoder in enumerate(decoders): # loop over different kind of decoders
                
                if decoder == 'LDA':
                    clf = LinearDiscriminantAnalysis(solver='lsqr', shrinkage='auto')
                elif decoder == 'SVM':
                    clf = svm.SVC(kernel='linear')
                else:
                    sys.exit('decoder not properly defined!')
                    
                # get rid of runs that do not have all the task labels
                data_decode, labels_decode, groups_decode = del_group_misslabel(data_ROI, labels_all, groups, full_task_labs)
                
                # if the number of runs below 3, then skip the decoding since cross-validation is impossible
                if np.unique(groups_decode).size < 3:
                    
                    results = results.append({
                        'subject' : subj, 
                        'decoder': decoder, 
                        'smoothing': file_substrings[1], 
                        'measure': file_substrings[2], 
                        'condition': file_substrings[3],
                        'ROI' : ROI, 
                        'ROI_size' : data_ROI.shape[1], 
                        'n_folds': np.unique(groups_decode).size, 
                        'mean_accuracy': float("nan")}, ignore_index=True)
                    
                    continue
                
                # running leave-one-run out decoding
                accuracies = logo_decoding(data_decode, labels_decode, groups_decode, clf, confusion=False)
                acc_average = np.mean(accuracies)
                
                results = results.append({
                    'subject' : subj, 
                    'decoder': decoder, 
                    'smoothing': file_substrings[1], 
                    'measure': file_substrings[2], 
                    'condition': file_substrings[3],
                    'ROI' : ROI, 
                    'ROI_size' : data_ROI.shape[1], 
                    'n_folds': accuracies.size, 
                    'mean_accuracy': acc_average}, ignore_index=True)
            
            
#%% Saving the data            
results.to_csv('/Users/mengqiao/Documents/fMRI_task_transform/Results/decoding_results_smth8_TF-all-c1.csv')      
            
#%% Using a brain atlas        
atlas_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas'    
atlas_name = 'schaefer_2018'
atlas_dir = os.path.join(atlas_root_dir, atlas_name)

brain_atlas_sch = datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=17, resolution_mm=2, data_dir=atlas_dir)

brain_atlas_sch_img = image.load_img(brain_atlas_sch['maps'])
brain_atlas_sch_data = brain_atlas_sch_img.get_fdata()
brain_atlas_sch_labels = [str(label, encoding='UTF-8') for label in brain_atlas_sch['labels']]
brain_atlas_sch_labels.insert(0, "Background")

# resampling the atlas
brain_atlas_sch_img_resamp = image.resample_to_img(brain_atlas_sch_img, big_4D, interpolation='nearest')
brain_atlas_sch_img_resamp_data = brain_atlas_sch_img_resamp.get_fdata()
np.sum(brain_atlas_sch_img_resamp_data == 0)/brain_atlas_sch_img_resamp_data.size

# save the resampled image
nib.nifti1.save(brain_atlas_sch_img_resamp, 'Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2.5mm.nii.gz')

# plotting the atlas
plotting.plot_roi(brain_atlas_sch_img_resamp, colorbar=True)

fsaverage = datasets.fetch_surf_fsaverage()
mesh = surface.load_surf_mesh(fsaverage.pial_right)
brain_atlas_sch_surface = surface.vol_to_surf(brain_atlas_sch_img, mesh)

plotting.plot_surf_roi(fsaverage['infl_right'], roi_map=brain_atlas_sch_surface,
                       hemi='right', view='lateral',
                       bg_map=fsaverage['sulc_right'], bg_on_data=True)

data111 = pd.read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas/schaefer_2018/schaefer_2018/Schaefer2018_400Parcels_17Networks_order.txt', delimiter = "\t", header = None)


