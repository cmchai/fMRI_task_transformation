#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:15:23 2024

Functions that used for decoding analysis

@author: mengqiao
"""

import os
import glob
import sys

import numpy as np
import numpy_indexed as npi
# import pandas as pd

from nilearn import image
from nilearn import masking
# import nibabel as nib

# from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
# import sklearn.svm as svm
from sklearn.preprocessing import LabelEncoder
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
    file_names.sort() # sorting the list in ascending order based on the letter
    
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
    conf_mats = np.zeros((n_labels, n_labels, n_groups))
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
            conf_mats[:,:,idx] = conf_mat
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

def extract_data_roi_byatlas(d4_data, atlas_map, ROI_idx):
    """
    Function to extract data ready for decoding from a certain ROI derived from an atlas

    Parameters
    ----------
    d4_data : 4d np array, the 4th dimension are the samples
    atlas_map : atlas image map in 3d np array, with different integers denoting different ROIs
    ROI_idx: the index of ROI that you want to extract data from
    
    the d4_data and atlas_map should have the same dimension. if not, please resample the atlas image first
    
    Returns
    -------
    data_roi : the ROI data ready for decoding, which is a 2-d np array [samples * features within the ROI]
    """    
        
    ROI_mask_array = atlas_map == ROI_idx # return a boolean 3-d array
    data_roi = np.transpose(d4_data[ROI_mask_array])
    
    # get rid of voxels that only contains zero values across samples
    non_zero_columns = np.any(data_roi != 0, axis=0)
    data_roi = data_roi[:,non_zero_columns]
    
    return data_roi


def extract_data_combined_roi_byatlas(d4_data, atlas_map, ROI_idxes):
    """
    Function to extract data ready for decoding from a group of certain ROIs (merged ROI) derived from an atlas

    Parameters
    ----------
    d4_data : 4d np array, the 4th dimension are the samples
    atlas_map : atlas image map in 3d np array, with different integers denoting different ROIs
    ROI_idxes: np array including indexes of ROIs that you want to extract data from as a merged (combined) region of interest
    
    the d4_data and atlas_map should have the same dimension. if not, please resample the atlas image first
    
    Returns
    -------
    data_roi : the ROI data ready for decoding, which is a 2-d np array [samples * features within the ROI]
    """    
        
    ROI_mask_array = np.isin(atlas_map, ROI_idxes) # return a boolean 3-d array
    data_roi = np.transpose(d4_data[ROI_mask_array])
    
    # get rid of voxels that only contains zero values across samples
    non_zero_columns = np.any(data_roi != 0, axis=0)
    data_roi = data_roi[:,non_zero_columns]
    
    return data_roi


def unique_preserve_order(groups):
    """
    Function to return the unique elements of a vector in the orginal order without reordering them 
    
    Parameters
    ----------
    groups : 1-d np array containing group identity of all samples

    Returns
    -------
    unique elements from this vector

    """
    unique_groups, index = np.unique(groups, return_index=True)
    sorted_index = np.argsort(index)
    return unique_groups[sorted_index]


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
    
    groups_idx = unique_preserve_order(groups)
    labels_bygroup = [labels[groups == group] for group in groups_idx]    # return a list of arrays
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

def is_numerical(array):
    try:
        # Try converting the array to a numerical type
        np.asarray(array, dtype=np.float64)
        return True
    except ValueError:
        return False
    

def from_4runs_to_8groups(labels):
    """
    the function of reassigning data from 4 runs into 8 folds for cross-validation
    mainly used whe both short and long trials c1 regressors were modeled  

    Parameters
    ----------
    labels : a numpy array with task ascendingly sorted for each fold
    
    Returns
    -------
    groups : a numpy array with corresponding group id(from 1 to 8)

    """
    groups = np.zeros(labels.size, dtype=int)
    
    if not is_numerical(labels):
        le = LabelEncoder()
        le.fit(labels)
        labels = le.transform(labels)
        
    labels_diffs = np.concatenate(([1],np.diff(labels))) # difference between neigboring elements of the labels array
    
    group = 1 # starting value of assinging group
    for label_idx, labels_diff in enumerate(labels_diffs):
        if labels_diff > 0:
            groups[label_idx] = group
        else:
            group += 1
            groups[label_idx] = group
            
    return groups

def split_groups(df):
    df['new_group'] = df.groupby(['runs', 'folds'], sort = False).ngroup() # the new group starts with 0
    return df['new_group'].to_numpy()