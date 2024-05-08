#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:36:26 2024

checking the decoding results systematically

@author: mengqiao
"""

#%% import packages

import os
import sys
import glob
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as multi
import pandas as pd

from nilearn import image, plotting
import nibabel as nib

#%% functions

def find_files_with_prefix_and_suffix(root_folder, prefix, suffix):
    file_paths = []
    for folder_name, _, filenames in os.walk(root_folder):
        for filename in filenames:
            if filename.startswith(prefix) and filename.endswith(suffix):
                file_paths.append(os.path.join(folder_name, filename))
    return file_paths

def extract_substring_between(string, substring_a, substring_b):
    start_index = string.find(substring_a)
    if start_index == -1:
        return None  # Substring A not found
    end_index = string.find(substring_b, start_index + len(substring_a))
    if end_index == -1:
        return None  # Substring B not found after substring A
    return string[start_index + len(substring_a):end_index]

# Define a function to check if a substring exists in a value
def has_substring(value, substring):
    return substring in value

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
atlas_labels_idx = np.arange(1,atlas_labels.shape[0]+1)

#%% import all the decoding results

root_folder = "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/compos/atlas/Schaefer_2018"
prefix = "decodeAcc"
suffix = "SchPar.npy"

matching_files = find_files_with_prefix_and_suffix(root_folder, prefix, suffix)

decoding_results = np.zeros((43, 400, len(matching_files)))
decoding_schemes = []


for file_ind, file_path in enumerate(matching_files):
    print(file_path)
    file_name = os.path.basename(file_path)
    decoding_scheme = extract_substring_between(file_name, 'decodeAcc_', '_SchPar')
    decoding_schemes.append(decoding_scheme)
    
    decode_result = np.squeeze(np.load(file_path)) # should be 43 * 400
    decoding_results[:,:,file_ind] = decode_result
    
decoding_schemes2 = np.array(decoding_schemes, dtype=np.str_)

#%% running one sample t test on the decoding accuracy across participants for each parcel

chance_level = 0.3333 # 0.1111 if 9 tasks, 0.3333 if 3 tasks
t_statistics, p_values = stats.ttest_1samp(a=decoding_results, popmean=chance_level, axis=0, nan_policy='omit', alternative='greater')

#%% adjusting p value using FDR
sig_level = 0.05
p_values_fdr = np.zeros(p_values.shape)

for idx in range(p_values.shape[1]):
    p_fdr = multi.fdrcorrection(p_values[:,idx])[1]
    p_values_fdr[:,idx] = p_fdr

p_values_fdr_sigs = p_values_fdr < sig_level

#%% compare decoding accuracy of different decoding schemes

sig_parcels_byscheme = np.sum(p_values_fdr_sigs, axis=0)

sig_parcels_byscheme_df = pd.DataFrame({'decode_scheme':decoding_schemes2,
                                        'nr_sig_parccels':sig_parcels_byscheme})

sig_parcels_bydim = sig_parcels_byscheme_df.groupby(lambda x: has_substring(sig_parcels_byscheme_df['decode_scheme'].iloc[x], 'stim')).mean()
sig_parcels_bysmth = sig_parcels_byscheme_df.groupby(lambda x: has_substring(sig_parcels_byscheme_df['decode_scheme'].iloc[x], 'smthN')).mean()
sig_parcels_byfold = sig_parcels_byscheme_df.groupby(lambda x: has_substring(sig_parcels_byscheme_df['decode_scheme'].iloc[x], 'RG_short')).mean()
sig_parcels_8fold = sig_parcels_byscheme_df[~sig_parcels_byscheme_df['decode_scheme'].str.contains('RG_short')].groupby(lambda x: has_substring(sig_parcels_byscheme_df['decode_scheme'].iloc[x], 'RG')).mean()

#%% compare different parcels
sig_parcels = np.mean(p_values_fdr_sigs, axis=1)  # the proportion of finding a significant decodable regionsm across decoding schemes
sig_parcels_df = pd.DataFrame({'parcel_name':atlas_labels,
                               'parcel_index':atlas_labels_idx,
                               'mean_sig':sig_parcels})

sig_parcels_img_data = np.zeros(atlas_image_map.shape)

for parcel_list_idx, parcel_idx in enumerate(atlas_labels_idx):
    sig_parcels_img_data[atlas_image_map == parcel_idx] =  sig_parcels[parcel_list_idx]

sig_parcels_img = nib.Nifti1Image(sig_parcels_img_data, atlas_image.affine, atlas_image.header)
plotting.plot_glass_brain(sig_parcels_img, threshold=0)
nib.save(sig_parcels_img, os.path.join(root_folder, 'summary', 'sig_parcels.nii'))

sig_parcels_df_sorted = sig_parcels_df[sig_parcels_df['mean_sig'] > 0].sort_values('mean_sig', ascending = False)

#%% creating nifti images of significant parcels for each decoding scheme

decoding_results_submean = np.mean(decoding_results, axis = 0)

decoding_results_submean_sig = decoding_results_submean
decoding_results_submean_sig[~p_values_fdr_sigs] = 0

for scheme_idx, decoding_scheme in enumerate(decoding_schemes2):
    decoding_acc = decoding_results_submean_sig[:,scheme_idx]
    
    img_data = np.zeros(atlas_image_map.shape)
    for parcel_list_idx, parcel_idx in enumerate(atlas_labels_idx):
        img_data[atlas_image_map == parcel_idx] =  decoding_acc[parcel_list_idx]
    img = nib.Nifti1Image(img_data, atlas_image.affine, atlas_image.header)
    img_name = 'decodeAcc_' + decoding_scheme + '_SchPar_Sig.nii'
    plotting.plot_glass_brain(img, threshold=0.001, title='decoding acc for ' + decoding_scheme, vmin=0.33, vmax=0.46)
    nib.save(img, os.path.join(root_folder, 'summary', img_name))
    
        


