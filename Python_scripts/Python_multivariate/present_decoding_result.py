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

def find_files_with_prefix_and_suffix_2(root_folder, prefix, suffix, midfix, nofix):
    file_paths = []
    for folder_name, _, filenames in os.walk(root_folder):
        for filename in filenames:
            if filename.startswith(prefix) and filename.endswith(suffix) and (midfix in filename) and (nofix not in filename):
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

root_folder = "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding"
prefix = "decodeAcc_smthN"
suffix = "SchPar.npy"
midfix = "short"
nofix = "4to8"

matching_files = find_files_with_prefix_and_suffix_2(root_folder, prefix, suffix, midfix, nofix)

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

comp_decoding_schemes2 = np.array([("stim" in decoding_scheme) or ("rule" in decoding_scheme) for decoding_scheme in decoding_schemes], dtype=bool)

#%% running one sample t test on the decoding accuracy across participants for each parcel

p_values = np.zeros((400, decoding_schemes2.shape[0]))

for scheme_ind, scheme in enumerate(decoding_schemes2):
    
    decoding_result = np.squeeze(decoding_results[:,:,scheme_ind])
    
    if comp_decoding_schemes2[scheme_ind]:
        chance_level = 0.3333
    else:
        chance_level = 0.1111
        
    t_statistic, p_value = stats.ttest_1samp(a=decoding_result, popmean=chance_level, axis=0, nan_policy='omit', alternative='greater')
    p_values[:,scheme_ind] = p_value.copy()

p_values_sigs = p_values < sig_level

#%% adjusting p value using FDR
sig_level = 0.05
p_values_fdr = np.zeros(p_values.shape)

for idx in range(p_values.shape[1]):
    p_fdr = multi.fdrcorrection(p_values[:,idx])[1]
    p_values_fdr[:,idx] = p_fdr.copy()

p_values_fdr_sigs = p_values_fdr < sig_level

#%% See which ROIs consistently show significant results

# results from 8-fold decoding
fold8_decoding_schemes2 = np.array([("ll_" in decoding_scheme) for decoding_scheme in decoding_schemes], dtype=bool)

fold8_decoding_all_sig_rois = np.sum(p_values_fdr_sigs[:,fold8_decoding_schemes2], axis=1) >= 3
fold8_decoding_all_sig_roi_idx = atlas_labels_idx[fold8_decoding_all_sig_rois]
fold8_decoding_all_sig_roi_labels = atlas_labels[fold8_decoding_all_sig_rois]

img_data = np.zeros(atlas_image_map.shape)
for parcel_list_idx, parcel_idx in enumerate(atlas_labels_idx):
    img_data[atlas_image_map == parcel_idx] =  fold8_decoding_all_sig_rois[parcel_list_idx]
img = nib.Nifti1Image(img_data, atlas_image.affine, atlas_image.header)
img_name = 'Sig_rois_8fold_SchPar.nii'
plotting.plot_glass_brain(img, threshold=0.001, title='sig ROIs for 8 fold decoding')

# results from 4-fold decoding (un-corrected result)
fold4_decoding_all_sig_rois = np.sum(p_values_sigs[:,~fold8_decoding_schemes2], axis=1) >= 4
fold4_decoding_all_sig_roi_idx = atlas_labels_idx[fold4_decoding_all_sig_rois]
fold4_decoding_all_sig_roi_labels = atlas_labels[fold4_decoding_all_sig_rois]

img_data2 = np.zeros(atlas_image_map.shape)
for parcel_list_idx, parcel_idx in enumerate(atlas_labels_idx):
    img_data2[atlas_image_map == parcel_idx] =  fold4_decoding_all_sig_rois[parcel_list_idx]
img2 = nib.Nifti1Image(img_data2, atlas_image.affine, atlas_image.header)
img_name2 = 'Sig_rois_4fold_SchPar.nii'
plotting.plot_glass_brain(img2, threshold=0.001, title='Un-corrected sig ROIs for 4 fold decoding')

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

#%% compare p-values of un-corrected results  
scheme_bool = np.array([has_substring(scheme, 'RG_short') for scheme in decoding_schemes])
short_schemes = decoding_schemes2[scheme_bool]

p_values_RGshort =  p_values[:,scheme_bool]
p_values_sigs = p_values_RGshort < 0.05
p_values_sigs_avg = np.sum(p_values_sigs, axis = 0)

decoding_results_RGshort_submean = decoding_results_submean[:,scheme_bool]
decoding_results_RGshort_submean_sig_uncor = decoding_results_RGshort_submean
decoding_results_RGshort_submean_sig_uncor[~p_values_sigs] = 0

for scheme_idx, decoding_scheme in enumerate(short_schemes):
    decoding_acc = decoding_results_RGshort_submean_sig_uncor[:,scheme_idx]
    
    img_data = np.zeros(atlas_image_map.shape)
    for parcel_list_idx, parcel_idx in enumerate(atlas_labels_idx):
        img_data[atlas_image_map == parcel_idx] =  decoding_acc[parcel_list_idx]
    img = nib.Nifti1Image(img_data, atlas_image.affine, atlas_image.header)
    img_name = 'decodeAcc_' + decoding_scheme + '_SchPar_Sig_Uncor.nii'
    plotting.plot_glass_brain(img, threshold=0.001, title='Uncorrected decoding acc for ' + decoding_scheme, vmin=0.33, vmax=0.46)
    nib.save(img, os.path.join(root_folder, 'summary', img_name))

#%% check conjunctive uncorrected result

conjunc_decod_acc = np.squeeze(np.load('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/conjunc/atlas/Schaefer_2018/smth8/RG_short/decodeAcc_smth8_spmT_RG_short_SchPar.npy'))
conjunc_decod_acc_submean = np.nanmean(conjunc_decod_acc, axis=0)

t_statistics, p_values = stats.ttest_1samp(a=conjunc_decod_acc, popmean=0.1111, axis=0, nan_policy='omit', alternative='greater')
np.sum(p_values < 0.05)

conjunc_decod_acc_submean_sig = conjunc_decod_acc_submean
conjunc_decod_acc_submean_sig[p_values >= 0.05] = 0

img_data = np.zeros(atlas_image_map.shape)
for parcel_list_idx, parcel_idx in enumerate(atlas_labels_idx):
    img_data[atlas_image_map == parcel_idx] =  conjunc_decod_acc_submean_sig[parcel_list_idx]
img = nib.Nifti1Image(img_data, atlas_image.affine, atlas_image.header)
plotting.plot_glass_brain(img, threshold=0.001, title='Uncorrected decoding acc for smth8 conjunc RG-short', vmin=0.11, vmax=0.15)

