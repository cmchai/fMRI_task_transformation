#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 17:13:00 2024

Script to smooth images

@author: mengqiao
"""
#%% load necessary packages

import os
import glob
import nibabel as nib
from nilearn import image

#%% create a list of images to be smoothed

folder_path = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/conjunc/searchlight/smthN'

img_file_pattern = '*SL_acc.nii'
search_pattern = os.path.join(folder_path, img_file_pattern)
imgs_path = glob.glob(search_pattern)


#%% conduct the smoothing

fwhm = 8 # in mm

smoothed_imgs = image.smooth_img(imgs_path, fwhm)

#%% Saving images

subfolder = 'after_smth8'

for idx, smth_img in enumerate(smoothed_imgs):
    img_name = os.path.basename(imgs_path[idx])
    smth_img_name = img_name.replace('SL_acc', 'SL_acc_smth8')
    nib.save(smth_img, os.path.join(folder_path, subfolder, smth_img_name))

