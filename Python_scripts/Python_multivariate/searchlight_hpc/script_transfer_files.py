#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 11:35:41 2024

@author: mengqiao
"""

import os
import paramiko
import numpy as np

#%% functions

def generate_source_folder_list(base_path, prefix, numbers):
    folders = [os.path.join(base_path, f"{prefix}{num}") for num in numbers]
    for folder in folders:
        if not os.path.exists(folder):
            print(f"Source folder '{folder}' does not exist.")
            folders.remove(folder)  # if the folder does not exsit, then remove the corresponding element from the folders list
    return folders

def search_files_by_suffixes(folder_path, suffixes):
    matching_files = []

    # Iterate through files in the specified folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        # Check if the file name ends with any of the specified suffixes
        for suffix in suffixes:
            if file_name.endswith(suffix):
                matching_files.append(file_path)
                # No need to check other suffixes if this file matches
                break

    return matching_files

#%% define the subj folder that need to be uploaded on HPC for all participants

# Define root folder and all the sub-folders for both source and destination
source_root_path = "/Volumes/extdrive/Task_Transform_GLM/GLM-02M/results"
hpc_root_path = "/kyukon/scratch/gent/vo/001/gvo00170/vsc43896/task_transform_fmri/FLM_4d_data/GLM-02M"

subj_prefix = "sub-00"
subj_all_numbers = np.concatenate((np.arange(2,6), np.arange(7,12), np.arange(13,18), np.arange(20,45), np.arange(46,50))).tolist()
subj_numbers = subj_all_numbers[4:]

source_folders = generate_source_folder_list(source_root_path, subj_prefix, subj_numbers)

#%% define the file's suffixes that we want to move

suffixes_to_copy = ["_4D.nii", "_labels.txt", "mask.nii", "-crop.nii"]

#%% Getting access to the HPC
ssh_client = paramiko.SSHClient()
ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

# Connect to HPC cluster
hostname = 'login.hpc.ugent.be' # 'gligar08.gastly.os'
port = 22
username = 'vsc43896'
password = 'Qqqnn467952.'
private_key_path = '/Users/mengqiao/.ssh/id_rsa_vsc'
private_key = paramiko.RSAKey.from_private_key_file(private_key_path, password=password)

# ssh_client.connect(hostname, username=username, password=password)
ssh_client.connect(hostname, username=username, pkey=private_key)

# Create SFTP client
sftp = ssh_client.open_sftp()

#%% upload the files on HPC

for source_folder in source_folders:
    
    # define the to-be uploaded files based on the suffixes
    
    local_files = search_files_by_suffixes(source_folder, suffixes_to_copy)
    hpc_folder = os.path.join(hpc_root_path, os.path.basename(source_folder))
    
    # Check if the destination directory exists on HPC
    try:
        sftp.chdir(hpc_folder)
    except FileNotFoundError:
        # Destination directory doesn't exist, create it
        sftp.mkdir(hpc_folder)
        sftp.chdir(hpc_folder)
        print(f"Created directory: {hpc_folder}")
        
    # Upload files
    for local_file in local_files:
        # Extract filename from the full path
        filename = os.path.basename(local_file)
    
        # Upload file to the HPC cluster
        sftp.put(local_file, filename)
    
        print(f"Uploaded {filename} to {hpc_folder}")

#%% Close SFTP and SSH connection
sftp.close()
ssh_client.close()

print("All files uploaded successfully.")