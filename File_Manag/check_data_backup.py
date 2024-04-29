#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 11:01:39 2023

Script that checks if all the fMRI subject folders on local machine are also uploaded on the share drive

@author: mengqiao
"""

import os

# Replace these with the paths to your two folders
local_path = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/dcm'
share_path = '/Volumes/mchai/shares/pp02_labbraem/Share/Task_Transform_fMRI/data/MRI/dcm'


#%% check if all the participants in the local folder are also present in the share drive

# List files in the first folder
files_in_local_path = set(os.listdir(local_path))

# List files in the second folder
files_in_share_path = set(os.listdir(share_path))

# Files present in local but not in share
files_only_in_local_path = files_in_local_path - files_in_share_path

# Files present in folder2 but not in folder1
files_only_in_share_path = files_in_share_path - files_in_local_path

# Files present in both folders
common_files = files_in_local_path.intersection(files_in_share_path)
print(len(common_files))

# Print the results or perform further actions as needed
print(f"Files only in local folder: {files_only_in_local_path}")
print(f"Files only in share drive: {files_only_in_share_path}")
print(f"Common files in both folders: {common_files}")

#%% check if the SIZE of all the participants's folder in the local folder are equivalent to the ones in the share drive

def get_folder_size(folder_path):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(folder_path):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            total_size += os.path.getsize(file_path)
    return total_size

# Get a list of folder names in local_path
folder_names1 = set(os.listdir(local_path))

# Get a list of folder names in share_path
folder_names2 = set(os.listdir(share_path))

# Find common folder names between the two parent directories
common_folder_names = folder_names1.intersection(folder_names2)

for folder_name in common_folder_names:
    folder_path1 = os.path.join(local_path, folder_name)
    folder_path2 = os.path.join(share_path, folder_name)

    size1 = get_folder_size(folder_path1)
    size2 = get_folder_size(folder_path2)

    size_discrepancy = size1 - size2

    if size_discrepancy > 0:
        print(f"{folder_name}: local_path is {size_discrepancy} bytes larger than share_path")
    elif size_discrepancy < 0:
        print(f"{folder_name}: share_path is {abs(size_discrepancy)} bytes larger than local_path")
    else:
        print(f"{folder_name}: Both folders have the same size")

# Check for folders that exist in one parent directory but not the other
for folder_name in folder_names1 - folder_names2:
    print(f"{folder_name} only exists in {local_path}")

for folder_name in folder_names2 - folder_names1:
    print(f"{folder_name} only exists in {share_path}")
    
    
    
#%% check if all the files in one folder are also present in another folder? this process could take QUITE SOME TIME !!!!

import os
import filecmp

# Replace these with the paths to the two root directories you want to compare
dir1 = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/dcm/sub-003'
dir2 = '/Volumes/mchai/shares/pp02_labbraem/Share/Task_Transform_fMRI/data/MRI/dcm/sub-003'

# Lists to store the results
files_only_in_dir1 = []
files_only_in_dir2 = []

def compare_directories(dir1, dir2):
    dir_comp = filecmp.dircmp(dir1, dir2)

    # Compare common subdirectories
    for common_dir in dir_comp.common_dirs:
        compare_subdirectories(os.path.join(dir1, common_dir), os.path.join(dir2, common_dir))

    # List files in dir1 but not in dir2
    files_only_in_dir1.extend(file for file in os.listdir(dir1) if file not in os.listdir(dir2))

    # List files in dir2 but not in dir1
    files_only_in_dir2.extend(file for file in os.listdir(dir2) if file not in os.listdir(dir1))

    # Recursively compare sub-subdirectories
    for common_subdir in dir_comp.common_dirs:
        compare_subdirectories(os.path.join(dir1, common_subdir), os.path.join(dir2, common_subdir))

def compare_subdirectories(subdir1, subdir2):
    dir_comp = filecmp.dircmp(subdir1, subdir2)

    # Compare files in the current subdirectory
    for common_file in dir_comp.common_files:
        file1 = os.path.join(subdir1, common_file)
        file2 = os.path.join(subdir2, common_file)

        if filecmp.cmp(file1, file2, shallow=False):
            pass
        else:
            print(f"Files in {common_file} are different.")

    # Recursively compare sub-subdirectories
    for common_subdir in dir_comp.common_dirs:
        compare_subdirectories(os.path.join(subdir1, common_subdir), os.path.join(subdir2, common_subdir))

# Start the comparison
compare_directories(dir1, dir2)

# Print the lists of files only in dir1 and dir2
print("Files only in dir1:")
for file in files_only_in_dir1:
    print(file)

print("Files only in dir2:")
for file in files_only_in_dir2:
    print(file)














