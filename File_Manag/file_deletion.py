#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 10:57:06 2024

Scripts that delete local files AFTER backing up them somewhere else

@author: mengqiao
"""
import os

#%% All the useful functions

def generate_folders(base_path, prefix, numbers):
    folders = [os.path.join(base_path, f"{prefix}{num}") for num in numbers]
    for folder in folders:
        if not os.path.exists(folder):
            print(f"Source folder '{folder}' does not exist.")
            folders.remove(folder)  # if the folder does not exsit, then remove the corresponding element from the folders list
    return folders


def delete_files_with_prefix(source_folders, prefix):
    try:
        for source_folder in source_folders:
            # Make sure the source folder exists
            if not os.path.exists(source_folder):
                print(f"Source folder '{source_folder}' does not exist.")
                continue

            # List all files in the source folder
            files = os.listdir(source_folder)

            # Delete files with the specified prefix
            for file_name in files:
                if file_name.startswith(prefix):
                    file_path = os.path.join(source_folder, file_name)

                    os.remove(file_path)

            print(f"Deleted files from {source_folder}")

        print(f"All files with prefix '{prefix}' deleted successfully.")

    except Exception as e:
        print(f"Error: {e}")

#%% Define root folder and all the sub-folders before deletion
base_source_path = "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/FIR-01B/results"
folder_prefix = "sub-00"
folder_numbers = list(range(31,45)) +list(range(46,50)) # Need to define clearly and accurately !!!

source_folders = generate_folders(base_source_path, folder_prefix, folder_numbers)

#%% Deleting files
delete_prefix = "beta"
delete_files_with_prefix(source_folders, delete_prefix)