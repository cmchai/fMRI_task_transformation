 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 11:50:24 2024

Scripts that copy local files from one folder to another
Project: Task Transformation Paradigm

@author: mengqiao
"""
import shutil
import os

#%% All the useful functions

def generate_source_folder_list(base_path, prefix, numbers):
    folders = [os.path.join(base_path, f"{prefix}{num}") for num in numbers]
    for folder in folders:
        if not os.path.exists(folder):
            print(f"Source folder '{folder}' does not exist.")
            folders.remove(folder)  # if the folder does not exsit, then remove the corresponding element from the folders list
    return folders

def generate_dest_folders(new_base_path, source_folders):
    folders = [os.path.join(new_base_path, os.path.basename(folder)) for folder in source_folders]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder) # if the folder does not exist in the destination based path, then create it
    return folders


def copy_files_with_prefixes(source_folders, destination_folders, prefixes):
    try:
        if len(source_folders) != len(destination_folders):
            print("Error: Number of source and destination folders must be the same.")
            return

        for source_folder, destination_folder in zip(source_folders, destination_folders):
            # Make sure the source folder exists
            if not os.path.exists(source_folder):
                print(f"Source folder '{source_folder}' does not exist.")
                continue

            # Make sure the destination folder exists, create it if not
            if not os.path.exists(destination_folder):
                os.makedirs(destination_folder)

            # List all files in the source folder
            files = os.listdir(source_folder)

            # Copy files with the specified prefixes to the destination folder
            for file_name in files:
                for prefix in prefixes:
                    if file_name.startswith(prefix):
                        source_path = os.path.join(source_folder, file_name)
                        destination_path = os.path.join(destination_folder, file_name)

                        shutil.copy(source_path, destination_path) # IMPORTANT, if the files already exist in the destination folder, they will be overwritten

                        print(f"Copied: {file_name} from {source_folder} to {destination_folder}")

        print(f"All files with prefixes {prefixes} copied successfully.")

    except Exception as e:
        print(f"Error: {e}")

#%% Define root folder and all the sub-folders for both source and destination

base_source_path = "/Volumes/extdrive/Task_Transform_GLM/FIR-01B/results"
base_destination_path = "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/FIR-01B/just_test"
folder_prefix = "sub-00"
folder_numbers = list(range(2,6)) + list(range(7,12)) +list(range(13,18)) +list(range(20,31))

source_folders = generate_source_folder_list(base_source_path, folder_prefix, folder_numbers)
destination_folders = generate_dest_folders(base_destination_path, source_folders)

#%% Copy files with specific prefixes from the source folder to the destination folder
prefixes_to_copy = ["con_", "spmT_"]

copy_files_with_prefixes(source_folders, destination_folders, prefixes_to_copy)

