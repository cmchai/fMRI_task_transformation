# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 11:23:56 2023

Script to generate all the block_trial_iterator for each and every participants
The generated block_trial_iterators are saved as a JSON file for each participants

@author: chaim
"""
import sys
import os
import random
import numpy as np
import json

#%% set working directory
cwd = 'D:\\Ghent_Braem\\Experiment\\1st_fMRI\\Experiment_script_forScan'
sys.path.append(cwd)

from imageList import stimulusList, taskList, subj_random
import genCon

#%% define the response modality(either behav_keyboard, scan_keyboard, or scan_cedrus)
responseMod = "behav_keyboard"

#%% set some important global variables
globalVars = {}
globalVars["taskDims"] = ["stim", "rule"]
globalVars["nTrialsBlock"] = 54
globalVars["CTIs"] = np.array([1250, 
             	   			2125, 2250, 2375, 2500, 2625, 2750, 
             	   			2808, 2865, 2923, 2981, 3038, 3096, 3154, 3212, 3269, 3327, 3385, 3442, 3500, 
             	   			3625, 3750, 3875, 4000, 4125, 4250, 
             	   			5000, 
             	   			5001, 
             	   			5875, 6000, 6125, 6250, 6375, 6500, 
             	   			6558, 6615, 6673, 6731, 6788, 6846, 6904, 6962, 7019, 7077, 7135, 7192, 7250, 
             	   			7375, 7500, 7625, 7750, 7875, 8000, 
             	   			8750])
globalVars["n_CTIs_long"] = int(len(globalVars["CTIs"]) * 0.5) # fMRI pilot setting, half of them are long, 27 trials

globalVars["numTran"] = 13 # or 14
globalVars["n_subj"] = 50

if responseMod == "behav_keyboard":
    globalVars["leftKey"] = "f"
    globalVars["rightKey"] = "j"
    sub_folder = "\\iterators\\behav_keyboard"
elif responseMod == "scan_keyboard":
    globalVars["leftKey"] = "c"
    globalVars["rightKey"] = "d"
    sub_folder = "\\iterators\\scan_keyboard"
elif responseMod == "scan_cedrus":
    globalVars["leftKey"] = 2
    globalVars["rightKey"] = 3
    sub_folder = "\\iterators\\scan_cedrus"

#%% generate and store the blockTrialIterator in JSON file for each participant
for i in range(globalVars["n_subj"] - 40): # generate the iterators for 10 participants first
    
    # define subject number
    subjNo = i+1
    subj_param = genCon.fromInttoDict(subjNo, subj_random, "subjNo")
    
    # define block IDs
    if subj_param["start_block"] == "RG":
        blockIDs = genCon.genFullBlockID2(["RG","RG","RG","RG","TF","TF","TF","TF"])
    else:
        blockIDs = genCon.genFullBlockID2(["TF","TF","TF","TF","RG","RG","RG","RG"])
    
    # define resp mappings
    responseMappings = tuple(map(int, subj_param["respMaps"].split(', ')))

    # define task IDs and CTIs
    TaskIDs_8blk, CTIs_8blk = genCon.genTaskIDs_CTIs_8blk(globalVars["CTIs"])
    
    # define stimuli
    stack_targets_4blk = genCon.genTargetImgIDs_4blk()
    stims8blk = genCon.genStim8blk(stack_targets_4blk)
    
    # make the iterator
    blockTrialIterator = genCon.gen8blocks(blockIDs, TaskIDs_8blk, CTIs_8blk, stims8blk, 
                                            globalVars["taskDims"], globalVars["nTrialsBlock"],
                                            globalVars["CTIs"], globalVars["n_CTIs_long"], globalVars["numTran"], 
                                            globalVars["leftKey"], globalVars["rightKey"], responseMappings)
    
    # create and save the blocktrialiterator as a json file   
    with open(cwd + sub_folder + "\\exp_iterator_subj" + str(subjNo) + ".json", "w") as out_file:
        json.dump(blockTrialIterator, out_file, indent = 4)
