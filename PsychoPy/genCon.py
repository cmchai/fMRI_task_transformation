# -*- coding: utf-8 -*-
"""
Condition_Creation for the Task Transfrom Paradigm

Created on Wed Oct 26 17:13:37 2022

@author: chaim
"""

Run = False

import sys
import os
import random
import numpy as np
import pandas as pd
from datetime import datetime, date
from numpy import random as rdm

if Run == True:
    cwd = os.path.abspath(os.path.dirname("__file__")) # the path of the current working directory
    sys.path.append(cwd)
else:
    sys.path.append('D:\\Ghent_Braem\\Experiment\\1st_fMRI\\Experiment_script_forScan')

from imageList import stimulusList, taskList
#%% Functions

def genFullBlockID(fullBlockSeq):
    # block type alternate
    seqArray = ["1","1","2","2","3","3","4","4"] # 4 blocks of each block type interleaved
    fullBlockID = []
    for i in range(len(seqArray)):
        blockID = fullBlockSeq[i] + "_" + seqArray[i]
        fullBlockID.append(blockID)
    return fullBlockID

def genFullBlockID2(fullBlockSeq):
    # block type NOT alternate
    seqArray = ["1","2","3","4","1","2","3","4"] # 4 blocks of each block type NO interleaving
    fullBlockID = []
    for i in range(len(seqArray)):
        blockID = fullBlockSeq[i] + "_" + seqArray[i]
        fullBlockID.append(blockID)
    return fullBlockID

def genPrompt(rule, responseMappings):
    ''' generate the prompt during stimulus familiarization phase '''
    leftPrompt_0 = ["young", "small", "water"]
    rightPrompt_0 = ["old", "big", "land"]
    
    leftPrompt = [None] * 3
    rightPrompt = [None] * 3
    
    for i in range(len(responseMappings)):
        if responseMappings[i] == 0:
            leftPrompt[i] = leftPrompt_0[i]
            rightPrompt[i] = rightPrompt_0[i]
        else:
            leftPrompt[i] = rightPrompt_0[i]
            rightPrompt[i] = leftPrompt_0[i]
    
    prompDict = {}
    if rule == "age":
        prompDict["leftPrompt"] = leftPrompt[0]
        prompDict["rightPrompt"] = rightPrompt[0]
    elif rule == "size":
        prompDict["leftPrompt"] = leftPrompt[1]
        prompDict["rightPrompt"] = rightPrompt[1]
    elif rule == "location":
        prompDict["leftPrompt"] = leftPrompt[2]
        prompDict["rightPrompt"] = rightPrompt[2]
    return prompDict
# genPrompt("size")

def trialtype(logicalCTI):
    # if input is True, then transform trial
    trialType = "RG"
    if logicalCTI == True:
        trialType = "TF"
    return trialType
# trialtype(True)

def genAnswer(rule, imageID):
    answer = ""
    if rule == "age":
        answer = "young"
        if "_o" in imageID:
            answer = "old"
    elif rule == "size":
        answer = "small"
        if "_b" in imageID:
            answer = "big"
    elif rule == "location":
        answer = "water"
        if "_l" in imageID:
            answer = "land"
    return answer
# genAnswer("location", "animal_y_b_l")

def genCorrectKey(answer, responseMappings, leftKey, rightKey):
    responseMapping = ""
    correctKey = ""
    if answer == "young":
        responseMapping = responseMappings[0]
        correctKey = leftKey
        if responseMapping == 1:
            correctKey = rightKey
    elif answer == "old":
        responseMapping = responseMappings[0]
        correctKey = rightKey
        if responseMapping == 1:
            correctKey = leftKey
    elif answer == "small":
        responseMapping = responseMappings[1]
        correctKey = leftKey
        if responseMapping == 1:
            correctKey = rightKey
    elif answer == "big":
        responseMapping = responseMappings[1]
        correctKey = rightKey
        if responseMapping == 1:
            correctKey = leftKey
    elif answer == "water":
        responseMapping = responseMappings[2]
        correctKey = leftKey
        if responseMapping == 1:
            correctKey = rightKey
    elif answer == "land":
        responseMapping = responseMappings[2]
        correctKey = rightKey
        if responseMapping == 1:
            correctKey = leftKey
    return correctKey

# genCorrectKey("land", responseMappings)

def genRuleCue(rule, responseMappings):
    responseMapping = ""
    ruleCue = ""
    if rule == "age":
        responseMapping = responseMappings[0]
        ruleCue = "young | old"
        if responseMapping == 1:
            ruleCue = "old | young"
    elif rule == "size":
        responseMapping = responseMappings[1]
        ruleCue = "small | big"
        if responseMapping == 1:
            ruleCue = "big | small"
    elif rule == "location":
        responseMapping = responseMappings[2]
        ruleCue = "water | land"
        if responseMapping == 1:
            ruleCue = "land | water"
    return ruleCue

# genRuleCue("location", responseMappings)


def shuffleNoRepeat(Alist):
    ''' shuffle a list without repetition between consecutive elements,
    the list has to be a list of integers.'''
    taskRep = False
    listOK = False
    resultList = Alist.copy()

    while listOK == False:
        for i in range(len(resultList)-1):
            if resultList[i] == resultList[i+1]:
                taskRep = True
                break
        if taskRep == True:
            resultList = np.random.permutation(resultList) # result in a np array
            # print("reshuffling list continues...")
            taskRep = False
        else:
            listOK = True
            
    if type(resultList) == np.array: # convert to a list if the result is an numpy array
        resultList = resultList.tolist()
    return resultList # return a list

def shuffleArrayNoRepeat(AArray):
    ''' shuffle a 1d array without repetition between consecutive elements.'''
    taskRep = False
    listOK = False
    resultArray = AArray.copy()

    while listOK == False:
        for i in range(AArray.shape[0] - 1):
            if resultArray[i] == resultArray[i+1]:
                taskRep = True
                break
        if taskRep == True:
            resultArray = np.random.permutation(resultArray)
            # print("reshuffling array continues...")
            taskRep = False
        else:
            listOK = True
    return resultArray
# norman = np.array([3,5,6,6,7,9,10,10,14,16,76])
# shuf_norman = shuffleArrayNoRepeat(norman)

def shuffleArrayNoRepeatWithLabel(AArray, labels):
    ''' shuffle a 1d array without repetition between consecutive elements. Also return accordingly sorted labels,
        labels should also be an numpy array'''
    taskRep = False
    listOK = False
    resultArray = AArray.copy()
    resultLabels = labels.copy()
    index = np.arange(AArray.shape[0])

    while listOK == False:
        for i in range(AArray.shape[0] - 1):
            if resultArray[i] == resultArray[i+1]:
                taskRep = True
                break
        if taskRep == True:
            permuteIndex = np.random.permutation(index)
            resultArray = AArray[permuteIndex]
            resultLabels = labels[permuteIndex]
            # print("reshuffling array continues...")
            taskRep = False
        else:
            listOK = True
    return resultArray, resultLabels

# norman = np.array([3,5,6,6,7,9,10,10,14,16,76])
# norman2 = np.array([5,8,9,10,77,9,72,36,5,10,9])
# shortlong = np.arange(norman.shape[0])
# shuf_norman, shuf_labels = shuffleArrayNoRepeatWithLabel(norman2, shortlong)

def shuffleNestList(ANestList):
    '''shuffle a nested list with each sublist shuffled without repetition'''
    resultNestList = np.random.permutation(ANestList).tolist()
    
    for i in range(len(ANestList)):
        resultNestList[i] = shuffleNoRepeat(np.random.permutation(resultNestList[i]).tolist())
        
    return resultNestList
    
# a = [[1,2,3,4,4,2],[5,6,7,7,7,8],[9,10,11,12,9,11]]
# b = shuffleNestList(a)

def shuffle2dArray(A2dArray):
    '''shuffle a 2d-array with each sub-array shuffled without repetition'''
    resultArray = np.random.permutation(A2dArray)
    
    for i in range(A2dArray.shape[0]):
        resultArray[i] = shuffleArrayNoRepeat(np.random.permutation(resultArray[i]))
        
    return resultArray

def fromInttoDict(aInt, dictList, mappingID):
    '''transform a integer into a dictionary by a mapping id,
    the mapping id is the key of the dictionary, shall be a string '''
    resultDict = list(filter(lambda x:x[mappingID] == aInt, dictList))[0]
    return resultDict

def fromInttoDict_list(intList, dictList, mappingID):
    '''transform a integer list into a dictionary list by a mapping id,
    the mapping id is the key of the dictionary, shall be a string '''
    resultDictList = []
    for i in range(len(intList)):
        elementInt = intList[i]
        elementDict = list(filter(lambda x:x[mappingID] == elementInt, dictList))[0]
        resultDictList.append(elementDict)
    return resultDictList


def bet_sub_randomization(participant_num, reg_color):
    # the assignment of fixation color for each block type
    # seven block type alternation
    fixationColors = [(178, 0, 237), (59, 177, 67)] # violet and green for regular and transform block respectively
    if reg_color == "#3BB143":
        fixationColors = [(59, 177, 67), (178, 0, 237)] # green and violet for regular and transform block respectively
    
    # start either with RG or TF block
    if participant_num % 2 == 1: # if participant number is a odd number
        blockSeq = ["RG","TF","RG","TF","RG","TF","RG","TF"]
    else: # if participant number is a even number
        blockSeq = ["TF","RG","TF","RG","TF","RG","TF","RG"]
    
    fullBlockID = genFullBlockID(blockSeq)
    
    # determine the response mappings
    # responseMappings = random.choices([0,1], k = 3) # the responseMapping should be the same as the online pilot before
    
    return fixationColors, fullBlockID

def bet_sub_randomization2(participant_num, reg_color):
    # NO block type alternation
    # the assignment of fixation color for each block type
    fixationColors = [(178, 0, 237), (59, 177, 67)] # violet and green for regular and transform block respectively
    if reg_color == "#3BB143":
        fixationColors = [(59, 177, 67), (178, 0, 237)] # green and violet for regular and transform block respectively
    
    # start either with RG or TF block
    if participant_num % 2 == 1: # if participant number is a odd number
        blockSeq = ["RG","RG","RG","RG","TF","TF","TF","TF"]
    else: # if participant number is a even number
        blockSeq = ["TF","TF","TF","TF","RG","RG","RG","RG"]
    
    fullBlockID = genFullBlockID2(blockSeq)
    
    # determine the response mappings
    # responseMappings = random.choices([0,1], k = 3) # the responseMapping should be the same as the online pilot before
    
    return fixationColors, fullBlockID

def bet_sub_randomization3(participant_num, reg_color):
    # less block type alternation: 4 times block alternation
    # the assignment of fixation color for each block type
    fixationColors = [(178, 0, 237), (59, 177, 67)] # violet and green for regular and transform block respectively
    if reg_color == "#3BB143":
        fixationColors = [(59, 177, 67), (178, 0, 237)] # green and violet for regular and transform block respectively
    
    # start either with RG or TF block
    if participant_num % 2 == 1: # if participant number is a odd number
        blockSeq = ["RG","TF","TF","RG","RG","TF","TF","RG"]
    else: # if participant number is a even number
        blockSeq = ["TF","RG","RG","TF","TF","RG","RG","TF"]
    
    fullBlockID = genFullBlockID(blockSeq)
    
    # determine the response mappings
    # responseMappings = random.choices([0,1], k = 3) # the responseMapping should be the same as the online pilot before
    
    return fixationColors, fullBlockID

# fixationColors, fullBlockID, responseMappings = bet_sub_randomization(22)

#%% defining the task sequence for each block with their corresponding CTI

def genTaskIDs_CTIs(CTIs):
    '''generate the task IDs for a block with corresponding CTIs,
       the result is a randomized and balanced task sequence without repeating tasks between consecutive trials,
       and the tasks are ALSO balanced between short and long CTI intervals. '''
    blockTaskIDs_order = np.array(list(range(1,10)) * 6)
    is_CTIs_short_order = CTIs <= 5000
    index = np.arange(blockTaskIDs_order.shape[0]) # data type is 'int32'
    
    # first shuffle once before calling the function
    permuteIndex_start = np.random.permutation(index)
    blockTaskIDs = blockTaskIDs_order[permuteIndex_start]
    is_CTIs_short = is_CTIs_short_order[permuteIndex_start]
    
    resultTaskIDs, resultCTILabels = shuffleArrayNoRepeatWithLabel(blockTaskIDs, is_CTIs_short)
    
    # transform CTI labels into CTIs
    resultCTIs = np.zeros(CTIs.shape)
    shuffleShortCTIs = np.random.permutation(CTIs[is_CTIs_short_order])
    shuffleLongCTIs = np.random.permutation(CTIs[~is_CTIs_short_order])
    
    resultCTIs[resultCTILabels] = shuffleShortCTIs
    resultCTIs[~resultCTILabels] = shuffleLongCTIs
    
    return resultTaskIDs, resultCTIs

# bob, cindy = genTaskIDs_CTIs(globalVars["CTIs"])
# long_bob = bob[cindy > 5000]
# sort_longbob = np.sort(long_bob) # should be 3 repetition of each task
# short_bob = bob[cindy <= 5000]
# sort_shortbob = np.sort(short_bob) # should be 3 repetition of each task
# np.array_equal(sort_longbob, sort_shortbob) # should be true

def genTaskIDs_CTIs_8blk(CTIs):
    TaskIDs_8blk = np.zeros((8,54), dtype=int) # dtype=int
    CTIs_8blk = np.zeros((8,54), dtype=int) # dtype=int
    
    for i in range(8):
        blockTaskIDs, blockCTIs = genTaskIDs_CTIs(CTIs)
        TaskIDs_8blk[i,] = blockTaskIDs
        CTIs_8blk[i,] = blockCTIs
    return TaskIDs_8blk, CTIs_8blk

 # TaskIDs_8blk, CTIs_8blk = genTaskIDs_CTIs_8blk(CTIs)

#%% Defining all the stimulus pool for each semantic rule

animalImgList = list(filter(lambda x:x["stim_type"] == "animal", stimulusList))
placeImgList = list(filter(lambda x:x["stim_type"] == "place", stimulusList))
vehicleImgList = list(filter(lambda x:x["stim_type"] == "vehicle", stimulusList))

#%% define the stimlus for the 4 experimental blocks

def genTargetImgIDs_4blk():
    '''generate target images for 4 blocks,
    with each image repeat 9 times across 4 blocks to have balanced target images,
    the result is a 3*4*18(stimulus type * block * trials) matrix'''
    animalImgIDs = np.array(range(1,9))           # animal image IDs: 1-8;
    placeImgIDs = np.array(range(1+8, 9+8))       # place image IDs: 9-16;
    vehicleImgIDs = np.array(range(1+8+8, 9+8+8)) # vehicle image IDs: 17-24.
    
    target_ani_4blk = np.resize(animalImgIDs, 8*9) # 72, here numpy array is different operation on list
    target_pla_4blk = np.resize(placeImgIDs, 8*9)
    target_veh_4blk = np.resize(vehicleImgIDs, 8*9)
    
    target_ani_4blk_nest = np.array(np.split(target_ani_4blk, 4)) # split into 4 blocks, resulting in a 2-d array
    target_pla_4blk_nest = np.array(np.split(target_pla_4blk, 4))
    target_veh_4blk_nest = np.array(np.split(target_veh_4blk, 4))
    
    stack_targets_4blk_nest = np.stack((target_ani_4blk_nest, target_pla_4blk_nest, target_veh_4blk_nest)) # the first stimulus type(either the animal, place, or vehicle) is the target dimension
    
    return stack_targets_4blk_nest

# stack_targets_4blk = genTargetImgIDs_4blk()

# shuf_ani_4blk_nest = shuffle2dArray(target_ani_4blk_nest)
# shuf_pla_4blk_nest = shuffle2dArray(target_pla_4blk_nest)
# shuf_veh_4blk_nest = shuffle2dArray(target_veh_4blk_nest)

# b = np.sort(a)
# c = a.tolist()

#%% generate balanced stimulus for 4 blocks

def genStim4blk(stacked_targets_4blk):
    '''generate compound stimulus for 4 blocks,
    the dimension of result matrix is 4*54*3(block*trials*image),
    from 1-18 trials the target imgaes are animals;
    from 19-36 trials the target images are places;
    from 37-54 trials the target images are vehicles,
    the 1st element of the image dimension is the ID of the target image'''
    
    target_ani_4blk_nest = stacked_targets_4blk[0,:,:]
    target_pla_4blk_nest = stacked_targets_4blk[1,:,:]
    target_veh_4blk_nest = stacked_targets_4blk[2,:,:]
    
    result3dArray = np.zeros([4, 54, 3], dtype = int) # dtype = int
    targetStimTypeList = ['animal', 'place', 'vehicle']
    
    for targetStim in targetStimTypeList:
        
        shuf_ani_4blk_nest = shuffle2dArray(target_ani_4blk_nest)
        shuf_pla_4blk_nest = shuffle2dArray(target_pla_4blk_nest)
        shuf_veh_4blk_nest = shuffle2dArray(target_veh_4blk_nest)
        
        if targetStim == "animal":
            stim_4blk = np.stack((shuf_ani_4blk_nest, shuf_pla_4blk_nest, shuf_veh_4blk_nest))
            stim_4blk = np.transpose(stim_4blk, (1,2,0)) # dim: 4 blocks *18 trials * 3 images 
            result3dArray[:,0:18,:] = stim_4blk
        elif targetStim == "place":
            stim_4blk = np.stack((shuf_pla_4blk_nest, shuf_ani_4blk_nest, shuf_veh_4blk_nest))
            stim_4blk = np.transpose(stim_4blk, (1,2,0)) # dim: 4*18*3
            result3dArray[:,18:36,:] = stim_4blk
        else: 
            stim_4blk = np.stack((shuf_veh_4blk_nest, shuf_ani_4blk_nest, shuf_pla_4blk_nest))
            stim_4blk = np.transpose(stim_4blk, (1,2,0)) # dim: 4*18*3
            result3dArray[:,36:54,:] = stim_4blk
    return result3dArray
            
def genStim8blk(stack_targets_4blk):
    '''generate compound stimulus for 8 blocks,
    the dimension of result matrix is 8*54*3(block*trials*image),
    from 1-18 trials the target imgaes are animals;
    from 19-36 trials the target images are places;
    from 37-54 trials the target images are vehicles,
    the 1st element of the image dimension is the ID of the target image'''
    first_stim4blk = genStim4blk(stack_targets_4blk)
    second_stim4blk = genStim4blk(stack_targets_4blk)
    stim8blk = np.concatenate((first_stim4blk, second_stim4blk), axis=0)
    return stim8blk

# stim8blk = genStim8blk(stack_targets_4blk)
# imageID_array = targets_ani_4blk[:,1,2]
# imageDict_list = fromInttoDict(imageID_array.tolist(), stimulusList, "image_no")

#%% generate a block

def gen8blocks(blockIDs, TaskIDs_8blk, CTIs_8blk, stimFor8blk, taskDimensions, nTrialsBlock, CTIs, n_CTIs_long, numTran, leftKey, rightKey, responseMappings):
    
    blockTrialList = []
    
    for blockIndex, blockID in enumerate(blockIDs):
        blockType = blockID[0:2] # either TF or RG

        # generate task sequence for the block
        blockTaskIDs = TaskIDs_8blk[blockIndex]
        blockCTIs = CTIs_8blk[blockIndex]
        taskDict_block = fromInttoDict_list(blockTaskIDs, taskList, "task_id")
        
        blockStimIDs_unmatch = stimFor8blk[blockIndex,:,:] # shape:54*3
        blockStimIDs = np.zeros(blockStimIDs_unmatch.shape, dtype = int) # shape:54*3, dtype = int
        blockStimIDs[blockTaskIDs <= 3,:] =  blockStimIDs_unmatch[0:18, :] # the animal target trials
        blockStimIDs[blockTaskIDs > 6,:] =  blockStimIDs_unmatch[36:54, :] # the vehicle target trials
        blockStimIDs[blockStimIDs[:,0] == 0, :] = blockStimIDs_unmatch[18:36, :] # the place target trials
              
        if blockType == "TF":
            tran_CTIs = np.random.choice(CTIs[-n_CTIs_long:], numTran, replace=False) # select transform trials
            bool_tran = np.isin(blockCTIs, tran_CTIs)
        
        trialList = [] # creating the trials list which is a list of dictionaries
        
        for i in range(nTrialsBlock):
            trial_dict = {} # create the trial dictionary
            trial_dict["block"] = blockID
            trial_dict["block_type"] = blockType
            trial_dict["task"] = taskDict_block[i] # a dictionary within a dictionary
            trial_dict["rule_cue"] = genRuleCue(trial_dict["task"]["rule"], responseMappings)
            
            trial_stimIDs = blockStimIDs[i,:] # which is a np array of 3 elements, first one is the target image
            trial_dict["img_target"] = fromInttoDict(trial_stimIDs[0], stimulusList, "image_no")
            
            trial_stimPos = np.random.permutation(trial_stimIDs) # first is up, second is mid, last is low
            trial_stim = fromInttoDict_list(trial_stimPos, stimulusList, "image_no")
            trial_dict["img_up"] = trial_stim[0]
            trial_dict["img_mid"] = trial_stim[1]
            trial_dict["img_low"] = trial_stim[2]
            
            trial_dict["CTI"] = blockCTIs[i].item()
            
            if blockType == "TF": # in the case of a transform block
                trial_dict["trial_type"] = trialtype(bool_tran[i])
            else: # in the case of a regular block
                trial_dict["trial_type"] = "RG"
            
            if trial_dict["trial_type"] == "TF":# in the case of transform trial
                trial_dict["dim_tran"] = str(np.random.choice(taskDimensions, 1, replace=False)[0]) # either stim or rule
                if trial_dict["dim_tran"] == "stim": # if stim transform
                    alt_tasks = list(filter(lambda x:(x["stim"] != trial_dict["task"]["stim"] and
                                                      x["rule"] == trial_dict["task"]["rule"]), taskList))
                    trial_dict["task_tran"] = np.random.choice(alt_tasks, 1, replace=False)[0] # the new comfound task
                    trial_dict["rule_cue_tran"] = genRuleCue(trial_dict["task_tran"]["rule"], responseMappings) # should be the same as orginal cue
                    trial_dict["img_target_tran"] = list(filter(lambda x:x["stim_type"] == trial_dict["task_tran"]["stim"], trial_stim))[0]
                    trial_dict["answer"] = genAnswer(trial_dict["task_tran"]["rule"], trial_dict["img_target_tran"]["id"])
                else: # rule transform
                    alt_tasks = list(filter(lambda x:(x["stim"] == trial_dict["task"]["stim"] and
                                                      x["rule"] != trial_dict["task"]["rule"]), taskList))
                    trial_dict["task_tran"] = np.random.choice(alt_tasks, 1, replace=False)[0] # the new comfound task
                    trial_dict["rule_cue_tran"] = genRuleCue(trial_dict["task_tran"]["rule"], responseMappings) # should be the same as orginal cue
                    trial_dict["img_target_tran"] = trial_dict["img_target"]
                    trial_dict["answer"] = genAnswer(trial_dict["task_tran"]["rule"], trial_dict["img_target_tran"]["id"])
            else: # in the case of a regular trial
                trial_dict["answer"] = genAnswer(trial_dict["task"]["rule"], trial_dict["img_target"]["id"])
            
            trial_dict["correct_key"] = genCorrectKey(trial_dict["answer"], responseMappings, leftKey, rightKey)
            trialList.append(trial_dict)        
        blockTrialList.append(trialList)
    return blockTrialList

# blockTrial_full = gen8blocks(fullBlockID, stim8blk)
# len(list(filter(lambda x:x["trial_type"] == "TF", trialList))) == numTran # should be true 





