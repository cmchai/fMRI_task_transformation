# -*- coding: utf-8 -*-
"""
main script for the Task Transfrom Paradigm

Created on Wed Oct 26 17:13:37 2022

@author: chaim
"""

debugging = True
running = True

#%% import all the modules and other script
import sys
import os
import random
import numpy as np
import pandas as pd
from datetime import datetime, date
from numpy import random as rdm

if running == True: # if running the experiment
    from psychopy import visual, core, event, gui
    cwd = os.path.abspath(os.path.dirname("__file__")) # the path of the current working directory

else:
    cwd = 'D:\\Ghent_Braem\\Experiment\\1st_fMRI\\Behaviour_Pilot\\2nd_InLab\\experiment_psychopy'

# print(cwd)
sys.path.append(cwd)
from imageList import stimulusList, taskList
import genCon

# b = os.path.abspath(getsourcefile(lambda:0))
# print(sys.argv)
# a = os.path.dirname(os.path.realpath("__file__"))
# b = os.path.realpath(os.path.dirname("__file__"))
#%% defining all the important global variables
globalVars = {}
globalVars["taskDims"] = ["stim", "rule"]
globalVars["stimFullList"] = ["animal", "place", "vehicle"]
globalVars["ruleFullList"] = ["age", "size", "location"]
globalVars["leftKey"] = "f"
globalVars["rightKey"] = "j"
globalVars["nTrialsBlock"] = 54
globalVars["nBlocks"] = 8
globalVars["intervalTaskCue"] = 1000 # 1000 ms
globalVars["intervalTaskCueTran"] = 750
globalVars["intervalStim"] = 2000
globalVars["ddlResponse"] = 2000
globalVars["intervalTrialFb"] = 750
globalVars["threBlockFb"] = 0.8 # the accuracy threshold, below which the block-wise feedback will show
globalVars["CTIs"] = np.array([1250, 
             	   			2125, 2250, 2375, 2500, 2625, 2750, 
             	   			2808, 2865, 2923, 2981, 3038, 3096, 3154, 3212, 3269, 3327, 3385, 3442, 3500, 
             	   			3625, 3750, 3875, 4000, 4125, 4250, 
             	   			5000, 
             	   			5001, 
             	   			5875, 6000, 6125, 6250, 6375, 6500, 
             	   			6558, 6615, 6673, 6731, 6788, 6846, 6904, 6962, 7019, 7077, 7135, 7192, 7250, 
             	   			7375, 7500, 7625, 7750, 7875, 8000, 
             	   			8750]) # 54 CTIs, mean is 5041 ms.
globalVars["n_CTIs_long"] = int(len(globalVars["CTIs"]) * 0.5) # fMRI pilot setting, half of them are long, 27 trials
globalVars["numTran"] = 13 # or 14
globalVars["ITIs"] = np.array([1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
                                1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,
                                1800,1800,1800,1800,1800,1800,1800,
                                2200,2200,2200,2200,2200,2200,
                                2600,2600,2600,2600,2600,2600,
                                3000,3000,3000,3000,3000,
                                3400,3400,3400,
                                3800,3800,
                                4200,
                                4600]) # 54 intervals
#%% collect demographic info
# Dialogue box
info = {'Participant ID': 'NaN', # the subject ID that is generated in the online session 
        'Counter balance No.': '999', # which is an number, used for counterbalancing the starting block
        'Age': '', 
        'Gender': ['Female','Male','Others','Prefer not to say'], 
        'Handedness': ['right-handed','left-handed'],
        'ResponseMappings': [(0,0,0),
                            (0,0,1),(0,1,0),(1,0,0),
                            (0,1,1),(1,1,0),(1,0,1),
                            (1,1,1)],
        'reg_color':["#B200ED", "#3BB143"]}
dlg = gui.DlgFromDict(dictionary = info, title = 'Experiment Setup', order = ['Participant ID', 'Counter balance No.', 'ResponseMappings', 'reg_color', 'Age', 'Gender', 'Handedness'])
if not dlg.OK: core.quit()
if info['Participant ID'] == '': core.quit()

if debugging == True:
    subjID = 'fakeID2'
    CB_no = 22
else:
    subjID = str(info['Participant ID'])
    CB_no = int(info['Counter balance No.'])

#%% defining all the between subject variables

fixationColors, blockIDs = genCon.bet_sub_randomization(CB_no, info['reg_color']) # if CB_no is an even number, then start with transform block
# importantly, the fixation color assignment should be consistent with the previous online training session

if fixationColors[0] == (178, 0, 237): # regular block is violet
    CTI_img_reg = 'images/red_cross.jpg'
    CTI_img_tra = 'images/green_cross.jpg'
else:
    CTI_img_reg = 'images/green_cross.jpg'
    CTI_img_tra = 'images/red_cross.jpg'

if running == True:    
    responseMappings_str = info['ResponseMappings']
    respMap_age = int(responseMappings_str[1])
    respMap_size = int(responseMappings_str[4])
    respMap_location = int(responseMappings_str[7])
    responseMappings = tuple([respMap_age, respMap_size, respMap_location])
else:
    responseMappings = tuple([0,1,0]) # as an example
# print("the response mappings tuple format are", responseMappings, ", the type of it is", type(responseMappings))

if responseMappings == (0,0,0):
    respMap_inst = "resp_map_0_0_0.jpg"
elif responseMappings == (0,0,1):
    respMap_inst = "resp_map_0_0_1.jpg"
elif responseMappings == (0,1,0):
    respMap_inst = "resp_map_0_1_0.jpg"
elif responseMappings == (1,0,0):
    respMap_inst = "resp_map_1_0_0.jpg"
elif responseMappings == (0,1,1):
    respMap_inst = "resp_map_0_1_1.jpg"
elif responseMappings == (1,1,0):
    respMap_inst = "resp_map_1_1_0.jpg"
elif responseMappings == (1,0,1):
    respMap_inst = "resp_map_1_0_1.jpg"
else:
    respMap_inst = "resp_map_1_1_1.jpg"
    
#%% create data frame using panda

# define columns
columns = ['subject', 'CB_no', 'age', 'gender', 'handedness', 'date',
           'reg_color', 'start_block', 'resp_map_age', 'resp_map_size', 'resp_map_location',
           'block_type', 'block_id',
           'trial_num', 'trial_id', 'trial_nature', 'stim', 'rule', 'image_up_id', 'image_mid_id', 'image_low_id', 'image_target_id', 'CTI',
           'dim_tran', 'task_tran', 'image_target_tran',
           'answer', 'correct_key', 
           'response', 'rt', 'accuracy', 'time_elapsed']
data = pd.DataFrame(columns = columns)
dateTimeNow = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")

# data file directory
datafile_path = cwd + '\\data\\'
print(datafile_path)

# datafile name
datafile_name = "{}data_{}.csv".format(datafile_path, str(subjID))

#%% create big block trial nested list for the looping
stack_targets_4blk = genCon.genTargetImgIDs_4blk()
stims8blk = genCon.genStim8blk(stack_targets_4blk)
blockTrialIterator = genCon.gen8blocks(blockIDs, stims8blk, globalVars["taskDims"], globalVars["nTrialsBlock"], 
                                       globalVars["CTIs"], globalVars["n_CTIs_long"], globalVars["numTran"], 
                                       globalVars["leftKey"], globalVars["rightKey"], responseMappings)


if debugging == True:
    TranTrials = list(map(lambda y:list(filter(lambda x:x["trial_type"] == "TF", y)), blockTrialIterator))
    CTIs_TranTrials = list(map(lambda y:[x['CTI'] for x in y], TranTrials))
    AllTaskDicts = list(map(lambda y:[x['task'] for x in y], blockTrialIterator))
    AllTaskIDs = list(map(lambda y:[x['task']['task_id'] for x in y], blockTrialIterator))

#%% create all the experiment psychopy objects

win = visual.Window(color=(1,1,1), fullscr = True) 
win.setMouseVisible(False)

instr_obj = visual.TextStim(win, text = '', wrapWidth = 1, font = 'monospace', color=(-1,-1,-1))
instr_img = visual.ImageStim(win, image = None, pos = (0, 0), units = "height")

taskCue_frame = visual.Rect(win, width=0.30, height=0.27, lineColor=(-1,-1,-1))
taskCue_obj = visual.TextStim(win, text='', height = 0.06, pos = (0,0),  bold=True, color=(-1,-1,-1))

CTIup_obj_img = visual.ImageStim(win, image = 'images/error.jpg', size = (0.005, 0.005), pos = (0, 0.25), units = "height")
CTImid_obj_img  = visual.ImageStim(win, image = 'images/error.jpg', size = (0.005, 0.005), pos = (0, 0), units = "height")
CTIlow_obj_img  = visual.ImageStim(win, image = 'images/error.jpg', size = (0.005, 0.005), pos = (0, -0.25), units = "height")

stimUp_obj = visual.ImageStim(win, image = None, size = (0.3, 0.2), pos = (0, 0.25), units = "height")
stimMid_obj = visual.ImageStim(win, image = None, size = (0.3, 0.2), pos = (0, 0), units = "height")
stimLow_obj = visual.ImageStim(win, image = None, size = (0.3, 0.2), pos = (0, -0.25), units = "height")

fb_obj = visual.TextStim(win, text='', height = 0.07, pos = (0,0), bold=True, color=(-1,-1,-1))

ITI_obj = visual.TextStim(win, text='+', height=0.05, pos = (0, 0),  bold=False, color=(-1,-1,-1))

#%% loop over each block and trial, present stimulus, record and save response for each trial

clock_global = core.Clock()
clock_cti = core.Clock()
clock_stim = core.Clock()

if debugging == True:
    blockIDs = blockIDs[0:2] # only run 2 blocks for debugging purpose

for blockIdx, blockID in enumerate(blockIDs):
    trialIterator = blockTrialIterator[blockIdx]
    if debugging == True:
        trialIterator = trialIterator[0:10] # only run 10 trials per block for debugging purpose
    random.shuffle(globalVars["ITIs"])
    #print(blockIdx, blockID)
    
    # show block instruction
    if blockID[0:2] == "RG":
        instr_img.setImage("regular_block_inst.jpg")
        CTIup_obj_img.setImage(CTI_img_reg)
        CTImid_obj_img.setImage(CTI_img_reg)
        CTIlow_obj_img.setImage(CTI_img_reg)
    else:
        instr_img.setImage("transform_block_inst.jpg")
        CTIup_obj_img.setImage(CTI_img_tra)
        CTImid_obj_img.setImage(CTI_img_tra)
        CTIlow_obj_img.setImage(CTI_img_tra)        
    
    instr_img.draw(); win.flip(); instr_resp = event.waitKeys(keyList=["space", "escape"])
    
    instr_img.setImage(respMap_inst); 
    instr_img.draw(); win.flip(); instr_resp = event.waitKeys(keyList=["space", "escape"])
    
    for trialIdx, trialDict in enumerate(trialIterator):
        # print(trialIdx)
        trialNum = trialIdx + 1 # since the index start from zero
        trialID = str(subjID) + "_" + str(blockID) + "_" + str(trialNum)
        
        # intertrial interval (ITI)
        ITI_obj.draw(); win.flip(); core.wait(globalVars["ITIs"][trialIdx]/1000)     
        
        # update and show task cue
        taskCue = "{}\n{}".format(trialDict["task"]["stim"], trialDict["rule_cue"])
        taskCue_obj.setText(taskCue); taskCue_frame.draw(); taskCue_obj.draw(); win.flip(); core.wait(globalVars["intervalTaskCue"]/1000)    
            
        clock_cti.reset()
        CTIup_obj_img.setSize(0.005); CTImid_obj_img.setSize(0.005); CTIlow_obj_img.setSize(0.005)
        while clock_cti.getTime() < round(trialDict["CTI"]/1000, 3):
            CTIup_obj_img.size += 0.00012
            CTImid_obj_img.size += 0.00012
            CTIlow_obj_img.size += 0.00012
            CTIup_obj_img.draw()
            CTImid_obj_img.draw()
            CTIlow_obj_img.draw()
            win.flip()
            
        # show transformed task cue if it's a transform trial
        if trialDict["trial_type"] == "TF":
            tranCue = "{}\n{}".format(trialDict["task_tran"]["stim"], trialDict["rule_cue_tran"])
            taskCue_obj.setText(tranCue); taskCue_frame.draw(); taskCue_obj.draw(); win.flip(); core.wait(globalVars["intervalTaskCueTran"]/1000)
        
        # show stimulus
        stimUp_obj.setImage(trialDict["img_up"]["full_name"])
        stimMid_obj.setImage(trialDict["img_mid"]["full_name"])
        stimLow_obj.setImage(trialDict["img_low"]["full_name"])
        stimUp_obj.draw(); stimMid_obj.draw(); stimLow_obj.draw(); win.flip()
        
        # time and record participants' response
        clock_stim.reset()
        respond = event.waitKeys(maxWait = globalVars["ddlResponse"]/1000, keyList=['f','j','escape'], timeStamped=clock_stim)
        
        # retrieve response, evaluate accuracy, and show feedback      
        if respond:
            response, rt = respond[0]
            if response == 'escape': win.close(); core.quit()
            if response == trialDict['correct_key']:
                acc = True
                fb_obj.setText('well done!')
            else:
                acc = False
                fb_obj.setText('Incorrect!')
        else:
            response = "null"
            rt = "null"
            acc = False
            fb_obj.setText('too slow!')
        
        fb_obj.draw(); win.flip(); core.wait(globalVars["intervalTrialFb"]/1000)
        
        # save trial data
        data = data.append({
            'subject': subjID, 
            'CB_no': CB_no, 
            'age': info['Age'], 
            'gender': info['Gender'],
            'handedness': info['Handedness'], 
            'date': dateTimeNow,
            'reg_color': fixationColors[0], 
            'start_block': blockIDs[0][0:2], 
            'resp_map_age': responseMappings[0], 
            'resp_map_size': responseMappings[1], 
            'resp_map_location': responseMappings[2],
            'block_type': blockID[0:2], 
            'block_id': blockID,
            'trial_num': trialNum,
            'trial_id': trialID, 
            'trial_nature': trialDict["trial_type"], 
            'stim': trialDict["task"]["stim"], 
            'rule': trialDict["task"]["rule"], 
            'image_up_id': trialDict["img_up"]["id"], 
            'image_mid_id': trialDict["img_mid"]["id"], 
            'image_low_id': trialDict["img_low"]["id"], 
            'image_target_id': trialDict["img_target"]["id"], 
            'CTI': trialDict["CTI"],
            'dim_tran': trialDict.get("dim_tran", np.nan), 
            'task_tran': trialDict.get("task_tran", np.nan), 
            'image_target_tran': trialDict.get("img_target_tran", np.nan),
            'answer': trialDict["answer"], 
            'correct_key': trialDict["correct_key"], 
            'response': response, 
            'rt': rt, 
            'accuracy': acc, 
            'time_elapsed' : clock_global.getTime()}, ignore_index=True)
    
    # show block-wise feedback based on the accuracy rate of both regular and transform trials
    blockFb = ""
    ACC_reg = data.loc[(data["block_id"] == blockID) & (data["trial_nature"] == "RG")]["accuracy"].mean()
    ACC_tran = data.loc[(data["block_id"] == blockID) & (data["trial_nature"] == "TF")]["accuracy"].mean()    
    if ACC_reg < globalVars["threBlockFb"]:
        ACC_reg_percent = int(ACC_reg*100)
        blockFb += "In the previous block, you have achieved {} percent of accuracy on the regular trials, which is, unfortunately, lower than other participants\'performance, please try harder to perform better.\n\n".format(ACC_reg_percent)
    if (np.isnan(ACC_tran) == False) & (ACC_tran < globalVars["threBlockFb"]):
        ACC_tran_percent = int(ACC_tran*100)
        blockFb += "In the previous block, you have achieved {} percent of accuracy on the transform trials, which is, unfortunately, lower than other participants\'performance, please try harder to perform better.\n\n".format(ACC_tran_percent)
    blockFb += "Please press the space bar to start the next block."          

    instr_obj.setText(blockFb); instr_obj.draw(); win.flip(); instr_resp = event.waitKeys(keyList=["space", "escape"])
            
#%% save data on the computer            
data.to_csv(datafile_name)
instr_datasave = "the data is being saved at the moment, please wait patiently for 15 seconds before the data storage is finished :)"
instr_obj.setText(instr_datasave); instr_obj.draw(); win.flip(); core.wait(15)

#%% close the whole experiment 
instr_goodbye = "Now you've finished the experiment ! \n\n Thanks a lot for your participation, now you can press spacebar to close the experiment interface."
instr_obj.setText(instr_goodbye); instr_obj.draw(); win.flip()
final_resp = event.waitKeys(keyList=["space", "escape"])
win.close()
