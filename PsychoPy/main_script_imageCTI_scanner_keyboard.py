# -*- coding: utf-8 -*-
"""
main script for the Task Transfrom Paradigm

Created on Wed Oct 26 17:13:37 2022

@author: chaim
"""
debugging = True
running = True
scanning = True

#%% import all the modules and other script
import sys
import os
import random
import numpy as np
import pandas as pd
import json
from datetime import datetime, date
from numpy import random as rdm

if running == True: # if running the experiment
    from psychopy import visual, core, event, gui, logging
    cwd = os.path.abspath(os.path.dirname("__file__")) # the path of the current working directory
else:
    cwd = 'D:\\Ghent_Braem\\Experiment\\1st_fMRI\\Experiment_script'

# print(cwd)
sys.path.append(cwd)
from imageList import stimulusList, taskList, subj_random
import genCon

if scanning == True:
    import scannertrigger as s # the module of 'pyxid' might not be installed

# b = os.path.abspath(getsourcefile(lambda:0))
# print(sys.argv)
# a = os.path.dirname(os.path.realpath("__file__"))
# b = os.path.realpath(os.path.dirname("__file__"))
#%% defining all the important global variables
globalVars = {}
globalVars["taskDims"] = ["stim", "rule"]
globalVars["stimFullList"] = ["animal", "place", "vehicle"]
globalVars["ruleFullList"] = ["age", "size", "location"]
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
globalVars["ITIs"] = np.array([1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 1750, 
                               2150, 2150, 2150, 2150, 2150, 2150, 2150, 2150, 2150, 2150, 2150, 
                               2550, 2550, 2550, 2550, 2550, 2550, 2550, 
                               2950, 2950, 2950, 2950, 2950, 2950, 
                               3350, 3350, 3350, 3350, 3350, 3350, 
                               3750, 3750, 3750, 3750, 3750, 
                               4150, 4150, 4150, 
                               4550, 4550, 
                               4950, 
                               5350]) # 54 intervals

if scanning == True: 
    globalVars["leftKey"] = "c"
    globalVars["rightKey"] = "d"
    keyList = ['c','d','escape']
else:
    globalVars["leftKey"] = "f"
    globalVars["rightKey"] = "j"
    keyList = ['f','j','escape']
   
#%% collect demographic info
# Dialogue box
info = {'Subject No.':  'NaN', # the subject number
        'Run No. Start':'1',
        'Run No. End':  '8'} 

if running == True:
    dlg = gui.DlgFromDict(dictionary = info, title = 'Experiment Setup', order = ['Subject No.', 'Run No. Start', 'Run No. End'])
    if not dlg.OK: core.quit()
    if info['Subject No.'] == '': core.quit()
    subjNo = int(info['Subject No.'])
    runNo_start = int(info['Run No. Start'])
    runNo_end = int(info['Run No. End'])
else:
    subjNo = random.randrange(1, 10)
    runNo_start = 3
    runNo_end = 4
    
subjID = "subj_" + str(info['Subject No.'])
n_runs = runNo_end - runNo_start + 1 # number of runs

#%% defining all the between subject variables
subj_param = genCon.fromInttoDict(subjNo, subj_random, "subjNo")

# fixation colors of 2 block types
fixationColors = [tuple(map(int, subj_param["reg_color"].split(', '))), # regular block color
                  tuple(map(int, subj_param["tran_color"].split(', ')))] # transform block color

if fixationColors[0] == (178, 0, 237): # regular block is violet
    CTI_img_reg = 'images/red_cross.jpg'
    CTI_img_tra = 'images/green_cross.jpg'
else:
    CTI_img_reg = 'images/green_cross.jpg'
    CTI_img_tra = 'images/red_cross.jpg'

# starting block and block sequence
if subj_param["start_block"] == "RG":
    blockIDs = genCon.genFullBlockID2(["RG","RG","RG","RG","TF","TF","TF","TF"])
else:
    blockIDs = genCon.genFullBlockID2(["TF","TF","TF","TF","RG","RG","RG","RG"])

# response mappings
responseMappings = tuple(map(int, subj_param["respMaps"].split(', ')))
respMap_age = int(responseMappings[0])
respMap_size = int(responseMappings[1])
respMap_location = int(responseMappings[2])

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
    
#%% Retrieve the corresponding block trial iterator based on the subject number
# read the json file and reconstruct the blocktrialiterator
with open(cwd +"\\iterators\\scan_keyboard\\exp_iterator_subj" + str(subjNo) + ".json", "r") as read_it:
      blockTrialIterator = json.load(read_it)

#%% create all the experiment psychopy objects
win = visual.Window(color=(1,1,1), fullscr = True) 
win.setMouseVisible(False)

instr_obj = visual.TextStim(win, text = '', wrapWidth = 1, font = 'monospace', color=(-1,-1,-1))
instr_img = visual.ImageStim(win, image = None, pos = (0, 0), size = 2, units = "norm")
instr_respMap = visual.ImageStim(win, image = respMap_inst, pos = (0, -0.4), size = (0.55, 0.45), units = "norm")

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
block_end = visual.TextStim(win, text='+', height=0.05, pos = (0, 0),  bold=False, color=(-1,-1,-1))

# Close if 'Escape' is pressed at any moment:
escKeyMessage = visual.TextStim(win,
                                text='ESC was pressed.\nClosing the experiment.',
                                height=0.07,
                                color=(-1, -1, -1))
def EscExit():
    print('Experiment terminated by the user (ESC)')
    # present a message for 2 seconds (120 frames)
    for frameN in range(120):
        escKeyMessage.draw()
        win.flip()
    win.close()
    core.quit()
event.globalKeys.add(key='escape', func=EscExit)

#%% detect the Cedrus button-box
if scanning == True:
    try: 
       import pyxid2 as pyxid
    except ImportError:
       import pyxid
    cedrusBox = None
    
    for n in range(10):  # doesn't always work first time!
        try:
            devices = pyxid.get_xid_devices()
            core.wait(0.1)
            cedrusBox = devices[0]
            break  # found the device so can break the loop
        except Exception:
            pass
    if not cedrusBox:
        logging.error('Could not find a Cedrus device.')
        core.quit()
    
    print('Cedrus device found: ' + cedrusBox.device_name)

#%% loop over each block and trial, present stimulus, record and save response for each trial
clock_global = core.Clock() # the time since the start of the script
clock_trigger = core.Clock()# the time since the trigger start --> which should be the benchmark time
clock_cti = core.Clock()
clock_stim = core.Clock()

# defining all the blocks to be run
if debugging == True:
    blockIDs = blockIDs[0:3] # only run 3 blocks for debugging purpose
else:
    blockIDs = blockIDs[runNo_start - 1:runNo_end] # including all the runs defined at the top of script

# loop over each block and each trial
for blockIdx, blockID in enumerate(blockIDs, runNo_start - 1):
    runID = "run_" + str(blockIdx + 1)
    trialIterator = blockTrialIterator[blockIdx]
    
    if debugging == True:
        trialIterator = trialIterator[0:5] # only run 5 trials per block for debugging purpose
        
    ############## VERY IMPORTANT !! ###################
    ################### shuffle the ITIs ###############    
    random.shuffle(globalVars["ITIs"]) # shuffle the ITI
    #print(blockIdx, blockID)

    #############################################
    ####### create data frame using panda #######
    #############################################
    # define columns #
    columns = ['subject', 'date',
               'reg_color', 'start_block', 'resp_map_age', 'resp_map_size', 'resp_map_location',
               'run_id', 'block_type', 'block_id',
               'trial_num', 'trial_id', 'trial_nature', 
               'task_id','stim', 'rule', 'image_up_id', 'image_mid_id', 'image_low_id', 'image_target_id', 'CTI',
               'dim_tran', 'task_tran', 'image_target_tran',
               'answer', 'correct_key', 
               'response', 'rt', 'accuracy', 
               ######### the following information all are the timing information for the fMRI GLM modelling ########
               ######################## all the timing information are the post-trigger time ########################
               'ITI_onset','ITI_duration','ITI_offset',
               'cue_onset','cue_duration','cue_offset',
               'CTI_onset','CTI_duration','CTI_offset',
               'tranCue_onset','tranCue_duration','tranCue_offset',
               'stim_onset','stim_duration','stim_offset',
               'resp_onset',
               'fb_onset','fb_duration','fb_offset',
               'run_time', 'time_elapsed']
    data = pd.DataFrame(columns = columns)
    dateTimeNow = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    
    # data file directory
    datafile_path = cwd + '\\data\\'
    print(datafile_path)
    
    # datafile name
    datafile_name = "{}behavioral_data_{}_{}.csv".format(datafile_path, str(subjID), str(runID))

    ######################################################################
    ############# block instruction and researcher start the #############
    ######################################################################
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
       
    instr_img.draw(); instr_respMap.draw(); win.flip()
    
    ############################################
    ######### synchronize the scanner ##########
    ############################################
    
    # Experimenter start the block by pressing the SPACE key
    start_block = event.waitKeys(keyList=["space", "escape"])
    
    if scanning == True:
        # scanner_port = 53
        stimScanner = visual.TextStim(
            win=win,
            name='text',
            text='Waiting for scanner...',
            font='Arial',
            pos=(0, 0),
            height=0.1,
            wrapWidth=None,
            ori=0,
            color=(-1,-1,-1),
            colorSpace='rgb',
            opacity=1,
            depth=0.0
            )
        stimScanner.draw()
        win.flip()
 
        portType = "cedrus"
        # Cedrus specific settings
        MR_settings = {
            'devicenr': 0,
            'sync': 4,
            }  
        
        # Open scanner trigger
        st = s.ScannerTrigger.create(
                win,
                clock_trigger,
                portType,
                portConfig=MR_settings,
                timeout=999999,
                #logLevel = logging.DATA,
                esc_key='escape')
        st.open()
        waitingTriggerOnset = clock_global.getTime()
        try:
            print("----------Waiting for Scanner...")
            triggered = st.waitForTrigger(skip=5) # It's True when it receives the 6th trigger(also the 6th TR)
            print('----------Trigger OK')
            #logging.flush()
        except Exception as e:
            # In case of errors
            #logging.flush()
            print("ScannerTrigger Error: {0}".format(e))
            core.quit()
            
    ############## VERY IMPORTANT !! ###################
    ######### reset the trigger clock !!!!!! ###########
    clock_trigger.reset() # reset the trigger time
    
    for trialIdx, trialDict in enumerate(trialIterator):
        # print(trialIdx)
        trialNum = trialIdx + 1 # since the index start from zero
        trialID = str(subjID) + "_" + str(blockID) + "_" + str(trialNum)
        
        # intertrial interval (ITI)
        ITI_obj.draw(); win.flip(); 
        ITI_onset = clock_trigger.getTime()
        core.wait(globalVars["ITIs"][trialIdx]/1000)
        ITI_duration = globalVars["ITIs"][trialIdx]/1000; ITI_offset = clock_trigger.getTime()
        
        # update and show task cue
        taskCue = "{}\n{}".format(trialDict["task"]["stim"], trialDict["rule_cue"])
        taskCue_obj.setText(taskCue); taskCue_frame.draw(); taskCue_obj.draw(); win.flip();
        cue_onset = clock_trigger.getTime()
        core.wait(globalVars["intervalTaskCue"]/1000)
        cue_duration = globalVars["intervalTaskCue"]/1000; cue_offset = clock_trigger.getTime()
            
        clock_cti.reset()
        CTIup_obj_img.setSize(0.005); CTImid_obj_img.setSize(0.005); CTIlow_obj_img.setSize(0.005)
        CTI_onset = clock_trigger.getTime()
        while clock_cti.getTime() < round(trialDict["CTI"]/1000, 3):
            CTIup_obj_img.size += 0.00012
            CTImid_obj_img.size += 0.00012
            CTIlow_obj_img.size += 0.00012
            CTIup_obj_img.draw()
            CTImid_obj_img.draw()
            CTIlow_obj_img.draw()
            win.flip()
        CTI_duration = trialDict["CTI"]/1000;  CTI_offset = clock_trigger.getTime()
        # show transformed task cue if it's a transform trial
        if trialDict["trial_type"] == "TF":
            tranCue = "{}\n{}".format(trialDict["task_tran"]["stim"], trialDict["rule_cue_tran"])
            taskCue_obj.setText(tranCue); taskCue_frame.draw(); taskCue_obj.draw(); win.flip() 
            tranCue_onset = clock_trigger.getTime()
            core.wait(globalVars["intervalTaskCueTran"]/1000)
            tranCue_duration = globalVars["intervalTaskCueTran"]/1000; tranCue_offset = clock_trigger.getTime()
        else:
            tranCue_onset = tranCue_duration = tranCue_offset = np.nan
        
        # show stimulus
        stimUp_obj.setImage(trialDict["img_up"]["full_name"])
        stimMid_obj.setImage(trialDict["img_mid"]["full_name"])
        stimLow_obj.setImage(trialDict["img_low"]["full_name"])
        stimUp_obj.draw(); stimMid_obj.draw(); stimLow_obj.draw(); win.flip()
        stim_onset = clock_trigger.getTime()
        
        # time and record participants' response
        clock_stim.reset()
        respond = event.waitKeys(maxWait = globalVars["ddlResponse"]/1000, keyList=keyList, timeStamped=clock_stim)
        stim_offset = clock_trigger.getTime()
        
        # retrieve response, evaluate accuracy, and show feedback
        if respond:
            resp_onset = clock_trigger.getTime()
            response, rt = respond[0]
            stim_duration = rt
            if response == 'escape': win.close(); core.quit()
            if response == trialDict['correct_key']:
                acc = True
                fb_obj.setText(' ')
            else:
                acc = False
                fb_obj.setText('Incorrect!')
        else:
            resp_onset = np.nan
            stim_duration = globalVars["ddlResponse"]/1000
            response = "null"
            rt = "null"
            acc = False
            fb_obj.setText('too slow!')
        
        fb_obj.draw(); win.flip(); 
        fb_onset = clock_trigger.getTime()
        core.wait(globalVars["intervalTrialFb"]/1000)
        fb_duration = globalVars["intervalTrialFb"]/1000; fb_offset = clock_trigger.getTime()
        
        # save trial data
        data = data.append({
            'subject': subjID, 
            'date': dateTimeNow,
            'reg_color': fixationColors[0], 
            'start_block': blockIDs[0][0:2], 
            'resp_map_age': responseMappings[0], 
            'resp_map_size': responseMappings[1], 
            'resp_map_location': responseMappings[2],
            'run_id': runID,
            'block_type': blockID[0:2], 
            'block_id': blockID,
            'trial_num': trialNum,
            'trial_id': trialID, 
            'trial_nature': trialDict["trial_type"],
            'task_id' : trialDict["task"]["task_id"],
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
            ########## all the timing information as follows #########
            'ITI_onset':ITI_onset, 'ITI_duration':ITI_duration, 'ITI_offset':ITI_offset,
            'cue_onset':cue_onset, 'cue_duration':cue_duration, 'cue_offset':cue_offset,
            'CTI_onset':CTI_onset, 'CTI_duration':CTI_duration, 'CTI_offset':CTI_offset,
            'tranCue_onset':tranCue_onset, 'tranCue_duration':tranCue_duration, 'tranCue_offset':tranCue_offset,
            'stim_onset':stim_onset, 'stim_duration':stim_duration, 'stim_offset':stim_offset,
            'resp_onset':resp_onset,
            'fb_onset':fb_onset, 'fb_duration':fb_duration, 'fb_offset':fb_offset,
            'run_time':clock_trigger.getTime(),
            'time_elapsed': clock_global.getTime()}, ignore_index=True)
 
    ######################################
    ########## save data !!!!!!! #########
    ######################################    
    data.to_csv(datafile_name)    
    
    ################################################################################################
    ### show block-wise feedback based on the accuracy rate of both regular and transform trials ###
    ################################################################################################
    if blockIdx < 7:
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
    
        instr_obj.setText(blockFb); instr_obj.draw(); win.flip()
    else:
        instr_goodbye = "Now you've finished the experiment ! \n\n Thanks a lot for your participation, now please wait for the researcher to come to you."
        instr_obj.setText(instr_goodbye); instr_obj.draw(); win.flip()
               
    core.wait(16)
    ######################################
    ########## close Cedrus port #########
    ######################################
    if scanning == True:
        st.close()
    print('Duration of ' + str(blockID) + ': ' + str(clock_trigger.getTime()) + "\n")

           
#%% close the whole experiment 
win.close()
