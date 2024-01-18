#################################################
### Generate all events for the fMRI analysis ###
#################### GLM-01 #####################
### Author: Mengqiao Chai, chaimengqiao@gmail.com

### import all packages ----
library(here)
library(tidyverse)

### functions to import ----

fun_RespToDbl <- function(dataframe) {
  # in the case of missing trial, the response and rt would be null, which need to changed into NA double
  if (typeof(dataframe$rt) != "double"){
    dataframe$rt <- parse_double(dataframe$rt, na = "null")
    dataframe$response <- parse_double(dataframe$response, na = "null")
  }
  return(dataframe)
}

### specify input/output directories and the subjects to run ----
root_dir <- here("scanner")
root_o_dir <- "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-01B/events"

all_subjs <- c(2:5,7:11,13:17,20:44,46:49)
subjs <- c(3:5,7:11,13:17,20:44,46:49) # subjects to run

sum_ses <- 2 # number of scanning sessions
sum_run <- 4 # number of functional runs per session

#### for one subject, load the data ----

for (subj in subjs) { # loop over each participant
  
  bids_subj <- paste0("sub-00", subj)
  
  input_subj_dir <- file.path(root_dir, bids_subj)
  output_subj_dir <- file.path(root_o_dir, bids_subj)
  
  # if the output directory does not exist, then create one
  if (!dir.exists(output_subj_dir)) {
    dir.create(output_subj_dir)
  }  
  
  list_rawdata <- list.files(path = input_subj_dir,
                             recursive = TRUE,
                             pattern = "*.csv",
                             full.name = TRUE)
  
  
  bigdf <- list_rawdata %>%
    lapply(read_csv) %>%
    lapply(fun_RespToDbl) %>%
    bind_rows() %>%
    mutate(CTI_length = if_else(CTI <= 5000, "short", "long"))
  
  # glimpse(bigdf)
  
  ##################  Generate and Save event(*.txt) files ##############-----
  ### Events of Interests ----
  
  # first select the trials that we want to look at 
  # in practice, the accuracy and the trial type(either regular or transform trials) would be relevant to pick up the trials of interest
  trials_int <- bigdf %>%
    filter(accuracy == TRUE & trial_nature == "RG")
  
  ####-- Event: Task cue + CTI ("qc") ----
  qc_dfs <- trials_int %>%
    select(subject, run_id, block_type, CTI_length, cue_onset, CTI_offset) %>%
    mutate(condition = paste0(block_type, "-", CTI_length, "-", "qc"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = CTI_offset - cue_onset) %>%
    group_by(run_id, block_type, CTI_length) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(qc_dfs)) {
    qc_df <- qc_dfs[[i]]
    run_id <- unique(qc_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(qc_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(qc_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  
  #####-- Event: Stimulus ("s") ----
  # Only for the regular trials, but not transform trials since the in transform trials, the CTI is followed by the transform task cue instead of stimulus
  s_dfs <- trials_int %>%
    select(subject, run_id, block_type, CTI_length, stim_onset, stim_offset) %>%
    mutate(condition = paste0(block_type, "-", CTI_length, "-", "s"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = stim_offset - stim_onset) %>%
    group_by(run_id, block_type, CTI_length) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(s_dfs)) {
    s_df <- s_dfs[[i]]
    run_id <- unique(s_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(s_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(s_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  
  ### Events of No interests ----
  
  # first select the trials that we are not interested in but need to account for in GLM
  # in practice, the accuracy and the trial type(either regular or transform trials) would be relevant to pick up the trials of interest
  
  #### Error trials ----
  trials_error <- bigdf %>%
    filter(accuracy == FALSE)
  
  ####-- Event: Task cue + CTI ("qc") for Error trials ----
  qc_error_dfs <- trials_error %>%
    select(subject, run_id, cue_onset, CTI_offset) %>%
    mutate(condition = "error-qc", .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = CTI_offset - cue_onset) %>%
    group_by(run_id) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(qc_error_dfs)) {
    qc_error_df <- qc_error_dfs[[i]]
    run_id <- unique(qc_error_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(qc_error_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(qc_error_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  
  ####-- Event: Stimulus("s") for Error trials ----
  s_error_dfs <- trials_error %>%
    select(subject, run_id, stim_onset, stim_offset) %>%
    mutate(condition = "error-s", .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = stim_offset - stim_onset) %>%
    group_by(run_id) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(s_error_dfs)) {
    s_error_df <- s_error_dfs[[i]]
    run_id <- unique(s_error_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(s_error_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(s_error_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  
  #### Accurate transform trials ----
  trials_tran <- bigdf %>%
    filter(accuracy == TRUE & trial_nature == "TF")
  
  ####-- Event: Task cue + CTI ("qc") ----
  qc_tran_dfs <- trials_tran %>%
    select(subject, run_id, block_type, CTI_length, cue_onset, CTI_offset) %>%
    mutate(condition = paste0(block_type, "-tran-", CTI_length, "-", "qc"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = CTI_offset - cue_onset) %>%
    group_by(run_id, block_type, CTI_length) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(qc_tran_dfs)) {
    qc_tran_df <- qc_tran_dfs[[i]]
    run_id <- unique(qc_tran_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(qc_tran_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(qc_tran_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
 
  #####-- Event: Transform Cues ("q2") ----
  # Only for transform trials, to account for the BOLD response right after the CTI offset
  q2_dfs <- trials_tran %>%
    select(subject, run_id, block_type, CTI_length, tranCue_onset, tranCue_offset) %>%
    mutate(condition = paste0(block_type, "-", CTI_length, "-", "q2"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = tranCue_offset - tranCue_onset) %>%
    group_by(run_id) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(q2_dfs)) {
    q2_df <- q2_dfs[[i]]
    run_id <- unique(q2_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(q2_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(q2_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }

  ## remove the data frame at the end of the each loop when this subject is done
  rm(bigdf, trials_int, trials_error, trials_tran)
}



