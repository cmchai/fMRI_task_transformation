#################################################
### Generate all events for the fMRI analysis ###
#################### GLM-02 #####################
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
root_o_dir <- "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/GLM-02B/events"

all_subjs <- c(2:5,7:11,13:17,20:44,46:49)
subjs <- c(3:5,7:11,13:17,20:44,46:49) # subjects to run

sum_ses <- 2 # number of scanning sessions
sum_run <- 4 # number of functional runs per session

CTI_mid <- 5 # in second, IMPORTANT to define in order to put the long CTI onset regressor

########################################
#### for one subject, load the data ----
########################################

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
  # Here we picked only (1) Accurate, (2) Regular, (3) Long trials
  trials_int <- bigdf %>%
    filter(accuracy == TRUE & CTI_length == "long" & trial_nature == "RG")
    
  ####-- Event: CTI onset("c1") ----
  c1_dfs <- trials_int %>%
    select(subject, run_id, block_type, CTI_onset) %>%
    mutate(condition = paste0(block_type, "-long", "-", "c1"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(offset = CTI_onset + CTI_mid,
           duration = CTI_mid) %>%
    group_by(run_id, block_type) %>% # here, the block type argument is redundant since there is only one block type per block
    group_split()
  
  
  ## save as txt file ###
  for (i in 1:length(c1_dfs)) {
    c1_df <- c1_dfs[[i]]
    run_id <- unique(c1_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(c1_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(c1_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }

  ####-- Event: CTI midpoint("c2") ----
  c2_dfs <- trials_int %>%
    select(subject, run_id, block_type, CTI_onset, CTI_duration, CTI_offset) %>%
    mutate(condition = paste0(block_type, "-long", "-", "c2"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(onset = CTI_onset + CTI_mid, .after = block_type) %>%
    mutate(duration = CTI_duration - CTI_mid) %>%
    select(-c(CTI_onset, CTI_duration)) %>%
    group_by(run_id, block_type) %>% # here, the block type argument is redundant since there is only one block type per block
    group_split()
  
  
  ## save as txt file ###
  for (i in 1:length(c2_dfs)) {
    c2_df <- c2_dfs[[i]]
    run_id <- unique(c2_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(c2_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(c2_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }  
  
  #####-- Event: Stimulus ("s") ----
  s_dfs <- trials_int %>%
    filter(trial_nature == "RG") %>%
    select(subject, run_id, block_type, stim_onset, stim_offset) %>%
    mutate(condition = paste0(block_type, "-long", "-", "s"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = stim_offset - stim_onset) %>%
    group_by(run_id, block_type) %>%
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
  
  ## Events of No interests ----
  
  # first select the trials that we are not interested in but need to account for in GLM
  # in practice, the accuracy and the trial type(either regular or transform trials) would be relevant to pick up the trials of interest
  
  ### Events of No interest 1: error trials ----
  trials_error <- bigdf %>%
    filter(accuracy == FALSE)
  
  ####-- Event: CTI onset ("c1") for Error trials ----
  c1_error_dfs <- trials_error %>%
    select(subject, run_id, CTI_onset, CTI_offset) %>%
    mutate(condition = "error-c1", .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = CTI_offset - CTI_onset) %>%
    group_by(run_id) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(c1_error_dfs)) {
    c1_error_df <- c1_error_dfs[[i]]
    run_id <- unique(c1_error_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(c1_error_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(c1_error_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
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
  
  ### Events of No interest 2: accurate short trials ----
  trials_short <- bigdf %>%
    filter(accuracy == TRUE & CTI_length == "short")  
  
  ####-- Event: TI onset ("c1") for short trials ----
  c1_short_dfs <- trials_short %>%
    select(subject, run_id, block_type, CTI_onset, CTI_offset) %>%
    mutate(condition = paste0(block_type, "-short", "-", "c1"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = CTI_offset - CTI_onset) %>%
    group_by(run_id, block_type) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(c1_short_dfs)) {
    c1_short_df <- c1_short_dfs[[i]]
    run_id <- unique(c1_short_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(c1_short_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(c1_short_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }  
  
  ####-- Event: Stimulus("s") for short trials ----
  s_short_dfs <- trials_short %>%
    select(subject, run_id, block_type, stim_onset, stim_offset) %>%
    mutate(condition = paste0(block_type, "-short", "-", "s"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = stim_offset - stim_onset) %>%
    group_by(run_id, block_type) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(s_short_dfs)) {
    s_short_df <- s_short_dfs[[i]]
    run_id <- unique(s_short_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(s_short_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(s_short_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  
  
  ### Events of No interest 3: accurate transform trials ----
  trials_tran <- bigdf %>%
    filter(accuracy == TRUE & trial_nature == "TF")
  
  ####-- Event: CTI onset("c") ----
  c_tran_dfs <- trials_tran %>%
    select(subject, run_id, block_type, CTI_onset, CTI_offset) %>%
    mutate(condition = paste0(block_type, "-tran-long", "-", "c"), .after = run_id) %>% # for the subsequent SPM GLM
    mutate(duration = CTI_offset - CTI_onset) %>%
    group_by(run_id, block_type) %>% # here, the block type argument is redundant since there is only one block type per block
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(c_tran_dfs)) {
    c_tran_df <- c_tran_dfs[[i]]
    run_id <- unique(c_tran_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(c_tran_df$condition)
    name_txt = paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(c_tran_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }  
  
  #####-- Event: Transform Cues ("q2") ----
  q2_dfs <- trials_tran %>%
    select(subject, run_id, block_type, tranCue_onset, tranCue_offset) %>%
    mutate(condition = paste0(block_type, "-long", "-", "q2"), .after = run_id) %>% # for the subsequent SPM GLM
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
  rm(bigdf, trials_int, trials_error, trials_short, trials_tran)
}



