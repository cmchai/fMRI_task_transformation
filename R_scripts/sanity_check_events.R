##########################################################
### Generate sanity check events for the fMRI analysis ###
### Author: Mengqiao Chai, chaimengqiao@gmail.com ###

### load package ###
library(here)
library(tidyverse)

## useful functions ----
fun_RespToDbl <- function(dataframe) {
  # in the case of missing trial, the response and rt would be null, which need to changed into double
  if (typeof(dataframe$rt) != "double"){
    dataframe$rt <- parse_double(dataframe$rt, na = "null")
    dataframe$response <- parse_double(dataframe$response, na = "null")
  }
  return(dataframe)
}

### specify input/output directories and the subjects to run ----
root_dir <- here("scanner")
root_o_dir <- "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/sanity_check/events"

all_subjs <- c(2:5,7:11,13:17,20:44,46:49)
subjs <- c(2,4,5) # subjects to run

for (subj in subjs) { # loop over each participant
  
  bids_subj <- paste0("sub-00", subj)
  
  input_subj_dir <- file.path(root_dir, bids_subj)
  output_subj_dir <- file.path(root_o_dir, bids_subj)
  
  # if the output directory does not exist, then create one
  if (!dir.exists(output_subj_dir)) {
    dir.create(output_subj_dir)
  }
  
  list_rawdata <- list.files(path = here(input_subj_dir),
                             recursive = TRUE,
                             pattern = "*.csv",
                             full.name = TRUE)
  
  
  bigdf <- list_rawdata %>%
    lapply(read_csv) %>%
    lapply(fun_RespToDbl) %>%
    bind_rows()
  
  ##################  generate key press timing txt files ##############-----
  key_run_dfs <- bigdf %>%
    select(subject, run_id, response, resp_onset) %>%
    filter(!is.na(response)) %>%
    group_by(run_id, response) %>%
    group_split()
  
  ## save as txt file ###
  
  for (i in 1:length(key_run_dfs)) {
    name_txt <- paste0(bids_subj, "_" ,unique(key_run_dfs[[i]]$run_id), "_keypress-", unique(key_run_dfs[[i]]$response), ".txt")
    write.table(key_run_dfs[[i]], file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  
  ##################  Generate encoding timing txt files ##############-----
  cue_run_dfs <- bigdf %>%
    select(subject, run_id, cue_onset) %>%
    group_by(run_id) %>%
    group_split()
  
  ## save as txt file ###
  
  for (i in 1:length(cue_run_dfs)) {
    name_txt <- paste0(bids_subj, "_" ,unique(cue_run_dfs[[i]]$run_id), "_cueOnset.txt")
    write.table(cue_run_dfs[[i]], file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  
  
  ##################  Generate stim onset timing txt files ##############-----
  stim_run_dfs <- bigdf %>%
    select(subject, run_id, stim_onset) %>%
    group_by(run_id) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(stim_run_dfs)) {
    name_txt <- paste0(bids_subj, "_" ,unique(stim_run_dfs[[i]]$run_id), "_stimOnset.txt")
    write.table(stim_run_dfs[[i]], file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  ## remove the data frame at the end of the each loop when this subject is done
  rm(bigdf) 
}



