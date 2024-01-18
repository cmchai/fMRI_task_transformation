##########################################################
### Generate FIR-01 events for the fMRI analysis ###
### Author: Mengqiao Chai, chaimengqiao@gmail.com ###
## model: FIR-01, using finite impulse function to compare the BOLD dynamics after cue onset between regular and transform blocks

#### remove all the variables in the global environment
rm(list = ls())

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
root_o_dir <- "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/FIR-01/events"

all_subjs <- c(2:5,7:11,13:17,20:44,46:49)
subjs <- c(3:5,7:11,13:17,20:44,46:49) # subjects to run

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
    bind_rows() %>%
    mutate(CTI_length = if_else(CTI <= 5000, "short", "long"))
  
  ##################  Generate and Save event(*.txt) files ##############-----
  ### Events of Interests ----
  
  # first select the trials that we want to look at 
  # in practice, the accuracy and the trial type(either regular or transform trials) would be relevant to pick up the trials of interest
  trials_int <- bigdf %>%
    filter(accuracy == TRUE)
  
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

  ### Events of No interests ----
  
  # first select the trials that we are not interested in but need to account for in GLM
  # in practice, the accuracy and the trial type(either regular or transform trials) would be relevant to pick up the trials of interest
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

  ## remove the data frame at the end of the each loop when this subject is done
  rm(bigdf) 
}
