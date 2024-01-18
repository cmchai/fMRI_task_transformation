##########################################################
### Generate FIR-00 events for the fMRI analysis ###
### Author: Mengqiao Chai, chaimengqiao@gmail.com ###
## model: FIR-00, using finite impulse function to do sanity check on stimulus onset 

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
root_o_dir <- "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/first_level/univariate/FIR-00/events"

all_subjs <- c(2:5,7:11,13:17,20:44,46:49)
subjs <- c(4,5,7:11,13:17,20:44,46:49) # subjects to run

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
  
##################  Generate stim onset timing txt files ##############-----
  stim_run_dfs <- bigdf %>%
    select(subject, run_id, stim_onset, stim_offset) %>%
    mutate(duration = stim_offset - stim_onset) %>%
    mutate(condition = paste0("all-", "s"), .after = run_id) %>%
    group_by(run_id) %>%
    group_split()
  
  ## save as txt file ###
  for (i in 1:length(stim_run_dfs)) {
    stim_run_df <- stim_run_dfs[[i]]
    run_id <- unique(stim_run_df$run_id)
    run_no <- str_sub(run_id, -1) # which is a string but NOT number
    condition <- unique(stim_run_df$condition)
    name_txt <- paste0(bids_subj, "_run-", run_no, "_condition-", condition, ".txt")
    write.table(stim_run_df, file = file.path(output_subj_dir, name_txt), sep = "\t", row.names = FALSE)
  }
  ## remove the data frame at the end of the each loop when this subject is done
  rm(bigdf) 
}
