### Analysis on the Behave data of fMRI study ###

library(here)
library(tidyverse)
library(rstatix)

library(ez)
library(feather)
library(wesanderson)
library(ggsignif)
library(ggnewscale)
library(papaja)
library(rempsyc)

library(lmerTest)
library(lme4)
library(emmeans)

### all the functions to be used ---

fun_RespToDbl <- function(dataframe) {
  # in the case of missing trial, the response and rt would be null, which need to changed into double
  if (typeof(dataframe$rt) != "double"){
    dataframe$rt <- parse_double(dataframe$rt, na = "null")
    dataframe$response <- parse_double(dataframe$response, na = "null")
  }
  return(dataframe)
}


### import the csv data and concatenate them into an Dataframe----
subj_exclude <- c("subj_1") # subjects to be excluded from analysis

list_rawdata <- list.files(path = here("data"),
                           recursive = TRUE,
                           pattern = "*.csv",
                           full.name = TRUE)

bigdf <- list_rawdata %>%
  lapply(read_csv) %>%
  lapply(fun_RespToDbl) %>%
  bind_rows() %>%
  filter(!subject %in% subj_exclude) %>%
  mutate(CTI_length = if_else(CTI > 5000, "long", "short"))

unique(bigdf$subject)
glimpse(bigdf)

### Check all the balancing design ----

bigdf %>%
  group_by(CTI_length) %>%
  summarise(count = n())

bigdf %>%
  convert_as_factor(task_id) %>%
  group_by(CTI_length, task_id) %>%
  summarise(count = n())

img_balance <- bigdf %>%
  filter(block_type == "RG") %>%
  group_by(image_target_id) %>%
  summarise(count = n())
  
img_balance2 <- bigdf %>%
  filter(block_type == "TF") %>%
  group_by(image_target_id) %>%
  summarise(count = n())

### saving the data for the next stage of pre-processing ----
save(bigdf, file = paste0(here(), "/results/bigdf_fmri.Rdata"))








