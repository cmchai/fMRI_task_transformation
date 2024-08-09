#################################################
######## Analysis of the Decoding Results #######
#################### GLM-02M and GLM-02M-Comp ####################
######### FOR MUltivairate task decoding ########
### Author: Mengqiao Chai, chaimengqiao@gmail.com

library(tidyverse)
library(ggh4x)
library(Hmisc)
library(ggprism)
library(rstatix)
library(ggpubr)
library(ez)

library(lmerTest)
library(lme4)
library(emmeans)

options(contrasts = c("contr.sum","contr.poly"))
options(scipen=999)

#####################################################
#################### importing data #################----
#####################################################

### for Schaefer atlas 

col_names = c("ROI_index", "ROI_label", "col_a", "col_b", "col_c", "col_d")
roi_labels_df = read.table('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas/schaefer_2018/schaefer_2018/Schaefer2018_400Parcels_17Networks_order.txt',
                           sep='\t', 
                           header=FALSE,col.names = col_names) %>%
  tibble() %>%
  select("ROI_index", "ROI_label")

setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/independent_roi/conA_conB_defB')
results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/independent_roi/decodeAcc_smthN_spmT_rois_contAB_defB.csv') %>%
  mutate(task_dim = case_when(str_detect(condition, "stim") ~ "stim",
                              str_detect(condition, "rule") ~ "rule",
                              .default = "conjunc"),
         block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
         CTI_window = if_else(str_detect(condition, "c1"), "short", "long")) %>%
  left_join(roi_labels_df, by = join_by(roi == ROI_index)) %>%
  mutate(roi_label = ROI_label) %>%
  separate(ROI_label, into = c("network_amount", "hemisphere", "network_label", "parcel_label", "parcel_extra_nr"), sep = "_", extra = "merge")
# filter(grepl("\\[", roi)) # add if the ROIs are combined together before decoding for Schaefer atlas

### for Glasser atlas
setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/Glasser')

Glasser_SupParcels <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/Glasser/decodeAcc_smthN_spmT_rois_FPN_Glasser.csv') %>%
  mutate(task_dim = case_when(str_detect(condition, "stim") ~ "stim",
                              str_detect(condition, "rule") ~ "rule",
                              .default = "conjunc"),
         block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
         CTI_window = if_else(str_detect(condition, "c1"), "short", "long")) %>%
  mutate(
    roi_glas = str_extract(roi, "^[^_]+"),
    hemisphere = str_extract(roi, "(?<=_)[^.]+")
  ) %>%
  left_join(Glasser_SupParcels, by = join_by(roi_glas == fpn_labels))

grouping_check <- results %>%
  group_by(task_dim, block_type, CTI_window, condition) %>%
  summarise(count = n())

acc_byDimRegion <- results %>%
  group_by(task_dim, roi) %>%
  summarise(mean_roi_size = mean(roi_size),
            mean_acc = mean(mean_accuracy, na.rm = TRUE))

### Visualizing the result of conjunctive decoding ###----
roi_mapping <- results %>%
  group_by(roi, hemisphere, network_label, parcel_label) %>%
  summarise(count = n())

ROIs <- unique(results$roi)

task_dim <- c()
rois <- c()
p_block <- c()
p_cti <- c()
p_inter <- c()

for (ROI in ROIs) {
  inter_data <- results %>%
    filter(roi == ROI, task_dim == "conjunc")
  
  # statistical analysis
  inter_anova <- ezANOVA(data = inter_data, 
                         dv = .(mean_accuracy), 
                         wid = .(subject), 
                         within = .(block_type, CTI_window))
  
  p_vals <- round(inter_anova$ANOVA$p, 3)
  
  task_dim <- c(task_dim, "conjunc")
  rois <- c(rois, ROI)
  p_block <- c(p_block, p_vals[1])
  p_cti <- c(p_cti, p_vals[2])
  p_inter <- c(p_inter, p_vals[3])
  
  # summary of the data for plotting
  inter_summary <- inter_data %>%
    group_by(block_type, CTI_window) %>%
    summarise(submean = mean(mean_accuracy, na.rm = TRUE),
              subsd = sd(mean_accuracy, na.rm = TRUE),
              se = subsd/((43)^.5)) %>%   # standard error
    ungroup()
  
  
  p4 <- ggplot(inter_summary, aes(x = CTI_window, y = submean, group = block_type, color = block_type)) +
    scale_x_discrete(name ="CTI window", 
                     limits=c("short", "long"))+
    geom_point(size = 5) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                  width = .2) +
    scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                       name = "block type:",
                       breaks = c("RG", "TF")) +
    geom_hline(aes(yintercept=0.1111),linetype="dashed", size=0.5) +
    coord_cartesian(ylim=c(0.10, 0.16)) +
    labs(y = "mean decoding acc") +
    ggtitle(paste0("conjunctive task ", "roi:", ROI)) +
    theme_bw()
    
  ggsave(paste0("conjunctive task_", "roi_", ROI, ".png"), width = 3.3, height = 3)
}


### Result of compositional decoding --- stimulus type ###----
for (ROI in ROIs) {
  inter_data <- results %>%
    filter(roi == ROI, task_dim == "stim")
  
  # statistical analysis
  inter_anova <- ezANOVA(data = inter_data,
                         dv = .(mean_accuracy),
                         wid = .(subject),
                         within = .(block_type, CTI_window))

  p_vals <- round(inter_anova$ANOVA$p, 3)

  task_dim <- c(task_dim, "stim")
  rois <- c(rois, ROI)
  p_block <- c(p_block, p_vals[1])
  p_cti <- c(p_cti, p_vals[2])
  p_inter <- c(p_inter, p_vals[3])
  
  # summary of the data for plotting
  inter_summary <- inter_data %>%
    group_by(block_type, CTI_window) %>%
    summarise(submean = mean(mean_accuracy, na.rm = TRUE),
              subsd = sd(mean_accuracy, na.rm = TRUE),
              se = subsd/((43)^.5)) %>%   # standard error
    ungroup()
  
  
  p4 <- ggplot(inter_summary, aes(x = CTI_window, y = submean, group = block_type, color = block_type)) +
    scale_x_discrete(name ="CTI window", 
                     limits=c("short", "long"))+
    geom_point(size = 5) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                  width = .2) +
    scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                       name = "block type:",
                       breaks = c("RG", "TF")) +
    geom_hline(aes(yintercept=0.3333),linetype="dashed", size=0.5) +
    coord_cartesian(ylim=c(0.3, 0.45)) +
    labs(y = "mean decoding acc") +
    ggtitle(paste0("stimulus task ", "roi:", ROI)) +
    theme_bw()
  
  # ggsave(paste0("stim-task_", "roi_", ROI, ".png"), width = 3.3, height = 3)
}


### Result of compositional decoding --- task rule ###----
for (ROI in ROIs) {
  inter_data <- results %>%
    filter(roi == ROI, task_dim == "rule")
  
  # statistical analysis
  inter_anova <- ezANOVA(data = inter_data, 
                         dv = .(mean_accuracy), 
                         wid = .(subject), 
                         within = .(block_type, CTI_window))
  
  p_vals <- round(inter_anova$ANOVA$p, 3)
  
  task_dim <- c(task_dim, "rule")
  rois <- c(rois, ROI)
  p_block <- c(p_block, p_vals[1])
  p_cti <- c(p_cti, p_vals[2])
  p_inter <- c(p_inter, p_vals[3])
  
  # summary of the data for plotting
  inter_summary <- inter_data %>%
    group_by(block_type, CTI_window) %>%
    summarise(submean = mean(mean_accuracy, na.rm = TRUE),
              subsd = sd(mean_accuracy, na.rm = TRUE),
              se = subsd/((43)^.5)) %>%   # standard error
    ungroup()
  
  
  p4 <- ggplot(inter_summary, aes(x = CTI_window, y = submean, group = block_type, color = block_type)) +
    scale_x_discrete(name ="CTI window", 
                     limits=c("short", "long"))+
    geom_point(size = 5) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                  width = .2) +
    scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                       name = "block type:",
                       breaks = c("RG", "TF")) +
    geom_hline(aes(yintercept=0.3333),linetype="dashed", size=0.5) +
    coord_cartesian(ylim=c(0.3, 0.45)) +
    labs(y = "mean decoding acc") +
    ggtitle(paste0("rule task ", "roi:", ROI)) +
    theme_bw()
  
  # ggsave(paste0("rule-task_", "roi_", ROI, ".png"), width = 3.3, height = 3)
}

decoding_anova <- tibble(task_dim,rois,p_block,p_cti,p_inter)
save(decoding_anova, file = "decoding_anova_results_rois_FLM.Rdata") 
load(file = "decoding_anova_results_rois_selectFromShort.Rdata")

sig_block <- decoding_anova %>%
  filter(p_block < 0.05)

sig_inter <- decoding_anova %>%
  filter(p_inter < 0.05)

############################################################
###### statistical analysis on Schaefer atlas results ######----
############################################################

### only include parcels from AFC, dlPFC, vlPFC, and PPC

parcels_network_fpn <- c("ContA_IPS", "ContB_IPL","ContA_PFClv", "ContA_PFCl", "ContB_PFClv", "DefaultB_PFCv")
fpn_SupParcels <- c("iPL", "iPL", "dlPF", "dlPF", "aPF", "vlPF")
unique_fpn_SupParcels <- c("iPL", "dlPF", "aPF", "vlPF")
hemispheres <- c("LH", "RH")

schaefer_parcel_to_fpn_sup <- tibble(parcels_network_fpn, fpn_SupParcels)

parcel_retain <- c("IPS", "IPL","PFClv", "PFCl", "PFCv") # control A PFClv belong to dlPFC, control B PFClv belong to AFC

results_schaefer_fpn <- results %>%
  filter(parcel_label %in% parcel_retain) %>%
  mutate(network_parcel = paste(network_label, parcel_label, sep = "_")) %>%
  left_join(schaefer_parcel_to_fpn_sup, by = join_by(network_parcel == parcels_network_fpn))

# check super parcel mapping

sup_parcel_mappings <- results_schaefer_fpn %>%
  group_by(fpn_SupParcels, hemisphere, network_parcel) %>%
  summarise(count = n())

# run analysis 
task_dims = c("conjunc","stim","rule")

t_dim <- c()
t_SupParcels <- c()
p_block <- c()
p_cti <- c()
p_hem <- c()
p_block_by_cti <- c()
p_block_by_hem <- c()
p_cti_by_hem <- c()
p_3way <- c()

for (i in 1:length(task_dims)) { # loop over each task dimension
  for (ii in 1:length(unique_fpn_SupParcels)) { # loop over each super parcel
    
    data_interest <- results_schaefer_fpn %>%
      filter(task_dim == task_dims[i],
             fpn_SupParcels == unique_fpn_SupParcels[ii]) %>%
      convert_as_factor(subject, block_type, CTI_window, hemisphere, roi)
    
    ## statistical analysis
    # -- ANOVA, which results from sub parcels are averaged 
    # inter_anova <- ezANOVA(data = data_interest, 
    #                        dv = .(mean_accuracy), 
    #                        wid = .(subject), 
    #                        within = .(block_type, CTI_window, hemisphere))
    # 
    # p_vals <- round(inter_anova$ANOVA$p, 3)
    # p_block <- c(p_block, p_vals[1])
    # p_cti <- c(p_cti, p_vals[2])
    # p_inter2 <- c(p_inter2, p_vals[4])
    # p_inter3 <- c(p_inter3, p_vals[7])
    
    # using mixed effect model
    inter_lme <- lmer(formula = mean_accuracy ~ block_type * CTI_window * hemisphere + (1 | subject),
                     data = data_interest)

    table_summary <- round((summary(inter_lme, correlation= FALSE))$coefficients, 4)
    
    t_dim <- c(t_dim, task_dims[i])
    t_SupParcels <- c(t_SupParcels, unique_fpn_SupParcels[ii])
    p_block <- c(p_block, table_summary[2,5])
    p_cti <- c(p_cti, table_summary[3,5])
    p_hem <- c(p_hem, table_summary[4,5])
    p_block_by_cti <- c(p_block_by_cti, table_summary[5,5])
    p_block_by_hem <- c(p_block_by_hem, table_summary[6,5])
    p_cti_by_hem <- c(p_cti_by_hem, table_summary[7,5])
    p_3way <- c(p_3way, table_summary[8,5])
    
    ### visualizing the interaction for each super parcel for each hemisphere separately
    # for (iii in 1:length(hemispheres)) {
    #   
    #   hem = hemispheres[iii]
    #   
    #   inter_summary <- data_interest %>%
    #     filter(hemisphere == hem) %>%
    #     group_by(block_type, CTI_window) %>%
    #     summarise(submean = mean(mean_accuracy, na.rm = TRUE),
    #               subsd = sd(mean_accuracy, na.rm = TRUE),
    #               se = subsd/((43)^.5)) %>%   # standard error
    #     ungroup()
    #   
    #   if (task_dims[i] == "conjunc") {
    #     ymin = 0.10
    #     ymax = 0.16
    #     chance = 0.1111
    #   } else {
    #     ymin = 0.28
    #     ymax = 0.47
    #     chance = 0.3333
    #   }
    #   
    #   p4 <- ggplot(inter_summary, aes(x = CTI_window, y = submean, group = block_type, color = block_type)) +
    #     scale_x_discrete(name ="CTI window", 
    #                      limits=c("short", "long"))+
    #     geom_point(size = 5) +
    #     geom_line(size = 1) +
    #     geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
    #                   width = .2) +
    #     scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
    #                        name = "block type:",
    #                        breaks = c("RG", "TF")) +
    #     geom_hline(aes(yintercept=chance),linetype="dashed", size=0.5) +
    #     coord_cartesian(ylim=c(ymin, ymax)) +
    #     labs(y = "mean decoding acc") +
    #     ggtitle(paste0("decoding acc in ", task_dims[i], " for ", hemispheres[iii], " ", unique_fpn_SupParcels[ii])) +
    #     theme_bw()
    #   
    #   ggsave(paste0(task_dims[i], "_SupParcel_", hemispheres[iii], "_", unique_fpn_SupParcels[ii], ".png"), width = 3.7, height = 3.2)
    #}
  }
}

decoding_scheafer_fpn_stats <- tibble(t_dim, t_SupParcels, p_block, p_cti, p_hem, p_block_by_cti, p_block_by_hem, p_cti_by_hem, p_3way)
save(decoding_scheafer_fpn_stats, file = "decoding_lme_results_rois_FPN_Schaefer.Rdata")
load(file = "decoding_lme_results_rois_FPN_Schaefer.Rdata")

decoding_scheafer_fpn_stats_oneroi <- decoding_scheafer_fpn_stats %>%
  filter(t_SupParcels == "aPF")

## conjunc decoding analysis of diff super parcels in one model

conjunc_data_fpn <- results_schaefer_fpn %>%
  filter(task_dim == "conjunc",
         fpn_SupParcels %in% c("iPL", "dlPFC", "aPF")) %>%
  convert_as_factor(subject, block_type, CTI_window, hemisphere, fpn_SupParcels)

## using mixed effect model ----

# re-level CTI before fitting the model
conjunc_data_fpn$CTI_window <- relevel(conjunc_data_fpn$CTI_window, ref = "short")

big_lme <- lmer(formula = mean_accuracy ~ block_type * CTI_window * hemisphere * fpn_SupParcels + (1 | subject),
                data = conjunc_data_fpn)

emmip(big_lme, block_type ~ CTI_window | fpn_SupParcels)

summary(big_lme, correlation=FALSE)

levels(conjunc_data_fpn$fpn_SupParcels)
levels(conjunc_data_fpn$CTI_window)
levels(conjunc_data_fpn$block_type)

#####################################################
#################### FIR—2M results##################----
#####################################################

TR_sec = 1.78 # in second(s)
CTI_mid = 5+1 # in second(s),including the task cue
CTI_longest = 8.75+1 # in second(s),including the task cue

setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/roi_short_trials')
results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/roi_short_trials/decodeAcc_smthN_spmT_combined_rois_selectFromShort.csv') %>%
  mutate(task_dim = case_when(str_detect(condition, "stim") ~ "stim",
                              str_detect(condition, "rule") ~ "rule",
                              .default = "conjunc"),
         block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
         TR = as.integer(str_extract(condition, "\\d+$")))

TR_check <- results %>%
  group_by(TR, condition) %>%
  summarise(count=n())

## visualizing result --- stim decoding
result_stim_summary <- results %>%
  filter(task_dim == "stim") %>%
  group_by(roi, block_type, TR) %>%
  summarise(submean = mean(mean_accuracy, na.rm = TRUE),
            subsd = sd(mean_accuracy, na.rm = TRUE),
            se = subsd/((43)^.5)) %>%   # standard error
  ungroup()

(
p_exp <- ggplot(result_stim_summary, aes(x = TR, y = submean, group = block_type, color = block_type)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                width = .2) +
  scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                     name = "block type:",
                     breaks = c("RG", "TF")) +
  scale_x_continuous(breaks=seq(0,10,1)) +
  geom_hline(aes(yintercept=0.3333),linetype="dashed", size=0.5) +
  geom_vline(aes(xintercept=CTI_mid/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
  geom_vline(aes(xintercept=CTI_longest/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
  coord_cartesian(ylim=c(0.28, 0.47)) +
  labs(y = "mean decoding acc") +
  ggtitle("stim decoding") +
  theme_bw() +
  facet_wrap(~ roi, nrow = 4)
)

## visualizing result --- rule decoding
result_rule_summary <- results %>%
  filter(task_dim == "rule") %>%
  group_by(roi, block_type, TR) %>%
  summarise(submean = mean(mean_accuracy, na.rm = TRUE),
            subsd = sd(mean_accuracy, na.rm = TRUE),
            se = subsd/((43)^.5)) %>%   # standard error
  ungroup()

(
  p_exp <- ggplot(result_rule_summary, aes(x = TR, y = submean, group = block_type, color = block_type)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                  width = .2) +
    scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                       name = "block type:",
                       breaks = c("RG", "TF")) +
    scale_x_continuous(breaks=seq(0,10,1)) +
    geom_hline(aes(yintercept=0.3333),linetype="dashed", size=0.5) +
    geom_vline(aes(xintercept=CTI_mid/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
    geom_vline(aes(xintercept=CTI_longest/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
    coord_cartesian(ylim=c(0.28, 0.43)) +
    labs(y = "mean decoding acc") +
    ggtitle("rule decoding") +
    theme_bw() +
    facet_wrap(~ roi, nrow = 4)
)

## saving the plot and doing analysis for each ROI

ROIs <- unique(results$roi)
task_dims <- unique(results$task_dim)

task_dim <- c() 
rois <- c()
p_block <- c()
p_tr <- c()
p_inter <- c()

for (dim in task_dims) {
  for (ROI in ROIs) {
    inter_data <- results %>%
      filter(roi == ROI, task_dim == dim)
    
    # statistical analysis
    inter_anova <- ezANOVA(data = inter_data, 
                           dv = .(mean_accuracy), 
                           wid = .(subject), 
                           within = .(block_type, TR))
    
    p_vals <- round(inter_anova$ANOVA$p, 3)
    
    task_dim <- c(task_dim, dim)
    rois <- c(rois, ROI)
    p_block <- c(p_block, p_vals[1])
    p_tr <- c(p_tr, p_vals[2])
    p_inter <- c(p_inter, p_vals[3])
    
    # summary of the data for plotting
    inter_summary <- inter_data %>%
      group_by(block_type, TR) %>%
      summarise(submean = mean(mean_accuracy, na.rm = TRUE),
                subsd = sd(mean_accuracy, na.rm = TRUE),
                se = subsd/((43)^.5)) %>%   # standard error
      ungroup()
    
    
    p_exp <- ggplot(inter_summary, aes(x = TR, y = submean, group = block_type, color = block_type)) +
      geom_point(size = 5) +
      geom_line(size = 1.7) +
      geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                    width = .2) +
      scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                         name = "block type:",
                         breaks = c("RG", "TF")) +
      scale_x_continuous(breaks=seq(0,10,1)) +
      geom_hline(aes(yintercept=0.3333),linetype="dashed", size=0.5) +
      geom_vline(aes(xintercept=CTI_mid/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
      geom_vline(aes(xintercept=CTI_longest/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
      coord_cartesian(ylim=c(0.28, 0.47)) +
      labs(y = "mean decoding acc") +
      ggtitle(paste0(dim, " decoding for ", "roi:", ROI)) +
      theme_bw()
    
    ggsave(paste0(dim, "_decoding_FIR_", "roi_", ROI, ".png"), width = 7, height = 3)
  }
}


############################
###### HCP-HMM atlas #######----
############################

mmp <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/mmp.csv')

fpn_labels <- c("SCEF", "8BM", "8C", "IFJp", "p9-46v", "a9-46v", "i6-8", "AVI", "IP1", "IP2", "PFm",
                "LIPd", "MIP", "AIP", "POS2", "p47r", "a32pr", "6r", "s6-8", "a10p", "11l", "TE1m",
                "TE1p", "FOP5", "p10p", "d32", "a47r", "PGs")

fpn_SupParcels <- c("mPF", "mPF", "dlPF", "dlPF", "dlPF", "aPF", "dPF", "Ins", "iPL","iPL", "iPL",
                    "iPL", "iPL", "iPL", "PreCu", "aPF", "mPF", "dlPF", "dPF", "aPF", "aPF", "mT",
                    "mT", "Ins", "aPF", "mPF", "aPF", "iPL")

# mPF: medial prefrontal; dlPF: dorsal lateral prefrontal; aPF: anterior prefrontal; dPF: dorsal prefrontal
# Ins: insular; iPL: inferior parietal lobe; PreCu: precuneus; mT: medial temporal

fpn_nrs <- mmp %>%
  filter(roi %in% fpn_labels) %>%
  pull(num.roi)

fpn_nrs

MD_SupParcels <- tibble(fpn_labels, fpn_SupParcels)
write_csv(MD_SupParcels, file = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

###########################################################
###### statistical analysis on MMP (Glasser) results ######----
###########################################################

# check the super parcel mapping
supparcel_map <- results %>%
  group_by(fpn_SupParcels, roi_glas) %>%
  summarise(count = n())

# only select fpn related ROIs
results_only_fpn <- results %>%
  filter(fpn_SupParcels %in% c("aPF", "dlPF", "iPL"))

conditions <- results_only_fpn %>%
  group_by(task_dim, fpn_SupParcels) %>%
  summarise(count = n())

task_dims <- conditions$task_dim
sup_parcels <- conditions$fpn_SupParcels
hemispheres <- c("left", "right")

t_dim <- c()
t_SupParcels <- c()
p_block <- c()
p_cti <- c()
p_hem <- c()
p_block_by_cti <- c()
p_block_by_hem <- c()
p_cti_by_hem <- c()
p_3way <- c()

for (i in 1:length(task_dims)) {
  
  data_interest <- results_only_fpn %>%
    filter(task_dim == task_dims[i],
           fpn_SupParcels == sup_parcels[i]) %>%
    convert_as_factor(subject, block_type, CTI_window, hemisphere, roi)
  
  # using mixed effect model
  inter_lme <- lmer(formula = mean_accuracy ~ block_type * CTI_window * hemisphere + (1 | subject),
                    data = data_interest)
  
  table_summary <- round((summary(inter_lme, correlation= FALSE))$coefficients, 4)
  
  t_dim <- c(t_dim, task_dims[i])
  t_SupParcels <- c(t_SupParcels, sup_parcels[i])
  p_block <- c(p_block, table_summary[2,5])
  p_cti <- c(p_cti, table_summary[3,5])
  p_hem <- c(p_hem, table_summary[4,5])
  p_block_by_cti <- c(p_block_by_cti, table_summary[5,5])
  p_block_by_hem <- c(p_block_by_hem, table_summary[6,5])
  p_cti_by_hem <- c(p_cti_by_hem, table_summary[7,5])
  p_3way <- c(p_3way, table_summary[8,5])
  
  ### visualizing the interaction for each super parcel separately for each hemisphere
  # for (ii in 1:length(hemispheres)) {
  #   
  #   inter_summary <- data_interest %>%
  #     filter(hemisphere == hemispheres[ii]) %>%
  #     group_by(block_type, CTI_window) %>%
  #     summarise(submean = mean(mean_accuracy, na.rm = TRUE),
  #               subsd = sd(mean_accuracy, na.rm = TRUE),
  #               se = subsd/((43)^.5)) %>%   # standard error
  #     ungroup()
  #   
  #   if (task_dims[i] == "conjunc") {
  #     ymin = 0.10
  #     ymax = 0.16
  #     chance = 0.1111
  #   } else {
  #     ymin = 0.28
  #     ymax = 0.47
  #     chance = 0.3333
  #   }
  #   
  #   p4 <- ggplot(inter_summary, aes(x = CTI_window, y = submean, group = block_type, color = block_type)) +
  #     scale_x_discrete(name ="CTI window", 
  #                      limits=c("short", "long"))+
  #     geom_point(size = 5) +
  #     geom_line(size = 1) +
  #     geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
  #                   width = .2) +
  #     scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
  #                        name = "block type:",
  #                        breaks = c("RG", "TF")) +
  #     geom_hline(aes(yintercept=chance),linetype="dashed", size=0.5) +
  #     coord_cartesian(ylim=c(ymin, ymax)) +
  #     labs(y = "mean decoding acc") +
  #     ggtitle(paste0("decoding acc in ", task_dims[i], " for ", hemispheres[ii], " ", sup_parcels[i])) +
  #     theme_bw()
  #   
  #   ggsave(paste0(task_dims[i], "_SupParcel_", hemispheres[ii], " ", sup_parcels[i], ".png"), width = 3.3, height = 3)
  #}
}

decoding_glasser_fpn_stats <- tibble(t_dim, t_SupParcels, p_block, p_cti, p_hem, p_block_by_cti, p_block_by_hem, p_cti_by_hem, p_3way)
save(decoding_glasser_fpn_stats, file = "decoding_lme_results_rois_FPN_glasser.Rdata")
load(file = "decoding_lme_results_rois_FPN_glasser.Rdata")
load(file = "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/independent_roi/conA_conB_defB/decoding_lme_results_rois_FPN_Schaefer.Rdata")

## conjunc decoding analysis of diff super parcels in one model

conjunc_data_fpn <- results_only_fpn %>%
  filter(task_dim == "conjunc") %>%
  convert_as_factor(subject, block_type, CTI_window, hemisphere, fpn_SupParcels)

## using mixed effect model

# re-level CTI before fitting the model
conjunc_data_fpn$CTI_window <- relevel(conjunc_data_fpn$CTI_window, ref = "short")

big_lme <- lmer(formula = mean_accuracy ~ block_type * CTI_window * hemisphere * fpn_SupParcels + (1 | subject),
                  data = conjunc_data_fpn)

emmip(big_lme, block_type ~ CTI_window | fpn_SupParcels)

summary(big_lme, correlation=FALSE)

levels(conjunc_data_fpn$fpn_SupParcels)
levels(conjunc_data_fpn$CTI_window)
levels(conjunc_data_fpn$block_type)

############# garbage code #############

folds <- results %>%
  filter(task_dim == "conjunc") %>%
  group_by(subject, n_folds) %>%
  summarise(count = n())

##################################################################
########### Cross-decoding results --- cross block type ##########
##################################################################

setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/Glasser/cross_decoding')

Glasser_SupParcels <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/Glasser/cross_decodeAcc_smthN_spmT_rois_FPN_Glasser.csv') %>%
  mutate(task_dim = case_when(str_detect(training, "stim") ~ "stim",
                              str_detect(training, "rule") ~ "rule",
                              .default = "conjunc"),
         block_type_train = if_else(str_detect(training, "RG-"), "RG", "TF"),
         block_type_test = if_else(str_detect(testing, "RG-"), "RG", "TF"),
         CTI_window_train = if_else(str_detect(training, "c1"), "short", "long"),
         CTI_window_test = if_else(str_detect(testing, "c1"), "short", "long"),
         cross_block = paste(block_type_train, block_type_test, sep = "-"),
         cross_cti = paste(CTI_window_train, CTI_window_test, sep = "-")) %>%
  mutate(
    roi_glas = str_extract(roi, "^[^_]+"),
    hemisphere = str_extract(roi, "(?<=_)[^.]+")
  ) %>%
  left_join(Glasser_SupParcels, by = join_by(roi_glas == fpn_labels))

# only select fpn related ROIs
results_only_fpn <- results %>%
  filter(fpn_SupParcels %in% c("aPF", "dlPF", "iPL"))

sup_parcels <- c("aPF", "dlPF", "iPL")
hemispheres <- c("left", "right")

t_SupParcels <- c()
p_cross_block <- c()
p_cross_cti <- c()
p_hem <- c()
p_block_by_cti <- c()
p_block_by_hem <- c()
p_cti_by_hem <- c()
p_3way <- c()

for (i in 1:length(sup_parcels)) {
  
  data_interest <- results_only_fpn %>%
    filter(fpn_SupParcels == sup_parcels[i]) %>%
    convert_as_factor(subject, cross_block, cross_cti, hemisphere, roi)
  
  # using mixed effect model
  inter_lme <- lmer(formula = mean_accuracy ~ cross_block * cross_cti * hemisphere + (1 | subject),
                    data = data_interest)
  
  table_summary <- round((summary(inter_lme, correlation= FALSE))$coefficients, 4)
  
  t_SupParcels <- c(t_SupParcels, sup_parcels[i])
  p_cross_block <- c(p_cross_block, table_summary[2,5])
  p_cross_cti <- c(p_cross_cti, table_summary[3,5])
  p_hem <- c(p_hem, table_summary[4,5])
  p_block_by_cti <- c(p_block_by_cti, table_summary[5,5])
  p_block_by_hem <- c(p_block_by_hem, table_summary[6,5])
  p_cti_by_hem <- c(p_cti_by_hem, table_summary[7,5])
  p_3way <- c(p_3way, table_summary[8,5])
  
  ### visualizing the interaction for each super parcel separately for each hemisphere
  for (ii in 1:length(hemispheres)) {

    inter_summary <- data_interest %>%
      filter(hemisphere == hemispheres[ii]) %>%
      group_by(cross_block, cross_cti) %>%
      summarise(submean = mean(mean_accuracy, na.rm = TRUE),
                subsd = sd(mean_accuracy, na.rm = TRUE),
                se = subsd/((43)^.5)) %>%   # standard error
      ungroup()

      ymin = 0.10
      ymax = 0.15
      chance = 0.1111

    p4 <- ggplot(inter_summary, aes(x = cross_cti, y = submean, group = cross_block, color = cross_block)) +
      scale_x_discrete(name ="CTI window",
                       limits=c("short-short", "long-long"))+
      geom_point(size = 5) +
      geom_line(size = 1) +
      geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                    width = .2) +
      scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                         name = "cross type:",
                         breaks = c("RG-TF", "TF-RG")) +
      geom_hline(aes(yintercept=chance),linetype="dashed", size=0.5) +
      coord_cartesian(ylim=c(ymin, ymax)) +
      labs(y = "mean cross-decoding acc") +
      ggtitle(paste0("Cross decoding acc in conjunc", " for ", hemispheres[ii], " ", sup_parcels[i])) +
      theme_bw()

    ggsave(paste0("CrossDecode_conjunc_SupParcel_", hemispheres[ii], " ", sup_parcels[i], ".png"), width = 3.7, height = 3.2)
  }
}

cross_decoding_glasser_fpn_stats <- tibble(t_SupParcels, p_cross_block, p_cross_cti, p_hem, p_block_by_cti, p_block_by_hem, p_cti_by_hem, p_3way)
save(cross_decoding_glasser_fpn_stats, file = "cross_decoding_lme_results_rois_FPN_glasser.Rdata")
load(file = "cross_decoding_lme_results_rois_FPN_glasser.Rdata")

## conjunc decoding analysis of diff super parcels in one model
## using mixed effect model

# re-level CTI before fitting the model
cross_decode_glasser_lme <- results_only_fpn %>%
  convert_as_factor(subject, cross_block, cross_cti, hemisphere, roi)

cross_decode_glasser_lme$cross_cti <- relevel(cross_decode_glasser_lme$cross_cti, ref = "short-short")

big_lme <- lmer(formula = mean_accuracy ~ cross_block * cross_cti * hemisphere * fpn_SupParcels + (1 | subject),
                data = cross_decode_glasser_lme)

summary(big_lme, correlation=FALSE)

##################################################################
############# Cross-decoding results --- cross CTI ###############
##################################################################

setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/Glasser/cross_decoding')

Glasser_SupParcels <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/Glasser/cross_decoding/cross_cti_decodeAcc_smthN_spmT_rois_FPN_Glasser.csv') %>%
  mutate(task_dim = case_when(str_detect(training, "stim") ~ "stim",
                              str_detect(training, "rule") ~ "rule",
                              .default = "conjunc"),
         block_type_train = if_else(str_detect(training, "RG-"), "RG", "TF"),
         block_type_test = if_else(str_detect(testing, "RG-"), "RG", "TF"),
         CTI_window_test = if_else(str_detect(testing, "c1"), "short", "long")) %>%
  mutate(
    roi_glas = str_extract(roi, "^[^_]+"),
    hemisphere = str_extract(roi, "(?<=_)[^.]+")
  ) %>%
  left_join(Glasser_SupParcels, by = join_by(roi_glas == fpn_labels))

# only select fpn related ROIs
results_only_fpn <- results %>%
  filter(fpn_SupParcels %in% c("aPF", "dlPF", "iPL"))

sup_parcels <- c("aPF", "dlPF", "iPL")
hemispheres <- c("left", "right")

t_SupParcels <- c()
p_block <- c()
p_cti <- c()
p_hem <- c()
p_block_by_cti <- c()
p_block_by_hem <- c()
p_cti_by_hem <- c()
p_3way <- c()

for (i in 1:length(sup_parcels)) {
  
  data_interest <- results_only_fpn %>%
    filter(fpn_SupParcels == sup_parcels[i]) %>%
    convert_as_factor(subject, block_type_test, CTI_window_test, hemisphere, roi)
  
  # using mixed effect model
  inter_lme <- lmer(formula = mean_accuracy ~ block_type_test * CTI_window_test * hemisphere + (1 | subject),
                    data = data_interest)
  
  table_summary <- round((summary(inter_lme, correlation= FALSE))$coefficients, 4)
  
  t_SupParcels <- c(t_SupParcels, sup_parcels[i])
  p_block <- c(p_block, table_summary[2,5])
  p_cti <- c(p_cti, table_summary[3,5])
  p_hem <- c(p_hem, table_summary[4,5])
  p_block_by_cti <- c(p_block_by_cti, table_summary[5,5])
  p_block_by_hem <- c(p_block_by_hem, table_summary[6,5])
  p_cti_by_hem <- c(p_cti_by_hem, table_summary[7,5])
  p_3way <- c(p_3way, table_summary[8,5])
  
  ### visualizing the interaction for each super parcel separately for each hemisphere
  for (ii in 1:length(hemispheres)) {
    
    inter_summary <- data_interest %>%
      filter(hemisphere == hemispheres[ii]) %>%
      group_by(block_type_test, CTI_window_test) %>%
      summarise(submean = mean(mean_accuracy, na.rm = TRUE),
                subsd = sd(mean_accuracy, na.rm = TRUE),
                se = subsd/((43)^.5)) %>%   # standard error
      ungroup()
    
    ymin = 0.10
    ymax = 0.15
    chance = 0.1111
    
    p4 <- ggplot(inter_summary, aes(x = CTI_window_test, y = submean, group = block_type_test, color = block_type_test)) +
      scale_x_discrete(name ="CTI window",
                       limits=c("short", "long"))+
      geom_point(size = 5) +
      geom_line(size = 1) +
      geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                    width = .2) +
      scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                         name = "block type:",
                         breaks = c("RG", "TF")) +
      geom_hline(aes(yintercept=chance),linetype="dashed", size=0.5) +
      coord_cartesian(ylim=c(ymin, ymax)) +
      labs(y = "mean cross-decoding acc") +
      ggtitle(paste0("Cross CTI decode conjunc for ", hemispheres[ii], " ", sup_parcels[i])) +
      theme_bw()
    
    ggsave(paste0("CrossCTI_Decode_conjunc_SupParcel_", hemispheres[ii], " ", sup_parcels[i], ".png"), width = 3.9, height = 3.4)
  }
}

cross_cti_decoding_glasser_fpn_stats <- tibble(t_SupParcels, p_block, p_cti, p_hem, p_block_by_cti, p_block_by_hem, p_cti_by_hem, p_3way)
save(cross_cti_decoding_glasser_fpn_stats, file = "cross_cti_decoding_lme_results_rois_FPN_glasser.Rdata")
load(file = "cross_cti_decoding_lme_results_rois_FPN_glasser.Rdata")

## conjunc decoding analysis of diff super parcels in one model
## using mixed effect model

# re-level CTI before fitting the model
cross_cti_decode_glasser_lme <- results_only_fpn %>%
  convert_as_factor(subject, block_type_test, CTI_window_test, hemisphere, roi)

cross_cti_decode_glasser_lme$CTI_window_test <- relevel(cross_cti_decode_glasser_lme$CTI_window_test, ref = "short")

big_lme <- lmer(formula = mean_accuracy ~ block_type_test * CTI_window_test * hemisphere * fpn_SupParcels + (1 | subject),
                data = cross_cti_decode_glasser_lme)

summary(big_lme, correlation=FALSE)

######################################################
####### cross cti decoding use Schaefer atlas ########
######################################################

col_names = c("ROI_index", "ROI_label", "col_a", "col_b", "col_c", "col_d")
roi_labels_df = read.table('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas/schaefer_2018/schaefer_2018/Schaefer2018_400Parcels_17Networks_order.txt',
                           sep='\t', 
                           header=FALSE,col.names = col_names) %>%
  tibble() %>%
  select("ROI_index", "ROI_label")

setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/independent_roi/conA_conB_defB/cross_decode')

results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/independent_roi/conA_conB_defB/cross_decode/cross_cti_decodeAcc_smthN_spmT_rois_FPN_Schaefer.csv') %>%
  mutate(block_type = if_else(str_detect(testing, "RG-"), "RG", "TF"),
         CTI_window = if_else(str_detect(testing, "c1"), "short", "long")) %>%
  left_join(roi_labels_df, by = join_by(roi == ROI_index)) %>%
  mutate(roi_label = ROI_label) %>%
  separate(ROI_label, into = c("network_amount", "hemisphere", "network_label", "parcel_label", "parcel_extra_nr"), sep = "_", extra = "merge")

### only include parcels from AFC, dlPFC, vlPFC, and PPC

parcels_network_fpn <- c("ContA_IPS", "ContB_IPL","ContA_PFClv", "ContA_PFCl", "ContB_PFClv", "DefaultB_PFCv")
fpn_SupParcels <- c("iPL", "iPL", "dlPF", "dlPF", "aPF", "vlPF")
unique_fpn_SupParcels <- c("iPL", "dlPF", "aPF", "vlPF")
hemispheres <- c("LH", "RH")

schaefer_parcel_to_fpn_sup <- tibble(parcels_network_fpn, fpn_SupParcels)

parcel_retain <- c("IPS", "IPL","PFClv", "PFCl", "PFCv") # control A PFClv belong to dlPFC, control B PFClv belong to AFC

results_schaefer_fpn <- results %>%
  filter(parcel_label %in% parcel_retain) %>%
  mutate(network_parcel = paste(network_label, parcel_label, sep = "_")) %>%
  left_join(schaefer_parcel_to_fpn_sup, by = join_by(network_parcel == parcels_network_fpn))

# run analysis 

t_SupParcels <- c()
p_block <- c()
p_cti <- c()
p_hem <- c()
p_block_by_cti <- c()
p_block_by_hem <- c()
p_cti_by_hem <- c()
p_3way <- c()


for (ii in 1:length(unique_fpn_SupParcels)) { # loop over each super parcel
  
  data_interest <- results_schaefer_fpn %>%
    filter(fpn_SupParcels == unique_fpn_SupParcels[ii]) %>%
    convert_as_factor(subject, block_type, CTI_window, hemisphere, roi)
  
  ## statistical analysis
  inter_lme <- lmer(formula = mean_accuracy ~ block_type * CTI_window * hemisphere + (1 | subject),
                    data = data_interest)
  
  table_summary <- round((summary(inter_lme, correlation= FALSE))$coefficients, 4)
  
  t_SupParcels <- c(t_SupParcels, unique_fpn_SupParcels[ii])
  p_block <- c(p_block, table_summary[2,5])
  p_cti <- c(p_cti, table_summary[3,5])
  p_hem <- c(p_hem, table_summary[4,5])
  p_block_by_cti <- c(p_block_by_cti, table_summary[5,5])
  p_block_by_hem <- c(p_block_by_hem, table_summary[6,5])
  p_cti_by_hem <- c(p_cti_by_hem, table_summary[7,5])
  p_3way <- c(p_3way, table_summary[8,5])
  
  ### visualizing the interaction for each super parcel for each hemisphere separately
  for (iii in 1:length(hemispheres)) {

    hem = hemispheres[iii]

    inter_summary <- data_interest %>%
      filter(hemisphere == hem) %>%
      group_by(block_type, CTI_window) %>%
      summarise(submean = mean(mean_accuracy, na.rm = TRUE),
                subsd = sd(mean_accuracy, na.rm = TRUE),
                se = subsd/((43)^.5)) %>%   # standard error
      ungroup()

    ymin = 0.10
    ymax = 0.15
    chance = 0.1111

    p4 <- ggplot(inter_summary, aes(x = CTI_window, y = submean, group = block_type, color = block_type)) +
      scale_x_discrete(name ="CTI window",
                       limits=c("short", "long"))+
      geom_point(size = 5) +
      geom_line(size = 1) +
      geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                    width = .2) +
      scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                         name = "block type:",
                         breaks = c("RG", "TF")) +
      geom_hline(aes(yintercept=chance),linetype="dashed", size=0.5) +
      coord_cartesian(ylim=c(ymin, ymax)) +
      labs(y = "mean cross-decoding acc") +
      ggtitle(paste0("cross cti decode for ", hemispheres[iii], " ", unique_fpn_SupParcels[ii])) +
      theme_bw()

    ggsave(paste0("cross_cti_decode_SupParcel_", hemispheres[iii], "_", unique_fpn_SupParcels[ii], ".png"), width = 3.9, height = 3.4)
  }
}

cross_cti_decoding_scheafer_fpn_stats <- tibble(t_SupParcels, p_block, p_cti, p_hem, p_block_by_cti, p_block_by_hem, p_cti_by_hem, p_3way)
save(cross_cti_decoding_scheafer_fpn_stats, file = "cross_cti_decoding_lme_results_rois_FPN_Schaefer.Rdata")
load(file = "cross_cti_decoding_lme_results_rois_FPN_Schaefer.Rdata")

#####################################################
#################### FIR—M results ##################----
#####################################################

############# Glasser atlas ###############----
Glasser_SupParcels <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

TR_sec = 1.78 # in second(s)
CTI_mid = 5+1 # in second(s),including the task cue
CTI_longest = 8.75+1 # in second(s),including the task cue

setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/Glasser')
results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/Glasser/decodeAcc_FIR_Conjunc_smthN_beta_rois_Glasser_FPN.csv') %>%
  mutate(task_dim = case_when(str_detect(condition, "stim") ~ "stim",
                              str_detect(condition, "rule") ~ "rule",
                              .default = "conjunc"),
         block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
         TR = as.integer(str_extract(condition, "\\d+$")),
         roi_glas = str_extract(roi, "^[^_]+"),
         hemisphere = str_extract(roi, "(?<=_)[^.]+")) %>%
  left_join(Glasser_SupParcels, by = join_by(roi_glas == fpn_labels))

results_only_fpn <- results %>%
  filter(fpn_SupParcels %in% c("aPF", "dlPF", "iPL"))

TR_check <- results %>%
  group_by(TR, condition) %>%
  summarise(count=n())


## visualizing result --- conjunc decoding
result_conjunc_summary <- results_only_fpn %>%
  group_by(block_type, TR, roi, hemisphere, fpn_SupParcels) %>% # average decoding acc between-subject first for each roi
  summarise(submean = mean(mean_accuracy, na.rm = TRUE),
            subsd = sd(mean_accuracy, na.rm = TRUE),
            se = subsd/((43)^.5)) %>%   # standard error
  ungroup() %>%
  group_by(block_type, TR, hemisphere, fpn_SupParcels) %>% # average across parcels within a super parcel
  summarise(sup_mean = mean(submean, na.rm = TRUE),
            sup_se = mean(se, na.rm = TRUE)) %>%
  ungroup()

(
  p_exp <- ggplot(result_conjunc_summary, aes(x = TR, y = sup_mean, group = block_type, color = block_type)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = sup_mean - sup_se, ymax = sup_mean + sup_se),
                  width = .2) +
    scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                       name = "block type:",
                       breaks = c("RG", "TF")) +
    scale_x_continuous(breaks=seq(0,10,1)) +
    geom_hline(aes(yintercept=0.1111),linetype="dashed", size=0.5) +
    geom_vline(aes(xintercept=CTI_mid/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
    geom_vline(aes(xintercept=CTI_longest/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
    coord_cartesian(ylim=c(0.10, 0.13)) +
    labs(y = "mean decoding acc") +
    ggtitle("conjunc decoding") +
    theme_bw() +
    facet_wrap(vars(fpn_SupParcels, hemisphere), nrow = 3)
)

############# ROIs from univariate contrast ###############----
TR_sec = 1.78 # in second(s)
CTI_mid = 5+1 # in second(s),including the task cue
CTI_longest = 8.75+1 # in second(s),including the task cue

setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/roi_FLM')
results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/decoding/roi_approach/w:o_feat_select/roi_FLM/decodeAcc_FIR_Conjunc_smthN_beta_rois_FLM.csv') %>%
  mutate(task_dim = case_when(str_detect(condition, "stim") ~ "stim",
                              str_detect(condition, "rule") ~ "rule",
                              .default = "conjunc"),
         block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
         TR = as.integer(str_extract(condition, "\\d+$")))

results_fpn <- results %>%
  filter(roi %in% c("Clus_R_AFC.nii", "Clus_L_IFG.nii", "Clus_R_IFG.nii", "Clus_R_IPL.nii"))

## visualizing result --- conjunc decoding
result_conjunc_summary2 <- results_fpn %>%
  group_by(block_type, TR, roi) %>% # average decoding acc between-subject first for each roi
  summarise(submean = mean(mean_accuracy, na.rm = TRUE),
            subsd = sd(mean_accuracy, na.rm = TRUE),
            se = subsd/((43)^.5)) %>%   # standard error
  ungroup()

(
  p_exp <- ggplot(result_conjunc_summary2, aes(x = TR, y = submean, group = block_type, color = block_type)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = submean - se, ymax = submean + se),
                  width = .2) +
    scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .6),
                       name = "block type:",
                       breaks = c("RG", "TF")) +
    scale_x_continuous(breaks=seq(0,10,1)) +
    geom_hline(aes(yintercept=0.1111),linetype="dashed", size=0.5) +
    geom_vline(aes(xintercept=CTI_mid/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
    geom_vline(aes(xintercept=CTI_longest/TR_sec),linetype="dashed", size=0.5, alpha = 0.5) +
    coord_cartesian(ylim=c(0.09, 0.14)) +
    labs(y = "mean decoding acc") +
    ggtitle("conjunc decoding") +
    theme_bw() +
    facet_wrap(vars(roi), nrow = 2)
)


