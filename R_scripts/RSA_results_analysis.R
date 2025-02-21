#################################################
########## Analysis of the RSA Results ##########
#################### GLM-02M ####################
############## FOR MUltivairate MVPA ############
### Author: Mengqiao Chai, chaimengqiao@gmail.com

library(tidyverse)
library(ggh4x)
library(Hmisc)
library(ggprism)
library(rstatix)
library(ggpubr)
library(ez)
library(papaja)

library(lmerTest)
library(lme4)
library(emmeans)

library(flextable)
library(magick)
library(cowplot)
library(wesanderson)
options(contrasts = c("contr.sum","contr.poly"))
options(scipen=999)

#########################################
############# Glasser ROIs ##############
setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA')
parcel_scheme <- 'glasser' # or 'glasser'

## for Glasser atlas
Glasser_SupParcels <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

## for Schaefer atlas
# import all the labels and their corresponding name of the ROI
col_names = c("ROI_index", "ROI_label", "col_a", "col_b", "col_c", "col_d")
roi_labels_df = read.table('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas/schaefer_2018/schaefer_2018/Schaefer2018_400Parcels_17Networks_order.txt',
                           sep='\t', 
                           header=FALSE,col.names = col_names) %>%
  tibble() %>%
  select("ROI_index", "ROI_label")

# include all the super parcel labels
parcels_network_fpn <- c("ContA_IPS", "ContB_IPL","ContA_PFClv", "ContA_PFCl", "ContB_PFClv", "DefaultB_PFCv") 
fpn_SupParcels <- c("iPL", "iPL", "dlPF", "dlPF", "aPF", "vlPF") 
unique_fpn_SupParcels <- c("iPL", "dlPF", "aPF", "vlPF") 
hemispheres <- c("LH", "RH")

schaefer_parcel_to_fpn_sup <- tibble(parcels_network_fpn, fpn_SupParcels)

## Create the big result data frame

if (parcel_scheme == 'glasser') {
  
  results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA/results_RDM_glasser_sorted.csv') %>%
    mutate(block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
           CTI_window = if_else(str_detect(condition, "c1"), "short", "long"),
           bet_or_within = if_else(str_detect(task_relation, "_conjunc"), "within", "between")) %>%
    mutate(
      spear_coef = 1 - distance,
      roi = str_extract(ROI, "^[^_]+"),
      hemisphere = str_extract(ROI, "(?<=_)[^.]+")
    ) %>%
    left_join(Glasser_SupParcels, by = join_by(roi == fpn_labels)) %>%
    convert_as_factor(subject, block_type, CTI_window, bet_or_within, fpn_SupParcels, hemisphere, roi)
  
} else if (parcel_scheme == 'schaefer') {
  
  results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA/results_RDM_schaefer_sorted.csv') %>%
    mutate(roi = ROI,
           block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
           CTI_window = if_else(str_detect(condition, "c1"), "short", "long"),
           bet_or_within = if_else(str_detect(task_relation, "_conjunc"), "within", "between"),
           spear_coef = 1 - distance) %>%
    left_join(roi_labels_df, by = join_by(roi == ROI_index)) %>%
    mutate(roi_label = ROI_label) %>%
    separate(ROI_label, into = c("network_amount", "hemisphere", "network_label", "parcel_label", "parcel_extra_nr"), sep = "_", extra = "merge") %>%
    mutate(network_parcel = paste(network_label, parcel_label, sep = "_")) %>%
    left_join(schaefer_parcel_to_fpn_sup, by = join_by(network_parcel == parcels_network_fpn)) %>%
    convert_as_factor(subject, block_type, CTI_window, bet_or_within, fpn_SupParcels, hemisphere, roi)
}

head(results)
n_sub = length(unique(results$subject))

# saving on object in RData format
save(results, file = "results_RDM_Glasser_sorted_preprocessed.RData")

# check grouping factors
results %>%
  group_by(bet_or_within, task_relation) %>%
  summarise(count = n())

############ First analysis : Compare within-task distance and between-task distance -----
####### hypothesis: within task distance should be smaller than between task distance ########

# mixed effect model (have to be done on HPC)
model_dist <- lmer(formula = distance ~ bet_or_within * fpn_SupParcels * hemisphere * block_type * CTI_window + (1 | subject),
                     data = results,
                     REML = TRUE,
                     control = lmerControl(optimizer = "bobyqa",
                                           calc.derivs = FALSE,
                                           optCtrl = list(maxfun = 2e5)))


bet_wit_lme <- readRDS("/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA/HPC/lme_model.rds")
summary(bet_wit_lme, correlation= FALSE)
emmeans_bet_or_within <- emmeans(bet_wit_lme, ~ bet_or_within)


# compare the distance across super parcels (aPF, dlPF, iPL)
within_between_dist <- results %>%
  filter(fpn_SupParcels %in% c("aPF", "dlPF", "iPL")) %>%
  group_by(hemisphere, fpn_SupParcels, bet_or_within, block_type, CTI_window, subject) %>%
  summarise(meanDis = mean(distance, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(hemisphere, fpn_SupParcels, bet_or_within, block_type, CTI_window) %>%
  summarise(subDis = mean(meanDis, na.rm = TRUE),
            subsd = sd(meanDis, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%   # standard error
  ungroup()
  
(
  p <- ggplot(within_between_dist, aes(x = factor(CTI_window, levels = c("short", "long")), y = subDis, color = block_type, shape = bet_or_within, group = block_type)) +
    geom_point(size = 3, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin = subDis - se, ymax = subDis + se),
                  width = .2,
                  position = position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    theme_apa(base_size = 14) +
    labs(x = "CTI", y = "distance") +
    facet_wrap(vars(fpn_SupParcels, hemisphere), nrow = 3)
)
  
############ Second analysis : Compare distance as a function of overlapping task component -----

############ Third analysis : Compare within-task between-run pattern similarity (equivalent to the opposite of pattern distance) between different conditions ----
load(file = "/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA/results_RDM_Glasser_sorted_preprocessed.RData")

n_sub <- 43

results_within <- results %>%
  filter(bet_or_within == "within", fpn_SupParcels %in% c("aPF", "dlPF", "iPL"))

within_similar <- results %>%
  filter(bet_or_within == "within", fpn_SupParcels %in% c("aPF", "dlPF", "iPL")) %>%
  group_by(fpn_SupParcels, hemisphere, subject, block_type, CTI_window) %>%
  summarise(meanSpearCorr = mean(spear_coef, na.rm = TRUE)) %>% # average within subject
  ungroup() %>%
  group_by(fpn_SupParcels, hemisphere, block_type, CTI_window) %>%
  summarise(subCorr = mean(meanSpearCorr, na.rm = TRUE), # average between subjects
            subsd = sd(meanSpearCorr, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%   # standard error
  ungroup()
  
(
  p3 <- ggplot(within_similar, aes(x = factor(CTI_window, levels = c("short", "long")), y = subCorr, color = block_type, group = block_type)) +
    geom_line(aes(group = block_type), size = 1) +
    geom_pointrange(aes(ymin = subCorr - se, ymax = subCorr + se),
                    shape=15, size = 1, linewidth = 1) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    theme_bw() +
    labs(x = "CTI", y = "Spearman's rho") +
    ggtitle("between-run pattern consistency of the same task") + 
    facet_wrap(vars(fpn_SupParcels, hemisphere), nrow = 3)
)

### plotting for publication

sup_parcels <- c("aPF", "dlPF", "iPL")
hemispheres <- c("left", "right") # "LH" and "RH" for the Schaefer atlas

# for saving plots
setwd("/Users/mengqiao/Documents/fMRI_task_transform/writing-up/reporting/multivariate")
slice_names <- c("aPF_L", "aPF_R", "dlPF_L", "dlPF_R", "iPL_L", "iPL_R")
counter <- 1

# for saving stats
vec_supPar <- c()
vec_hem <- c()

t_intercept <- c()
p_intercept <- c()
F_block <- c()
p_block <- c()
F_cti <- c()
p_cti <- c()
F_inter <- c()
p_inter <- c()


for (i in 1:length(sup_parcels)) {
  for (ii in 1:length(hemispheres)) {
    
    # ## run lme on each hemisphere separately
    # data_hem <- results_within %>%
    #   filter(fpn_SupParcels == sup_parcels[i],
    #          hemisphere == hemispheres[ii]) %>%
    #   group_by(subject, roi, block_type, CTI_window) %>%
    #   summarise(mean_spear_coef = mean(spear_coef, na.rm = TRUE)) %>%
    #   convert_as_factor(subject, block_type, CTI_window, roi)
    # 
    # inter_lme <- lmer(formula = mean_spear_coef ~ block_type * CTI_window + (1 | subject),
    #                   data = data_hem)
    # 
    # table_anova <- round(anova(inter_lme), 4)
    # intercept_summary <- coef(summary(inter_lme))["(Intercept)", ]
    # 
    # vec_supPar <- c(vec_supPar, sup_parcels[i])
    # vec_hem <- c(vec_hem, hemispheres[ii])
    # 
    # t_intercept <- c(t_intercept, round(intercept_summary["t value"],4))
    # p_intercept <- c(p_intercept, round(intercept_summary["Pr(>|t|)"],4))
    # F_block <- c(F_block, table_anova$`F value`[1])
    # p_block <- c(p_block, table_anova$`Pr(>F)`[1])
    # F_cti <- c(F_cti, table_anova$`F value`[2])
    # p_cti <- c(p_cti, table_anova$`Pr(>F)`[2])
    # F_inter <- c(F_inter, table_anova$`F value`[3])
    # p_inter <- c(p_inter, table_anova$`Pr(>F)`[3])
    
    # plotting
    data_aux <- within_similar %>%
      filter(fpn_SupParcels == sup_parcels[i], hemisphere == hemispheres[ii])

    p4 <- ggplot(data_aux, aes(x = factor(CTI_window, levels = c("short", "long")), y = subCorr, color = block_type, group = block_type)) +
      scale_x_discrete(name ="CTI window",
                       limits=c("short", "long"))+
      geom_line(size = 1.3, position = position_dodge(0.08)) +
      geom_pointrange(aes(ymin = subCorr - se, ymax = subCorr + se),
                      shape=22, size = 1.1, linewidth = 1.3, fill = "white", stroke = 1.4, linetype = 1,
                      position = position_dodge(0.08)) +
      scale_color_manual(values = c("#00008B", "#DC143C"),
                         name = "block type:",
                         breaks = c("RG", "TF")) +
      labs(y = "Spearman's rho") +
      coord_cartesian(ylim=c(0.03, 0.22)) +
      theme_half_open() +
      theme(axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            legend.position="none",
            axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11))

    # print(p4)
    slice_path = paste0("/Users/mengqiao/Documents/fMRI_task_transform/writing-up/reporting/multivariate/ROIs/Glasser/", slice_names[counter],".png")
    ggdraw(p4) +
      draw_image(slice_path, x = 1, y = 1, hjust = 1, vjust = 1, width = 0.24, height = 0.24)
    ggsave(paste0("pat_consist_SupParcel_", hemispheres[ii], " ", sup_parcels[i], ".png"), width = 2.6, height = 3)

    counter <- counter + 1
  }
}

pat_consist_glasser_fpn_stats <- tibble(vec_supPar, vec_hem, t_intercept, p_intercept, F_block, p_block, F_cti, p_cti, F_inter, p_inter)
pat_consist_glasser_fpn_stats_table <- flextable(pat_consist_glasser_fpn_stats)
print(pat_consist_glasser_fpn_stats_table)

##### run lme on the overal results

# average across diff kind of between-run distance or consistency(run1-run2, run1-run3, run2-run4....) anda cross diff tasks(9 of them)
data_lme_consist <- results %>%
  filter(bet_or_within == "within", fpn_SupParcels %in% c("aPF", "dlPF", "iPL")) %>%
  group_by(subject, fpn_SupParcels, hemisphere, roi, block_type, CTI_window) %>%
  summarise(mean_spear_coef = mean(spear_coef, na.rm = TRUE),
            mean_dist = mean(distance, na.rm = TRUE),
            count = n()) %>%
  ungroup() %>%
  mutate(CTI_window = factor(CTI_window, levels = c("short", "long")))

# random intercept model
model_consist <- lmer(formula = mean_spear_coef ~ fpn_SupParcels * hemisphere * block_type * CTI_window + (1 | subject),
                   data = data_lme_consist,
                   REML = TRUE,
                   control = lmerControl(optimizer = "bobyqa",
                                         calc.derivs = FALSE,
                                         optCtrl = list(maxfun = 2e5)))

summary(model_consist, correlation= FALSE)
anova(model_consist)

## post-hoc check on the results
marginal_SupParcels <- emmeans(model_consist, "fpn_SupParcels")
pairs(marginal_SupParcels)

marginal_hem <- emmeans(model_consist, "hemisphere")

marginal_block <- emmeans(model_consist, "block_type")

marginal_CTI <- emmeans(model_consist, "CTI_window")

marginal_inter_SupBlock <- emmeans(model_consist, ~ fpn_SupParcels * block_type)
print(marginal_inter_SupBlock)
pairs(marginal_inter_SupBlock, simple = "each")
emmip(model_consist, fpn_SupParcels ~ block_type)

marginal_inter_SupHem <- emmeans(model_consist, ~ fpn_SupParcels * hemisphere)
print(marginal_inter_SupHem)
pairs(marginal_inter_SupHem, simple = "each")
emmip(model_consist, fpn_SupParcels ~ hemisphere)

marginal_inter_HemBlock <- emmeans(model_consist, ~ hemisphere * block_type)
print(marginal_inter_HemBlock)
pairs(marginal_inter_HemBlock, simple = "each")
emmip(model_consist,  block_type ~ hemisphere)

marginal_inter_HemCTI <- emmeans(model_consist, ~ hemisphere * CTI_window)
print(marginal_inter_HemCTI)
pairs(marginal_inter_HemCTI, simple = "each")
emmip(model_consist, hemisphere ~ CTI_window)

emmip(model_consist, hemisphere ~ CTI_window | fpn_SupParcels)

## Post-hoc plotting on the results for publication

# the block type effect
pat_consist_block_stats<- results_within %>%
  convert_as_factor(subject, block_type, CTI_window, hemisphere, fpn_SupParcels) %>%
  group_by(block_type, subject) %>%
  summarise(pat_consist = mean(spear_coef, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(block_type) %>%
  summarise(submean = mean(pat_consist, na.rm = TRUE),
            subsd = sd(pat_consist, na.rm = TRUE),
            se = subsd/((43)^.5)) %>%   # standard error
  ungroup()

(
  p_bar_1 <- ggplot(pat_consist_block_stats, aes(x = block_type, y = submean, color = block_type)) +
    geom_bar(stat="identity", width = 0.5, fill = "white", linewidth = 1) +
    geom_errorbar(aes(ymin=submean-se, ymax=submean+se), width=.1, linewidth=0.5) +
    coord_cartesian(ylim=c(0.08, 0.145)) +
    labs(y = "Spearman's rho", x = "block type") +
    scale_x_discrete(limits=c("RG", "TF"), labels = c("Regular", "Transform")) +
    scale_color_manual(values=c("#00008B", "#DC143C")) + theme_half_open() +
    theme(legend.position = "none")
)

# the CTI effect
pat_consist_cti_stats<- results_within %>%
  convert_as_factor(subject, block_type, CTI_window, hemisphere, fpn_SupParcels) %>%
  group_by(CTI_window, subject) %>%
  summarise(pat_consist = mean(spear_coef, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(CTI_window) %>%
  summarise(submean = mean(pat_consist, na.rm = TRUE),
            subsd = sd(pat_consist, na.rm = TRUE),
            se = subsd/((43)^.5)) %>%   # standard error
  ungroup()

(
  p_bar_2 <- ggplot(pat_consist_cti_stats, aes(x = CTI_window, y = submean, color = CTI_window)) +
    geom_bar(stat="identity", width = 0.5, fill = "white", linewidth = 1) +
    geom_errorbar(aes(ymin=submean-se, ymax=submean+se), width=.1, linewidth=0.5) +
    coord_cartesian(ylim=c(0.08, 0.145)) +
    labs(y = "Spearman's rho", x = "CTI window") +
    scale_x_discrete(limits=c("short", "long")) +
    scale_color_manual(values=wes_palette(n=2, name="Darjeeling1")) + theme_half_open() +
    theme(legend.position = "none")
)

# the hemisphere effect
pat_consist_hem_stats<- results_within %>%
  convert_as_factor(subject, block_type, CTI_window, hemisphere, fpn_SupParcels) %>%
  group_by(hemisphere, subject) %>%
  summarise(pat_consist = mean(spear_coef, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(hemisphere) %>%
  summarise(submean = mean(pat_consist, na.rm = TRUE),
            subsd = sd(pat_consist, na.rm = TRUE),
            se = subsd/((43)^.5)) %>%   # standard error
  ungroup()

(
  p_bar_3 <- ggplot(pat_consist_hem_stats, aes(x = hemisphere, y = submean, color = hemisphere)) +
    geom_bar(stat="identity", width = 0.5, fill = "white", linewidth = 1) +
    geom_errorbar(aes(ymin=submean-se, ymax=submean+se), width=.1, linewidth=0.5) +
    coord_cartesian(ylim=c(0.08, 0.145)) +
    labs(y = "Spearman's rho", x = "hemisphere") +
    scale_x_discrete(limits=c("left", "right")) +
    scale_color_manual(values=wes_palette(n=2, name="Royal1")) + theme_half_open() +
    theme(legend.position = "none")
)

############ Fourth analysis :the correlation between empirical and model RDM ----
setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA')
parcel_scheme <- 'glasser' # or 'glasser'

## for Glasser atlas
Glasser_SupParcels <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

## for Schaefer atlas
# import all the labels and their corresponding name of the ROI
col_names = c("ROI_index", "ROI_label", "col_a", "col_b", "col_c", "col_d")
roi_labels_df = read.table('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/Atlas/schaefer_2018/schaefer_2018/Schaefer2018_400Parcels_17Networks_order.txt',
                           sep='\t', 
                           header=FALSE,col.names = col_names) %>%
  tibble() %>%
  select("ROI_index", "ROI_label")

# include all the super parcel labels
parcels_network_fpn <- c("ContA_IPS", "ContB_IPL","ContA_PFClv", "ContA_PFCl", "ContB_PFClv", "DefaultB_PFCv") 
fpn_SupParcels <- c("iPL", "iPL", "dlPF", "dlPF", "aPF", "vlPF") 
unique_fpn_SupParcels <- c("iPL", "dlPF", "aPF", "vlPF") 
hemispheres <- c("LH", "RH")

schaefer_parcel_to_fpn_sup <- tibble(parcels_network_fpn, fpn_SupParcels)

## Create the big result data frame
model_RDM <- "compo" # or "compo"

if (model_RDM == "conjunc") {
  data_file_name <- 'spear_empirical_conjunc_model_glasser_sorted.csv'
  plot_folder <- "/Users/mengqiao/Documents/fMRI_task_transform/writing-up/reporting/multivariate/conjunc_rsa"
  plot_y_limit <-  c(-0.02, 0.03)
  point_shape <- 23
} else if (model_RDM == "compo") {
  data_file_name <- 'spear_empirical_model_glasser_sorted.csv'
  plot_folder <- "/Users/mengqiao/Documents/fMRI_task_transform/writing-up/reporting/multivariate/compositional_rsa"
  plot_y_limit <- c(-0.04, 0.05)
  point_shape <- 24
} else {
  print("not a valid model RDM")
}

if (parcel_scheme == 'glasser') {
  
  results_model <- read_csv(paste0('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA/', data_file_name)) %>%
    mutate(block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
           CTI_window = if_else(str_detect(condition, "c1"), "short", "long")) %>%
    mutate(
      roi = str_extract(ROI, "^[^_]+"),
      hemisphere = str_extract(ROI, "(?<=_)[^.]+")
    ) %>%
    left_join(Glasser_SupParcels, by = join_by(roi == fpn_labels)) %>%
    convert_as_factor(subject, block_type, CTI_window, fpn_SupParcels, hemisphere, roi)
  
} else if (parcel_scheme == 'schaefer') {
  
  results_model <- read_csv(paste0('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA/', data_file_name)) %>%
    mutate(roi = ROI,
           block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
           CTI_window = if_else(str_detect(condition, "c1"), "short", "long")) %>%
    left_join(roi_labels_df, by = join_by(roi == ROI_index)) %>%
    mutate(roi_label = ROI_label) %>%
    separate(ROI_label, into = c("network_amount", "hemisphere", "network_label", "parcel_label", "parcel_extra_nr"), sep = "_", extra = "merge") %>%
    mutate(network_parcel = paste(network_label, parcel_label, sep = "_")) %>%
    left_join(schaefer_parcel_to_fpn_sup, by = join_by(network_parcel == parcels_network_fpn)) %>%
    convert_as_factor(subject, block_type, CTI_window, fpn_SupParcels, hemisphere, roi)
}

head(results_model)
n_sub = length(unique(results_model$subject))

# check some mapping
inter_mapping <- results_model %>%
  group_by(block_type, CTI_window, condition) %>%
  summarise(count = n())

# plotting the result
results_model_summary <- results_model %>%
  filter(fpn_SupParcels %in% c("aPF", "dlPF", "iPL")) %>%
  group_by(fpn_SupParcels, hemisphere, subject, block_type, CTI_window) %>%
  summarise(meanSpearCorr = mean(spear_rho_empirical_model, na.rm = TRUE)) %>% # average within subject
  ungroup() %>%
  group_by(fpn_SupParcels, hemisphere, block_type, CTI_window) %>%
  summarise(subCorr = mean(meanSpearCorr, na.rm = TRUE), # average between subjects
            subsd = sd(meanSpearCorr, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%   # standard error
  ungroup()

(
  p4 <- ggplot(results_model_summary, aes(x = factor(CTI_window, levels = c("short", "long")), y = subCorr, color = block_type, group = block_type)) +
    geom_line(aes(group = block_type), size = 1) +
    geom_pointrange(aes(ymin = subCorr - se, ymax = subCorr + se),
                    shape = 17, size = 1, linewidth = 1) +
    geom_hline(aes(yintercept=0),linetype="dashed", size=0.5) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    theme_bw() +
    labs(x = "CTI", y = "Spearman's rho") +
    ggtitle("correlation between empirical and model task RDM") + 
    facet_wrap(vars(fpn_SupParcels, hemisphere), nrow = 3)
)

### run analysis per super parcel per hemisphere and plot them separately
sup_parcels <- c("aPF", "dlPF", "iPL")
hemispheres <- c("left", "right")

# for saving plots
slice_names <- c("aPF_L", "aPF_R", "dlPF_L", "dlPF_R", "iPL_L", "iPL_R")
counter <- 1

# for saving the result of analysis
vec_supPar <- c()
vec_hem <- c()

t_intercept <- c()
p_intercept <- c()
F_block <- c()
p_block <- c()
F_cti <- c()
p_cti <- c()
F_inter <- c()
p_inter <- c()

# for saving the plots for publication

for (i in 1:length(sup_parcels)) {
  for (ii in 1:length(hemispheres)) {
    
    ## run lme on each hemisphere separately
    # data_hem <- results_model %>%
    #   filter(fpn_SupParcels == sup_parcels[i],
    #          hemisphere == hemispheres[ii]) %>%
    #   convert_as_factor(subject, block_type, CTI_window, hemisphere, ROI)
    # 
    # inter_lme <- lmer(formula = spear_rho_empirical_model ~ block_type * CTI_window + (1 | subject),
    #                   data = data_hem)
    # 
    # table_anova <- round(anova(inter_lme), 4)
    # intercept_summary <- coef(summary(inter_lme))["(Intercept)", ]
    # 
    # vec_supPar <- c(vec_supPar, sup_parcels[i])
    # vec_hem <- c(vec_hem, hemispheres[ii])
    # 
    # t_intercept <- c(t_intercept, round(intercept_summary["t value"],4))
    # p_intercept <- c(p_intercept, round(intercept_summary["Pr(>|t|)"],4))
    # F_block <- c(F_block, table_anova$`F value`[1])
    # p_block <- c(p_block, table_anova$`Pr(>F)`[1])
    # F_cti <- c(F_cti, table_anova$`F value`[2])
    # p_cti <- c(p_cti, table_anova$`Pr(>F)`[2])
    # F_inter <- c(F_inter, table_anova$`F value`[3])
    # p_inter <- c(p_inter, table_anova$`Pr(>F)`[3])
    
    ## plotting
    data_aux <- results_model_summary %>%
      filter(fpn_SupParcels == sup_parcels[i], hemisphere == hemispheres[ii])

    p4 <- ggplot(data_aux, aes(x = factor(CTI_window, levels = c("short", "long")), y = subCorr, color = block_type, group = block_type)) +
      geom_hline(aes(yintercept=0),linetype="dashed", size=0.5, color = 'grey47') +
      geom_line(size = 1.3, position = position_dodge(0.06)) +
      scale_x_discrete(name ="CTI window",
                       limits=c("short", "long"))+
      geom_pointrange(aes(ymin = subCorr - se, ymax = subCorr + se), position = position_dodge(0.06),
                      shape=point_shape, size = 1.1, linewidth = 1.3, fill = "white", stroke = 1.4, linetype = 1) +
      scale_color_manual(values = c("#00008B", "#DC143C"),
                         name = "block type:",
                         breaks = c("RG", "TF")) +
      labs(y = "Spearman's rho") +
      coord_cartesian(ylim=plot_y_limit) +
      theme_half_open() +
      theme(axis.title.x=element_text(size=11),
            axis.title.y=element_text(size=11),
            legend.position="none",
            axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11))
    
    # print(p4)
    slice_path = paste0("/Users/mengqiao/Documents/fMRI_task_transform/writing-up/reporting/multivariate/ROIs/Glasser/", slice_names[counter],".png")
    ggdraw(p4) +
      draw_image(slice_path, x = 1, y = 1, hjust = 1, vjust = 1, width = 0.24, height = 0.24)
    ggsave(paste0(plot_folder, "/",model_RDM, "_rsa_SupParcel_", hemispheres[ii], " ", sup_parcels[i], ".png"), width = 2.6, height = 3)
    counter <- counter + 1
  }
}

task_rsa_glasser_fpn_stats <- tibble(vec_supPar, vec_hem, t_intercept, p_intercept, F_block, p_block, F_cti, p_cti, F_inter, p_inter)
task_rsa_glasser_fpn_stats_table <- flextable(task_rsa_glasser_fpn_stats)
print(task_rsa_glasser_fpn_stats_table)

##### run lme on the results

# average across diff kind of b
data_lme_compositional <- results_model %>%
  filter(fpn_SupParcels %in% c("aPF", "dlPF", "iPL"))

# random intercept model
model_compositional <- lmer(formula = spear_rho_empirical_model ~ fpn_SupParcels * hemisphere * block_type * CTI_window + (1 | subject),
                      data = data_lme_compositional,
                      REML = TRUE,
                      control = lmerControl(optimizer = "bobyqa",
                                            calc.derivs = FALSE,
                                            optCtrl = list(maxfun = 2e5)))

summary(model_compositional, correlation= FALSE)
anova(model_compositional)
