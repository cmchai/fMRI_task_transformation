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

options(contrasts = c("contr.sum","contr.poly"))
options(scipen=999)

#########################################
############# Glasser ROIs ##############

### for Glasser atlas
setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA')

Glasser_SupParcels <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

results <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA/results_RDM_glasser_sorted.csv') %>%
  mutate(block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
         CTI_window = if_else(str_detect(condition, "c1"), "short", "long"),
         bet_or_within = if_else(str_detect(task_relation, "_conjunc"), "within", "between")) %>%
  mutate(
    spear_coef = 1 - distance,
    roi_glas = str_extract(ROI, "^[^_]+"),
    hemisphere = str_extract(ROI, "(?<=_)[^.]+")
  ) %>%
  left_join(Glasser_SupParcels, by = join_by(roi_glas == fpn_labels)) %>%
  convert_as_factor(subject, block_type, CTI_window, bet_or_within, fpn_SupParcels, hemisphere)

head(results)
n_sub = length(unique(results$subject))

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
                    size = 1, linewidth = 1) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    theme_classic() +
    labs(x = "CTI", y = "Spearman's rho") +
    ggtitle("between-run pattern consistency of the same task") + 
    facet_wrap(vars(fpn_SupParcels, hemisphere), nrow = 3)
)

### plotting for psychonomics

sup_parcels <- c("aPF", "dlPF", "iPL")
hemispheres <- c("left", "right")

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
    
    ## run lme on each hemisphere separately
    data_hem <- results_within %>%
      filter(fpn_SupParcels == sup_parcels[i],
             hemisphere == hemispheres[ii]) %>%
      group_by(subject, roi_glas, block_type, CTI_window) %>%
      summarise(mean_spear_coef = mean(spear_coef, na.rm = TRUE)) %>%
      convert_as_factor(subject, block_type, CTI_window, roi_glas)
    
    inter_lme <- lmer(formula = mean_spear_coef ~ block_type * CTI_window + (1 | subject),
                      data = data_hem)
    
    table_anova <- round(anova(inter_lme), 4)
    intercept_summary <- coef(summary(inter_lme))["(Intercept)", ]
    
    vec_supPar <- c(vec_supPar, sup_parcels[i])
    vec_hem <- c(vec_hem, hemispheres[ii])
    
    t_intercept <- c(t_intercept, round(intercept_summary["t value"],4))
    p_intercept <- c(p_intercept, round(intercept_summary["Pr(>|t|)"],4))
    F_block <- c(F_block, table_anova$`F value`[1])
    p_block <- c(p_block, table_anova$`Pr(>F)`[1])
    F_cti <- c(F_cti, table_anova$`F value`[2])
    p_cti <- c(p_cti, table_anova$`Pr(>F)`[2])
    F_inter <- c(F_inter, table_anova$`F value`[3])
    p_inter <- c(p_inter, table_anova$`Pr(>F)`[3])
    
    # plotting
    # data_aux <- within_similar %>%
    #   filter(fpn_SupParcels == sup_parcels[i], hemisphere == hemispheres[ii])
    # 
    # p4 <- ggplot(data_aux, aes(x = factor(CTI_window, levels = c("short", "long")), y = subCorr, color = block_type, group = block_type)) +
    #   scale_x_discrete(name ="CTI window",
    #                    limits=c("short", "long"))+
    #   geom_pointrange(aes(ymin = subCorr - se, ymax = subCorr + se),
    #                   shape=15, size = 1.4, linewidth = 0.8) +
    #   geom_line(size = 1.3) +
    #   scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .85),
    #                      name = "block type:",
    #                      breaks = c("RG", "TF")) +
    #   coord_cartesian(ylim=c(0.03, 0.21)) +
    #   theme_bw() +
    #   theme(axis.title.x=element_blank(),
    #         axis.title.y=element_blank(),
    #         legend.position="none",
    #         axis.text.x = element_text(size=14),
    #         axis.text.y = element_text(size=14))
    # # print(p4)
    # ggsave(paste0("pat_consist_SupParcel_", hemispheres[ii], " ", sup_parcels[i], "_psychonomics.png"), width = 2.6, height = 3)
  }
}

pat_consist_glasser_fpn_stats <- tibble(vec_supPar, vec_hem, t_intercept, p_intercept, F_block, p_block, F_cti, p_cti, F_inter, p_inter)
pat_consist_glasser_fpn_stats_table <- flextable(pat_consist_glasser_fpn_stats)
print(pat_consist_glasser_fpn_stats_table)

##### run lme on the overal results

# average across diff kind of between-run distance or consistency(run1-run2, run1-run3, run2-run4....) anda cross diff tasks(9 of them)
data_lme_consist <- results %>%
  filter(bet_or_within == "within", fpn_SupParcels %in% c("aPF", "dlPF", "iPL")) %>%
  group_by(subject, fpn_SupParcels, hemisphere, roi_glas, block_type, CTI_window) %>%
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

# post-hoc check on the results
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
emmip(model_consist, hemisphere ~ block_type)

marginal_inter_HemCTI <- emmeans(model_consist, ~ hemisphere * CTI_window)
print(marginal_inter_HemCTI)
pairs(marginal_inter_HemCTI, simple = "each")
emmip(model_consist, hemisphere ~ CTI_window)

emmip(model_consist, hemisphere ~ CTI_window | fpn_SupParcels)

############ Fourth analysis :the correlation between empirical and model RDM ----

### for Glasser atlas
setwd('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA')

Glasser_SupParcels <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/resources/fmri-extract-HCP-mask-main/fpn_SupParcels.csv')

results_model <- read_csv('/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/RSA/spear_empirical_model_glasser_sorted.csv') %>%
      mutate(block_type = if_else(str_detect(condition, "RG-"), "RG", "TF"),
             CTI_window = if_else(str_detect(condition, "c1"), "short", "long")) %>%
      mutate(
        roi_glas = str_extract(ROI, "^[^_]+"),
        hemisphere = str_extract(ROI, "(?<=_)[^.]+")
      ) %>%
      left_join(Glasser_SupParcels, by = join_by(roi_glas == fpn_labels)) %>%
      convert_as_factor(subject, block_type, CTI_window, fpn_SupParcels, hemisphere)

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
                    size = 1, linewidth = 1) +
    geom_hline(aes(yintercept=0),linetype="dashed", size=0.5) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    theme_classic() +
    labs(x = "CTI", y = "Spearman's rho") +
    ggtitle("correlation between empirical and model task RDM") + 
    facet_wrap(vars(fpn_SupParcels, hemisphere), nrow = 3)
)

### run analysis per super parcel per hemisphere and plot them separately
sup_parcels <- c("aPF", "dlPF", "iPL")
hemispheres <- c("left", "right")

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
    
    ## run lme on each hemisphere separately
    data_hem <- results_model %>%
      filter(fpn_SupParcels == sup_parcels[i],
             hemisphere == hemispheres[ii]) %>%
      convert_as_factor(subject, block_type, CTI_window, hemisphere, roi_glas)
    
    inter_lme <- lmer(formula = spear_rho_empirical_model ~ block_type * CTI_window + (1 | subject),
                      data = data_hem)
    
    table_anova <- round(anova(inter_lme), 4)
    intercept_summary <- coef(summary(inter_lme))["(Intercept)", ]
    
    vec_supPar <- c(vec_supPar, sup_parcels[i])
    vec_hem <- c(vec_hem, hemispheres[ii])
    
    t_intercept <- c(t_intercept, round(intercept_summary["t value"],4))
    p_intercept <- c(p_intercept, round(intercept_summary["Pr(>|t|)"],4))
    F_block <- c(F_block, table_anova$`F value`[1])
    p_block <- c(p_block, table_anova$`Pr(>F)`[1])
    F_cti <- c(F_cti, table_anova$`F value`[2])
    p_cti <- c(p_cti, table_anova$`Pr(>F)`[2])
    F_inter <- c(F_inter, table_anova$`F value`[3])
    p_inter <- c(p_inter, table_anova$`Pr(>F)`[3])
    
    ## plotting 
    data_aux <- results_model_summary %>%
      filter(fpn_SupParcels == sup_parcels[i], hemisphere == hemispheres[ii])
    
    p4 <- ggplot(data_aux, aes(x = factor(CTI_window, levels = c("short", "long")), y = subCorr, color = block_type, group = block_type)) +
      scale_x_discrete(name ="CTI window",
                       limits=c("short", "long"))+
      geom_pointrange(aes(ymin = subCorr - se, ymax = subCorr + se),
                      shape=17, size = 1.4, linewidth = 0.8) +
      geom_hline(aes(yintercept=0),linetype="dashed", size=0.5) +
      geom_line(size = 1.3) +
      scale_color_manual(values = alpha(c("#4E84C4", "#FC4E07"), .85),
                         name = "block type:",
                         breaks = c("RG", "TF")) +
      coord_cartesian(ylim=c(-0.04, 0.05)) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14))
    # print(p4)
    ggsave(paste0("task_rsa_SupParcel_", hemispheres[ii], " ", sup_parcels[i], ".png"), width = 2.6, height = 3)
  }
}

task_rsa_glasser_fpn_stats <- tibble(vec_supPar, vec_hem, t_intercept, p_intercept, F_block, p_block, F_cti, p_cti, F_inter, p_inter)
task_rsa_glasser_fpn_stats_table <- flextable(task_rsa_glasser_fpn_stats)
print(task_rsa_glasser_fpn_stats_table)

##### run lme on the results

# average across diff kind of between-run distance or consistency(run1-run2, run1-run3, run2-run4....) anda cross diff tasks(9 of them)
data_lme_consist <- results %>%
  filter(bet_or_within == "within", fpn_SupParcels %in% c("aPF", "dlPF", "iPL")) %>%
  group_by(subject, fpn_SupParcels, hemisphere, roi_glas, block_type, CTI_window) %>%
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
