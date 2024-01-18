################# statistical analysis ########################
############ behavioral analysis of fMRI experiment ###########
library(here)
library(tidyverse)
library(ez)
library(feather)
library(rstatix)
library(wesanderson)
library(ggsignif)
library(ggnewscale)
library(papaja)
library(rempsyc)

library(lmerTest)
library(lme4)
library(emmeans)

options(contrasts = c("contr.sum","contr.poly"))
options(scipen=999)

scale_this <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

###########################################################
########### the Analysis on the Reaction time #############
###########################################################
load(file = paste0(here(), "/results/bigdf_fmri_clean.Rdata"))
glimpse(bigdf_clean)
n_sub <- length(unique(bigdf_clean$subject))
# 
# start_block_df <- bigdf_clean %>%
#   filter(run_id == "run_1") %>%
#   group_by(subject, block_id) %>%
#   summarise(count = n()) %>%
#   mutate(start_block = str_sub(block_id, 1, 2)) %>%
#   select(-c(block_id, count)) %>%
#   ungroup()
# 
# bigdf_clean_aux <- bigdf_clean %>%
#   select(-start_block) %>%
#   left_join(start_block_df, by = "subject") %>%
#   relocate(start_block, .after = run_id)
# 
# bigdf_clean <- bigdf_clean_aux
# save(bigdf_clean, file = paste0(here(), "/results/bigdf_fmri_clean.Rdata"))
# save(bigdf_clean_error, file = paste0(here(), "/results/bigdf_fmri_clean_error.Rdata"))

length(unique(bigdf_clean$subject))

bigdf_clean <- bigdf_clean %>%
  mutate(rt = rt*1000)

start_block_balance <- bigdf_clean %>%
  group_by(start_block, subject) %>%
  summarise(count_trials = n()) %>%
  ungroup() %>%
  group_by(start_block) %>%
  summarise(count_subj = n())

### step1. only include the regular trials for analysis ###
regular_trials <- bigdf_clean %>%
  filter(trial_nature == "RG") %>%
  select(CTI_length, rt, start_block, block_type, block_id, stim, rule, CTI, resp_map_age,
         resp_map_size, resp_map_location, subject, accuracy, incongruence_score_stim, incongruence_score_rule)

### step1. only include the transform trials for analysis ###
transform_trials <- bigdf_clean %>%
  filter(trial_nature == "TF") %>%
  select(CTI_length, rt, start_block, block_type, block_id, stim, rule, CTI, resp_map_age,
         resp_map_size, resp_map_location, subject, accuracy, incongruence_score_stim, incongruence_score_rule)

###########################################################
################## the Analysis on the RT #################
###########################################################

# - RT summary and plotting of regular trials ========
rt_descrip <- regular_trials %>%
  convert_as_factor(CTI_length, subject, start_block, block_type) %>%
  group_by(subject, start_block, block_type, CTI_length) %>%
  summarise(count = n(),
            meanCTI = mean(CTI),
            meanRT = mean(rt, na.rm = TRUE)) %>%
  ungroup()


# -- subplots for each subject =======================
(
  p_indiv <- ggplot(rt_descrip, aes(x = CTI_length, y = meanRT, color = block_type)) +
    geom_point(size = 2.5, position = position_dodge(0.2)) +
    geom_line(aes(group = block_type), 
              size = 0.8, position = position_dodge(0.2), alpha = .5) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw() +
    facet_wrap(~ subject, nrow = 5)
)

# -- subplot for each participant depend on the starting block type
(
  p_indiv_RGst <- ggplot(subset(rt_descrip, start_block == "RG"), aes(x = CTI_length, y = meanRT, color = block_type)) +
    geom_point(size = 2.5, position = position_dodge(0.2)) +
    geom_line(aes(group = block_type), 
              size = 0.8, position = position_dodge(0.2), alpha = .5) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw() +
    facet_wrap(~ subject, nrow = 5)
)

(
  p_indiv_TFst <- ggplot(subset(rt_descrip, start_block == "TF"), aes(x = CTI_length, y = meanRT, color = block_type)) +
    geom_point(size = 2.5, position = position_dodge(0.2)) +
    geom_line(aes(group = block_type), 
              size = 0.8, position = position_dodge(0.2), alpha = .5) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw() +
    facet_wrap(~ subject, nrow = 5)
)

# -- group level plot ================================
rt_descrip2 <- rt_descrip %>%
  group_by(block_type, CTI_length) %>%
  summarise(submeanCTI = mean(meanCTI),
            subRT = mean(meanRT, na.rm = TRUE),
            subsd = sd(meanRT, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>% # standard error
  ungroup()

group_vars(rt_descrip) # need to check the grouping variables after using "summarize" since it will drop one grouping variable 

(
  p_groupRT <- ggplot(rt_descrip2, aes(x = CTI_length, y = subRT, color = block_type)) +
    geom_point(size = 7, position = position_dodge(0.1)) +
    geom_line(aes(group = block_type), size = 1.5,
              position = position_dodge(0.1),
              linetype="dotted") +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .1, size = 1,
                  position = position_dodge(0.1)) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw()
)

# group level plot depends on the starting block type

# Start with RG block
rt_descrip2_RGsta <- rt_descrip %>%
  filter(start_block == "RG") %>%
  group_by(block_type, CTI_length) %>%
  summarise(submeanCTI = mean(meanCTI),
            subRT = mean(meanRT, na.rm = TRUE),
            subsd = sd(meanRT, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>% # standard error
  ungroup()

(
  p_groupRT_RGst <- ggplot(rt_descrip2_RGsta, aes(x = CTI_length, y = subRT, color = block_type)) +
    geom_point(size = 7, position = position_dodge(0.1)) +
    geom_line(aes(group = block_type), size = 1.5,
              position = position_dodge(0.1),
              linetype="dotted") +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .1, size = 1,
                  position = position_dodge(0.1)) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw()
)

# Start with TF block
rt_descrip2_TFsta <- rt_descrip %>%
  filter(start_block == "TF") %>%
  group_by(block_type, CTI_length) %>%
  summarise(submeanCTI = mean(meanCTI),
            subRT = mean(meanRT, na.rm = TRUE),
            subsd = sd(meanRT, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>% # standard error
  ungroup()

(
  p_groupRT_TFst <- ggplot(rt_descrip2_TFsta, aes(x = CTI_length, y = subRT, color = block_type)) +
    geom_point(size = 7, position = position_dodge(0.1)) +
    geom_line(aes(group = block_type), size = 1.5,
              position = position_dodge(0.1),
              linetype="dotted") +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .1, size = 1,
                  position = position_dodge(0.1)) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw()
)

# -- plot both individual and group level ============
(
  p <- ggplot(rt_descrip, aes(x = CTI_length, y = meanRT, color = block_type)) +
    geom_point(size = 4, alpha = .1) +
    geom_line(aes(group = interaction(subject, block_type)),
              size = 2, alpha = .1) +
    geom_point(data = rt_descrip2, aes(y = subRT), size = 5) +
    geom_line(data = rt_descrip2, aes(y = subRT, group = block_type), size = 1) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw()
)

# -- visualize the data in another way using CTI instead of CTI bins #######
rt_descrip3 <- regular_trials %>%
  convert_as_factor(subject, block_type) %>%
  group_by(subject, block_type, CTI) %>%
  summarise(count = n(),
            meanRT = mean(rt, na.rm = TRUE)) %>%
  ungroup()

rt_descrip4 <- rt_descrip3 %>%
  group_by(block_type, CTI) %>%
  summarise(subRT = mean(meanRT, na.rm = TRUE),
            subsd = sd(meanRT, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%
  ungroup()  

# plot individually
(
  p <- ggplot(rt_descrip3, aes(x = CTI, y = meanRT, color = block_type)) +
    geom_point(size = 1,position=position_dodge(0.2)) +
    geom_line(aes(group = block_type), size = 0.5, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw() +
    facet_wrap(~ subject, nrow = 4)
)

# plot group level
(
  p <- ggplot(rt_descrip4, aes(x = CTI, y = subRT, color = block_type)) +
    geom_point(size = 5,position=position_dodge(0.2)) +
    geom_line(aes(group = block_type), size = 1) +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .1, size = 1,
                  position = position_dodge(0.2)) +
    geom_vline(aes(xintercept=5000),linetype="dashed", size=1) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    coord_cartesian(ylim=c(1050, 1300)) +
    theme_bw()
)

# add an smooth line for publication
(
  p1 <- ggplot(rt_descrip4, aes(x = CTI, y = subRT, color = block_type)) +
    geom_point(size = 5,position=position_dodge(0.2), alpha = .3) +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .1, size = 1,
                  position = position_dodge(0.2), alpha = .3) +
    geom_smooth(method = "loess", span = 0.75, size = 2, linetype="dashed") +
    geom_vline(aes(xintercept=5000),linetype="dashed", size=1) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    coord_cartesian(ylim=c(1050, 1300)) +
    theme_apa(base_size = 14) +
    theme(axis.title.x = element_blank()) +
    labs(x = "CTI in ms", y = "Reaction time in ms") +
    ggtitle("fMRI Experiment")
)

# -- ANOVA testing =======================
rt_reg.aov <- anova_test(
  data = rt_descrip,
  dv = meanRT,
  wid = subject,
  within = c(block_type, CTI_length)
)

get_anova_table(rt_reg.aov)

# - RT congruency effect ============================
# -- stimulus congruency effect ======================
rt_cgc_stim <- regular_trials %>%
  convert_as_factor(subject, block_type) %>%
  group_by(subject, block_type, incongruence_score_stim) %>%
  summarise(meanRT = mean(rt, na.rm = TRUE)) %>%
  ungroup()

rt_cgc_stim2 <- rt_cgc_stim%>%
  group_by(block_type, incongruence_score_stim) %>%
  summarise(subRT = mean(meanRT, na.rm = TRUE)) %>%
  ungroup()

(
  p <- ggplot(rt_cgc_stim, aes(x = incongruence_score_stim, y = meanRT, color = block_type)) +
    geom_point(size = 4, position = position_dodge(0.2), alpha = .08) +
    geom_line(aes(group = interaction(subject, block_type)), 
              size = 1, position = position_dodge(0.2), alpha = .2) +
    geom_point(data = rt_cgc_stim2, aes(y = subRT), size = 7,
               position=position_dodge(0.2)) +
    geom_line(data = rt_cgc_stim2, aes(y = subRT, group = block_type), 
              size = 1.6, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw()
)
# -- task rule congruency effect ======================
rt_cgc_rule <- regular_trials %>%
  convert_as_factor(subject, block_type) %>%
  group_by(subject, block_type, incongruence_score_rule) %>%
  summarise(meanRT = mean(rt, na.rm = TRUE)) %>%
  ungroup()

rt_cgc_rule2 <- rt_cgc_rule %>%
  group_by(block_type, incongruence_score_rule) %>%
  summarise(subRT = mean(meanRT, na.rm = TRUE)) %>%
  ungroup()  

(
  p <- ggplot(rt_cgc_rule, aes(x = incongruence_score_rule, y = meanRT, color = block_type)) +
    geom_point(size = 4, position = position_dodge(0.2), alpha = .08) +
    geom_line(aes(group = interaction(subject, block_type)), 
              size = 1, position = position_dodge(0.2), alpha = .2) +
    geom_point(data = rt_cgc_rule2, aes(y = subRT), size = 7) +
    geom_line(data = rt_cgc_rule2, aes(y = subRT, group = block_type), 
              size = 1.6, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw()
)

# - Mixed Effect Model on RT =========================
regular_trials_lme <- regular_trials %>%
  mutate(CTI_sec = CTI/1000) %>%
  convert_as_factor(block_type, subject, CTI_length) %>%
  ungroup()

contrasts(regular_trials_lme$block_type)
contrasts(regular_trials_lme$CTI_length)

## specify the model using continuous CTI ---- 
## - random intercept 
model_conRT <- lmer(formula = rt ~ (CTI_sec + incongruence_score_stim + incongruence_score_rule) * block_type
                 + (1 | subject),
                 data = regular_trials_lme)

anova(model_conRT)
summary(model_conRT, correlation= FALSE)

## - random slope (maximal random structure)
model_conRT_max <- lmer(formula = rt ~ (CTI_sec + incongruence_score_stim + incongruence_score_rule) * block_type
                     + ((CTI_sec + incongruence_score_stim + incongruence_score_rule) * block_type | subject),
                     data = regular_trials_lme,
                     REML = TRUE,
                     control = lmerControl(optimizer = "nloptwrap",
                                           calc.derivs = FALSE,
                                           optCtrl = list(maxfun = 2e5)))

anova(model_conRT_max)
summary(model_conRT_max, correlation= FALSE)

# specify the model with using CTI length(short vs. long) ----
## - random intercept 
model_facRT <- lmer(formula = rt ~ (CTI_length + incongruence_score_stim + incongruence_score_rule) * block_type
                        + (1 | subject),
                        data = regular_trials_lme)

anova(model_facRT)
summary(model_facRT, correlation= FALSE)

regular_trials_lme$fitted_RT_Bi <- fitted(model_RT2)

## - random slope (maximal random structure)
model_RT_max <- lmer(formula = rt ~ (CTI_length + incongruence_score_stim + incongruence_score_rule) * block_type
                  + ((CTI_length + incongruence_score_stim + incongruence_score_rule) * block_type | subject),
                  data = regular_trials_lme,
                  REML = TRUE,
                  control = lmerControl(optimizer = "nloptwrap",
                                        calc.derivs = FALSE,
                                        optCtrl = list(maxfun = 2e5)))

anova(model_RT_max)
summary(model_RT_max, correlation= FALSE)

# plot the fitted value
rt_fitted_descrip2 <- regular_trials_lme %>%
  group_by(block_type, CTI_length) %>%
  summarise(subRT_fit = mean(fitted_RT_Bi, na.rm = TRUE)) %>%
  ungroup()

(
  p <- ggplot(rt_fitted_descrip2, aes(x = CTI_length, y = subRT_fit, color = block_type)) +
    geom_point(size = 7, position = position_dodge(0.1)) +
    geom_line(aes(group = block_type), size = 2,
              position = position_dodge(0.1),
              linetype=2) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw()
)

#########################################################
# - analyzing RT on transform trials =======================
#########################################################
rt_tran_descrip1 <- transform_trials %>%
  convert_as_factor(subject) %>%
  group_by(subject, CTI) %>%
  summarise(count = n(),
            meanRT = mean(rt, na.rm = TRUE)) %>%
  ungroup()

rt_tran_descrip2 <- rt_tran_descrip1 %>%
  group_by(CTI) %>%
  summarise(subRT = mean(meanRT, na.rm = TRUE),
            subsd = sd(meanRT, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%
  ungroup()

(
  p <- ggplot(rt_tran_descrip2, aes(x = CTI, y = subRT)) +
    geom_point(size = 5,position=position_dodge(0.2), color = "#FC4E07") +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .1, size = 1,
                  position = position_dodge(0.2), alpha = .3) +
    geom_line(size = 1, alpha = .3)  +
    coord_cartesian(ylim=c(950, 1200)) +
    theme_bw()
)

# - Mixed Effect Model on transform trial RT =========================
tran_trials_lme <- transform_trials %>%
  mutate(CTI_sec = CTI/1000) %>%
  convert_as_factor(subject)

## specify the model using continuous CTI ---- 
## - random intercept 
model_tran_conRT <- lmer(formula = rt ~ CTI_sec + (CTI_sec | subject),
                    data = tran_trials_lme,
                    REML = TRUE,
                    control = lmerControl(optimizer = "nloptwrap",
                                          calc.derivs = FALSE,
                                          optCtrl = list(maxfun = 2e5)))

anova(model_tran_conRT)
summary(model_tran_conRT, correlation= FALSE)

#### the across-subject correlation between block-type RT diff in the regular trials and the RT of transform trials

rt_regular_long_diff <- regular_trials %>% 
  filter(CTI_length == "long") %>%
  group_by(subject, block_type) %>%
  summarise(meanRT = mean(rt, na.rm = TRUE)) %>%
  pivot_wider(names_from = block_type, values_from = meanRT) %>%
  mutate(long_rg_diff = TF - RG) %>% # should be > 0 if people for prepare for transform, the larger it is, the more prepared  ppl are
  ungroup() %>%
  select(subject, long_rg_diff)

subject_TFshort_mean <- regular_trials %>%
  filter(block_type == "TF", CTI_length == "short") %>%
  group_by(subject) %>%
  summarise(TF_reg_short = mean(rt, na.rm = TRUE)) %>% # which serve as baseline for the tranform trial RT
  ungroup()

subject_tranTrial_mean <- transform_trials %>%
  group_by(subject) %>%
  summarise(meanRT_tran = mean(rt, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(subject_TFshort_mean, by = "subject") %>%
  mutate(tranRT = meanRT_tran - TF_reg_short) %>% # the smaller it is, the better ppl perform on the transform task
  left_join(rt_regular_long_diff, by = "subject")
  
corr1 <- cor.test(subject_tranTrial_mean$tranRT, subject_tranTrial_mean$long_rg_diff)

## visualize this trend
corr_coef <- round(corr1$estimate, 3)
corr_p <- round(corr1$p.value, 3)

(
p_corr <- ggplot(subject_tranTrial_mean, aes(x = tranRT, y = long_rg_diff)) +
  geom_point(size = 3, alpha = .3) +
  geom_smooth(method = "lm", color = "black") +
  theme_classic() +
  labs(x = "RT of transform trials(compared to short reg trials in TF blocks)", y = "The RT diff. of long trials between TF and RG blocks") +
  annotate("text", x=0.5, y=120, label=paste0("r = ", corr_coef), size = 6) +
  annotate("text", x=0.5, y=95, label=paste0("italic(p) == ", corr_p), size = 6, parse = TRUE)
)

###########################################################
########### the Analysis on the Accuracy rate #############
###########################################################

load(file = paste0(here(), "/results/bigdf_fmri_clean_error.Rdata"))
bigdf_clean_error <- bigdf_clean_error %>%
  filter(!is.na(rt)) # get rid of missing trials before analyzing accuracy

# bigdf_clean_error_aux <- bigdf_clean_error %>%
#   select(-start_block) %>%
#   left_join(start_block_df, by = "subject") %>%
#   relocate(start_block, .after = run_id)
# 
# bigdf_clean_error <- bigdf_clean_error_aux

### step1. only include the regular trials for analysis ###
regular_trials_error <- bigdf_clean_error %>%
  filter(trial_nature == "RG") %>%
  select(rt, start_block, block_type, block_id, stim, rule, CTI, resp_map_age,
         resp_map_size, resp_map_location, subject, accuracy, CTI_length,
         incongruence_score_stim, incongruence_score_rule)

tran_trials_error <- bigdf_clean_error %>%
  filter(trial_nature == "TF") %>%
  select(rt, start_block, block_type, block_id, stim, rule, CTI, resp_map_age,
         resp_map_size, resp_map_location, subject, accuracy, CTI_length,
         incongruence_score_stim, incongruence_score_rule) %>%
  mutate(ACC_bi = if_else(accuracy == TRUE, 1, 0))

### step2. compare the subject level accuracy rate between regular and transform block
accRate_descrip <- regular_trials_error %>%
  group_by(subject, block_type, CTI_length) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

accRate_descrip2 <- accRate_descrip %>%
  group_by(block_type, CTI_length) %>%
  summarise(subaccRate = mean(accRate, na.rm = TRUE),
            subsd = sd(accRate, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%
  ungroup()

(
  p <- ggplot(accRate_descrip2, aes(x = CTI_length, y = subaccRate, color = block_type)) +
    geom_point(size = 5, position=position_dodge(0.2)) +
    geom_line(aes(group = block_type), size = 1, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin = subaccRate - se, ymax = subaccRate + se),
                  width = .1, size = 1,
                  position = position_dodge(0.2)) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    coord_cartesian(ylim=c(0.87, 0.915)) +
    labs(y = "accuracy") +
    theme_bw()
)

# -- visualize the data in another way using CTI instead of CTI length #######
accRate_descrip3 <- regular_trials_error %>%
  group_by(subject, block_type, CTI) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

accRate_descrip4 <- accRate_descrip3 %>%
  group_by(block_type, CTI) %>%
  summarise(subaccRate = mean(accRate, na.rm = TRUE),
            subsd = sd(accRate, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%
  ungroup() 

# plot individually
(
  p <- ggplot(accRate_descrip3, aes(x = CTI, y = accRate, color = block_type)) +
    geom_point(size = 1,position=position_dodge(0.2)) +
    geom_line(aes(group = block_type), size = 0.5, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw() +
    facet_wrap(~ subject, nrow = 4)
)

# plot group level
(
  p <- ggplot(accRate_descrip4, aes(x = CTI, y = subaccRate, color = block_type)) +
    geom_point(size = 5,position=position_dodge(0.2)) +
    geom_line(aes(group = block_type), size = 1) +
    geom_errorbar(aes(ymin = subaccRate - se, ymax = subaccRate + se),
                  width = .1, size = 1,
                  position = position_dodge(0.2)) +
    geom_vline(aes(xintercept=5000),linetype="dashed", size=1) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    coord_cartesian(ylim=c(0.75, 1)) +
    theme_bw()
)

# add an smooth line for publication
(
  p1 <- ggplot(accRate_descrip4, aes(x = CTI, y = subaccRate, color = block_type)) +
    geom_point(size = 5,position=position_dodge(0.2), alpha = .3) +
    geom_smooth(method = "loess", span = 0.75, size = 2, linetype="dashed") +
    geom_vline(aes(xintercept=5000),linetype="dashed", size=1) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    coord_cartesian(ylim=c(0.78, 1)) +
    theme_apa(base_size = 14) +
    labs(x = "CTI in ms", y = "accuracy") +
    ggtitle("fMRI Experiment")
)


### the ANOVA analysis on ACC ----
acc_reg.aov <- anova_test(
  data = accRate_descrip,
  dv = accRate,
  wid = subject,
  within = c(block_type, CTI_length)
)

get_anova_table(acc_reg.aov)

### the mixed effect logistic regression model on ACC ----
options(contrasts = c("contr.sum","contr.poly"))
options(scipen = FALSE)

regular_trials_error_lme <- regular_trials_error %>%
  mutate(CTI_sec = CTI/1000,
         ACC_bi = if_else(accuracy == TRUE, 1, 0)) %>%
  convert_as_factor(block_type, subject) %>%
  ungroup()

contrasts(regular_trials_error_lme$block_type)

### -- random intercept model
model_acc <- glmer(formula = ACC_bi ~ (CTI_sec + incongruence_score_stim + incongruence_score_rule) * block_type 
                   + (1|subject),
                   data = regular_trials_error_lme,
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa"))

anova(model_acc)
summary(model_acc, correlation= FALSE)

### -- maximal model
model_ac_max <- glmer(formula = ACC_bi ~ (CTI_sec + incongruence_score_stim + incongruence_score_rule) * block_type 
                    + ((CTI_sec + incongruence_score_stim + incongruence_score_rule) * block_type | subject),
                    data = regular_trials_error_lme,
                    family = binomial,
                    REML = TRUE,
                    control = glmerControl(optimizer = "nloptwrap",
                                         calc.derivs = FALSE,
                                         optCtrl = list(maxfun = 2e5)))

anova(model_acc_max)
summary(model_ac_max, correlation= FALSE)
saveRDS(model_ac_max, file = paste0(here(), "/results/lme_max_acc.rda"))

emtrends(model_ac_max, ~ block_type, var = 'incongruence_score_stim')

############################################################
# - Analyzing accuracy on transform trials =================
############################################################
acc_tran_descrip1 <- tran_trials_error %>%
  convert_as_factor(subject) %>%
  group_by(subject, CTI) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

acc_tran_descrip2 <- acc_tran_descrip1 %>%
  group_by(CTI) %>%
  summarise(subaccRate = mean(accRate, na.rm = TRUE),
            subsd = sd(accRate, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%
  ungroup()

(
  p <- ggplot(acc_tran_descrip2, aes(x = CTI, y = subaccRate)) +
    geom_point(size = 5,position=position_dodge(0.2), color = "#FC4E07") +
    geom_errorbar(aes(ymin = subaccRate - se, ymax = subaccRate + se),
                  width = .1, size = 1,
                  position = position_dodge(0.2), alpha = .3) +
    geom_line(size = 1, alpha = .3)  +
    coord_cartesian(ylim=c(0.80, 1)) +
    theme_bw()
)

# - Mixed Effect Model on transform trial RT =====================
tran_trials_lme <- transform_trials %>%
  mutate(CTI_sec = CTI/1000) %>%
  convert_as_factor(subject)

## specify the model using continuous CTI ---- 
## - random intercept 
model_tran_conRT <- lmer(formula = rt ~ CTI_sec + (CTI_sec | subject),
                         data = tran_trials_lme,
                         REML = TRUE,
                         control = lmerControl(optimizer = "nloptwrap",
                                               calc.derivs = FALSE,
                                               optCtrl = list(maxfun = 2e5)))

anova(model_tran_conRT)
summary(model_tran_conRT, correlation= FALSE)

### logistic regression on ACC of Transform trials ###
tran_trials_error_lme <- tran_trials_error %>%
  mutate(CTI_sec = CTI/1000) %>%
  convert_as_factor(subject)

model_tranAcc <- glmer(formula = ACC_bi ~ CTI_sec + (CTI_sec|subject),
                       data = tran_trials_error_lme,
                       family = binomial,
                       control = glmerControl(optimizer = "bobyqa",
                                              calc.derivs = FALSE,
                                              optCtrl = list(maxfun = 2e5)))

summary(model_tranAcc, correlation = FALSE)


### Explore the between-subject correlation between the RT main effect and the ACC in transform trials
acc_subject_TFshort_mean <- regular_trials_error %>%
  filter(block_type == "TF", CTI_length == "short") %>%
  group_by(subject) %>%
  summarise(acc_TF_reg_short = sum(accuracy)/length(accuracy)) %>% # which serve as baseline for the Transform trial accuracy
  ungroup()

acc_subject_tranTrial_mean <- tran_trials_error %>%
  group_by(subject) %>%
  summarise(subaccRate_tran = sum(accuracy)/length(accuracy)) %>%
  ungroup() %>%
  left_join(acc_subject_TFshort_mean, by = "subject") %>%
  mutate(tranACC = subaccRate_tran - acc_TF_reg_short) %>% # the higher it is, the better ppl perform on the transform task
  left_join(rt_regular_long_diff, by = "subject")

corr2 <- cor.test(acc_subject_tranTrial_mean$subaccRate_tran, acc_subject_tranTrial_mean$long_rg_diff)

## visualize this trend
corr_coef2 <- round(corr2$estimate, 3)
corr_p2 <- round(corr2$p.value, 3)

(
  p_corr2 <- ggplot(acc_subject_tranTrial_mean, aes(x = tranACC, y = long_rg_diff)) +
    geom_point(size = 3, alpha = .3) +
    geom_smooth(method = "lm", color = "black") +
    theme_classic() +
    labs(x = "ACC of transform trials(compared to short reg trials in TF blocks)", y = "The RT diff. of long trials between TF and RG blocks") +
    annotate("text", x=0, y=120, label=paste0("r = ", corr_coef2), size = 6) +
    annotate("text", x=0, y=95, label=paste0("italic(p) == ", corr_p2), size = 6, parse = TRUE)
)

###########################################################
############# Congruence effect on accRate ################
###########################################################

## visualize the congruence effect

## the stimulus congruence effect
acc_cgc_stim <- regular_trials_error %>%
  convert_as_factor(subject, block_type, incongruence_score_stim) %>%
  group_by(block_type, subject, incongruence_score_stim) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

acc_cgc_stim2 <- acc_cgc_stim %>%
  group_by(block_type, incongruence_score_stim) %>%
  summarise(subaccRate = mean(accRate, na.rm = TRUE),
            subsd = sd(accRate, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%
  ungroup()

(
  p <- ggplot(acc_cgc_stim2, aes(x = incongruence_score_stim, y = subaccRate, color = block_type)) +
    geom_point(size = 6, position=position_dodge(0.2)) +
    geom_line(aes(group = block_type), size = 1.7, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin = subaccRate - se, ymax = subaccRate + se),
                  width = .1, size = 1,
                  position = position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(0.85, 0.93) +
    labs(x = "stimulus incongruency", y = "average accuracy") +
    theme_bw()
)

# the rule congruence effect

acc_cgc_rule <- regular_trials_error %>%
  convert_as_factor(subject, block_type) %>%
  group_by(block_type, subject, incongruence_score_rule) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

acc_cgc_rule2 <- acc_cgc_rule %>%
  group_by(block_type, incongruence_score_rule) %>%
  summarise(subaccRate = mean(accRate, na.rm = TRUE),
            subsd = sd(accRate, na.rm = TRUE),
            se = subsd/((n_sub)^.5)) %>%
  ungroup()

(
  p <- ggplot(acc_cgc_rule2, aes(x = incongruence_score_rule, y = subaccRate, color = block_type)) +
    geom_point(size = 6, position=position_dodge(0.2)) +
    geom_line(aes(group = block_type), size = 1.7, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin = subaccRate - se, ymax = subaccRate + se),
                  width = .1, size = 1,
                  position = position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(0.80, 0.96) +
    labs(x = "task rule incongruency", y = "average accuracy") +
    theme_bw()
)

#################################################################################
############# Congruency effect on IES(inverse efficiency score) ################
#################################################################################

######### compute and visualize the IES of stim congruency effect ###########
ies_cgc_stim <- left_join(rt_cgc_stim, acc_cgc_stim, by = c("subject", "block_type", "incongruence_score_stim")) %>%
  mutate(ies = meanRT / accRate) # compute the IES score

ies_cgc_stim2 <- ies_cgc_stim %>%
  group_by(block_type, incongruence_score_stim) %>%
  summarise(sub_ies = mean(ies, na.rm = TRUE)) %>%
  ungroup()

(
  p <- ggplot(ies_cgc_stim, aes(x = incongruence_score_stim, y = ies, color = block_type)) +
    geom_point(size = 2, position=position_dodge(0.2), alpha = .08) +
    geom_line(aes(group = interaction(subject, block_type)), 
              size = 1, position=position_dodge(0.2), alpha = .2) +
    geom_point(data = ies_cgc_stim2, aes(y = sub_ies), 
               size = 6, position=position_dodge(0.2)) +
    geom_line(data = ies_cgc_stim2, aes(y = sub_ies, group = block_type),
              size = 1.7, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(1000, 1500) +
    theme_bw()
)

### the anova test
ies_cgc_stim.aov <- anova_test(
  data = ies_cgc_stim,
  dv = ies,
  wid = subject,
  within = c(block_type, incongruence_score_stim)
)

get_anova_table(ies_cgc_stim.aov)


######### compute and visualize the IES of rule congruency effect ###########
ies_cgc_rule <- left_join(rt_cgc_rule, acc_cgc_rule, by = c("subject", "block_type", "incongruence_score_rule")) %>%
  mutate(ies = meanRT / accRate) # compute the IES score

ies_cgc_rule2 <- ies_cgc_rule %>%
  group_by(block_type, incongruence_score_rule) %>%
  summarise(sub_ies = mean(ies, na.rm = TRUE)) %>%
  ungroup()

(
  p <- ggplot(ies_cgc_rule, aes(x = incongruence_score_rule, y = ies, color = block_type)) +
    geom_point(size = 2, position=position_dodge(0.2), alpha = .08) +
    geom_line(aes(group = interaction(subject, block_type)), 
              size = 1, position=position_dodge(0.2), alpha = .2) +
    geom_point(data = ies_cgc_rule2, aes(y = sub_ies), 
               size = 6, position=position_dodge(0.2)) +
    geom_line(data = ies_cgc_rule2, aes(y = sub_ies, group = block_type),
              size = 1.7, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(1000, 1500) +
    theme_bw()
)

### the anova test
ies_cgc_rule.aov <- anova_test(
  data = ies_cgc_rule,
  dv = ies,
  wid = subject,
  within = c(block_type, incongruence_score_rule)
)

get_anova_table(ies_cgc_rule.aov)

### the comparison between the rule congruency and the stim congruency effect ###
ies_cgc_stim_comb <- ies_cgc_stim %>%
  rename(incongruence_score = incongruence_score_stim) %>%
  select(!c(meanRT:accRate)) %>%
  mutate(dim = "stim")

ies_cgc_rule_comb <- ies_cgc_rule %>%
  rename(incongruence_score = incongruence_score_rule) %>%
  select(!c(meanRT:accRate)) %>%
  mutate(dim = "rule")

ies_cgc_comb <- bind_rows(ies_cgc_stim_comb, ies_cgc_rule_comb)

ies_cgc_comb.aov <- anova_test(
  data = ies_cgc_comb,
  dv = ies,
  wid = subject,
  within = c(block_type, incongruence_score, dim)
)

get_anova_table(ies_cgc_comb.aov)

################################################################################
############ post-hoc check on the accuracy on the tasks and stimulus ##########
################################################################################

bigdf_clean_error <- bigdf_clean_error %>%
  mutate(image_tar_id = str_extract(image_target, "(?<=id\":\").+(?=\",\"stim_type)"))

length(unique(bigdf_clean_error$image_tar_id)) # should be 24

acc_stim <- bigdf_clean_error %>%
  filter(trial_nature == "RG") %>%
  group_by(image_tar_id) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

glimpse(acc_stim)

ggplot(acc_stim, aes(x = as.factor(image_tar_id), y = accRate)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim=c(0.5,1)) +
  coord_flip(ylim=c(0.5,1))
  

acc_place_y_s_w <-  bigdf_clean_error %>%
  filter(trial_nature == "RG",
         image_tar_id == "place_y_s_w") %>%
  group_by(rule) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()  





