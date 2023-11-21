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

length(unique(bigdf_clean$subject))

bigdf_clean <- bigdf_clean %>%
  mutate(rt = rt*1000)

### step1. only include the regular trials for analysis ###
regular_trials <- bigdf_clean %>%
  filter(trial_nature == "RG") %>%
  select(CTI_length, rt, block_type, block_id, stim, rule, CTI, resp_map_age,
         resp_map_size, resp_map_location, subject, accuracy, incongruence_score_stim, incongruence_score_rule)

### step1. only include the transform trials for analysis ###
transform_trials <- bigdf_clean %>%
  filter(trial_nature == "TF") %>%
  select(CTI_length, rt, block_type, block_id, stim, rule, CTI, resp_map_age,
         resp_map_size, resp_map_location, subject, accuracy, incongruence_score_stim, incongruence_score_rule)

###########################################################
################## the Analysis on the RT #################
###########################################################

# - RT summary and plotting of regular trials ========
rt_descrip <- regular_trials %>%
  convert_as_factor(CTI_length, subject, block_type) %>%
  group_by(subject, block_type, CTI_length) %>%
  summarise(count = n(),
            meanCTI = mean(CTI),
            meanRT = mean(rt, na.rm = TRUE)) %>%
  ungroup()

n_sub <- length(unique(rt_descrip$subject))

# -- subplots for each subject =======================
(
  p <- ggplot(rt_descrip, aes(x = CTI_length, y = meanRT, color = block_type)) +
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
  p <- ggplot(rt_descrip2, aes(x = CTI_length, y = subRT, color = block_type)) +
    geom_point(size = 5, position = position_dodge(0.1)) +
    geom_line(aes(group = block_type), size = 1.5,
              position = position_dodge(0.1),
              linetype="dotted") +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .1, size = 1,
                  position = position_dodge(0.1)) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    labs(y = "RT in ms") +
    theme_bw()
)

ggsave("RT_interaction.png", width = 4.8, height = 3.5)

# -- plot both individual and group level ============
(
  p <- ggplot(rt_descrip, aes(x = CTI_length, y = meanRT, color = block_type)) +
    geom_point(size = 4, alpha = .1) +
    geom_line(aes(group = interaction(subject, block_type)),
              size = 2, alpha = .1) +
    geom_point(data = rt_descrip2, aes(y = subRT), size = 5) +
    geom_line(data = rt_descrip2, aes(y = subRT, group = block_type), size = 1) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    theme_bw()
)

# -- plotting for publication ========================
# plot the CTI bins as categorical variables with equal distance between bins
(
  p <- ggplot(rt_descrip2, aes(x = CTI_bins, y = subRT, group = block_type, color = block_type)) +
    geom_point(size = 5, position=position_dodge(0.2)) +
    geom_line(size = 1, position=position_dodge(0.2)) +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .2,
                  position = position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    ylim(1050, 1200) +
    theme_apa(base_size = 14) +
    labs(x = "CTI bin", y = "Reaction time in ms")
)

# plot the CTI bins as continuous variables -- using the mean CTI as distance between bins
cti_mean <- regular_trials %>%
  filter(session == "online")
  convert_as_factor(CTI_bins) %>%
  group_by(CTI_bins) %>%
  summarise(meanCTI = mean(CTI)) %>%
  ungroup()

(
  p <- ggplot(rt_descrip2, aes(x = submeanCTI, y = subRT, group = block_type, color = block_type)) +
    geom_point(size = 5) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = subRT - se, ymax = subRT + se),
                  width = .2) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07"),
                       name = "block type:",
                       breaks = c("RG", "TF"),
                       labels = c("Regular", "Transform")) +
    scale_x_continuous(breaks = cti_mean$meanCTI,
                       labels = seq(1,6)) +    
    ylim(1050, 1200) +
    theme_apa(base_size = 14) +
    labs(x = "CTI bin", y = "Reaction time in ms")
)

# add error bar
(
  p <- ggplot(rt_descrip, aes(x = CTI_bins, y = meanRT, color = block_type)) +
    geom_point(size = 4, alpha = .06) +
    geom_point(data = rt_descrip2, aes(y = subRT), size = 5,
               position=position_dodge(0.2)) +
    geom_pointrange(data = rt_descrip2, aes(y = subRT, ymin = subRT - subsd, ymax = subRT + subsd), size =.8,
                    position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(900, 1350) +
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
            subsd = sd(meanRT, na.rm = TRUE)) %>%
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
    geom_vline(aes(xintercept=5000),linetype="dashed", size=1) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    theme_bw()
)

# -- ANOVA testing =======================
rt_online.aov <- anova_test(
  data = rt_descrip,
  dv = meanRT,
  wid = subject,
  within = c(block_type, CTI_length)
)

get_anova_table(rt_online.aov)

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

# specify the model using continuous CTI #
model_RT <- lmer(formula = rt ~ (CTI_sec + incongruence_score_stim + incongruence_score_rule) * block_type
                 + (1 | subject),
                 data = regular_trials_lme)

anova(model_RT)
summary(model_RT, correlation= FALSE)

# specify the model with using CTI length(short vs. long) #
model_RT2 <- lmer(formula = rt ~ (CTI_length + incongruence_score_stim + incongruence_score_rule) * block_type
                        + (1 | subject),
                        data = regular_trials_lme)

anova(model_RT2)
summary(model_RT2, correlation= FALSE)

regular_trials_lme$fitted_RT_Bi <- fitted(model_RT2)

# the maximum model:
model_RT_max <- lmer(formula = rt ~ (CTI_length + incongruence_score_stim + incongruence_score_rule) * block_type
                  + ((CTI_length + incongruence_score_stim + incongruence_score_rule) * block_type || subject),
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

# - analyzing on transform trials =======================
rt_tran_descrip1 <- transform_trials %>%
  convert_as_factor(subject) %>%
  group_by(subject, CTI) %>%
  summarise(count = n(),
            meanRT = mean(rt, na.rm = TRUE)) %>%
  ungroup()

rt_tran_descrip2 <- rt_tran_descrip1 %>%
  group_by(CTI) %>%
  summarise(subRT = mean(meanRT, na.rm = TRUE),
            subsd = sd(meanRT, na.rm = TRUE)) %>%
  ungroup()

(
  p <- ggplot(rt_tran_descrip2, aes(x = CTI, y = subRT)) +
    geom_point(size = 5,position=position_dodge(0.2)) +
    geom_line(size = 1)  +
    scale_color_manual(values = "#4E84C4") +
    ylim(1000, 1250) +
    theme_bw()
)

#### the across-subject correlation between block-type diff in the regular trials and
#### the performance of transform trials

rt_regular_6bin <- regular_trials %>%
  filter(CTI_bins != 6) %>%
  group_by(subject, block_type) %>%
  summarise(meanRT = mean(rt, na.rm = TRUE)) %>%
  pivot_wider(names_from = block_type, values_from = meanRT) %>%
  mutate(sixbin_effect = TF - RG) %>%
  ungroup() %>%
  select(subject, sixbin_effect)

subject_regularblock_mean <- regular_trials %>%
  filter(block_type == "TF", CTI_bins == 6) %>%
  group_by(subject) %>%
  summarise(meanRT = mean(rt, na.rm = TRUE)) %>%
  ungroup()

subject_tranTrial_mean <- transform_trials %>%
  group_by(subject) %>%
  summarise(meanRT_tran = mean(rt, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(subject_regularblock_mean, by = "subject") %>%
  mutate(tranRT = meanRT_tran - meanRT) %>%
  left_join(rt_regular_6bin, by = "subject")
  
cor.test(subject_tranTrial_mean$tranRT, subject_tranTrial_mean$sixbin_effect)
cor.test(subject_tranTrial_mean$meanRT_tran, subject_tranTrial_mean$sixbin_effect)


# pairwise t tests
rt_descrip %>%
  pairwise_t_test(meanRT ~ CTI_bins,
                  paired = TRUE,
                  p.adjust.method = "bonferroni")

rt_descrip %>%
  pairwise_t_test(meanRT ~ block_type,
                  paired = TRUE,
                  p.adjust.method = "bonferroni")


###########################################################
########### the Analysis on the Accuracy rate #############
###########################################################

load(file = paste0(here(), "/results/bigdf_fmri_clean_error.Rdata"))
bigdf_clean_error <- bigdf_clean_error %>%
  filter(!is.na(rt)) # get rid of missing trials before analyzing accuracy

### step1. only include the regular trials for analysis ###
regular_trials_error <- bigdf_clean_error %>%
  filter(trial_nature == "RG") %>%
  select(rt, block_type, block_id, stim, rule, CTI, resp_map_age,
         resp_map_size, resp_map_location, subject, accuracy, CTI_length,
         incongruence_score_stim, incongruence_score_rule)

tran_trials_error <- bigdf_clean_error %>%
  filter(trial_nature == "TF") %>%
  select(rt, block_type, block_id, stim, rule, CTI, resp_map_age,
         resp_map_size, resp_map_location, subject, accuracy, CTI_length,
         incongruence_score_stim, incongruence_score_rule) %>%
  mutate(ACC_bi = if_else(accuracy == TRUE, 1, 0))

### step2. compare the subject level accuracy rate between regular and transform block

accRate_descrip <- regular_trials_error %>%
  group_by(block_type, subject, CTI_length) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

accRate_descrip2 <- accRate_descrip %>%
  group_by(block_type, CTI_length) %>%
  summarise(subaccRate = mean(accRate, na.rm = TRUE)) %>%
  ungroup()

(
  p <- ggplot(accRate_descrip, aes(x = CTI_length, y = accRate, color = block_type)) +
    geom_point(size = 4, alpha = .05) +
    geom_point(data = accRate_descrip2, aes(y = subaccRate), 
               size = 5, position=position_dodge(0.2)) +
    geom_line(data = accRate_descrip2, aes(y = subaccRate, group = block_type), 
              size = 1, position=position_dodge(0.2)) +
    scale_x_discrete(name = "CTI length", limits = c("short", "long")) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(0.87, 0.95) +
    theme_bw()
)

### step 4. visualize the data w/o the bins ###
accRate_descrip4 <- regular_trials_error %>%
  group_by(block_type, CTI) %>%
  summarise(grandAccRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

(
  p <- ggplot(accRate_descrip4, aes(x = CTI, y = grandAccRate, color = block_type)) +
    geom_point(size = 5, position=position_dodge(0.2)) +
    geom_line(aes(group = block_type), 
              size = 1, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(0.75, 0.99) +
    theme_bw()
)


### the mixed effect logistic regression model on ACC ###

options(contrasts = c("contr.sum","contr.poly"))
options(scipen = FALSE)

regular_trials_error_lme <- regular_trials_error %>%
  mutate(CTI_sec = CTI/1000,
         ACC_bi = if_else(accuracy == TRUE, 1, 0)) %>%
  convert_as_factor(block_type, subject) %>%
  ungroup()

contrasts(regular_trials_error_lme$block_type)

model_acc <- glmer(formula = ACC_bi ~ (CTI_sec + incongruence_score_stim + incongruence_score_rule) * block_type 
                   + (1|subject),
                   data = regular_trials_error_lme,
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa"))

anova(model_acc)
summary(model_acc, correlation= FALSE)

# for reporting #
statsACC.table <- as.data.frame(summary(model_acc)$coefficients)
names(statsACC.table) <- c("estimate", "SE", "z value", "p")

statsACC.table2 <- as.data.frame(matrix(ncol = 4, nrow = 8))
colnames(statsACC.table2) <- c('Coefficient', 'Estimate', 'SE', 'p value')

statsACC.table2$Coefficient <- c("(Intercept)",
                                  "CTI",
                                  "incongruency score of stimulus type",
                                  "incongruency score of task rule",
                                  "block type",
                                  "CTI:block type",
                                  "incongruency score of stimulus type:block type",
                                  "incongruency score of task rule:block type")
statsACC.table2$Estimate <- round(statsACC.table$estimate, digit= 3)
statsACC.table2$SE <- round(statsACC.table$SE, digit= 3)
statsACC.table2$"p value" <- round(statsACC.table$p, digit= 4)

fun_round <- function(x) {formatC(x, format = "f", digits = 3)}

tableACC_ready <- nice_table(statsACC.table2, 
                             col.format.p = 4, 
                             col.format.custom = 2:3,
                             format.custom = "fun_round",
                             title = "Table \nThe result of mixed effect regression model on accuracy in rule transform experiment", 
                             footnote = "SE = standard error.\n* denotes p < .05,
                          ** denotes p < .01, *** denates p < .001")
tableACC_ready
save_as_docx(table_ready, path = paste0(getwd(), "/nice_tablehere.docx"))

# merging RT and ACC data #
colnames(stats.table2) <- c('Coefficient', 'RT.Estimate', 'RT.SE', 'RT.p value')
colnames(statsACC.table2) <- c('Coefficient', 'ACC.Estimate', 'ACC.SE', 'ACC.p value')

statsALL.table <- merge(x=stats.table2, y=statsACC.table2,
                        by='Coefficient', all.x = TRUE)

tableALL_ready <- nice_table(statsALL.table, 
                             col.format.p = c(4,7), 
                             col.format.custom = c(2:3,5:6),
                             format.custom = "fun_round",
                             separate.header = TRUE,
                             italics = seq(statsALL.table),
                             title = "Table 1 \nThe result of mixed effect regression model on RT and accuray in rule transform experiment", 
                             footnote = "RT = reaction time, ACC = accuracy, SE = standard error.\n* denotes p < .05, ** denotes p < .01, *** denotes p < .001")

tableALL_ready
save_as_docx(tableALL_ready, path = paste0(here(), "/APA_table_rule_tran.docx"))

# visualize the 2-way interaction between CTI and block type
emmip(model_acc, block_type ~ CTI_sec, cov.reduce = range, type = "response")
emtrends(model_acc, ~ block_type, var = 'CTI_sec', type = "response")
test(emtrends(model_acc, ~ block_type, var = 'CTI_sec'))
# 3-way interaction #
model_acc2 <- glmer(formula = ACC_bi ~ (incongruence_score_stim + incongruence_score_rule) * CTI_sec * block_type
                    + (1 | subject),
                   data = regular_trials_error_lme,
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa"))
anova(model_acc2)
summary(model_acc2, correlation= FALSE)

regular_trials_error_lme %>%
  convert_as_factor(incongruence_score_rule) %>%
  mutate(CTI_bins = as.numeric(CTI_bins)) %>%
  group_by(incongruence_score_rule, CTI_bins) %>%
  summarise(meanAcc = sum(accuracy)/length(accuracy)) %>%
  ggplot(aes(x = incongruence_score_rule, y = meanAcc, color = CTI_bins)) +
  geom_point(size = 4) +
  geom_line(aes(group = CTI_bins), size = 1) +
  theme_bw()


### logistic regression on acc of transform trials ###

tran_trials_error_lme <- tran_trials_error %>%
  mutate(CTI_sec = CTI/1000)

model_tranAcc <- glmer(formula = ACC_bi ~ CTI_sec + (1|subject),
                       data = tran_trials_error_lme,
                       family = binomial,
                       control = glmerControl(optimizer = "bobyqa"))

summary(model_tranAcc)

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
  summarise(subaccRate = mean(accRate, na.rm = TRUE)) %>%
  ungroup()

(
  p <- ggplot(acc_cgc_stim, aes(x = incongruence_score_stim, y = accRate, color = block_type)) +
    geom_point(size = 2, position=position_dodge(0.2), alpha = .08) +
    geom_line(aes(group = interaction(subject, block_type)), 
              size = 1, position=position_dodge(0.2), alpha = .2) +
    geom_point(data = acc_cgc_stim2, aes(y = subaccRate), 
               size = 6, position=position_dodge(0.2)) +
    geom_line(data = acc_cgc_stim2, aes(y = subaccRate, group = block_type),
              size = 1.7, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(0.7, 1.0) +
    theme_bw()
)


acc_cgc_stim.aov <- anova_test(
  data = acc_cgc_stim,
  dv = accRate,
  wid = subject,
  within = c(block_type, incongruence_score_stim)
)

get_anova_table(acc_cgc_stim.aov)


# the rule congruence effect

acc_cgc_rule <- regular_trials_error %>%
  convert_as_factor(subject, block_type) %>%
  group_by(block_type, subject, incongruence_score_rule) %>%
  summarise(count = n(),
            accRate = sum(accuracy)/length(accuracy)) %>%
  ungroup()

acc_cgc_rule2 <- acc_cgc_rule %>%
  group_by(block_type, incongruence_score_rule) %>%
  summarise(subaccRate = mean(accRate, na.rm = TRUE)) %>%
  ungroup()

(
  p <- ggplot(acc_cgc_rule, aes(x = incongruence_score_rule, y = accRate, color = block_type)) +
    geom_point(size = 2, position=position_dodge(0.2), alpha = .08) +
    geom_line(aes(group = interaction(subject, block_type)), 
              size = 1, position=position_dodge(0.2), alpha = .2) +
    geom_point(data = acc_cgc_rule2, aes(y = subaccRate), 
               size = 6, position=position_dodge(0.2)) +
    geom_line(data = acc_cgc_rule2, aes(y = subaccRate, group = block_type),
              size = 1.7, position=position_dodge(0.2)) +
    scale_color_manual(values = c("#4E84C4", "#FC4E07")) +
    ylim(0.85, 0.95) +
    theme_bw()
)

acc_cgc_rule.aov <- anova_test(
  data = acc_cgc_rule,
  dv = accRate,
  wid = subject,
  within = c(block_type, incongruence_score_rule)
)

get_anova_table(acc_cgc_rule.aov)


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





