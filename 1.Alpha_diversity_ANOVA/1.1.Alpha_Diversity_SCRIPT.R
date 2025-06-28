## Author: 
#  Panayiotis Louca 

## Purpose of script: 
#  

## Date Created: 
#  19 March 2024 ------------------------------


## Notes: 
#  

# Clear environment 
rm(list = ls()) 

## Set functions: 
#  

## set working directory 
setwd("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/")

## load up packages: 

### core 
library(tidyverse)
library(skimr)
library(broom)

library(readxl)
library(ggpubr)


library(lme4)
library(broom.mixed)
library(lmerTest)

library(betareg)
library(emmeans)

  # ------------------------------------------------------------------------------------------------------------------- # 

##   IMPORT DATA ---- 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/dataset/CC_gm_polyp_dataset.rds")
dat <- read_rds(path) %>% as.data.frame(.)

# check first few names with outcomes  
names(dat)[1:15]

# setup variable denoting all the outcomes 
outcomes_to_test =
  dat %>%
  select(matches("^outcome\\d|^SA\\d")) %>%
  names(.)

outcomes_to_test

sapply(dat[ ,outcomes_to_test], class)

dat %>%
  select(all_of(outcomes_to_test)) %>%
  summary(.)

dat %>% 
  glimpse()  

###   check classes ---- 

dat %>% # all gut microbiome vars are numeric 
  select(
    starts_with("species_"),
    UNCLASSIFIED, 
    alpha_div_shannon, 
    alpha_div_taxonomic_richness) %>%
  sapply(class) %>%
  table() 

dat %>% #
  select(age, 
         sex_imputed, 
         batch,
         BMI_imputed, 
         ethnicity_imputed, 
         site,
         aspirin_imputed,
         ppi,
         years_education_imputed,
         symptom_or_screening_imputed, 
         pack_years_smoking,
         energy_kcal_imputed, 
         alcohol_g_adj_imputed, 
         calcium_mg_adj_imputed,
         nsp_g_adj_imputed, 
         red_meat_g_adj_imputed,
         HDI_imputed) %>%
  skim() 


 
# ------------------------------------------------------------------------------------------------------------------- # 

        # -------------------------------------------------------------------------- # 
        
        # ************************* # 

        ####   DIFFERENCES WITH AND WITHOUT NEOPLASIA  ---- 

        # ************************* # 
        
        # -------------------------------------------------------------------------- #  
  

df_analyse <- dat %>%
  filter(!is.na(outcome1_any_neoplasia_detected))

format(round(mean(df_analyse$alpha_div_taxonomic_richness), 2), nsmall = 2)
format(round(mean(df_analyse$alpha_div_shannon), 2), nsmall = 2)

format(round(sd(df_analyse$alpha_div_taxonomic_richness), 2), nsmall = 2)
format(round(sd(df_analyse$alpha_div_shannon), 2), nsmall = 2)

# no neoplasia detected 
format(round(mean(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 0, 'alpha_div_taxonomic_richness']), 2), nsmall = 2)
format(round(mean(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 0, 'alpha_div_shannon']), 2), nsmall = 2)

format(round(sd(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 0, 'alpha_div_taxonomic_richness']), 2), nsmall = 2)
format(round(sd(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 0, 'alpha_div_shannon']), 2), nsmall = 2)

# neoplasia detected 
format(round(mean(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 1, 'alpha_div_shannon']), 2), nsmall = 2)
format(round(sd(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 1, 'alpha_div_shannon']), 2), nsmall = 2)

format(round(mean(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 1, 'alpha_div_taxonomic_richness']), 2), nsmall = 2)
format(round(sd(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 1, 'alpha_div_taxonomic_richness']), 2), nsmall = 2)


# create object of div measures to test 
alpha_div_metrics <- c("alpha_div_shannon", "alpha_div_taxonomic_richness")

  tmp_results <- data.frame()  # Create an empty dataframe to store the results 
  
  # Loop through each metric 
  for (alpha_div_metric in alpha_div_metrics) {
    
  metric <- alpha_div_metric

      # run unadjusted model 
  
  # Dynamically create formula 
  formula <- formula(paste(
    alpha_div_metric,
    paste("~ age + sex_imputed + BMI_imputed + ethnicity_imputed + site + batch + outcome1_any_neoplasia_detected", 
          sep = "")
  ))
  

  
  shannon_model <- aov(data = df_analyse,
                       formula)
  
  
  # extract p value 
  pval = tidy(shannon_model) %>%
    filter(term == 'outcome1_any_neoplasia_detected') %>%
    pull(p.value)

    
    # calculate means for each group 
    level_0_mean <- mean(unlist(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 0, alpha_div_metric]))
    level_1_mean <- mean(unlist(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 1, alpha_div_metric]))
  
  
  difference <- ifelse(!is.na(level_0_mean) & level_1_mean > level_0_mean, yes = "higher",
                       no =
                         ifelse(!is.na(level_0_mean) & level_1_mean < level_0_mean, yes = "lower",
                                no = NA))
  
  
  tmp <- data.frame(outcome = "outcome1_any_neoplasia_detected", metric = metric,
                    model = 'base_model (age, sex, BMI, ethnicity, site and batch)',
                    pval = pval, level_0_mean = level_0_mean, level_1_mean = level_1_mean, 
                    difference = difference,
                    model = paste(deparse(formula), collapse = ""))
  
  tmp_results <- rbind(tmp_results, tmp)
  
      # -------------------------------------------------------------------------- #    
  
  # run adjusted model 
  
  
  # Dynamically create formula 
  formula <- formula(paste(
    alpha_div_metric,
    paste("~ age + sex_imputed + BMI_imputed + ethnicity_imputed + site + batch + energy_kcal_imputed + HDI_imputed + alcohol_g_adj_imputed + pack_years_smoking +  ppi + aspirin_imputed + outcome1_any_neoplasia_detected", 
          sep = "")
  ))
  
  
  shannon_model <- aov(data = df_analyse,
                       formula)
  
  # extract p value 
  pval = tidy(shannon_model) %>%
    filter(term == 'outcome1_any_neoplasia_detected') %>%
    pull(p.value)
  
  
  # calculate means for each group 
  level_0_mean <- mean(unlist(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 0, alpha_div_metric]))
  level_1_mean <- mean(unlist(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 1, alpha_div_metric]))
  
  
  difference <- ifelse(!is.na(level_0_mean) & level_1_mean > level_0_mean, yes = "higher",
                       no =
                         ifelse(!is.na(level_0_mean) & level_1_mean < level_0_mean, yes = "lower",
                                no = NA))
  
  
  tmp <- data.frame(outcome = "outcome1_any_neoplasia_detected", metric = metric,
                    model = 'adjusted',
                    pval = pval, level_0_mean = level_0_mean, level_1_mean = level_1_mean, 
                    difference = difference,
                    model = paste(deparse(formula), collapse = ""))
  
  tmp_results <- rbind(tmp_results, tmp)
  
  
      # -------------------------------------------------------------------------- #  
  
  # run LMER adjusted model using random effects for site and batch 
  
  
  # Dynamically create formula 
  formula <- formula(paste(
    alpha_div_metric,
    paste("~ age + sex_imputed + BMI_imputed + ethnicity_imputed + (1 | site) + (1 | batch) + energy_kcal_imputed + HDI_imputed + alcohol_g_adj_imputed + pack_years_smoking +  ppi + aspirin_imputed + outcome1_any_neoplasia_detected",
          sep = "")
  ))
  
  
  lmer_model <- lmer(data = df_analyse,
                       formula)
  
  
  # extract p value 
  pval = tidy(lmer_model) %>%
    filter(term == 'outcome1_any_neoplasia_detected1') %>%
    pull(p.value)
  
  
  # calculate means for each group 
  level_0_mean <- mean(unlist(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 0, alpha_div_metric]))
  level_1_mean <- mean(unlist(df_analyse[df_analyse$outcome1_any_neoplasia_detected == 1, alpha_div_metric]))
  
  
  difference <- ifelse(!is.na(level_0_mean) & level_1_mean > level_0_mean, yes = "higher",
                       no =
                         ifelse(!is.na(level_0_mean) & level_1_mean < level_0_mean, yes = "lower",
                                no = NA))
  
  
  tmp <- data.frame(outcome = "outcome1_any_neoplasia_detected", metric = metric,
                    model = 'LMER_adjusted',
                    pval = pval, level_0_mean = level_0_mean, level_1_mean = level_1_mean, 
                    difference = difference,
                    model = paste(deparse(formula), collapse = ""))
  
  tmp_results <- rbind(tmp_results, tmp)
  
}

  with_and_without_neoplasia_results <- tmp_results
  
  
  write.csv(with_and_without_neoplasia_results,
            "/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_alpha_diversity/with_and_without_neoplasia_alpha_div_RESULTS.csv",
            row.names = FALSE)
  
      # -------------------------------------------------------------------------- #  
  
  # ************************* # 
  
  ####   DIFFERENCES BY NEOPLASIA TYPE  ---- 
  
  # ************************* # 
  
  # -------------------------------------------------------------------------- #  
  
  # define polyp type outcomes 
  outcomes = c("outcome1_cancer","outcome2_advanced_neoplasia_exc_cancer","outcome3_non_advanced_neoplasia",
               "outcome4_only_adenomatous","outcome5_only_serrated_lesions","outcome6_only_sessile_serrated_lesions",
               'SA1_advanced_neoplasia_inc_cancer')
  
  # Create an empty dataframe to store the results 
  By_neoplasia_type_result_list <- list()
  
  # Loop through each factor and outcome 
  for (outcome in outcomes) {
    
  df_analyse <- dat %>%
    filter(!is.na(!!sym(outcome))) 
  
  # create object of div measures to test 
  alpha_div_metrics <- c("alpha_div_shannon", "alpha_div_taxonomic_richness")
  
  tmp_results <- data.frame()  # Create an empty dataframe to store the results 
  
  # Loop through each metric 
  for (alpha_div_metric in alpha_div_metrics) {
    
    metric <- alpha_div_metric
    
    # run unadjusted model 
    
    # Dynamically create formula 
    formula <- formula(paste(
      alpha_div_metric,
      paste("~ age + sex_imputed + BMI_imputed + ethnicity_imputed + site + batch + ", 
            outcome,
            sep = "")
    ))
    
    
    shannon_model <- aov(data = df_analyse,
                         formula)
    
    
    # extract p value 
    pval = tidy(shannon_model) %>%
      filter(term == outcome) %>%
      pull(p.value)
    
    
    # calculate means for each group 
    level_0_mean <- mean(unlist(df_analyse[df_analyse[ ,outcome] == 0, alpha_div_metric]))
    level_1_mean <- mean(unlist(df_analyse[df_analyse[ ,outcome] == 1, alpha_div_metric]))
    
    
    difference <- ifelse(!is.na(level_0_mean) & level_1_mean > level_0_mean, yes = "higher",
                         no =
                           ifelse(!is.na(level_0_mean) & level_1_mean < level_0_mean, yes = "lower",
                                  no = NA))
    
    
    tmp <- data.frame(outcome = outcome, metric = metric,
                      model = 'unadjusted',
                      pval = pval, level_0_mean = level_0_mean, level_1_mean = level_1_mean, 
                      difference = difference,
                      model_formula = paste(deparse(formula), collapse = ""))
    
    tmp_results <- rbind(tmp_results, tmp)
    
    # -------------------------------------------------------------------------- #    
    
    # run adjusted model 
    
    
    # Dynamically create formula 
    formula <- formula(paste(
      alpha_div_metric,
      paste("~ age + sex_imputed + BMI_imputed + ethnicity_imputed + site + batch + energy_kcal_imputed + HDI_imputed + alcohol_g_adj_imputed + pack_years_smoking +  ppi + aspirin_imputed + ", 
            outcome,
            sep = "")
    ))
    
    
    shannon_model <- aov(data = df_analyse,
                         formula)
    
    # extract p value 
    pval = tidy(shannon_model) %>%
      filter(term == outcome) %>%
      pull(p.value)
    
    
    # calculate means for each group 
    level_0_mean <- mean(unlist(df_analyse[df_analyse[ ,outcome] == 0, alpha_div_metric]))
    level_1_mean <- mean(unlist(df_analyse[df_analyse[ ,outcome] == 1, alpha_div_metric]))
    
    
    difference <- ifelse(!is.na(level_0_mean) & level_1_mean > level_0_mean, yes = "higher",
                         no =
                           ifelse(!is.na(level_0_mean) & level_1_mean < level_0_mean, yes = "lower",
                                  no = NA))
    
    
    tmp <- data.frame(outcome = outcome, metric = metric,
                      model = 'adjusted',
                      pval = pval, level_0_mean = level_0_mean, level_1_mean = level_1_mean, 
                      difference = difference,
                      model_formula = paste(deparse(formula), collapse = ""))
    
    tmp_results <- rbind(tmp_results, tmp)
    
    
    # -------------------------------------------------------------------------- #  
    
    # run LMER adjusted model using random effects for site and batch 
    
    
    # Dynamically create formula 
    formula <- formula(paste(
      alpha_div_metric,
      paste("~ ",
            outcome, 
            " + age + sex_imputed + BMI_imputed + ethnicity_imputed + (1 | site) + (1 | batch) + energy_kcal_imputed + HDI_imputed + alcohol_g_adj_imputed + pack_years_smoking +  ppi + aspirin_imputed", 
            sep = "")
    ))
    
    
    
    lmer_model <- lmer(data = df_analyse,
                       formula)
    
    
    
    # extract p value 
    pval = tidy(lmer_model) %>%
      filter(term == paste0(outcome, "1")) %>%
      pull(p.value)

    
    # calculate means for each group 
    level_0_mean <- mean(unlist(df_analyse[df_analyse[ ,outcome] == 0, alpha_div_metric]))
    level_1_mean <- mean(unlist(df_analyse[df_analyse[ ,outcome] == 1, alpha_div_metric]))
    
    
    difference <- ifelse(!is.na(level_0_mean) & level_1_mean > level_0_mean, yes = "higher",
                         no =
                           ifelse(!is.na(level_0_mean) & level_1_mean < level_0_mean, yes = "lower",
                                  no = NA))
    
    
    tmp <- data.frame(outcome = outcome, metric = metric,
                      model = 'LMER_adjusted',
                      pval = pval, level_0_mean = level_0_mean, level_1_mean = level_1_mean, 
                      difference = difference,
                      model_formula = paste(deparse(formula), collapse = ""))
    
    tmp_results <- rbind(tmp_results, tmp)
    
  }
  
  By_neoplasia_type_result_list[[outcome]] <- tmp_results
  }
  
    # -------------------------------------------------------------------------- #  

########   Check results  

names(By_neoplasia_type_result_list)
any(By_neoplasia_type_result_list$outcome1_cancer$pval < 0.05)
any(By_neoplasia_type_result_list$outcome2_advanced_neoplasia_exc_cancer$pval < 0.05)
any(By_neoplasia_type_result_list$outcome3_non_advanced_neoplasia$pval < 0.05)
any(By_neoplasia_type_result_list$outcome4_only_adenomatous$pval < 0.05)
any(By_neoplasia_type_result_list$outcome5_only_serrated_lesions$pval < 0.05)
any(By_neoplasia_type_result_list$outcome6_only_sessile_serrated_lesions$pval < 0.05)
any(By_neoplasia_type_result_list$outcome7_mixed_polyp_profile$pval < 0.05)
any(By_neoplasia_type_result_list$outcome8_complex_polyp_profile$pval < 0.05)
any(By_neoplasia_type_result_list$outcome9_polyp_burden$pval < 0.05)
any(By_neoplasia_type_result_list$SA1_advanced_neoplasia_inc_cancer$pval < 0.05)

########   Collapse results  
By_neoplasia_type_results <- do.call(rbind, By_neoplasia_type_result_list)

########   Save results  
write.csv(By_neoplasia_type_results,
          '/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_alpha_diversity/by_neoplasia_type_alpha_div_RESULTS.csv',
          row.names = FALSE) 

    