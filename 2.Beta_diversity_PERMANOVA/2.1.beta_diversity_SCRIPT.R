## Author: 
#  Panayiotis Louca 

## Purpose of script: 
#  

## Date Created: 
#  26 March 2024 ------------------------------


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
library(vegan)

#   IMPORT DATA ---- 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/dataset/CC_gm_polyp_dataset.rds")
dat <- read_rds(path) %>% as.data.frame(.)

# setup variable denoting all the outcomes 
outcomes_to_test =
  dat %>%
  select(matches("^outcome\\d|^SA\\d")) %>%
  names(.)

outcomes_to_test

dat %>% 
  glimpse()  

dat %>%
  select(all_of(outcomes_to_test)) %>%
  summary(.)

###   check classes ---- 

dat %>% # all gut microbiome vars are numeric 
  select(
    starts_with("species_"),
    UNCLASSIFIED, 
    alpha_div_shannon, 
    alpha_div_simpson,
    alpha_div_taxonomic_richness) %>%
  sapply(class) %>%
  table() 

dat %>% 
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

# setup unadjusted and adjusted formula (without outcomes) 
base_formula_unadjusted <- paste(
  "df_gm ~ age + sex_imputed + BMI_imputed + ethnicity_imputed + batch", 
  sep = "")

base_formula_adjusted <- paste(
  "df_gm ~ age + sex_imputed + BMI_imputed + ethnicity_imputed + energy_kcal_imputed + HDI_imputed + alcohol_g_adj_imputed + pack_years_smoking + ppi + aspirin_imputed + batch",
  sep = "")

all_vars <- c('age', 
              'sex_imputed',
              'BMI_imputed', 
              'ethnicity_imputed', 
              'energy_kcal_imputed', 
              'HDI_imputed', 
              'alcohol_g_adj_imputed',
              'pack_years_smoking', 
              'ppi',
              'aspirin_imputed',
              'site', 
              'batch')
    
# -------------------------------------------------------------------------- # 

# ************************* # 

####   DIFFERENCES WITH AND WITHOUT NEOPLASIA  ---- 

# ************************* # 

# -------------------------------------------------------------------------- #  
 
set.seed(1234)

outcome_results <- data.frame()

# split data 

df_gm <-  # just gut microbiome data and iid coerced to rownmaes 
  dat %>%
  filter(!is.na(!!sym('outcome1_any_neoplasia_detected')) ) %>%
  select(record_id,
         starts_with("species_") ) %>%
  column_to_rownames('record_id') %>%
  as.matrix()


df_meta <- # just metadata 
  dat %>% 
  filter(!is.na(!!sym('outcome1_any_neoplasia_detected')) ) %>%
  select(record_id,
         outcome1_any_neoplasia_detected,
         all_of(all_vars)) %>%
  column_to_rownames('record_id')

# ------------------------------------------------------------------------------------------------------------------- # 

#####   Permanova ---- 

######   Bray-Curtis  ---- 

res <- data.frame()
#########   unadjusted model ---- 


formula_unadjusted <- as.formula(paste(base_formula_unadjusted, " + outcome1_any_neoplasia_detected", sep = ""))

 perm_res <- adonis2(
     formula = formula_unadjusted,
       data = df_meta,
       permutations = 1000,
       method = "bray",
       by = "margin",  # alt: "terms"
       strata = df_meta$site,
       parallel = parallel::detectCores() - 1
   )
 
 perm_res
 
 # R2 
 r_sq <- perm_res[grep("outcome1_any_neoplasia_detected", rownames(perm_res)), "R2"] 
 # P-value 
 pval <- perm_res[grep("outcome1_any_neoplasia_detected", rownames(perm_res)), 'Pr(>F)'] 
 
 unadj <- data.frame(outcome = 'outcome1_any_neoplasia_detected', metric = "Bray-Curtis", r2_basic_adj = r_sq, pvalue_basic_adj = pval) 
 unadj


# ------------------------------------------------------------------------------------------------------------------- # 

######   full model ---- 

formula_adjusted <- as.formula(paste(base_formula_adjusted, " + outcome1_any_neoplasia_detected", sep = ""))

perm_res <- adonis2(
  formula = formula_adjusted,
  data = df_meta,
  permutations = 1000,
  method = "bray",
  by = "margin",  # or "terms"
  strata = df_meta$site,
  parallel = parallel::detectCores() - 1
)

perm_res
# R2 
r_sq <- perm_res[grep("outcome1_any_neoplasia_detected", rownames(perm_res)), "R2"] 
# P-value 
pval <- perm_res[grep("outcome1_any_neoplasia_detected", rownames(perm_res)), 'Pr(>F)'] 

adj <- data.frame(r2_full_adj = r_sq, pvalue_full_adj = pval) 
adj

res <- cbind(unadj,adj)
res

outcome_res <- res

######   Jaccard  ---- 

res <- data.frame()

#########   unadjusted model ---- 

perm_res <- adonis2(
  formula = formula_unadjusted,
  data = df_meta,
  permutations = 1000,
  method = "jaccard",
  by = "margin",  # or "terms"
  strata = df_meta$site,
  parallel = parallel::detectCores() - 1
)

perm_res
# R2 
r_sq <- perm_res[grep("outcome1_any_neoplasia_detected", rownames(perm_res)), "R2"] 
# P-value 
pval <- perm_res[grep("outcome1_any_neoplasia_detected", rownames(perm_res)), 'Pr(>F)'] 

unadj <- data.frame(outcome = 'outcome1_any_neoplasia_detected', metric = "Jaccard", r2_basic_adj = r_sq, pvalue_basic_adj = pval) 
unadj

# ------------------------------------------------------------------------------------------------------------------- # 

######   full model ---- 

perm_res <- adonis2(
  formula = formula_adjusted,
  data = df_meta,
  permutations = 1000,
  method = "jaccard",
  by = "margin",  # or "terms"
  strata = df_meta$site,
  parallel = parallel::detectCores() - 1
)

perm_res
# R2 
r_sq <- perm_res[grep("outcome1_any_neoplasia_detected", rownames(perm_res)), "R2"] 
# P-value 
pval <- perm_res[grep("outcome1_any_neoplasia_detected", rownames(perm_res)), 'Pr(>F)'] 


adj <- data.frame(r2_full_adj = r_sq, pvalue_full_adj = pval) 
adj



res <- cbind(unadj,adj)
res

# rbind bray and jaccard results 

outcome_res <-  rbind(outcome_res, res)

# ------------------------------------------------------------------------------------------------------------------- # 

# view results 
outcome_res

# save results 
path <- file.path(
  "/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_beta_diversity/with_and_without_neoplasia_beta_div_RESULTS.csv"
  )

write.csv(outcome_res, 
          path,
          row.names = FALSE)

# -------------------------------------------------------------------------- # 

# ************************* # 

####   DIFFERENCES BY POLYP TYPE  ---- 

# ************************* # 

# -------------------------------------------------------------------------- #  

set.seed(1234)

# define polyp type outcomes 
outcomes = c("outcome1_cancer","outcome2_advanced_neoplasia_exc_cancer","outcome3_non_advanced_neoplasia",
             "outcome4_only_adenomatous","outcome5_only_serrated_lesions","outcome6_only_sessile_serrated_lesions",
             "outcome7_mixed_polyp_profile", 'SA1_advanced_neoplasia_inc_cancer')

outcome_results <- data.frame()

for(outcome in outcomes) {
  
  # split data 
  
  df_gm <-  # just gut microbiome data and iid coerced to rownmaes 
    dat %>%
    filter(!is.na(!!sym(outcome)) ) %>%
    select(record_id,
           starts_with("species_") ) %>%
    column_to_rownames('record_id') %>%
    as.matrix()
  
  
  df_meta <- # just metadata 
    dat %>% 
    filter(!is.na(!!sym(outcome)) ) %>%
    select(record_id,
           outcome,
           all_vars) %>%
    column_to_rownames('record_id')
  
  # ------------------------------------------------------------------------------------------------------------------- # 
  
  #####   Permanova ---- 
  
  ######   Bray-Curtis  ---- 
  
  res <- data.frame()
  #########   unadjusted model ---- 
  
  formula_unadjusted <- as.formula(paste(base_formula_unadjusted, " + ",outcome, sep = ""))
  
  perm_res <- adonis2(
    formula = formula_unadjusted,
    data = df_meta,
    permutations = 1000,
    method = "bray",
    by = "margin",  # alt: "terms"
    strata = df_meta$site,
    parallel = parallel::detectCores() - 1
  )
  
  perm_res
  
  # R2 
  r_sq <- perm_res[grep(outcome, rownames(perm_res)), "R2"] 
  # P-value 
  pval <- perm_res[grep(outcome, rownames(perm_res)), 'Pr(>F)'] 
  
  unadj <- data.frame(outcome = outcome, metric = "Bray-Curtis", r2_basic_adj = r_sq, pvalue_basic_adj = pval) 
  unadj
  
  # ------------------------------------------------------------------------------------------------------------------- # 
  
  ######   full model ---- 
  
  formula_adjusted <- as.formula(paste(base_formula_adjusted, " + ", outcome, sep = ""))
  
  perm_res <- adonis2(
    formula = formula_adjusted,
    data = df_meta,
    permutations = 1000,
    method = "bray",
    by = "margin",  # alt: "terms"
    strata = df_meta$site,
    parallel = parallel::detectCores() - 1
  )
  
  perm_res
  # R2 
  r_sq <- perm_res[grep(outcome, rownames(perm_res)), "R2"] 
  # P-value 
  pval <- perm_res[grep(outcome, rownames(perm_res)), 'Pr(>F)'] 
  
  adj <- data.frame(r2_full_adj = r_sq, pvalue_full_adj = pval) 
  adj
  
  res <- cbind(unadj,adj)
  res
  
  outcome_res <- res
  
  ######   Jaccard  ---- 
  
  res <- data.frame()
  #########   unadjusted model ---- 

  perm_res <- adonis2(
    formula = formula_unadjusted,
    data = df_meta,
    permutations = 1000,
    method = "jaccard",
    by = "margin",  # or "terms"
    strata = df_meta$site,
    parallel = parallel::detectCores() - 1
  )
  
  perm_res
  
  # R2 
  r_sq <- perm_res[grep(outcome, rownames(perm_res)), "R2"] 
  # P-value 
  pval <- perm_res[grep(outcome, rownames(perm_res)), 'Pr(>F)'] 
  
  unadj <- data.frame(outcome = outcome, metric = "Jaccard", r2_basic_adj = r_sq, pvalue_basic_adj = pval) 
  unadj
  
  # ------------------------------------------------------------------------------------------------------------------- # 
  
  ######   full model ---- 
  
  perm_res <- adonis2(
    formula = formula_adjusted,
    data = df_meta,
    permutations = 1000,
    method = "jaccard",
    by = "margin",  # or "terms"
    strata = df_meta$site,
    parallel = parallel::detectCores() - 1
  )
  
  perm_res
  # R2 
  r_sq <- perm_res[grep(outcome, rownames(perm_res)), "R2"] 
  # P-value 
  pval <- perm_res[grep(outcome, rownames(perm_res)), 'Pr(>F)'] 
  
  adj <- data.frame(r2_full_adj = r_sq, pvalue_full_adj = pval) 
  adj
  
  res <- cbind(unadj,adj)
  res
  
  # rbind bray and jaccard results 
  outcome_res <-  rbind(outcome_res, res)
  
  outcome_results <- rbind(outcome_results, outcome_res)
}

# ------------------------------------------------------------------------------------------------------------------- # 

# view results 
outcome_results

# save results 
path <- file.path(
  "/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_beta_diversity/by_polyp_type_beta_div_RESULTS.csv"
)

write.csv(outcome_results, 
          path,
          row.names = FALSE)

