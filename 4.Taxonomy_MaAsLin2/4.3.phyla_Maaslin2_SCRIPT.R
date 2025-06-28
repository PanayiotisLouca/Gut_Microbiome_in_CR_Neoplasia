## Author: 
#  Panayiotis Louca 

## Purpose of script: 
#  

## Date Created: 
#  01 August 2023 ------------------------------


## Notes: 
#  

## set working directory
setwd("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/")

## load up packages: 
library(tidyverse)
library(readxl)

library(Maaslin2)

library(skimr)

##                                  

rm(list = ls())

# ------------------------------------------------------------------------------------------------------------------- # 

##   IMPORT DATA ---- 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/dataset/CC_gm_polyp_dataset.rds")
dat <- read_rds(path) %>% as.data.frame(.)


# setup variable denoting all the outcomes 
outcomes_to_test =
  c(grep("outcome", names(dat), value = TRUE), # main outcomes 
    grep("SA", names(dat), value = TRUE) # secondary outcomes/sensitivity analysis 
    ) 
outcomes_to_test

# setup variable denoting all models to test 
standard_model <- c('age', 'sex_imputed', 'BMI_imputed', 'ethnicity_imputed')
full_model <- c(standard_model,
                 'energy_kcal_imputed', 'HDI_imputed', # diet 
                  'alcohol_g_adj_imputed','pack_years_smoking', # lifestyle 
                  'ppi', 'aspirin_imputed') # medication 


    # -------------------------------------------------------------------------- #  

# setup combined models  
full_model_any_neoplasia <- c(full_model,
                              'carbohydrate_component_score_imputed',
                              'height_imputed',
                              'previous_bowel_cancer',
                              'previous_polyps_imputed',
                              'protein_g_adj_imputed',
                              'red_meat_g_adj_imputed',
                              'waist_circumference_imputed',
                              'weight_imputed',
                              'employment',
                              'fish_oil_imputed',
                              'metformin_imputed',
                              'statin_imputed',
                              'hypertension_imputed',
                              'non_milk_extrinsic_sugars_component_score_imputed',
                              'nsp_g_adj_imputed',
                              'pulses_and_nuts_component_score_imputed',
                              'years_education_imputed',
                              "fibre_component_score_imputed",                    
                              'calcium_component_score_imputed',
                              'fh_crc_imputed',
                              'fish_component_score_imputed',
                              'cholesterol_component_score_imputed',
                              'other_cancers_imputed',
                              'protein_component_score_imputed',
                              'calcium_mg_adj_imputed',
                              'red_meat_component_score_imputed',
                              'bowel_prep_imputed',
                              'sat_fat_component_score_imputed',
                              'fruit_and_veg_component_score_imputed',
                              'pufa_component_score_imputed',
                              'antibiotic_use_imputed',
                              'calcium_vitd_imputed')
                              
                              
full_model_cancer_detected <- c(full_model,
                                'weight_imputed',
                                'nsp_g_adj_imputed')
                              
                              


full_model_advanced_neoplasia_exc_cancer <- c(full_model,
                                              'nsp_g_adj_imputed',
                                              'protein_g_adj_imputed',
                                              'waist_circumference_imputed',
                                              'weight_imputed',
                                              'carbohydrate_component_score_imputed',
                                              'non_milk_extrinsic_sugars_component_score_imputed',
                                              'red_meat_g_adj_imputed',
                                              'height_imputed',
                                              'pulses_and_nuts_component_score_imputed',
                                              'fish_component_score_imputed',
                                              'metformin_imputed',
                                              'other_cancers_imputed')


full_model_non_advanced_neoplasia <- c(full_model,
                                       'carbohydrate_component_score_imputed',
                                       'diabetes',
                                       'fish_oil_imputed',
                                       'hypertension_imputed',
                                       "metformin_imputed",
                                       "nsp_g_adj_imputed",
                                       "previous_bowel_cancer",
                                       "previous_polyps_imputed",
                                       "pulses_and_nuts_component_score_imputed",
                                       'red_meat_g_adj_imputed',
                                       "waist_circumference_imputed",
                                       "weight_imputed",
                                       "years_education_imputed",
                                       "calcium_component_score_imputed",
                                       "employment",
                                       "fh_crc_imputed",
                                       "fibre_component_score_imputed",
                                       "fish_component_score_imputed",
                                       "height_imputed",
                                       "non_milk_extrinsic_sugars_component_score_imputed",
                                       "protein_component_score_imputed",
                                       "protein_g_adj_imputed",
                                       "red_meat_component_score_imputed",
                                       "sat_fat_component_score_imputed",
                                       "statin_imputed",
                                       "bowel_prep_imputed",
                                       "cholesterol_component_score_imputed",
                                       "fruit_and_veg_component_score_imputed",
                                       "liver_disease",
                                       "antibiotic_use_imputed",
                                       "other_cancers_imputed",
                                       "calcium_mg_adj_imputed",
                                       "calcium_vitd_imputed",
                                       "pufa_component_score_imputed",
                                       "anti_obesity_surgery_imputed")


full_model_only_adenoma <- c(full_model,
                             "carbohydrate_component_score_imputed",
                             "employment",
                             "height_imputed",
                             "previous_bowel_cancer",
                             "waist_circumference_imputed",
                             "weight_imputed",
                             "fibre_component_score_imputed",
                             "metformin_imputed",
                             "previous_polyps_imputed",
                             "red_meat_g_adj_imputed",
                             "statin_imputed",
                             "nsp_g_adj_imputed",
                             "fish_oil_imputed",
                             "hypertension_imputed",
                             "protein_g_adj_imputed",
                             "other_cancers_imputed",
                             "pulses_and_nuts_component_score_imputed",
                             "non_milk_extrinsic_sugars_component_score_imputed",
                             "years_education_imputed",
                             "calcium_component_score_imputed",
                             "fish_component_score_imputed",
                             "protein_component_score_imputed",
                             "red_meat_component_score_imputed",
                             "sat_fat_component_score_imputed",
                             "fh_crc_imputed",
                             "cholesterol_component_score_imputed",
                             "fruit_and_veg_component_score_imputed",
                             "bowel_prep_imputed",
                             "antibiotic_use_imputed",
                             "calcium_mg_adj_imputed",
                             "pufa_component_score_imputed")          
                             
full_model_only_serrated_lesions <- c(full_model,
                                      'employment',
                                      "nsp_g_adj_imputed",
                                      "waist_circumference_imputed",
                                      "weight_imputed",
                                      "protein_g_adj_imputed",
                                      "previous_polyps_imputed",
                                      "fruit_and_veg_component_score_imputed",
                                      "pulses_and_nuts_component_score_imputed",
                                      "red_meat_g_adj_imputed",
                                      "cholesterol_component_score_imputed",
                                      "height_imputed",
                                      "hypertension_imputed",
                                      "fibre_component_score_imputed",
                                      "calcium_component_score_imputed",
                                      "years_education_imputed")

                                      

full_model_only_ssl <- c(full_model,
                         'weight_imputed',
                         "nsp_g_adj_imputed",
                         "protein_g_adj_imputed")


full_model_advanced_neoplasia_inc_cancer <- c(full_model,
                                              'protein_g_adj_imputed',
                                              'waist_circumference_imputed',
                                              'weight_imputed',
                                              "nsp_g_adj_imputed",
                                              "metformin_imputed",
                                              "carbohydrate_component_score_imputed",
                                              "non_milk_extrinsic_sugars_component_score_imputed",
                                              "height_imputed",
                                              "red_meat_g_adj_imputed",
                                              "other_cancers_imputed",
                                              "fibre_component_score_imputed",
                                              "pulses_and_nuts_component_score_imputed")


# vectorise all unique vars across models 
all_vars = unique(c(full_model_any_neoplasia,
                    full_model_cancer_detected,
                    full_model_advanced_neoplasia_exc_cancer,
                    full_model_non_advanced_neoplasia,
                    full_model_only_adenoma,
                    full_model_only_serrated_lesions,
                    full_model_only_ssl,
                    full_model_advanced_neoplasia_inc_cancer))

    # -------------------------------------------------------------------------- #  

###   check classes ---- 

dat %>% # all gut microbiome vars are numeric 
  select(
    starts_with("genus_"),
    UNCLASSIFIED, 
    alpha_div_shannon, 
    alpha_div_simpson,
    alpha_div_taxonomic_richness) %>%
  sapply(class) %>%
  table() 


dat %>% 
  select(all_of(all_vars)) %>%
  skim() 

dat %>% #  factors 
  select(all_of(outcomes_to_test)
         ) %>%
  glimpse()  


# ------------------------------------------------------------------------------------------------------------------- # 

names(dat)

folder_names <- c("Standard_model_Age_sex_BMI_ethnic", 
                  "Full_model")

model_vars <- list(
  standard_model = standard_model,
  full_model = full_model)

    # -------------------------------------------------------------------------- #  

# make folders if not already present 

# Define outcomes that should have the combined_model folder
combined_model_outcomes <- c(
  "outcome1_any_neoplasia_detected",
  "outcome1_cancer",
  "outcome2_advanced_neoplasia_exc_cancer",
  "outcome3_non_advanced_neoplasia",
  "outcome4_only_adenomatous",
  "outcome5_only_serrated_lesions",
  "outcome6_only_sessile_serrated_lesions",
  "SA1_advanced_neoplasia_inc_cancer"
)

# # Make folders if not already present 
 for (outcome in outcomes_to_test) {

   base_dir <- "/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_Taxonomy/Phyla_results/5%_prev"


#   # Define the base directory path for each outcome
   outcome_path <- file.path(base_dir, outcome)

   # Create the outcome directory if it doesn't exist
   if (!dir.exists(outcome_path)) {
     dir.create(outcome_path)
     cat("Created outcome folder:", outcome_path, "\n")
   }

#   # Loop through each folder in folder_names to create subdirectories
   for (folder in folder_names) {
     subfolder_path <- file.path(outcome_path, folder)
     if (!dir.exists(subfolder_path)) {
       dir.create(subfolder_path)
       cat("Created subfolder:", subfolder_path, "\n")
     }
   }

   # Create combined_model folder only for specified outcomes
   if (outcome %in% combined_model_outcomes) {
     combined_folder_path <- file.path(outcome_path, "combined_model")
     if (!dir.exists(combined_folder_path)) {
       dir.create(combined_folder_path)
       cat("Created combined_model folder:", combined_folder_path, "\n")
     }
   }
 }


# *********************************************************************************************** # 
#           Run loop to analyse all outcomes and all models sequentially                       ---- 
# *********************************************************************************************** # 

for (outcome in outcomes_to_test) {
  
  for (folder in 1:length(folder_names)) {

  ###   split data ---- 
  
  # gut microbiome data only 
  dat_microb <- dat %>% 
    filter(!is.na(!!sym(outcome))
    ) %>%
    select(record_id,
           starts_with("phyla_")
    )
  
  # convert microbe data to decimal proportions 
  dat_microb[ ,-1] <- sapply(dat_microb[ ,-1], function(x) (x/100))
  
  dat_microb <- # move iid to rownames 
    dat_microb %>%
    tibble::column_to_rownames('record_id')
  
  # set prevalence threshold (5%) 
  prevalence_threshold <- 0.05
  
  # filter microbes present in at least prevalence_threshold of samples
  microb_prevalence <- colMeans(dat_microb > 0)  # proportion of samples with non-zero abundance per microbe
  dat_microb <- dat_microb[ , microb_prevalence >= prevalence_threshold]  # keep microbes above threshold
  
  
  # filter to metadata and outcomes 
  dat_meta <- dat %>%
    filter(!is.na(!!sym(outcome))) %>%
    select(record_id,
           site,
           batch,
           all_of(all_vars),
           !!sym(outcome)) 
  
  dat_meta <-  # move iid to rownames 
    dat_meta %>%
    tibble::column_to_rownames('record_id')
  
  # ------------------------------------------------------------------------------------------------------------------- # 
  
  ###   Run Maaslin2 ---- 
  
  #####  standard model (site + age, sex, BMI, ethnicity) ---- 
  
  # Construct the fixed effects vector
  
  # run maaslin2 
  maaslin_results <-
    Maaslin2(
      input_data = dat_microb,
      input_metadata = dat_meta,
      output = paste(
        "/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_Taxonomy/Phyla_results/5%_prev/",
        outcome,
        "/",
        folder_names[folder],
        sep = ""),
      fixed_effects = c(outcome, model_vars[[folder]]),
      random_effects = c('site', 'batch'),
      plot_scatter = FALSE,
      plot_heatmap = FALSE,
      normalization = 'CLR', # CLR 
      transform = "NONE", # NONE 
      min_abundance = 0, # default is also 0 
      max_significance = 0.1, # fdr < 0.1 is saved in sep file 
      min_prevalence = 0,
      cores = parallel::detectCores()-1
    )
  
  }
}

# ------------------------------------------------------------------------------------------------------------------- # 

  ######   RUN AGAIN ON THOSE THAT HAVE COMBINED MODELS   ---- 

outcome_model_map <- list(
  "outcome1_any_neoplasia_detected" = full_model_any_neoplasia,
  "outcome1_cancer" = full_model_cancer_detected,
  "outcome2_advanced_neoplasia_exc_cancer" = full_model_advanced_neoplasia_exc_cancer, 
  "outcome3_non_advanced_neoplasia" = full_model_non_advanced_neoplasia,
  "outcome4_only_adenomatous" = full_model_only_adenoma,
  "outcome5_only_serrated_lesions" = full_model_only_serrated_lesions,
  "outcome6_only_sessile_serrated_lesions" = full_model_only_ssl,
  "SA1_advanced_neoplasia_inc_cancer" = full_model_advanced_neoplasia_inc_cancer)
  


for (outcome in combined_model_outcomes) {
  
    ###   split data ---- 
    
    # gut microbiome data only 
    dat_microb <- dat %>% 
      filter(!is.na(!!sym(outcome))
      ) %>%
      select(record_id,
             starts_with("phyla_")
      )
    
    # convert microbe data to decimal proportions 
    dat_microb[ ,-1] <- sapply(dat_microb[ ,-1], function(x) (x/100))
    
    dat_microb <- # move iid to rownames 
      dat_microb %>%
      tibble::column_to_rownames('record_id')
    
    
    # set prevalence threshold (5%) 
    prevalence_threshold <- 0.05
    
    # filter microbes present in at least prevalence_threshold of samples
    microb_prevalence <- colMeans(dat_microb > 0)  # proportion of samples with non-zero abundance per microbe
    dat_microb <- dat_microb[ , microb_prevalence >= prevalence_threshold]  # keep microbes above threshold
    
    
    # filter to metadata and outcomes 
    dat_meta <- dat %>%
      filter(!is.na(!!sym(outcome))) %>%
      select(record_id,
             site,
             batch,
             all_of(all_vars),
             !!sym(outcome)) 
    
    dat_meta <-  # move iid to rownames 
      dat_meta %>%
      tibble::column_to_rownames('record_id')
    
    # ------------------------------------------------------------------------------------------------------------------- # 
    
    ###   Run Maaslin2 ---- 
    
    #####  standard model (site + age, sex, BMI, ethnicity) ---- 

    
    model = outcome_model_map[[outcome]]
    
    
    # run maaslin2 
    maaslin_results <-
      Maaslin2(
        input_data = dat_microb,
        input_metadata = dat_meta,
        output = paste(
          "/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_Taxonomy/Phyla_results/5%_prev/",
          outcome,
          "/combined_model",
          sep = ""),
        fixed_effects = c(outcome, model),
        random_effects = c('site', 'batch'),
        plot_scatter = FALSE,
        plot_heatmap = FALSE,
        normalization = 'CLR', # CLR 
        transform = "NONE", # NONE 
        min_abundance = 0, # default is also 0 
        max_significance = 0.1, # fdr < 0.1 is saved in sep file 
        min_prevalence = 0,
        cores = parallel::detectCores()-1
      )
    
  }


# ------------------------------------------------------------------------------------------------------------------- # 

