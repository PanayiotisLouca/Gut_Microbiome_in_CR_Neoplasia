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
require(tidyverse)
require(readxl)
#library(microbiomics)
library(Maaslin2)
#library(pkgconfig)
#require(skimr)
##                                  

rm(list = ls())

# ------------------------------------------------------------------------------------------------------------------- # 

##   IMPORT DATA ---- 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_functionality/P1_Pathways/dataset/P1_pathways_DATASET.rds")
dat <- read_rds(path)

# setup variable denoting all the outcomes 
outcomes_to_test <- names(dat)[
  grepl("^outcome", names(dat)) | grepl("^SA[0-9]", names(dat))
]

outcomes_to_test

names(dat)

# vectorise pathway names 
pathway_names = names(dat)[65:ncol(dat)] # all pathway names 


    # -------------------------------------------------------------------------- #    

###   check classes ---- 

# all pathways are numeric 
dat %>% 
  select(all_of(pathway_names)) %>%
  sapply(class) %>%
  table(useNA = 'ifany')


dat %>% 
  select(age, sex_imputed, BMI_imputed, ethnicity_imputed, site, years_education_imputed,symptom_or_screening_imputed, pack_years_smoking, alcohol_g_adj_imputed, calcium_mg_adj_imputed, nsp_g_adj_imputed, red_meat_g_adj_imputed) %>%
  glimpse() 

dat %>% 
  select(all_of(outcomes_to_test)) %>%
  glimpse()  

# ------------------------------------------------------------------------------------------------------------------- # 

# setup covar vector  

covars <- c('age', 'sex_imputed', 'BMI_imputed', 'ethnicity_imputed',
            'energy_kcal_imputed', 'HDI_imputed', # diet 
            'alcohol_g_adj_imputed','pack_years_smoking', # lifestyle 
            'ppi', 'aspirin_imputed') # medication 

names(dat)

###   Split data ---- 

# gut microbiome pathway data only 
dat_pathways <- dat %>% 
  select(record_id, all_of(pathway_names))

dat_pathways <- dat_pathways %>%
  as.data.frame() %>%
  `rownames<-`(NULL) %>%  # forcefully drop all row names
  tibble::column_to_rownames('record_id')

# filter to metadata and outcomes 
dat_meta <- dat %>%
  select(record_id, site, batch, outcomes_to_test, 
         all_of(covars))

dat_meta <-  # move iid to rownames 
  dat_meta %>%
  as.data.frame() %>%
  `rownames<-`(NULL) %>%  # forcefully drop all row names
  tibble::column_to_rownames('record_id')

# ------------------------------------------------------------------------------------------------------------------- # 

# run loop for outcomes and GLM maaslin2 models 
for(i in outcomes_to_test) {

dir.create(paste("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_functionality/P1_Pathways/Maaslin2_results/",
                 i, sep = ""),
           showWarnings = FALSE)
  
# run maaslin2 
maaslin_results <- 
  Maaslin2(input_data = dat_pathways, 
           input_metadata = dat_meta,
           output = paste("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_functionality/P1_Pathways/Maaslin2_results/",
                          i, sep = ""),
           fixed_effects = c(i, covars),
           random_effects = c('site', 'batch'),
           plot_scatter = FALSE,
           plot_heatmap = FALSE,
           normalization = 'CLR', # CLR 
           transform = "NONE", # NONE 
           min_abundance = 0, # default is also 0 
           max_significance = 0.1, # fdr < 0.1 is saved in sep file 
           min_prevalence = 0.05, # 5% prevalence filter 
           cores = parallel::detectCores()-1
  )

}








