## Author: 
    # Panayiotis Louca 

## Purpose of script: 
    #  

## Date Created: 
    # 05 May 2025 

## Notes: 
    #  

## Clear environment 
    rm(list = ls()) 

## Set seed 
    set.seed(1234)

## Set functions: 
      #  
      
## load up packages: 

    ### core 
    library(tidyverse)
    library(readxl)
    
    # mixed effect modelling 
    library(glmmTMB)
    library(lmerTest)
    library(broom.mixed)
    
    library(parallel)
    library(future.apply)
    
    library(Maaslin2)

# -------------------------------------------------------------------------- # 
  
    # ************************* # 
    #   IMPORT & PREP DATA   ---- 
    # ************************* # 
   
# -------------------------------------------------------------------------- #  

      ##   Dataset ---- 
    
    path <- file.path("/scratch/users/k2480753/NCL/P1_gene_families__DATASET.rds")
    df <- readRDS(path) %>% as.data.frame(.)
    
    # vectorise pathway column names 
    pathway_cols = names(df)[65:ncol(df)]
    
    # setup variable denoting all the outcomes 
    outcomes_to_test <- names(df)[
      grepl("^outcome", names(df)) | grepl("^SA[0-9]", names(df))
    ]
    
        # -------------------------------------------------------------------------- #  
    
    ###   check classes ---- 
    
    # all pathways are numeric 
    df %>% 
      select(all_of(pathway_cols)) %>%
      sapply(class) %>%
      table(useNA = 'ifany')
    
    df %>% 
      select(age, sex_imputed, BMI_imputed, ethnicity_imputed, site, years_education_imputed,symptom_or_screening_imputed, pack_years_smoking, alcohol_g_adj_imputed, calcium_mg_adj_imputed, nsp_g_adj_imputed, red_meat_g_adj_imputed) %>%
      glimpse() 
    
    df %>% 
      select(all_of(outcomes_to_test)) %>%
      glimpse()  

        # -------------------------------------------------------------------------- # 
    
    # setup covar vector  
    covars <- c('age', 'sex_imputed', 'BMI_imputed', 'ethnicity_imputed',
                'energy_kcal_imputed', 'HDI_imputed', # diet 
                'alcohol_g_adj_imputed','pack_years_smoking', # lifestyle 
                'ppi', 'aspirin_imputed') # medication 
    
    
    ###   Split data ---- 
    
    # gut microbiome pathway data only 
    dat_pathways <- df %>% 
      select(record_id, all_of(pathway_cols))
    
    dat_pathways <- dat_pathways %>%
      as.data.frame() %>%
      `rownames<-`(NULL) %>%  # forcefully drop all row names 
      tibble::column_to_rownames('record_id')
    
    # filter to metadata and outcomes 
    dat_meta <- df %>%
      select(record_id, site, batch, outcomes_to_test, 
             all_of(covars))
    
    dat_meta <-  # move iid to rownames 
      dat_meta %>%
      as.data.frame() %>%
      `rownames<-`(NULL) %>%  # forcefully drop all row names
      tibble::column_to_rownames('record_id')
    
    # ------------------------------------------------------------------------------------------------------------------- # 
    
    # Detect number of cores from SLURM 
    cores <- as.integer(Sys.getenv("SLURM_NTASKS"))
    
    # setup loop 
    base_dir <- "/scratch/users/k2480753/NCL/"
    
    # Loop through each variable 
     for (i in outcomes_to_test) {
     
       dir.create(paste(base_dir, 
                        i, sep = ""),
                   showWarnings = FALSE)
       
       out_dir <- paste(base_dir, 
                        i, sep = "")
      
      # Run Maaslin2 
      maaslin_results <- 
        Maaslin2(input_data = dat_pathways, 
                 input_metadata = dat_meta,
                 output = out_dir,
                 fixed_effects = c(i, covars),
                 random_effects = c('site', 'batch'),
                 plot_scatter = FALSE,
                 plot_heatmap = FALSE,
                 normalization = 'CLR', # CLR 
                 transform = "NONE", # NONE 
                 min_abundance = 0, # default is also 0 
                 max_significance = 0.1, # fdr < 0.1 is saved in sep file 
                 min_prevalence = 0.05, # 5% prevalence filter 
                 cores = cores
        )
     }    
    
    
    
    
    
    
    
    
    
    
    
    
    
    