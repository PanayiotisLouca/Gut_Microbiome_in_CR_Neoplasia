## Author: 
    # Panayiotis Louca 

## Purpose of script: 
    #  

## Date Created: 
    # 08 June 2025 

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
    
    # for matching cases and controls 
    library(MatchIt)
    
    # for random forest & parallelisation 
    library(caret)
    library(ranger)
    library(pROC)
    
# -------------------------------------------------------------------------- # 
  
    # ************************* # 
    #   IMPORT & PREP DATA   ---- 
    # ************************* # 
   
# -------------------------------------------------------------------------- #  

      ##   Dataset ---- 
    path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/dataset/CC_gm_polyp_dataset.rds")
    dat <- read_rds(path) %>% as.data.frame(.)
    
# import risk  data 
    path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_high_risk_prediction/polyp_risk_data.csv")
    dat.risk <- read.csv(path)

    table(dat.risk$high_risk)        

     # subset dataset 
    dat.risk <- dat.risk %>%
      select(record_id, high_risk)

    table(dat$record_id %in% dat.risk$record_id)

    # merge  
    dat <- merge(dat, dat.risk,
                 by = "record_id",
                 all.x = TRUE)  

    table(dat$high_risk, useNA = 'ifany')    

    
        # -------------------------------------------------------------------------- #  
    
    # Create a new variable for matching pool: high_risk + neoplasia status 
    dat$eligible_control <- ifelse(dat$high_risk == FALSE & dat$outcome1_any_neoplasia_detected == 0, TRUE, FALSE)
    
    # Subset to only high_risk == TRUE and eligible controls 
    dat_match <- dat[dat$high_risk == TRUE | dat$eligible_control == TRUE, ]
    
    # Run propensity score matching 
    m.out <- matchit(
      high_risk ~ age + sex_imputed + BMI_imputed, # matching on age, sex, and BMI 
      data = dat_match,
      method = "nearest",         # Nearest neighbor matching 
      ratio = 1,                  # 1:1 matching 
      caliper = 0.2,               #restrict matches to 0.2 * SD of PS 
      std.caliper = TRUE 
    )
    
    # Examine matching summary 
    summary(m.out)
    
    # Extract matched dataset 
    matched_data <- match.data(m.out) 
    
    # Check match counts 
    table(matched_data$high_risk)

    write.csv(matched_data, 
              "/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_high_risk_prediction/P1_matched_high_risk_n_414.csv",
              row.names = FALSE)
    
        # -------------------------------------------------------------------------- #  
    
# *********************************************************** # 
##   Select data for gut microbiome and covariates         ---- 
# *********************************************************** # 

    # non gut microbiome variables 
  covars <- c( #'age', # remove age, sex, and BMI as these are matched upon 
               # 'sex_imputed',
               # 'BMI_imputed',
              'ethnicity_imputed',
              'years_education_imputed',
              'weight_imputed',
              'height_imputed',
              'waist_circumference_imputed',
              'bowel_prep_imputed',
              'previous_polyps_imputed',
              'previous_bowel_cancer',
              'other_cancers_imputed',
              'fh_crc_imputed',
              'symptom_or_screening_imputed',
              'diabetes',
              'angina_imputed',
              'hypertension_imputed',
              'liver_disease',
              'anti_obesity_surgery_imputed',
              'aspirin_imputed',
              'metformin_imputed',
              'ppi',
              'statin_imputed',
              'antibiotic_use_imputed',
              'fish_oil_imputed',
              'calcium_vitd_imputed',
              'alcohol_g_adj_imputed',
              'pack_years_smoking',
              'energy_kcal_imputed',
              'HDI_imputed',
              'nsp_g_imputed',
              'protein_g_adj_imputed',
              'calcium_mg_adj_imputed',
              'red_meat_g_adj_imputed')
    
    
    # select non gut microbiome data alongside outcome 
    dat_non_gm <- matched_data %>% # all gut microbiome vars are numeric 
      select(record_id,
             high_risk,
             all_of(covars),
             contains("component_score"))
    
    names(dat_non_gm)         
        
    # Select gut microbiome data 
    dat_gm <- matched_data %>% 
      select(record_id,
             high_risk,
             starts_with("species_"))

    names(dat_gm)    

    # merge datasets 
    dat <- merge(dat_gm, dat_non_gm,
                 by = c("record_id",
                        "high_risk"))
    
    # drop record id 
    dat <- dat %>%
      select(-record_id)

    # -------------------------------------------------------------------------- #  
    
    # ************************* # 
    ##   Prepare data        ---- 
    # ************************* # 
    
    # Rename and convert response to factor
    data <- dat %>%
      rename(response = high_risk) %>%
      mutate(response = case_when(response == TRUE ~ "high",
                                  .default = "low"),
             response = factor(response, levels = c("low", "high")))
    
    # Remove zero variance predictors 
    vars_no_outcome <- names(data)[names(data) != "response"] # Get the names except response 
    
    # Run nearZeroVar on subset without response 
    nzv <- nearZeroVar(data %>% select(all_of(vars_no_outcome)))
    
    # get names of variables with near zero variance 
    nzv_names <- vars_no_outcome[nzv]
    nzv_names
    
    if (length(nzv_names) > 0) data <- data %>% select(-all_of(nzv_names))
    
    # -------------------------------------------------------------------------- #
    
    # ************************* #
    ##   5-fold CV with 80/20 split ----
    # ************************* #
    
    set.seed(1234)
    
    cat("dataset size:", nrow(data), "\n")
    cat("Class distribution:\n", table(data$response))
    
    
    # Set up stratified 5-fold CV on data 
    set.seed(1234)
    cv_folds <- createFolds(data$response, k = 5, list = TRUE, returnTrain = FALSE)
    
   
    cv_results <- data.frame(
      fold = integer(),
      auc = numeric(),
      sensitivity = numeric(),
      specificity = numeric(),
      accuracy = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Store predictions for each fold
    all_predictions <- data.frame(
      obs = factor(),
      pred = factor(),
      high = numeric(),
      low = numeric(),
      fold = integer()
    )
    
    # Store variable importance across folds 
    var_importance_list <- list()
    
    # Perform 5-fold CV 
    for (i in 1:5) { 
      
      cat("Processing fold", i, "of 5...\n")
      
      # Split data for this fold 
      test_fold <- data[cv_folds[[i]], ] 
      train_fold <- data[-cv_folds[[i]], ]
      
      # Train model on this fold 
      rf_model <- ranger(response ~ ., data = train_fold, num.trees = 500, mtry = floor(sqrt(ncol(train_fold) - 1)), 
                         importance = "impurity", probability = TRUE, seed = 1234 + i )
      
            # Extract variable importance 
            imp <- rf_model$variable.importance
            imp_df <- data.frame(variable = names(imp), importance = as.numeric(imp), fold = i)
            
            # Store 
            var_importance_list[[i]] <- imp_df
            
      # Predict on test fold 
      fold_pred <- predict(rf_model, test_fold) 
      fold_probs <- fold_pred$predictions
      
      # Calculate ROC and AUC 
      roc_obj <- roc(test_fold$response, fold_probs[, "high"], quiet = TRUE) 
      fold_auc <- as.numeric(auc(roc_obj))
      
      # Get optimal threshold 
      best_thresh <- coords(roc_obj, "best", ret = "threshold")$threshold 
      fold_pred_class <- ifelse(fold_probs[, "high"] > best_thresh, "high", "low") 
      fold_pred_class <- factor(fold_pred_class, levels = c("low", "high"))
      
      # Calculate performance metrics 
      cm <- confusionMatrix(fold_pred_class, test_fold$response)
      
      # Store results 
      cv_results <- rbind(cv_results, data.frame(fold = i, auc = fold_auc, 
                                                 sensitivity = cm$byClass["Sensitivity"],
                                                 specificity = cm$byClass["Specificity"], 
                                                 accuracy = cm$overall["Accuracy"] ))
      
      # Store predictions for aggregated analysis 
      fold_preds <- data.frame(obs = test_fold$response,
                               pred = fold_pred_class, 
                               high = fold_probs[, "high"], 
                               low = fold_probs[, "low"], fold = i ) 
      
      all_predictions <- rbind(all_predictions, fold_preds)
      
      }
    
    # Combine all importance data frames
    all_importance <- bind_rows(var_importance_list)
    
    # Average importance across folds 
    mean_importance <- all_importance %>%
      group_by(variable) %>%
      summarize(mean_importance = mean(importance), .groups = "drop") %>%
      arrange(desc(mean_importance)) %>%
      slice(1:20)  # Top 20 variables
    
        # -------------------------------------------------------------------------- #  
    
    # clean feature names 
    mean_importance = mean_importance %>%
      mutate(variable = variable %>%
               str_replace_all("_", " ") %>% # Replace underscores 
               str_remove_all("\\bimputed\\b") %>% # Remove the word "imputed"
               str_squish() %>% # Remove extra whitespace 
               str_replace("^species (.+)", "(species) \\1") %>% # Wrap species names 
               str_replace_all("(\\b\\w+ meat|alcohol|nsp|fibre|fiber|carbohydrate|energy|red meat|protein|processed meat|vegetables|fruit) g", "\\1 (g)") %>% 
               str_replace_all("(fat|sugar|salt) mg", "\\1 (mg)") %>%  # Add (mg) 
               str_to_sentence() %>%  # Apply sentence case to everything first 
               str_replace("\\(species\\)", "(Species)") # Fix the species wrapper 
      ) %>%
      # Handle species case formatting 
      mutate(variable = case_when(
        # GGB taxa should be uppercase 
        str_detect(variable, "\\(Species\\) ggb") ~ 
          str_replace(variable, "\\(Species\\) (.+)", function(x) {
            species_part <- str_extract(x, "(?<=\\(Species\\) ).*")
            paste0("(Species) ", str_to_upper(species_part))
          }),
        # Other species should be title case 
        str_detect(variable, "\\(Species\\)") ~ 
          str_replace(variable, "\\(Species\\) (.+)", function(x) {
            species_part <- str_extract(x, "(?<=\\(Species\\) ).*")
            paste0("(Species) ", str_to_title(species_part))
          }),
        # Everything else stays as is 
        TRUE ~ variable
      ))
    
    # -------------------------------------------------------------------------- #
    
    # *************************************** #
    ##   Enhanced Feature Importance Plot  ----
    # *************************************** #
    
    # Plot top variables 
    feat_imp_gut_microbiome_and_covars <- ggplot(mean_importance, aes(x = reorder(variable, mean_importance), y = mean_importance)) +
      geom_col(fill = "#1F77B4", width = 0.7) +
      coord_flip() +
      labs(
        title = "Variable Importance – Gut microbiome & covariates",
        x = NULL,
        y = "Importance (Gini decrease)"
      ) +
      theme_classic() +
      theme(
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.grid.minor = element_line(color = "gray95", size = 0.1),
        plot.margin = margin(10, 20, 10, 10)
      )
    
    feat_imp_gut_microbiome_and_covars
    
    # Summarise CV results
    cat("\n=== 5-Fold Cross-Validation Results ===\n")
    cat("Mean AUC:", round(mean(cv_results$auc), 3), "±", round(sd(cv_results$auc), 3), "\n")
    cat("AUC Range:", round(min(cv_results$auc), 3), "-", round(max(cv_results$auc), 3), "\n")
    cat("Mean Sensitivity:", round(mean(cv_results$sensitivity, na.rm = TRUE), 3), "±", round(sd(cv_results$sensitivity, na.rm = TRUE), 3), "\n")
    cat("Mean Specificity:", round(mean(cv_results$specificity, na.rm = TRUE), 3), "±", round(sd(cv_results$specificity, na.rm = TRUE), 3), "\n")
    cat("Mean Accuracy:", round(mean(cv_results$accuracy), 3), "±", round(sd(cv_results$accuracy), 3), "\n")
    
    
    # -------------------------------------------------------------------------- #
    
    # ************************* #
    ##   AUC Plot   ----
    # ************************* #
    
    # Get ROC for each fold
    roc_list <- all_predictions %>%
      group_split(fold) %>%
      map(~ roc(response = .x$obs, predictor = .x$high, quiet = TRUE))
    
    # Interpolate sensitivity at common FPR points
    fpr_grid <- seq(0, 1, length.out = 100)
    
    roc_interp <- map_dfr(seq_along(roc_list), function(i) {
      r <- roc_list[[i]]
      sens <- r$sensitivities
      fpr <- 1 - r$specificities
      
      # Create a data frame to ensure alignment
      roc_df <- tibble(fpr = fpr, sens = sens) %>%
        filter(!is.na(fpr), !is.na(sens)) %>%
        distinct(fpr, .keep_all = TRUE) %>%
        arrange(fpr)
      
      # Ensure at least two points for interpolation
      if (nrow(roc_df) < 2) {
        return(tibble(fpr = fpr_grid, sens = NA_real_, fold = as.character(i)))
      }
      
      # Interpolate sensitivities at common FPR grid
      tibble(
        fpr = fpr_grid,
        sens = approx(roc_df$fpr, roc_df$sens, xout = fpr_grid, rule = 2)$y,
        fold = as.character(i)
      )
    }, .id = NULL) %>%
      group_by(fpr) %>%
      summarise(
        mean_sens = mean(sens, na.rm = TRUE),
        min_sens = min(sens, na.rm = TRUE),
        max_sens = max(sens, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Create AUC plot 
    mean_auc <- mean(cv_results$auc)
    range_auc <- range(cv_results$auc)
    
    p_gut_microbiome_and_covars <- ggplot(roc_interp, aes(x = fpr)) +
      geom_ribbon(aes(ymin = min_sens, ymax = max_sens), fill = "gray80", alpha = 0.3) +
      geom_line(aes(y = mean_sens), color = "#1F77B4", size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
      labs(
        x = "False Positive Rate (1 - Specificity)",
        y = "Sensitivity",
        title = "ROC Curve – Gut microbiome & covariates",
        subtitle = sprintf("5-fold Cross-Validation: Mean AUC: %.3f (Range: %.3f–%.3f)", mean_auc, range_auc[1], range_auc[2])
      )  +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      theme_classic() +
      theme(
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.grid.minor = element_line(color = "gray95", size = 0.1),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      coord_fixed(ratio = 1)
    
    print(p_gut_microbiome_and_covars)
    
    
    # ----------------------------------------------------------------------------------------------------------------------- #      
    
    # ************************************************** # 
    ##   Select data for only gut microbiome          ---- 
    # ************************************************** # 
    
   
    # Select gut microbiome data 
    dat <- matched_data %>% 
      select(record_id,
             high_risk,
             starts_with("species_"))
    
    names(dat)    
    
    # drop record id 
    dat <- dat %>%
      select(-record_id)
    
    # -------------------------------------------------------------------------- #  
  
    # ************************* # 
    ##   Prepare data        ---- 
    # ************************* # 
    
    # Rename and convert response to factor
    data <- dat %>%
      rename(response = high_risk) %>%
      mutate(response = case_when(response == TRUE ~ "high",
                                  .default = "low"),
             response = factor(response, levels = c("low", "high")))
    
    # Remove zero variance predictors 
    vars_no_outcome <- names(data)[names(data) != "response"] # Get the names except response 
    
    # Run nearZeroVar on subset without response
    nzv <- nearZeroVar(data %>% select(all_of(vars_no_outcome)))
    
    # get names of variables with near zero variance 
    nzv_names <- vars_no_outcome[nzv]
    nzv_names
    
    if (length(nzv_names) > 0) data <- data %>% select(-all_of(nzv_names))
    
    # -------------------------------------------------------------------------- #
    
    # ************************* #
    ##   5-fold CV with 80/20 split ----
    # ************************* #
    
    # Create 80/20 split for final evaluation
    set.seed(1234)
    
    cat("dataset size:", nrow(data), "\n")
    cat("Class distribution:\n", table(data$response))
    
    
    set.seed(1234)
    cv_folds <- createFolds(data$response, k = 5, list = TRUE, returnTrain = FALSE)
    
    cv_results <- data.frame(
      fold = integer(),
      auc = numeric(),
      sensitivity = numeric(),
      specificity = numeric(),
      accuracy = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Store predictions for each fold
    all_predictions <- data.frame(
      obs = factor(),
      pred = factor(),
      high = numeric(),
      low = numeric(),
      fold = integer()
    )
    
    
    # Store variable importance across folds
    var_importance_list <- list()
    
    # Perform 5-fold CV 
    for (i in 1:5) { 
      
      cat("Processing fold", i, "of 5...\n")
      
      # Split data for this fold 
      test_fold <- data[cv_folds[[i]], ] 
      train_fold <- data[-cv_folds[[i]], ]
      
      # Train model on this fold 
      rf_model <- ranger(response ~ ., data = train_fold, num.trees = 500, mtry = floor(sqrt(ncol(train_fold) - 1)), 
                         importance = "impurity", probability = TRUE, seed = 1234 + i )
      
      # Extract variable importance
      imp <- rf_model$variable.importance
      imp_df <- data.frame(variable = names(imp), importance = as.numeric(imp), fold = i)
      
      # Store 
      var_importance_list[[i]] <- imp_df
      
      # Predict on test fold 
      fold_pred <- predict(rf_model, test_fold) 
      fold_probs <- fold_pred$predictions
      
      # Calculate ROC and AUC 
      roc_obj <- roc(test_fold$response, fold_probs[, "high"], quiet = TRUE) 
      fold_auc <- as.numeric(auc(roc_obj))
      
      # Get optimal threshold 
      best_thresh <- coords(roc_obj, "best", ret = "threshold")$threshold 
      fold_pred_class <- ifelse(fold_probs[, "high"] > best_thresh, "high", "low") 
      fold_pred_class <- factor(fold_pred_class, levels = c("low", "high"))
      
      # Calculate performance metrics 
      cm <- confusionMatrix(fold_pred_class, test_fold$response)
      
      # Store results 
      cv_results <- rbind(cv_results, data.frame(fold = i, auc = fold_auc, 
                                                 sensitivity = cm$byClass["Sensitivity"],
                                                 specificity = cm$byClass["Specificity"], 
                                                 accuracy = cm$overall["Accuracy"] ))
      
      # Store predictions for aggregated analysis 
      fold_preds <- data.frame(obs = test_fold$response,
                               pred = fold_pred_class, 
                               high = fold_probs[, "high"], 
                               low = fold_probs[, "low"], fold = i ) 
      
      all_predictions <- rbind(all_predictions, fold_preds)
      
    }
    
    # Combine all importance data frames
    all_importance <- bind_rows(var_importance_list)
    
    # Average importance across folds 
    mean_importance <- all_importance %>%
      group_by(variable) %>%
      summarize(mean_importance = mean(importance), .groups = "drop") %>%
      arrange(desc(mean_importance)) %>%
      slice(1:20)  # Top 20 variables 
    
    # -------------------------------------------------------------------------- #  
    
    # clean feature names 
    mean_importance = mean_importance %>%
      mutate(variable = variable %>%
               str_replace_all("_", " ") %>% # Replace underscores 
               str_remove_all("\\bimputed\\b") %>% # Remove the word "imputed" 
               str_squish() %>% # Remove extra whitespace 
               str_replace("^species (.+)", "(species) \\1") %>% # Wrap species names 
               str_replace_all("(\\b\\w+ meat|alcohol|nsp|fibre|fiber|carbohydrate|energy|red meat|protein|processed meat|vegetables|fruit) g", "\\1 (g)") %>% 
               str_replace_all("(fat|sugar|salt) mg", "\\1 (mg)") %>% # Add (mg) 
               str_to_sentence() %>% # Apply sentence case to everything first 
               str_replace("\\(species\\)", "(Species)") # Fix the species wrapper 
      ) %>%
      mutate(variable = case_when(
        # GGB taxa should be uppercase 
        str_detect(variable, "\\(Species\\) ggb") ~ 
          str_replace(variable, "\\(Species\\) (.+)", function(x) {
            species_part <- str_extract(x, "(?<=\\(Species\\) ).*")
            paste0("(Species) ", str_to_upper(species_part))
          }),
        # Other species should be title case 
        str_detect(variable, "\\(Species\\)") ~ 
          str_replace(variable, "\\(Species\\) (.+)", function(x) {
            species_part <- str_extract(x, "(?<=\\(Species\\) ).*")
            paste0("(Species) ", str_to_title(species_part))
          }),
        # Everything else stays as is 
        TRUE ~ variable
      ))
    
    # -------------------------------------------------------------------------- #
    
    # *************************************** #
    ##   Enhanced Feature Importance Plot  ----
    # *************************************** #
    
    # Plot top variables 
    feat_imp_gut_microbiome_only <- ggplot(mean_importance, aes(x = reorder(variable, mean_importance), y = mean_importance)) +
      geom_col(fill = "#1F77B4", width = 0.7) +
      coord_flip() +
      labs(
        title = "Variable Importance – Gut microbiome data",
        x = NULL,
        y = "Importance (Gini decrease)"
      ) +
      theme_classic() +
      theme(
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.grid.minor = element_line(color = "gray95", size = 0.1),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    feat_imp_gut_microbiome_only
    
    # Summarise CV results
    cat("\n=== 5-Fold Cross-Validation Results ===\n")
    cat("Mean AUC:", round(mean(cv_results$auc), 3), "±", round(sd(cv_results$auc), 3), "\n")
    cat("AUC Range:", round(min(cv_results$auc), 3), "-", round(max(cv_results$auc), 3), "\n")
    cat("Mean Sensitivity:", round(mean(cv_results$sensitivity, na.rm = TRUE), 3), "±", round(sd(cv_results$sensitivity, na.rm = TRUE), 3), "\n")
    cat("Mean Specificity:", round(mean(cv_results$specificity, na.rm = TRUE), 3), "±", round(sd(cv_results$specificity, na.rm = TRUE), 3), "\n")
    cat("Mean Accuracy:", round(mean(cv_results$accuracy), 3), "±", round(sd(cv_results$accuracy), 3), "\n")
    
    
    # -------------------------------------------------------------------------- #
    
    # ************************* #
    ##   Enhanced AUC Plot   ----
    # ************************* #
    
    # Get ROC for each fold
    roc_list <- all_predictions %>%
      group_split(fold) %>%
      map(~ roc(response = .x$obs, predictor = .x$high, quiet = TRUE))
    
    # Interpolate sensitivity at common FPR points
    fpr_grid <- seq(0, 1, length.out = 100)
    
    roc_interp <- map_dfr(seq_along(roc_list), function(i) {
      r <- roc_list[[i]]
      sens <- r$sensitivities
      fpr <- 1 - r$specificities
      
      # Create a data frame to ensure alignment
      roc_df <- tibble(fpr = fpr, sens = sens) %>%
        filter(!is.na(fpr), !is.na(sens)) %>%
        distinct(fpr, .keep_all = TRUE) %>%
        arrange(fpr)
      
      # Ensure at least two points for interpolation
      if (nrow(roc_df) < 2) {
        return(tibble(fpr = fpr_grid, sens = NA_real_, fold = as.character(i)))
      }
      
      # Interpolate sensitivities at common FPR grid
      tibble(
        fpr = fpr_grid,
        sens = approx(roc_df$fpr, roc_df$sens, xout = fpr_grid, rule = 2)$y,
        fold = as.character(i)
      )
    }, .id = NULL) %>%
      group_by(fpr) %>%
      summarise(
        mean_sens = mean(sens, na.rm = TRUE),
        min_sens = min(sens, na.rm = TRUE),
        max_sens = max(sens, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Create AUC plot 
    mean_auc <- mean(cv_results$auc)
    range_auc <- range(cv_results$auc)
    
    p_gut_microbiome_only <- ggplot(roc_interp, aes(x = fpr)) +
      geom_ribbon(aes(ymin = min_sens, ymax = max_sens), fill = "gray80", alpha = 0.3) +
      geom_line(aes(y = mean_sens), color = "#1F77B4", size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
      labs(
        x = "False Positive Rate (1 - Specificity)",
        y = "Sensitivity",
        title = "ROC Curve – Gut microbiome data",
        subtitle = sprintf("5-fold Cross-Validation: Mean AUC: %.3f (Range: %.3f–%.3f)", mean_auc, range_auc[1], range_auc[2])
      ) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      theme_classic() +
      theme(
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.grid.minor = element_line(color = "gray95", size = 0.1),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      coord_fixed(ratio = 1)
    
    print(p_gut_microbiome_only)
    
    # ----------------------------------------------------------------------------------------------------------------------- #      
    
    # ************************************************** # 
    ##   Select data for only Covar data             ---- 
    # ************************************************** # 
    
    # Select gut microbiome data 
    dat <- matched_data %>% 
      select(record_id,
             high_risk,
             all_of(covars))
    
    names(dat)    
    
    # drop record id 
    dat <- dat %>%
      select(-record_id)
    
    # -------------------------------------------------------------------------- #  
    
    # ************************* # 
    ##   Prepare data        ---- 
    # ************************* # 
    
    # Rename and convert response to factor
    data <- dat %>%
      rename(response = high_risk) %>%
      mutate(response = case_when(response == TRUE ~ "high",
                                  .default = "low"),
             response = factor(response, levels = c("low", "high")))
    
    # Remove zero variance predictors 
    vars_no_outcome <- names(data)[names(data) != "response"] # Get the names except response 
    
    # Run nearZeroVar on subset without response
    nzv <- nearZeroVar(data %>% select(all_of(vars_no_outcome)))
    
    # get names of variables with near zero variance 
    nzv_names <- vars_no_outcome[nzv]
    nzv_names
    
    if (length(nzv_names) > 0) data <- data %>% select(-all_of(nzv_names))
    
    # -------------------------------------------------------------------------- #
    
    # ************************* #
    ##   5-fold CV  ----
    # ************************* #
    
    set.seed(1234)
    
    cat("dataset size:", nrow(data), "\n")
    cat("Class distribution:\n", table(data$response))
    
    cv_folds <- createFolds(data$response, k = 5, list = TRUE, returnTrain = FALSE)
    
    # Manual 5-fold CV 
    cv_results <- data.frame(
      fold = integer(),
      auc = numeric(),
      sensitivity = numeric(),
      specificity = numeric(),
      accuracy = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Store predictions for each fold
    all_predictions <- data.frame(
      obs = factor(),
      pred = factor(),
      high = numeric(),
      low = numeric(),
      fold = integer()
    )
    
    # Perform 5-fold CV 
    for (i in 1:5) { 
      
      cat("Processing fold", i, "of 5...\n")
      
      # Split data for this fold 
      test_fold <- data[cv_folds[[i]], ] 
      train_fold <- data[-cv_folds[[i]], ]
      
      # Train model on this fold 
      rf_model <- ranger(response ~ ., data = train_fold, num.trees = 500, mtry = floor(sqrt(ncol(train_fold) - 1)), 
                         importance = "impurity", probability = TRUE, seed = 1234 + i )
      
      # Extract variable importance
      imp <- rf_model$variable.importance
      imp_df <- data.frame(variable = names(imp), importance = as.numeric(imp), fold = i)
      
      # Store 
      var_importance_list[[i]] <- imp_df
      
      # Predict on test fold 
      fold_pred <- predict(rf_model, test_fold) 
      fold_probs <- fold_pred$predictions
      
      # Calculate ROC and AUC 
      roc_obj <- roc(test_fold$response, fold_probs[, "high"], quiet = TRUE) 
      fold_auc <- as.numeric(auc(roc_obj))
      
      # Get optimal threshold 
      best_thresh <- coords(roc_obj, "best", ret = "threshold")$threshold 
      fold_pred_class <- ifelse(fold_probs[, "high"] > best_thresh, "high", "low") 
      fold_pred_class <- factor(fold_pred_class, levels = c("low", "high"))
      
      # Calculate performance metrics 
      cm <- confusionMatrix(fold_pred_class, test_fold$response)
      
      # Store results 
      cv_results <- rbind(cv_results, data.frame(fold = i, auc = fold_auc, 
                                                 sensitivity = cm$byClass["Sensitivity"],
                                                 specificity = cm$byClass["Specificity"], 
                                                 accuracy = cm$overall["Accuracy"] ))
      
      # Store predictions for aggregated analysis 
      fold_preds <- data.frame(obs = test_fold$response,
                               pred = fold_pred_class, 
                               high = fold_probs[, "high"], 
                               low = fold_probs[, "low"], fold = i ) 
      
      all_predictions <- rbind(all_predictions, fold_preds)
      
    }
    
    # Combine all importance data frames
    all_importance <- bind_rows(var_importance_list)
    
    # Average importance across folds
    mean_importance <- all_importance %>%
      group_by(variable) %>%
      summarize(mean_importance = mean(importance), .groups = "drop") %>%
      arrange(desc(mean_importance)) %>%
      slice(1:20)  # Top 20 variables
    
    
    # -------------------------------------------------------------------------- #  
    
    # clean feature names 
    mean_importance = mean_importance %>%
      mutate(variable = variable %>%
               str_replace_all("_", " ") %>% # Replace underscores 
               str_remove_all("\\bimputed\\b") %>% # Remove the word "imputed" 
               str_squish() %>% # Remove extra whitespace 
               str_replace("^species (.+)", "(species) \\1") %>% # Wrap species names 
               str_replace_all("(\\b\\w+ meat|alcohol|nsp|fibre|fiber|carbohydrate|energy|red meat|protein|processed meat|vegetables|fruit) g", "\\1 (g)") %>% 
               str_replace_all("(fat|sugar|salt) mg", "\\1 (mg)") %>% # Add (mg) 
               str_to_sentence() %>% # Apply sentence case to everything first 
               str_replace("\\(species\\)", "(Species)") # Fix the species wrapper 
      ) %>%
      mutate(variable = case_when(
        # GGB taxa should be uppercase 
        str_detect(variable, "\\(Species\\) ggb") ~ 
          str_replace(variable, "\\(Species\\) (.+)", function(x) {
            species_part <- str_extract(x, "(?<=\\(Species\\) ).*")
            paste0("(Species) ", str_to_upper(species_part))
          }),
        # Other species should be title case 
        str_detect(variable, "\\(Species\\)") ~ 
          str_replace(variable, "\\(Species\\) (.+)", function(x) {
            species_part <- str_extract(x, "(?<=\\(Species\\) ).*")
            paste0("(Species) ", str_to_title(species_part))
          }),
        # Everything else stays as is 
        TRUE ~ variable
      ))
    
    # -------------------------------------------------------------------------- #
    
    # *************************************** #
    ##   Feature Importance Plot  ----
    # *************************************** #
    
    # Plot top variables 
    feat_imp_covar_data_only <- ggplot(mean_importance, aes(x = reorder(variable, mean_importance), y = mean_importance)) +
      geom_col(fill = "#1F77B4", width = 0.7) +
      coord_flip() +
      labs(
        title = "Variable Importance – Covariate data",
        x = NULL,
        y = "Importance (Gini decrease)"
      ) +
      theme_classic() +
      theme(
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.grid.minor = element_line(color = "gray95", size = 0.1),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    feat_imp_covar_data_only
    
    # Summarise CV results
    cat("\n=== 5-Fold Cross-Validation Results ===\n")
    cat("Mean AUC:", round(mean(cv_results$auc), 3), "±", round(sd(cv_results$auc), 3), "\n")
    cat("AUC Range:", round(min(cv_results$auc), 3), "-", round(max(cv_results$auc), 3), "\n")
    cat("Mean Sensitivity:", round(mean(cv_results$sensitivity, na.rm = TRUE), 3), "±", round(sd(cv_results$sensitivity, na.rm = TRUE), 3), "\n")
    cat("Mean Specificity:", round(mean(cv_results$specificity, na.rm = TRUE), 3), "±", round(sd(cv_results$specificity, na.rm = TRUE), 3), "\n")
    cat("Mean Accuracy:", round(mean(cv_results$accuracy), 3), "±", round(sd(cv_results$accuracy), 3), "\n")
    
    
    # -------------------------------------------------------------------------- #
    
    # ************************* #
    ##  AUC Plot   ----
    # ************************* #
    
    # Get ROC for each fold
    roc_list <- all_predictions %>%
      group_split(fold) %>%
      map(~ roc(response = .x$obs, predictor = .x$high, quiet = TRUE))
    
    # Interpolate sensitivity at common FPR points
    fpr_grid <- seq(0, 1, length.out = 100)
    
    roc_interp <- map_dfr(seq_along(roc_list), function(i) {
      r <- roc_list[[i]]
      sens <- r$sensitivities
      fpr <- 1 - r$specificities
      
      # Create a data frame to ensure alignment
      roc_df <- tibble(fpr = fpr, sens = sens) %>%
        filter(!is.na(fpr), !is.na(sens)) %>%
        distinct(fpr, .keep_all = TRUE) %>%
        arrange(fpr)
      
      # Ensure at least two points for interpolation
      if (nrow(roc_df) < 2) {
        return(tibble(fpr = fpr_grid, sens = NA_real_, fold = as.character(i)))
      }
      
      # Interpolate sensitivities at common FPR grid
      tibble(
        fpr = fpr_grid,
        sens = approx(roc_df$fpr, roc_df$sens, xout = fpr_grid, rule = 2)$y,
        fold = as.character(i)
      )
    }, .id = NULL) %>%
      group_by(fpr) %>%
      summarise(
        mean_sens = mean(sens, na.rm = TRUE),
        min_sens = min(sens, na.rm = TRUE),
        max_sens = max(sens, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Create enhanced AUC plot
    mean_auc <- mean(cv_results$auc)
    range_auc <- range(cv_results$auc)
    
    p_covar_data_only <- ggplot(roc_interp, aes(x = fpr)) +
      geom_ribbon(aes(ymin = min_sens, ymax = max_sens), fill = "gray80", alpha = 0.3) +
      geom_line(aes(y = mean_sens), color = "#1F77B4", size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
      labs(
        x = "False Positive Rate (1 - Specificity)",
        y = "Sensitivity",
        title = "ROC Curve – Covariate data",
        subtitle = sprintf("5-fold Cross-Validation: Mean AUC: %.3f (Range: %.3f–%.3f)", mean_auc, range_auc[1], range_auc[2])
      ) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      theme_classic() +
      theme(
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.grid.minor = element_line(color = "gray95", size = 0.1),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      coord_fixed(ratio = 1)
    
    print(p_covar_data_only)

    
    # ----------------------------------------------------------------------------------------------------------------------- #  
    
    # combine plots 
    
    library(patchwork)

    # Combine all plots
    wrap_plots(p_gut_microbiome_only, p_covar_data_only, p_gut_microbiome_and_covars, 
               feat_imp_gut_microbiome_only, feat_imp_covar_data_only, feat_imp_gut_microbiome_and_covars,
               ncol = 3, nrow = 2) +
      plot_annotation(
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)),
        tag_levels = 'A',     # Adds tags A, B, C, ...
        tag_prefix = "(",     # Makes them look like (A), (B), ...
        tag_suffix = ")"
      ) &
      theme(plot.tag.position = c(0.2, 1), plot.tag = element_text(size = 12, face = "bold"))    
    
    # save plot 
    path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_high_risk_prediction/P1_high_risk_prediction_feature_importance_sensitivity_analysis.png")

    png(
      path,
      width = 21,
      height = 12,
      units = 'in',
      res = 300
    )
    
    wrap_plots(p_gut_microbiome_only, p_covar_data_only, p_gut_microbiome_and_covars, 
               feat_imp_gut_microbiome_only, feat_imp_covar_data_only, feat_imp_gut_microbiome_and_covars,
               ncol = 3, nrow = 2) +
      plot_annotation(
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)),
        tag_levels = 'A',     # Adds tags A, B, C, ...
        tag_prefix = "(",     # Makes them look like (A), (B), ...
        tag_suffix = ")"
      ) &
      theme(plot.tag.position = c(0.2, 1), plot.tag = element_text(size = 12, face = "bold"))    
    
    dev.off()  
    