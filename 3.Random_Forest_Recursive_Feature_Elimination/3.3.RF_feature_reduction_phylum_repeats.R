## Author: 
#  Panayiotis Louca 

## Purpose of script: 
#  

## Date Created: 
#  27 February 2024 ------------------------------


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

# tidymodels collection 
library(tidymodels) 

# machine learning 
library(caret)
library(randomForest) # random forest 

# parallelisation 
library(doParallel)
library(foreach)

# interpretation 
library(vip)

# other 
library(doBy) # sortBy 

# ------------------------------------------------------------------------------------------------------------------- # 

##   IMPORT DATA ---- 
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

sapply(dat[ ,outcomes_to_test], function(x) table(x, useNA = 'ifany'))

###   check classes ---- 

# gut microbiome species 
dat %>% # all gut microbiome vars are numeric 
  select(
    starts_with("phyla_"),
    UNCLASSIFIED) %>%
  sapply(class) %>%
  table() 

# pull out gm data 
dat_gm <- dat %>% 
  select(
    starts_with("phyla_"),
    UNCLASSIFIED) 

# covariates 
dat %>% 
  select(
    -all_of(outcomes_to_test),
    -starts_with("phyla_"),
    -starts_with("species_"),
    -starts_with("genus_"),
    -UNCLASSIFIED, 
    -starts_with("alpha_div_")) %>%
  names(.)
  
dat_covars <- dat %>% 
  select(
    -all_of(outcomes_to_test),
    -starts_with("phyla_"),
    -starts_with("species_"),
    -starts_with("genus_"),
    -UNCLASSIFIED, 
    -starts_with("alpha_div_"))
  
# remove unadjusted fibre intake and endoscopist 
dat_covars <- dat_covars %>% select(-nsp_g_imputed,
                                    -endoscopist_imputed,
                                    -ffq_returned)

dat_covars %>% # Some characters instead of factor 
  skim() 

# ------------------------------------------------------------------------------------------------------------------- # 

# outcomes 
dat_outcomes <- dat %>%
  select(all_of(outcomes_to_test) )

dat_outcomes %>% 
  glimpse()  


# ------------------------------------------------------------------------------------------------------------------- # 

# merge back together 

dat <- cbind(dat_outcomes, dat_covars, dat_gm)

# ------------------------------------------------------------------------------------------------------------------- # 

#                                                 #
# Run on outcome: Any neoplasia detected       ----
#                                                 #

outcome  <- outcomes_to_test[2]

df <- 
  dat %>%
  filter(!is.na(!!sym(outcome) ) ) %>%
  select(record_id,
         outcome, 
         all_of(names(dat_covars) ),
         all_of(names(dat_gm) ) )

table(df[ ,outcome], useNA = 'ifany')

data <-
  df %>%
  select(-record_id) %>%
  dplyr::rename(response = outcome) # rename response to 'response' 

data <- as.data.frame(data)

    # -------------------------------------------------------------------------- #  

# subset taxa to only 5% most prevalent 
phylum_names <- grep("phyla_", names(data), value = TRUE)
presence_counts <- colSums(data[, phylum_names] > 0)
cols_to_keep <- names(presence_counts)[presence_counts >= (0.05 * nrow(data))]

cat("Total phylum columns:", length(phylum_names), "\n")
cat("Phyla present in ≥5% of samples:", length(cols_to_keep), "\n")

data <- data %>%
  select(-all_of(phylum_names)) %>%
  bind_cols(select(data, all_of(cols_to_keep)))

# Final NA and zero variance check
data <- data[complete.cases(data), ]
zero_vars <- nearZeroVar(data[, -which(names(data) == "response")])
if (length(zero_vars) > 0) data <- data[, -zero_vars]
data$response <- as.factor(data$response)

# --- Parallel setup ---
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Settings --- 
set.seed(1234)
n_repeats <- 100
seeds <- sample(1:10000, n_repeats)

# --- Results container ---
results_list <- list()

# --- Run RFE and evaluate on 100 seeds ---
for (i in seq_len(n_repeats)) {
  cat("Running seed", i, "/", n_repeats, "\n")
  set.seed(seeds[i])
  
  # Split train/test
  train_index <- createDataPartition(data$response, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test  <- data[-train_index, ]
  
  # RFE control
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verbose = FALSE,
    allowParallel = TRUE
  )
  
  # Run RFE
  rfe_fit <- rfe(
    x = train[, setdiff(names(train), "response")],
    y = train$response,
    sizes = c(5, 10, 20, 30, 50, 75, 100),
    rfeControl = ctrl,
    tuneLength = 5
  )
  
  selected_features <- predictors(rfe_fit)
  
  # Train model on training set using selected features
  rf_final <- randomForest(
    response ~ ., data = train[, c("response", selected_features)],
    ntree = 1000
  )
  
  # Predict on test set
  test_pred <- predict(rf_final, test[, selected_features])
  
  # Confusion matrix + metrics
  cm <- confusionMatrix(test_pred, test$response)
  acc <- cm$overall["Accuracy"]
  ci <- cm$overall[c("AccuracyLower", "AccuracyUpper")]
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Store results
  results_list[[i]] <- list(
    seed = seeds[i],
    selected_features = selected_features,
    accuracy = acc,
    ci_lower = ci["AccuracyLower"],
    ci_upper = ci["AccuracyUpper"],
    sensitivity = sens,
    specificity = spec
  )
}

# Stop parallel
stopCluster(cl)

# --- Format output ---
results_df <- map_dfr(results_list, ~{
  tibble(
    seed = .x$seed,
    accuracy = .x$accuracy,
    ci_lower = .x$ci_lower,
    ci_upper = .x$ci_upper,
    sensitivity = .x$sensitivity,
    specificity = .x$specificity,
    features = paste(.x$selected_features, collapse = ",")
  )
})

# View summary
head(results_df)

# pull out top features 
feature_lists <- strsplit(results_df$features, ",") # Split comma-separated features into a list-column

# Unlist all features from all seeds 
all_features <- unlist(feature_lists)

# Count frequency of each feature
feature_summary_df <- as.data.frame(table(all_features))
colnames(feature_summary_df) <- c("feature", "times_selected")

# Add percentage selection 
n_models <- nrow(results_df)
feature_summary_df <- feature_summary_df %>%
  mutate(selection_percent = (times_selected / n_models) * 100) %>%
  arrange(desc(times_selected))

# View top selected features
print(head(feature_summary_df, 20))

# list results 
final_results <- list(results_df,
                      feature_summary_df)

# save results 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_covariate_feature_reduction/phyla_imp_features/Any_neoplasia_species_RF_results.xlsx")
writexl::write_xlsx(final_results,
                    path)

    # -------------------------------------------------------------------------- #  
    # -------------------------------------------------------------------------- #  

# ------------------------------------------------------------------------------------------------------------------- # 

#                                                 #
# Run on outcome: Cancer detected              ----
#                                                 #

outcome  <- outcomes_to_test[3]

df <- 
  dat %>%
  filter(!is.na(!!sym(outcome) ) ) %>%
  select(record_id,
         outcome, 
         all_of(names(dat_covars) ),
         all_of(names(dat_gm) ) )

table(df[ ,outcome], useNA = 'ifany')

data <-
  df %>%
  select(-record_id) %>%
  dplyr::rename(response = outcome) # rename response to 'response' 

data <- as.data.frame(data)

# -------------------------------------------------------------------------- #  

# subset taxa to only 5% most prevalent 
phylum_names <- grep("phyla_", names(data), value = TRUE)
presence_counts <- colSums(data[, phylum_names] > 0)
cols_to_keep <- names(presence_counts)[presence_counts >= (0.05 * nrow(data))]

cat("Total genus columns:", length(phylum_names), "\n")
cat("Genera present in ≥5% of samples:", length(cols_to_keep), "\n")

data <- data %>%
  select(-all_of(phylum_names)) %>%
  bind_cols(select(data, all_of(cols_to_keep)))

# Final NA and zero variance check
data <- data[complete.cases(data), ]
zero_vars <- nearZeroVar(data[, -which(names(data) == "response")])
if (length(zero_vars) > 0) data <- data[, -zero_vars]
data$response <- as.factor(data$response)

# --- Parallel setup ---
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Settings --- 
set.seed(1234)
n_repeats <- 100
seeds <- sample(1:10000, n_repeats)

# --- Results container ---
results_list <- list()

# --- Run RFE and evaluate on 100 seeds ---
for (i in seq_len(n_repeats)) {
  cat("Running seed", i, "/", n_repeats, "\n")
  set.seed(seeds[i])
  
  # Split train/test
  train_index <- createDataPartition(data$response, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test  <- data[-train_index, ]
  
  # RFE control
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verbose = FALSE,
    allowParallel = TRUE
  )
  
  # Run RFE
  rfe_fit <- rfe(
    x = train[, setdiff(names(train), "response")],
    y = train$response,
    sizes = c(5, 10, 20, 30, 50, 75, 100),
    rfeControl = ctrl,
    tuneLength = 5
  )
  
  selected_features <- predictors(rfe_fit)
  
  # Train model on training set using selected features
  rf_final <- randomForest(
    response ~ ., data = train[, c("response", selected_features)],
    ntree = 1000
  )
  
  # Predict on test set
  test_pred <- predict(rf_final, test[, selected_features])
  
  # Confusion matrix + metrics
  cm <- confusionMatrix(test_pred, test$response)
  acc <- cm$overall["Accuracy"]
  ci <- cm$overall[c("AccuracyLower", "AccuracyUpper")]
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Store results
  results_list[[i]] <- list(
    seed = seeds[i],
    selected_features = selected_features,
    accuracy = acc,
    ci_lower = ci["AccuracyLower"],
    ci_upper = ci["AccuracyUpper"],
    sensitivity = sens,
    specificity = spec
  )
}

# Stop parallel
stopCluster(cl)

# --- Format output ---
results_df <- map_dfr(results_list, ~{
  tibble(
    seed = .x$seed,
    accuracy = .x$accuracy,
    ci_lower = .x$ci_lower,
    ci_upper = .x$ci_upper,
    sensitivity = .x$sensitivity,
    specificity = .x$specificity,
    features = paste(.x$selected_features, collapse = ",")
  )
})

# View summary
head(results_df)

# pull out top features 
feature_lists <- strsplit(results_df$features, ",") # Split comma-separated features into a list-column

# Unlist all features from all seeds 
all_features <- unlist(feature_lists)

# Count frequency of each feature
feature_summary_df <- as.data.frame(table(all_features))
colnames(feature_summary_df) <- c("feature", "times_selected")

# Add percentage selection 
n_models <- nrow(results_df)
feature_summary_df <- feature_summary_df %>%
  mutate(selection_percent = (times_selected / n_models) * 100) %>%
  arrange(desc(times_selected))

# View top selected features
print(head(feature_summary_df, 20))

# list results 
final_results <- list(results_df,
                      feature_summary_df)

# save results 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_covariate_feature_reduction/phyla_imp_features/Cancer_detected_species_RF_results.xlsx")
writexl::write_xlsx(final_results,
                    path)

# -------------------------------------------------------------------------- #  
# -------------------------------------------------------------------------- #  

# ------------------------------------------------------------------------------------------------------------------- # 

#                                                              #
# Run on outcome: Advanced neoplasia (exc cancer) detected  ----
#                                                              #

outcome  <- outcomes_to_test[4]

df <- 
  dat %>%
  filter(!is.na(!!sym(outcome) ) ) %>%
  select(record_id,
         outcome, 
         all_of(names(dat_covars) ),
         all_of(names(dat_gm) ) )

table(df[ ,outcome], useNA = 'ifany')

data <-
  df %>%
  select(-record_id) %>%
  dplyr::rename(response = outcome) # rename response to 'response' 

data <- as.data.frame(data)

# -------------------------------------------------------------------------- #  

# subset taxa to only 5% most prevalent 
phylum_names <- grep("phyla_", names(data), value = TRUE)
presence_counts <- colSums(data[, phylum_names] > 0)
cols_to_keep <- names(presence_counts)[presence_counts >= (0.05 * nrow(data))]

cat("Total phyla columns:", length(phylum_names), "\n")
cat("Phylum present in ≥5% of samples:", length(cols_to_keep), "\n")

data <- data %>%
  select(-all_of(phylum_names)) %>%
  bind_cols(select(data, all_of(cols_to_keep)))

# Final NA and zero variance check
data <- data[complete.cases(data), ]
zero_vars <- nearZeroVar(data[, -which(names(data) == "response")])
if (length(zero_vars) > 0) data <- data[, -zero_vars]
data$response <- as.factor(data$response)

# --- Parallel setup ---
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Settings --- 
set.seed(1234)
n_repeats <- 100
seeds <- sample(1:10000, n_repeats)

# --- Results container ---
results_list <- list()

# --- Run RFE and evaluate on 100 seeds ---
for (i in seq_len(n_repeats)) {
  cat("Running seed", i, "/", n_repeats, "\n")
  set.seed(seeds[i])
  
  # Split train/test
  train_index <- createDataPartition(data$response, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test  <- data[-train_index, ]
  
  # RFE control
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verbose = FALSE,
    allowParallel = TRUE
  )
  
  # Run RFE
  rfe_fit <- rfe(
    x = train[, setdiff(names(train), "response")],
    y = train$response,
    sizes = c(5, 10, 20, 30, 50, 75, 100),
    rfeControl = ctrl,
    tuneLength = 5
  )
  
  selected_features <- predictors(rfe_fit)
  
  # Train model on training set using selected features
  rf_final <- randomForest(
    response ~ ., data = train[, c("response", selected_features)],
    ntree = 1000
  )
  
  # Predict on test set
  test_pred <- predict(rf_final, test[, selected_features])
  
  # Confusion matrix + metrics
  cm <- confusionMatrix(test_pred, test$response)
  acc <- cm$overall["Accuracy"]
  ci <- cm$overall[c("AccuracyLower", "AccuracyUpper")]
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Store results
  results_list[[i]] <- list(
    seed = seeds[i],
    selected_features = selected_features,
    accuracy = acc,
    ci_lower = ci["AccuracyLower"],
    ci_upper = ci["AccuracyUpper"],
    sensitivity = sens,
    specificity = spec
  )
}

# Stop parallel
stopCluster(cl)

# --- Format output ---
results_df <- map_dfr(results_list, ~{
  tibble(
    seed = .x$seed,
    accuracy = .x$accuracy,
    ci_lower = .x$ci_lower,
    ci_upper = .x$ci_upper,
    sensitivity = .x$sensitivity,
    specificity = .x$specificity,
    features = paste(.x$selected_features, collapse = ",")
  )
})

# View summary
head(results_df)

# pull out top features 
feature_lists <- strsplit(results_df$features, ",") # Split comma-separated features into a list-column

# Unlist all features from all seeds 
all_features <- unlist(feature_lists)

# Count frequency of each feature
feature_summary_df <- as.data.frame(table(all_features))
colnames(feature_summary_df) <- c("feature", "times_selected")

# Add percentage selection 
n_models <- nrow(results_df)
feature_summary_df <- feature_summary_df %>%
  mutate(selection_percent = (times_selected / n_models) * 100) %>%
  arrange(desc(times_selected))

# View top selected features
print(head(feature_summary_df, 20))

# list results 
final_results <- list(results_df,
                      feature_summary_df)

# save results 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_covariate_feature_reduction/phyla_imp_features/Advanced_neoplasia_exc_cancer_species_RF_results.xlsx")
writexl::write_xlsx(final_results,
                    path)

# -------------------------------------------------------------------------- #  
# -------------------------------------------------------------------------- #  

# ------------------------------------------------------------------------------------------------------------------- # 

#                                                      #
# Run on outcome: None advanced neoplasia detected  ----
#                                                      #

outcome  <- outcomes_to_test[5]

df <- 
  dat %>%
  filter(!is.na(!!sym(outcome) ) ) %>%
  select(record_id,
         outcome, 
         all_of(names(dat_covars) ),
         all_of(names(dat_gm) ) )

table(df[ ,outcome], useNA = 'ifany')

data <-
  df %>%
  select(-record_id) %>%
  dplyr::rename(response = outcome) # rename response to 'response' 

data <- as.data.frame(data)

# -------------------------------------------------------------------------- #  

# subset taxa to only 5% most prevalent 
phylum_names <- grep("phyla_", names(data), value = TRUE)
presence_counts <- colSums(data[, phylum_names] > 0)
cols_to_keep <- names(presence_counts)[presence_counts >= (0.05 * nrow(data))]

cat("Total phyla columns:", length(phylum_names), "\n")
cat("Phylum present in ≥5% of samples:", length(cols_to_keep), "\n")

data <- data %>%
  select(-all_of(phylum_names)) %>%
  bind_cols(select(data, all_of(cols_to_keep)))

# Final NA and zero variance check
data <- data[complete.cases(data), ]
zero_vars <- nearZeroVar(data[, -which(names(data) == "response")])
if (length(zero_vars) > 0) data <- data[, -zero_vars]
data$response <- as.factor(data$response)

# --- Parallel setup ---
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Settings --- 
set.seed(1234)
n_repeats <- 100
seeds <- sample(1:10000, n_repeats)

# --- Results container ---
results_list <- list()

# --- Run RFE and evaluate on 100 seeds ---
for (i in seq_len(n_repeats)) {
  cat("Running seed", i, "/", n_repeats, "\n")
  set.seed(seeds[i])
  
  # Split train/test
  train_index <- createDataPartition(data$response, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test  <- data[-train_index, ]
  
  # RFE control
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verbose = FALSE,
    allowParallel = TRUE
  )
  
  # Run RFE
  rfe_fit <- rfe(
    x = train[, setdiff(names(train), "response")],
    y = train$response,
    sizes = c(5, 10, 20, 30, 50, 75, 100),
    rfeControl = ctrl,
    tuneLength = 5
  )
  
  selected_features <- predictors(rfe_fit)
  
  # Train model on training set using selected features
  rf_final <- randomForest(
    response ~ ., data = train[, c("response", selected_features)],
    ntree = 1000
  )
  
  # Predict on test set
  test_pred <- predict(rf_final, test[, selected_features])
  
  # Confusion matrix + metrics
  cm <- confusionMatrix(test_pred, test$response)
  acc <- cm$overall["Accuracy"]
  ci <- cm$overall[c("AccuracyLower", "AccuracyUpper")]
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Store results
  results_list[[i]] <- list(
    seed = seeds[i],
    selected_features = selected_features,
    accuracy = acc,
    ci_lower = ci["AccuracyLower"],
    ci_upper = ci["AccuracyUpper"],
    sensitivity = sens,
    specificity = spec
  )
}

# Stop parallel
stopCluster(cl)

# --- Format output ---
results_df <- map_dfr(results_list, ~{
  tibble(
    seed = .x$seed,
    accuracy = .x$accuracy,
    ci_lower = .x$ci_lower,
    ci_upper = .x$ci_upper,
    sensitivity = .x$sensitivity,
    specificity = .x$specificity,
    features = paste(.x$selected_features, collapse = ",")
  )
})

# View summary
head(results_df)

# pull out top features 
feature_lists <- strsplit(results_df$features, ",") # Split comma-separated features into a list-column

# Unlist all features from all seeds 
all_features <- unlist(feature_lists)

# Count frequency of each feature
feature_summary_df <- as.data.frame(table(all_features))
colnames(feature_summary_df) <- c("feature", "times_selected")

# Add percentage selection 
n_models <- nrow(results_df)
feature_summary_df <- feature_summary_df %>%
  mutate(selection_percent = (times_selected / n_models) * 100) %>%
  arrange(desc(times_selected))

# View top selected features
print(head(feature_summary_df, 20))

# list results 
final_results <- list(results_df,
                      feature_summary_df)

# save results 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_covariate_feature_reduction/phyla_imp_features/Non_advanced_neoplasia_species_RF_results.xlsx")
writexl::write_xlsx(final_results,
                    path)

# -------------------------------------------------------------------------- #  
# -------------------------------------------------------------------------- #  


# ------------------------------------------------------------------------------------------------------------------- # 

#                                                      #
# Run on outcome: Only adenomatous polyps detected  ----
#                                                      #

outcome  <- outcomes_to_test[6]

df <- 
  dat %>%
  filter(!is.na(!!sym(outcome) ) ) %>%
  select(record_id,
         outcome, 
         all_of(names(dat_covars) ),
         all_of(names(dat_gm) ) )

table(df[ ,outcome], useNA = 'ifany')

data <-
  df %>%
  select(-record_id) %>%
  dplyr::rename(response = outcome) # rename response to 'response' 

data <- as.data.frame(data)

# -------------------------------------------------------------------------- #  

# subset taxa to only 5% most prevalent 
phylum_names <- grep("phyla_", names(data), value = TRUE)
presence_counts <- colSums(data[, phylum_names] > 0)
cols_to_keep <- names(presence_counts)[presence_counts >= (0.05 * nrow(data))]

cat("Total phylum columns:", length(phylum_names), "\n")
cat("Phylum present in ≥5% of samples:", length(cols_to_keep), "\n")

data <- data %>%
  select(-all_of(phylum_names)) %>%
  bind_cols(select(data, all_of(cols_to_keep)))

# Final NA and zero variance check
data <- data[complete.cases(data), ]
zero_vars <- nearZeroVar(data[, -which(names(data) == "response")])
if (length(zero_vars) > 0) data <- data[, -zero_vars]
data$response <- as.factor(data$response)

# --- Parallel setup ---
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Settings --- 
set.seed(1234)
n_repeats <- 100
seeds <- sample(1:10000, n_repeats)

# --- Results container ---
results_list <- list()

# --- Run RFE and evaluate on 100 seeds ---
for (i in seq_len(n_repeats)) {
  cat("Running seed", i, "/", n_repeats, "\n")
  set.seed(seeds[i])
  
  # Split train/test
  train_index <- createDataPartition(data$response, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test  <- data[-train_index, ]
  
  # RFE control
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verbose = FALSE,
    allowParallel = TRUE
  )
  
  # Run RFE
  rfe_fit <- rfe(
    x = train[, setdiff(names(train), "response")],
    y = train$response,
    sizes = c(5, 10, 20, 30, 50, 75, 100),
    rfeControl = ctrl,
    tuneLength = 5
  )
  
  selected_features <- predictors(rfe_fit)
  
  # Train model on training set using selected features
  rf_final <- randomForest(
    response ~ ., data = train[, c("response", selected_features)],
    ntree = 1000
  )
  
  # Predict on test set
  test_pred <- predict(rf_final, test[, selected_features])
  
  # Confusion matrix + metrics
  cm <- confusionMatrix(test_pred, test$response)
  acc <- cm$overall["Accuracy"]
  ci <- cm$overall[c("AccuracyLower", "AccuracyUpper")]
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Store results
  results_list[[i]] <- list(
    seed = seeds[i],
    selected_features = selected_features,
    accuracy = acc,
    ci_lower = ci["AccuracyLower"],
    ci_upper = ci["AccuracyUpper"],
    sensitivity = sens,
    specificity = spec
  )
}

# Stop parallel
stopCluster(cl)

# --- Format output ---
results_df <- map_dfr(results_list, ~{
  tibble(
    seed = .x$seed,
    accuracy = .x$accuracy,
    ci_lower = .x$ci_lower,
    ci_upper = .x$ci_upper,
    sensitivity = .x$sensitivity,
    specificity = .x$specificity,
    features = paste(.x$selected_features, collapse = ",")
  )
})

# View summary
head(results_df)

# pull out top features 
feature_lists <- strsplit(results_df$features, ",") # Split comma-separated features into a list-column

# Unlist all features from all seeds 
all_features <- unlist(feature_lists)

# Count frequency of each feature
feature_summary_df <- as.data.frame(table(all_features))
colnames(feature_summary_df) <- c("feature", "times_selected")

# Add percentage selection 
n_models <- nrow(results_df)
feature_summary_df <- feature_summary_df %>%
  mutate(selection_percent = (times_selected / n_models) * 100) %>%
  arrange(desc(times_selected))

# View top selected features
print(head(feature_summary_df, 20))

# list results 
final_results <- list(results_df,
                      feature_summary_df)

# save results 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_covariate_feature_reduction/phyla_imp_features/Only_adenoma_species_RF_results.xlsx")
writexl::write_xlsx(final_results,
                    path)

# -------------------------------------------------------------------------- #  
# -------------------------------------------------------------------------- #  

# ------------------------------------------------------------------------------------------------------------------- # 

#                                                    #
# Run on outcome: Only serrated lesions detected  ----
#                                                    #

outcome  <- outcomes_to_test[7]

df <- 
  dat %>%
  filter(!is.na(!!sym(outcome) ) ) %>%
  select(record_id,
         outcome, 
         all_of(names(dat_covars) ),
         all_of(names(dat_gm) ) )

table(df[ ,outcome], useNA = 'ifany')

data <-
  df %>%
  select(-record_id) %>%
  dplyr::rename(response = outcome) # rename response to 'response' 

data <- as.data.frame(data)

# -------------------------------------------------------------------------- #  

# subset taxa to only 5% most prevalent 
phylum_names <- grep("phyla_", names(data), value = TRUE)
presence_counts <- colSums(data[, phylum_names] > 0)
cols_to_keep <- names(presence_counts)[presence_counts >= (0.05 * nrow(data))]

cat("Total phylum columns:", length(phylum_names), "\n")
cat("Phylum present in ≥5% of samples:", length(cols_to_keep), "\n")

data <- data %>%
  select(-all_of(phylum_names)) %>%
  bind_cols(select(data, all_of(cols_to_keep)))

# Final NA and zero variance check
data <- data[complete.cases(data), ]
zero_vars <- nearZeroVar(data[, -which(names(data) == "response")])
if (length(zero_vars) > 0) data <- data[, -zero_vars]
data$response <- as.factor(data$response)

# --- Parallel setup ---
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Settings --- 
set.seed(1234)
n_repeats <- 100
seeds <- sample(1:10000, n_repeats)

# --- Results container ---
results_list <- list()

# --- Run RFE and evaluate on 100 seeds ---
for (i in seq_len(n_repeats)) {
  cat("Running seed", i, "/", n_repeats, "\n")
  set.seed(seeds[i])
  
  # Split train/test
  train_index <- createDataPartition(data$response, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test  <- data[-train_index, ]
  
  # RFE control
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verbose = FALSE,
    allowParallel = TRUE
  )
  
  # Run RFE
  rfe_fit <- rfe(
    x = train[, setdiff(names(train), "response")],
    y = train$response,
    sizes = c(5, 10, 20, 30, 50, 75, 100),
    rfeControl = ctrl,
    tuneLength = 5
  )
  
  selected_features <- predictors(rfe_fit)
  
  # Train model on training set using selected features
  rf_final <- randomForest(
    response ~ ., data = train[, c("response", selected_features)],
    ntree = 1000
  )
  
  # Predict on test set
  test_pred <- predict(rf_final, test[, selected_features])
  
  # Confusion matrix + metrics
  cm <- confusionMatrix(test_pred, test$response)
  acc <- cm$overall["Accuracy"]
  ci <- cm$overall[c("AccuracyLower", "AccuracyUpper")]
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Store results
  results_list[[i]] <- list(
    seed = seeds[i],
    selected_features = selected_features,
    accuracy = acc,
    ci_lower = ci["AccuracyLower"],
    ci_upper = ci["AccuracyUpper"],
    sensitivity = sens,
    specificity = spec
  )
}

# Stop parallel
stopCluster(cl)

# --- Format output ---
results_df <- map_dfr(results_list, ~{
  tibble(
    seed = .x$seed,
    accuracy = .x$accuracy,
    ci_lower = .x$ci_lower,
    ci_upper = .x$ci_upper,
    sensitivity = .x$sensitivity,
    specificity = .x$specificity,
    features = paste(.x$selected_features, collapse = ",")
  )
})

# View summary
head(results_df)

# pull out top features 
feature_lists <- strsplit(results_df$features, ",") # Split comma-separated features into a list-column

# Unlist all features from all seeds 
all_features <- unlist(feature_lists)

# Count frequency of each feature
feature_summary_df <- as.data.frame(table(all_features))
colnames(feature_summary_df) <- c("feature", "times_selected")

# Add percentage selection 
n_models <- nrow(results_df)
feature_summary_df <- feature_summary_df %>%
  mutate(selection_percent = (times_selected / n_models) * 100) %>%
  arrange(desc(times_selected))

# View top selected features
print(head(feature_summary_df, 20))

# list results 
final_results <- list(results_df,
                      feature_summary_df)

# save results 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_covariate_feature_reduction/phyla_imp_features/Only_serrated_lesions_species_RF_results.xlsx")
writexl::write_xlsx(final_results,
                    path)

# -------------------------------------------------------------------------- #  
# -------------------------------------------------------------------------- #  

# ------------------------------------------------------------------------------------------------------------------- # 
  
#                                                    #
# Run on outcome: Only SSL detected  ----
#                                                    #

outcome  <- outcomes_to_test[8]

df <- 
  dat %>%
  filter(!is.na(!!sym(outcome) ) ) %>%
  select(record_id,
         outcome, 
         all_of(names(dat_covars) ),
         all_of(names(dat_gm) ) )

table(df[ ,outcome], useNA = 'ifany')

data <-
  df %>%
  select(-record_id) %>%
  dplyr::rename(response = outcome) # rename response to 'response' 

data <- as.data.frame(data)

# -------------------------------------------------------------------------- #  

# subset taxa to only 5% most prevalent 
phylum_names <- grep("phyla_", names(data), value = TRUE)
presence_counts <- colSums(data[, phylum_names] > 0)
cols_to_keep <- names(presence_counts)[presence_counts >= (0.05 * nrow(data))]

cat("Total Phylum columns:", length(phylum_names), "\n")
cat("Phyla present in ≥5% of samples:", length(cols_to_keep), "\n")

data <- data %>%
  select(-all_of(phylum_names)) %>%
  bind_cols(select(data, all_of(cols_to_keep)))

# Final NA and zero variance check
data <- data[complete.cases(data), ]
zero_vars <- nearZeroVar(data[, -which(names(data) == "response")])
if (length(zero_vars) > 0) data <- data[, -zero_vars]
data$response <- as.factor(data$response)

# --- Parallel setup ---
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Settings --- 
set.seed(1234)
n_repeats <- 100
seeds <- sample(1:10000, n_repeats)

# --- Results container ---
results_list <- list()

# --- Run RFE and evaluate on 100 seeds ---
for (i in seq_len(n_repeats)) {
  cat("Running seed", i, "/", n_repeats, "\n")
  set.seed(seeds[i])
  
  # Split train/test
  train_index <- createDataPartition(data$response, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test  <- data[-train_index, ]
  
  # RFE control
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verbose = FALSE,
    allowParallel = TRUE
  )
  
  # Run RFE
  rfe_fit <- rfe(
    x = train[, setdiff(names(train), "response")],
    y = train$response,
    sizes = c(5, 10, 20, 30, 50, 75, 100),
    rfeControl = ctrl,
    tuneLength = 5
  )
  
  selected_features <- predictors(rfe_fit)
  
  # Train model on training set using selected features
  rf_final <- randomForest(
    response ~ ., data = train[, c("response", selected_features)],
    ntree = 1000
  )
  
  # Predict on test set
  test_pred <- predict(rf_final, test[, selected_features])
  
  # Confusion matrix + metrics
  cm <- confusionMatrix(test_pred, test$response)
  acc <- cm$overall["Accuracy"]
  ci <- cm$overall[c("AccuracyLower", "AccuracyUpper")]
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Store results
  results_list[[i]] <- list(
    seed = seeds[i],
    selected_features = selected_features,
    accuracy = acc,
    ci_lower = ci["AccuracyLower"],
    ci_upper = ci["AccuracyUpper"],
    sensitivity = sens,
    specificity = spec
  )
}

# Stop parallel
stopCluster(cl)

# --- Format output ---
results_df <- map_dfr(results_list, ~{
  tibble(
    seed = .x$seed,
    accuracy = .x$accuracy,
    ci_lower = .x$ci_lower,
    ci_upper = .x$ci_upper,
    sensitivity = .x$sensitivity,
    specificity = .x$specificity,
    features = paste(.x$selected_features, collapse = ",")
  )
})

# View summary
head(results_df)

# pull out top features 
feature_lists <- strsplit(results_df$features, ",") # Split comma-separated features into a list-column

# Unlist all features from all seeds 
all_features <- unlist(feature_lists)

# Count frequency of each feature
feature_summary_df <- as.data.frame(table(all_features))
colnames(feature_summary_df) <- c("feature", "times_selected")

# Add percentage selection 
n_models <- nrow(results_df)
feature_summary_df <- feature_summary_df %>%
  mutate(selection_percent = (times_selected / n_models) * 100) %>%
  arrange(desc(times_selected))

# View top selected features
print(head(feature_summary_df, 20))

# list results 
final_results <- list(results_df,
                      feature_summary_df)

# save results 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_covariate_feature_reduction/phyla_imp_features/Only_ssl_species_RF_results.xlsx")
writexl::write_xlsx(final_results,
                    path)

# -------------------------------------------------------------------------- #  
# -------------------------------------------------------------------------- #  


# ------------------------------------------------------------------------------------------------------------------- # 

#                                                              #
# Run on outcome: Advanced neoplasia (inc cancer) detected  ----
#                                                              #

outcome  <- outcomes_to_test[11]

df <- 
  dat %>%
  filter(!is.na(!!sym(outcome) ) ) %>%
  select(record_id,
         outcome, 
         all_of(names(dat_covars) ),
         all_of(names(dat_gm) ) )

table(df[ ,outcome], useNA = 'ifany')

data <-
  df %>%
  select(-record_id) %>%
  dplyr::rename(response = outcome) # rename response to 'response' 

data <- as.data.frame(data)

# -------------------------------------------------------------------------- #  

# subset taxa to only 5% most prevalent 
phylum_names <- grep("phyla_", names(data), value = TRUE)
presence_counts <- colSums(data[, phylum_names] > 0)
cols_to_keep <- names(presence_counts)[presence_counts >= (0.05 * nrow(data))]

cat("Total phyla columns:", length(phylum_names), "\n")
cat("Phylum present in ≥5% of samples:", length(cols_to_keep), "\n")

data <- data %>%
  select(-all_of(phylum_names)) %>%
  bind_cols(select(data, all_of(cols_to_keep)))

# Final NA and zero variance check
data <- data[complete.cases(data), ]
zero_vars <- nearZeroVar(data[, -which(names(data) == "response")])
if (length(zero_vars) > 0) data <- data[, -zero_vars]
data$response <- as.factor(data$response)

# --- Parallel setup ---
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Settings --- 
set.seed(1234)
n_repeats <- 100
seeds <- sample(1:10000, n_repeats)

# --- Results container ---
results_list <- list()

# --- Run RFE and evaluate on 100 seeds ---
for (i in seq_len(n_repeats)) {
  cat("Running seed", i, "/", n_repeats, "\n")
  set.seed(seeds[i])
  
  # Split train/test
  train_index <- createDataPartition(data$response, p = 0.8, list = FALSE)
  train <- data[train_index, ]
  test  <- data[-train_index, ]
  
  # RFE control
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    verbose = FALSE,
    allowParallel = TRUE
  )
  
  # Run RFE
  rfe_fit <- rfe(
    x = train[, setdiff(names(train), "response")],
    y = train$response,
    sizes = c(5, 10, 20, 30, 50, 75, 100),
    rfeControl = ctrl,
    tuneLength = 5
  )
  
  selected_features <- predictors(rfe_fit)
  
  # Train model on training set using selected features
  rf_final <- randomForest(
    response ~ ., data = train[, c("response", selected_features)],
    ntree = 1000
  )
  
  # Predict on test set
  test_pred <- predict(rf_final, test[, selected_features])
  
  # Confusion matrix + metrics
  cm <- confusionMatrix(test_pred, test$response)
  acc <- cm$overall["Accuracy"]
  ci <- cm$overall[c("AccuracyLower", "AccuracyUpper")]
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Store results
  results_list[[i]] <- list(
    seed = seeds[i],
    selected_features = selected_features,
    accuracy = acc,
    ci_lower = ci["AccuracyLower"],
    ci_upper = ci["AccuracyUpper"],
    sensitivity = sens,
    specificity = spec
  )
}

# Stop parallel
stopCluster(cl)

# --- Format output ---
results_df <- map_dfr(results_list, ~{
  tibble(
    seed = .x$seed,
    accuracy = .x$accuracy,
    ci_lower = .x$ci_lower,
    ci_upper = .x$ci_upper,
    sensitivity = .x$sensitivity,
    specificity = .x$specificity,
    features = paste(.x$selected_features, collapse = ",")
  )
})

# View summary
head(results_df)

# pull out top features 
feature_lists <- strsplit(results_df$features, ",") # Split comma-separated features into a list-column

# Unlist all features from all seeds 
all_features <- unlist(feature_lists)

# Count frequency of each feature
feature_summary_df <- as.data.frame(table(all_features))
colnames(feature_summary_df) <- c("feature", "times_selected")

# Add percentage selection 
n_models <- nrow(results_df)
feature_summary_df <- feature_summary_df %>%
  mutate(selection_percent = (times_selected / n_models) * 100) %>%
  arrange(desc(times_selected))

# View top selected features
print(head(feature_summary_df, 20))

# list results 
final_results <- list(results_df,
                      feature_summary_df)

# save results 
path <- file.path("/Users/panayiotislouca/Documents/NCL_files/COLO-COHORT/Gut_microbiome_work/phase1_analysis/P1_covariate_feature_reduction/phyla_imp_features/Advanced_neoplasia_(inc_cancer)_species_RF_results.xlsx")
writexl::write_xlsx(final_results,
                    path)

# -------------------------------------------------------------------------- #  
# -------------------------------------------------------------------------- #  

