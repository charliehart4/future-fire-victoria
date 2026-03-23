# Setup 
library(PresenceAbsence)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

#### FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### True Skill Statistic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate TSS across thresholds for one fold

calc_tss_list <- function(df){
  # TSS Across a set of thresholds
  thresholds <- seq(0, 1, by = 0.01)
  calc_tss <- function(df, thresholds) {
    sapply(thresholds, function(thresh) {
      conf_mat <- try(cmx(DATA = df[, 1:3], threshold = thresh), silent = TRUE)
      if (inherits(conf_mat, "try-error")) return(NA)
      
      sens <- try(sensitivity(conf_mat)[1, 1], silent = TRUE)
      spec <- try(specificity(conf_mat)[1, 1], silent = TRUE)
      
      if (inherits(sens, "try-error") || inherits(spec, "try-error")) return(NA)
      
      tss = (sens + spec - 1)
      
      return(tss)
    })
  }
  
  # Apply to each fold and return a data frame
  tss_list <- map2(df, seq_along(df), ~{
    tss_vals <- calc_tss(.x, thresholds)
    tibble(fold = .y, threshold = thresholds, TSS = tss_vals)
  })
  
  tss_df <- bind_rows(tss_list)
  return(tss_df)
}


calc_tss<- function(df, model){
  df_tss = lapply(df, FUN = function(x){
    optimal_threshold <- optimal.thresholds(DATA = x[,1:3], threshold = 100)[2,2]
    conf_mat <- cmx(DATA = x[,1:3], threshold = optimal_threshold)
    tss <- sensitivity(conf_mat)[1,1] + specificity(conf_mat)[1,1] - 1
    tss
    return(data.frame(composite_fold_no = unique(x$composite_fold_no),
                      max_tss = tss))})
  
  df_tss = do.call("rbind", df_tss)
  df_tss$model <- model
  return(df_tss)
}

calc_brier <- function(x, model){
  brier_scores <- lapply(x, function(df) {
    brier = mean((df$pred - df$obs)^2)
    composite_fold_no = unique(df$composite_fold_no)
    return(data.frame(brier = brier, composite_fold_no=composite_fold_no))
  })
  
 df_brier = bind_rows(brier_scores)
 df_brier$model <- model
 return(df_brier)
}


#### Read in the data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cv_base_path <- here("datan", "cv_results")

read_cv_folder <- function(folder_name) {
  path <- file.path(cv_base_path, folder_name)
  files <- list.files(path, pattern = ".csv", full.names = TRUE)
  lapply(files, read.csv)
}

no_tsf_no_rac_data <- read_cv_folder("no_tsf_no_rac")
tsf_no_rac_data    <- read_cv_folder("tsf_no_rac")
no_tsf_rac_data    <- read_cv_folder("no_tsf_rac")
tsf_rac_data       <- read_cv_folder("tsf_rac")

# Update the master CV fold file
cv_data <- read.csv(here("data", "fire_modelling_rac_ready_forCV.csv"))
cv_folds = cv_data %>% group_by(composite_fold_no, composite_fold, temporal_fold, spatial_folds) %>% summarise(min_year = min(year), max_year = max(year))

#### Area Under the Receiver Operating Curve ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fold_numbers = sapply(no_tsf_no_rac_data, FUN = function(x){composite_fold_no = unique(x$composite_fold_no)})

##### No TSF No RAC
no_tsf_no_rac_auc = sapply(no_tsf_no_rac_data, FUN = function(x){pROC::roc(x$obs, x$pred)$auc[[1]]})
no_tsf_no_rac_auc = data.frame(composite_fold_no = fold_numbers, 
                             AUC = no_tsf_no_rac_auc, 
                             model = "no_tsf_no_rac")
##### TSF No RAC
tsf_no_rac_auc = sapply(tsf_no_rac_data, FUN = function(x){pROC::roc(x$obs, x$pred)$auc[[1]]})
tsf_no_rac_auc = data.frame(composite_fold_no = fold_numbers, 
                               AUC = tsf_no_rac_auc, 
                            model = "tsf_no_rac")
##### No TSF  RAC
no_tsf_rac_auc = sapply(no_tsf_rac_data, FUN = function(x){pROC::roc(x$obs, x$pred)$auc[[1]]})
no_tsf_rac_auc = data.frame(composite_fold_no = fold_numbers, 
                               AUC = no_tsf_rac_auc, 
                            model = "no_tsf_rac")
##### TSF  RAC
tsf_rac_auc = sapply(tsf_rac_data, FUN = function(x){pROC::roc(x$obs, x$pred)$auc[[1]]})
tsf_rac_auc = data.frame(composite_fold_no = fold_numbers, 
                               AUC = tsf_rac_auc, 
                         model = "tsf_rac")


auc_fold_statistics = rbind(no_tsf_no_rac_auc, tsf_no_rac_auc, no_tsf_rac_auc, tsf_rac_auc)


#### True Skill Statistic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
no_tsf_no_rac_tss = calc_tss(no_tsf_no_rac_data, model = "no_tsf_no_rac")
tsf_no_rac_tss = calc_tss(tsf_no_rac_data, model = "tsf_no_rac")
no_tsf_rac_tss = calc_tss(no_tsf_rac_data, model = "no_tsf_rac")
tsf_rac_tss = calc_tss(tsf_rac_data, model = "tsf_rac")

tss_fold_statistics = rbind(no_tsf_no_rac_tss, tsf_no_rac_tss, no_tsf_rac_tss, tsf_rac_tss)

#### BRIER SCORES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
no_tsf_no_rac_brier = calc_brier(no_tsf_no_rac_data, model = "no_tsf_no_rac")
tsf_no_rac_brier = calc_brier(tsf_no_rac_data, model = "tsf_no_rac")
no_tsf_rac_brier = calc_brier(no_tsf_rac_data, model = "no_tsf_rac")
tsf_rac_brier = calc_brier(tsf_rac_data, model = "tsf_rac")

brier_fold_statistics = rbind(no_tsf_no_rac_brier, tsf_no_rac_brier, no_tsf_rac_brier, tsf_rac_brier)



#### Compile Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compiled_statistics = left_join(auc_fold_statistics, cv_folds)
compiled_statistics = left_join(compiled_statistics, tss_fold_statistics)
compiled_statistics = left_join(compiled_statistics, brier_fold_statistics)

compiled_statistics$temporal_fold_lab = paste0(compiled_statistics$min_year, "_", compiled_statistics$max_year)

#### Plot summaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
auc.plot = ggplot(compiled_statistics) + geom_tile(aes(x = temporal_fold_lab, y = as.factor(spatial_folds), fill = AUC)) + scale_fill_viridis_c() + theme_minimal() + labs(x = "Temporal Fold", y = "Spatial Fold") + facet_wrap(~model, ncol=4)
tss.plot = ggplot(compiled_statistics) + geom_tile(aes(x = temporal_fold_lab, y = as.factor(spatial_folds), fill = max_tss)) + scale_fill_viridis_c() + theme_minimal() + labs(x = "Temporal Fold", y = "Spatial Fold") + facet_wrap(~model, ncol=4)
brier.plot = ggplot(compiled_statistics) + geom_tile(aes(x = temporal_fold_lab, y = as.factor(spatial_folds), fill = brier)) + scale_fill_viridis_c() + theme_minimal() + labs(x = "Temporal Fold", y = "Spatial Fold")  + facet_wrap(~model, ncol=4)

library(patchwork)
combined.plot = auc.plot/tss.plot/brier.plot + plot_layout(axes = "collect", guides = "collect")

ggsave(plot = combined.plot, 
       filename = here("outputs", "crossvalidation_results_summary.pdf"), 
       width = 6, height = 8, scale = 1.5)

#### Boxplots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boxplots = compiled_statistics %>% 
  pivot_longer(cols = c(AUC, max_tss, brier), names_to = "Variable", values_to = "Value") %>%
  ggplot() + geom_boxplot(aes(x = model, y = Value, fill=model)) + scale_fill_viridis_d(option="B") + 
  facet_wrap(~Variable, scales="free_y", ncol=1) + xlab("") + theme_bw() + theme(legend.position = "none")
boxplots





# Average TSS at a Threshold

mean(fold_statistics$max_tss) # Mean



# Combine into a single tidy data frame
no_rac_tss_df <- bind_rows(no_rac_tss_list)

# Compute average TSS across folds
mean_tss_df <- no_rac_tss_df %>%
  group_by(threshold) %>%
  summarise(mean_TSS = mean(TSS, na.rm = TRUE), .groups = "drop") %>%
  mutate(fold = "mean")

# Plot using ggplot
ggplot(no_rac_tss_df, aes(x = threshold, y = TSS, group = as.factor(fold))) +
  geom_line(alpha = 0.4, colour = "grey50") +
  geom_line(data = mean_tss_df, aes(x = threshold, y = mean_TSS,group = as.factor(fold)), color = "red", linewidth = 1.2) +
  labs(x = "Threshold", y = "True Skill Statistic") +
  theme_minimal()

#### BRIER SCORES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Assuming `data` is a list of data frames with 'obs' and 'pred' columns
brier_scores <- lapply(no_rac_data, function(df) {
  brier = mean((df$pred - df$obs)^2)
  composite_fold_no = unique(df$composite_fold_no)
  return(data.frame(brier = brier, composite_fold_no=composite_fold_no))
})

brier_scores = bind_rows(brier_scores)

fold_statistics = left_join(fold_statistics, brier_scores)


mean(brier_scores$brier)  # Average across folds

brier_df <- tibble(
  fold = brier_scores$composite_fold_no,
  brier_score = brier_scores$brier
)

# Plot
ggplot(brier_df, aes(x = fold, y = brier_score)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = mean(brier_scores), linetype = "dashed", color = "red") +
  labs(title = "Brier Score per Fold",
       y = "Brier Score",
       x = "Fold") +
  theme_minimal()

#### Compile Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fold_statistics = left_join(fold_statistics, cv_folds)
fold_statistics$Autocorrelation = "No RAC"

#### Read in the data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rac_data <- lapply(
  list.files(here("data", "cv_results", "tsf_rac"), 
             pattern = ".csv", 
             full.names = TRUE), 
  read.csv
)

#### Area Under the Receiver Operating Curve ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fold_numbers = sapply(rac_data, FUN = function(x){composite_fold_no = unique(x$composite_fold_no)})
rac_auc = sapply(rac_data, FUN = function(x){pROC::roc(x$obs, x$pred)$auc[[1]]})
fold_statistics_with_rac = data.frame(composite_fold_no = fold_numbers, 
                             AUC = rac_auc)
mean(fold_statistics_with_rac$AUC)



#### True Skill Statistic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate TSS across thresholds for one fold
calc_tss <- function(df, thresholds) {
  sapply(thresholds, function(thresh) {
    conf_mat <- try(cmx(DATA = df[, 1:3], threshold = thresh), silent = TRUE)
    if (inherits(conf_mat, "try-error")) return(NA)
    
    sens <- try(sensitivity(conf_mat)[1, 1], silent = TRUE)
    spec <- try(specificity(conf_mat)[1, 1], silent = TRUE)
    
    if (inherits(sens, "try-error") || inherits(spec, "try-error")) return(NA)
    
    tss = (sens + spec - 1)
    
    return(tss)
  })
}

# Average TSS at a Threshold
rac_tss = lapply(rac_data, FUN = function(x){
  optimal_threshold <- optimal.thresholds(DATA = x[,1:3], threshold = 100)[2,2]
  conf_mat <- cmx(DATA = x[,1:3], threshold = optimal_threshold)
  tss <- sensitivity(conf_mat)[1,1] + specificity(conf_mat)[1,1] - 1
  tss
  return(data.frame(composite_fold_no = unique(x$composite_fold_no),
                    max_tss = tss))})

rac_tss = do.call("rbind", rac_tss)

fold_statistics_with_rac = left_join(fold_statistics_with_rac, rac_tss)

mean(fold_statistics_with_rac$max_tss) # Mean

# TSS Across a set of thresholds
thresholds <- seq(0, 1, by = 0.01)

# Apply to each fold and return a data frame
rac_tss_list <- map2(rac_data, seq_along(rac_data), ~{
  tss_vals <- calc_tss(.x, thresholds)
  tibble(fold = .y, threshold = thresholds, TSS = tss_vals)
})

# Combine into a single tidy data frame
rac_tss_df <- bind_rows(rac_tss_list)

# Compute average TSS across folds
mean_tss_df <- rac_tss_df %>%
  group_by(threshold) %>%
  summarise(mean_TSS = mean(TSS, na.rm = TRUE), .groups = "drop") %>%
  mutate(fold = "mean")

# Plot using ggplot
ggplot(rac_tss_df, aes(x = threshold, y = TSS, group = as.factor(fold))) +
  geom_line(alpha = 0.4, colour = "grey50") +
  geom_line(data = mean_tss_df, aes(x = threshold, y = mean_TSS,group = as.factor(fold)), color = "red", linewidth = 1.2) +
  labs(x = "Threshold", y = "True Skill Statistic") +
  theme_minimal()

#### BRIER SCORES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Assuming `data` is a list of data frames with 'obs' and 'pred' columns
brier_scores <- lapply(rac_data, function(df) {
  brier = mean((df$pred - df$obs)^2)
  composite_fold_no = unique(df$composite_fold_no)
  return(data.frame(brier = brier, composite_fold_no=composite_fold_no))
})

brier_scores = bind_rows(brier_scores)

fold_statistics_with_rac = left_join(fold_statistics_with_rac, brier_scores)


mean(brier_scores$brier)  # Average across folds

brier_df <- tibble(
  fold = brier_scores$composite_fold_no,
  brier_score = brier_scores$brier
)

# Plot
ggplot(brier_df, aes(x = fold, y = brier_score)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = mean(brier_scores), linetype = "dashed", color = "red") +
  labs(title = "Brier Score per Fold",
       y = "Brier Score",
       x = "Fold") +
  theme_minimal()

