
library(dplyr)
library(ggplot2)
library(tidyr)
cv.files = list.files("F:/vic_fire_mapping/output_data/cross_validation/rolling_cv", full=T)
cvs = lapply(cv.files, read.csv)

sens_spec = function(df){
  thresholds_test <- seq(0, 1, by = 0.01)
  tss_values <- numeric(length(thresholds_test))
  for (i in seq_along(thresholds_test)) {
    pred <- ifelse(df$prediction >= thresholds_test[i], 1, 0)
    tp <- sum(df$burnt == 1 & pred == 1)
    tn <- sum(df$burnt == 0 & pred == 0)
    fp <- sum(df$burnt == 0 & pred == 1)
    fn <- sum(df$burnt == 1 & pred == 0)
    
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    tss_values[i] <- sensitivity + specificity - 1
  }
  
  # Find threshold that maximizes TSS
  optimal_idx <- which.max(tss_values)
  threshold <- thresholds_test[optimal_idx]
  max_tss <- tss_values[optimal_idx]
  
  cat(sprintf("Optimal threshold (max TSS): %.3f\n", threshold))
  cat(sprintf("Maximum TSS value: %.4f\n\n", max_tss))
  
  # Convert probabilities to binary predictions
  predicted_class <- ifelse(df$prediction >= threshold, 1, 0)
  
  # Create confusion matrix components
  actual_positive <- df$burnt == 1
  actual_negative <- df$burnt == 0
  predicted_positive <- predicted_class == 1
  predicted_negative <- predicted_class == 0
  
  # Calculate confusion matrix values
  TP <- sum(actual_positive & predicted_positive)  # True Positives
  TN <- sum(actual_negative & predicted_negative)  # True Negatives
  FP <- sum(actual_negative & predicted_positive)  # False Positives
  FN <- sum(actual_positive & predicted_negative)  # False Negatives
  
  # Calculate rates
  FPR <- FP / (FP + TN)  # False Positive Rate
  FNR <- FN / (FN + TP)  # False Negative Rate
  
  out = data.frame(FPR = FPR, 
                   FNR = FNR, 
                   burnt.points = sum(df$burnt),
                   year = unique(df$year))
  return(out)
}

sens.spec.out = lapply(cvs, sens_spec)
sens.spec.out = do.call("rbind", sens.spec.out)

sens.spec.out = sens.spec.out %>% pivot_longer(cols = c(FNR, FPR), names_to = "metric", values_to = "values")

ggplot(sens.spec.out) + 
  geom_line(aes(x = year, y = values, colour = metric))

ggplot(sens.spec.out, aes(x = burnt.points, y = values, colour = metric)) + geom_point() + geom_smooth(method="lm")

# Regional
library(sf)
library(terra)
random_sample = vect("F:/vic_fire_mapping/covariates/random_points.shp")
random_sample = random_sample %>% st_as_sf()

fd_path <- "F:/vic_fire_mapping/output_data/fire_risk_maps/Fire districts/LF_DISTRICT.shp"

regions_sf <- st_read(fd_path, quiet = TRUE) |> 
  st_make_valid() %>% st_transform(crs(random_sample))

random_sample_regions = st_join(random_sample, regions_sf, largest=TRUE)
random_sample_regions_df = cbind(st_drop_geometry(random_sample_regions), st_coordinates(random_sample_regions))
random_sample_regions_df = random_sample_regions_df %>% select(x="X",y="Y", REGIONNAME, DSTRCTNAME)

sens_spec_reg = function(df, regions){
  # Calculate False Positive Rate and False Negative Rate per region
  # Join regions to dataframe
  df$coords = paste0(df$x,"_",df$y)
  regions$coords = paste0(regions$x,"_",regions$y)
  
  df <- left_join(df, regions, by = "coords")
  
  # Find optimal threshold that maximizes True Skill Statistic (TSS) globally
  thresholds_test <- seq(0, 1, by = 0.01)
  tss_values <- numeric(length(thresholds_test))
  
  for (i in seq_along(thresholds_test)) {
    pred <- ifelse(df$prediction >= thresholds_test[i], 1, 0)
    tp <- sum(df$burnt == 1 & pred == 1)
    tn <- sum(df$burnt == 0 & pred == 0)
    fp <- sum(df$burnt == 0 & pred == 1)
    fn <- sum(df$burnt == 1 & pred == 0)
    
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    tss_values[i] <- sensitivity + specificity - 1
  }
  
  # Find threshold that maximizes TSS
  optimal_idx <- which.max(tss_values)
  threshold <- thresholds_test[optimal_idx]
  max_tss <- tss_values[optimal_idx]
  
  cat(sprintf("Optimal threshold (max TSS): %.3f\n", threshold))
  cat(sprintf("Maximum TSS value: %.4f\n\n", max_tss))
  
  # Calculate FPR and FNR for each region
  regions_list <- unique(df$REGIONNAME)
  region_results <- data.frame(
    Region = character(),
    N = integer(),
    TP = integer(),
    TN = integer(),
    FP = integer(),
    FN = integer(),
    FPR = numeric(),
    FNR = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    TSS = numeric(),
    AUC = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (region in regions_list) {
    df_region <- df[df$REGIONNAME == region, ]
    
    # Convert probabilities to binary predictions using optimal threshold
    predicted_class <- ifelse(df_region$prediction >= threshold, 1, 0)
    
    # Calculate confusion matrix values
    tp <- sum(df_region$burnt == 1 & predicted_class == 1)
    tn <- sum(df_region$burnt == 0 & predicted_class == 0)
    fp <- sum(df_region$burnt == 0 & predicted_class == 1)
    fn <- sum(df_region$burnt == 1 & predicted_class == 0)
    
    # Calculate rates
    fpr <- fp / (fp + tn)
    fnr <- fn / (fn + tp)
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    tss <- sensitivity + specificity - 1
    
    if (length(unique(df_region$burnt)) <2){
      auc_value = NA
    } else{
      roc_obj <- pROC::roc(df_region$burnt, df_region$prediction, quiet=TRUE)
      auc_value <- as.numeric(pROC::auc(roc_obj))
    }
    
    # Add to results
    region_results <- rbind(region_results, data.frame(
      Region = region,
      year = unique(df_region$year),
      N = nrow(df_region),
      TP = tp,
      TN = tn,
      FP = fp,
      FN = fn,
      FPR = fpr,
      FNR = fnr,
      Sensitivity = sensitivity,
      Specificity = specificity,
      TSS = tss,
      AUC = auc_value
    ))
  }
  
  cat("\n--- Results by Region ---\n")
  print(region_results, row.names = FALSE)
  
  # Convert probabilities to binary predictions
  predicted_class <- ifelse(df$prediction >= threshold, 1, 0)
  
  # Create confusion matrix components
  actual_positive <- df$burnt == 1
  actual_negative <- df$burnt == 0
  predicted_positive <- predicted_class == 1
  predicted_negative <- predicted_class == 0
  
  # Calculate confusion matrix values
  TP <- sum(actual_positive & predicted_positive)  # True Positives
  TN <- sum(actual_negative & predicted_negative)  # True Negatives
  FP <- sum(actual_negative & predicted_positive)  # False Positives
  FN <- sum(actual_positive & predicted_negative)  # False Negatives
  
  # Calculate rates
  FPR <- FP / (FP + TN)  # False Positive Rate
  FNR <- FN / (FN + TP)  # False Negative Rate
  
  # Display results
  cat("Confusion Matrix:\n")
  cat(sprintf("True Positives (TP): %d\n", TP))
  cat(sprintf("True Negatives (TN): %d\n", TN))
  cat(sprintf("False Positives (FP): %d\n", FP))
  cat(sprintf("False Negatives (FN): %d\n", FN))
  cat("\n")
  cat(sprintf("False Positive Rate (FPR): %.4f\n", FPR))
  cat(sprintf("False Negative Rate (FNR): %.4f\n", FNR))
  cat("\n")
  cat(sprintf("Sensitivity (True Positive Rate): %.4f\n", TP / (TP + FN)))
  cat(sprintf("Specificity (True Negative Rate): %.4f\n", TN / (TN + FP)))
  
  return(region_results)
  
}

sens.spec.out = lapply(cvs, sens_spec_reg, regions = random_sample_regions_df)

sens.spec.out = do.call("rbind", sens.spec.out)


ggplot(sens.spec.out) + geom_line(aes(x=year, y = FNR)) + facet_wrap(~Region)

ggplot(sens.spec.out) + geom_line(aes(x=year, y = AUC)) + facet_wrap(~Region)


#### Temporal Fold PRediction Plots ####
library(terra)
library(tidyterra)
library(ggplot2)
library(patchwork)
map.files = list.files("F:/vic_fire_mapping/output_data/temporal_fold_prediction/prediction", pattern=".tif", full=T)

temporal.maps = rast(map.files)

map.plot = function(ras){
  ggplot() + geom_spatraster(data=ras) +
    scale_fill_viridis_c(option = "B", na.value = "transparent",limits = c(0,0.5)) + 
    theme_void()
}

maps = lapply(temporal.maps, map.plot)

wrap_plots(maps) + plot_layout(guides="collect")
