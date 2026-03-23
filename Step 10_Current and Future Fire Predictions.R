library(dplyr)
library(tidyr)
library(gbm)
library(terra)
library(readr)
library(tidyterra)
library(viridis)
library(scico)
library(patchwork)
library(tmap)
library(R.utils)

# Script 10: Generate Fire Predictions, Ensemble Mean, and Visualizations
# This script applies the trained BRT model to the 10 covariate stacks to generate Pr(Fire) 
# maps for observed, GCM baseline, and GCM future periods. It calculates and saves the 
# Ensemble Mean predictions and the Change Maps (Future - Baseline), and generates a 
# comprehensive, interactive tmap visualization.

# === CRITICAL MEMORY/DISK MANAGEMENT ===
options(scipen=999) # Prevents scientific notation in numbers
output_dir <- "./outputs"
output_stack_dir <- file.path(output_dir, "prediction_stacks")
output_pred_dir <- file.path(output_dir, "predictions")

if (!dir.exists(output_pred_dir)) dir.create(output_pred_dir, recursive = TRUE)

# ======================================================================
# GLOBAL SETUP & MODEL LOADING
# ======================================================================

# Define all climate models to process
CLIMATE_MODELS <- c(
  "ACCESS1-0 RCP85",
  "GFDL-ESM2M RCP85", 
  "CNRM-CM5 RCP85", 
  "MIROC5 RCP85"
)

# Load the trained BRT model
mod <- readRDS(file.path(output_dir, "brt_base_lr0.005_ffdi_95.RDS"))


# Helper function to convert the long model name (e.g., "ACCESS1-0 RCP85") 
# into a short, file-safe format used for loading GCM stacks (e.g., ACCESS1-0_RCP85)
# FIX: The original function was replacing '-' with '_', which caused the file loading failure.
# This version only replaces spaces and periods, preserving the hyphen.
clean_name_for_loading <- function(model) { gsub(" |\\.", "_", model) }

# Helper function to convert the long model name (e.g., "ACCESS1-0 RCP85") 
# into a safe format used for internal lists and prediction output (e.g., ACCESS1_0_RCP85).
# Retain the original function's name for internal consistency (Phases 3/4/5).
clean_name <- function(model) { gsub(" |-|\\.", "_", model) }


# ======================================================================
# FUNCTION: MAKE PREDICTION (REUSABLE)
# This function applies the BRT model to the long-term covariate stack
# to generate a single raster of the predicted fire probability (Pr(Fire))
# for the entire time period represented by the stack (e.g., 1985-2005).
# ======================================================================

make_prediction <- function(covariate.stack, brt.model){
  
  cat("	-> Making predictions...\n")
  
  # Ensure the stack has the expected number of layers (Confirmed 12 static + 6 dynamic = 18)
  if (terra::nlyr(covariate.stack) != 18) {
    warning(paste("Stack has", terra::nlyr(covariate.stack), "layers. Expected 18. Check input data."))
  }
  
  # Make the stack a data frame (xy=TRUE is critical for rasterization later)
  prediction.xy <- terra::as.data.frame(covariate.stack, xy=TRUE)
  
  # Ensure categorical variables are treated correctly 
  if ("fuel_management_zones" %in% names(prediction.xy)) {
    prediction.xy$fuel_management_zones <- as.factor(prediction.xy$fuel_management_zones)
  }  
  
  # Make the predictions
  predicted <- gbm::predict.gbm(brt.model, prediction.xy, type="response")
  
  prediction.xy$prediction = predicted
  predicted.xy = prediction.xy %>% dplyr::select(x, y, prediction)
  
  # Rasterise
  predicted.map <- terra::rast(predicted.xy, type="xyz", crs = crs(covariate.stack))
  
  rm(prediction.xy, predicted.xy); gc()
  return(predicted.map)
  }


# ======================================================================
# PHASE 1: EXPLICITLY LOAD ALL 10 PREDICTION STACKS (Verification Check)
# ======================================================================

cat("\n--- Phase 1: Loading All Prediction Stacks for Verification ---\n")

# --- Observed Baselines ---
cat("\n[OBSERVED BASELINES]\n")
tryCatch({
  # Match exact file name from Script 9c output
  stack_obs_hist <- rast(file.path(output_stack_dir, "baseline_1985_2005_observed_predstack.tif"))
  cat("Loaded stack_obs_hist. Layers:\n"); print(names(stack_obs_hist))
}, error = function(e) { stop(paste("Error loading Historical Observed stack:", e$message)) })

tryCatch({
  # Match exact file name from Script 9c output
  stack_obs_recent <- rast(file.path(output_stack_dir, "baseline_2006_2023_observed_predstack.tif"))
  cat("Loaded stack_obs_recent. Layers:\n"); print(names(stack_obs_recent))
}, error = function(e) { stop(paste("Error loading Recent Observed stack:", e$message)) })


# --- GCM Modelled Baselines and Futures ---
GCM_STACKS <- list()
for (model in CLIMATE_MODELS) {
  # Use the special cleaning function to match the filename structure
  model_clean_load <- clean_name_for_loading(model)
  # Use the original cleaning function for the internal list names
  model_clean_internal <- clean_name(model)
  cat(paste0("\n[", model, "]\n"))
  
  # Baseline Modelled (1985-2005)
  # FIX: Use 'model_clean_load' to match the file structure (hyphen retained)
  baseline_filename <- paste0("baseline_1985_2005_", model_clean_load, "_modelled_predstack.tif")
  tryCatch({
    stack_baseline <- rast(file.path(output_stack_dir, baseline_filename))
    GCM_STACKS[[paste0(model_clean_internal, "_baseline")]] <- stack_baseline
    cat(paste0("Loaded ", model, " Baseline Stack. Layers:\n")); print(names(stack_baseline))
  }, error = function(e) { warning(paste("Error loading", model, "Baseline stack:", e$message)) })
  
  # Future Modelled (2081-2099)
  # FIX: Use 'model_clean_load' to match the file structure (hyphen retained)
  future_filename <- paste0("future_2081_2099_", model_clean_load, "_modelled_predstack.tif")
  tryCatch({
    stack_future <- rast(file.path(output_stack_dir, future_filename))
    GCM_STACKS[[paste0(model_clean_internal, "_future")]] <- stack_future
    cat(paste0("Loaded ", model, " Future Stack. Layers:\n")); print(names(stack_future))
  }, error = function(e) { warning(paste("Error loading", model, "Future stack:", e$message)) })
}

# Remove the initial input stacks once they are loaded into GCM_STACKS and validated
rm(stack_obs_hist, stack_obs_recent); gc() 


# ======================================================================
# PHASE 1B: ALIGN AND VALIDATE VARIABLE ORDER IN ALL GCM STACKS
# ======================================================================

cat("\n--- Phase 1B: Aligning and validating variable names across all stacks ---\n")

expected_vars <- mod$var.names
cat("Model expects", length(expected_vars), "predictors:\n")
print(expected_vars)

for (nm in names(GCM_STACKS)) {
  stk <- GCM_STACKS[[nm]]
  cat("\nChecking:", nm, "\n")
  
  # Check for missing or extra layers
  missing_vars <- setdiff(expected_vars, names(stk))
  extra_vars   <- setdiff(names(stk), expected_vars)
  
  if (length(missing_vars) > 0) {
    warning(paste("⚠ Stack", nm, "is missing:", paste(missing_vars, collapse = ", ")))
  }
  if (length(extra_vars) > 0) {
    cat("Note: Stack", nm, "has extra layers that will be dropped:", paste(extra_vars, collapse = ", "), "\n")
  }
  
  # Subset and reorder to match model variable order
  present_vars <- intersect(expected_vars, names(stk))
  stk_aligned  <- stk[[present_vars]]
  stk_aligned  <- stk_aligned[[expected_vars]]  # ensures correct order (even if missing some)
  
  # Validate
  if (identical(names(stk_aligned), expected_vars)) {
    cat("✅ Alignment verified for", nm, "\n")
  } else {
    cat("⚠ Alignment incomplete for", nm, "\n")
  }
  
  # Replace original entry in list
  GCM_STACKS[[nm]] <- stk_aligned
}

cat("\nAll GCM stacks aligned to model variable order.\n")
gc()

mod$var.type
names(mod$var.names)


# ======================================================================
# PHASE 2: OBSERVED BASELINE PREDICTIONS
# ======================================================================

cat("\n--- Phase 2: Generating Observed Baseline Predictions ---\n")

# NOTE: Since the Observed Baseline stacks were removed above, we reload them here
# as temporary objects for prediction.
# Match exact file name from Script 9c output
stack_obs_hist <- rast(file.path(output_stack_dir, "baseline_1985_2005_observed_predstack.tif"))
# Match exact file name from Script 9c output
stack_obs_recent <- rast(file.path(output_stack_dir, "baseline_2006_2023_observed_predstack.tif"))


# 1. Historical Observed Baseline (1985-2005)
baseline_observed_fire_hist <- make_prediction(stack_obs_hist, mod)
writeRaster(baseline_observed_fire_hist, 
            file.path(output_pred_dir, "fire_observed_baseline_1985_2005.tif"), overwrite=TRUE)
cat("✅ Saved Historical Observed (1985-2005) prediction.\n")

# Cleanup
rm(baseline_observed_fire_hist); gc()

# 2. Recent Observed Baseline (2006-2023)
baseline_observed_fire_recent <- make_prediction(stack_obs_recent, mod)
writeRaster(baseline_observed_fire_recent, 
            file.path(output_pred_dir, "fire_observed_baseline_2006_2023.tif"), overwrite=TRUE)
cat("✅ Saved Recent Observed (2006-2023) prediction.\n")

# Cleanup
rm(baseline_observed_fire_recent, stack_obs_hist, stack_obs_recent); gc()


# ======================================================================
# PHASE 3: GCM BASELINE, FUTURE, AND CHANGE PREDICTIONS (Looping)
# ======================================================================

cat("\n--- Phase 3: Generating GCM Modelled Predictions and Change Maps ---\n")

for (model in CLIMATE_MODELS) {
  
  model_clean <- clean_name(model) # Used for output files and internal names
  cat(paste0("\n*** Processing Model: ", model, " ***\n"))
  
  # Access the pre-loaded stacks from the GCM_STACKS list (using the internal name)
  stack_baseline <- GCM_STACKS[[paste0(model_clean, "_baseline")]]
  stack_future <- GCM_STACKS[[paste0(model_clean, "_future")]]
  
  if (is.null(stack_baseline) || is.null(stack_future)) {
    cat(paste0("	> Skipping ", model, ": One or more prediction stacks failed to load.\n"))
    next
  }
  
  # Define Output File Names
  output_base_name <- file.path(output_pred_dir, paste0("fire_", model_clean, "_baseline_1985_2005.tif"))
  output_future_name <- file.path(output_pred_dir, paste0("fire_", model_clean, "_future_2081_2099.tif"))
  output_diff_name <- file.path(output_pred_dir, paste0("fire_change_", model_clean, "_2081_2099.tif"))
  
  # 1. Baseline Prediction
  baseline_fire_map <- make_prediction(stack_baseline, mod)
  writeRaster(baseline_fire_map, output_base_name, overwrite=TRUE)
  cat(paste0("✅ Saved Baseline fire prediction.\n"))
  
  # 2. Future Prediction
  future_fire_map <- make_prediction(stack_future, mod)
  writeRaster(future_fire_map, output_future_name, overwrite=TRUE)
  cat(paste0("✅ Saved Future fire prediction.\n"))
  
  # 3. Calculate and Save Change Map (Future - Baseline)
  change_map <- future_fire_map - baseline_fire_map
  names(change_map) <- paste0("fire_change_", model_clean)
  writeRaster(change_map, output_diff_name, overwrite=TRUE)
  cat(paste0("✅ Saved Change map (Future - Baseline).\n"))
  
  # Store the individual GCM predictions for Phase 4 ensemble calculation
  GCM_STACKS[[paste0(model_clean, "_pred_baseline")]] <- baseline_fire_map
  GCM_STACKS[[paste0(model_clean, "_pred_future")]] <- future_fire_map
  
  # Cleanup: The input stacks are no longer needed, and the temporary change map is saved.
  rm(stack_baseline, stack_future, change_map); gc()
}

# ======================================================================
# PHASE 4: ENSEMBLE PREDICTION (Change from Ensemble Mean Predictions)
# ======================================================================

cat("\n--- Phase 4: Generating Ensemble Mean Predictions and Change ---\n")

# Load individual GCM predictions from the GCM_STACKS list created in Phase 3
baseline_preds_list <- list()
future_preds_list <- list()

for (model in CLIMATE_MODELS) {
  model_clean <- clean_name(model)
  
  # Get baseline and future prediction maps
  base_pred <- GCM_STACKS[[paste0(model_clean, "_pred_baseline")]]
  fut_pred <- GCM_STACKS[[paste0(model_clean, "_pred_future")]]
  
  # Safety check (in case the loop was skipped in Phase 3)
  if (!is.null(base_pred) && !is.null(fut_pred)) {
    baseline_preds_list[[model_clean]] <- base_pred
    future_preds_list[[model_clean]] <- fut_pred
  } else {
    warning(paste("Missing individual prediction map for ensemble calculation:", model))
  }
}

# Cleanup: Clear the individual prediction rasters now that they are copied to lists
# and the GCM_STACKS list is no longer needed.
rm(GCM_STACKS); gc()

ensemble_future_map <- NULL
ensemble_change_map <- NULL

if (length(baseline_preds_list) == length(CLIMATE_MODELS)) {
  
  # 1. Calculate Ensemble Baseline Prediction (Mean across 4 GCM Baselines)
  stacked_baseline_preds <- rast(baseline_preds_list)
  ensemble_baseline_map <- mean(stacked_baseline_preds, na.rm = TRUE)
  names(ensemble_baseline_map) <- "fire_Ensemble_Mean_Baseline"
  output_ensemble_base_name <- file.path(output_pred_dir, "fire_Ensemble_Mean_Baseline_1985_2005.tif")
  writeRaster(ensemble_baseline_map, output_ensemble_base_name, overwrite = TRUE)
  cat("✅ Saved Ensemble Mean Baseline prediction.\n")
  
  # 2. Calculate Ensemble Future Prediction (Mean across 4 GCM Futures)
  stacked_future_preds <- rast(future_preds_list)
  ensemble_future_map <- mean(stacked_future_preds, na.rm = TRUE)
  names(ensemble_future_map) <- "fire_Ensemble_Mean_Future"
  output_ensemble_future_name <- file.path(output_pred_dir, "fire_Ensemble_Mean_Future_2081_2099.tif")
  writeRaster(ensemble_future_map, output_ensemble_future_name, overwrite = TRUE)
  cat("✅ Saved Ensemble Mean Future prediction.\n")
  
  # 3. Calculate Ensemble Change Map (Ensemble Future - Ensemble Baseline)
  ensemble_change_map <- ensemble_future_map - ensemble_baseline_map
  names(ensemble_change_map) <- "fire_change_Ensemble_Mean"
  output_ensemble_name <- file.path(output_pred_dir, "fire_change_Ensemble_Mean_2081_2099.tif")
  writeRaster(ensemble_change_map, output_ensemble_name, overwrite = TRUE)
  cat("✅ Saved Ensemble Mean Change map (Future Ensemble - Baseline Ensemble).\n")
  
  # Cleanup: Remove the large, temporary stacked rasters and the ensemble maps now they are saved.
  rm(stacked_baseline_preds, stacked_future_preds, baseline_preds_list, future_preds_list, 
     ensemble_baseline_map, ensemble_future_map); gc()
  
} else {
  warning("Not all 4 GCM prediction maps were found for ensemble calculation. Skipping ensemble map creation.")
}

gc()



# ======================================================================
# PHASE 5: SIMPLIFIED VISUALIZATION SCRIPT (Final Fix)
# ======================================================================
# =========================
# Multi-layer tmap exports
# =========================

# Libraries (same set as your example)
library(tmap)
library(viridis)
library(scico)
library(terra)

tmap_mode("view")

# Root output dir

# ----------------------------
# 1) PREDICTION MAPS (10 LAYERS)
# ----------------------------

# Files in the order you provided
pred_files <- c(
  "fire_ACCESS1_0_RCP85_baseline_1985_2005.tif",
  "fire_ACCESS1_0_RCP85_future_2081_2099.tif",
  "fire_GFDL_ESM2M_RCP85_baseline_1985_2005.tif",
  "fire_GFDL_ESM2M_RCP85_future_2081_2099.tif",
  "fire_CNRM_CM5_RCP85_baseline_1985_2005.tif",
  "fire_CNRM_CM5_RCP85_future_2081_2099.tif",
  "fire_MIROC5_RCP85_baseline_1985_2005.tif",
  "fire_MIROC5_RCP85_future_2081_2099.tif",
  "fire_Ensemble_Mean_Baseline_1985_2005.tif",
  "fire_Ensemble_Mean_Future_2081_2099.tif"
)

pred_labels <- c(
  "ACCESS1-0 Baseline",
  "ACCESS1-0 Future",
  "GFDL-ESM2M Baseline",
  "GFDL-ESM2M Future",
  "CNRM-CM5 Baseline",
  "CNRM-CM5 Future",
  "MIROC5 Baseline",
  "MIROC5 Future",
  "Ensemble Mean Baseline",
  "Ensemble Mean Future"
)

# Shared legend settings for predictions
pred_breaks  <- c(0, 0.01, 0.05, 0.1, 0.15, 0.25, 0.35, 0.5, 0.75, 1)
pred_palette <- viridis(10, option = "inferno")

# Build layered map with a single legend (shown on the first layer only)
pred_map <- NULL
for (i in seq_along(pred_files)) {
  fpath <- file.path(out_dir, pred_files[i])
  r     <- rast(fpath)
  
  layer_i <- tm_shape(r, name = pred_labels[i]) +
    tm_raster(
      palette = pred_palette,
      breaks  = pred_breaks,
      n       = 10,
      style   = "fixed",
      title   = "Pr(Fire)",
      legend.show = (i == 1)  # show legend only on the first layer
    ) +
    tm_layout(
      title = "Fire Probability (Predictions)",
      title.position = c("left", "top"),
      legend.outside = TRUE
    )
  
  pred_map <- if (is.null(pred_map)) layer_i else (pred_map + layer_i)
}

# Save
tmap_save(
  pred_map,
  filename = file.path(out_dir, "01_Predictions_All_Maps_Interactive.html")
)

# ------------------------
# 2) CHANGE MAPS (5 LAYERS)
# ------------------------

change_files <- c(
  "fire_change_ACCESS1_0_RCP85_2081_2099.tif",
  "fire_change_GFDL_ESM2M_RCP85_2081_2099.tif",
  "fire_change_CNRM_CM5_RCP85_2081_2099.tif",
  "fire_change_MIROC5_RCP85_2081_2099.tif",
  "fire_change_Ensemble_Mean_2081_2099.tif"
)

change_labels <- c(
  "ACCESS1-0 Change",
  "GFDL-ESM2M Change",
  "CNRM-CM5 Change",
  "MIROC5 Change",
  "Ensemble Mean Change"
)

# Shared legend settings for changes
chg_breaks  <- seq(-0.5, 0.5, by = 0.1)
chg_palette <- scico(11, palette = "vik")

change_map <- NULL
for (i in seq_along(change_files)) {
  fpath <- file.path(out_dir, change_files[i])
  r     <- rast(fpath)
  
  layer_i <- tm_shape(r, name = change_labels[i]) +
    tm_raster(
      palette = chg_palette,
      breaks  = chg_breaks,
      style   = "fixed",
      title   = "Change in Pr(Fire)",
      midpoint = 0,
      legend.show = (i == 1)  # single legend
    ) +
    tm_layout(
      title = "Change in Fire Probability (Future – Baseline)",
      title.position = c("left", "top"),
      legend.outside = TRUE
    )
  
  change_map <- if (is.null(change_map)) layer_i else (change_map + layer_i)
}

# Save
tmap_save(
  change_map,
  filename = file.path(out_dir, "02_Change_All_Maps_Interactive.html")
)









# --- Minimal reload to inspect predictors ---

library(terra)

# Path to your prediction stacks
stack_path <- file.path(output_stack_dir, "future_2081_2099_ACCESS1-0_RCP85_modelled_predstack.tif")

# Load the stack for that GCM
s <- rast(stack_path)

# Check what layers exist
names(s)

# Pick the key variables most likely to show artefacts
key_vars <- c("spei12_mean", "thunderstorm_days", "ffdi_95_days")

# Subset
stack_test <- s[[key_vars]]

# Plot to look for horizontal/vertical striping
plot(stack_test)

# Optional: also check grid alignment & resolution
res(s)
ext(s)





library(terra)
stack_path <- "/prediction_stacks/future_2081_2099_ACCESS1-0_RCP85_modelled_predstack.tif"
s <- rast(stack_path)
names(s)

# Numeric summaries for all predictors
summary(values(s, na.rm = TRUE))

# If you have training data (or min/max saved earlier), compare:
training_ranges <- data.frame(
  variable = mod$var.names,
  min_train = mod$data$min,  # if stored; otherwise skip
  max_train = mod$data$max
)

unique_fmz <- unique(values(s$fuel_management_zones))
unique_thunder <- unique(values(s$thunderstorm_days))

cat("Fuel Management Zones:\n"); print(unique_fmz)
cat("Thunderstorm Days (sample):\n"); print(head(sort(unique_thunder)))


mod$var.type
data.frame(variable = mod$var.names, type = mod$var.type)


levels(mod)
# or, if it doesn’t show, try:
attr(mod, "var.levels")

plot(s[["fuel_management_zones"]])
plot(s[["thunderstorm_days"]])







# ============================================================
# Diagnose future thunderstorm_days (MIROC5 RCP85)
# ============================================================

library(terra)

# Path to the future prediction stack for MIROC5
stack_path <- "/prediction_stacks/future_2081_2099_MIROC5_RCP85_modelled_predstack.tif"

# Load the stack
s_miroc5 <- rast(stack_path)

# Confirm layer names
names(s_miroc5)

# Extract just the thunderstorm_days layer
r_thunder <- s_miroc5[["thunderstorm_days"]]

# --- Quick numeric checks ---
summary(values(r_thunder, na.rm = TRUE))
cat("Unique sample values:\n")
print(head(sort(unique(values(r_thunder))), 20))

# --- Plot to check for banding or tile artefacts ---
plot(r_thunder, main = "Thunderstorm Days (Future 2081–2099, MIROC5 RCP85)")

# Optional: rescale for better contrast (helps reveal subtle bands)
r_thunder_rel <- (r_thunder - global(r_thunder, "min", na.rm = TRUE)[1]) /
  (global(r_thunder, "max", na.rm = TRUE)[1] - global(r_thunder, "min", na.rm = TRUE)[1])
plot(r_thunder_rel, main = "Thunderstorm Days (Normalised)")

# ------------------------------------------------------------
# FAST diagnostic using a random subset (~10,000 pixels)
# ------------------------------------------------------------
set.seed(42)
vals_sample <- sample(values(r_thunder, na.rm = TRUE), 10000)

summary(vals_sample)
cat("\nMin:", min(vals_sample), "\nMax:", max(vals_sample), "\n")

hist(vals_sample,
     breaks = 50,
     main = "Thunderstorm Days (Sampled, MIROC5 Future)",
     xlab = "Thunderstorm Days per Year")

# Quick visual check on a cropped extent (e.g., central Victoria)
r_thunder_crop <- crop(r_thunder, ext(2300000, 2400000, 5200000, 5300000))
plot(r_thunder_crop, main = "Thunderstorm Days (Cropped MIROC5 Future)")



# ------------------------------------------------------------
# Load and check thunderstorm_days from model input data
# ------------------------------------------------------------
INPUT_CSV <- "data/full_fire_covariates_no_folds.csv"

dat <- read.csv(INPUT_CSV)

# Basic summary
summary(dat$thunderstorm_days)

# Range
range(dat$thunderstorm_days, na.rm = TRUE)

# Quick distribution plot
hist(dat$thunderstorm_days,
     breaks = 50,
     col = "grey",
     main = "Thunderstorm Days in Model Training Data",
     xlab = "Thunderstorm Days per Year")

# Optionally overlay density curve
lines(density(dat$thunderstorm_days, na.rm = TRUE), col = "red", lwd = 2)



# ------------------------------------------------------------
# Filter thunderstorm raster to show only high-end values (>50)
# ------------------------------------------------------------
r_thunder_high <- r_thunder
r_thunder_high[r_thunder_high < 50] <- NA

# Optional: give a small buffer on the upper range for colour contrast
high_min <- 50
high_max <- global(r_thunder_high, "max", na.rm = TRUE)[1]

# Plot
plot(r_thunder_high,
     col = viridis::viridis(100, option = "plasma"),
     main = "Thunderstorm Days > 50 (Future 2081–2099, MIROC5 RCP85)",
     zlim = c(high_min, high_max))  # <-- use zlim, not range

# Add legend label
mtext(paste0("Values range: ", round(high_min,1), "–", round(high_max,1)), side = 1, line = 2)







# ======================================================================
# NEW SECTION: RAC MODEL PREDICTIONS
# ======================================================================

cat("\n==============================\n")
cat("   RAC MODEL PREDICTIONS\n")
cat("==============================\n\n")

# Load RAC model
rac_mod <- readRDS("outputs/brt_full_rac_lr0.005_lr0.005_ffdi_95.rds")

expected_vars_rac <- rac_mod$var.names
cat("RAC model expects predictors:\n")
print(expected_vars_rac)

# List to store RAC predictions for ensemble
RAC_BASE <- list()
RAC_FUT  <- list()

for (model in CLIMATE_MODELS) {
  
  model_clean <- clean_name(model)
  cat(paste0("\n*** RAC Processing Model: ", model, " ***\n"))
  
  # Load stacks again (baseline + future)
  model_clean_load <- clean_name_for_loading(model)
  
  stack_baseline <- rast(file.path(
    output_stack_dir,
    paste0("baseline_1985_2005_", model_clean_load, "_modelled_predstack.tif")
  ))
  
  stack_future <- rast(file.path(
    output_stack_dir,
    paste0("future_2081_2099_", model_clean_load, "_modelled_predstack.tif")
  ))
  
  # ---------------------------
  # ADD RAC LAYER (all zeros)
  # ---------------------------
  rac_layer_base  <- stack_baseline[[1]] * 0
  rac_layer_fut   <- stack_future[[1]]  * 0
  names(rac_layer_base) <- "rac"
  names(rac_layer_fut)  <- "rac"
  
  # Append to stacks
  stack_baseline <- c(stack_baseline, rac_layer_base)
  stack_future   <- c(stack_future,   rac_layer_fut)
  
  # --------------------------------------------
  # Reorder to match RAC model variable order
  # --------------------------------------------
  missing_vars <- setdiff(expected_vars_rac, names(stack_baseline))
  if (length(missing_vars) > 0) {
    stop(paste("Missing predictors:", paste(missing_vars, collapse=", ")))
  }
  
  stack_baseline <- stack_baseline[[expected_vars_rac]]
  stack_future   <- stack_future[[expected_vars_rac]]
  
  # ---------------------------
  # MAKE PREDICTIONS
  # ---------------------------
  out_base  <- make_prediction(stack_baseline, rac_mod)
  out_fut   <- make_prediction(stack_future, rac_mod)
  out_delta <- out_fut - out_base
  
  # ---------------------------
  # SAVE OUTPUTS
  # ---------------------------
  fn_base  <- file.path(output_pred_dir, paste0("fire_", model_clean, "_baseline_1985_2005_RAC.tif"))
  fn_fut   <- file.path(output_pred_dir, paste0("fire_", model_clean, "_future_2081_2099_RAC.tif"))
  fn_delta <- file.path(output_pred_dir, paste0("fire_change_", model_clean, "_2081_2099_RAC.tif"))
  
  writeRaster(out_base,  fn_base,  overwrite=TRUE)
  writeRaster(out_fut,   fn_fut,   overwrite=TRUE)
  writeRaster(out_delta, fn_delta, overwrite=TRUE)
  
  cat("✓ Saved RAC baseline, future and change maps\n")
  
  # Store for ensemble later
  RAC_BASE[[model_clean]] <- out_base
  RAC_FUT[[model_clean]]  <- out_fut
  
  rm(out_base, out_fut, out_delta, stack_baseline, stack_future); gc()
}

# ======================================================================
# RAC ENSEMBLE
# ======================================================================

cat("\n--- RAC ENSEMBLE ---\n")

stacked_base <- rast(RAC_BASE)
stacked_fut  <- rast(RAC_FUT)

ens_base_rac <- mean(stacked_base, na.rm=TRUE)
ens_fut_rac  <- mean(stacked_fut,  na.rm=TRUE)
ens_delta_rac <- ens_fut_rac - ens_base_rac

writeRaster(ens_base_rac,  file.path(output_pred_dir, "fire_Ensemble_Mean_Baseline_1985_2005_RAC.tif"), overwrite=TRUE)
writeRaster(ens_fut_rac,   file.path(output_pred_dir, "fire_Ensemble_Mean_Future_2081_2099_RAC.tif"), overwrite=TRUE)
writeRaster(ens_delta_rac, file.path(output_pred_dir, "fire_change_Ensemble_Mean_2081_2099_RAC.tif"), overwrite=TRUE)

cat("✓ RAC Ensemble predictions saved\n")



# ============================================================
# RAC OBSERVED BASELINE (1985–2005)
# ============================================================

cat("\n--- RAC Observed Baseline (1985–2005) ---\n")

# Load observed historical covariate stack
stack_obs_hist <- rast(
  file.path(output_stack_dir, "baseline_1985_2005_observed_predstack.tif")
)

# ------------------------------------------------------------
# Add RAC layer (all zeros)
# ------------------------------------------------------------
rac_layer_obs <- stack_obs_hist[[1]] * 0
names(rac_layer_obs) <- "rac"

stack_obs_hist_rac <- c(stack_obs_hist, rac_layer_obs)

# ------------------------------------------------------------
# Reorder to match RAC model predictor order
# ------------------------------------------------------------
missing_vars <- setdiff(rac_mod$var.names, names(stack_obs_hist_rac))
if (length(missing_vars) > 0) {
  stop(paste("Observed RAC stack missing predictors:",
             paste(missing_vars, collapse = ", ")))
}

stack_obs_hist_rac <- stack_obs_hist_rac[[rac_mod$var.names]]

# ------------------------------------------------------------
# Predict
# ------------------------------------------------------------
fire_observed_baseline_1985_2005_RAC <-
  make_prediction(stack_obs_hist_rac, rac_mod)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
writeRaster(
  fire_observed_baseline_1985_2005_RAC,
  file.path(output_pred_dir,
            "fire_observed_baseline_1985_2005_RAC.tif"),
  overwrite = TRUE
)

cat("✓ Saved RAC observed baseline (1985–2005)\n")

rm(stack_obs_hist, stack_obs_hist_rac,
   fire_observed_baseline_1985_2005_RAC); gc()


