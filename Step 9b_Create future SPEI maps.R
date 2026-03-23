
# Script 9b: Calculate SPEI for Bias-Corrected Future Water Balance
# This script loads the bias-corrected Water Balance stacks produced by Script 9a,
# combines them into a single time series, and calculates the 12-month (SPEI-12)
# and 24-month (SPEI-24) SPEI indices for the future projection period.

library(terra)
library(SPEI)

# === CRITICAL MEMORY/DISK MANAGEMENT (COPIED FROM 9a) ===
options(scipen=999) # Prevents scientific notation in numbers
# 1. Set the temp directory to the dedicated F: drive path.
terraOptions(tempdir = "F:/vic_fire_mapping/terra_temp")
# 2. Tell terra to use a HUGE amount of memory (1 Terabyte) before trying to load
# raster blocks into RAM. This effectively forces it to use the disk for all operations.
terraOptions(memmax = 1e+12) 
gc()
# =========================================================


# ======================================================================
# GLOBAL SETUP
# ======================================================================

# Define all climate models to process (MUST match the folder names on disk)
CLIMATE_MODELS <- c(
  "ACCESS1-0 RCP85",
  "GFDL-ESM2M RCP85", 
  "CNRM-CM5 RCP85", 
  "MIROC5 RCP85"
)

# Define global paths
future.climate = here("data", "future_predictions")

# Define the reference period used for standardization in SPEI.
# *** CRITICAL FIX: Set to 1980-2005 to match historical training data baseline. ***
SPEI_REF_START <- c(1980, 1) # Start year and month
SPEI_REF_END <- c(2005, 12) # End year and month

# ======================================================================
# FUNCTION: CALCULATE AND SAVE SPEEI
# ======================================================================

# Define a function to calculate SPEI for a specific scale
calculate_and_save_spei <- function(input_stack, scale, model_name_clean, output_dir_spei) {
  
  output_file_prefix <- paste0("spei_", scale, "m_", model_name_clean)
  
  # ----------------------------------------------------------------------
  # --- NEW ROBUST FILE EXISTENCE CHECK (checks all expected months) ---
  # ----------------------------------------------------------------------
  
  # 1. Determine the dates for all valid SPEI output layers
  all_input_dates <- time(input_stack)
  
  if (length(all_input_dates) <= scale) {
    cat(paste0("Skipping SPEI-", scale, "m: Input stack has only ", length(all_input_dates), " months, insufficient for scale ", scale, ".\n"))
    return(NULL)
  }
  
  # Valid output starts AFTER the accumulation period (scale months)
  valid_output_dates <- all_input_dates[(scale + 1):length(all_input_dates)]
  
  # 2. Construct the full list of expected file paths using standardized date format
  expected_files <- file.path(output_dir_spei, paste0(output_file_prefix, "_", format(valid_output_dates, "%Y-%m-%d"), ".tif"))
  
  # 3. Check for missing files
  missing_files <- expected_files[!file.exists(expected_files)]
  
  if (length(missing_files) == 0) {
    cat(paste0("Skipping SPEI-", scale, "m: All ", length(expected_files), 
               " monthly files (", format(valid_output_dates[1], "%Y-%m-%d"), 
               " to ", format(valid_output_dates[length(valid_output_dates)], "%Y-%m-%d"), 
               ") already exist. ✅\n"))
    return(NULL) # Exit the function if all files are present
  } else {
    # Proceed with calculation, informing the user of the gaps
    total_files_expected <- length(expected_files)
    cat(paste0("Calculating SPEI-", scale, "m... (Total expected: ", total_files_expected, 
               ", Missing: ", length(missing_files), ". Recalculating all layers, but only saving missing ones.)\n")) # Updated log
    cat(paste0(" Earliest expected but missing file: ", basename(missing_files[1]), "\n"))
  }
  # ----------------------------------------------------------------------
  # --- END ROBUST FILE EXISTENCE CHECK ---
  # ----------------------------------------------------------------------
  
  # Extract the actual start date of the input raster stack to define the time series start
  start_year <- as.numeric(format(time(input_stack)[[1]], "%Y"))
  start_month <- as.numeric(format(time(input_stack)[[1]], "%m"))
  ts_stack_start <- c(start_year, start_month)
  
  # SPEI Calculation (Must run on the full stack for correct time series fitting)
  spei_stack <- terra::app(input_stack, fun = function(x, scale, ref_start, ref_end, ts_stack_start) {
    
    # The ts object must start at the beginning of the combined stack
    ts_data <- ts(x, start = ts_stack_start, frequency = 12)
    
    res <- tryCatch({
      # Calculate SPEI using the consistent reference period (1980-2005)
      spei_out <- SPEI::spei(ts_data, scale = scale,
                             ref.start = ref_start, ref.end = ref_end,
                             na.rm = TRUE, verbose = FALSE)
      as.numeric(spei_out$fitted)
    }, error = function(e) {
      # On error, return NA for that pixel
      rep(NA, length(x))
    })
    return(res)
  }, scale = scale, ref_start = SPEI_REF_START, ref_end = SPEI_REF_END, ts_stack_start = ts_stack_start)
  
  # Assign time information back
  time(spei_stack) <- time(input_stack)
  
  # Save individual layers (crucial for prediction input)
  cat("Saving individual SPEI layers (overwriting missing files to ensure completion)...\n")
  for (i in 1:nlyr(spei_stack)) {
    t = time(spei_stack[[i]])
    
    # Only process layers starting after the accumulation period
    if (i > scale) {	
      # CRITICAL: Construct the filename using the standardized date format
      output_file <- file.path(output_dir_spei, paste0(output_file_prefix, "_", format(t, "%Y-%m-%d"), ".tif"))
      
      # --- NEW: SKIP EXISTING FILES TO SAVE TIME ---
      if (file.exists(output_file)) {
        next # Skip to the next layer iteration
      }
      # ---------------------------------------------
      
      writeRaster(spei_stack[[i]],
                  filename = output_file,
                  overwrite = TRUE) # overwrite=TRUE is critical here
    }
  }
  cat(paste0("SPEI-", scale, "m layers saved/checked for ", model_name_clean, ". \n"))
  if (length(missing_files) > 0) {
    cat(paste0("   -> Successfully filled ", length(missing_files), " missing layers.\n"))
  }
}

# ======================================================================
# EXECUTION LOOP
# ======================================================================

for (model in CLIMATE_MODELS) {
  
  cat(paste0("\n--- Starting SPEI Calculation for Model: ", model, " ---\n"))
  
  # Use the model folder name to create the clean name used in the file path
  model_name_clean <- gsub(" |-|\\.", "_", model)
  
  # Define directory paths
  output_dir_wb <- file.path(future.climate, model, "wb_output")
  output_dir_spei <- file.path(future.climate, model, "spei_output")
  if (!dir.exists(output_dir_spei)) dir.create(output_dir_spei, recursive = TRUE)
  
  
  # 1. ROBUSTLY IDENTIFY AND LOAD THE WATER BALANCE STACKS
  wb.files <- list.files(output_dir_wb, full.names = TRUE, pattern = "\\.tif$")
  
  # Reference File: Ends in '_ref_bias_corrected.tif' AND contains the clean model name
  ref_file <- wb.files[
    grepl(model_name_clean, wb.files, ignore.case = TRUE) &
      grepl("_ref_bias_corrected\\.tif$", wb.files, ignore.case = TRUE)
  ]
  
  # Future File: Ends in '_future_bias_corrected.tif' AND contains the clean model name
  future_file <- wb.files[
    grepl(model_name_clean, wb.files, ignore.case = TRUE) &
      grepl("_future_bias_corrected\\.tif$", wb.files, ignore.case = TRUE)
  ]
  
  # Check if we found exactly one of each
  if (length(ref_file) != 1 || length(future_file) != 1) {
    cat(paste0("!!! ERROR: Could not find both bias-corrected files for ", model, ". Skipping SPEI calculation.\n"))
    cat(paste0(" Please check the '", output_dir_wb, "' folder.\n"))
    cat(paste0(" Expected internal fragment: '", model_name_clean, "'\n"))
    cat(paste0(" Expected file suffixes: '_ref_bias_corrected.tif' and '_future_bias_corrected.tif'\n"))
    next
  }
  
  tryCatch({
    cat("Loading and combining bias-corrected reference and future stacks...\n")
    ref_wb <- rast(ref_file)
    future_wb <- rast(future_file)
    
    # Combine them into one continuous time series (Reference + Future)
    water_balance <- c(ref_wb, future_wb)
    
    # Check time integrity
    if (is.null(time(water_balance))) {
      cat("!!! WARNING: Time attributes missing. SPEI calculation may fail. Skipping.\n")
      next
    }
    
    
    # 2. CALCULATE SPEI-12
    calculate_and_save_spei(
      input_stack = water_balance,
      scale = 12,
      model_name_clean = model_name_clean,
      output_dir_spei = output_dir_spei
    )
    
    # 3. CALCULATE SPEI-24
    calculate_and_save_spei(
      input_stack = water_balance,
      scale = 24,
      model_name_clean = model_name_clean,
      output_dir_spei = output_dir_spei
    )
    
    # Explicit cleanup
    rm(ref_wb, future_wb, water_balance)
    gc()
    
    cat(paste0("--- All SPEI Calculations Complete for: ", model, " ---\n"))
    
  }, error = function(e) {
    cat(paste0("\n!!! FATAL ERROR during SPEEI calculation for model ", model, ": ", e$message, "\n"))
    # Attempt to clear memory in case of crash
    gc()
  })
}

cat("\n\n*** Future SPEI (12m and 24m) calculation complete for all models. ***\n")



wb_file <- "F:/future_predictions/MIROC5 RCP85/wb_output/water_balance_MIROC5_RCP85_ref_bias_corrected.tif"
wb_stack <- rast(wb_file)
time(wb_stack)








library(terra)
library(SPEI)

# === SETTINGS ===
model <- "CNRM-CM5 RCP85"
base_dir <- "F:/future_predictions"
wb_dir   <- file.path(base_dir, model, "wb_output")
spei_dir <- file.path(base_dir, model, "diagnostic_spei_test")
if (!dir.exists(spei_dir)) dir.create(spei_dir, recursive = TRUE)

# === Load bias-corrected WB stack ===
ref_file  <- list.files(wb_dir, pattern = "_ref_bias_corrected\\.tif$", full.names = TRUE)
fut_file  <- list.files(wb_dir, pattern = "_future_bias_corrected\\.tif$", full.names = TRUE)
wb_stack  <- c(rast(ref_file), rast(fut_file))

# Limit to first ~100 months for speed
if (nlyr(wb_stack) > 100) wb_stack <- wb_stack[[1:100]]

# === Function to calculate SPEI-24 with custom reference period ===
calc_spei_test <- function(stack, ref_start, ref_end) {
  terra::app(stack, fun = function(x) {
    ts_data <- ts(x, start = c(1980, 1), frequency = 12)
    out <- tryCatch({
      spei_fit <- SPEI::spei(ts_data, scale = 24,
                             ref.start = ref_start, ref.end = ref_end,
                             na.rm = TRUE, verbose = FALSE)
      as.numeric(spei_fit$fitted)
    }, error = function(e) rep(NA, length(x)))
    out
  })
}

# === Calculate two versions ===
cat("Calculating SPEI-24 with 1980–2005 baseline...\n")
spei_1980 <- calc_spei_test(wb_stack, c(1980, 1), c(2005, 12))

cat("Calculating SPEI-24 with 1985–2005 baseline...\n")
spei_1985 <- calc_spei_test(wb_stack, c(1985, 1), c(2005, 12))

# === Compare distributions ===
vals_1980 <- values(spei_1980, na.rm = TRUE)
vals_1985 <- values(spei_1985, na.rm = TRUE)

par(mfrow = c(1, 2))
hist(vals_1980, breaks = 40, col = "grey70", main = "SPEI-24 (1980–2005 ref)",
     xlab = "SPEI value", xlim = c(-3, 3))
hist(vals_1985, breaks = 40, col = adjustcolor("blue", 0.5), main = "SPEI-24 (1985–2005 ref)",
     xlab = "SPEI value", xlim = c(-3, 3))

# === Summary stats ===
cat("\n=== Variance Comparison ===\n")
cat("1980–2005 SD:", round(sd(vals_1980, na.rm = TRUE), 3), "\n")
cat("1985–2005 SD:", round(sd(vals_1985, na.rm = TRUE), 3), "\n")

# === Quick spatial check (one month) ===
par(mfrow = c(1, 2))
plot(spei_1980[[1]], main = "SPEI-24 (1980–2005, layer 1)",
     col = rev(terrain.colors(100)))
plot(spei_1985[[1]], main = "SPEI-24 (1985–2005, layer 1)",
     col = rev(terrain.colors(100)))

