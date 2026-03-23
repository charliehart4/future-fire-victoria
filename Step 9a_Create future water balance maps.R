# Script 9a: Water Balance Calculation and Bias Correction for Future Projections
# This script performs three main tasks:
# 1. Calculates Water Balance (P-E) and SPEI for historical observed data.
# 2. Calculates raw Water Balance (P-E) for all four GCMs.
# 3. Applies quantile mapping bias correction to the GCM Water Balance data.

library(terra)
library(R.utils)
library(SPEI)
library(qmap)

# === CRITICAL MEMORY/DISK MANAGEMENT ===
options(scipen=999) # Prevents scientific notation in numbers
# 1. Set the temp directory to the dedicated F: drive path.
terraOptions(tempdir = "/vic_fire_mapping/terra_temp")
# 2. Tell terra to use a HUGE amount of memory (1 Terabyte) before trying to load
# raster blocks into RAM. This effectively forces it to use the disk for all operations.
terraOptions(memmax = 1e+12) 
gc()
# =======================================


# ======================================================================
# GLOBAL SETUP
# ======================================================================

# Define all climate models to process
CLIMATE_MODELS <- c(
  "ACCESS1-0 RCP85", # Corrected model name to match the file path F:/future_predictions/ACCESS1-0 RCP85
  "GFDL-ESM2M RCP85", 
  "CNRM-CM5 RCP85", 
  "MIROC5 RCP85"
)

# Define global paths
future.climate = here("data", "future_predictions")
temp_dir = "F:/vic_fire_mapping/terra_temp"
# Load a single raster to serve as a template for projection and cropping
baserast = rast(here("data", "baseline_observed_predstack.tif"))

# Helper function to sort files based on the date embedded in the name
sort_files_by_date <- function(files, pattern_prefix) {
  # Extract the YYYYMM part from the filename
  dates_str <- gsub(paste0(".*", pattern_prefix, "_?(\\d{6}).*"), "\\1", basename(files))
  # Create an ordered index
  return(files[order(dates_str)])
}

# Define observed data file locations (Used by Task 1: calculate_spei_raster)
rainfall_files <- list.files(here("data", "agcd_climate_data", "precip"),pattern = "\\.nc$", 
                             full.names = TRUE)

rainfall_files <- rainfall_files[grepl("r005", rainfall_files)]

et_files <- list.files(here("data", "awo_climate_data", "etot"),
                       pattern = "\\.nc$", 
                       full.names = TRUE)


# ======================================================================
#### 1) CALCULATE OBSERVED WATER BALANCE AND SPEI 
# ======================================================================

calculate_spei_raster <- function(rainfall_files, et_files, output_dir, temp_directory = temp_dir) {
  
  output_wb_file <- here("data", "water_balance_obs_vic.tif")
  
  # Check if observed WB file already exists
  if (file.exists(output_wb_file)) {
    cat(paste0("✅ Observed Water Balance file exists. Loading and skipping calculation.\n"))
    water_balance <- rast(output_wb_file)
  } else {
    # Read raster stacks
    cat("Reading observed rainfall rasters...\n")
    rainfall_stack <- rast(rainfall_files)
    cat("Reading observed evapotranspiration rasters...\n")
    et_stack <- rast(et_files)
    
    # Make ET monthly (sum daily data)
    cat("Aggregating ET to monthly sums...\n")
    layer_dates <- terra::time(et_stack)
    month_index <- format(layer_dates, "%Y-%m")
    # The following tapp command ensures ET is aggregated to monthly sums
    monthly_et_sum <- tapp(et_stack, index = month_index, fun = sum, na.rm = TRUE)
    et_stack <- monthly_et_sum
    
    # Ensure rasters line up
    rainfall_stack <- project(rainfall_stack, et_stack) 
    
    # Calculate Water Balance
    cat("Calculating Water Balance (P - E) \n")
    
    water_balance = rainfall_stack - et_stack
    water_balance = project(water_balance, y = crs(baserast))
    water_balance = crop(water_balance, baserast)
    
    month_index <- format(time(rainfall_stack), "%Y-%m")
    time(water_balance) <- as.Date(paste0(month_index, "-01"))
    
    writeRaster(water_balance, filename = output_wb_file, overwrite=TRUE)
    cat(paste("Observed Water Balance saved to:", output_wb_file, "\n"))
  }
  
  # --- START OF SPEI CHECK ---
  # Check if SPEI outputs exist by looking for any file starting with 'spei_12m_'
  spei_files_exist <- length(list.files(output_dir, pattern = "^spei_12m_.*\\.tif$", full.names = FALSE)) > 0
  
  if (spei_files_exist) {
    cat(paste0("✅ Observed SPEI files appear to exist in ", output_dir, ". Skipping SPEI calculation.\n"))
    return(water_balance) # Return WB for use in bias correction
  }
  # --- END OF SPEI CHECK ---
  
  # Create SPEI output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Calculate SPEI-12
  cat("Calculating SPEI-12m \n")
  spei_stack <- terra::app(water_balance, fun = function(x, scale = 12, na.rm = TRUE, ...) {
    res <- tryCatch({
      # Start date must match the actual data start date (1980, 1 for AGCD data)
      ts_data <- ts(x, start = c(1980, 1), frequency = 12)	
      spei_out <- SPEI::spei(ts_data, scale = 12, 
                             ref.start = c(1981, 1), ref.end = c(2005, 12), 
                             na.rm = TRUE, verbose = FALSE, ...)
      as.numeric(spei_out$fitted)
    }, error = function(e) rep(NA, length(x)))
    res
  })
  
  # Calculate SPEI-24
  cat("Calculating SPEI-24m \n")
  spei24_stack <- terra::app(water_balance, fun = function(x, scale = 24, na.rm = TRUE, ...) {
    res <- tryCatch({
      ts_data <- ts(x, start = c(1980, 1), frequency = 12)
      spei_out <- SPEI::spei(ts_data, scale = 24, 
                             ref.start = c(1981, 1), ref.end = c(2005, 12), 
                             na.rm = TRUE, verbose = FALSE, ...)
      as.numeric(spei_out$fitted)
    }, error = function(e) rep(NA, length(x)))
    res
  })
  
  # Assign time information
  if (!is.null(time(water_balance))) {
    time(spei_stack) <- time(water_balance)
    time(spei24_stack) <- time(water_balance)
  }
  
  # Save individual SPEI layers
  cat("Saving SPEI results...\n")
  
  for (i in 1:nlyr(spei_stack)) {
    t = time(spei_stack[[i]])
    writeRaster(spei_stack[[i]], 
                filename = file.path(output_dir, paste0("spei_12m_", t, ".tif")),
                overwrite = TRUE)
    writeRaster(spei24_stack[[i]], 
                filename = file.path(output_dir, paste0("spei_24m_", t, ".tif")),
                overwrite = TRUE)
  }
  
  cat("SPEI calculation complete!\n")
  return(water_balance) # Return WB for use in bias correction
}


# ======================================================================
#### 2) GCM WATER BALANCE CALCULATION
# ======================================================================

calculate_water_balance <- function(model_name, baserast, temp_directory = temp_dir) {
  
  cat(paste0("\n=== Processing Raw Water Balance for model: ", model_name, " ===\n"))
  
  model_folder = file.path(future.climate, model_name)
  # Sanitize model name for filename (replace spaces, hyphens, dots with underscores)
  model_name_clean <- gsub(" |-|\\.", "_", model_name)
  output_wb_file <- file.path(model_folder, "processed", paste0("water_balance_", model_name_clean, ".tif"))
  
  # Check if the output file already exists - USE A RETRY MECHANISM HERE
  max_retries <- 3
  for (attempt in 1:max_retries) {
    if (file.exists(output_wb_file)) {
      tryCatch({
        wb_rast <- rast(output_wb_file)
        cat(paste0("✅ Raw WB file already exists at ", output_wb_file, ". Loading and skipping calculation.\n"))
        return(wb_rast)
      }, error = function(e) {
        # If reading fails (e.g., corrupt file from previous failed run), delete it and retry calculation
        cat(paste0("⚠️ Failed to read existing file (Attempt ", attempt, "). Deleting and re-calculating.\n"))
        file.remove(output_wb_file)
      })
    }
    break # Exit loop if file check or deletion was successful
  }
  
  
  # Create necessary directories
  output_dir_proc <- file.path(model_folder, "processed")
  if (!dir.exists(output_dir_proc)) dir.create(output_dir_proc, recursive = TRUE)
  if (!dir.exists(temp_directory)) dir.create(temp_directory, recursive = TRUE)
  
  
  # 1. Get all files for the current model
  all_files <- list.files(model_folder, pattern = "\\.nc\\.gz$", full.names = TRUE)
  
  rainfall_files <- all_files[grepl("rr", all_files)]
  et_files <- all_files[grepl("etot", all_files)]
  
  # Sort files chronologically
  rainfall_files <- sort_files_by_date(rainfall_files, "rr")
  et_files <- sort_files_by_date(et_files, "etot")
  
  if (length(rainfall_files) != length(et_files) || length(rainfall_files) == 0) {
    stop(paste0("File mismatch for ", model_name, ". Found ", length(rainfall_files), " rr and ", length(et_files), " etot files."))
  }
  
  # 2. Decompress and Read Data
  cat("Decompressing and reading rasters...\n")
  
  temp_nc_files <- c()
  
  # Decompression helper function
  decompress_file <- function(file_path) {
    temp_path <- file.path(temp_directory, paste0(tools::file_path_sans_ext(basename(file_path), compression = TRUE), "_temp.nc"))
    R.utils::gunzip(file_path, destname = temp_path, remove = FALSE, overwrite = TRUE)
    temp_nc_files <<- c(temp_nc_files, temp_path) # Store for cleanup
    return(temp_path)
  }
  
  temp_rainfall_files <- sapply(rainfall_files, decompress_file)
  temp_et_files <- sapply(et_files, decompress_file)
  
  # Load the massive raster stacks (rely on terra's disk-backed implementation)
  rainfall_stack <- rast(temp_rainfall_files)
  et_stack <- rast(temp_et_files)
  
  # Apply filters and set time
  rainfall_stack[rainfall_stack < 0] <- NA
  et_stack[et_stack > 2000] <- NA
  
  # Extract dates to set time for the terra stacks	
  rr_dates_str <- gsub(".*rr_?(\\d{6}).*", "\\1", basename(rainfall_files))
  rr_formatted_dates <- as.Date(paste0(rr_dates_str, "01"), format = "%Y%m%d")
  
  terra::time(rainfall_stack) <- rr_formatted_dates
  terra::time(et_stack) <- rr_formatted_dates	
  
  # 3. Calculate Water Balance
  cat("Calculating Water Balance and reprojecting...\n")
  water_balance = rainfall_stack - et_stack
  water_balance = project(water_balance, y = crs(baserast))
  water_balance = crop(water_balance, baserast)
  
  # 4. Save and Cleanup (Retry mechanism applied here)
  for (attempt in 1:max_retries) {
    tryCatch({
      cat(paste0("Attempting to save raw WB file (Attempt ", attempt, ")...\n"))
      # Use overwrite = TRUE to prevent errors if a corrupt file exists from a previous failed run
      writeRaster(water_balance, filename = output_wb_file, overwrite = TRUE)
      
      # Explicitly reload the file to ensure the disk operation succeeded
      wb_rast_out <- rast(output_wb_file)	
      
      cat("Cleaning up temporary files...\n")
      if (length(temp_nc_files) > 0) {
        unlink(temp_nc_files)
      }
      
      cat(paste0("WB calculation complete for ", model_name, "!\n"))
      return(wb_rast_out) # Return the successfully saved and reloaded file
      
    }, error = function(e) {
      # Handle I/O failure
      cat(paste0("!!! Error during save/reload for ", model_name, " on attempt ", attempt, ": ", e$message, "\n"))
      if (attempt < max_retries) {
        sleep_time <- 2 ^ attempt
        cat(paste0("Retrying in ", sleep_time, " seconds...\n"))
        Sys.sleep(sleep_time)
      } else {
        stop(paste0("Failed to save/reload raw WB after ", max_retries, " attempts for ", model_name))
      }
    })
  }
}

# ======================================================================
#### 3) BIAS CORRECTION FUNCTIONS
# ======================================================================

# Function to perform quantile mapping bias correction on the reference period
wb_bias_correction <- function(model_wb_stack,
                               observed_wb_stack,
                               output_dir,
                               name,
                               model_name_clean,
                               method = "QUANT") {
  
  # Output filename includes the clean model name for clarity
  output_filename <- file.path(output_dir, paste0(name, "_", model_name_clean, "_ref_bias_corrected.tif"))
  if (file.exists(output_filename)) {
    cat(paste0("✅ Bias-corrected reference WB already exists at ", output_filename, ". Loading and skipping.\n"))
    return(rast(output_filename))
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Resample observed data to match model grid if needed
  if (!identical(ext(model_wb_stack), ext(observed_wb_stack)) ||
      !identical(res(model_wb_stack), res(observed_wb_stack))) {
    cat("Resampling observed reference data to match model grid...\n")
    observed_wb_stack <- terra::project(observed_wb_stack, model_wb_stack, method = "bilinear")
  }
  
  cat("Performing quantile mapping bias correction on reference period...\n")
  
  n_time <- nlyr(model_wb_stack)
  
  # Define the correction function to apply to each pixel
  correct_pixel <- function(values) {
    model_ts <- values[1:n_time]
    observed_ts <- values[(n_time + 1):(2 * n_time)]
    
    # Check for sufficient non-NA data
    if (sum(!is.na(model_ts)) < 10 || sum(!is.na(observed_ts)) < 10) {
      return(model_ts)
    }
    
    tryCatch({
      # Fit correction using reference period
      qm_fit <- fitQmap(obs = observed_ts, mod = model_ts, method = method,
                        wet.day = FALSE, qstep = 0.01)
      # Apply correction to the reference period (model data)
      return(doQmap(model_ts, qm_fit))
    }, error = function(e) {
      # warning(paste("Error", e)) # Using warning instead of message to be visible
      return(model_ts) # Return original if correction fails
    })
  }
  
  # Combine model and observed stacks
  combined_stack <- c(model_wb_stack, observed_wb_stack)
  
  # Apply correction function using terra::app()
  corrected_stack <- app(combined_stack, fun = correct_pixel)
  
  names(corrected_stack) <- paste0("WB_corrected_", time(model_wb_stack))
  
  if (!is.null(time(model_wb_stack))) {
    time(corrected_stack) <- time(model_wb_stack)
  }
  
  cat("Saving bias-corrected reference WB data...\n")
  
  writeRaster(corrected_stack,
              filename = output_filename,
              overwrite = TRUE)
  
  cat("Bias correction for reference period complete!\n")
  return(corrected_stack)
}


# Function to apply correction to future projections
apply_wb_correction <- function(future_wb_stack,
                                reference_model_stack,
                                reference_observed_stack,
                                output_dir,
                                name,
                                model_name_clean,
                                method = "QUANT") {
  
  # Output filename includes the clean model name
  output_filename <- file.path(output_dir, paste0(name, "_", model_name_clean, "_future_bias_corrected.tif"))
  if (file.exists(output_filename)) {
    cat(paste0("✅ Final future bias-corrected WB already exists at ", output_filename, ". Loading and skipping.\n"))
    return(rast(output_filename))
  }
  
  cat("Applying bias correction to future projections...\n")
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 1. Ensure all reference stacks align
  if (!identical(ext(reference_model_stack), ext(reference_observed_stack)) ||
      !identical(res(reference_model_stack), res(reference_observed_stack))) {
    cat("Resampling observed reference data to match model grid...\n")
    reference_observed_stack <- terra::project(reference_observed_stack, reference_model_stack, method = "bilinear")
  }
  
  # 2. Ensure future stack aligns with reference model grid
  if (!identical(ext(future_wb_stack), ext(reference_model_stack)) ||
      !identical(res(future_wb_stack), res(reference_model_stack))) {
    cat("Resampling future data to match reference grid...\n")
    future_wb_stack <- terra::project(future_wb_stack, reference_model_stack, method = "bilinear")
  }
  
  n_future <- nlyr(future_wb_stack)
  n_ref <- nlyr(reference_model_stack)
  
  cat("Processing", n_future, "future layers with", n_ref, "reference layers...\n")
  
  correct_future_pixel <- function(values) {
    # values structure: [future_ts, ref_model_ts, ref_obs_ts]
    future_ts <- values[1:n_future]
    ref_model_ts <- values[(n_future + 1):(n_future + n_ref)]
    ref_obs_ts <- values[(n_future + n_ref + 1):(n_future + 2*n_ref)]
    
    # Check for sufficient data in reference period
    if (sum(!is.na(ref_model_ts)) < 10 || sum(!is.na(ref_obs_ts)) < 10) {
      return(future_ts)	# Return original future values if insufficient reference data
    }
    
    # Skip if future data is all NA
    if (all(is.na(future_ts))) {
      return(future_ts)
    }
    
    tryCatch({
      # Fit correction using reference period
      qm_fit <- fitQmap(obs = ref_obs_ts, mod = ref_model_ts, method = method,
                        wet.day = FALSE, qstep = 0.01)
      
      # Apply correction to future data
      corrected_future <- doQmap(future_ts, qm_fit)
      return(corrected_future)
      
    }, error = function(e) {
      # If correction fails, return original future values
      # warning(paste("Correction failed for pixel:", e$message))
      return(future_ts)
    })
  }
  
  # Combine all stacks in the correct order: [FUTURE, REF_MODEL, REF_OBS]
  combined_stack <- c(future_wb_stack, reference_model_stack, reference_observed_stack)
  
  cat("Applying bias correction...\n")
  
  # Apply correction using terra::app()
  corrected_future <- app(combined_stack, fun = correct_future_pixel)
  
  # Set proper names and time
  names(corrected_future) <- paste0("wb_future_corrected_", time(future_wb_stack))
  if (!is.null(time(future_wb_stack))) {
    time(corrected_future) <- time(future_wb_stack)
  }
  
  cat("Saving results to:", output_filename, "\n")
  
  # Save results
  writeRaster(corrected_future,
              filename = output_filename,
              overwrite = TRUE)
  
  cat("Future projection correction complete!\n")
  return(corrected_future)
}

# ======================================================================
#### 4) EXECUTION LOOP FOR ALL GCMS
# ======================================================================

# Define reference period for bias correction
ref_start_date <- as.Date("1980-01-01")
ref_end_date <- as.Date("2006-01-01")
# Note: future_start_date is defined to pick up layers starting from 2050, 
# you may want to adjust this if your future period starts earlier (e.g., "2006-01-01")
future_start_date <- as.Date("2050-01-01") 

# Extract the observed reference period data once
cat("\n--- Task 1: Calculating/Loading Observed Water Balance and SPEI ---\n")
wb_observed_stack_full <- calculate_spei_raster(
  rainfall_files, 
  et_files, 
  output_dir = here("data", "spei_output"), 
  temp_directory = temp_dir
)

wb_observed_stack_reference <- subset(wb_observed_stack_full, 
                                      time(wb_observed_stack_full) >= ref_start_date &
                                        time(wb_observed_stack_full) < ref_end_date)


for (model in CLIMATE_MODELS) {
  tryCatch({
    
    cat(paste0("\n\n#####################################################\n"))
    cat(paste0("### STARTING MODEL: ", model, " ###\n"))
    cat(paste0("#####################################################\n"))
    
    model_name_clean <- gsub(" |-|\\.", "_", model)
    output_dir_wb <- file.path(future.climate, model, "wb_output")
    if (!dir.exists(output_dir_wb)) dir.create(output_dir_wb, recursive = TRUE)
    
    
    # 4.1) CALCULATE RAW GCM WATER BALANCE (P-E)
    # This function handles decompression, calculation, and disk-backed storage.
    wb_access_stack_full <- calculate_water_balance(model, baserast, temp_dir)
    
    # 4.2) SUBSET GCM DATA INTO REFERENCE AND FUTURE PERIODS
    wb_access_stack_reference <- subset(wb_access_stack_full, 
                                        time(wb_access_stack_full) >= ref_start_date &
                                          time(wb_access_stack_full) < ref_end_date)
    
    wb_access_stack_future <- subset(wb_access_stack_full, 
                                     time(wb_access_stack_full) >= future_start_date)
    
    # Free up the combined stack memory
    rm(wb_access_stack_full)
    gc()
    
    # 4.3) BIAS CORRECT REFERENCE PERIOD (FITS THE QMAP MODEL)
    cat(paste0("\n--- Task 3a: Bias Correcting Reference Period for ", model, " ---\n"))
    corrected_reference <- wb_bias_correction(
      model_wb_stack = wb_access_stack_reference,
      observed_wb_stack = wb_observed_stack_reference,
      output_dir = output_dir_wb,
      name = "water_balance",
      model_name_clean = model_name_clean,
      method = "QUANT"
    )
    
    # 4.4) APPLY BIAS CORRECTION TO FUTURE PROJECTIONS
    cat(paste0("\n--- Task 3b: Applying Bias Correction to Future Projections for ", model, " ---\n"))
    corrected_future <- apply_wb_correction(
      future_wb_stack = wb_access_stack_future,
      reference_model_stack = wb_access_stack_reference, # GCM Ref data (used for fitting)
      reference_observed_stack = wb_observed_stack_reference, # Observed Ref data (used for fitting)
      output_dir = output_dir_wb,
      name = "water_balance",
      model_name_clean = model_name_clean,
      method = "QUANT"
    )
    
    # Clean up model-specific stacks
    rm(wb_access_stack_reference, wb_access_stack_future, corrected_reference, corrected_future)
    gc()
    
  }, error = function(e) {
    cat(paste0("\n!!! FATAL ERROR processing model ", model, ": ", e$message, "\n"))
  })
}

cat("\n\n*** Full multi-model water balance and bias correction workflow complete! ***\n")






