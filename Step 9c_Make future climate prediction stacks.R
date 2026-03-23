# Script 9c: Prepare Final Future Prediction Stacks (FINAL VERSION - ENHANCED ERROR REPORTING & SPEI SKIP)
library(terra)
library(R.utils)

# === GLOBAL CONFIGURATION ===

# Define all climate models to process (MUST match the folder names on disk)
CLIMATE_MODELS <- c(
  "ACCESS1-0 RCP85",
  "GFDL-ESM2M RCP85",
  "CNRM-CM5 RCP85",
  "MIROC5 RCP85"
)

# Define all years to process
PERIODS <- list(
  baseline_modelled = 1985:2005,
  future_modelled = 2081:2099,
  # We use the union of years for the annual processing phase (Phase 1)
  all_modelled_years = unique(c(1985:2005, 2081:2099))
)

# Define file paths
here("data", "future_predictions")
here("data", "covariates")
here("data", "output_data")
temp_dir = "F:/vic_fire_mapping/terra_temp"

# Ensure temp directory exists
if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)

# Load the mask and static covariates once for Phase 2/3
cat("\nLoading study region mask and static covariate stack...\n")
cov.stack.static <- rast(file.path(cov.path, "masked/masked_static_covariate_stack.tif"))
mask.ras <- cov.stack.static[[1]] # Use the first layer as the mask
cat(paste0("✅ Static stack (", nlyr(cov.stack.static), " layers) and mask loaded successfully.\n"))

# ======================================================================
# CORE FUNCTIONS - PHASE 1: ANNUAL AGGREGATION
# (These functions are kept separate but modified to take model_name)
# ======================================================================

# --- FFDI Mean ---
mean_summer_rasters_gz <- function(gz_folder, year, model_id, temp_directory = temp_dir, pattern = "FFDI", output_folder = NULL) {
  prev_year <- year - 1
  # File pattern for FFDI mean is unique (not month-specific) and includes 'QME'
  months_to_match <- c(paste0(prev_year, c("11", "12")), paste0(year, c("01", "02", "03")))
  
  # List all .gz files in the folder
  all_gz_files <- list.files(gz_folder, pattern = "\\.nc\\.gz$", full.names = TRUE)
  
  summer_gz_files <- character(5)
  for (i in 1:5) {
    month_str <- months_to_match[i]
    # CRITICAL FIX: Match the FFDI QME files
    regex_pattern <- paste0(pattern, "_", month_str, ".*QME\\.nc\\.gz$")
    found_file <- grep(regex_pattern, all_gz_files, value = TRUE)
    
    if (length(found_file) != 1) {
      warning(paste("FFDI Mean: No file found for month", month_str))
      return(NULL)
    }
    summer_gz_files[i] <- found_file
  }
  
  if (is.null(output_folder)) output_folder <- file.path(gz_folder, "processed")
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  output_file <- file.path(output_folder, sprintf("FFDI_summer_mean_%d.nc", year))
  
  if (file.exists(output_file)) {
    cat("Skipping FFDI Mean:", basename(output_file), "already exists.\n")
    return(output_file)
  }
  
  # Decompress files
  temp_nc_files <- character(length(summer_gz_files))
  for(i in seq_along(summer_gz_files)) {
    temp_nc_files[i] <- file.path(temp_directory, paste0(tools::file_path_sans_ext(basename(summer_gz_files[i]), compression = TRUE), "_ffdi_mean_", year, "_", i, ".nc"))
    R.utils::gunzip(summer_gz_files[i], destname = temp_nc_files[i], remove = FALSE, overwrite = TRUE)
  }
  
  # Read rasters and average them
  summer_rasters <- tryCatch({
    rast(temp_nc_files)
  }, error = function(e) {
    warning(paste("Error reading rasters for FFDI Mean year", year, ":", e$message))
    return(NULL)
  })
  if(is.null(summer_rasters)) return(NULL)
  
  summer_mean <- mean(summer_rasters)
  summer_mean[summer_mean < 0] <- NA
  
  writeRaster(summer_mean, output_file, overwrite = TRUE)
  unlink(temp_nc_files)
  return(output_file)
}

# --- FFDI 95 / KBDI 95 (Same logic, different pattern) ---
sum_summer_rasters_gz <- function(gz_folder, year, model_id, temp_directory = temp_dir, pattern = "FFDI_gt95perc_", output_file_base = "FFDI_95_summer_sum", output_folder = NULL) {
  prev_year <- year - 1
  months <- c(sprintf("%d11", prev_year), sprintf("%d12", prev_year),
              sprintf("%d01", year), sprintf("%d02", year), sprintf("%d03", year))
  
  all_gz_files <- list.files(gz_folder, pattern = "\\.nc\\.gz$", full.names = TRUE)
  
  # CRITICAL FIX: Generalize the file pattern match
  summer_gz_files <- c()
  for (month_str in months) {
    # Pattern must be model-agnostic, checking for month and pattern (e.g., FFDI_gt95perc_198411)
    regex_pattern <- paste0(pattern, month_str, ".*\\.nc\\.gz$")
    found_file <- grep(regex_pattern, all_gz_files, value = TRUE)
    if (length(found_file) == 1) {
      summer_gz_files <- c(summer_gz_files, found_file)
    }
  }
  
  if (length(summer_gz_files) != 5) {
    warning(paste(output_file_base, "No 5 matching monthly files found for summer of", year, ". Found:", length(summer_gz_files)))
    return(NULL)
  }
  
  if (is.null(output_folder)) output_folder <- file.path(gz_folder, "processed")
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  output_file <- file.path(output_folder, sprintf("%s_%d.nc", output_file_base, year))
  
  if (file.exists(output_file)) {
    cat(paste("Skipping", output_file_base, ":", basename(output_file), "already exists.\n"))
    return(output_file)
  }
  
  # Decompress files
  temp_nc_files <- character(length(summer_gz_files))
  for(i in seq_along(summer_gz_files)) {
    temp_nc_files[i] <- file.path(temp_directory, paste0(tools::file_path_sans_ext(basename(summer_gz_files[i]), compression = TRUE), "_", output_file_base, "_", year, "_", i, ".nc"))
    R.utils::gunzip(summer_gz_files[i], destname = temp_nc_files[i], remove = FALSE, overwrite = TRUE)
  }
  
  # Read rasters and sum them
  summer_sum <- tryCatch({
    summer_rasters <- rast(temp_nc_files)
    sum(summer_rasters, na.rm = TRUE)
  }, error = function(e) {
    warning(paste("Error reading/summing rasters for", output_file_base, year, ":", e$message))
    return(NULL)
  }, finally = {
    unlink(temp_nc_files)
  })
  
  if(is.null(summer_sum)) return(NULL)
  
  writeRaster(summer_sum, output_file, overwrite = TRUE)
  return(output_file)
}

# --- SPEI Mean (TIF stack) ---
mean_spei_summer <- function(spei_lag_stack, lag_length, year, output_folder){
  prev_year <- year - 1
  # CRITICAL FIX: The SPEI files should be read from TIFs generated separately
  dates <- c(as.Date(paste0("01-11-", prev_year), format = c("%d-%m-%Y")),
             as.Date(paste0("01-12-", prev_year), format = c("%d-%m-%Y")),
             as.Date(paste0("01-01-", year), format = c("%d-%m-%Y")),
             as.Date(paste0("01-02-", year), format = c("%d-%m-%Y")),
             as.Date(paste0("01-03-", year), format = c("%d-%m-%Y")))
  
  if (is.null(output_folder)) stop("Output folder must be defined for SPEI.")
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  output_file <- file.path(output_folder, paste0("spei", lag_length,"_summer_mean_", year,".tif"))
  
  if (file.exists(output_file)) {
    cat(paste("Skipping SPEI", lag_length, ":", basename(output_file), "already exists.\n"))
    return(output_file)
  }
  
  if (is.null(spei_lag_stack)) {
    warning(paste("SPEI", lag_length, "stack is NULL, cannot process year", year))
    return(NULL)
  }
  
  # Filter files matching the required summer months
  year_stack <- tryCatch({
    subset(spei_lag_stack, time(spei_lag_stack) %in% dates)
  }, error = function(e) {
    warning(paste("Error subsetting SPEI", lag_length, "stack for year", year, ":", e$message))
    return(NULL)
  })
  
  if (is.null(year_stack) || nlyr(year_stack) != 5) {
    warning(paste("SPEI", lag_length, "Found", nlyr(year_stack), "months, expected 5. Skipping year", year))
    return(NULL)
  }
  
  # Create average for summer period
  spei_summer_mean = mean(year_stack, na.rm=TRUE)
  
  writeRaster(spei_summer_mean, output_file, overwrite = TRUE)
  return(output_file)
}

# --- Thunderstorm (BTE) ---
is_leap_year <- function(year) {
  return((year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0))
}

sum_summer_thunderstorm <- function(gz_folder, year, model_id, temp_directory = temp_dir, output_folder = NULL) {
  prev_year <- year - 1
  model_base_name <- gsub(" RCP85", "", model_id) # e.g. ACCESS1-0
  pattern <- paste0("BTE_", model_base_name, "_")
  
  if (is.null(output_folder)) output_folder <- file.path(gz_folder, "processed")
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  output_file <- file.path(output_folder, sprintf("thunderstorm_summer_days_%d.nc", year))
  
  if (file.exists(output_file)) {
    cat("Skipping Thunderstorm:", basename(output_file), "already exists.\n")
    return(output_file)
  }
  
  # Find the files for both years
  all_gz_files <- list.files(gz_folder, pattern = "\\.nc\\.gz$", full.names = TRUE)
  
  # CRITICAL FIX: The BTE files are generally year-based (e.g., BTE_ACCESS1-0_1984.nc.gz)
  summer_gz_files <- c(
    grep(paste0(pattern, prev_year, ".*\\.nc\\.gz$"), all_gz_files, value = TRUE),
    grep(paste0(pattern, year, ".*\\.nc\\.gz$"), all_gz_files, value = TRUE)
  )
  
  # Decompression (must be done BEFORE checking layer count as we need the layer metadata)
  temp_nc_files <- character(length(summer_gz_files))
  if (length(summer_gz_files) != 2) {
    warning(paste("Thunderstorm: Expected 2 raw files (prev/current year). Found:", length(summer_gz_files), ". Skipping year."))
    return(NULL)
  }
  
  for(i in seq_along(summer_gz_files)) {
    temp_nc_files[i] <- file.path(temp_directory, paste0(tools::file_path_sans_ext(basename(summer_gz_files[i]), compression = TRUE), "_bte_", year, "_", i, ".nc"))
    R.utils::gunzip(summer_gz_files[i], destname = temp_nc_files[i], remove = FALSE, overwrite=TRUE)
  }
  
  # Load rasters
  r_prev <- rast(temp_nc_files[1])
  r_current <- rast(temp_nc_files[2]) # CRITICAL FIX: Load current year from the second file
  
  # Define date range (Nov 1 of prev year to Mar 31 of current year)
  start_date <- as.Date(paste0(prev_year, "-11-01"))
  end_date <- as.Date(paste0(year, "-03-31"))
  
  # Subset layers by time
  summer_prev <- subset(r_prev, time(r_prev) >= start_date & time(r_prev) <= as.Date(paste0(prev_year, "-12-31")))
  summer_curr <- subset(r_current, time(r_current) >= as.Date(paste0(year, "-01-01")) & time(r_current) <= end_date)
  
  # Merge Nov-Mar layers
  summer_rasters <- c(summer_prev, summer_curr)
  layers_per_day <- 4 # Assumed sub-daily data
  
  if (nlyr(summer_rasters) == 0) {
    warning(paste("Thunderstorm: Subset resulted in 0 layers for summer", year, ". Skipping."))
    unlink(temp_nc_files)
    return(NULL)
  }
  
  # Convert to daily presence/absence (max over 4 layers per day)
  num_days <- nlyr(summer_rasters) / layers_per_day
  if (num_days != round(num_days)) {
    warning(paste("Thunderstorm: Layer count is not a multiple of 4 (daily layers). Found:", nlyr(summer_rasters), ". Skipping."))
    unlink(temp_nc_files)
    return(NULL)
  }
  
  daily_presence <- rast(lapply(seq_len(num_days), function(d) {
    max(summer_rasters[[((d - 1) * layers_per_day + 1):(d * layers_per_day)]], na.rm = TRUE)
  }))
  
  # Sum across summer days
  summer_sum <- sum(daily_presence, na.rm = TRUE)
  
  writeRaster(summer_sum, output_file, overwrite = TRUE)
  unlink(temp_nc_files)
  return(output_file)
}

# ======================================================================
# CORE FUNCTION - PHASE 2: PROJECT AND MASK
# ======================================================================

align_and_mask_raster <- function(f, mask_ras, model_name, output_folder_name, cov_path) {
  if (grepl("\\.nc\\.aux\\.xml$", f) || grepl("\\.tif\\.aux\\.xml$", f)) return(NULL)	
  
  output_path_dir <- file.path(cov_path, "modelled_masked", model_name, output_folder_name)
  if (!dir.exists(output_path_dir)) dir.create(output_path_dir, recursive = TRUE)
  
  filename <- gsub("\\.(nc|tif)$", ".tif", basename(f)) # Standardize all output to .tif
  output_file <- file.path(output_path_dir, filename)
  
  if (file.exists(output_file)) {
    cat(paste0("	> SKIPPING Alignment: Masked output already saved to ", basename(output_file), ".\n"))
    return(NULL)
  }
  
  tryCatch({
    input_ras <- rast(f)
    if (nlyr(input_ras) == 0) return(NULL)
    
    input_ras <- terra::project(input_ras, crs(mask_ras))
    input_ras <- terra::resample(input_ras, mask_ras, method = "bilinear")
    input_ras <- terra::mask(input_ras, mask_ras)
    
    writeRaster(input_ras, output_file, overwrite = TRUE)
    cat(paste0("	✅ PHASE 2 SAVED: ", basename(output_file), "\n"))
    
  }, error = function(e) {
    warning(paste("Error processing", basename(f), ":", e$message))
  }, finally = {
    gc()
  })
}




# ======================================================================
# CORE FUNCTION - PHASE 3: CREATE PREDICTION STACKS (AGGRESSIVE MEMORY)
# ======================================================================

make_prediction_stack <- function(start_year, end_year, time_period_label, model_name, prediction_stack_static, cov_path, output_dir){
  cat(paste0("\n--- PHASE 3: Generating FINAL STACK for ", model_name, " (", time_period_label, ") ---\n"))
  
  # Define the appropriate base path for reading annual files
  if (model_name == "Observed") { # Special case for observed baseline stacks
    base_path <- file.path(cov_path, "masked")
  } else {
    base_path <- file.path(cov_path, "modelled_masked", model_name)
  }
  
  year_pattern <- paste0("(", paste(start_year:end_year, collapse = "|"), ")")
  output_stack_dir <- file.path(output_dir, "prediction_stacks")
  if (!dir.exists(output_stack_dir)) dir.create(output_stack_dir, recursive = TRUE)
  
  # ----------------------------------------------------------------------
  # 1. DEFINE FILENAME AND PATTERNS (MUST HAPPEN BEFORE FILE.EXISTS CHECK)
  # ----------------------------------------------------------------------
  if (model_name == "Observed") {
    # Observed Patterns (using your specific file names)
    ffdi_mean_pattern <- "FFDI_summer_mean_" 
    ffdi95_sum_pattern <- "FFDI_95_summer_sum_" 
    kbdi95_sum_pattern <- "KBDI_summer_sum_"     #  KBDI FIX
    spei12_pattern <- "SPEI12_summer_avg_"      #  SPEI FIX
    spei24_pattern <- "SPEI24_summer_avg_"      #  SPEI FIX
    bte_pattern <- "thunderstorm_summer_days_" 
    
    # Set final output filename for Observed
    filename <- file.path(output_stack_dir, paste0("baseline_", time_period_label, "_observed_predstack.tif"))
    
  } else {
    # Modelled Patterns (original, standard patterns)
    ffdi_mean_pattern <- "FFDI_summer_mean_"
    ffdi95_sum_pattern <- "FFDI_95_summer_sum_"
    kbdi95_sum_pattern <- "KBDI_95_summer_sum_"
    spei12_pattern <- "spei12_summer_mean_"
    spei24_pattern <- "spei24_summer_mean_"
    bte_pattern <- "thunderstorm_summer_days_"
    
    # Set final output filename for Modelled
    filename <- file.path(output_stack_dir, paste0(time_period_label, "_", gsub(" ", "_", model_name), "_modelled_predstack.tif"))
  }
  # ----------------------------------------------------------------------
  
  # 2. FILE EXISTENCE CHECK (Now that 'filename' is defined)
  if (file.exists(filename)) {
    cat(paste0("	> SKIPPING: Final Prediction Stack already saved to ", basename(filename), ".\n"))
    return(NULL)	
  }
  
  prediction.stack <- prediction_stack_static
  
  covariate_specs <- list(
    # Use the dynamically set patterns
    list(name = "ffdi_mean", folder = "ffdi_mean", pattern = ffdi_mean_pattern),
    list(name = "ffdi_95_days", folder = "ffdi95_sum", pattern = ffdi95_sum_pattern),
    list(name = "kbdi_95_days", folder = "kbdi95_sum", pattern = kbdi95_sum_pattern),
    list(name = "spei12_mean", folder = "spei12", pattern = spei12_pattern),
    list(name = "spei24_mean", folder = "spei24", pattern = spei24_pattern),
    list(name = "thunderstorm_days", folder = "thunderstorm", pattern = bte_pattern)
  )
  
  for (spec in covariate_specs) {
    # Searches for annual TIFs created in Phase 2/other scripts
    cov_files <- list.files(file.path(base_path, spec$folder),	
                            pattern = paste0(spec$pattern, year_pattern, "\\.tif$"),	
                            full.names = TRUE)
    
    if (length(cov_files) == 0) {
      warning(paste("No annual TIF files found for", spec$name, "in", file.path(base_path, spec$folder), ". Skipping."))
      next
    }
    
    # Load the full stack for the current covariate
    cov_stack <- rast(cov_files)
    if (grepl("spei", spec$name)) {
      cov_stack[is.infinite(cov_stack)] <- NA	
    }
    
    # Calculate the long-term mean
    cov_mean <- mean(cov_stack, na.rm = TRUE)
    names(cov_mean) <- spec$name
    
    # Append the new mean layer to the final stack
    prediction.stack <- c(prediction.stack, cov_mean)
    
    cat(paste0("✅ Added long-term mean for: ", spec$name, " (", length(cov_files), " years)\n"))
    
    # --- AGGRESSIVE MEMORY MANAGEMENT ---
    rm(cov_stack, cov_mean)
    gc(verbose = FALSE, full = TRUE)
    terra::tmpFiles(remove=TRUE)
    # ------------------------------------
  }
  
  writeRaster(prediction.stack, filename = filename, overwrite = TRUE)
  cat(paste0("--- ✅ PHASE 3 SAVED: FINAL STACK (", nlyr(prediction.stack), " layers): ", basename(filename), " ---\n"))
  
  return(prediction.stack)
}

# ======================================================================
# MASTER EXECUTION LOOP (Iterates over all four climate models)
# ======================================================================

# --- PHASE 0: Initialize Stacks (Only need to initialize once outside the loop) ---
# NOTE: The stacks will be re-defined inside the loop for each model
# spei12_stack <- NULL # Remove this line (and the spei24 one) from here
# spei24_stack <- NULL # Remove this line (and the spei12 one) from here


for (model in CLIMATE_MODELS) {
  
  cat(paste0("\n=================================================================\n"))
  cat(paste0("STARTING MODEL PROCESSING: ", model, "\n"))
  cat(paste0("=================================================================\n"))
  
  source_folder <- file.path(future.climate, model)
  model_id <- model
  
  
  # --------------------------------------------------------------------------
  # --- NEW LOCATION for PHASE 0: Pre-load SPEI stacks (Now Model-Specific) ---
  # --------------------------------------------------------------------------
  
  # 1. Define the model-specific SPEI path
  spei_path <- file.path(source_folder, "spei_output") 
  
  # 2. List all SPEI files in the current model's directory
  spei_files = list.files(spei_path, pattern = "\\.tif$", full.names = TRUE)
  
  # 3. Filter files for 12-month and 24-month SPEI
  # NOTE: We keep the model-specific files (e.g., 'spei_12m_GFDL-ESM2M_...')
  spei12_files = spei_files[grepl("spei_12m", spei_files)]
  spei24_files = spei_files[grepl("spei_24m", spei_files)]
  
  # 4. Load the raster stacks (These variables MUST be inside the loop)
  spei12_stack <- NULL
  spei24_stack <- NULL
  
  if (length(spei12_files) > 0) spei12_stack <- rast(spei12_files)
  if (length(spei24_files) > 0) spei24_stack <- rast(spei24_files)
  
  # 5. Set layer names
  if (!is.null(spei12_stack)) names(spei12_stack) <- time(spei12_stack)
  if (!is.null(spei24_stack)) names(spei24_stack) <- time(spei24_stack)
  

  
  # ----------------------------------------------------------------------
  # PHASE 1: PROCESS RAW GZ DATA (Annual Aggregation)
  # ----------------------------------------------------------------------
  
  cat("\n--- PHASE 1: RAW DATA PROCESSING AND ANNUAL AGGREGATION ---\n")
  
  for (year in PERIODS$all_modelled_years) {
    cat(paste0("\n*** Processing Year: ", year, " ***\n"))
    
    # 1. FFDI Mean
    ffdi_mean_file <- mean_summer_rasters_gz(
      gz_folder = source_folder, year = year, model_id = model_id, output_folder = file.path(source_folder, "processed")
    )
    gc()
    
    # 2. FFDI 95 Sum
    ffdi95_file <- sum_summer_rasters_gz(
      gz_folder = source_folder, year = year, model_id = model_id, pattern = "FFDI_gt95perc_", 
      output_file_base = "FFDI_95_summer_sum", output_folder = file.path(source_folder, "processed")
    )
    gc()
    
    # 3. KBDI 95 Sum
    kbdi95_file <- sum_summer_rasters_gz(
      gz_folder = source_folder, year = year, model_id = model_id, pattern = "KBDI_gt95perc_", 
      output_file_base = "KBDI_95_summer_sum", output_folder = file.path(source_folder, "processed")
    )
    gc()
    
    # 4. Thunderstorm Sum (BTE)
    # Note: Thunderstorm files may be in a different subfolder path, this assumes the main model folder
    # is the correct one containing the BTE files (BTE_ACCESS1-0_1984.nc.gz, etc.)
    bte_file <- sum_summer_thunderstorm(
      gz_folder = source_folder, year = year, model_id = model_id, output_folder = file.path(source_folder, "processed")
    )
    gc()
    
    # 5. SPEI Mean
    spei12_file <- mean_spei_summer(
      spei_lag_stack = spei12_stack, lag_length = "12", year = year, output_folder = file.path(source_folder, "processed")
    )
    gc()
    spei24_file <- mean_spei_summer(
      spei_lag_stack = spei24_stack, lag_length = "24", year = year, output_folder = file.path(source_folder, "processed")
    )
    gc()
  }
  
  # ----------------------------------------------------------------------
  # PHASE 2: ALIGN AND MASK ALL ANNUAL LAYERS
  # ----------------------------------------------------------------------
  
  cat("\n--- PHASE 2: ALIGNING AND MASKING ANNUAL RASTERS ---\n")
  processed_dir <- file.path(source_folder, "processed")
  
  annual_data_specs <- list(
    list(pattern = "FFDI_summer_mean", out_folder = "ffdi_mean"),	
    list(pattern = "FFDI_95_summer_sum", out_folder = "ffdi95_sum"),
    list(pattern = "KBDI_95_summer_sum", out_folder = "kbdi95_sum"),
    list(pattern = "thunderstorm_summer_days", out_folder = "thunderstorm"),
    list(pattern = "spei12_summer_mean", out_folder = "spei12"),
    list(pattern = "spei24_summer_mean", out_folder = "spei24")
  )
  
  for (spec in annual_data_specs) {
    # FFDI/KBDI/BTE are NC files; SPEI are TIF files
    file_extension_pattern <- ifelse(grepl("spei", spec$pattern), "\\.tif$", "\\.nc$")
    annual_files <- list.files(processed_dir, pattern = paste0(spec$pattern, ".*", file_extension_pattern), full.names = TRUE)
    
    for (f in annual_files) {
      align_and_mask_raster(f, mask.ras, model, spec$out_folder, cov.path)
    }
  }
  
  # ----------------------------------------------------------------------
  # PHASE 3: GENERATE FINAL PREDICTION STACKS (LONG-TERM AVERAGES)
  # ----------------------------------------------------------------------
  
  # 1. Baseline Modelled Historical (1985-2005)
  make_prediction_stack(
    start_year = PERIODS$baseline_modelled[1],
    end_year = PERIODS$baseline_modelled[length(PERIODS$baseline_modelled)],
    time_period_label = paste0("baseline_", PERIODS$baseline_modelled[1], "_", PERIODS$baseline_modelled[length(PERIODS$baseline_modelled)]),
    model_name = model,
    prediction_stack_static = cov.stack.static,
    cov_path = cov.path,
    output_dir = output_dir
  )
  
  # 2. Future Modelled (2081-2099)
  make_prediction_stack(
    start_year = PERIODS$future_modelled[1],
    end_year = PERIODS$future_modelled[length(PERIODS$future_modelled)],
    time_period_label = paste0("future_", PERIODS$future_modelled[1], "_", PERIODS$future_modelled[length(PERIODS$future_modelled)]),
    model_name = model,
    prediction_stack_static = cov.stack.static,
    cov_path = cov.path,
    output_dir = output_dir
  )
  
  cat(paste0("\n=================================================================\n"))
  cat(paste0("COMPLETED MODEL PROCESSING: ", model, "\n"))
  cat(paste0("=================================================================\n"))
  
  gc()
}

# ======================================================================
# FINAL STACKS: OBSERVED DATA (Independent of the model loop)
# ======================================================================

cat(paste0("\n=================================================================\n"))
cat(paste0("STARTING OBSERVED STACK GENERATION (No Phase 1/2)\n"))
cat(paste0("=================================================================\n"))

# Baseline - Observed Data (1985-2005)
make_prediction_stack(
  start_year = 1985, end_year = 2005,
  time_period_label = "1985_2005",
  model_name = "Observed", # Special identifier
  prediction_stack_static = cov.stack.static,
  cov_path = cov.path,
  output_dir = output_dir
)

# Baseline - Observed Data (2006-2023)
make_prediction_stack(
  start_year = 2006, end_year = 2023,
  time_period_label = "2006_2023",
  model_name = "Observed",
  prediction_stack_static = cov.stack.static,
  cov_path = cov.path,
  output_dir = output_dir
)

cat(paste0("\n=================================================================\n"))
cat(paste0("ALL PROCESSING COMPLETE\n"))
cat(paste0("=================================================================\n"))









