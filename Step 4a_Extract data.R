# =============================================================================
# 4a) Covariate Extraction and Merging (All Paths Updated)
# =============================================================================
library(terra)
library(dplyr)
library(readr)

message("🔄 Starting Covariate Extraction and Merging...")

random_sample = vect(here("covariates", "random_points.shp"))
cov.path <- here("covariates")
output_file_path <- here("output_data", "full_fire_covariates_no_folds.csv")

## Extract covariates
cov_stack = rast(file.path(cov.path, "masked/masked_static_covariate_stack.tif"))
cov_extracted = terra::extract(cov_stack, random_sample, xy=TRUE)

# Check correlation between variables
# NOTE: Ensure 'corrplot' is installed if you run this line outside of a working environment
# cors = cor(cov_extracted[,2:ncol(cov_extracted)], method="pearson", use="pairwise.complete.obs")
# corrplot::corrplot(cors)

## Extract FFDI values - mean
ffdi.maps = list.files(file.path(cov.path, "masked", "ffdi_mean"), pattern=".tif", full=T)
extract_ffdi_data = function(ffdi.ras.path, random.points){
  ffdi.ras = rast(ffdi.ras.path)
  ffdi.data = terra::extract(ffdi.ras, random.points, xy=TRUE)
  names(ffdi.data) <- c("ID", "ffdi_mean", "x","y")
  ffdi.data$year = readr::parse_number(names(ffdi.ras))
  return(ffdi.data)
}

ffdi_extract = lapply(ffdi.maps, FUN = extract_ffdi_data, random_sample)
ffdi_extract_long = do.call("rbind", ffdi_extract)




## Extract kbdi values (UPDATED PATH: kbdi95_sum)
kbdi.maps = list.files(file.path(cov.path, "masked", "kbdi95_sum"), pattern=".tif", full=T)
extract_kbdi_data = function(kbdi.ras.path, random.points){
  kbdi.ras = rast(kbdi.ras.path)
  kbdi.data = terra::extract(kbdi.ras, random.points, xy=TRUE)
  names(kbdi.data) <- c("ID", "kbdi_95_days", "x","y")
  kbdi.data$year = readr::parse_number(names(kbdi.ras))
  return(kbdi.data)
}

# --- CRITICAL SAFETY CHECK FOR KBDI (Still in place for robustness) ---
if (length(kbdi.maps) == 0) {
  message("⚠️ WARNING: No KBDI files found in the 'kbdi95_sum' directory. Creating an empty placeholder data frame to prevent crash.")
  # Create a data frame with zero rows but the correct column structure
  kbdi_extract_long = data.frame(ID = integer(), x = numeric(), y = numeric(), year = integer(), kbdi_95_days = numeric())
} else {
  kbdi_extract = lapply(kbdi.maps, FUN = extract_kbdi_data, random_sample)
  kbdi_extract_long = do.call("rbind", kbdi_extract)
}


## Extract spei values (UPDATED PATH: spei12)
spei.maps = list.files(file.path(cov.path, "masked", "spei12"), pattern=".tif", full=T)
extract_spei_data = function(spei.ras.path, random.points){
  spei.ras = rast(spei.ras.path)
  spei.data = terra::extract(spei.ras, random.points, xy=TRUE)
  names(spei.data) <- c("ID", "spei12_mean", "x","y")
  # Correctly reading year from name that might contain other text
  spei.data$year = readr::parse_number(gsub("SPEI12","", names(spei.ras))) 
  return(spei.data)
}

spei_extract = lapply(spei.maps, FUN = extract_spei_data, random_sample)
spei12_extract_long = do.call("rbind", spei_extract)

## Extract spei values (UPDATED PATH: spei24)
spei.maps = list.files(file.path(cov.path, "masked", "spei24"), pattern=".tif", full=T)
extract_spei_data = function(spei.ras.path, random.points){
  spei.ras = rast(spei.ras.path)
  spei.data = terra::extract(spei.ras, random.points, xy=TRUE)
  names(spei.data) <- c("ID", "spei24_mean", "x","y")
  # Correctly reading year from name that might contain other text
  spei.data$year = readr::parse_number(gsub("SPEI24","", names(spei.ras))) 
  return(spei.data)
}

spei_extract = lapply(spei.maps, FUN = extract_spei_data, random_sample)
spei24_extract_long = do.call("rbind", spei_extract)


## Extract TSF (Time Since Fire) values - ALIGNED (CRITICAL STEP)
tsf.maps <- list.files(here("covariates", "raw", "tsf"), pattern = "\\.tif$", full.names = TRUE)

extract_tsf_data <- function(tsf.ras.path, random.points){
  tsf.ras  <- terra::rast(tsf.ras.path)
  tsf.data <- terra::extract(tsf.ras, random.points, xy = TRUE)
  names(tsf.data) <- c("ID", "tsf", "x", "y")
  tsf.data$year <- readr::parse_number(names(tsf.ras))
  return(tsf.data)
}

tsf_extract     <- lapply(tsf.maps, FUN = extract_tsf_data, random.points = random_sample)
tsf_extract_long <- do.call("rbind", tsf_extract)

# The TSF value for year Y (tsf_Y.tif) is the correct predictor for the fire in year Y (burnt_Y.tif).
# NO LAG needed here. We just rename it for the join.
tsf_extract_aligned <- tsf_extract_long 
message("TSF extracted and aligned (No lag applied).")

## Extract Thunderstorm values
thunder.maps = list.files(file.path(cov.path, "masked", "thunderstorm"), pattern=".tif", full=T)
extract_thunder_data = function(thunder.ras.path, random.points){
  thunder.ras = rast(thunder.ras.path)
  thunder.data = terra::extract(thunder.ras, random.points, xy=TRUE)
  names(thunder.data) <- c("ID", "thunderstorm_days", "x","y")
  thunder.data$year = readr::parse_number(names(thunder.ras))
  return(thunder.data)
}

thunder_extract = lapply(thunder.maps, FUN = extract_thunder_data, random_sample)
thunder_extract_long = do.call("rbind", thunder_extract)




# =============================================================================
# 2. Extract FFDI 95th Percentile (FIXED)
# =============================================================================

# Helper function to correctly extract the year from the filename path
extract_ffdi95_data = function(ffdi.ras.path, random.points){
  # 1. Load raster and extract data
  ffdi.ras = rast(ffdi.ras.path)
  ffdi.data = terra::extract(ffdi.ras, random.points, xy=TRUE)
  names(ffdi.data) <- c("ID", "ffdi_95_days", "x","y")
  
  # 2. CRITICAL FIX: Extract the year from the filename (e.g., '1981.tif')
  # This avoids the error when using layer names which might return '1' or '95'.
  year_string = gsub(".*_sum_", "", basename(ffdi.ras.path))
  ffdi.data$year = readr::parse_number(year_string)
  
  return(ffdi.data)
}

# 3. List the files
ffdi.maps = list.files(file.path(cov.path, "masked", "ffdi95_sum"), pattern=".tif", full=T)

# 4. Run the extraction using the fixed function 'extract_ffdi95_data'
ffdi95_extract = lapply(ffdi.maps, FUN = extract_ffdi95_data, random_sample)
ffdi95_extract_long = do.call("rbind", ffdi95_extract)

message("--- FFDI 95th Percentile Head Check (Year should be 4 digits) ---")
print(head(ffdi95_extract_long))
message("---------------------------------------------------------------")




## Extract fire presences and absences. 
fire.maps = list.files(here("data", "fire_data"), pattern = ".tif", full=TRUE)

extract_fire_data = function(fire.ras.path, random.points){
  fire.ras = rast(fire.ras.path)
  fire.data = terra::extract(fire.ras, random.points, xy=TRUE)
  names(fire.data) <- c("ID", "burnt", "x","y")
  fire.data$year = readr::parse_number(names(fire.ras))
  return(fire.data)
}

burn_extract = lapply(fire.maps, FUN = extract_fire_data, random_sample)

burn_extract_long = do.call("rbind", burn_extract)


# Merge all extracted covariates into one dataset
joined_data <- burn_extract_long %>%
  left_join(ffdi_extract_long, by = c("ID", "x", "y", "year")) %>%
  left_join(ffdi95_extract_long, by = c("ID", "x", "y", "year")) %>%
  left_join(thunder_extract_long, by = c("ID", "x", "y", "year")) %>%
  left_join(spei12_extract_long, by = c("ID", "x", "y", "year")) %>%
  left_join(spei24_extract_long, by = c("ID", "x", "y", "year")) %>%
  left_join(kbdi_extract_long, by = c("ID", "x", "y", "year")) %>%
  # *** TSF ALIGNMENT FIX IN PLACE ***
  left_join(tsf_extract_aligned,   by = c("ID", "x", "y", "year")) %>% 
  left_join(cov_extracted, by = c("ID", "x", "y"))


write.csv(joined_data, output_file_path, row.names = FALSE)

# Save full dataset including all covariates 
summary(joined_data)