library(terra)
library(dplyr)
library(readr)
library(fs) 

# --- A. Define Paths and Constants ---
# Use the paths from your last successful setup
INPUT_CSV <- "/outputs/models/full_fire_covariates_no_folds.csv"
RASTER_BASE_DIR <- "/outputs/models/fire_data" 
OUTPUT_CSV_PATH <- "/outputs/final_duplicate_review.csv"

# --- B. Load Data and Filter Problematic Points ---

message("1. Loading data and filtering problematic points...")
data <- read_csv(INPUT_CSV, show_col_types = FALSE)

if (file.exists(OUTPUT_CSV_PATH)) {
  problem_points_df <- read_csv(OUTPUT_CSV_PATH, show_col_types = FALSE)
  message("Resuming review from partially completed CSV.")
} else {
  problem_points_df <- data %>%
    filter(tsf == 1, burnt == 1) %>%
    dplyr::select(ID, x, y, year, burnt, tsf) %>%
    mutate(is_duplicate = NA_integer_)
  message(paste("Starting fresh review for", nrow(problem_points_df), "points."))
}

# --- C. Utility Function for Plotting ONLY (No interactive input) ---

view_contradiction <- function(p) {
  
  year_before <- p$year - 1
  file_before <- path(RASTER_BASE_DIR, paste0("burned_area_", year_before, "_75m.tif"))
  file_of <- path(RASTER_BASE_DIR, paste0("burned_area_", p$year, "_75m.tif"))
  
  if (!file.exists(file_before) || !file.exists(file_of)) {
    warning(paste("Missing raster files for year", p$year, ". Skipping point ID:", p$ID))
    return(NA) 
  }
  
  # Load Rasters
  burn_rast_before <- rast(file_before)
  burn_rast_of <- rast(file_of)
  
  # Define Plotting Extent
  zoom_size <- 7500 
  plot_extent <- ext(p$x - zoom_size, p$x + zoom_size, p$y - zoom_size, p$y + zoom_size)
  
  # Set up plotting area
  par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
  burn_col_map <- c("0" = "lightgray", "1" = "darkred")
  
  # Plot 1: Year Before
  plot(burn_rast_before, main = paste("ID:", p$ID, " | Year:", year_before),
       col = burn_col_map, type = "classes", ext = plot_extent, legend = FALSE)
  points(p$x, p$y, pch = 3, col = "blue", cex = 1.5, lwd = 2) 
  
  # Plot 2: Year Of
  plot(burn_rast_of, main = paste("ID:", p$ID, " | Year:", p$year),
       col = burn_col_map, type = "classes", ext = plot_extent, legend = FALSE)
  points(p$x, p$y, pch = 3, col = "blue", cex = 1.5, lwd = 2) 
  
  # Return the index for external processing
  return(NA) # The function no longer handles input
}

# --- D. Iteration and Review Loop (External Input Required) ---

message("\n2. Starting manual review loop. Follow the prompts in the console.")

# Find the next unreviewed point
start_index <- which(is.na(problem_points_df$is_duplicate))[1]
if (is.na(start_index)) start_index <- 1

for (i in start_index:nrow(problem_points_df)) {
  
  p <- problem_points_df[i, ]
  
  # If the point has already been reviewed, skip it (shouldn't be needed with start_index but safer)
  if (!is.na(p$is_duplicate)) next
  
  # 1. Plot the rasters
  view_contradiction(p) 
  
  # 2. Print instructions and pause the script
  message(paste0("\nREVIEWING ID ", p$ID, " (", p$year - 1, " vs ", p$year, ")."))
  message("Action Required: Please type your result (1 for Duplicate, 0 for Not Duplicate, or 's' to skip) and press Enter.")
  
  # Use standard console input (a much safer alternative to readline inside the function)
  input <- readline(prompt="Enter 1, 0, or s: ")
  
  # 3. Process Input and Save
  if (tolower(input) == '1') {
    problem_points_df$is_duplicate[i] <- 1L
  } else if (tolower(input) == '0') {
    problem_points_df$is_duplicate[i] <- 0L
  } else {
    problem_points_df$is_duplicate[i] <- NA_integer_
    message("Point skipped.")
  }
  
  par(mfrow = c(1, 1)) # Reset plot layout after input
  
  # Save progress every 10 points
  if (i %% 10 == 0) {
    write_csv(problem_points_df, OUTPUT_CSV_PATH)
    message("Progress saved.")
  }
}

# --- E. Final Save ---
message("\n3. Review complete! Saving final result...")
write_csv(problem_points_df, OUTPUT_CSV_PATH)
message(paste("Final review CSV saved to:", OUTPUT_CSV_PATH))

# --- F. Next Step Guidance ---
message("The next step is to use the saved CSV to filter your main data and retrain your BRT model.")






#summary
library(dplyr)
library(readr)

# Define the path to the completed review file (REPLACE THIS LINE IF NEEDED)
REVIEW_CSV_PATH <- "C:/Users/charliehart/The University of Melbourne/Billy Geary - vic_fire_model/outputs/final_duplicate_review.csv"

# 1. Load the completed review file
review_results <- read_csv(REVIEW_CSV_PATH, show_col_types = FALSE)

# 2. Summarize duplicates by year
duplicates_by_year <- review_results %>%
  filter(is_duplicate == 1) %>%
  group_by(year) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  arrange(year)

# 3. Calculate total count
total_duplicates <- sum(duplicates_by_year$Count)

# 4. Print results
cat("--- Summary of Confirmed Duplicates (TSF=1 & Burnt=1) ---\n")
cat(paste("Total Duplicates Confirmed:", total_duplicates, "rows\n\n"))
cat("Distribution by Year (Chronological):\n")
print(duplicates_by_year)