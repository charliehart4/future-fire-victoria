library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(patchwork)
library(ozmaps)
library(stringr)
library(dplyr)
library(ggspatial) 

# --- 1. SET PATHS ---
fd_path   <- "data/LF_DISTRICT.shp"
fire_path <- "data/vic_fire_mapping/fire_data"
out_path  <- "outputs/plots/study_area_map_FINAL.png"

# --- 2. VECTOR DATA: REGIONS & BACKGROUND ---
fdist <- st_read(fd_path, quiet = TRUE) |> st_make_valid()

# Dissolve into 6 regions
regions <- fdist |>
  group_by(REGIONNAME) |>
  summarise(.groups = "drop") |>
  mutate(REGIONNAME = str_to_title(REGIONNAME)) |>
  arrange(REGIONNAME) |> 
  mutate(ID = row_number())

region_labels <- st_point_on_surface(regions)
vic_bg <- ozmap_states |> filter(NAME == "Victoria")

# --- 3. RASTER DATA: SUM FIRE FREQUENCY (1980-2023) ---
all_files <- list.files(fire_path, pattern = "(?i)burned_area_.*\\.tif$", full.names = TRUE)
file_years <- as.numeric(str_extract(basename(all_files), "[0-9]{4}"))

# Filter for the modelling window
valid_files <- all_files[file_years >= 1980 & file_years <= 2023]

# Load, sum, and mask
fire_stack <- rast(valid_files)
fire_sum   <- sum(fire_stack, na.rm = TRUE)

# Match CRS and mask to regions only
raster_crs <- crs(fire_sum)
regions_v  <- vect(st_transform(regions, raster_crs))
fire_freq  <- crop(fire_sum, regions_v, mask = TRUE)


##map with codes rather than numbers
# --- 2. VECTOR DATA: REGIONS & BACKGROUND ---
fdist <- st_read(fd_path, quiet = TRUE) |> st_make_valid()

# Dissolve into 6 regions and add Abbreviations
regions <- fdist |>
  group_by(REGIONNAME) |>
  summarise(.groups = "drop") |>
  mutate(REGIONNAME = str_to_title(REGIONNAME)) |>
  arrange(REGIONNAME) |> 
  mutate(
    ID = row_number(),
    # Adding your requested codes based on the alphabetical order (1-6)
    CODE = case_match(ID,
                      1 ~ "BSW",
                      2 ~ "GI",
                      3 ~ "GR",
                      4 ~ "H",
                      5 ~ "LM",
                      6 ~ "PP")
  )

# Regenerate labels to include the new CODE column
region_labels <- st_point_on_surface(regions)
vic_bg <- ozmap_states |> filter(NAME == "Victoria")

##new plot
# --- 4. CREATE THE PLOTS ---

# A. CONTEXT MAP (AUSTRALIA)
plot_aus_context <- ggplot() +
  geom_sf(data = ozmap_states, fill = "#d1d1d1", color = "white", linewidth = 0.1) +
  geom_sf(data = st_transform(vic_bg, st_crs(ozmap_states)), fill = "#e63946", color = "white", linewidth = 0.1) + 
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5))

# B. MAIN STUDY AREA MAP
plot_vic_main <- ggplot() +
  geom_sf(data = st_transform(vic_bg, raster_crs), fill = "#f9f9f9", color = "grey80", linewidth = 0.3) +
  
  geom_spatraster(data = fire_freq) +
  scale_fill_viridis_c(
    option = "inferno", 
    direction = -1, 
    na.value = "transparent",
    name = "Times Burnt\n(1980-2023)",
    # FORCE INTEGERS: Set breaks at whole numbers (0 to 8 or max found)
    breaks = seq(0, 10, by = 2), 
    labels = seq(0, 10, by = 2)
  ) +
  
  geom_sf(data = st_transform(regions, raster_crs), fill = NA, color = "black", linewidth = 0.5) +
  
  # 15% LARGER REGION LABELS: Increased from 5.0 to 5.75
  geom_sf_text(data = st_transform(region_labels, raster_crs), aes(label = CODE), 
               size = 5.75, fontface = "bold", color = "black") +
  
  # LARGER SCALE BAR: Increased width_hint and text size
  annotation_scale(
    location = "br", 
    width_hint = 0.25, 
    unit_category = "metric",
    text_size = 5, # Enlarged text
    pad_x = unit(2.5, "cm"), 
    pad_y = unit(0.5, "cm")
  ) +
  
  # LARGER NORTH ARROW: Increased height/width from 0.8 to 0.95
  annotation_north_arrow(
    location = "br", 
    which_north = "true", 
    pad_x = unit(3.1, "cm"), 
    pad_y = unit(1.2, "cm"), 
    height = unit(0.95, "cm"), 
    width = unit(0.95, "cm"),
    style = north_arrow_orienteering(text_size = 0, fill = c("black", "white"))
  ) +
  
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.margin = margin(r = 40),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.height = unit(1.2, "cm")
  )

# --- 5. FINAL ASSEMBLY & SAVE ---
out_path_coded <- "outputs/plots/study_area_map_CODED.png"

# CONTEXT MAP: Shifted 'bottom' and 'top' up
final_map <- plot_vic_main + 
  inset_element(plot_aus_context, 
                left = 0.55,   
                bottom = 0.68, # Moved up from 0.60
                right = 0.95, 
                top = 1.0)     # Moved up to the very top

ggsave(out_path_coded, final_map, width = 12, height = 9, dpi = 300)

print(final_map)



##CLEAN MAP FOR POWERPOINT

# --- 1. SET PATHS ---
fd_path     <- "outputs/fire_risk_maps/Fire districts/LF_DISTRICT.shp"
fire_path   <- "data/fire_data"
out_path_hd <- "outputs/plots/study_area_map_HD_CLEAN.png"

# --- 2. VECTOR DATA: DISSOLVE REGIONS ---
fdist <- st_read(fd_path, quiet = TRUE) |> st_make_valid()
regions <- fdist |> group_by(REGIONNAME) |> summarise(.groups = "drop")
vic_bg <- ozmap_states |> filter(NAME == "Victoria")

# --- 3. RASTER DATA: SUM FIRE FREQUENCY ---
all_files <- list.files(fire_path, pattern = "(?i)burned_area_.*\\.tif$", full.names = TRUE)
file_years <- as.numeric(str_extract(basename(all_files), "[0-9]{4}"))
valid_files <- all_files[file_years >= 1980 & file_years <= 2023]

fire_stack <- rast(valid_files)
fire_sum   <- sum(fire_stack, na.rm = TRUE)

raster_crs <- crs(fire_sum)
regions_v  <- vect(st_transform(regions, raster_crs))
fire_freq  <- crop(fire_sum, regions_v, mask = TRUE)

# --- 4. CREATE THE PLOTS ---

# A. CONTEXT MAP (AUSTRALIA)
plot_aus_context <- ggplot() +
  geom_sf(data = ozmap_states, fill = "#d1d1d1", color = "white", linewidth = 0.05) +
  geom_sf(data = st_transform(vic_bg, st_crs(ozmap_states)), fill = "#e63946", color = "white", linewidth = 0.05) + 
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5))

# B. MAIN STUDY AREA MAP
plot_vic_main <- ggplot() +
  geom_sf(data = st_transform(vic_bg, raster_crs), fill = "#f9f9f9", color = "grey80", linewidth = 0.3) +
  geom_spatraster(data = fire_freq) +
  scale_fill_viridis_c(
    option = "inferno", 
    direction = -1, 
    na.value = "transparent",
    name = "Times Burnt\n(1980-2023)",
    breaks = seq(0, 10, by = 2), 
    labels = seq(0, 10, by = 2)
  ) +
  geom_sf(data = st_transform(regions, raster_crs), fill = NA, color = "black", linewidth = 0.6) +
  
  # ALIGNED SCALE BAR (Bottom-Left)
  annotation_scale(
    location = "bl", 
    width_hint = 0.3, 
    unit_category = "metric",
    text_size = 22, 
    height = unit(0.4, "cm"),
    pad_x = unit(1.0, "cm"),   # Positioned 1cm from the left
    pad_y = unit(0.8, "cm")
  ) +
  
  # NORTH ARROW (Now positioned immediately to the right of the bar)
  annotation_north_arrow(
    location = "bl",           # CHANGED: Anchored to the left
    which_north = "true", 
    pad_x = unit(8.0, "cm"),   # ADJUST THIS: Pushes the arrow just past the end of the bar
    pad_y = unit(0.8, "cm"),   # MATCHED: Level with the bar
    height = unit(1.2, "cm"), 
    width = unit(1.2, "cm"),
    style = north_arrow_orienteering(text_size = 0, fill = c("black", "white"))
  ) +
  
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.margin = margin(r = 50),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.height = unit(1.5, "cm")
  )

# --- 5. FINAL ASSEMBLY & SAVE ---
final_map_hd <- plot_vic_main + 
  inset_element(plot_aus_context, 
                left = 0.55, bottom = 0.72, right = 0.95, top = 1.0)

ggsave(out_path_hd, final_map_hd, width = 14, height = 10, dpi = 600, scale=.8)


print(final_map_hd)

#HD PDF 
out_path_pdf <- "outputs/plots/study_area_map_HD_CLEAN.pdf"
ggsave(out_path_pdf, final_map_hd, width = 14, height = 10, device = "pdf", scale = 0.8)




##part 2. temporal summary
# --- 6. OPTIMISED TEMPORAL EXTRACTION ---
cat("Starting optimised extraction of annual hectares burned (1980-2023)...\n")

# Transform regions once to match the raster CRS
first_r <- rast(valid_files[1])
regions_vect <- vect(st_transform(regions, crs(first_r)))

stats_list <- list()

for(i in seq_along(valid_files)){
  year_val <- file_years[file_years >= 1980 & file_years <= 2023][i]
  
  # Load the annual raster
  r <- rast(valid_files[i])
  
  # Optimization: Crop the raster to the regions' extent to ignore the rest of the state
  r_cropped <- crop(r, regions_vect)
  
  # Extract sum of pixels (1 = burnt) per region
  # zonal() returns a dataframe with the sums
  ex <- terra::zonal(r_cropped, regions_vect, fun = "sum", na.rm = TRUE)
  
  # Calculate Hectares: (Count * 5625m2) / 10000m2
  # ex[,1] contains the sums for the burnt area layer
  hectares <- (ex[, 1] * 5625) / 10000
  
  stats_list[[i]] <- data.frame(
    Year = year_val,
    Region = regions$REGIONNAME,
    ID = regions$ID,
    Hectares = hectares
  )
  
  if(i %% 5 == 0) cat(paste0("Processed year: ", year_val, "\n"))
}

# Combine and Save
burn_history <- bind_rows(stats_list)

csv_path <- "data/annual_burn_stats.csv"
write.csv(burn_history, csv_path, row.names = FALSE)

cat("✅ Data extracted and saved to CSV.\n")




# --- 7. CREATE CUMULATIVE STACKED AREA CHART ---

# 1. Prepare and reorder the data
burn_history_cumul <- burn_history |>
  mutate(
    CODE = case_match(as.numeric(as.character(ID)),
                      1 ~ "BSW",
                      2 ~ "GI",
                      3 ~ "GR",
                      4 ~ "H",
                      5 ~ "LM",
                      6 ~ "PP")
  ) |>
  mutate(CODE = factor(CODE, levels = c("BSW", "GR", "LM", "PP", "H", "GI"))) |>
  arrange(CODE, Year) |>
  group_by(CODE) |>
  mutate(Cumulative_Hectares = cumsum(Hectares)) |>
  ungroup()

# 2. Plotting
plot_temporal_cumul <- ggplot(burn_history_cumul, aes(x = Year, y = Cumulative_Hectares, fill = CODE)) +
  geom_area(alpha = 0.85, color = "white", linewidth = 0.1) +
  scale_fill_viridis_d(
    option = "inferno", 
    direction = -1, 
    name = "Region",
    guide = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 1), # Makes the legend colors solid/easier to see
      keywidth = unit(1.5, "cm")      # Makes the legend color boxes wider
    )
  ) +
  scale_x_continuous(breaks = seq(1980, 2023, by = 5), expand = c(0,0)) +
  scale_y_continuous(labels = scales::label_comma(), expand = c(0,0)) +
  # Title removed as requested
  labs(
    x = "Year",
    y = "Cumulative Hectares Burned"
  ) +
  theme_minimal(base_size = 14) + # Boosts all text sizes globally
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    # Emphasize Axis Titles
    axis.title = element_text(face = "bold", size = 16), 
    axis.text = element_text(size = 12, color = "black"),
    # Legend Text size
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    # Adjust plot margins to make the actual plot appear slightly smaller/more framed
    plot.margin = margin(t = 20, r = 40, b = 10, l = 20),
    legend.box.spacing = unit(1, "cm")
  )

# Save the plot
plot_out_cumul <- "outputs/plots/cumulative_burn_chart_FINAL.png"
ggsave(plot_out_cumul, plot_temporal_cumul, width = 10, height = 7, dpi = 300, bg = "white")

print(plot_temporal_cumul)


head(burn_history_cumul)





# --- 7b. CREATE FACETED LINE GRAPH PER REGION ---

plot_regional_lines <- ggplot(burn_history_cumul, aes(x = Year, y = Hectares, color = CODE)) +
  # Using a thicker line for the cumulative trend
  geom_line(linewidth = 1.2) +
  # Removing scales = "free_y" forces the same scale across all regions
  facet_wrap(~CODE, nrow = 2) + 
  scale_color_viridis_d(option = "inferno", direction = -1) +
  scale_x_continuous(breaks = seq(1980, 2023, by = 10)) +
  # Fixed scale with comma formatting
  scale_y_continuous(labels = scales::label_comma(), expand = c(0, 0.1)) +
  labs(
    title = NULL,
    x = "Year",
    y = "Area Burned (Ha)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey90", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

# Save the faceted plot
line_out_cumul <- "outputs/plots/regional_lines_FIXED.png"
ggsave(line_out_cumul, plot_regional_lines, width = 12, height = 7, dpi = 300, bg = "white")

print(plot_regional_lines)




# --- 7c. LOOP TO SAVE ALL REGIONAL PLOTS (MAXIMUM TEXT SIZE - LINE ONLY) ---

region_list <- levels(burn_history_cumul$CODE)

for (reg in region_list) {
  # 1. Filter data for specific region
  reg_data <- burn_history_cumul %>% filter(CODE == reg)
  
  # 2. Create plot
  p <- ggplot(reg_data, aes(x = Year, y = Hectares)) +
    # Points removed, only the line remains
    geom_line(color = "black", linewidth = 2) + 
    
    facet_wrap(~CODE) + 
    
    scale_x_continuous(breaks = seq(1980, 2023, by = 10)) +
    scale_y_continuous(
      labels = scales::label_comma(), 
      expand = expansion(mult = c(0, 0.25)) 
    ) +
    
    labs(x = "Year", y = "Area Burned (Ha)") +
    
    # --- GLOBAL SIZE BOOST ---
    theme_minimal(base_size = 32) + 
    theme(
      # --- BANNER (STRIP) - GIANT ---
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 2),
      strip.text = element_text(face = "bold", size = 45, margin = margin(t=10, b=10)),
      
      # --- AXIS TITLES ---
      axis.title.x = element_text(face = "bold", margin = margin(t = 20)),
      axis.title.y = element_text(face = "bold", margin = margin(r = 20)),
      
      # --- AXIS TICK LABELS ---
      axis.text = element_text(color = "black", face = "bold", size = 28),
      axis.text.x = element_text(angle = 45, hjust = 1),
      
      # --- PANEL & BORDER ---
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(30, 30, 30, 30) 
    )
  
  # 3. Save individual file
  file_name <- paste0("annual_fire_CLEAN_LINE_", reg, ".png")
  ggsave(
    filename = file.path(out_dir_plots, file_name), 
    plot = p, 
    width = 9, 
    height = 7, 
    dpi = 600
  )
}

message("Individual plots saved with clean lines and giant labels.")



#y = not free
# --- 1. CALCULATE GLOBAL Y-AXIS MAXIMUM ---
# Find the highest annual hectares recorded in any region to set the limit
global_max_h <- max(burn_history_cumul$Hectares, na.rm = TRUE)

# --- 2. LOOP TO SAVE ALL REGIONAL PLOTS (STANDARDISED Y-AXIS) ---
region_list <- levels(burn_history_cumul$CODE)

for (reg in region_list) {
  # Filter data for specific region
  reg_data <- burn_history_cumul %>% filter(CODE == reg)
  
  # Create plot
  p <- ggplot(reg_data, aes(x = Year, y = Hectares)) +
    geom_line(color = "black", linewidth = 2) + 
    
    facet_wrap(~CODE) + 
    
    scale_x_continuous(breaks = seq(1980, 2023, by = 10)) +
    
    # Standardised Y-axis using the global maximum
    scale_y_continuous(
      labels = scales::label_comma(),
      limits = c(0, global_max_h),
      expand = expansion(mult = c(0, 0.25)) # Room for the banner
    ) +
    
    labs(x = "Year", y = "Area Burned (Ha)") +
    
    # --- GLOBAL SIZE BOOST ---
    theme_minimal(base_size = 32) + 
    theme(
      # --- BANNER (STRIP) - GIANT ---
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 2),
      strip.text = element_text(face = "bold", size = 45, margin = margin(t=10, b=10)),
      
      # --- AXIS TITLES ---
      axis.title.x = element_text(face = "bold", margin = margin(t = 20)),
      axis.title.y = element_text(face = "bold", margin = margin(r = 20)),
      
      # --- AXIS TICK LABELS ---
      axis.text = element_text(color = "black", face = "bold", size = 28),
      axis.text.x = element_text(angle = 45, hjust = 1),
      
      # --- PANEL & BORDER ---
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(30, 30, 30, 30) 
    )
  
  # 3. Save individual file
  file_name <- paste0("annual_fire_STANDARDISED_", reg, ".png")
  ggsave(
    filename = file.path(out_dir_plots, file_name), 
    plot = p, 
    width = 9, 
    height = 7, 
    dpi = 600
  )
}

message("Individual plots saved with a standardised Y-axis.")



##log plots
# --- 1. CALCULATE TRUE DYNAMIC LIMITS ---
# Find the actual maximum hectare value across all regions
actual_max <- max(burn_history_cumul$Hectares, na.rm = TRUE)

# Set the ceiling to either 150,000 or slightly above the actual max
# This ensures even massive years (like 200k+ Ha) are captured
global_log_max <- max(150000, actual_max * 1.2)

# --- 2. LOOP TO SAVE ALL REGIONAL PLOTS ---
region_list <- levels(burn_history_cumul$CODE)

for (reg in region_list) {
  # Filter data for specific region
  reg_data <- burn_history_cumul %>% filter(CODE == reg)
  
  # Create plot
  p <- ggplot(reg_data, aes(x = Year, y = Hectares + 1)) +
    geom_line(color = "black", linewidth = 2.2) + 
    facet_wrap(~CODE) + 
    
    scale_x_continuous(breaks = seq(1980, 2023, by = 10)) +
    
    scale_y_log10(
      # Standard breaks + the next power of 10 if necessary
      breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
      labels = c("0", "10", "100", "1,000", "10,000", "100,000", "1,000,000"),
      limits = c(1, global_log_max),
      expand = expansion(mult = c(0, 0.2)) 
    ) +
    
    labs(x = "Year", y = "Area Burned (Hectares)") +
    
    theme_minimal(base_size = 32) + 
    theme(
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 2),
      strip.text = element_text(face = "bold", size = 45, margin = margin(t=10, b=10)),
      axis.title.x = element_text(face = "bold", margin = margin(t = 20)),
      axis.title.y = element_text(face = "bold", margin = margin(r = 20), size = 26),
      axis.text = element_text(color = "black", face = "bold", size = 28),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(30, 30, 30, 30) 
    )
  
  # Save individual file
  file_name <- paste0("annual_fire_LOG_FINAL_", reg, ".png")
  ggsave(
    filename = file.path(out_dir_plots, file_name), 
    plot = p, 
    width = 9, 
    height = 8, 
    dpi = 600,
    scale=0.9
  )
}

message("Success: All regional log-plots saved with dynamic limits. Actual max found was: ", actual_max)












##PATCHWORK PLOT

# 1. Prepare the Map (Panel A)
# Setting all margins to 0 and using an aggressive negative bottom margin
p1 <- final_map + 
  labs(tag = NULL) +
  theme(
    plot.margin = margin(t = 0, r = 0, b = -80, l = 0), # Pulls the plot up significantly
    # Ensure no padding within the panel itself
    panel.spacing = unit(0, "cm")
  )

# 2. Prepare the Cumulative Plot (Panel B)
p2 <- plot_temporal_cumul + 
  labs(tag = NULL) +
  theme(
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )

# 3. Combine with a dominant Map height
# We use 'guides = "collect"' just in case, and 'heights' to maintain the map's dominance
combined_figure <- (p1 / p2) + 
  plot_layout(
    heights = c(3, 1),
    guides = "keep"
  ) & 
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    # This removes the space between the patchwork panels
    panel.spacing = unit(0, "lines") 
  )

# 4. Save with high resolution
final_out_path <- "outputs/plots/Fig1_Final_MAX_MAP.png"

# Shifting to a 12x16 aspect ratio to tighten the overall look
ggsave(
  final_out_path, 
  combined_figure, 
  width = 12, 
  height = 16, 
  dpi = 300, 
  bg = "white"
)

# Show it
print(combined_figure)






