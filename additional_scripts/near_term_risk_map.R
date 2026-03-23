### Near-term Fire Risk Forecast Mapping
library(terra)
library(ggplot2)
library(tidyterra)
library(viridis)
library(dplyr)
library(sf)

# === 1. Load Data ===
pred_path <- "outputs/predictions_nt/fire_prediction_observed_2021to2024.tif"
r_nt <- rast(pred_path)

# Native Vegetation Raster for masking
veg_path <- "dataNVR2017_CONDITION.tif"
r_veg <- rast(veg_path)

# === 2. Processing & Masking ===

# Align veg raster to prediction geometry
r_veg_proj <- project(r_veg, r_nt, method = "near")
r_veg_mask <- r_veg_proj > 0

# Mask the prediction raster to native vegetation
r_nt_final <- mask(r_nt, r_veg_mask, maskvalues = 0)

# === 3. Plotting (Continuous Scale) ===

p_nt_continuous <- ggplot() +
  # Use geom_spatraster for the continuous variable
  geom_spatraster(data = r_nt_final) +
  scale_fill_viridis_c(
    option = "inferno",
    name = "Pr(Fire)",
    na.value = "transparent",
    # Adjust limits if you want to standardise across different maps (e.g., 0 to 0.5)
    limits = c(0, max(values(r_nt_final), na.rm = TRUE)),
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black",
      barwidth = 1.2,
      barheight = 15
    )
  ) +
  # Clean theme with no title as requested
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_nt_continuous)

# === 3. Summary Statistics ===
# Extract values while removing NAs (non-veg and missing data)
stats_vals <- values(r_nt_final, na.rm = TRUE)

summary_table <- tibble(
  Variable = "Near-term Pr(Fire) 2021-2024",
  Min = min(stats_vals),
  Max = max(stats_vals),
  Mean = mean(stats_vals),
  SD = sd(stats_vals)
)

cat("\n--- Summary Statistics for Near-term Forecast ---\n")
print(summary_table)

# === 4. Export (Optional) ===
# ggsave("NearTerm_Forecast_Continuous.png", p_nt_continuous, width = 10, height = 8, dpi = 300, bg = "white")