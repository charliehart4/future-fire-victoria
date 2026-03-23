# ============================================================
# Script 11: Violin Plots of Dynamic Predictors (Baseline / Observed / Future)
# ============================================================

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
stack_dir <- "F:/vic_fire_mapping/output_data/prediction_stacks"
plot_dir  <- "F:/vic_fire_mapping/output_data/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Climate models
models <- c("ACCESS1-0", "GFDL-ESM2M", "CNRM-CM5", "MIROC5")

# Dynamic predictors only
dynamic_vars <- c("spei12_mean", "spei24_mean", "ffdi_95_days", "kbdi_95_days", "thunderstorm_days")

# ------------------------------------------------------------
# Helper function: sample values from stack and reshape
# ------------------------------------------------------------
sample_stack <- function(path, model_name, source_label, vars = dynamic_vars, n = 5000) {
  s <- rast(path)[[vars]]  # subset dynamic layers
  total_cells <- terra::ncell(s)
  
  # Random cell indices
  idx <- sample.int(total_cells, min(n, total_cells))
  
  # Extract values (no df argument)
  vals <- terra::extract(s, idx)
  
  # Drop the first column if it's the cell index (depends on terra version)
  if ("ID" %in% names(vals)) vals <- vals[ , -1, drop = FALSE]
  
  # Reshape and add metadata
  vals_long <- vals %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(model = model_name,
           source = source_label)
  
  return(vals_long)
}




# ------------------------------------------------------------
# 1️⃣  BASELINE (Modelled)
# ------------------------------------------------------------
baseline_list <- lapply(models, function(m) {
  f <- file.path(stack_dir, paste0("baseline_1985_2005_", m, "_RCP85_modelled_predstack.tif"))
  sample_stack(f, m, "Baseline (Modelled)")
})
df_baseline <- bind_rows(baseline_list)

# ------------------------------------------------------------
# 2️⃣  BASELINE OBSERVED (1985–2005)
# ------------------------------------------------------------
obs_file <- file.path(stack_dir, "baseline_1985_2005_observed_predstack.tif")
df_obs <- sample_stack(obs_file, "Observed", "Baseline (Observed)")

# ------------------------------------------------------------
# 3️⃣  FUTURE (Modelled)
# ------------------------------------------------------------
future_list <- lapply(models, function(m) {
  f <- file.path(stack_dir, paste0("future_2081_2099_", m, "_RCP85_modelled_predstack.tif"))
  sample_stack(f, m, "Future (Modelled)")
})
df_future <- bind_rows(future_list)

# ------------------------------------------------------------
# Helper function: make violin plot
# ------------------------------------------------------------
plot_violins <- function(df, title) {
  ggplot(df, aes(x = model, y = value, fill = model)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, color = "black", size = 0.2) +
    geom_boxplot(width = 0.12, fill = "white", color = "black", outlier.size = 0.5, outlier.alpha = 0.6) +
    facet_wrap(~variable, scales = "free_y", ncol = 3) +
    scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.9) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", size = 10)
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL
    )
}

# ------------------------------------------------------------
# Generate and save plots
# ------------------------------------------------------------
p_base <- plot_violins(df_baseline, "Dynamic Predictors — Baseline (1985–2005, Modelled)")
p_obs  <- plot_violins(df_obs, "Dynamic Predictors — Baseline (1985–2005, Observed)")
p_fut  <- plot_violins(df_future, "Dynamic Predictors — Future (2081–2099, Modelled)")

ggsave(file.path(plot_dir, "violin_dynamic_baseline_modelled.png"), p_base, width = 10, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "violin_dynamic_baseline_observed.png"), p_obs, width = 10, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "violin_dynamic_future_modelled.png"), p_fut, width = 10, height = 6, dpi = 300)

cat("\n✅ Violin plots saved to:\n", plot_dir, "\n")
