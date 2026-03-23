
# Visualise the model outputs
library(dismo)
library(gbm)

# RAC MODEL
rac.model = readRDS("outputs/brt_full_rac_lr0.005_lr0.005_ffdi_95.rds")

summary(rac.model)

gbm.interactions(rac.model)

gbm.perspec(rac.model, x = "rac", y = "spei12_mean", z.range = c(0,0.03))

##cur
plot.gbm(rac.model, "rac", type = "response")
plot.gbm(rac.model, "fuel_management_zones", type = "response") # THis looks good
plot.gbm(rac.model, "bio5", type = "response") # THis looks good
plot.gbm(rac.model, "bio18", type = "response") # THis looks good
plot.gbm(rac.model, "broad_refuges", type = "response") # THis looks good
plot.gbm(rac.model, "local_refuges", type = "response") # THis looks good
plot.gbm(rac.model, "spei12_mean", type = "response", smooth=TRUE) # high parts in right spot but a bit funky
plot.gbm(rac.model, "ffdi_95_days", type = "response", smooth=TRUE) # High parts in right spot but a bit funky
plot.gbm(rac.model, "spei24_mean", type = "response", smooth=TRUE) # High parts in two ends for 24 - fuel?
plot.gbm(rac.model, "thunderstorm_days", type = "response", smooth=TRUE) # High parts in right spot
plot.gbm(rac.model, "tsf", type = "response", smooth=TRUE) # High parts in right spot

plot.gbm(rac.model, c("fuel_management_zones", "ffdi_95_days"), type = "response") # THis looks good
plot.gbm(rac.model, c("fuel_management_zones", "spei12_mean"), type = "response") # THis looks good
plot.gbm(rac.model, c("fuel_management_zones", "spei24_mean"), type = "response") # THis looks good
plot.gbm(rac.model, c("fuel_management_zones", "thunderstorm_days"), type = "response") # THis looks good

plot.gbm(rac.model, c("spei12_mean", "ffdi_95_days"), type = "response") # THis looks good
plot.gbm(rac.model, c("spei24_mean", "ffdi_95_days"), type = "response") # THis looks good
plot.gbm(rac.model, c("rac", "ffdi_95_days"), type = "response") # THis looks good

# No RAC MODEL
model = readRDS("outputs/brt_base_lr0.005_ffdi_95.RDS")
summary(model)







# PDP interaction plots ALIGNED 
library(gbm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales) 

# --- 1. SETUP & CONFIGURATION ---
output_dir <- "outputs/plots"

var_labels <- list(
  ffdi_95_days          = "Days FFDI >95th percentile",
  spei12_mean           = "SPEI (12-month)",
  thunderstorm_days     = "Thunderstorm Days",
  kbdi_95_days          = "Days KBDI >95th percentile",
  spei24_mean           = "SPEI (24-month)",
  rac                   = "Residual Autocovariate",
  fuel_management_zones = "Fuel Management Zone",
  distance_roads        = "Distance from Roads"
)

# --- 2. CORE FUNCTIONS ---

get_gbm_2d <- function(model, vars) {
  pd <- plot.gbm(model, i.var = vars, type = "response", return.grid = TRUE)
  resp_col <- setdiff(names(pd), vars)[1]
  pd <- pd[, c(vars, resp_col)]
  names(pd) <- c("x", "y", "z")
  return(pd)
}

make_2d_plot <- function(df, xlab, ylab) {
  p <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma") +
    labs(x = xlab, y = ylab, fill = "Pr(Fire)") +
    theme_minimal(base_size = 28) + 
    theme(
      aspect.ratio = 1,
      panel.border = element_blank(), 
      panel.grid.major = element_line(color = "grey92"), 
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),    
      axis.title = element_text(size = 24, face = "plain"), 
      axis.text = element_text(size = 18), 
      legend.title = element_text(size = 22), 
      legend.text = element_text(size = 18),  
      legend.key.height = unit(1.2, "cm"),    
      plot.margin = margin(20, 20, 20, 20) 
    )
  
  # Targeted Logic for X-axis
  if(is.numeric(df$x)) {
    if(xlab == "Residual Autocovariate") {
      # RAC gets fewer breaks and rounded decimals
      p <- p + scale_x_continuous(n.breaks = 4, labels = label_number(accuracy = 0.01))
    } else {
      # Everything else stays default with overlap check
      p <- p + scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
    }
  }
  
  # Targeted Logic for Y-axis
  if(is.numeric(df$y)) {
    if(ylab == "Residual Autocovariate") {
      p <- p + scale_y_continuous(n.breaks = 4, labels = label_number(accuracy = 0.01))
    } else {
      p <- p + scale_y_continuous(guide = guide_axis(check.overlap = TRUE))
    }
  }
  
  return(p)
}

build_interaction_plots <- function(model, interaction_list) {
  lapply(interaction_list, function(vars) {
    df <- get_gbm_2d(model, vars)
    xlab <- var_labels[[vars[1]]]
    ylab <- var_labels[[vars[2]]]
    make_2d_plot(df, xlab, ylab)
  })
}

# --- 3. BASE MODEL WORKFLOW ---
top_interactions_base <- list(
  c("ffdi_95_days", "spei12_mean"),
  c("thunderstorm_days", "ffdi_95_days"),
  c("thunderstorm_days", "kbdi_95_days"),
  c("kbdi_95_days", "spei24_mean"),
  c("spei24_mean", "ffdi_95_days"),
  c("thunderstorm_days", "spei12_mean")
)

plots_base <- build_interaction_plots(model, top_interactions_base)
interaction_grid_base <- wrap_plots(plots_base, ncol = 3) & theme(legend.position = "right")

print(interaction_grid_base)

ggsave(filename = file.path(output_dir, "interaction_pdp_top6_BASE_final.pdf"),
       plot = interaction_grid_base, width = 20, height = 15, units = "in")

ggsave(filename = file.path(output_dir, "interaction_pdp_top6_BASE_final.png"),
       plot = interaction_grid_base, width = 20, height = 15, units = "in", dpi = 300)

# --- 4. RAC MODEL WORKFLOW ---
top_interactions_rac <- list(
  c("rac", "thunderstorm_days"),
  c("rac", "spei24_mean"),
  c("rac", "fuel_management_zones"),
  c("thunderstorm_days", "ffdi_95_days"),
  c("rac", "kbdi_95_days"),
  c("rac", "ffdi_95_days"),
  c("thunderstorm_days", "spei24_mean"),
  c("rac", "distance_roads")
)

plots_rac <- build_interaction_plots(rac.model, top_interactions_rac)
interaction_grid_rac <- wrap_plots(plots_rac, ncol = 3) & theme(legend.position = "right")

print(interaction_grid_rac)

ggsave(filename = file.path(output_dir, "interaction_pdp_top8_RAC_final.pdf"),
       plot = interaction_grid_rac, width = 20, height = 22, units = "in")

ggsave(filename = file.path(output_dir, "interaction_pdp_top8_RAC_final.png"),
       plot = interaction_grid_rac, width = 20, height = 22, units = "in", dpi = 300)






##Standarised scale

# --- 1. GLOBAL LIMIT CALCULATION ---

# Function to extract all interaction data to find the global Z range
get_global_range <- function(model, interaction_list) {
  all_responses <- sapply(interaction_list, function(vars) {
    df <- get_gbm_2d(model, vars)
    return(df$z)
  })
  return(range(all_responses, na.rm = TRUE))
}

# Calculate ranges for both models
z_range_base <- get_global_range(model, top_interactions_base)
z_range_rac <- get_global_range(rac.model, top_interactions_rac)

# --- 2. UPDATED PLOTTING FUNCTION (Standard Linear) ---

make_2d_plot_linear <- function(df, xlab, ylab, z_limits) {
  p <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    # Standard Linear Scale with fixed limits
    scale_fill_viridis_c(
      option = "magma",
      limits = z_limits,
      labels = label_number(accuracy = 0.001),
      name = "Pr(Fire)"
    ) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 28) + 
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_line(color = "grey92"), 
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 22), 
      axis.text = element_text(size = 16), 
      legend.title = element_text(size = 22), 
      legend.text = element_text(size = 18),  
      legend.key.height = unit(1.5, "cm"),    
      plot.margin = margin(10, 10, 10, 10) 
    )
  
  # Handle X-axis specifically for categorical variables (FMZ)
  if(!is.numeric(df$x)) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))
  } else {
    p <- p + scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
  }
  
  if(is.numeric(df$y)) {
    p <- p + scale_y_continuous(guide = guide_axis(check.overlap = TRUE))
  }
  
  return(p)
}

# --- 3. ASSEMBLY ---

# BASE MODEL
plots_base <- lapply(top_interactions_base, function(vars) {
  df <- get_gbm_2d(model, vars)
  make_2d_plot_linear(df, var_labels[[vars[1]]], var_labels[[vars[2]]], z_range_base)
})

interaction_grid_base <- wrap_plots(plots_base, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# RAC MODEL
plots_rac <- lapply(top_interactions_rac, function(vars) {
  df <- get_gbm_2d(rac.model, vars)
  make_2d_plot_linear(df, var_labels[[vars[1]]], var_labels[[vars[2]]], z_range_rac)
})

interaction_grid_rac <- wrap_plots(plots_rac, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# --- 4. SAVE ---
ggsave(file.path(output_dir, "interaction_BASE_standard_scale.png"), interaction_grid_base, width = 20, height = 15, dpi = 300)
ggsave(file.path(output_dir, "interaction_RAC_standard_scale.png"), interaction_grid_rac, width = 20, height = 22, dpi = 300)





## logit scale interaction plots - Final Polish
# --- 1. SET HARD-CODED LIMITS & BREAKS ---

# Base Model: 0.0004 to 0.779
z_lims_base   <- c(0.0004, 0.8)
z_breaks_base <- c(0.001, 0.01, 0.1, 0.4, 0.8)

# RAC Model: 0.0003 to 0.989 (Added 0.7 intermediate break)
z_lims_rac    <- c(0.0003, 0.99)
z_breaks_rac  <- c(0.001, 0.01, 0.1, 0.5, 0.7, 0.99)

# --- 2. UPDATED PLOTTING FUNCTION ---

make_interaction_plot_final <- function(df, xlab, ylab, limits, breaks) {
  
  p <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c(
      option = "magma",
      trans = "logit",
      limits = limits,
      breaks = breaks,
      labels = label_number(accuracy = 0.001),
      name = "logit(Pr(Fire))"
    ) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 28) + 
    theme(
      aspect.ratio = 1,
      legend.key.height = unit(2.5, "cm"), 
      # Larger title and increased bottom margin for spacing
      legend.title = element_text(size = 24, face = "plain", margin = margin(b = 20)),
      legend.text = element_text(size = 18),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 16),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Handle X-axis for Fuel Management Zones
  if(!is.numeric(df$x)) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))
  } else {
    p <- p + scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
  }
  
  return(p)
}

# --- 3. ASSEMBLY ---

# Base Model
plots_base <- lapply(top_interactions_base, function(vars) {
  df <- get_gbm_2d(model, vars)
  make_interaction_plot_final(df, var_labels[[vars[1]]], var_labels[[vars[2]]], z_lims_base, z_breaks_base)
})

final_base <- wrap_plots(plots_base, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# RAC Model
plots_rac <- lapply(top_interactions_rac, function(vars) {
  df <- get_gbm_2d(rac.model, vars)
  make_interaction_plot_final(df, var_labels[[vars[1]]], var_labels[[vars[2]]], z_lims_rac, z_breaks_rac)
})

final_rac <- wrap_plots(plots_rac, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")


# Base Model Exports
ggsave(file.path(output_dir, "interaction_BASE_logit_final_polished.png"), 
       final_base, width = 20, height = 15, dpi = 300)

ggsave(file.path(output_dir, "interaction_BASE_logit_final_polished.pdf"), 
       final_base, width = 20, height = 15, device = cairo_pdf)

# RAC Model Exports
ggsave(file.path(output_dir, "interaction_RAC_logit_final_polished.png"), 
       final_rac, width = 20, height = 22, dpi = 300)

ggsave(file.path(output_dir, "interaction_RAC_logit_final_polished.pdf"), 
       final_rac, width = 20, height = 22, device = cairo_pdf)



## pdp facet plot for report
library(ggplot2)
library(gbm)

# Simpler method using plot.gbm directly
get_gbm_plot_data <- function(gbm_model, var_name) {
  
  # Use plot.gbm to generate the plot data
  plot_data <- plot.gbm(gbm_model, i.var = var_name, 
                        return.grid = TRUE, type = "response")
  
  # Extract the data
  data.frame(
    variable = plot_data[[var_name]],  # Use variable name as column name
    response = plot_data$y
  )
}

# Create Response Plots
ffdi_plot <- get_gbm_plot_data(model, "ffdi_95_days") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Days FFDI>95th Percentile", 
       y = "Pr(Fire)") +
  theme_minimal()

spei12_plot <- get_gbm_plot_data(model, "spei12_mean") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Mean SPEI (12 Months)", 
       y = "Pr(Fire)") +
  theme_minimal()

spei24_plot <- get_gbm_plot_data(model, "spei24_mean") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Mean SPEI (24 Months)", 
       y = "Pr(Fire)") +
  theme_minimal()

bio5_plot <- get_gbm_plot_data(model, "bio5") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Max Temperature of Warmest Month", 
       y = "Pr(Fire)") +
  theme_minimal()

kbdi_plot <- get_gbm_plot_data(model, "kbdi_95_days") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Days KBDI >95th Percentile", 
       y = "Pr(Fire)") +
  theme_minimal()


fmz_plot <- get_gbm_plot_data(model, "fuel_management_zones") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_bar(color = "black", size = 1, stat="identity") +
  labs(x = "Fuel Management Zone", 
       y = "Pr(Fire)") +
  theme_minimal()

thunderstorm_plot <- get_gbm_plot_data(model, "thunderstorm_days") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Days of Thunderstorm Occurrence", 
       y = "Pr(Fire)") +
  theme_minimal()

bio18_plot <- get_gbm_plot_data(model, "bio18") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Precipitation of Warmest Quarter", 
       y = "Pr(Fire)") +
  theme_minimal()

twi_plot <- get_gbm_plot_data(model, "twi") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Topographic Wetness Index", 
       y = "Pr(Fire)") +
  theme_minimal()

broad_refuges_plot <- get_gbm_plot_data(model, "broad_refuges") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Landscape Refuges", 
       y = "Pr(Fire)") +
  theme_minimal()

local_refuges_plot <- get_gbm_plot_data(model, "local_refuges") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Local Refuges", 
       y = "Pr(Fire)") +
  theme_minimal()

distance_roads_plot <- get_gbm_plot_data(model, "distance_roads") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Distance from Roads", 
       y = "Pr(Fire)") +
  theme_minimal()


#theme_set(theme_minimal(base_size = 16))

final_plot <-
  ffdi_plot + spei12_plot + spei24_plot + bio5_plot + kbdi_plot + fmz_plot +
  thunderstorm_plot + bio18_plot + twi_plot + broad_refuges_plot +
  distance_roads_plot + local_refuges_plot +
  plot_layout(axes = "collect", ncol = 3) &
  theme(
    axis.title.x = element_text(
      size = 22,
      margin = margin(t = 10)   
    ),
    axis.title.y = element_text(
      size = 24,
      margin = margin(r = 10)   
    ),
    axis.text = element_text(size = 23)
  )
ggsave(
  filename = file.path(output_dir, "pdp_all_covariates_highres.png"),
  plot = final_plot,
  width = 18,
  height = 15,
  dpi = 600,
  units = "in"
)
ggsave(
  filename = file.path(output_dir, "pdp_all_covariates_highres.pdf"),
  plot = final_plot,
  width = 18,
  height = 15,
  units = "in",
  device = cairo_pdf
)


## BASE MODEL pdp facet plot - Version 3.1 (Logit Scale + Fixed FMZ + Expanded Labels)
library(ggplot2)
library(gbm)
library(dplyr)
library(patchwork)
library(scales) 

# 1. Load data and create a sample for the rug
raw_data <- read.csv("outputs/fire_modelling_rac_ready_forCV.csv")
set.seed(123)
rug_sample <- raw_data[sample(nrow(raw_data), 1000), ] 

# 2. Extract plot data and find global Y-axis limits
vars_to_plot <- c("ffdi_95_days", "spei12_mean", "spei24_mean", "bio5", "kbdi_95_days", 
                  "fuel_management_zones", "thunderstorm_days", "bio18", "twi", 
                  "broad_refuges", "local_refuges", "distance_roads")

all_pdp_data <- lapply(vars_to_plot, function(v) get_gbm_plot_data(model, v))

# Standardise Y limits across all plots
y_min <- min(sapply(all_pdp_data, function(d) min(d$response)))
y_max <- max(sapply(all_pdp_data, function(d) max(d$response)))
y_limits <- c(max(0.0001, y_min), min(0.9999, y_max))

# 3. Create Response Plots
# Expanded breaks for better readability on logit scale
logit_breaks <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)

apply_logit_style <- function(p, x_lab) {
  p + 
    scale_y_continuous(
      trans = "logit", 
      limits = y_limits,
      breaks = logit_breaks,
      labels = label_number(accuracy = 0.001) # Forces decimal display
    ) +
    labs(x = x_lab, y = "logit(Pr(Fire))") +
    theme_minimal()
}

# FFDI
ffdi_plot <- get_gbm_plot_data(model, "ffdi_95_days") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days FFDI > 95th Percentile") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = ffdi_95_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# SPEI 12
spei12_plot <- get_gbm_plot_data(model, "spei12_mean") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Mean SPEI (12 Months)") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = spei12_mean), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# SPEI 24
spei24_plot <- get_gbm_plot_data(model, "spei24_mean") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Mean SPEI (24 Months)") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = spei24_mean), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# BIO5
bio5_plot <- get_gbm_plot_data(model, "bio5") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Max Temp Warmest Month") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = bio5), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# KBDI
kbdi_plot <- get_gbm_plot_data(model, "kbdi_95_days") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days KBDI > 95th Percentile") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = kbdi_95_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# FIXED FMZ Plot 
fmz_plot <- fmz_data %>%
  mutate(variable = as.factor(variable)) %>% 
  ggplot(aes(x = variable, y = response)) +
  geom_segment(aes(x=variable, xend=variable, y=y_limits[1], yend=response), color="grey") +
  geom_point(size=4, color="black") +
  scale_y_continuous(
    trans = "logit", 
    limits = y_limits, 
    breaks = logit_breaks,
    labels = label_number(accuracy = 0.001)
  ) +
  labs(x = "Fuel Management Zone", y = "Pr(Fire)") +
  theme_minimal()
fmz_plot

# Thunderstorm
thunderstorm_plot <- get_gbm_plot_data(model, "thunderstorm_days") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days of Thunderstorm") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = thunderstorm_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# BIO18
bio18_plot <- get_gbm_plot_data(model, "bio18") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Precip Warmest Quarter") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = bio18), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# TWI
twi_plot <- get_gbm_plot_data(model, "twi") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Topographic Wetness Index") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = twi), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# Broad Refuges
broad_refuges_plot <- get_gbm_plot_data(model, "broad_refuges") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Landscape Refuges") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = broad_refuges), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# Distance Roads
distance_roads_plot <- get_gbm_plot_data(model, "distance_roads") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Distance from Roads") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = distance_roads), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# Local Refuges
local_refuges_plot <- get_gbm_plot_data(model, "local_refuges") %>%
  ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Local Refuges") +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = local_refuges), inherit.aes = FALSE, alpha = 0.3, sides = "b")


final_plot <-
  ffdi_plot + spei12_plot + spei24_plot + bio5_plot + kbdi_plot + fmz_plot +
  thunderstorm_plot + bio18_plot + twi_plot + broad_refuges_plot +
  distance_roads_plot + local_refuges_plot +
  plot_layout(axes = "collect", ncol = 3) &
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.text.y = element_text(size = 14), # Smaller text to fit more labels
    axis.text.x = element_text(size = 16)
  )

ggsave(filename = file.path(output_dir, "pdp_all_covariates_logit_final.png"),
       plot = final_plot, width = 18, height = 15, dpi = 600, units = "in")

ggsave(filename = file.path(output_dir, "pdp_all_covariates_logit_final.pdf"),
       plot = final_plot, width = 18, height = 15, units = "in", device = cairo_pdf)





## BASE MODEL pdp facet plot - Version 3.3 (Logit Scale + Lollipop FMZ + 1k Rugs)
library(ggplot2)
library(gbm)
library(dplyr)
library(patchwork)
library(scales) 

# 1. Load data and create rug sample
raw_data <- read.csv("outputs/fire_modelling_rac_ready_forCV.csv")
set.seed(123)
rug_sample <- raw_data[sample(nrow(raw_data), 1000), ] 

# 2. Extract plot data and find global Y-axis limits
vars_to_plot <- c("ffdi_95_days", "spei12_mean", "spei24_mean", "bio5", "kbdi_95_days", 
                  "fuel_management_zones", "thunderstorm_days", "bio18", "twi", 
                  "broad_refuges", "local_refuges", "distance_roads")

all_pdp_data <- lapply(vars_to_plot, function(v) get_gbm_plot_data(model, v))

# Standardise Y limits
y_min <- min(sapply(all_pdp_data, function(d) min(d$response)))
y_max <- max(sapply(all_pdp_data, function(d) max(d$response)))
y_limits <- c(max(0.0001, y_min), min(0.9999, y_max))
logit_breaks <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)

# Helper function for continuous plots
apply_logit_style <- function(p, x_lab) {
  p + 
    scale_y_continuous(
      trans = "logit", 
      limits = y_limits,
      breaks = logit_breaks,
      labels = label_number(accuracy = 0.001)
    ) +
    labs(x = x_lab, y = "logit(Pr(Fire))") +
    theme_minimal()
}

# 3. Individual Plots
ffdi_plot <- get_gbm_plot_data(model, "ffdi_95_days") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days FFDI > 95th Percentile") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = ffdi_95_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

spei12_plot <- get_gbm_plot_data(model, "spei12_mean") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Mean SPEI (12 Months)") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = spei12_mean), inherit.aes = FALSE, alpha = 0.3, sides = "b")

spei24_plot <- get_gbm_plot_data(model, "spei24_mean") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Mean SPEI (24 Months)") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = spei24_mean), inherit.aes = FALSE, alpha = 0.3, sides = "b")

bio5_plot <- get_gbm_plot_data(model, "bio5") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Max Temp Warmest Month") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = bio5), inherit.aes = FALSE, alpha = 0.3, sides = "b")

kbdi_plot <- get_gbm_plot_data(model, "kbdi_95_days") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days KBDI > 95th Percentile") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = kbdi_95_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# FIXED FMZ Plot (Lollipop)
fmz_plot <- get_gbm_plot_data(model, "fuel_management_zones") %>%
  mutate(variable = as.factor(variable)) %>% 
  ggplot(aes(x = variable, y = response)) +
  geom_segment(aes(x=variable, xend=variable, y=y_limits[1], yend=response), color="grey50") +
  geom_point(size=4, color="black") +
  scale_y_continuous(
    trans = "logit", 
    limits = y_limits, 
    breaks = logit_breaks,
    labels = label_number(accuracy = 0.001)
  ) +
  labs(x = "Fuel Management Zone", y = "logit(Pr(Fire))") +
  theme_minimal()

thunderstorm_plot <- get_gbm_plot_data(model, "thunderstorm_days") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days of Thunderstorm") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = thunderstorm_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

bio18_plot <- get_gbm_plot_data(model, "bio18") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Precip Warmest Quarter") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = bio18), inherit.aes = FALSE, alpha = 0.3, sides = "b")

twi_plot <- get_gbm_plot_data(model, "twi") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Topographic Wetness Index") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = twi), inherit.aes = FALSE, alpha = 0.3, sides = "b")

broad_refuges_plot <- get_gbm_plot_data(model, "broad_refuges") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Landscape Refuges") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = broad_refuges), inherit.aes = FALSE, alpha = 0.3, sides = "b")

distance_roads_plot <- get_gbm_plot_data(model, "distance_roads") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Distance from Roads") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = distance_roads), inherit.aes = FALSE, alpha = 0.3, sides = "b")

local_refuges_plot <- get_gbm_plot_data(model, "local_refuges") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Local Refuges") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = local_refuges), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# 4. Final Assembly

final_plot <-
  ffdi_plot + spei12_plot + spei24_plot + bio5_plot + kbdi_plot + fmz_plot +
  thunderstorm_plot + bio18_plot + twi_plot + broad_refuges_plot +
  distance_roads_plot + local_refuges_plot +
  plot_layout(axes = "collect", ncol = 3) &
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.text.y = element_text(size = 14), 
    axis.text.x = element_text(size = 16)
  )

ggsave(filename = file.path(output_dir, "pdp_all_covariates_logit_final.png"),
       plot = final_plot, width = 18, height = 15, dpi = 600, units = "in")

ggsave(filename = file.path(output_dir, "pdp_all_covariates_logit_final.pdf"),
       plot = final_plot, width = 18, height = 15, units = "in", device = cairo_pdf)







###same for RAC model 

## pdp facet plot for report

library(ggplot2)
library(gbm)
library(patchwork)

# Simpler method using plot.gbm directly
get_gbm_plot_data <- function(gbm_model, var_name) {
  
  # Use plot.gbm to generate the plot data
  plot_data <- plot.gbm(gbm_model, i.var = var_name, 
                        return.grid = TRUE, type = "response")
  
  # Extract the data
  data.frame(
    variable = plot_data[[var_name]],  # Use variable name as column name
    response = plot_data$y
  )
}

# Create Response Plots
ffdi_plot <- get_gbm_plot_data(rac.model, "ffdi_95_days") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Days FFDI>95th Percentile", 
       y = "Pr(Fire)") +
  theme_minimal()

spei12_plot <- get_gbm_plot_data(rac.model, "spei12_mean") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Mean SPEI (12 Months)", 
       y = "Pr(Fire)") +
  theme_minimal()

spei24_plot <- get_gbm_plot_data(rac.model, "spei24_mean") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Mean SPEI (24 Months)", 
       y = "Pr(Fire)") +
  theme_minimal()

bio5_plot <- get_gbm_plot_data(rac.model, "bio5") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Max Temperature of Warmest Month", 
       y = "Pr(Fire)") +
  theme_minimal()

kbdi_plot <- get_gbm_plot_data(rac.model, "kbdi_95_days") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Days KBDI >95th Percentile", 
       y = "Pr(Fire)") +
  theme_minimal()


fmz_plot <- get_gbm_plot_data(rac.model, "fuel_management_zones") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_bar(color = "black", size = 1, stat="identity") +
  labs(x = "Fuel Management Zone", 
       y = "Pr(Fire)") +
  theme_minimal()

thunderstorm_plot <- get_gbm_plot_data(rac.model, "thunderstorm_days") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Days of Thunderstorm Occurrence", 
       y = "Pr(Fire)") +
  theme_minimal()

bio18_plot <- get_gbm_plot_data(rac.model, "bio18") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Precipitation of Warmest Quarter", 
       y = "Pr(Fire)") +
  theme_minimal()

twi_plot <- get_gbm_plot_data(rac.model, "twi") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Topographic Wetness Index", 
       y = "Pr(Fire)") +
  theme_minimal()

broad_refuges_plot <- get_gbm_plot_data(rac.model, "broad_refuges") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Landscape Refuges", 
       y = "Pr(Fire)") +
  theme_minimal()

local_refuges_plot <- get_gbm_plot_data(rac.model, "local_refuges") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Local Refuges", 
       y = "Pr(Fire)") +
  theme_minimal()

distance_roads_plot <- get_gbm_plot_data(rac.model, "distance_roads") %>%
  ggplot(aes(x = variable, y = response)) +
  geom_line(color = "black", size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  labs(x = "Distance from Roads", 
       y = "Pr(Fire)") +
  theme_minimal()



#theme_set(theme_minimal(base_size = 16))

final_plot <-
  ffdi_plot + spei12_plot + spei24_plot + bio5_plot + kbdi_plot + fmz_plot +
  thunderstorm_plot + bio18_plot + twi_plot + broad_refuges_plot +
  distance_roads_plot + local_refuges_plot +
  plot_layout(axes = "collect", ncol = 3) &
  theme(
    axis.title.x = element_text(
      size = 22,
      margin = margin(t = 10)   
    ),
    axis.title.y = element_text(
      size = 24,
      margin = margin(r = 10)   
    ),
    axis.text = element_text(size = 23)
  )
ggsave(
  filename = file.path(output_dir, "RAC_pdp_all_covariates_highres.png"),
  plot = final_plot,
  width = 18,
  height = 15,
  dpi = 600,
  units = "in"
)
ggsave(
  filename = file.path(output_dir, "RAC_pdp_all_covariates_highres.pdf"),
  plot = final_plot,
  width = 18,
  height = 15,
  units = "in",
  device = cairo_pdf
)




## RAC MODEL pdp facet plot - Version 3.6 (Cleaned Axis + Lollipop FMZ)
library(ggplot2)
library(gbm)
library(dplyr)
library(patchwork)
library(scales) 

# 1. Load data and create rug sample
raw_data <- read.csv("outputs/fire_modelling_rac_ready_forCV.csv")
set.seed(123)
rug_sample <- raw_data[sample(nrow(raw_data), 1000), ] 

# 2. Extract plot data and define coarse limits/breaks
vars_to_plot <- c("ffdi_95_days", "spei12_mean", "spei24_mean", "bio5", "kbdi_95_days", 
                  "fuel_management_zones", "thunderstorm_days", "bio18", "twi", 
                  "broad_refuges", "local_refuges", "distance_roads")

all_pdp_data <- lapply(vars_to_plot, function(v) get_gbm_plot_data(rac.model, v))

# Limits with slight padding at the top for clarity
y_limits_rac <- c(0.001, 0.022) 
# Coarser breaks to prevent label overlapping
logit_breaks_rac <- c(0.001, 0.002, 0.005, 0.01, 0.015, 0.02)

# Helper function for continuous plots
apply_logit_style <- function(p, x_lab) {
  p + 
    scale_y_continuous(
      trans = "logit", 
      limits = y_limits_rac,
      breaks = logit_breaks_rac,
      labels = label_number(accuracy = 0.001)
    ) +
    labs(x = x_lab, y = "logit(Pr(Fire))") +
    theme_minimal()
}

# 3. Individual Plots
ffdi_plot <- get_gbm_plot_data(rac.model, "ffdi_95_days") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days FFDI > 95th Percentile") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = ffdi_95_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

spei12_plot <- get_gbm_plot_data(rac.model, "spei12_mean") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Mean SPEI (12 Months)") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = spei12_mean), inherit.aes = FALSE, alpha = 0.3, sides = "b")

spei24_plot <- get_gbm_plot_data(rac.model, "spei24_mean") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Mean SPEI (24 Months)") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = spei24_mean), inherit.aes = FALSE, alpha = 0.3, sides = "b")

bio5_plot <- get_gbm_plot_data(rac.model, "bio5") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Max Temp Warmest Month") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = bio5), inherit.aes = FALSE, alpha = 0.3, sides = "b")

kbdi_plot <- get_gbm_plot_data(rac.model, "kbdi_95_days") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days KBDI > 95th Percentile") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = kbdi_95_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

# RAC FMZ Plot (Lollipop)
fmz_plot <- get_gbm_plot_data(rac.model, "fuel_management_zones") %>%
  mutate(variable = as.factor(variable)) %>% 
  ggplot(aes(x = variable, y = response)) +
  geom_segment(aes(x=variable, xend=variable, y=y_limits_rac[1], yend=response), color="grey50") +
  geom_point(size=4, color="black") +
  scale_y_continuous(
    trans = "logit", 
    limits = y_limits_rac, 
    breaks = logit_breaks_rac, 
    labels = label_number(accuracy = 0.001)
  ) +
  labs(x = "Fuel Management Zone", y = "logit(Pr(Fire))") +
  theme_minimal()

thunderstorm_plot <- get_gbm_plot_data(rac.model, "thunderstorm_days") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Days of Thunderstorm") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = thunderstorm_days), inherit.aes = FALSE, alpha = 0.3, sides = "b")

bio18_plot <- get_gbm_plot_data(rac.model, "bio18") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Precip Warmest Quarter") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = bio18), inherit.aes = FALSE, alpha = 0.3, sides = "b")

twi_plot <- get_gbm_plot_data(rac.model, "twi") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Topographic Wetness Index") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = twi), inherit.aes = FALSE, alpha = 0.3, sides = "b")

broad_refuges_plot <- get_gbm_plot_data(rac.model, "broad_refuges") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Landscape Refuges") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = broad_refuges), inherit.aes = FALSE, alpha = 0.3, sides = "b")

distance_roads_plot <- get_gbm_plot_data(rac.model, "distance_roads") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Distance from Roads") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = distance_roads), inherit.aes = FALSE, alpha = 0.3, sides = "b")

local_refuges_plot <- get_gbm_plot_data(rac.model, "local_refuges") %>% ggplot(aes(x = variable, y = response)) %>%
  apply_logit_style("Local Refuges") + geom_line(size = 1) +
  geom_smooth(method="loess", span = 1, se = FALSE, color = "blue") +
  geom_rug(data = rug_sample, aes(x = local_refuges), inherit.aes = FALSE, alpha = 0.3, sides = "b")



final_plot_rac <-
  ffdi_plot + spei12_plot + spei24_plot + bio5_plot + kbdi_plot + fmz_plot +
  thunderstorm_plot + bio18_plot + twi_plot + broad_refuges_plot +
  distance_roads_plot + local_refuges_plot +
  plot_layout(axes = "collect", ncol = 3) &
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.text.y = element_text(size = 14), 
    axis.text.x = element_text(size = 16)
  )

ggsave(filename = file.path(output_dir, "pdp_rac_model_logit_final.png"),
       plot = final_plot_rac, width = 18, height = 15, dpi = 600, units = "in")

ggsave(filename = file.path(output_dir, "pdp_rac_model_logit_final.pdf"),
       plot = final_plot_rac, width = 18, height = 15, units = "in", device = cairo_pdf)





# Variable Importance
library(tidyverse)
labels = data.frame(var = c("ffdi_95_days","spei12_mean","spei24_mean","bio5","kbdi_95_days","fuel_management_zones", "thunderstorm_days",
                            "bio18","twi","broad_refuges","distance_roads","local_refuges"),
                    var.lab = c("Days FFDI>95th Percentile", "Mean SPEI (12 Months)", "Mean SPEI (24 Months)", "Max Temperature of Warmest Month", 
                                "Days KBDI >95th Percentile", "Fuel Management Zone", "Days of Thunderstorm Occurrence", 
                                "Precipitation of Warmest Quarter", "Topographic Wetness Index", 
                                "Landscape Refuges", "Distance from Road", "Local Refuges"))
vimp = data.frame(summary(model))
vimp = left_join(vimp, labels)
vimp = vimp %>% mutate(var.lab = fct_reorder(var.lab, rel.inf))
ggplot(vimp) + 
  geom_bar(aes(x=var.lab, y =rel.inf), stat="identity") +
  theme_bw() + labs(x="", y="Relative Importance (%)") + coord_flip() 



print(vimp)

# Variable Importance
library(tidyverse)

# Labels
labels <- data.frame(
  var = c(
    "ffdi_95_days", "spei12_mean", "spei24_mean", "bio5", "kbdi_95_days",
    "fuel_management_zones", "thunderstorm_days", "bio18", "twi",
    "broad_refuges", "distance_roads", "local_refuges"
  ),
  var.lab = c(
    "Days FFDI>95th Percentile", "Mean SPEI (12 Months)", "Mean SPEI (24 Months)",
    "Max Temperature of Warmest Month", "Days KBDI >95th Percentile",
    "Fuel Management Zone", "Days of Thunderstorm Occurrence",
    "Precipitation of Warmest Quarter", "Topographic Wetness Index",
    "Landscape Refuges", "Distance from Road", "Local Refuges"
  )
)

# Dynamic vs static classification
labels$predictor_type <- ifelse(
  labels$var %in% c(
    "ffdi_95_days", "spei12_mean", "spei24_mean", 
    "kbdi_95_days", "thunderstorm_days"
  ),
  "Dynamic",
  "Static"
)

# Import model importance
vimp <- data.frame(summary(model)) |>
  left_join(labels, by = "var") |>
  mutate(var.lab = fct_reorder(var.lab, rel.inf))

# Plot
p_vimp <- ggplot(vimp, aes(x = var.lab, y = rel.inf, fill = predictor_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Dynamic" = "#1B9E77", "Static" = "#756BB1")) +
  theme_bw(base_size = 14) +
  labs(x = "",
       y = "Relative Importance (%)",
       fill = "Predictor type") +
  coord_flip()

p_vimp


ggsave(
  filename = file.path(output_dir, "figure3_predictor_importance.pdf"),
  plot = p_vimp,
  width = 8,
  height = 6,
  units = "in"
)

# Optional PNG (high res)
ggsave(
  filename = file.path(output_dir, "figure3_predictor_importance.png"),
  plot = p_vimp,
  width = 8,
  height = 6,
  dpi = 600,
  units = "in"
)







##Faceted VIP plots 

library(tidyverse)
library(forcats)
library(patchwork)

# ---------------------------------------------------
# Labels + predictor types (add RAC)
# ---------------------------------------------------
labels <- data.frame(
  var = c(
    "ffdi_95_days","spei12_mean","spei24_mean","bio5","kbdi_95_days",
    "fuel_management_zones","thunderstorm_days","bio18","twi",
    "broad_refuges","distance_roads","local_refuges","rac"
  ),
  var.lab = c(
    "Days FFDI>95th Percentile", "Mean SPEI (12 Months)", "Mean SPEI (24 Months)",
    "Max Temp (Warmest Month)", "Days KBDI >95th Percentile",
    "Fuel Management Zone", "Thunderstorm Days",
    "Precip. Warmest Quarter", "Topographic Wetness Index",
    "Landscape Refuges", "Distance from Roads", "Local Refuges",
    "Residual Autocovariate"
  )
)

labels$predictor_type <- ifelse(
  labels$var %in% c(
    "ffdi_95_days","spei12_mean","spei24_mean",
    "kbdi_95_days","thunderstorm_days"
  ),
  "Dynamic", "Static"
)
labels$predictor_type[labels$var == "rac"] <- "Static"

# ---------------------------------------------------
# Extract variable importance
# ---------------------------------------------------
vimp_base <- data.frame(summary(model)) %>% mutate(model = "BASE")
vimp_rac  <- data.frame(summary(rac.model)) %>% mutate(model = "RAC")

# ---------------------------------------------------
# Join labels & remove NA rows
# ---------------------------------------------------
v_base <- vimp_base %>%
  left_join(labels, by = "var") %>%
  filter(!is.na(var.lab))

v_rac <- vimp_rac %>%
  left_join(labels, by = "var") %>%
  filter(!is.na(var.lab))

# ---------------------------------------------------
# Reorder each model independently (separate scales)
# ---------------------------------------------------
v_base <- v_base %>%
  mutate(var.lab = fct_reorder(var.lab, rel.inf))

v_rac <- v_rac %>%
  mutate(var.lab = fct_reorder(var.lab, rel.inf))

# ---------------------------------------------------
# Build RAC plot (panel a)
# ---------------------------------------------------
p_rac <- ggplot(v_rac, aes(x = var.lab, y = rel.inf, fill = predictor_type)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "",
    y = "Relative Importance (%)",
    fill = "Predictor type",
    title = "a)"
  ) +
  scale_fill_manual(values = c("Dynamic" = "#1B9E77", "Static" = "#756BB1")) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right"
  )

# ---------------------------------------------------
# Build BASE plot (panel b)
# ---------------------------------------------------
p_base <- ggplot(v_base, aes(x = var.lab, y = rel.inf, fill = predictor_type)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "",
    y = "Relative Importance (%)",
    fill = "Predictor type",
    title = "b)"
  ) +
  scale_fill_manual(values = c("Dynamic" = "#1B9E77", "Static" = "#756BB1")) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "none"  # shared legend from RAC plot
  )

# ---------------------------------------------------
#  Combine side by side with shared legend
# ---------------------------------------------------
p_final <- p_rac + p_base +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

p_final

# ---------------------------------------------------
# Save output
# ---------------------------------------------------

ggsave(
  file.path(output_dir, "VIP_RAC_vs_BASE_side_by_side.pdf"),
  p_final,
  width = 14, height = 7, units = "in"
)

ggsave(
  file.path(output_dir, "VIP_RAC_vs_BASE_side_by_side.png"),
  p_final,
  width = 14, height = 7, dpi = 600, units = "in"
)




tsf.model = readRDS("./outputs/brt_base_lr0.005_tsf_ffdi_95.rds")

plot.gbm(tsf.model, c("thunderstorm_days", ""), type = "response") # THis looks good
plot.gbm(tsf.model, c("fuel_management_zones", "tsf"), type = "response") # THis looks good
#tsf












##Forecasts



##plot of baseline fire risk

library(terra)
library(tidyterra)
library(ggplot2)
library(viridis)
library(sf)
library(ozmaps)

# --- Load raster ---------------------------------------------
pred_dir <- "outputs/predictions"
f_base <- file.path(pred_dir, "fire_observed_baseline_19851995.tif")
r_base <- rast(f_base); names(r_base) <- "prob"

# --- Clipping extents ----------------------------------------
r_ext <- ext(r_base)
ymin <- terra::ymin(r_ext)
ymax <- terra::ymax(r_ext)

x_buffer <- 150000  # 150 km buffer
xlim <- c(terra::xmin(r_ext) - x_buffer, terra::xmax(r_ext) + x_buffer)

# --- Binning setup -------------------------------------------
breaks_fixed <- c(0, 0.01, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50)
labels_fixed <- c("0–0.01", "0.01–0.05", "0.05–0.10", "0.10–0.15",
                  "0.15–0.25", "0.25–0.35", "0.35–0.50")

r_binned <- classify(r_base, rcl = matrix(c(breaks_fixed[-length(breaks_fixed)],
                                            breaks_fixed[-1],
                                            1:length(labels_fixed)), ncol = 3))

# Apply labels to factor
levels(r_binned) <- data.frame(
  ID    = 1:length(labels_fixed),
  label = labels_fixed
)

# --- Colour palette ------------------------------------------
palette_inf <- viridis(n = length(labels_fixed), option = "inferno", direction = 1)

# --- Plot ----------------------------------------------------
ggplot() +
  geom_spatraster(data = r_binned) +
  scale_fill_manual(
    values = palette_inf,
    name = "Probability of fire",
    na.value = "transparent"
  ) +
  guides(fill = guide_legend(na.translate = FALSE)) +
  coord_sf(xlim = xlim, ylim = c(ymin, ymax), expand = FALSE) +
  labs(title = "Baseline Fire Probability (Binned to 0.50)") +
  theme_void(base_size = 11) +  # ✅ removes grid and axes
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# --- Save plot -----------------------------------------------
out_path <- file.path(output_dir, "baseline_probability_map.png")

ggsave(
  filename = out_path,
  width = 8, height = 10, units = "in", dpi = 300,
  device = ragg::agg_png
)
message("✅ Saved: ", out_path)







## plot of changing future fire risk
# ================================
# Binned ΔPr(Fire) Plots — With & Without Fire Regions
# ================================

# --- Libraries --------------------------------------------------
libs <- c("terra", "ggplot2", "tidyterra", "sf", "dplyr", "stringr")
invisible(lapply(setdiff(libs, rownames(installed.packages())), install.packages))
invisible(lapply(libs, library, character.only = TRUE))

# --- Paths ------------------------------------------------------
pred_dir <- "./outputs/predictions"
fd_path  <- "./data/spatial/LF_DISTRICT.shp"
fig_dir  <- "./plots"
out_base <- file.path(fig_dir, "future_minus_baseline_binned.png")
out_regions <- file.path(fig_dir, "future_minus_baseline_binned_regions.png")

# --- Raster load & difference -----------------------------------
f_base <- file.path(pred_dir, "fire_observed_baseline_19851995.tif")
f_fut  <- file.path(pred_dir, "fire_ACCESS1-0_rcp85_20812099.tif")
r_base <- rast(f_base); names(r_base) <- "prob"
r_fut  <- rast(f_fut);  names(r_fut)  <- "prob"

if (!compareGeom(r_fut, r_base, stopOnError = FALSE)) {
  r_fut <- resample(r_fut, r_base, method = "bilinear")
}
r_diff <- clamp(r_fut - r_base, -0.5, 0.5, values = TRUE)

# --- Downsample for plotting ------------------------------------
target_maxcells <- 3e5
fact <- max(1L, ceiling(sqrt(ncell(r_diff) / target_maxcells)))
r_plot <- if (fact > 1L) aggregate(r_diff, fact = fact, fun = "mean", na.rm = TRUE) else r_diff

# --- Bin classification setup -----------------------------------
step <- 0.05; zero_tol <- 0.01; neg_min <- -0.15; pos_max <- 0.20
neg_from <- seq(-step, neg_min, by = -step); neg_to <- neg_from + step
pos_from <- seq(0.00, pos_max - step, by = step); pos_to <- pos_from + step
n_neg <- length(neg_from); n_pos <- length(pos_from); zero_id <- n_neg + 1L

bin_fun <- function(v) {
  out <- rep(NA_integer_, length(v))
  nz <- is.finite(v); vv <- v[nz]
  out[nz][vv >= -zero_tol & vv <= zero_tol] <- zero_id
  idx_neg <- which(vv < -zero_tol)
  if (length(idx_neg)) {
    k <- ceiling((-vv[idx_neg]) / step); k <- pmin(pmax(k, 1L), n_neg)
    out[nz][idx_neg] <- k
  }
  idx_pos <- which(vv > zero_tol)
  if (length(idx_pos)) {
    k <- ceiling(vv[idx_pos] / step); k <- pmin(pmax(k, 1L), n_pos)
    out[nz][idx_pos] <- zero_id + k
  }
  out
}

r_class <- app(r_plot, bin_fun)

# --- Labels & palette -------------------------------------------
neg_labs <- sprintf("%.2f to %.2f", neg_from, neg_to)
neg_labs[neg_to == 0] <- "-0.05 to 0"
pos_labs <- sprintf("%.2f to %.2f", pos_from, pos_to)
pos_labs[pos_from == 0] <- "0 to 0.05"
zero_lab <- "−0.01 to 0.01"
r_class <- as.factor(r_class)
levels(r_class) <- data.frame(
  ID = c(seq_len(n_neg), zero_id, zero_id + seq_len(n_pos)),
  bin = c(neg_labs, zero_lab, pos_labs)
)

neg_cols <- colorRampPalette(c("#08519C", "#6BAED6", "#DEEBF7"))(n_neg)
pos_cols <- colorRampPalette(c("#FEE5D9", "#FB6A4A", "#CB181D", "#67000D"))(n_pos)
zero_col <- "#BFBFBF"

val_cols <- c(
  setNames(neg_cols, rev(neg_labs)),
  setNames(zero_col, zero_lab),
  setNames(pos_cols, pos_labs)
)

breaks_ord <- c(rev(neg_labs), zero_lab, pos_labs)
ex <- ext(r_class)

# --- Fire region shapefile --------------------------------------
fdist <- st_read(fd_path, quiet = TRUE) |> st_make_valid()
regions <- fdist |>
  group_by(REGIONNAME) |>
  summarise(.groups = "drop") |>
  st_transform(st_crs(r_class))
regions$REGIONNAME <- str_to_title(regions$REGIONNAME)
region_labels <- st_point_on_surface(regions["REGIONNAME"])

# --- Common theme + layers --------------------------------------
base_map <- function(add_layers = NULL) {
  ggplot() +
    geom_spatraster(data = r_class) +
    scale_fill_manual(
      values = val_cols, breaks = breaks_ord,
      name = "Change in probability of fire",
      na.value = "transparent"
    ) +
    guides(fill = guide_legend(na.translate = FALSE)) +
    (add_layers %||% NULL) +
    coord_sf(
      xlim = c(xmin(ex) - 150000, xmax(ex) + 150000),
      ylim = c(ymin(ex), ymax(ex)),
      expand = FALSE
    ) +
    labs(title = "Future – Baseline: Change in Fire Probability (Binned)") +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

# --- Plot 1: No shapefile ---------------------------------------
p_clean <- base_map()

# --- Plot 2: Fire regions overlaid ------------------------------
p_regions <- base_map(
  list(
    geom_sf(data = regions, fill = NA, color = "black", linewidth = 0.3),
    geom_sf_text(data = region_labels, aes(label = REGIONNAME),
                 size = 2.5, fontface = "plain", color = "black", check_overlap = TRUE)
  )
)

# --- Save both --------------------------------------------------
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

ragg::agg_png(out_base, width = 9, height = 10, units = "in", res = 300)
print(p_clean)
dev.off()
message("✅ Saved: ", out_base)

ragg::agg_png(out_regions, width = 9, height = 10, units = "in", res = 300)
print(p_regions)
dev.off()
message("✅ Saved: ", out_regions)

















##Predictions BASE MODEL

##read all predictions 
library(terra)
library(ggplot2)
library(tidyterra)
library(viridis)

pred_dir <- "outputs/predictions"

# === Observed Baseline prediction (1985–2005) ===
r_base_OBS <- rast(
  file.path(pred_dir, "fire_observed_baseline_1985_2005.tif")
)

# === Future predictions (2081–2099) ===
r_fut_ACCESS1  <- rast(file.path(pred_dir, "fire_ACCESS1_0_RCP85_future_2081_2099.tif"))
r_fut_CNRM     <- rast(file.path(pred_dir, "fire_CNRM_CM5_RCP85_future_2081_2099.tif"))
r_fut_GFDL     <- rast(file.path(pred_dir, "fire_GFDL_ESM2M_RCP85_future_2081_2099.tif"))
r_fut_MIROC    <- rast(file.path(pred_dir, "fire_MIROC5_RCP85_future_2081_2099.tif"))
# Ensemble future
r_fut_ENSEMBLE <- rast(file.path(pred_dir, "fire_Ensemble_Mean_Future_2081_2099.tif"))

# === Difference maps (Future – Baseline) ===
r_diff_ACCESS1  <- rast(file.path(pred_dir, "fire_change_ACCESS1_0_RCP85_2081_2099.tif"))
r_diff_CNRM     <- rast(file.path(pred_dir, "fire_change_CNRM_CM5_RCP85_2081_2099.tif"))
r_diff_GFDL     <- rast(file.path(pred_dir, "fire_change_GFDL_ESM2M_RCP85_2081_2099.tif"))
r_diff_MIROC    <- rast(file.path(pred_dir, "fire_change_MIROC5_RCP85_2081_2099.tif"))
# Ensemble difference map
r_diff_ENSEMBLE <- rast(file.path(pred_dir, "fire_change_Ensemble_Mean_2081_2099.tif"))

summary(r_base_OBS)
plot(r_base_OBS, main = "Observed Baseline (Masked to Native Vegetation)")


# ============================================================
# ADD CLIPPING STEPS HERE - Native Vegetation Mask
# ============================================================

library(terra)

# --- Load Native Vegetation Raster ---
veg_path <- "./data/raw/NVR2017_CONDITION.tif"
r_veg <- rast(veg_path)

cat("\nLoaded veg raster:\n")
print(r_veg)
plot(r_veg, main = "Raw Native Vegetation Raster")

# ------------------------------------------------------------
# Reproject vegetaton raster → match prediction raster CRS
# (GDA2020 / VicGrid: use nearest neighbour for categorical data)
# ------------------------------------------------------------
r_veg_proj <- project(r_veg, r_base_OBS, method = "near")

# Quick inspect values
vals_sample <- unique(values(r_veg_proj)[1:10000])
cat("\nSample veg raster values (after reprojection):\n")
print(vals_sample)

# ------------------------------------------------------------
# Create native vegetation mask
# If NVR2017 is coded: 0 = cleared, 1–100 = native vegetation condition,
# then simply use >0 to keep native vegetation.
# ------------------------------------------------------------
r_veg_mask <- r_veg_proj > 0

plot(r_veg_mask, main = "Native Vegetation Mask (1 = native vegetation)")

# ============================================================
# ALIGN FUTURE RASTERS TO OBSERVED BASELINE GEOMETRY
# ============================================================

align_to_baseline <- function(r) {
  resample(r, r_base_OBS, method = "bilinear")
}

r_fut_ACCESS1  <- align_to_baseline(r_fut_ACCESS1)
r_fut_CNRM     <- align_to_baseline(r_fut_CNRM)
r_fut_GFDL     <- align_to_baseline(r_fut_GFDL)
r_fut_MIROC    <- align_to_baseline(r_fut_MIROC)
r_fut_ENSEMBLE <- align_to_baseline(r_fut_ENSEMBLE)


# ============================================================
# MASK ALL MODEL RASTERS (Baseline, Future, Delta)
# ============================================================

mask_pred <- function(r) {
  mask(r, r_veg_mask, maskvalues = 0)   # mask out non-native vegetation
}

# --- Observed Baseline ---
r_base_OBS <- mask_pred(r_base_OBS)

# --- Future ---
r_fut_ACCESS1   <- mask_pred(r_fut_ACCESS1)
r_fut_CNRM      <- mask_pred(r_fut_CNRM)
r_fut_GFDL      <- mask_pred(r_fut_GFDL)
r_fut_MIROC     <- mask_pred(r_fut_MIROC)
r_fut_ENSEMBLE  <- mask_pred(r_fut_ENSEMBLE)

# # --- Delta (future – baseline) ---
# r_diff_ACCESS1  <- mask_pred(r_diff_ACCESS1)
# r_diff_CNRM     <- mask_pred(r_diff_CNRM)
# r_diff_GFDL     <- mask_pred(r_diff_GFDL)
# r_diff_MIROC    <- mask_pred(r_diff_MIROC)
# r_diff_ENSEMBLE <- mask_pred(r_diff_ENSEMBLE)

# Visual confirmation
plot(r_base_ENSEMBLE, main = "Baseline Ensemble (Masked to Native Vegetation)")


##added code for resuts (Min/Max/Mean)
## ============================================================
## TABLE X: GLOBAL MIN / MAX / MEAN (BASE / non-RAC MODEL)
## ============================================================

library(dplyr)

summarise_pred <- function(r) {
  v <- values(r, na.rm = TRUE)
  tibble(
    min  = min(v),
    max  = max(v),
    mean = mean(v)
  )
}

table_x_base <- bind_rows(
  
  # --- Observed baseline (single) ---
  summarise_pred(r_base_OBS) %>%
    mutate(period = "Baseline (Observed)", model = "Observed"),
  
  # --- Future GCMs ---
  summarise_pred(r_fut_ACCESS1) %>%
    mutate(period = "Future", model = "ACCESS1"),
  
  summarise_pred(r_fut_GFDL) %>%
    mutate(period = "Future", model = "GFDL"),
  
  summarise_pred(r_fut_MIROC) %>%
    mutate(period = "Future", model = "MIROC"),
  
  summarise_pred(r_fut_CNRM) %>%
    mutate(period = "Future", model = "CNRM"),
  
  # --- Future Ensemble ---
  summarise_pred(r_fut_ENSEMBLE) %>%
    mutate(period = "Future", model = "Ensemble")
) %>%
  select(period, model, min, max, mean) %>%
  arrange(period, model)

print(table_x_base)




#plot ensembles
# --- Binning scheme (same as before) ---
PROB_BREAKS <- c(0, 0.01, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50)
PROB_LABELS <- c(
  "0–0.01", "0.01–0.05", "0.05–0.10", "0.10–0.15",
  "0.15–0.25", "0.25–0.35", "0.35–0.50"
)
PROB_PALETTE <- viridis(
  n = length(PROB_LABELS),
  option = "inferno",
  direction = 1
)

# --- Function to produce the binned map ---
plot_prob_map <- function(r_input, title_text) {
  
  # 1. Classify
  r_binned <- classify(
    r_input,
    rcl = matrix(c(
      PROB_BREAKS[-length(PROB_BREAKS)],
      PROB_BREAKS[-1],
      1:length(PROB_LABELS)
    ), ncol = 3)
  )
  
  # 2. Apply labels
  levels(r_binned) <- data.frame(
    ID = 1:length(PROB_LABELS),
    label = PROB_LABELS
  )
  
  # 3. Plot (no saving yet)
  ggplot() +
    geom_spatraster(data = r_binned) +
    scale_fill_manual(
      values = PROB_PALETTE,
      name = "Probability of fire",
      na.value = "transparent",
      na.translate = FALSE
    ) +
    guides(fill = guide_legend(na.translate = FALSE)) +
    labs(title = title_text) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 14),     
      legend.title = element_text(size = 15),    
      legend.key.size = unit(0.8, "cm")           
    )
}

# --- Make plots ---

p_base_ens <- plot_prob_map(r_base_ENSEMBLE, "")#title removed Ensemble Baseline Prediction
p_base_ens

p_fut_ens <- plot_prob_map(r_fut_ENSEMBLE, "")#removed Ensemble Future Prediction
p_fut_ens





##plot delta map 
names(r_diff_ENSEMBLE) <- "delta_pr"

# --- Define simple positive-only breaks ---
delta_breaks <- c(0, 0.01, 0.05, 0.10, 0.20, 0.40, 0.60, 0.70)

delta_labels <- c(
  "0–0.01",
  "0.01–0.05",
  "0.05–0.10",
  "0.10–0.20",
  "0.20–0.40",
  "0.40–0.60",
  ">0.60"
)

# Classification matrix
rcl <- cbind(
  delta_breaks[-length(delta_breaks)],
  delta_breaks[-1],
  1:(length(delta_breaks)-1)
)

r_delta_bin <- classify(r_diff_ENSEMBLE, rcl)

# Apply labels
levels(r_delta_bin) <- data.frame(
  ID = 1:length(delta_labels),
  label = delta_labels
)

# --- Colour palette (light → dark red) ---
pal <- c(
  "#CCCCCC",  # near-zero
  "#FEE5D9",
  "#FCBBA1",
  "#FC9272",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15"   # big increase
)

# --- Plot ---
p_delta <- ggplot() +
  geom_spatraster(data = r_delta_bin, maxcell = 2e6) +
  scale_fill_manual(
    values = pal,
    labels = delta_labels,
    name = expression(paste(Delta, "Pr(Fire)")),
    na.value = "transparent",
    na.translate = FALSE
  ) +
  labs(title = "") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 14),      
    legend.title = element_text(size = 15),     
    legend.key.size = unit(0.8, "cm")          
  )

p_delta





# --- Load fire district shapefile ---
fd_path  <- "data/spatial/LF_DISTRICT.shp"

regions <- st_read(fd_path, quiet = TRUE) |>
  st_make_valid() |>
  group_by(REGIONNAME) |>
  summarise(.groups = "drop")

# Transform to raster CRS
regions <- st_transform(regions, crs(r_delta_bin))

# Label points (so labels appear inside polygons)
region_labels <- st_point_on_surface(regions["REGIONNAME"])


p_delta_regions <- ggplot() +
  geom_spatraster(data = r_delta_bin, maxcell = 2e6) +
  geom_sf(data = regions, fill = NA, color = "black", linewidth = 0.3) +
  geom_sf_text(
    data = region_labels,
    aes(label = REGIONNAME),
    size = 3,
    color = "black",
    fontface = "plain",
    check_overlap = TRUE
  ) +
  scale_fill_manual(
    values = pal,
    labels = delta_labels,
    name = expression(paste(Delta, "Probability of fire")),
    na.value = "transparent",
    na.translate = FALSE
  ) +
  labs(title = "") +
  theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 14),      # bigger legend labels
      legend.title = element_text(size = 15),     # bigger legend title
      legend.key.size = unit(0.8, "cm")           # bigger colour boxes
    )

p_delta_regions



##plot all 3
p_composite <- 
  p_base_ens /
  p_fut_ens /
  p_delta_regions +
  plot_annotation(
    title = "",
    tag_levels = "a",       # adds a), b), c)
    tag_prefix = "",        # nothing before "a"
    tag_suffix = ") ",      # formats as a) b) c)
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(0, 0, 10, 0)
      ),
      plot.tag = element_text(size = 20, face = "bold", hjust = -0.2)  # makes labels larger & shifts left
    )
  )


p_composite

ggsave(
  filename = file.path(output_dir, "fire_risk_composite.jpg"),
  plot = p_composite,
  width = 10,      # adjust if needed
  height = 18,     # taller because 3 rows
  dpi = 400,
  units = "in",
  device = "jpg"
)








## ============================================================================
## Regional summaries (
## ============================================================================

library(terra)
library(dplyr)
library(readr)
library(janitor)
library(purrr)

## ------------------------------------------------------------
## 1. READ + FIX REGIONS
## ------------------------------------------------------------


regions_sf <- st_read(fd_path, quiet = TRUE) |> 
  st_make_valid()

regions_v <- vect(regions_sf)

# Project regions into raster CRS
regions_v <- terra::project(regions_v, crs(r_base_ENSEMBLE))

# Dissolve into 6 regions
regions_v <- aggregate(regions_v, by = "REGIONNAME")

print(regions_v)

## ------------------------------------------------------------
## 2. MODEL–RASTER LIST
## ------------------------------------------------------------

model_rasters <- list(
  ACCESS1 = list(base = r_base_ACCESS1,  fut = r_fut_ACCESS1,  delta = r_diff_ACCESS1),
  CNRM    = list(base = r_base_CNRM,     fut = r_fut_CNRM,     delta = r_diff_CNRM),
  GFDL    = list(base = r_base_GFDL,     fut = r_fut_GFDL,     delta = r_diff_GFDL),
  MIROC   = list(base = r_base_MIROC,    fut = r_fut_MIROC,    delta = r_diff_MIROC),
  Ensemble= list(base = r_base_ENSEMBLE, fut = r_fut_ENSEMBLE, delta = r_diff_ENSEMBLE)
)

ha_per_cell <- (75 * 75) / 10000    # 0.5625 ha

## ------------------------------------------------------------
## 3. NEW FUNCTION: clip each raster to the region FIRST
## ------------------------------------------------------------

summarise_raster <- function(r, region_poly, prefix) {
  
  # crop + mask raster to region
  r_clip <- crop(r, region_poly)
  r_clip <- mask(r_clip, region_poly)
  
  vals <- values(r_clip)
  vals <- vals[!is.na(vals)]
  
  n_cells <- length(vals)
  
  tibble(
    region = region_poly$REGIONNAME,
    !!paste0(prefix, "_n_cells") := n_cells,
    !!paste0(prefix, "_total_area_ha") := n_cells * ha_per_cell,
    !!paste0(prefix, "_mean") := ifelse(n_cells > 0, mean(vals), NA),
    !!paste0(prefix, "_median") := ifelse(n_cells > 0, median(vals), NA),
    !!paste0(prefix, "_prop_gt_0_10") := ifelse(n_cells > 0, mean(vals > 0.10), NA),
    !!paste0(prefix, "_area_gt_0_10_ha") := ifelse(n_cells > 0, sum(vals > 0.10) * ha_per_cell, NA)
  )
}

## ------------------------------------------------------------
## 4. LOOP THROUGH MODELS + REGIONS
## ------------------------------------------------------------

summary_list <- list()

for (m in names(model_rasters)) {
  cat("Processing model:", m, "\n")
  
  base_res  <- map_dfr(1:nrow(regions_v), ~ summarise_raster(model_rasters[[m]]$base,  regions_v[.x], "baseline"))
  fut_res   <- map_dfr(1:nrow(regions_v), ~ summarise_raster(model_rasters[[m]]$fut,   regions_v[.x], "future"))
  delta_res <- map_dfr(1:nrow(regions_v), ~ summarise_raster(model_rasters[[m]]$delta, regions_v[.x], "delta"))
  
  df <- base_res |>
    left_join(fut_res,   by = "region") |>
    left_join(delta_res, by = "region") |>
    mutate(
      model = m,
      relative_change_pct =
        100 * (future_mean - baseline_mean) / pmax(baseline_mean, 0.001)
    )
  
  summary_list[[m]] <- df
}

regional_summary <- bind_rows(summary_list) |>
  arrange(region, model) |>
  clean_names()

## ------------------------------------------------------------
## 5. SAVE OUTPUT
## ------------------------------------------------------------

out_csv <- file.path(output_dir, "regional_fire_risk_summary.csv")
write_csv(regional_summary, out_csv)

message("✅ Saved: ", out_csv)

head(regional_summary)











##PLOTS

## plot with inferno 
library(tidyverse)
library(viridis)

df <- read_csv(
  file.path(output_dir, "regional_fire_risk_summary.csv"),
  show_col_types = FALSE
)

gcm_models <- c("ACCESS1", "CNRM", "GFDL", "MIROC")
ensemble_model <- "Ensemble"

# --- Reshape to long format ---
plot_df <- df %>%
  dplyr::select(region, model, baseline_mean, future_mean) %>%
  pivot_longer(
    cols = c(baseline_mean, future_mean),
    names_to = "period",
    values_to = "mean_prob"
  ) %>%
  mutate(
    period = recode(period,
                    "baseline_mean" = "Baseline",
                    "future_mean"   = "Future"),
    period = factor(period, levels = c("Baseline", "Future"))
  )

# --- Boxplot data = only the GCM models (no ensemble) ---
box_df <- plot_df %>% filter(model %in% gcm_models)

# --- Inferno colours by region ---
region_list <- unique(plot_df$region)
region_cols <- viridis(length(region_list), option = "inferno")
names(region_cols) <- region_list

# --- Shapes and sizes for models ---
shape_map <- c(
  "ACCESS1" = 15,
  "CNRM"    = 17,
  "GFDL"    = 3,
  "MIROC"   = 8,
  "Ensemble" = 19
)

size_map <- c(
  "ACCESS1" = 2.0,
  "CNRM"    = 2.0,
  "GFDL"    = 2.0,
  "MIROC"   = 2.0,
  "Ensemble" = 4.0
)

# --- B/F labels positioned just above x-axis ---
label_df <- plot_df %>%
  filter(model == "Ensemble") %>%
  mutate(
    label = ifelse(period == "Baseline", "B", "F"),
    y = -0.02     # adjust if needed depending on your y scale
  )

p <- ggplot() +
  
  geom_boxplot(
    data = box_df,
    aes(x = region, y = mean_prob, fill = region,
        group = interaction(region, period)),
    position = position_dodge(width = 0.6),
    alpha = 0.75, outlier.shape = NA, colour = "black"
  ) +
  
  geom_point(
    data = plot_df,
    aes(x = region, y = mean_prob, shape = model, size = model,
        group = interaction(region, period)),
    position = position_dodge(width = 0.6),
    colour = "black", alpha = 0.9
  ) +
  
  # --- B & F LABELS near the x-axis ---
  geom_text(
    data = label_df,
    aes(x = region, y = y, label = label,
        group = interaction(region, period)),
    position = position_dodge(width = 0.6),
    fontface = "bold",
    size = 4,
    vjust = 1
  ) +
  
  scale_fill_manual(values = region_cols, name = "Region") +
  scale_shape_manual(values = shape_map, name = "Climate Model") +
  scale_size_manual(values = size_map, guide = "none") +
  labs(
    x = "Region",
    y = "Mean Probability of Fire",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

p

ggsave(
  filename = file.path(output_dir, "regional_fire_probabilities.png"),
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)



# ============================================================
# PRINT MEAN VALUES FOR TABLE / FIGURE CAPTION
# ============================================================
cat("\n================ MODEL MEANS ================\n")

mean_table <- plot_df %>%
  group_by(model, period) %>%
  summarise(mean_prob = mean(mean_prob), .groups = "drop")

print(mean_table)

cat("=============================================\n")

# --- Overall mean Baseline vs Future across all regions & models ---
overall_period_means <- plot_df %>%
  group_by(period) %>%
  summarise(overall_mean = mean(mean_prob))

cat("\n========= OVERALL PERIOD MEANS (B vs F) =========\n")
print(overall_period_means)
cat("=================================================\n")







#delta plot
delta_df <- df %>%
  mutate(delta = future_mean - baseline_mean)

p_delta <- ggplot(delta_df, aes(x = region, y = delta)) +
  geom_boxplot(fill = "#E69F00", alpha = 0.6, outlier.shape = NA) +
  geom_point(color = "black", size = 2,
             position = position_jitter(width = 0.15)) +
  labs(
    x = "Region",
    y = "Change in Mean Probability",
    title = "Change in Fire Probability (Future – Baseline)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

p_delta














##Predictions RAC MODEL

##read all predictions 
library(terra)
library(ggplot2)
library(tidyterra)
library(viridis)

pred_dir <- "F:/vic_fire_mapping/output_data/predictions"

# === Observed Baseline prediction (1985–2005, RAC) ===
r_base_OBS_RAC <- rast(
  file.path(pred_dir, "fire_observed_baseline_1985_2005_RAC.tif")
)

# === Future predictions (2081–2099) ===
r_fut_ACCESS1  <- rast(file.path(pred_dir, "fire_ACCESS1_0_RCP85_future_2081_2099_RAC.tif"))
r_fut_CNRM     <- rast(file.path(pred_dir, "fire_CNRM_CM5_RCP85_future_2081_2099_RAC.tif"))
r_fut_GFDL     <- rast(file.path(pred_dir, "fire_GFDL_ESM2M_RCP85_future_2081_2099_RAC.tif"))
r_fut_MIROC    <- rast(file.path(pred_dir, "fire_MIROC5_RCP85_future_2081_2099_RAC.tif"))
# Ensemble future
r_fut_ENSEMBLE <- rast(file.path(pred_dir, "fire_Ensemble_Mean_Future_2081_2099_RAC.tif"))

# === Difference maps (Future – Baseline) ===
r_diff_ACCESS1  <- rast(file.path(pred_dir, "fire_change_ACCESS1_0_RCP85_2081_2099_RAC.tif"))
r_diff_CNRM     <- rast(file.path(pred_dir, "fire_change_CNRM_CM5_RCP85_2081_2099_RAC.tif"))
r_diff_GFDL     <- rast(file.path(pred_dir, "fire_change_GFDL_ESM2M_RCP85_2081_2099_RAC.tif"))
r_diff_MIROC    <- rast(file.path(pred_dir, "fire_change_MIROC5_RCP85_2081_2099_RAC.tif"))
# Ensemble difference map
r_diff_ENSEMBLE <- rast(file.path(pred_dir, "fire_change_Ensemble_Mean_2081_2099_RAC.tif"))

summary(r_base_ENSEMBLE)
summary(r_fut_ENSEMBLE)


# ============================================================
# ADD CLIPPING STEPS HERE - Native Vegetation Mask
# ============================================================

library(terra)

# --- Load Native Vegetation Raster ---
veg_path <- "F:/deeca_data/smp/NVR2017_CONDITION.tif"
r_veg <- rast(veg_path)

cat("\nLoaded veg raster:\n")
print(r_veg)
plot(r_veg, main = "Raw Native Vegetation Raster")

# ------------------------------------------------------------
# Reproject vegetaton raster → match prediction raster CRS
# (GDA2020 / VicGrid: use nearest neighbour for categorical data)
# ------------------------------------------------------------
r_veg_proj <- project(r_veg, r_base_OBS_RAC, method = "near")

# Quick inspect values
vals_sample <- unique(values(r_veg_proj)[1:10000])
cat("\nSample veg raster values (after reprojection):\n")
print(vals_sample)

# ------------------------------------------------------------
# Create native vegetation mask
# If NVR2017 is coded: 0 = cleared, 1–100 = native vegetation condition,
# then simply use >0 to keep native vegetation.
# ------------------------------------------------------------
r_veg_mask <- r_veg_proj > 0

plot(r_veg_mask, main = "Native Vegetation Mask (1 = native vegetation)")

align_to_baseline <- function(r) {
  resample(r, r_base_OBS_RAC, method = "bilinear")
}

r_fut_ACCESS1  <- align_to_baseline(r_fut_ACCESS1)
r_fut_CNRM     <- align_to_baseline(r_fut_CNRM)
r_fut_GFDL     <- align_to_baseline(r_fut_GFDL)
r_fut_MIROC    <- align_to_baseline(r_fut_MIROC)
r_fut_ENSEMBLE <- align_to_baseline(r_fut_ENSEMBLE)

# ============================================================
# MASK ALL MODEL RASTERS (Baseline, Future, Delta)
# ============================================================

mask_pred <- function(r) {
  mask(r, r_veg_mask, maskvalues = 0)   # mask out non-native vegetation
}

# --- Observed Baseline (RAC) ---
r_base_OBS_RAC <- mask_pred(r_base_OBS_RAC)

# --- Future ---
r_fut_ACCESS1   <- mask_pred(r_fut_ACCESS1)
r_fut_CNRM      <- mask_pred(r_fut_CNRM)
r_fut_GFDL      <- mask_pred(r_fut_GFDL)
r_fut_MIROC     <- mask_pred(r_fut_MIROC)
r_fut_ENSEMBLE  <- mask_pred(r_fut_ENSEMBLE)

# --- Delta (future – baseline) ---
r_diff_ACCESS1  <- mask_pred(r_diff_ACCESS1)
r_diff_CNRM     <- mask_pred(r_diff_CNRM)
r_diff_GFDL     <- mask_pred(r_diff_GFDL)
r_diff_MIROC    <- mask_pred(r_diff_MIROC)
r_diff_ENSEMBLE <- mask_pred(r_diff_ENSEMBLE)

summary(r_base_OBS_RAC)
plot(r_base_OBS_RAC, main = "Observed Baseline (RAC, Masked to Native Vegetation)")


## ============================================================
## TABLE X: GLOBAL MIN / MAX / MEAN (RAC MODEL)
## ============================================================

library(dplyr)

summarise_pred <- function(r) {
  v <- values(r, na.rm = TRUE)
  tibble(
    min  = min(v),
    max  = max(v),
    mean = mean(v)
  )
}

table_x_rac <- bind_rows(
  
  # --- Observed baseline (RAC) ---
  summarise_pred(r_base_OBS_RAC) %>%
    mutate(period = "Baseline (Observed)", model = "Observed"),
  
  # --- Future GCMs ---
  summarise_pred(r_fut_ACCESS1) %>%
    mutate(period = "Future", model = "ACCESS1"),
  
  summarise_pred(r_fut_GFDL) %>%
    mutate(period = "Future", model = "GFDL"),
  
  summarise_pred(r_fut_MIROC) %>%
    mutate(period = "Future", model = "MIROC"),
  
  summarise_pred(r_fut_CNRM) %>%
    mutate(period = "Future", model = "CNRM"),
  
  # --- Future Ensemble ---
  summarise_pred(r_fut_ENSEMBLE) %>%
    mutate(period = "Future", model = "Ensemble")
) %>%
  select(period, model, min, max, mean) %>%
  mutate(
    period = factor(period, levels = c("Baseline (Observed)", "Future"))
  ) %>%
  arrange(period, model)

print(table_x_rac)


#plot ensembles
# --- Binning scheme (same as before) ---
PROB_BREAKS <- c(0, 0.01, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.75, 0.95)
PROB_LABELS <- c(
  "0–0.01",
  "0.01–0.05",
  "0.05–0.10",
  "0.10–0.15",
  "0.15–0.25",
  "0.25–0.35",
  "0.35–0.50",
  "0.50–0.75",
  "0.75–0.95"
)
PROB_PALETTE <- viridis(
  n = length(PROB_LABELS),
  option = "inferno",
  direction = 1
)

# --- Function to produce the binned map ---
plot_prob_map <- function(r_input, title_text) {
  
  # 1. Classify
  r_binned <- classify(
    r_input,
    rcl = matrix(c(
      PROB_BREAKS[-length(PROB_BREAKS)],
      PROB_BREAKS[-1],
      1:length(PROB_LABELS)
    ), ncol = 3)
  )
  
  # 2. Apply labels
  levels(r_binned) <- data.frame(
    ID = 1:length(PROB_LABELS),
    label = PROB_LABELS
  )
  
  # 3. Plot (no saving yet)
  ggplot() +
    geom_spatraster(data = r_binned) +
    scale_fill_manual(
      values = PROB_PALETTE,
      name = "Probability of fire",
      na.value = "transparent",
      na.translate = FALSE
    ) +
    guides(fill = guide_legend(na.translate = FALSE)) +
    labs(title = title_text) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 14),     
      legend.title = element_text(size = 15),    
      legend.key.size = unit(0.8, "cm")           
    )
}

# --- Make plots ---

p_base_ens <- plot_prob_map(r_base_ENSEMBLE, "")#title removed Ensemble Baseline Prediction
p_base_ens

p_fut_ens <- plot_prob_map(r_fut_ENSEMBLE, "")#removed Ensemble Future Prediction
p_fut_ens





##plot delta map 
names(r_diff_ENSEMBLE) <- "delta_pr"

# 1. Set negative values to zero
r_diff_ENSEMBLE[r_diff_ENSEMBLE < 0] <- 0

# 2. Updated break structure
delta_breaks <- c(0, 0.01, 0.05, 0.10, 0.20, 0.40, 0.60, 0.70, 0.85)

delta_labels <- c(
  "0–0.01",
  "0.01–0.05",
  "0.05–0.10",
  "0.10–0.20",
  "0.20–0.40",
  "0.40–0.60",
  "0.60–0.70",
  "0.70–0.85"
)

# Classification matrix
rcl <- cbind(
  delta_breaks[-length(delta_breaks)],
  delta_breaks[-1],
  1:(length(delta_breaks)-1)
)

r_delta_bin <- classify(r_diff_ENSEMBLE, rcl)

# Apply labels
levels(r_delta_bin) <- data.frame(
  ID = 1:length(delta_labels),
  label = delta_labels
)

# --- Colour palette (light → dark red) ---
pal <- c(
  "#f7f7f7",  # near-zero
  "#FEE5D9",
  "#FCBBA1",
  "#FC9272",
  "#FB6A4A",
  "#DE2D26",
  "#A50F15",
  "#67000D"   # darker for highest bin
)

# --- Plot ---
p_delta <- ggplot() +
  geom_spatraster(data = r_delta_bin, maxcell = 2e6) +
  scale_fill_manual(
    values = pal,
    labels = delta_labels,
    name = expression(paste(Delta, "Pr(Fire)")),
    na.value = "transparent",
    na.translate = FALSE
  ) +
  labs(title = "") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 14),      
    legend.title = element_text(size = 15),     
    legend.key.size = unit(0.8, "cm")          
  )

p_delta



# --- Load fire district shapefile ---
fd_path <- "data/LF_DISTRICT.shp"

regions <- st_read(fd_path, quiet = TRUE) |>
  st_make_valid() |>
  group_by(REGIONNAME) |>
  summarise(.groups = "drop")

# Transform to raster CRS
regions <- st_transform(regions, crs(r_delta_bin))

# Label points (so labels appear inside polygons)
region_labels <- st_point_on_surface(regions["REGIONNAME"])


p_delta_regions <- ggplot() +
  geom_spatraster(data = r_delta_bin, maxcell = 2e6) +
  geom_sf(data = regions, fill = NA, color = "black", linewidth = 0.3) +
  geom_sf_text(
    data = region_labels,
    aes(label = REGIONNAME),
    size = 3,
    color = "black",
    fontface = "plain",
    check_overlap = TRUE
  ) +
  scale_fill_manual(
    values = pal,
    labels = delta_labels,
    name = expression(paste(Delta, "Probability of fire")),
    na.value = "transparent",
    na.translate = FALSE
  ) +
  labs(title = "") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 14),      # bigger legend labels
    legend.title = element_text(size = 15),     # bigger legend title
    legend.key.size = unit(0.8, "cm")           # bigger colour boxes
  )

p_delta_regions



##plot all 3
p_composite <- 
  p_base_ens /
  p_fut_ens /
  p_delta_regions +
  plot_annotation(
    title = "",
    tag_levels = "a",       # adds a), b), c)
    tag_prefix = "",        # nothing before "a"
    tag_suffix = ") ",      # formats as a) b) c)
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(0, 0, 10, 0)
      ),
      plot.tag = element_text(size = 20, face = "bold", hjust = -0.2)  # makes labels larger & shifts left
    )
  )


p_composite

ggsave(
  filename = "outputs/plots/fire_risk_composite_RAC.jpg",
  plot = p_composite,
  width = 10,      # adjust if needed
  height = 18,     # taller because 3 rows
  dpi = 400,
  units = "in",
  device = "jpg"
)








## ============================================================================
## Regional summaries (
## ============================================================================

library(terra)
library(dplyr)
library(readr)
library(janitor)
library(purrr)

## ------------------------------------------------------------
## 1. READ + FIX REGIONS
## ------------------------------------------------------------

fd_path <- "data/LF_DISTRICT.shp"

regions_sf <- st_read(fd_path, quiet = TRUE) |> 
  st_make_valid()

regions_v <- vect(regions_sf)

# Project regions into raster CRS
regions_v <- terra::project(regions_v, crs(r_base_ENSEMBLE))

# Dissolve into 6 regions
regions_v <- aggregate(regions_v, by = "REGIONNAME")

print(regions_v)

## ------------------------------------------------------------
## 2. MODEL–RASTER LIST
## ------------------------------------------------------------

model_rasters <- list(
  ACCESS1 = list(base = r_base_ACCESS1,  fut = r_fut_ACCESS1,  delta = r_diff_ACCESS1),
  CNRM    = list(base = r_base_CNRM,     fut = r_fut_CNRM,     delta = r_diff_CNRM),
  GFDL    = list(base = r_base_GFDL,     fut = r_fut_GFDL,     delta = r_diff_GFDL),
  MIROC   = list(base = r_base_MIROC,    fut = r_fut_MIROC,    delta = r_diff_MIROC),
  Ensemble= list(base = r_base_ENSEMBLE, fut = r_fut_ENSEMBLE, delta = r_diff_ENSEMBLE)
)

ha_per_cell <- (75 * 75) / 10000    # 0.5625 ha

## ------------------------------------------------------------
## 3. NEW FUNCTION: clip each raster to the region FIRST
## ------------------------------------------------------------

summarise_raster <- function(r, region_poly, prefix) {
  
  # crop + mask raster to region
  r_clip <- crop(r, region_poly)
  r_clip <- mask(r_clip, region_poly)
  
  vals <- values(r_clip)
  vals <- vals[!is.na(vals)]
  
  n_cells <- length(vals)
  
  tibble(
    region = region_poly$REGIONNAME,
    !!paste0(prefix, "_n_cells") := n_cells,
    !!paste0(prefix, "_total_area_ha") := n_cells * ha_per_cell,
    !!paste0(prefix, "_mean") := ifelse(n_cells > 0, mean(vals), NA),
    !!paste0(prefix, "_median") := ifelse(n_cells > 0, median(vals), NA),
    !!paste0(prefix, "_prop_gt_0_10") := ifelse(n_cells > 0, mean(vals > 0.10), NA),
    !!paste0(prefix, "_area_gt_0_10_ha") := ifelse(n_cells > 0, sum(vals > 0.10) * ha_per_cell, NA)
  )
}

## ------------------------------------------------------------
## 4. LOOP THROUGH MODELS + REGIONS
## ------------------------------------------------------------

summary_list <- list()

for (m in names(model_rasters)) {
  cat("Processing model:", m, "\n")
  
  base_res  <- map_dfr(1:nrow(regions_v), ~ summarise_raster(model_rasters[[m]]$base,  regions_v[.x], "baseline"))
  fut_res   <- map_dfr(1:nrow(regions_v), ~ summarise_raster(model_rasters[[m]]$fut,   regions_v[.x], "future"))
  delta_res <- map_dfr(1:nrow(regions_v), ~ summarise_raster(model_rasters[[m]]$delta, regions_v[.x], "delta"))
  
  df <- base_res |>
    left_join(fut_res,   by = "region") |>
    left_join(delta_res, by = "region") |>
    mutate(
      model = m,
      relative_change_pct =
        100 * (future_mean - baseline_mean) / pmax(baseline_mean, 0.001)
    )
  
  summary_list[[m]] <- df
}

regional_summary <- bind_rows(summary_list) |>
  arrange(region, model) |>
  clean_names()

## ------------------------------------------------------------
## 5. SAVE OUTPUT
## ------------------------------------------------------------

out_csv <- "outputs/regional_fire_risk_summary_RAC.csv"
write_csv(regional_summary, out_csv)

message("✅ Saved: ", out_csv)

head(regional_summary)






##PLOTS

## plot with inferno 
library(tidyverse)
library(viridis)

df <- read_csv(
  "outputs/regional_fire_risk_summary_RAC.csv",
  show_col_types = FALSE
)

gcm_models <- c("ACCESS1", "CNRM", "GFDL", "MIROC")
ensemble_model <- "Ensemble"

# --- Reshape to long format ---
plot_df <- df %>%
  dplyr::select(region, model, baseline_mean, future_mean) %>%
  pivot_longer(
    cols = c(baseline_mean, future_mean),
    names_to = "period",
    values_to = "mean_prob"
  ) %>%
  mutate(
    period = recode(period,
                    "baseline_mean" = "Baseline",
                    "future_mean"   = "Future"),
    period = factor(period, levels = c("Baseline", "Future"))
  )

# --- Boxplot data = only the GCM models (no ensemble) ---
box_df <- plot_df %>% filter(model %in% gcm_models)

# --- Inferno colours by region ---
region_list <- unique(plot_df$region)
region_cols <- viridis(length(region_list), option = "inferno")
names(region_cols) <- region_list

# --- Shapes and sizes for models ---
shape_map <- c(
  "ACCESS1" = 15,
  "CNRM"    = 17,
  "GFDL"    = 3,
  "MIROC"   = 8,
  "Ensemble" = 19
)

size_map <- c(
  "ACCESS1" = 2.0,
  "CNRM"    = 2.0,
  "GFDL"    = 2.0,
  "MIROC"   = 2.0,
  "Ensemble" = 4.0
)

# --- B/F labels positioned just above x-axis ---
label_df <- plot_df %>%
  filter(model == "Ensemble") %>%
  mutate(
    label = ifelse(period == "Baseline", "B", "F"),
    y = -0.02     # adjust if needed depending on your y scale
  )

p <- ggplot() +
  
  geom_boxplot(
    data = box_df,
    aes(x = region, y = mean_prob, fill = region,
        group = interaction(region, period)),
    position = position_dodge(width = 0.6),
    alpha = 0.75, outlier.shape = NA, colour = "black"
  ) +
  
  geom_point(
    data = plot_df,
    aes(x = region, y = mean_prob, shape = model, size = model,
        group = interaction(region, period)),
    position = position_dodge(width = 0.6),
    colour = "black", alpha = 0.9
  ) +
  
  # --- B & F LABELS near the x-axis ---
  geom_text(
    data = label_df,
    aes(x = region, y = y, label = label,
        group = interaction(region, period)),
    position = position_dodge(width = 0.6),
    fontface = "bold",
    size = 4,
    vjust = 1
  ) +
  
  scale_fill_manual(values = region_cols, name = "Region") +
  scale_shape_manual(values = shape_map, name = "Climate Model") +
  scale_size_manual(values = size_map, guide = "none") +
  labs(
    x = "Region",
    y = "Mean Probability of Fire",
    title = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

p

ggsave(
  filename = "outputs/plots/regional_fire_probabilities_RAC.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)



# ============================================================
# PRINT MEAN VALUES FOR TABLE / FIGURE CAPTION
# ============================================================
cat("\n================ MODEL MEANS ================\n")

mean_table <- plot_df %>%
  group_by(model, period) %>%
  summarise(mean_prob = mean(mean_prob), .groups = "drop")

print(mean_table)

cat("=============================================\n")

# --- Overall mean Baseline vs Future across all regions & models ---
overall_period_means <- plot_df %>%
  group_by(period) %>%
  summarise(overall_mean = mean(mean_prob))

cat("\n========= OVERALL PERIOD MEANS (B vs F) =========\n")
print(overall_period_means)
cat("=================================================\n")







#delta plot
delta_df <- df %>%
  mutate(delta = future_mean - baseline_mean)

p_delta <- ggplot(delta_df, aes(x = region, y = delta)) +
  geom_boxplot(fill = "#E69F00", alpha = 0.6, outlier.shape = NA) +
  geom_point(color = "black", size = 2,
             position = position_jitter(width = 0.15)) +
  labs(
    x = "Region",
    y = "Change in Mean Probability",
    title = "Change in Fire Probability (Future – Baseline)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

p_delta



## ============================================================
## STACKED BOXPLOTS (Maximized Data Area, PNG/PDF Export)
## ============================================================
library(tidyverse)
library(viridis)
library(patchwork)

# --- 1. SETTINGS & PATHS ---
gcm_models     <- c("ACCESS1", "CNRM", "GFDL", "MIROC")
path_rac    <- "outputs/regional_fire_risk_summary_RAC.csv"
path_nonrac <- "outputs/regional_fire_risk_summary.csv"

out_dir <- "outputs/plots"
out_png <- file.path(out_dir, "regional_boxplots_max_size.png")
out_pdf <- file.path(out_dir, "regional_boxplots_max_size.pdf")

# --- 2. THEME & SCALES ---
shape_map <- c("ACCESS1" = 15, "CNRM" = 17, "GFDL" = 3, "MIROC" = 8, "Ensemble" = 21)
size_map  <- c("ACCESS1" = 2.5, "CNRM" = 2.5, "GFDL" = 2.5, "MIROC" = 2.5, "Ensemble" = 3.5)

# --- 3. PLOT FUNCTION ---
create_regional_boxplot <- function(file_path, show_x = TRUE) {
  
  df_raw <- read_csv(file_path, show_col_types = FALSE)
  if (!"model" %in% names(df_raw)) { df_raw <- df_raw %>% mutate(model = "Ensemble") }
  
  plot_df <- df_raw %>%
    dplyr::select(region, model, baseline_mean, future_mean) %>%
    pivot_longer(cols = c(baseline_mean, future_mean), names_to = "period", values_to = "mean_prob") %>%
    mutate(
      period = recode(period, "baseline_mean" = "Baseline", "future_mean" = "Future"),
      period = factor(period, levels = c("Baseline", "Future"))
    )
  
  box_df <- plot_df %>% filter(model %in% gcm_models)
  # B/F text placement
  label_df <- plot_df %>% filter(model == "Ensemble") %>% mutate(label = ifelse(period == "Baseline", "B", "F"), y = -0.008) 
  
  region_cols <- viridis(length(unique(plot_df$region)), option = "inferno")
  
  p <- ggplot() +
    {if(nrow(box_df) > 0) geom_boxplot(
      data = box_df,
      aes(x = region, y = mean_prob, fill = region, group = interaction(region, period)),
      position = position_dodge(width = 0.6),
      alpha = 0.7, outlier.shape = NA, colour = "black"
    )} +
    geom_point(
      data = plot_df,
      aes(x = region, y = mean_prob, shape = model, size = model, group = interaction(region, period)),
      position = position_dodge(width = 0.6),
      fill = "white", colour = "black", alpha = 0.9
    ) +
    geom_text(
      data = label_df,
      aes(x = region, y = y, label = label, group = interaction(region, period)),
      position = position_dodge(width = 0.6),
      fontface = "bold", size = 5.5, vjust = 1.2 # Tighter vjust
    ) +
    scale_fill_manual(values = region_cols, guide = "none") + 
    scale_shape_manual(values = shape_map, name = "Climate Model") +
    scale_size_manual(values = size_map, guide = "none") +
    # TIGHTER expansion: c(0.15, 0.05) stretches the boxplots more vertically
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    labs(x = if(show_x) "Region" else NULL, y = "Mean Pr(Fire)") +
    theme_bw(base_size = 18) +
    theme(
      plot.tag = element_text(face = "bold", size = 24, margin = margin(b = 10)),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(t = 5, r = 10, b = 5, l = 10)
    ) +
    guides(shape = guide_legend(
      title.position = "top", 
      title.hjust = 0, 
      nrow = 1,
      override.aes = list(size = 5)
    ))
  
  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
                   plot.margin = margin(t = 5, r = 10, b = -10, l = 10)) # Negative bottom margin
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
                   plot.margin = margin(t = -10, r = 10, b = 5, l = 10)) # Negative top margin
  }
  
  return(p)
}

# --- 4. CREATE ---
p_nonrac <- create_regional_boxplot(path_nonrac, show_x = FALSE)
p_rac    <- create_regional_boxplot(path_rac, show_x = TRUE)

# --- 5. ASSEMBLE ---
composite_plot <- wrap_plots(p_nonrac, p_rac, ncol = 1) +
  plot_layout(guides = "collect", heights = c(1, 1)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") & 
  theme(
    legend.position = "bottom",
    legend.box.just = "left",
    legend.justification = "left",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(t = 5, r = 10, b = 60, l = 10) 
  )

# --- 6. SAVE ---
# Width remains narrow (10), Height increased (15) to stretch plots further
ggsave(filename = out_png, plot = composite_plot, width = 10, height = 15, dpi = 300, scale = 0.8)
ggsave(filename = out_pdf, plot = composite_plot, width = 10, height = 15, device = cairo_pdf, scale = 0.8)

message("Success! Maximized boxplots saved to: ", out_dir)







##Cleaned Fig 6 and 7 Script

library(terra)
library(sf)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(tidyterra)
library(viridis)
library(janitor)

# ============================================================
# 1. PATHS
# ============================================================

pred_dir <- "outputs/predictions"
veg_path <- "data/NVR2017_CONDITION.tif"
fd_path  <- "data/LF_DISTRICT.shp"

out_dir  <- "outputs/plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. NATIVE VEGETATION MASK (ONCE)
# ============================================================

r_veg <- rast(veg_path)

make_veg_mask <- function(template_raster) {
  r_veg_proj <- project(r_veg, template_raster, method = "near")
  r_veg_proj > 0
}

mask_to_veg <- function(r, veg_mask) {
  mask(r, veg_mask, maskvalues = 0)
}

# ============================================================
# 3. REGIONS (ONCE)
# ============================================================

regions_sf <- st_read(fd_path, quiet = TRUE) |>
  st_make_valid() |>
  group_by(REGIONNAME) |>
  summarise(.groups = "drop")

regions_v <- vect(regions_sf)

# ============================================================
# 4. LOAD + PREP RASTERS
# ============================================================

load_model_stack <- function(type = c("nonRAC", "RAC")) {
  type <- match.arg(type)
  suf  <- ifelse(type == "RAC", "_RAC", "")
  
  base <- rast(file.path(pred_dir,
                         paste0("fire_observed_baseline_1985_2005", suf, ".tif")
  ))
  
  fut <- list(
    ACCESS1  = rast(file.path(pred_dir, paste0("fire_ACCESS1_0_RCP85_future_2081_2099", suf, ".tif"))),
    CNRM     = rast(file.path(pred_dir, paste0("fire_CNRM_CM5_RCP85_future_2081_2099", suf, ".tif"))),
    GFDL     = rast(file.path(pred_dir, paste0("fire_GFDL_ESM2M_RCP85_future_2081_2099", suf, ".tif"))),
    MIROC    = rast(file.path(pred_dir, paste0("fire_MIROC5_RCP85_future_2081_2099", suf, ".tif"))),
    Ensemble = rast(file.path(pred_dir, paste0("fire_Ensemble_Mean_Future_2081_2099", suf, ".tif")))
  )
  
  # Align
  fut <- map(fut, ~ resample(.x, base, method = "bilinear"))
  
  veg_mask <- make_veg_mask(base)
  
  base <- mask_to_veg(base, veg_mask)
  fut  <- map(fut, mask_to_veg, veg_mask = veg_mask)
  
  delta <- map(fut, ~ .x - base)
  
  list(
    type  = type,
    base  = base,
    fut   = fut,
    delta = delta,
    veg_mask = veg_mask
  )
}

# ============================================================
# 5. GLOBAL SUMMARY
# ============================================================

summarise_raster <- function(r) {
  v <- values(r, na.rm = TRUE)
  tibble(
    min  = min(v),
    max  = max(v),
    mean = mean(v)
  )
}

global_summary <- function(stack) {
  
  bind_rows(
    summarise_raster(stack$base) |>
      mutate(period = "Baseline (Observed)", model = "Observed"),
    
    imap_dfr(stack$fut, ~
               summarise_raster(.x) |>
               mutate(period = "Future", model = .y)
    )
  ) |>
    mutate(model_type = stack$type) |>
    select(model_type, period, model, min, max, mean)
}


# ============================================================
# 6. REGIONAL SUMMARY (MIN / MAX / MEAN)
#   - Baseline: Observed only
#   - Future: ALL models (ACCESS1, CNRM, GFDL, MIROC, Ensemble)
# ============================================================

regional_summary <- function(stack, regions_v) {
  
  regions_v <- terra::project(regions_v, crs(stack$base))
  
  summarise_region <- function(r, region) {
    r_clip <- crop(r, region) |> mask(region)
    v <- values(r_clip, na.rm = TRUE)
    
    tibble(
      min  = min(v),
      max  = max(v),
      mean = mean(v)
    )
  }
  
  map_dfr(seq_len(nrow(regions_v)), function(i) {
    
    reg <- regions_v[i]
    
    # --- Baseline (Observed only) ---
    baseline_df <- summarise_region(stack$base, reg) |>
      mutate(
        period = "Baseline",
        model  = "Observed"
      )
    
    # --- Future (ALL models) ---
    future_df <- imap_dfr(stack$fut, function(r_fut, model_name) {
      summarise_region(r_fut, reg) |>
        mutate(
          period = "Future",
          model  = model_name
        )
    })
    
    bind_rows(baseline_df, future_df) |>
      mutate(
        region     = reg$REGIONNAME,
        model_type = stack$type
      )
  })
}


# ============================================================
# 7. RUN PIPELINE
# ============================================================

stack_nonrac <- load_model_stack("nonRAC")
stack_rac    <- load_model_stack("RAC")

global_table <- bind_rows(
  global_summary(stack_nonrac),
  global_summary(stack_rac)
)

regional_table <- bind_rows(
  regional_summary(stack_nonrac, regions_v),
  regional_summary(stack_rac, regions_v)
)

# ============================================================
# 8. SAVE OUTPUTS
# ============================================================

write_csv(global_table,
          file.path(out_dir, "global_fire_probability_summary.csv")
)

write_csv(regional_table,
          file.path(out_dir, "regional_fire_probability_summary.csv")
)





## ============================================================
## STACKED BARPLOTS (Observed Baseline + Future GCMs) 
## ============================================================

library(tidyverse)
library(viridis)
library(patchwork)

# ------------------------------------------------------------
# 1. SETTINGS & PATHS
# ------------------------------------------------------------

gcm_models <- c("ACCESS1", "CNRM", "GFDL", "MIROC")

path_regional <- "outputs/regional_fire_probability_summary.csv"

out_dir <- "outputs/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(out_dir, "regional_barplots_observed_baseline.png")
out_pdf <- file.path(out_dir, "regional_barplots_observed_baseline.pdf")

# ------------------------------------------------------------
# 2. THEME & SCALES (UNCHANGED)
# ------------------------------------------------------------

shape_map <- c(
  "ACCESS1"  = 15,
  "CNRM"     = 17,
  "GFDL"     = 3,
  "MIROC"    = 8,
  "Ensemble" = 21
)

size_map <- c(
  "ACCESS1"  = 2.5,
  "CNRM"     = 2.5,
  "GFDL"     = 2.5,
  "MIROC"    = 2.5,
  "Ensemble" = 3.5
)

# ------------------------------------------------------------
# 3. PLOT FUNCTION (FIXED: explicit left/right bar positions)
# ------------------------------------------------------------

create_regional_barplot <- function(df_all, model_type_use, show_x = TRUE) {
  
  # --- Prep ---
  plot_df <- df_all %>%
    filter(model_type == model_type_use) %>%
    mutate(
      period    = factor(period, levels = c("Baseline", "Future")),
      mean_prob = mean,
      # numeric x so we can control left/right placement precisely
      region_f  = factor(region, levels = unique(region)),
      region_id = as.numeric(region_f)
    )
  
  # Offsets for left (Baseline) and right (Future)
  x_offset <- 0.22
  bar_width <- 0.38
  
  # --- Bars ---
  # Baseline bar: observed only
  bar_baseline <- plot_df %>%
    filter(period == "Baseline", model == "Observed") %>%
    mutate(x = region_id - x_offset)
  
  # Future bar: ensemble mean only
  bar_future <- plot_df %>%
    filter(period == "Future", model == "Ensemble") %>%
    mutate(x = region_id + x_offset)
  
  bar_df <- bind_rows(bar_baseline, bar_future)
  
  # --- Points ---
  # Future climate models (incl Ensemble) -> ON THE RIGHT BAR ONLY
  point_future <- plot_df %>%
    filter(period == "Future") %>%
    mutate(x = region_id + x_offset)
  
  # --- B / F labels (based on bars) ---
  label_df <- bar_df %>%
    mutate(
      label = ifelse(period == "Baseline", "B", "F"),
      y = -0.008
    )
  
  # Region colours (inferno)
  region_cols <- viridis(length(unique(plot_df$region_f)), option = "inferno")
  names(region_cols) <- levels(plot_df$region_f)
  
  # --- Plot ---
  p <- ggplot() +
    
    # Bars: left baseline + right future
    geom_col(
      data = bar_df,
      aes(x = x, y = mean_prob, fill = region_f),
      width = bar_width,
      colour = "black",
      alpha = 0.8
    ) +
    
    # Future model symbols: only on right bar
    geom_point(
      data = point_future,
      aes(x = x, y = mean_prob, shape = model, size = model),
      fill = "white",
      colour = "black",
      alpha = 0.9
    ) +
    
    # B/F labels beneath bars (kept identical style)
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      fontface = "bold",
      size = 5.5,
      vjust = 1.2
    ) +
    
    scale_fill_manual(values = region_cols, guide = "none") +
    scale_shape_manual(values = shape_map, name = "Climate Model") +
    scale_size_manual(values = size_map, guide = "none") +
    
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    
    # numeric x-axis but labelled as regions (matches old plot appearance)
    scale_x_continuous(
      breaks = unique(plot_df$region_id),
      labels = levels(plot_df$region_f)
    ) +
    
    labs(
      x = if (show_x) "Region" else NULL,
      y = "Mean Pr(Fire)"
    ) +
    
    theme_bw(base_size = 18) +
    theme(
      plot.tag = element_text(face = "bold", size = 24, margin = margin(b = 10)),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(t = 5, r = 10, b = 5, l = 10)
    ) +
    guides(
      shape = guide_legend(
        title.position = "top",
        title.hjust = 0,
        nrow = 1,
        override.aes = list(size = 5)
      )
    )
  
  if (!show_x) {
    p <- p +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin  = margin(t = 5, r = 10, b = -10, l = 10)
      )
  } else {
    p <- p +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.margin = margin(t = -10, r = 10, b = 5, l = 10)
      )
  }
  
  return(p)
}

# ------------------------------------------------------------
# 4. LOAD DATA
# ------------------------------------------------------------

regional_df <- read_csv(path_regional, show_col_types = FALSE)

# ------------------------------------------------------------
# 5. CREATE PLOTS
# ------------------------------------------------------------

p_nonrac <- create_regional_barplot(regional_df, "nonRAC", show_x = FALSE)
p_rac    <- create_regional_barplot(regional_df, "RAC",    show_x = TRUE)

# ------------------------------------------------------------
# 6. ASSEMBLE
# ------------------------------------------------------------

composite_plot <- wrap_plots(p_nonrac, p_rac, ncol = 1) +
  plot_layout(guides = "collect", heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(
    legend.position      = "bottom",
    legend.box.just      = "left",
    legend.justification = "left",
    legend.text          = element_text(size = 14),
    legend.title         = element_text(size = 16, face = "bold"),
    plot.margin          = margin(t = 5, r = 10, b = 60, l = 10)
  )

composite_plot

# ------------------------------------------------------------
# 7. SAVE
# ------------------------------------------------------------

ggsave(out_png, composite_plot, width = 10, height = 15, dpi = 300, scale = 0.8)
ggsave(out_pdf, composite_plot, width = 10, height = 15, device = cairo_pdf, scale = 0.8)

message("✅ Observed-baseline regional barplots saved to: ", out_dir)











### with numbers for regions
## ============================================================
## 0. PACKAGES
## ============================================================
library(terra)
library(ggplot2)
library(tidyterra)
library(sf)
library(patchwork)
library(dplyr)
library(scales)
library(grid)

## ============================================================
## 1. FILE PATHS
## ============================================================
pred_dir <- "outputs/predictions"
veg_path <- "data/NVR2017_CONDITION.tif"
fd_path  <- "data/LF_DISTRICT.shp"

out_dir  <- "C:/Users/charliehart/OneDrive - The University of Melbourne/fire_modelling/plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png  <- file.path(out_dir, "fire_risk_composite_continuous_numbered_v4.png")
out_pdf  <- file.path(out_dir, "fire_risk_composite_continuous_numbered_v4.pdf")

## ============================================================
## 2. LOAD RASTERS (OBSERVED BASELINE + FUTURE)
## ============================================================

# --- Observed baseline ---
nr_base  <- rast(file.path(pred_dir, "fire_observed_baseline_1985_2005.tif"))
rac_base <- rast(file.path(pred_dir, "fire_observed_baseline_1985_2005_RAC.tif"))

# --- Future ensemble ---
nr_fut  <- rast(file.path(pred_dir, "fire_Ensemble_Mean_Future_2081_2099.tif"))
rac_fut <- rast(file.path(pred_dir, "fire_Ensemble_Mean_Future_2081_2099_RAC.tif"))

## ============================================================
## 2b. CLIP TO SHARED EXTENT (NO RESAMPLING)  ✅ NEW
## ============================================================

# --- non-RAC ---
ext_nr <- terra::intersect(
  terra::ext(nr_base),
  terra::ext(nr_fut)
)

nr_base <- terra::crop(nr_base, ext_nr)
nr_fut  <- terra::crop(nr_fut,  ext_nr)

# --- RAC ---
ext_rac <- terra::intersect(
  terra::ext(rac_base),
  terra::ext(rac_fut)
)

rac_base <- terra::crop(rac_base, ext_rac)
rac_fut  <- terra::crop(rac_fut,  ext_rac)

## ============================================================
## 2c. DELTA: FUTURE − OBSERVED BASELINE  ✅ UPDATED
## ============================================================

nr_delta  <- nr_fut  - nr_base
rac_delta <- rac_fut - rac_base

## ============================================================
## 3. MASK TO NATIVE VEGETATION
## ============================================================
r_veg <- rast(veg_path)
r_veg_proj <- project(r_veg, nr_base, method = "near")
veg_mask   <- r_veg_proj > 0

mask_pred <- function(r) {
  mask(r, veg_mask, maskvalues = 0)
}

nr_base  <- mask_pred(nr_base)
nr_fut   <- mask_pred(nr_fut)
nr_delta <- mask_pred(nr_delta)

rac_base  <- mask_pred(rac_base)
rac_fut   <- mask_pred(rac_fut)
rac_delta <- mask_pred(rac_delta)

## ============================================================
## 4. REGIONS & CENTROIDS
## ============================================================
regions <- st_read(fd_path, quiet = TRUE) |>
  st_make_valid() |>
  group_by(REGIONNAME) |>
  summarise(.groups = "drop") |>
  mutate(region_id = row_number())

regions <- st_transform(regions, crs(nr_delta))
region_labels <- st_centroid(regions)

## ============================================================
## 5. LIMITS & BREAKS
## ============================================================
prob_limits <- c(0, 1)
prob_breaks <- seq(0, 1, by = 0.2)
prob_labels <- number_format(accuracy = 0.1)(prob_breaks)

nr_min <- min(
  terra::global(nr_delta, "min", na.rm = TRUE)[1,1],
  terra::global(rac_delta, "min", na.rm = TRUE)[1,1]
)

nr_max <- max(
  terra::global(nr_delta, "max", na.rm = TRUE)[1,1],
  terra::global(rac_delta, "max", na.rm = TRUE)[1,1]
)

delta_limits <- c(nr_min, nr_max)
delta_breaks <- pretty(delta_limits, n = 6)
delta_labels <- number_format(accuracy = 0.1)(delta_breaks)

## ============================================================
## 6. THEMES & GUIDES
## ============================================================
guide_prob <- guide_colorbar(
  direction = "horizontal",
  title.position = "top",
  barwidth = unit(18, "cm"),
  barheight = unit(0.8, "cm")
)

guide_delta <- guide_colorbar(
  direction = "horizontal",
  title.position = "top",
  barwidth = unit(18, "cm"),
  barheight = unit(0.8, "cm")
)

base_theme_map <- theme_void() +
  theme(
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 20)
  )

## ============================================================
## 7. PLOT FUNCTIONS
## ============================================================
plot_prob <- function(r) {
  ggplot() +
    geom_spatraster(data = r, na_color = "white") +
    scale_fill_viridis_c(
      option = "inferno",
      limits = prob_limits,
      oob = scales::squish,
      breaks = prob_breaks,
      labels = prob_labels,
      name = "Pr(Fire)",
      na.value = "white"
    ) +
    guides(fill = guide_prob) +
    base_theme_map
}

plot_delta <- function(r) {
  ggplot() +
    geom_spatraster(data = r, na_color = "white") +
    geom_sf(data = regions, fill = NA, colour = "black", linewidth = 0.3) +
    geom_sf_text(
      data = region_labels,
      aes(label = region_id),
      size = 6,
      fontface = "bold",
      colour = "black"
    ) +
    scale_fill_gradient(
      low = "#FEE5D9",
      high = "#A50F15",
      limits = delta_limits,
      oob = scales::squish,
      breaks = delta_breaks,
      labels = delta_labels,
      name = "Change in Pr(Fire)",
      na.value = "white"
    ) +
    guides(fill = guide_delta) +
    base_theme_map
}

## ============================================================
## 8. BUILD PANELS
## ============================================================
p_a <- plot_prob(nr_base)
p_b <- plot_prob(rac_base)
p_c <- plot_prob(nr_fut)
p_d <- plot_prob(rac_fut)
p_e <- plot_delta(nr_delta)
p_f <- plot_delta(rac_delta)

## ============================================================
## 9–11. ASSEMBLY (UNCHANGED)
## ============================================================
block_prob <- wrap_plots(p_a, p_b, p_c, p_d, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 30, face = "bold"),
    plot.margin = margin(t = 15, r = 10, b = 5, l = 10)
  )

block_delta <- wrap_plots(p_e, p_f, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 30, face = "bold"),
    plot.margin = margin(t = 15, r = 10, b = 5, l = 10)
  )

final_plot <- wrap_plots(block_prob, block_delta, ncol = 1) +
  plot_layout(heights = c(2.2, 1.2)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(aspect.ratio = NULL)

## ============================================================
## 12. SAVE
## ============================================================
ggsave(out_png, final_plot, width = 14, height = 18, units = "in", dpi = 300)
ggsave(out_pdf, final_plot, width = 14, height = 18, units = "in", device = cairo_pdf)

out_png
