library(terra)
library(sf)
library(dplyr)

#### Make analysis mask
mask_path <- here("data", "HDM_10253_Sooty_Owl_75m.tif")
vic_path  <- here("data", "STE_2021_AUST_GDA2020.shp")
mask = rasterize(vic, mask)

names(mask) <- "victoria_mask"
varnames(mask) <- "victoria_mask"

writeRaster(mask, here("covariates", "victoria_mask_75m.tif"))

#### Make random points for study region
points = terra::spatSample(mask, size = 10000, method = "stratified", as.points = TRUE)

writeVector(points, here("covariates", "random_points.shp"))