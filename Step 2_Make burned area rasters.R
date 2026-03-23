# Step 2 - Create annual burned areas
library(terra)
library(dplyr)
library(tidyr)
library(sf)
library(here)

#### Maksf#### Make Time Since Fire Map ####
fire_history_path <- here("data", "FIRE_HISTORY.shp")
mask_path <- here("covariates", "victoria_mask_75m.tif")

firehistoryT = firehistory %>% sf::st_transform(crs(mask))
## Season is the last year in the fire season. 
### Loop through each fire SEASON to create raster of fires
seasons = sort(unique(firehistory$SEASON))
for (i in 1:length(seasons)){
  print(i)
  s = seasons[i]
  fires = firehistoryT %>% filter(SEASON == s) %>% vect()
  fire.ras = rasterize(fires, mask, field = "SEASON", background = NA)
  fire.ras <- ifel(fire.ras > 1, 1, 0)
  fire.ras[is.na(fire.ras)] <-0
  names(fire.ras) <- paste0("burned_area_",s)
  out_path <- here("fire_data", paste0("burned_area_", s, "_75m.tif"))
  writeRaster(fire_ras, out_path, overwrite = TRUE)
}

# Create season start and end file
dates = firehistory %>% st_drop_geometry() %>%
  select(SEASON, START_DATE) %>% distinct() %>% 
  filter(SEASON > 1980) %>% drop_na() #%>%
  group_by(SEASON) %>% summarise(mindate = min(START_DATE),
                                 maxdate = max(START_DATE))

