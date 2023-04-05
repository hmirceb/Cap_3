library(climaemet) # meteorological data
library(mapSpain) # base maps of Spain
library(classInt) # classification
library(terra) # raster handling
library(sf) # spatial shape handling
library(gstat) # for spatial interpolation
library(geoR) # for spatial analysis
library(tidyverse) # collection of R packages designed for data science
library(tidyterra) # tidyverse methods for terra package
rm(list = ls())

stations <- aemet_stations()

# Select data
date_start <- 1981
date_end <- 2022
estaciones <- unique(stations$indicativo)
clim_data_list <- c()
for (est in seq_along(estaciones)) {
  tryCatch({
  clim_data_list[[est]] <- aemet_monthly_period(station = estaciones[[est]],
                       start = date_start,
                       end = date_end,
                       return_sf = T,
                       verbose = T)
  print(paste0(100*est/length(estaciones), '%'))},
  error=function(e){}
  )
}
save(clim_data_list, file = 'Datos_AEMET_1981_2022.RData')