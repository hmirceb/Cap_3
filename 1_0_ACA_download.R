library(tidyverse)
library(terra)
library(readxl)
library(sp)
rm(list=ls())
graphics.off()
setwd("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")
setwd("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")

# Cargamos funciones para descargar el Atlas Climatico de Aragon
source('Scripts/ACA_Aragon.R')

# Cargamos datos de MONITO
load("Datos/MONITO_2022.RData")

# Cargamos coordenadas de MUs y convertimos a formato SpatVector para usar luego con terra
MAIN_DETAILED <- read_excel("~/Dropbox/DATA__LAB/__MONITO/MONITO_ADOPTA/5_DATOS & ANÁLISIS/MONITO_FILEMAKER/FM_EXPORTED/Main_Tables_090123/MAIN_DETAILED.xlsx")
mus_coords<-terra::vect(SpatialPoints(cbind(MAIN_DETAILED$X30_DET, MAIN_DETAILED$Y30_DET),
                          proj4string = CRS("+init=epsg:25830")))

# Destino donde guardar los datos descargados
dest_promedio <- '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/Datos/Clima/Promedios/'

#########################################################################
# Empezamos por descargar los valores promedio anuales de cada variable #
#########################################################################
vars <- c('pcp', 'tmax', 'tmed', 'tmin')
for (var in seq_along(vars)) {
  ACA_promedios(periodo = 'anual', x = vars[[var]],
                dest = dest_promedio)
}

# Cargamos los archivos raster
pcp_raster_annual <- rast(paste0(dest_promedio, list.files(dest_promedio, pattern = 'pcp')))
tmax_raster_annual <- rast(paste0(dest_promedio, list.files(dest_promedio, pattern = 'tmax')))
tmed_raster_annual <- rast(paste0(dest_promedio, list.files(dest_promedio, pattern = 'tmed')))
tmin_raster_annual <- rast(paste0(dest_promedio, list.files(dest_promedio, pattern = 'tmin')))

# Extraemos los valores de cada MU y renombramos
mus_temps <- cbind(terra::extract(pcp_raster_annual, mus_coords, ID = F), 
                   terra::extract(tmax_raster_annual, mus_coords, ID = F),
                   terra::extract(tmed_raster_annual, mus_coords, ID = F),
                   terra::extract(tmin_raster_annual, mus_coords, ID = F),
                   as.data.frame(mus_coords, geom = 'xy'))
names(mus_temps)<- c('pcp_anual', 'tmax_anual', 'tmed_anual', 'tmin_anual', 'X30_DET', 'Y30_DET')
mus_temps <- MAIN_DETAILED %>%
  inner_join(mus_temps) %>%
  distinct()
save(mus_temps, file = 'Datos/Promedios_MU.RData')

#################################################
# Descargamos valores de cada anyo para cada MU #
#################################################

# Establecemos los anyos de cada MU y quitamos los posteriores a 2020
# porque no se incluyen en el ACA
dest_brutos <- '~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/Datos/Clima/Brutos/'

years_mu <- MONITO_DATA %>%
  dplyr::select(c(MU_DYN, Year)) %>%
  filter(Year <= 2020) %>%
  group_by(MU_DYN) %>%
  summarise(Inicio = min(Year),
            Final = max(Year))

MAIN_DETAILED <- MAIN_DETAILED %>%
  mutate(X30_ACA = as.integer(
    paste0(substr(as.character(round(X30_DET, 0)), 
                  start = 1, 
                  stop = nchar(as.character(round(X30_DET, 0)))-3), 
           '500')),
    Y30_ACA = as.integer(
      paste0(substr(as.character(round(Y30_DET, 0)), 
                            start = 1, 
                            stop = nchar(as.character(round(Y30_DET, 0)))-3), 
                     '500'))) 
mus_coords <- MAIN_DETAILED %>%
  dplyr::select(c(X30_ACA, Y30_ACA)) %>%
  distinct() 

errores <- c()
for (mu in seq_len(dim(mus_coords)[1])) {
  ACA_coord(lon = mus_coords$X30_ACA[[mu]], 
            lat = mus_coords$Y30_ACA[[mu]],
            dest = dest_brutos)
  print(paste(mus_coords$X30_ACA[[mu]], mus_coords$Y30_ACA[[mu]], sep = '_'))
  errores[[mu]] <- list.files(dest_brutos, pattern = "\\.zip")[mu]
}

zip_list <- list.files(dest_brutos, pattern = "\\.zip")
coord_abiotic <- c()
for (coord in seq_along(zip_list)) {
   t <- read.table(unz(paste0(dest_brutos, zip_list[coord]),
                      paste0('csv/', paste0(substr(zip_list[coord], 
                                                   start = 1, 
                                                   stop = nchar(zip_list[coord])-3), 'csv'))), 
                  header=T, dec = ',', sep=";")
   t$X30_ACA <- as.numeric(unlist(lapply(strsplit(zip_list[coord], '_'), function(x) x[1])))
   t$Y30_ACA <-  unlist(lapply(strsplit(zip_list[coord], '_'), function(x) x[2]))
   t$Y30_ACA <- as.numeric(substr(t$Y30_ACA, 1, nchar(t$Y30_ACA)-4)) # Quitamos el .zip
  coord_abiotic[[coord]] <- t
}
coord_abiotic <- do.call('rbind', coord_abiotic)
coord_abiotic$Year <- unlist(lapply(strsplit(coord_abiotic$fecha, '-'), function(x) x[1]))
coord_abiotic$Month <- unlist(lapply(strsplit(coord_abiotic$fecha, '-'), function(x) x[2]))
coord_abiotic$Day <- unlist(lapply(strsplit(coord_abiotic$fecha, '-'), function(x) x[3]))

coord_abiotic <- coord_abiotic %>%
  inner_join((MAIN_DETAILED %>%
  mutate(X30_ACA = as.numeric(
    paste0(substr(as.character(X30_DET), 
                  start = 1, 
                  stop = nchar(as.character(X30_DET))-3), 
           '500')),
    Y30_ACA = as.numeric(
      paste0(substr(as.character(Y30_DET), 
                    start = 1, 
                    stop = nchar(as.character(Y30_DET))-3), 
             '500'))) %>%
  dplyr::select(c(MU_DYN, X30_ACA, Y30_ACA))),
  by = c('X30_ACA' = 'X30_ACA', 'Y30_ACA' = 'Y30_ACA')) %>%
  distinct()

save(coord_abiotic, file = 'Datos/Diarios_MU.RData')

coord_abiotic %>%
  group_by(MU_DYN, Year) %>%
  summarise(pcp = mean(pcp),
            tmax = mean(tmedmax),
            tmin = mean(tmedmin)) %>%
  ggplot(aes(x = as.numeric(Year), y = pcp))+
  theme_classic()+
  geom_line()+
  geom_smooth()+
  facet_wrap(~MU_DYN)

