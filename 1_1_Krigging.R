library(climaemet)
library(mapSpain) # Base maps of Spain
library(sf) # spatial shape handling
library(terra) # Spatial raster handling
library(gstat) # for spatial interpolation
library(tidyterra)
library(geoR)
library(tidyverse)
library(foreach)
library(automap)
library(stars)
library(readxl)
library(sp)
rm(list = ls())

setwd("/home/hector/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")

#---------------
# Cargamos datos
#---------------
# Establecemos el EPSG general (para no tener que estar todo el rato repitiendo)
epsg <- "epsg:25830"

# Datos de AEMET
load("~/MONITO/Data/Datos_AEMET_1981_2022.RData")

# Posicion de las MUs
MAIN_DETAILED <- read_excel("~/MONITO/Data/MAIN_DETAILED.xlsx")
# Creamos una capa de puntos de las MUs en SF
mus_coords<-st_as_sf(MAIN_DETAILED, 
                     coords = c('X30_DET', 'Y30_DET')) %>%
  st_set_crs(epsg)
# Estaciones de la AEMET
stations <- aemet_stations()

# Juntamos datos de AEMET y sus estaciones 
# y nos quedamos solo con los meses (el mes 13 indica el anyo)
clim_data_all <- clim_data_list %>%
  bind_rows() %>%
  mutate(Year = unlist(lapply(strsplit(fecha, '-'), function(x) x[1])),
         Month = unlist(lapply(strsplit(fecha, '-'), function(x) x[2]))) %>%
  left_join(stations, by = 'indicativo') %>%
  filter(Month != 13) %>%
  dplyr::select(indicativo, nombre, provincia, altitud, longitud, latitud,
                Year, Month, p_mes, tm_mes, ta_min, tm_min, ta_max, tm_max,
                geometry) %>%
  mutate(ta_min_day = str_remove_all(unlist(lapply(strsplit(ta_min, '\\('), function(x) x[2])), '\\)'),
         ta_max_day = str_remove_all(unlist(lapply(strsplit(ta_max, '\\('), function(x) x[2])), '\\)'),
         ta_min = unlist(lapply(strsplit(ta_min, '\\('), function(x) x[1])),
         ta_max = unlist(lapply(strsplit(ta_max, '\\('), function(x) x[1]))) %>%
  mutate(across(p_mes:tm_max, as.numeric)) %>%
  filter(!provincia %in% c("LAS PALMAS", "STA. CRUZ DE TENERIFE")) %>%
  st_transform(crs = epsg)
rm(clim_data_list, stations)

#--------------
# Grid de 1x1km
#--------------
# Descargamos contorno de las comunidades autonomas
CCAA <- esp_get_ccaa(epsg = 4326) %>%
  # Quitamos Canarias y Baleares
  filter(!ine.ccaa.name %in% c("Canarias")) %>%
  st_transform(crs = epsg)

# Creamos el grid
side <- 1000 # Lado del cuadrado en las unidades de la capa CCAA
grd <- st_as_stars(rast(CCAA, res = c(side, side)))

#-------------------------------------
# ALTERNATIVA: Creamos un DEM de 1x1km
#-------------------------------------
# Se ha usado el DEM del IGN con resolucion de 200m (porque total los iba a juntar a 1km y asi ocupaban menos)
# Leemos los datos y nos quedamos solo con los que estan en el huso 30 (algunas imagenes estan reptidos en varios usos)
f <- list.files(path = "~/MONITO/Data/IGN",pattern = '.*\\.(tif|asc)$') 
f <- f[grepl(pattern = 'HU30', x = f)] # Huso 30
r <- lapply(lapply(f, function(x) paste0('~/MONITO/Data/IGN/', x)), rast)

# Como vienen sin CRS se lo ponemos y de paso reducimos la resolucion de los datos
# a 1x1km
z<-c()
for (ras in seq_along(r)) {
  crs(r[[ras]]) <- epsg
  z[[ras]]<-terra::aggregate(r[[ras]], fact = 5, fun = 'mean')
}
# Juntamos todo en una sola imagen
dem1000 <- do.call('merge', z)
names(dem1000)<-'altitud'
rm(r, z, f, ras)

# Transformamos a una malla tipo stars (para usar en autoKrige)
dem1000_grid <- st_as_stars(dem1000)
rm(dem1000)

#------------------------
# Preparamos el paralelo
#------------------------
# Numero de cores
n.cores <- 42
# Armamos el cluster (Usar FORK en Linux, PSOCK en otros SO)
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
doParallel::registerDoParallel(cl = my.cluster)

vars <- c('p_mes', 'tm_min', 'tm_max')
years <- unique(clim_data_all$Year)
months <- unique(clim_data_all$Month)

system.time( # Esto es solo para ver cuanto tarda
  # Bucle for en paralelo y anidado por variable, anyo y mes
  # combine = 'c' indica que los resultados de los bucles anidados se guarden en una sola lista
  krig <- foreach(var = seq_along(vars),
                .combine='c', .multicombine=TRUE)  %:% # Esto indica el anidamiento
  foreach(month = seq_along(months),
          .combine='c', .multicombine=TRUE)  %:%
  foreach(year = seq_along(years),
          .combine='c', .multicombine=TRUE,
          .packages = c('tidyverse', # Hay que especificar los paquetes que se usan
                        'stars',
                        'terra',
                        'gstat',
                        'automap',
                        'geoR',
                        'sf')) %dopar% {
            
            # Formula para el krigging       
            form <- as.formula(paste(vars[[var]], 'altitud', sep = ' ~ '))
            
            # Filtramos los datos para que correspondan al anyo y mes de interes
            test_year <- clim_data_all %>% 
              filter(Year == years[year])
            test_month <- test_year %>% 
              filter(Month == months[month]) %>%
              drop_na() %>%
              distinct(geometry, 
                       .keep_all = TRUE) %>%
              filter((ta_min-ta_max) != 0 | # Filtramos posibles datos raros
                       ta_max >= 50 |
                       ta_max <= -30 |
                       ta_min >= 40 |
                       ta_min <= -35)
            
            # Autokrigging con el paquete automap
            krg<-autoKrige(formula = form, # la formula del modelo
                           input_data = test_month, # datos de entrada
                           new_data = dem1000_grid, # grid para predecir
                           verbose = F) # que nos de la turra o no
                           # tipo de modelo (Se puede no especificar y la funcion lo elige de forma automatica)
            
            # Validacion cruzada con leave one out
            krg_cv<-autoKrige.cv(formula = form, 
                                 input_data = test_month)
            
            krg_extract <- st_extract(x = krg$krige_output, 
                                      at = mus_coords)
            rm(krg)
            krg_extract$Year <- years[[year]]
            krg_extract$Month <- months[[month]]
            krg_extract$Var <- vars[[var]]
            krg_extract$MU_DYN <- MAIN_DETAILED$MU_DYN
            krg_extract$MU_DET <- MAIN_DETAILED$MU_DET
            gc()
            # Guardamos el krigging y la validacion cruzada
            list(krg_extract, krg_cv)
            }
)
# Desarmamos el cluster para no consumir recursos
parallel::stopCluster(cl = my.cluster)
clim_mu <- krig[seq(from = 1, to = length(krig), by = 2)]
save(clim_mu, krig, file = 'Krigging.RData')
