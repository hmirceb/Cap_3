library(tidyverse)
library(readxl)
library(sf)

rm(list=ls())
graphics.off()
setwd("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")
setwd("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")

# Cargamos datos MONITO
data <- read_excel("Datos/Poblaciones/MONITO.xlsx",
                   sheet = 1,
                   col_types = 
                     c('text', 'text', 'text', 'text', 'text',
                       'numeric', 'numeric', 'numeric', 'text', 'text'))

# Quitamos espacios en blanco y seguimientos malos
MONITO_DATA <- data %>%
  mutate(across(MU_DYN:Year, trimws)) %>%
  mutate(across(Comentarios:Extension, trimws)) %>%
  filter(!is.na(TAXON)) %>%
  mutate(TAXON = str_remove_all(TAXON, '\\.')) %>% # Quitamos los puntos porque dan problemas luego
  group_by(MU_DYN, TAXON, Indicator) %>%
  mutate(Plot_unif = as.numeric(as.factor(Plot))) %>% # Unificamos los nombres de los plots entre MUs
  ungroup() %>%
  group_by(MU_DYN, Plot_unif, TAXON, Indicator) %>% # Quitamos los seguimientos con menos de dos anyos
  ungroup()  %>%
  group_by(l = cumsum(is.na(N) & is.na(N) != dplyr::lag(is.na(N))), # Esto cuenta los NAs
           MU_DYN, Plot_unif, TAXON, Indicator) %>% # previos, sale de stackoverflow
  mutate(lagsrtsz = replace_na(ifelse(lag(is.na(N)) & 
                                        !is.na(N), sum(is.na(N)), 0), 0)) %>%
  ungroup() %>% 
  dplyr::select(-l) %>% # quita la variable dummy que se usa para contar NAs
  mutate(lagsrtsz = lagsrtsz+1) %>% # suma uno para que los anyos buenos sean 1
  group_by(MU_DYN, Plot_unif, TAXON, Indicator) %>%
  mutate(lagsrtsz = ifelse(Year == min(Year[!is.na(N)]) & !is.na(N), -1, lagsrtsz)) %>% # pone el primer anyo de cada seguimienot como -1
  ungroup() %>% 
  mutate(lagsrtsz = ifelse(is.na(N), 0, lagsrtsz)) # Pone 0 en los NAs


# Cargamos y formateamos MONITO_TRENDS
MONITO_TRENDS_2022 <- read_excel("~/Dropbox/DATA__LAB/__MONITO/MONITO_ADOPTA/6_RESULTADOS/MONITO_TRENDS_2022.xlsx", 
                                 sheet = "data_2022", skip = 1)
MONITO_TRENDS_2022 <- read_excel("C:/Users/18172844S/Dropbox/DATA__LAB/__MONITO/MONITO_ADOPTA/6_RESULTADOS/MONITO_TRENDS_2022.xlsx", 
                                 sheet = "data_2022", skip = 1)

names(MONITO_TRENDS_2022)[28:40] <- paste('Lambda',
                                          unlist(lapply(strsplit(names(MONITO_TRENDS_2022)[28:40], '\\...'), function(x) x[1])),
                                          sep = ' ')
names(MONITO_TRENDS_2022)[45:58] <- paste('Lambda R',
                                          unlist(lapply(strsplit(names(MONITO_TRENDS_2022)[45:58], '\\...'), function(x) x[1])),
                                          sep = ' ')
names(MONITO_TRENDS_2022)[62:71] <- paste('Incidence',
                                          unlist(lapply(strsplit(names(MONITO_TRENDS_2022)[62:71], '\\...'), function(x) x[1])),
                                          sep = ' ')
names(MONITO_TRENDS_2022)[75:83] <- paste('Plant cover',
                                          unlist(lapply(strsplit(names(MONITO_TRENDS_2022)[75:83] , '\\...'), function(x) x[1])),
                                          sep = ' ')

MONITO_TRENDS_2022 <- MONITO_TRENDS_2022 %>%
  dplyr::select(-c(41:44, 58:61, 72:74, 84:86)) %>%
  dplyr::select(-c("Indicator":"annual", 'Indicator_obs')) %>%
  mutate(across("Lambda 10-11":"Plant cover 22-23", as.numeric)) %>%
  pivot_longer(names_to = "Dummy",
               cols = "Lambda 10-11":"Plant cover 22-23",
               values_to = 'Change') %>%
  filter(!is.na(Change)) %>%
  mutate(Indicator = ifelse(startsWith(Dummy, 'Lambda R') | startsWith(Dummy, 'Plant cover'), 
                            unlist(lapply(strsplit(Dummy, ' '), function(x) paste(x[1], x[2]))),
                            unlist(lapply(strsplit(Dummy, ' '), function(x) x[1]))),
         Years = ifelse(startsWith(Dummy, 'Lambda R') | startsWith(Dummy, 'Plant cover'), 
                        unlist(lapply(strsplit(Dummy, ' '), function(x) x[3])),
                        unlist(lapply(strsplit(Dummy, ' '), function(x) x[2])))) %>%
  dplyr::select(-Dummy) %>%
  rename('N3' = 'N3+/N4') %>%
  mutate(TAXON = str_remove_all(TAXON, '\\.'))

#------------------------------------------------------------
# Calculamos los valores medios de T y la precipitacion total
# para el anyo hidrologico 
#-------------------------------------------------------------
load("Datos/Clima/Krigging.RData")
clim_mu <- lapply(clim_mu, st_drop_geometry)
clim_mu <- do.call('rbind', clim_mu)
clim_mu <- clim_mu %>%
  group_by(Year, Month, Var, MU_DYN) %>%
  summarise(var1.pred = mean(var1.pred)) %>%
  ungroup() %>%
  pivot_wider(names_from = Var, 
              values_from = var1.pred) %>%
  mutate(Month = ifelse(nchar(Month) == 1, paste0('0', Month), Month))

start_month <- '09'
mus <- unique(clim_mu$MU_DYN)
results <- c()
for(mu in seq_along(mus)){
  dat<-clim_mu %>%
    filter(MU_DYN == mus[mu]) %>%
    arrange(Year, Month)
  
  start_row <- which(dat$Month == start_month, arr.ind = T)[1]
  NYears <- length(unique(dat$Year))
  result_period <- c()
  for(i in seq_len(NYears-1)){
    first_month <- 12*(i-1)+start_row
    last_month <- first_month+11
    t <- dat[first_month:last_month,]
    result_period[[i]] <- data.frame(MU_DYN = unique(t$MU_DYN)[1],
                                     start = unique(t$Year)[1], 
                                     end = unique(t$Year)[2],
                                     pcp_m = sum(t$p_mes),
                                     tmax_m = mean(t$tm_max, na.rm = T),
                                     tmin_m = mean(t$tm_min, na.rm = T))
  }
  results[[mu]] <- do.call('rbind', result_period)
}
temps_hydro <- do.call('rbind', results)

#---------------------------------------------------------------------------
# Calculamos la temperatura media para el periodo normal (ultimos 30 anyos)
#---------------------------------------------------------------------------

temps_normal <- temps_hydro %>%
  filter(start < 1991) %>%
  group_by(MU_DYN) %>%
  summarise(pcp_anual = mean(pcp_m),
            tmin_anual = mean(tmin_m),
            tmax_anual = mean(tmax_m))

# #------------------------------------------------------------
# # Calculamos los valores medios de T y la precipitacion total
# # para el anyo hidrologico 
# #-------------------------------------------------------------
# load("Datos/Diarios_MU.RData")
# 
# start_month <- '09'
# mus <- unique(coord_abiotic$MU_DYN)
# results <- c()
# for(mu in seq_along(mus)){
#   dat<-coord_abiotic %>%
#     filter(MU_DYN == mus[mu]) %>%
#     group_by(MU_DYN, Year, Month) %>%
#     summarise(pcp_m = sum(pcp),
#               tmax_m = mean(tmedmax, na.rm = T),
#               tmin_m = mean(tmedmin, na.rm = T)) %>%
#     ungroup()
#   
#   start_row <- which(dat$Month == start_month, arr.ind = T)[1]
#   NYears <- length(unique(dat$Year))
#   result_period <- c()
#   for(i in seq_len(NYears-1)){
#     first_month <- 12*(i-1)+start_row
#     last_month <- first_month+11
#     t <- dat[first_month:last_month,]
#     result_period[[i]] <- data.frame(MU_DYN = unique(t$MU_DYN)[1],
#                                      start = unique(t$Year)[1], 
#                                      end = unique(t$Year)[2],
#                                      pcp_m = sum(t$pcp_m),
#                                      tmax_m = mean(t$tmax_m, na.rm = T),
#                                      tmin_m = mean(t$tmin_m, na.rm = T))
#   }
#   results[[mu]] <- do.call('rbind', result_period)
# }
# temps_hydro <- do.call('rbind', results)
# 
# #---------------------------------------------------------------------------
# # Calculamos la temperatura media para el periodo normal (ultimos 30 anyos)
# #---------------------------------------------------------------------------
# 
# temps_normal <- temps_hydro %>%
#   filter(start < 1991) %>%
#   group_by(MU_DYN) %>%
#   summarise(pcp_anual = mean(pcp_m),
#             tmax_anual = mean(tmax_m),
#             tmin_anual = mean(tmin_m))
# 
# #------------------------
# # Juntamos todo en un DF
# #-----------------------

MONITO_TEMPS <- temps_hydro %>%
  left_join((temps_normal %>%
               dplyr::select(MU_DYN, tmax_anual, pcp_anual, tmin_anual)),
            'MU_DYN',
            multiple = 'all') %>%
  distinct()

save(MONITO_DATA, MONITO_TRENDS_2022, MONITO_TEMPS, file = 'Datos/MONITO_2022.RData')
