library(tidyverse)
library(readxl)
library(rjags)
library(runjags)
library(popbio)
rm(list = ls())
setwd('/home/hector/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/')

MAIN_GENERAL <- read_excel("~/Dropbox/DATA__LAB/__MONITO/MONITO_ADOPTA/5_DATOS & ANÁLISIS/MONITO_FILEMAKER/FM_EXPORTED/Main_Tables_090123/MAIN_GENERAL.xlsx")
load("Datos/JAGS_Mods.RData")
load("Datos/JAGS_Mods_Clim.RData")
load("Datos/MONITO_2022.RData")
load("Datos/Lambdas_clasicas.RData")

N_final <- MONITO_DATA %>%
  mutate(N = N+1) %>%
  group_by(MU_DYN, TAXON, Indicator) %>%
  filter(length(unique(Year)) > 4) %>%
  ungroup() %>%
  group_by(MU_DYN, TAXON, Indicator, Year) %>%
  summarise(N_f = sum(N, na.rm = T)) %>%
  filter(Year == max(Year)) %>%
  ungroup()

props <- matrix(NA, ncol = 5, nrow = length(jag.mod.list))
tmax <- 30
for (mu in seq_along(jag.mod.list)) {
  mu_dyn <- unlist(lapply(strsplit(names(jag.mod.list)[mu], '\\.'), function(x) x[1]))
  sp <- unlist(lapply(strsplit(names(jag.mod.list)[mu], '\\.'), function(x) x[2]))
  ind <- unlist(lapply(strsplit(names(jag.mod.list)[mu], '\\.'), function(x) x[3]))
  N1 <- N_final %>%
    filter(MU_DYN == mu_dyn &
             TAXON == sp &
             Indicator == ind) %>%
    pull(N_f)
  chains <- do.call('rbind',
                    lapply(jag.mod.list[[mu]]$mcmc, 
                           function(x) x[(dim(x)[1]-999):dim(x)[1], 1:2]))
  t <- apply(chains, 1, function(x) extCDF(mu = x[1], 
                                           sig2 = x[2]^2,
                                           Nc = N1,
                                           Ne = 0.9*N1,
                                           tmax = tmax)[tmax])
  props[mu,1] <- mu_dyn
  props[mu,2] <- sp
  props[mu,3] <- ind
  props[mu,4] <- mean(t, na.rm = T)
  props[mu,5] <- sd(t, na.rm = T)
}

ext_prob <- as.data.frame(props)
names(ext_prob) <- c('MU_DYN', 'TAXON', 'Indicator', 'MeanEXT', 'sdEXT')
ext_prob$Type <- 'Bayes'

#--------------------------------------
# Pseudoextinction con lambdas clasicas
#--------------------------------------
N_final <- MONITO_DATA %>%
  mutate(N = N+1) %>%
  group_by(MU_DYN, TAXON, Indicator) %>%
  filter(length(unique(Year)) > 1) %>%
  ungroup() %>%
  group_by(MU_DYN, TAXON, Indicator, Year) %>%
  summarise(N_f = sum(N, na.rm = T)) %>%
  filter(Year == max(Year)) %>%
  ungroup()

mus <- unique(resultados_confi$MU_DYN)
ppp <- list()
for (mu in seq_along(mus)) {
  t_mu <- resultados_confi %>%
    filter(MU_DYN == mus[mu])
  
  sps <- unique(t_mu$TAXON)
  for (sp in seq_along(sps)) {
    t_sps <- t_mu %>%
      filter(TAXON == sps[sp])
    ind <- unique(t_sps$Indicator)
    
    for (type in seq_along(ind)) {
      t_ind <- t_sps %>%
        filter(Indicator == ind[type])
      N1 <- N_final %>%
        filter(MU_DYN == mus[mu] &
                 TAXON == sps[sp] &
                 Indicator == ind[type]) %>%
        pull(N_f)
      t <- extCDF(mu = t_ind$coef, 
                  sig2 = t_ind$var,
                  Nc = N1,
                  Ne = 0.9*N1,
                  tmax = tmax)[tmax]
      z <- data.frame(MU_DYN = mus[mu], TAXON = sps[sp], Indicator = ind[type], t)
      ppp <- append(ppp,
                    z)
    }
  }
}
ext_prob_classic <- as.data.frame(matrix(unlist(ppp), ncol = 4, byrow = T))
names(ext_prob_classic) <- c('MU_DYN', 'TAXON', 'Indicator', 'MeanEXT')
ext_prob_classic$MeanEXT <- as.numeric(ext_prob_classic$MeanEXT)
ext_prob_classic$Type <- 'Classic'

ext_prob <- ext_prob %>%
  mutate(MeanEXT = as.numeric(MeanEXT)) %>%
  bind_rows(ext_prob_classic)

save(ext_prob, file = 'Datos/Extinction_Probabilites.RData')
