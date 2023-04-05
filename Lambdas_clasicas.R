library(tidyverse)

rm(list=ls())
graphics.off()
setwd("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")
setwd("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")
load('Datos/MONITO_2022.RData')

# Primero calculamos el promedio entre el primer y el doble conteo
MONITO_DATA$N_mean <- ifelse(is.na(MONITO_DATA$N2), MONITO_DATA$N, (MONITO_DATA$N+MONITO_DATA$N2)/2)+1
#MONITO_DATA$N_mean <- ifelse(MONITO_DATA$Indicator %in% c('Plant cover', 'Incidence'),
#100*MONITO_DATA$N_mean/MONITO_DATA$Max, MONITO_DATA$N_mean)

# Nos quedamos solo con los seguimientos buenos
MONITO_DATA <- MONITO_DATA %>%
  full_join((MONITO_TRENDS_2022 %>%
               dplyr::select(c(Analysis, MU_DYN, TAXON)) %>%
               distinct()), by = c('MU_DYN', 'TAXON')) 

#### NOTA IMPORTANTE ####
# La lambda o tasa de cambio calculada a partir de todos los seguimientos 
# juntos (por ejemplo sumando el numero de individuos de todos los plots) 
# es igual a la tasa de cambio del promedio de individuos en cada plot, no a la tasa de cambio 
# promedio de los plots (lo mismo se aplica a cberturas y presencias)


########################
####  Tasas por MU  ####
########################

MONITO_DATA_test <- MONITO_DATA %>% 
  filter(lagsrtsz != 0) %>%
  group_by(MU_DYN, TAXON, Indicator, Plot) %>%
  mutate(trans = n()) %>%
  ungroup() %>%
  filter(trans > 1) %>%
  filter(Analysis == 'good')

mus <- unique(MONITO_DATA_test$MU_DYN)
mu_df_list <- c()
for (mu in 1:length(mus)) {
  t_mu <- MONITO_DATA_test %>%
    filter(MU_DYN == mus[mu])
  
  species_df_list <- c()
  sps <- unique(t_mu$TAXON)
  for (sp in 1:length(sps)) {
    t_sps <- t_mu %>%
      filter(TAXON == sps[sp])
    
    types <- unique(t_sps$Indicator)
    ty <- c()
    types_df_list <- c()
    for (type in 1:length(types)) {
      t_type <- t_sps %>%
        filter(Indicator == types[type])
      
      lam_plot <- c()
      years <- sort(unique(t_type$Year))
      trans <- c()
      for (yy in 2:length(years)) {
        t_type_t <- t_type %>% filter(Year == years[yy] & !is.na(t_type$N_mean))
        t_type_t0 <- t_type[t_type$Year == years[yy-1] & !is.na(t_type$N_mean),]
        common_plots <- intersect(t_type_t$Plot, t_type_t0$Plot)
        t_type_t <- t_type_t %>% filter(Plot %in% common_plots) 
        t_type_t0 <- t_type_t0 %>% filter(Plot %in% common_plots)
        
        Nt <- sum(t_type_t$N_mean, na.rm = T)
          N0 <- sum(t_type_t0$N_mean, na.rm = T)
          pot <- (1/unique(t_type_t$lagsrtsz))
          lam_plot[[yy-1]] <- (Nt/N0)^pot
        
        trans[[yy-1]] <- paste(str_sub(unique(t_type_t0$Year), start = 3), str_sub(unique(t_type_t$Year), start = 3), sep = '-')
      }
      t <- data.frame(Years = unlist(trans), 
                      Lambda = unlist(lam_plot))
      t$Indicator <- rep(types[type], times = dim(t)[1])
      types_df_list[[type]]<- t
    }  
    t <- do.call('rbind', types_df_list)
    t$TAXON <- rep(sps[sp], times = dim(t)[1])
    species_df_list[[sp]] <- t
  }
  t <- do.call('rbind', species_df_list)
  t$MU_DYN <- rep(mus[mu], times = dim(t)[1])
  mu_df_list[[mu]] <- t
}

resultados_mu <- do.call('rbind', mu_df_list)


#### Calculamos la lambda media y los intervalos de confianza
# segun Morris y Doak
resultados_mu$t <- as.numeric(unlist(lapply(strsplit(resultados_mu$Years, '-'), function(x) x[1])))
resultados_mu$t2 <- as.numeric(unlist(lapply(strsplit(resultados_mu$Years, '-'), function(x) x[2])))
resultados_mu$x <- sqrt(resultados_mu$t2-resultados_mu$t)
resultados_mu$x <- ifelse(is.nan(resultados_mu$x), 1, resultados_mu$x)

resultados_mu$t <- as.numeric(unlist(lapply(strsplit(resultados_mu$Years, '-'), function(x) x[1])))
resultados_mu$t2 <- as.numeric(unlist(lapply(strsplit(resultados_mu$Years, '-'), function(x) x[2])))
resultados_mu$x <- sqrt(resultados_mu$t2-resultados_mu$t)
resultados_mu$x <- ifelse(is.nan(resultados_mu$x), 1, resultados_mu$x)
resultados_mu$y <- log(resultados_mu$Lambda)/resultados_mu$x

mus <- unique(resultados_mu$MU_DYN)
mu_df_list <- c()
for (mu in 1:length(mus)) {
  t_mu <- resultados_mu %>%
    filter(MU_DYN == mus[mu])
  sps <- unique(t_mu$TAXON)
  
  species_df_list <- c()
  for (sp in 1:length(sps)) {
    t_sps <- t_mu %>%
      filter(TAXON == sps[sp])
    indicators <- unique(t_sps$Indicator)
    mod <- c()
    for (ind in 1:length(indicators)) {
      t_ind <- t_sps %>%
        filter(Indicator == indicators[ind])
      t_m <- lm(y~0+x, data = t_ind)
      mod[[ind]] <- data.frame(coef = t_m$coefficients[1],
                               var = vcov(t_m)[1],
                               low95 = confint(t_m)[1], 
                               up95 = confint(t_m)[2],
                               Indicator = indicators[ind])
    }
    t <- do.call('rbind', mod)
    t$TAXON <- rep(sps[sp], times = dim(t)[1])
    species_df_list[[sp]] <- t
  }
  t <- do.call('rbind', species_df_list)
  t$MU_DYN <- rep(mus[mu], times = dim(t)[1])
  mu_df_list[[mu]] <- t
}
resultados_confi <- do.call('rbind', mu_df_list)

save(resultados_mu, resultados_confi, file = 'Lambdas_clasicas.RData')

resultados_confi %>%
  mutate(cosa = paste(MU_DYN, TAXON)) %>% 
  ggplot()+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', alpha = 0.5)+
  geom_errorbar(aes(x = cosa, ymax = up95, ymin = low95, color = Indicator))+
  geom_point(aes(x = cosa, y = coef, color = Indicator))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
