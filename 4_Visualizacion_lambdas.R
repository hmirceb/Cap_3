library(ggridges)
library(tidyverse); theme_set(theme_classic())
library(readxl)
library(rjags)
library(runjags)
rm(list = ls())
setwd('/home/hector/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/')

MAIN_GENERAL <- read_excel("~/Dropbox/DATA__LAB/__MONITO/MONITO_ADOPTA/5_DATOS & ANÁLISIS/MONITO_FILEMAKER/FM_EXPORTED/Main_Tables_090123/MAIN_GENERAL.xlsx")
load("Datos/JAGS_Mods.RData")
load("Datos/JAGS_Mods_Clim.RData")
load("Datos/Extinction_Probabilites.RData")
load("Datos/MONITO_2022.RData")
load("Datos/Lambdas_clasicas.RData")

#----------------------------------------------
# Creamos DF con los modelos sin datos de clima
#----------------------------------------------
mod_df_list<-lapply(jag.mod.list, function(x) as.data.frame(summary(x)))
mod_df<-do.call('rbind', mod_df_list)
mod_df$MU_DYN <- unlist(lapply(strsplit(rownames(mod_df), '\\.'), function(x) x[1]))
mod_df$TAXON <- unlist(lapply(strsplit(rownames(mod_df), '\\.'), function(x) x[2]))
mod_df$Indicator <- unlist(lapply(strsplit(rownames(mod_df), '\\.'), function(x) x[3]))
mod_df$Var <- unlist(lapply(strsplit(rownames(mod_df), '\\.'), function(x) paste(x[4], x[5], x[6], sep = '.')))
mod_df$Var <- str_remove(mod_df$Var, '.NA'); mod_df$Var <- str_remove(mod_df$Var, '.NA')
rownames(mod_df) <- NULL
mod_df$trans <- str_remove_all(unlist(lapply(strsplit(mod_df$Var, '\\['), function(x) x[2])), '\\]')
mod_df$trans <- ifelse(is.na(mod_df$trans), 0, mod_df$trans)
mod_df$Clim <- 'No'

#----------------------
# DF con datos de clima
#----------------------
mod_df_list_clim<-lapply(jag.mod.list.clim, function(x) as.data.frame(summary(x)))
mod_df_clim<-do.call('rbind', mod_df_list_clim)
mod_df_clim$MU_DYN <- unlist(lapply(strsplit(rownames(mod_df_clim), '\\.'), function(x) x[1]))
mod_df_clim$TAXON <- unlist(lapply(strsplit(rownames(mod_df_clim), '\\.'), function(x) x[2]))
mod_df_clim$Indicator <- unlist(lapply(strsplit(rownames(mod_df_clim), '\\.'), function(x) x[3]))
mod_df_clim$Var <- unlist(lapply(strsplit(rownames(mod_df_clim), '\\.'), function(x) paste(x[4], x[5], x[6], sep = '.')))
mod_df_clim$Var <- str_remove(mod_df_clim$Var, '.NA'); mod_df_clim$Var <- str_remove(mod_df_clim$Var, '.NA')
rownames(mod_df_clim) <- NULL
mod_df_clim$trans <- str_remove_all(unlist(lapply(strsplit(mod_df_clim$Var, '\\['), function(x) x[2])), '\\]')
mod_df_clim$trans <- ifelse(is.na(mod_df_clim$trans), 0, mod_df_clim$trans)
mod_df_clim$Clim <- 'Yes'

#--------------------------------------------
# Juntamos todos los resultados en un solo DF
#--------------------------------------------
RESULTS_DF <- MAIN_GENERAL %>%
  dplyr::select(c(MU_DYN, HABITAT, MU_ALT, LON, LAT, EUNIS_DESC)) %>%
  inner_join(bind_rows(mod_df, mod_df_clim), by = 'MU_DYN')  %>%
  inner_join(ext_prob, by = c('MU_DYN', 'TAXON', 'Indicator')) 

#---------
# Graficos
#---------
# Mean log-Lambda and confidence interval
RESULTS_DF %>%
  filter(startsWith(Var, 'annual')) %>%
  ggplot(aes(x = Median))+
  geom_vline(aes(xintercept = quantile(Median, 0.025)),
             linetype = 'dashed')+
  geom_vline(aes(xintercept = quantile(Median, 0.975)), 
             linetype = 'dashed')+
  geom_vline(aes(xintercept = mean(Median)), 
             color = 'red')+
  geom_density()

# log-Lambda por altitud sin tener en cuenta el clima
RESULTS_DF %>%
  mutate(Alt_group = cut_number(MU_ALT, 10)) %>%
  filter(Var == 'mn.log.lam' &
           Type == 'Bayes' &
           Clim == 'No') %>%
  ggplot(aes(x = Median, y = Alt_group))+
  geom_density_ridges_gradient(aes(fill = Alt_group), alpha = 0.7)+
  geom_vline(aes(xintercept = 0), 
             color = 'red')+
  ylab('Altitude m a.s.l.')+
  xlab('log(Lambda)')

# Log-Lambda promedio por habitat sin tener en cuenta el clima
RESULTS_DF %>%
  filter(Var == 'mn.log.lam' &
           Clim == 'No' &
           Type == 'Bayes') %>%
  ggplot(aes(x = Median, y = EUNIS_DESC))+
  geom_density_ridges_gradient(aes(fill = EUNIS_DESC), alpha = 0.7)+
  geom_vline(aes(xintercept = 0), 
             color = 'red')+
  ylab('Altitude m a.s.l.')+
  xlab('log(Lambda)')

# Probabilidad de pseudoextincion a 30 anyos por habitat sin tener en cuenta el clima
RESULTS_DF %>%
  filter(Var == 'mn.log.lam' &
           Clim == 'No' &
           Type == 'Bayes') %>%
  ggplot(aes(x = MeanEXT, y = EUNIS_DESC))+
  geom_density_ridges_gradient(aes(fill = EUNIS_DESC), alpha = 0.2)+
  geom_vline(aes(xintercept = mean(MeanEXT)), 
             color = 'red')+
  ylab('Habitat')+
  xlab('Probability of pseudoextiction')

# Probabilidad de pseudoextincion a 30 anyos por altitud sin tener en cuenta el clima
RESULTS_DF %>%
  mutate(Alt_group = cut_number(MU_ALT, 10)) %>%
  filter(Var == 'mn.log.lam' &
           Clim == 'No' &
           Type == 'Bayes') %>%
  ggplot(aes(x = MeanEXT, y = Alt_group))+
  geom_density_ridges_gradient(aes(fill = Alt_group), alpha = 0.1)+
  geom_vline(aes(xintercept = mean(MeanEXT)), 
             color = 'red')+
  ylab('Altitude m a.s.l.')+
  xlab('Probability of pseudoextiction')

# Boxplot de PE30 por habitats
RESULTS_DF %>%
  mutate(Alt_group = cut_number(MU_ALT, 20)) %>%
  filter(Var == 'mn.log.lam' &
           Clim == 'No' &
           Type == 'Bayes') %>%
  ggplot(aes(x = MeanEXT, y = EUNIS_DESC))+
  geom_boxplot(aes(fill = EUNIS_DESC), alpha = 0.7)+
  geom_vline(aes(xintercept = mean(MeanEXT)), 
             color = 'black', linetype = 'dashed', linewidth = 1)+
  ylab('Altitude m a.s.l.')+
  xlab('Probability of pseudoextiction')

# Relacion entre PE30 y log-Lambda
RESULTS_DF %>%
  filter(Var == 'mn.log.lam' &
           Clim == 'No' &
           Type == 'Bayes') %>%
  ggplot(aes(x = Median, y = MeanEXT))+
  geom_point()+
  geom_smooth()+
  geom_vline(xintercept = 0, linetype = 'dashed')+
  geom_hline(aes(yintercept = mean(MeanEXT)), color = 'red')+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))

# Coefficientes de pcp, tmax y tmin por habitat
RESULTS_DF %>%
  filter(Var %in% c('pcp_coefficient', 'tmax_coefficient', 'tmin_coefficient') &
           Clim == 'Yes' &
           Type == 'Bayes') %>%
  ggplot(aes(y = Median, x = EUNIS_DESC))+
  geom_boxplot(aes(fill = EUNIS_DESC), varwidth = T)+
  geom_hline(yintercept = 0)+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  facet_wrap(~Var)

# Proporcion de especies que suben o bajan
escala <- c('Growing' = '#228833', 'Decreasing' = '#EE6677', 'Stable' = '#DDAA33')

RESULTS_DF %>%
  filter(Var == 'mn.log.lam' &
           Clim == 'No' &
           Type == 'Bayes') %>%
  mutate(l_sig = ifelse(Median > 0, 'Growing',
                        ifelse(Median < 0,
                               'Decreasing', 'Stable'))) %>%
  ggplot(aes(x = EUNIS_DESC))+
  geom_bar(aes(fill = l_sig), color = 'black',
           position = 'fill')+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = escala)

# Proporcion de especies que suben o bajan teniendo en cuanta la incertidumbre
RESULTS_DF %>%
  filter(Var == 'mn.log.lam' &
           Clim == 'No' &
           Type == 'Bayes') %>%
  mutate(l_sig = ifelse(Upper95 > 0 & Lower95 > 0, 'Growing',
                        ifelse(Upper95 < 0 & Lower95 < 0,
                               'Decreasing', 'Stable'))) %>%
  ggplot(aes(x = EUNIS_DESC))+
  geom_bar(aes(fill = l_sig), color = 'black',
           position = 'fill')+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = escala)

# Proporcion de especies que suben o bajan lambda clasica por habitat
resultados_confi %>%
  left_join((RESULTS_DF %>%
               dplyr::select(c(MU_DYN, TAXON, Indicator, EUNIS_DESC, MU_ALT)) %>%
               distinct()),
            by = c('Indicator', 'TAXON', 'MU_DYN')) %>%
  mutate(l_sig = ifelse(coef > 0, 'Growing',
                        ifelse(coef < 0,
                               'Decreasing', 'Stable'))) %>%
  ggplot(aes(x = EUNIS_DESC))+
  geom_bar(aes(fill = l_sig), color = 'black',
           position = 'fill')+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = escala)

# Proporcion de especies que suben o bajan teniendo en cuanta la incertidumbre
# segun labda clasica por habitat
resultados_confi %>%
  left_join((RESULTS_DF %>%
               dplyr::select(c(MU_DYN, TAXON, Indicator, EUNIS_DESC, MU_ALT)) %>%
               distinct()),
            by = c('Indicator', 'TAXON', 'MU_DYN')) %>%
  mutate(l_sig = ifelse(up95 > 0 & low95 > 0, 'Growing',
                        ifelse(up95 < 0 & low95 < 0,
                               'Decreasing', 'Stable'))) %>%
  ggplot(aes(x = EUNIS_DESC))+
  geom_bar(aes(fill = l_sig), color = 'black',
           position = 'fill')+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = escala)

#
resultados_confi %>%
  left_join((RESULTS_DF %>%
               dplyr::select(c(MU_DYN, TAXON, Indicator, EUNIS_DESC, MU_ALT)) %>%
               distinct()),
            by = c('Indicator', 'TAXON', 'MU_DYN')) %>%
  mutate(l_sig = ifelse(coef > 0, 'Growing',
                        ifelse(coef < 0,
                               'Decreasing', 'Stable'))) %>%
  mutate(Alt_group = cut_number(MU_ALT, 10)) %>%
  ggplot(aes(x = Alt_group))+
  geom_bar(aes(fill = l_sig), color = 'black',
           position = 'fill')+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = escala)

#
resultados_confi %>%
  left_join((RESULTS_DF %>%
               dplyr::select(c(MU_DYN, TAXON, Indicator, EUNIS_DESC, MU_ALT)) %>%
               distinct()),
            by = c('Indicator', 'TAXON', 'MU_DYN')) %>%
  mutate(l_sig = ifelse(up95 > 0 & low95 > 0, 'Growing',
                        ifelse(up95 < 0 & low95 < 0,
                               'Decreasing', 'Stable'))) %>%
  mutate(Alt_group = cut_number(MU_ALT, 10)) %>%
  ggplot(aes(x = Alt_group))+
  geom_bar(aes(fill = l_sig), color = 'black',
           position = 'fill')+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = escala)


#
RESULTS_DF %>%
  filter(Var == 'mn.log.lam' &
           Clim == 'No' &
           Type == 'Bayes') %>%
  mutate(l_sig = ifelse(Upper95 > 0 & Lower95 > 0, 'Growing',
                        ifelse(Upper95 < 0 & Lower95 < 0,
                               'Decreasing', 'Stable'))) %>%
  ggplot(aes(x = Median, y = MU_ALT))+
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red')+
  geom_point()+
  geom_errorbar(aes(xmin = Lower95, xmax = Upper95))+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = escala)


#
resultados_confi %>%
  left_join((RESULTS_DF %>%
               dplyr::select(c(MU_DYN, TAXON, Indicator, EUNIS_DESC, MU_ALT)) %>%
               distinct()),
            by = c('Indicator', 'TAXON', 'MU_DYN')) %>%
  mutate(l_sig = ifelse(up95 > 0 & low95 > 0, 'Growing',
                        ifelse(up95 < 0 & low95 < 0,
                               'Decreasing', 'Stable'))) %>%
  ggplot(aes(x = coef, y = MU_ALT))+
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red')+
  geom_point()+
  geom_errorbar(aes(xmin = low95, xmax = up95))+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_fill_manual(values = escala)
