library(tidyverse)
library(lme4)
library(rjags)
library(runjags)
library(readxl)
rm(list=ls())
graphics.off()
setwd("~/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")
setwd("C:/Users/18172844S/Dropbox/DATA__LAB/Hector_tesis/Cap. 3 - Population trends/")

load("Datos/MONITO_2022.RData")
doubles <- read_excel("Datos/Poblaciones/MONITO.xlsx",
                      sheet = 2)

MONITO_DATA_FULL<-MONITO_DATA %>%
  mutate(N = N+1) %>%
  full_join((MONITO_TRENDS_2022 %>%
                dplyr::select(c(MU_DYN, TAXON, Indicator, N3, Analysis))),
             by = c('MU_DYN', 'TAXON', 'Indicator'),
            multiple = 'all')  %>%
  distinct() %>%
  filter(Analysis == 'good') %>%
  mutate(N = ifelse(is.na(Max), N, 10*N)) %>%
  left_join((MONITO_TEMPS %>% 
               rename(Year = end)), 
            by = c('MU_DYN', 'Year'),
            multiple = 'all') %>%
  filter(!is.na(pcp_m)) %>% # Esto quita los anyos 2021 y 2022
  group_by(MU_DYN, TAXON, Indicator) %>%
  filter(length(unique(Year)) > 4) %>%
  ungroup() %>%
  arrange(MU_DYN, TAXON, Indicator, Plot, Year)

doubles<- doubles %>%
  filter(is.na(Comentarios)) %>%
  filter(!is.na(N2))

monito_mus<-unique(MONITO_DATA_FULL$MU_DYN)
monito_mus <- 'andpyr_ses1'
mod.list.mu <- c()

# Parametros para el modelo
# final number of iterations n.chains*(sample-burnin)/thin
nc <- 3
nb <- 500 # burnin
nt <- 500 # thining
ns <- 1000 # samples
nad <- 1000 # adaptation

for (mu in seq_along(monito_mus)) {
  temp_mu<-MONITO_DATA_FULL %>% 
    filter(MU_DYN == monito_mus[mu])
  
  monito_sps <- unique(temp_mu$TAXON)
  mod.list.sps <- c()
  for (sp in seq_along(monito_sps)) {
    temp_sps <- temp_mu %>%
      filter(TAXON == monito_sps[sp])
    
    monito_ind <- unique(temp_sps$Indicator)
    mod.list.ind <- c()
    for (ind in seq_along(monito_ind)) {
      temp_ind <- temp_sps %>%
        filter(Indicator == monito_ind[ind]) %>%
        mutate(MU_DYN = as.factor(MU_DYN),
               Plot_unif = as.factor(Plot_unif),
               Year = as.factor(Year))
      
      logN=log(as.integer(round(temp_ind$N, 0)))
      allcases = length(logN)
      lags=temp_ind$lagsrtsz 
      goodendrows=which(lags>0) 
      Ngoodendrows=length(goodendrows)
      NPlots<-length(unique(temp_ind$Plot_unif))
      NYears <- length(unique(temp_ind$Year))
      pcp <- unique(scale(temp_ind$pcp_m-temp_ind$pcp_anual)[,1])
      pcp_all <- scale(temp_ind$pcp_m-temp_ind$pcp_anual)[,1]
      tmax <- unique(scale(temp_ind$tmax_m-temp_ind$tmax_anual)[,1])
      
      #------------------------------
      # Doubles censuses
      # ----------------------------
      dups_ind <- doubles %>%
        filter(MU_DYN == monito_mus[mu] & 
                 TAXON == monito_sps[sp] & 
                 Indicator == monito_ind[ind])
      dups_sps <- doubles %>%
        filter(TAXON == monito_sps[sp] & 
                 Indicator == monito_ind[ind])
      
      if(dim(dups_ind)[1] != 0){
        dups <- cbind(dups_ind$N, dups_ind$N2)
        Ndups <- dim(dups)[1]
        dupvars=NULL
        for (i in seq_len(Ndups)) {
          dupvars = c(dupvars, var(dups[i,]))
        }
      }
      
      if(dim(dups_ind)[1] == 0 & dim(dups_sps)[1] != 0){
        temp_dups <- dups_sps[dups_sps$Indicator == unique(temp_ind$Indicator)[1],]
        dups <- cbind(temp_dups$N, temp_dups$N2)
        Ndups <- dim(dups)[1]
        dupvars=NULL
        for (i in seq_len(Ndups)) {
          dupvars = c(dupvars, var(dups[i,]))
        }
      }
       if(dim(dups_ind)[1] == 0 & dim(dups_sps)[1] == 0){
         temp <- doubles %>%
           filter(Indicator == monito_ind[ind]) %>%
             select(c(N, N2))
         t_sd <- apply(temp, 1, sd)
         t_mean <- apply(temp, 1, mean)
         t_mod <- lm(t_sd~t_mean)
         sd_dup <- predict(t_mod, 
                           newdata = data.frame(
                             t_mean = mean(temp_ind$N, na.rm = T)))
         
         dupvars <- (sd_dup)^2
         Ndups<-length(dupvars) 
       }
      # NoDups <- (dim(dups_ind)[1] == 0 & dim(dups_sps)[1] == 0)
      
      #-----------------------------------
      # Valores iniciales
      #----------------------------------
      inits_func <- function(chain){
        gen_list <- function(chain = chain){
          list( # Esto hay que elegirlo para todos los priors
            mn.log.lam=rnorm(1, 0, 0.001),
            Mesp.precision = rgamma(1, 0.1, 0.1),
            loglam.esp.prec = rgamma(1, 0.1, 0.1),
            Plot_precision = rgamma(1, 0.1, 0.1),
            Year_precision = rgamma(1, 0.1, 0.1),
            .RNG.name = switch(chain,
                               "1" = "base::Wichmann-Hill",
                               "2" = "base::Marsaglia-Multicarry",
                               "3" = "base::Super-Duper",
                               "4" = "base::Mersenne-Twister",
                               "5" = "base::Wichmann-Hill",
                               "6" = "base::Marsaglia-Multicarry",
                               "7" = "base::Super-Duper",
                               "8" = "base::Mersenne-Twister"),
            .RNG.seed = sample(1:1e+06, 1)
          )
        }
        return(switch(chain,           
                      "1" = gen_list(chain),
                      "2" = gen_list(chain),
                      "3" = gen_list(chain),
                      "4" = gen_list(chain),
                      "5" = gen_list(chain),
                      "6" = gen_list(chain),
                      "7" = gen_list(chain),
                      "8" = gen_list(chain)
        )
        )
      }
      
      
      #------------------------------
      # Modelos
      #-------------------------------
      inits_func <- NULL
      # Si solo hay un plot, quitamos el efecto de los plots porque da problemas
      if(length(unique(temp_ind$Plot)) == 1){
        if(is.na(unique(temp_ind$N3))){# Si NO viene de un N3+ o N4, usamos el error de muestreo
          runjags.options(force.summary = T)
          mod.list.ind[[ind]] <- run.jags('Scripts/Modelos_base/Norm_NoPlot.r', 
                                          n.chains=nc,
                                          data=temp_ind,
                                          burnin=nb, 
                                          thin=nt, 
                                          sample=ns, 
                                          adapt=nad,
                                          method='parallel',
                                          inits = inits_func)}
        
        if(!is.na(unique(temp_ind$N3))){# Si viene de un N3+ o N4, NO usamos el error de muestreo
          runjags.options(force.summary = T)
          mod.list.ind[[ind]] <- run.jags('Scripts/Modelos_base/Norm_NoPlot_NoObsErr.r', 
                                          n.chains=nc,
                                          data=temp_ind,
                                          burnin=nb, 
                                          thin=nt, 
                                          sample=ns, 
                                          adapt=nad,
                                          method='parallel',
                                          inits = inits_func)}
      }
      # Si hay varios plots mantenemos el efecto aleatorio
      if(length(unique(temp_ind$Plot)) != 1){
        if(is.na(unique(temp_ind$N3))){# Si NO viene de un N3+ o N4, usamos el error de muestreo 
          runjags.options(force.summary = T)
          mod.list.ind[[ind]] <- run.jags('Scripts/Modelos_base/Norm.r', 
                                          n.chains=nc,
                                          data=temp_ind,
                                          burnin=nb, 
                                          thin=nt, 
                                          sample=ns, 
                                          adapt=nad,
                                          method='parallel',
                                          inits = inits_func)}
        
        if(!is.na(unique(temp_ind$N3))){# Si viene de un N3+ o N4, no usamos el error de muestreo
          runjags.options(force.summary = T)
          mod.list.ind[[ind]] <- run.jags('Scripts/Modelos_base/Norm_NoObsErr.r', 
                                          n.chains=nc,
                                          data=temp_ind,
                                          burnin=nb, 
                                          thin=nt, 
                                          sample=ns, 
                                          adapt=nad,
                                          method='parallel',
                                          inits = inits_func)}
      }
      
      names(mod.list.ind)[[ind]] <- monito_ind[ind]
    }
    mod.list.sps[[sp]] <- mod.list.ind
    names(mod.list.sps)[[sp]] <- monito_sps[sp]
  }
  mod.list.mu[[mu]] <- mod.list.sps
  names(mod.list.mu)[[mu]] <- monito_mus[mu]
}
jag.mod.list <- unlist(unlist(mod.list.mu, recursive = F), recursive = F)

a<-lapply(jag.mod.list, function(x) as.data.frame(summary(x)))
b<-do.call('rbind', a)
b$MU_DYN <- unlist(lapply(strsplit(rownames(b), '\\.'), function(x) x[1]))
b$TAXON <- unlist(lapply(strsplit(rownames(b), '\\.'), function(x) x[2]))
b$Indicator <- unlist(lapply(strsplit(rownames(b), '\\.'), function(x) x[3]))
b$Var <- unlist(lapply(strsplit(rownames(b), '\\.'), function(x) paste(x[4], x[5], x[6], sep = '.')))
b$Var <- str_remove(b$Var, '.NA'); b$Var <- str_remove(b$Var, '.NA')
rownames(b) <- NULL
b$trans <- str_remove_all(unlist(lapply(strsplit(b$Var, '\\['), function(x) x[2])), '\\]')
b$trans <- ifelse(is.na(b$trans), 0, b$trans)

b %>%
  group_by(MU_DYN, TAXON, Indicator) %>%
  filter(trans != max(as.numeric(trans), na.rm = T)) %>%
  ungroup() %>%
  filter(startsWith(Var, 'annuallamest'))  %>%
  ggplot(aes(x = as.numeric(trans), color = Indicator, group = interaction(MU_DYN, TAXON)))+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_hline(data = b %>%
               filter(startsWith(Var, 'mn.log.lam')),
             aes(yintercept = Median))+
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95),
                position = position_dodge(1),
                width = 0.2)+
  geom_point(aes(y = Median, shape = Indicator),
             position = position_dodge(1), size = 2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  facet_wrap(MU_DYN~TAXON+Indicator)

b %>%
  filter(startsWith(Var, 'annual')) %>%
  ggplot(aes(x = Median))+
  theme_classic()+
  geom_vline(aes(xintercept = quantile(Median, 0.025)),
             linetype = 'dashed')+
  geom_vline(aes(xintercept = quantile(Median, 0.975)), 
             linetype = 'dashed')+
  geom_vline(aes(xintercept = mean(Median)), 
             color = 'red')+
  geom_density()

b %>%
  filter(Var == 'mn.log.lam') %>%
  ggplot(aes(x = Median))+
  theme_classic()+
  geom_vline(aes(xintercept = quantile(Median, 0.025)),
             linetype = 'dashed')+
  geom_vline(aes(xintercept = quantile(Median, 0.975)), 
             linetype = 'dashed')+
  geom_vline(aes(xintercept = mean(Median)), 
             color = 'red')+
  geom_density()

b %>%
  filter(Var == 'pcp_coefficient') %>%
  ggplot(aes(x = Median))+
  theme_classic()+
  geom_vline(aes(xintercept = quantile(Median, 0.025)),
             linetype = 'dashed')+
  geom_vline(aes(xintercept = quantile(Median, 0.975)), 
             linetype = 'dashed')+
  geom_vline(aes(xintercept = median(Median)), 
             color = 'red')+
  geom_density()
