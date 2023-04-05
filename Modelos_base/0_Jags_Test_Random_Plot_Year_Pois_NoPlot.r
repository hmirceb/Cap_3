
# this version of the model uses poisson regression, so no estimate of obs variance,since variance is equal to the mean

model{
  gammarate <- 0.5*Mesp.precision
  for (ii in 1:Ndups) { 
    dupvars[ii] ~ dgamma(1/2,gammarate)}  # variance tienen que ser positivas, el 1/2 es el parametro de la forma de la distribucion es el "correct" sampling distirbution de la varianza y gamma esta entre 0 y algo positivo
  
  # this loop says that all annual mean log lambdas are equal to a mean estimate: this can be changed to make lambda dependent on covariates. also, it draws random meaurement error amounts for each year: this has to be done this way as these are used in both the ending and starting calcs for transitions. 
  for (ii in 1:allcases) {
    mn.loglams[ii] <- mn.log.lam + Year_randomeffect[Year[ii]] # esto es como decir que cada valor en un momento dado es igual a la media
    Mepsilon[ii] ~ dnorm(0,Mesp.precision)
  }
  
  # this loop says that the best mean estimate for ending log size for a transition is the the result of adding up the log lambdas and starting sizes + measurement errors. Then, you use these along with the estimated variance in the estimate that comes from the process error in loglambda to compare the answer to the prediction. 
  for (ii in 1:Ngoodendrows) {
    
    poissonlams[goodendrows[ii]] <- exp(sum(mn.loglams[(goodendrows[ii]-lags[goodendrows[ii]]):(goodendrows[ii]-1)])+ log(logN[goodendrows[ii]-lags[goodendrows[ii]]])+
                                         Mepsilon[goodendrows[ii]-lags[goodendrows[ii]]] - Mepsilon[goodendrows[ii]])
    
    logN[goodendrows[ii]] ~ dpois(poissonlams[goodendrows[ii]])
    
    
  }
  
  #This makes the annual lambda estimate, adding yr effect to mean lam estimate
  for(Year_iterator in 1:NYears) {
    
    annuallamest[Year_iterator] <- mn.log.lam + Year_randomeffect[Year_iterator]
  }
  
  
  #Prior distrs:
  Mesp.precision ~ dgamma(0.001, 0.001) # we could use a lognormal, uniform or gamma, se necesita algo que este entre 0 y algo positivo, da un poco mas igual
  obs.error.sd <- sqrt(1/Mesp.precision)

  mn.log.lam ~ dnorm(0, 0.001)
  Year_precision ~ dgamma(0.001, 0.001)
  year.sd <- sqrt(1/Year_precision)
  
  # for(Year_iterator in 1:allcases){ # Random effect por plot
  #   Year_randomeffect[Year_iterator] ~ dnorm(0, Year_precision)
  # }
  for(Year_iterator in 1:NYears){ # Random effect por plot
    Year_randomeffect[Year_iterator] ~ dnorm(0, Year_precision)
  }
  
  
}# end of model specification 

# run.jags allows you to specify the list of variables to monitor from the fit, and also the list of data to be imported from the home environment. There is also a list for fitted parameters you wish to get output on:

#monitor#  obs.error.sd, mn.log.lam, annuallamest

#data#  Ndups, dupvars, logN, goodendrows, Ngoodendrows, allcases, lags, Year, NYears

