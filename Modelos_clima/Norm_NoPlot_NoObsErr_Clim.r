
# this version of the model uses poisson regression, so no estimate of obs variance,since variance is equal to the mean

model{
  # this loop says that all annual mean log lambdas are equal to a mean estimate: this can be changed to make lambda dependent on covariates. also, it draws random meaurement error amounts for each year: this has to be done this way as these are used in both the ending and starting calcs for transitions. 
  for (ii in 1:allcases) {
    mn.loglams[ii] <- mn.log.lam + 
      tmin_coefficient*tmin_all[ii] + 
      tmax_coefficient*tmax_all[ii] +
      pcp_coefficient*pcp_all[ii] + 
      Year_randomeffect[Year[ii]] # esto es como decir que cada valor en un momento dado es igual a la media
  }
  
  # this loop says that the best mean estimate for ending log size for a transition is the the result of adding up the log lambdas and starting sizes + measurement errors. Then, you use these along with the estimated variance in the estimate that comes from the process error in loglambda to compare the answer to the prediction. 
  for (ii in 1:Ngoodendrows) {
    
    mnestimate[goodendrows[ii]] <- sum(mn.loglams[(goodendrows[ii]-lags[goodendrows[ii]]):(goodendrows[ii]-1)])+logN[goodendrows[ii]-lags[goodendrows[ii]]]
    
    logN[goodendrows[ii]] ~ dnorm(mnestimate[goodendrows[ii]],(loglam.esp.prec/lags[goodendrows[ii]]))
    
    
  }
  
  #This makes the annual lambda estimate, adding yr effect to mean lam estimate
  for(Year_iterator in 1:NYears) {
    
    annuallamest[Year_iterator] <- mn.log.lam + 
      tmin_coefficient*tmin[Year_iterator] + 
      tmax_coefficient*tmax[Year_iterator] + 
      pcp_coefficient*pcp[Year_iterator] + 
      Year_randomeffect[Year_iterator]
  }
  
  
  #Prior distrs:
  mn.log.lam ~ dnorm(0,0.001)
  pcp_coefficient ~ dnorm(0, 0.001)
  tmax_coefficient ~ dnorm(0, 0.001)
  tmin_coefficient ~ dnorm(0, 0.001)
  
  loglam.esp.prec ~ dgamma(0.001, 0.001) 
  process.loglam.sd <- sqrt(1/loglam.esp.prec )
  
  #Year random effect
  Year_precision ~ dgamma(0.001, 0.001) 
  year.sd <- sqrt(1/Year_precision)
  
  for(Year_iterator in 1:NYears){
    Year_randomeffect[Year_iterator] ~ dnorm(0, Year_precision)
  }
  
  
}# end of model specification 

# run.jags allows you to specify the list of variables to monitor from the fit, and also the list of data to be imported from the home environment. There is also a list for fitted parameters you wish to get output on:

#monitor#  mn.log.lam, process.loglam.sd, annuallamest, pcp_coefficient, tmax_coefficient, tmin_coefficient

#data#  logN, goodendrows, Ngoodendrows, allcases, lags, Year, NYears, pcp, pcp_all, tmax, tmax_all, tmin, tmin_all

