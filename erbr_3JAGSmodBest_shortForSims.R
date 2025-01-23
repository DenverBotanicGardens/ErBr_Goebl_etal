
model{

# In the BUGS/JAGS language we must use an explicit for loop:

#########################################################################
  ## Doing one loop to generate the estimated size for each plant in each year with a survival data point (all these have size, except the survival=0 years) 
  ## This loop also estimates the survival to each year: for unobserved years, this includes the chances of survival each unobserved year in the past 
  ## Cases in this loop are pairs of observation years, that can include years without any observations 
  for(i in 1:Ncases){
  
  ## Do first lag, which must use the starting size (RosNew)
  regression_mean[goodrows[i]-lagvals[i]+1] <- exp(grwth_intercept + grwth_RosCoef*log(RosNew[goodrows[i]-lagvals[i]]) + grwth_TempFallCoef*TempFall[goodrows[i]-lagvals[i]+1] 
                                               + grwth_TempSummerCoef*TempSummer[goodrows[i]-lagvals[i]+1] + grwth_TempWinterCoef*TempWinter[goodrows[i]-lagvals[i]+1]
                                               + grwth_PptFallCoef*PptFall[goodrows[i]-lagvals[i]+1] + grwth_PptSummerCoef*PptSummer[goodrows[i]-lagvals[i]+1]
                                               + grwth_PptWinterCoef*PptWinter[goodrows[i]-lagvals[i]+1] + grwth_Transect_randomeffect[TransectNew.num[goodrows[i]-lagvals[i]]]) 
      
  r.growth[goodrows[i]-lagvals[i]+1] <- exp(grwthvar_intercept + grwthvar_RosCoef*log(RosNew[goodrows[i]-lagvals[i]])) 

  ## Survival prob, based on the last years observed size 
  Surv_mu[goodrows[i]-lagvals[i]+1] <- 1/(1+exp(-(surv_intercept + surv_RosCoef*log(RosNew[goodrows[i]-lagvals[i]]) + surv_PptWinterCoef*PptWinter[goodrows[i]-lagvals[i]+1] 
                                       + surv_TempWinterCoef*TempWinter[goodrows[i]-lagvals[i]+1] + surv_TempFallCoef*TempFall[goodrows[i]-lagvals[i]+1] 
                                       + surv_TempSummerCoef*TempSummer[goodrows[i]-lagvals[i]+1] + surv_Transect_randomeffect[TransectNew.num[goodrows[i]-lagvals[i]]])))

  
  ## This is the loop of remaining lagged years for this row i: note that the : operator is different in jags than R & can't be decreasing, hence the use of negative lag: 
  ## also, if lag=1, then this will be from 0 to -1, & the loop will be skipped 
  ## Survival rates are based on the inferred previous year size, since it was not observed 
  for (j in (goodrows[i]-lagvals[i]+1):(goodrows[i]-1)) { 
    
     regression_mean[j+1] <- exp(grwth_intercept + grwth_RosCoef*log(regression_mean[j]) + grwth_TempFallCoef*TempFall[j+1] 
                             + grwth_TempSummerCoef*TempSummer[j+1] + grwth_TempWinterCoef*TempWinter[j+1]
                             + grwth_PptFallCoef*PptFall[j+1] + grwth_PptSummerCoef*PptSummer[j+1]
                             + grwth_PptWinterCoef*PptWinter[j+1] + grwth_Transect_randomeffect[TransectNew.num[j]])
    
     r.growth[j+1] <- exp(grwthvar_intercept + grwthvar_RosCoef*log(regression_mean[j]))

     Surv_mu[j+1] <- Surv_mu[j]*1/(1+exp(-(surv_intercept + surv_RosCoef*log(regression_mean[j]) + surv_PptWinterCoef*PptWinter[j+1] 
                                + surv_TempWinterCoef*TempWinter[j+1] + surv_TempFallCoef*TempFall[j+1] 
                                + surv_TempSummerCoef*TempSummer[j+1] + surv_Transect_randomeffect[TransectNew.num[j]])))

   } #End lags loop 
} #End of going through cases for growth and survival
#####################################################################
  
 
  
  
  #####################################################################
  ## Then a loop that matches the estimated sizes with the observed sizes, just for rows with observed sizes
  for(i in 1:Ngrowcases){
    ## These lines describe the response distribution and linear model terms
    
    regression_residual[i] <- RosNew[goodgrowrows[i]] - regression_mean[goodgrowrows[i]]
    
    p[goodgrowrows[i]] <- r.growth[goodgrowrows[i]]/(r.growth[goodgrowrows[i]]+regression_mean[goodgrowrows[i]])
    
    RosNew[goodgrowrows[i]] ~ dnegbin(p[goodgrowrows[i]], r.growth[goodgrowrows[i]])T(1,)
    
  } #End of going through cases for sizes
  #####################################################################

  #####################################################################
  ## Now a loop to do survival predictions, for the ending years with observations
  for(i in 1:Ncases){
    ## Fit the survival function:
    Survs[goodrows[i]] ~ dbern(Surv_mu[goodrows[i]])
  } #End loop to do survival predictions
#####################################################################
  

  
  
########################################
########################################    
## These lines give the prior distributions for the parameters to be estimated
	# Note: the prior for the dispersion parameter k is quite important for convergence

## GROWTH  
grwth_intercept ~ dunif(-5,5) 
grwth_RosCoef ~ dunif(-3,3) 
grwth_PptFallCoef ~ dunif(-2,2)
grwth_PptSummerCoef ~ dunif(-2,2)
grwth_PptWinterCoef ~ dunif(-2,2)
grwth_TempFallCoef ~ dunif(-2,2)
grwth_TempSummerCoef ~ dunif(-2,2)
grwth_TempWinterCoef ~ dunif(-2,2)

## Growth transect random effects
for(TransectNew.num_iterator in 1:numtrans){
  grwth_Transect_randomeffect[TransectNew.num_iterator] ~ dnorm(0, grwth_Transect_precision) }
grwth_Transect_precision ~ dunif(1,10) 


## VARIANCE IN GROWTH
grwthvar_intercept ~ dunif(-3,3) 
grwthvar_RosCoef ~ dunif(-3,3) 


## SURVIVAL
surv_intercept ~  dnorm(0, 10^-6)
surv_RosCoef ~  dnorm(0, 10^-6)
surv_PptWinterCoef ~ dnorm(0, 10^-6)
surv_TempFallCoef ~ dnorm(0, 10^-6)
surv_TempSummerCoef ~ dnorm(0, 10^-6) 
surv_TempWinterCoef ~ dnorm(0, 10^-6)

## Survival transect random effects
for(TransectNew.num_iterator in 1:12){
  surv_Transect_randomeffect[TransectNew.num_iterator] ~ dnorm(0, surv_Transect_precision) }
surv_Transect_precision ~ dunif(0,1) 





#resid.sum.sq <- sum(regression_residual^2)
} # end of model specification 


# These lines are hooks to be read by runjags (they are ignored by JAGS):
#monitor# deviance,grwth_Transect_randomeffect,surv_Transect_randomeffect,grwth_intercept,grwth_RosCoef,grwth_TempFallCoef,grwth_TempSummerCoef,grwth_TempWinterCoef,grwth_PptFallCoef,grwth_PptSummerCoef,grwth_PptWinterCoef,grwthvar_intercept,grwthvar_RosCoef,surv_intercept,surv_RosCoef,surv_PptWinterCoef,surv_TempFallCoef,surv_TempSummerCoef,surv_TempWinterCoef,surv_Transect_precision,grwth_Transect_precision
#modules# glm on
#response# RosNew
#residual# regression_residual
#fitted# regression_fitted
#data#  Ncases,Ngrowcases,goodrows,goodgrowrows,lagvals,TransectNew.num,RosNew,InflNew,InflYesNo,Survs,rows.w.sz,rows.wo.sz,Ndirectszcases,Nindirectszcases,numtrans,rows.w.inflors,Nrows.w.inflors,newPltlines,yrtranscombo,newplt.yrtranscombo,newplts,PptFall,PptWinter,PptSummer,TempWinter,TempFall,TempSummer,rows.wo.sz.alive



######################################################################################################
#### Initial values 
######################################################################################################
######################################################################################################

inits{
  "grwth_intercept" <- 2
  "grwth_RosCoef" <- 1
  "grwth_PptFallCoef" <- 0.01
  "grwth_PptSummerCoef" <- 0.01
  "grwth_PptWinterCoef" <- 0.01
  "grwth_TempFallCoef" <- 0.01
  "grwth_TempSummerCoef" <- 0.01
  "grwth_TempWinterCoef" <- 0.01
  "grwth_Transect_precision" <- 2
  
  "grwthvar_intercept" <- 0.01
  "grwthvar_RosCoef" <- 0.01
  
  "surv_intercept" <- -1
  "surv_RosCoef" <- 0.01
  "surv_PptWinterCoef" <- 0.01
  "surv_TempWinterCoef" <- 0.01
  "surv_TempFallCoef" <- 0.01
  "surv_TempSummerCoef" <- 0.01
  "surv_Transect_precision" <- 0.001


}
inits{
  "grwth_intercept" <- 1.5
  "grwth_RosCoef" <- 1 
  "grwth_PptFallCoef" <- 0.001
  "grwth_PptSummerCoef" <- 0.001
  "grwth_PptWinterCoef" <- 0.001
  "grwth_TempFallCoef" <- 0.001
  "grwth_TempSummerCoef" <- 0.001
  "grwth_TempWinterCoef" <- 0.001
  "grwth_Transect_precision" <- 2
  
  "grwthvar_intercept" <- 0.1
  "grwthvar_RosCoef" <- 0.1
  
  "surv_intercept" <- -10
  "surv_RosCoef" <- -0.01
  "surv_PptWinterCoef" <- 0.05 
  "surv_TempWinterCoef" <- 0.5 
  "surv_TempFallCoef" <- 0.5
  "surv_TempSummerCoef" <- 0.5
  "surv_Transect_precision" <- 0.001

 
}
inits{
  "grwth_intercept" <- 1.5
  "grwth_RosCoef" <- 1
  "grwth_PptFallCoef" <- 0.001
  "grwth_PptSummerCoef" <- 0.001
  "grwth_PptWinterCoef" <- 0.001
  "grwth_TempFallCoef" <- 0.001
  "grwth_TempSummerCoef" <- 0.001
  "grwth_TempWinterCoef" <- 0.001
  "grwth_Transect_precision" <- 2
  
  "grwthvar_intercept" <- 0.01
  "grwthvar_RosCoef" <- -0.01
  
  "surv_intercept" <- -20
  "surv_RosCoef" <- 0.01
  "surv_PptWinterCoef" <- 0.05
  "surv_TempWinterCoef" <- 1
  "surv_TempFallCoef" <- 0.5
  "surv_TempSummerCoef" <- 0.5
  "surv_Transect_precision" <- 0.001

 
}