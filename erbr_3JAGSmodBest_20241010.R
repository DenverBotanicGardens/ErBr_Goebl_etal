
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

   ## Survival probability, based on the last years observed size 
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
  ## A loop that makes each reproduction yes/no and reproduction amount estimate for years with observed sizes
  for(i in 1:Ndirectszcases){
    
    ## Do the comparison for the reproduction yes/no data in this loop, as all observations either reproduced or not
    InflYesNo[rows.w.sz[i]] ~ dbern(repro_prob[rows.w.sz[i]])
    
    logit(repro_prob[rows.w.sz[i]]) <- reproyesno_intercept + reproyesno_RosCoef*log(RosNew[rows.w.sz[i]]) + reproyesno_TempFallCoef*TempFall[rows.w.sz[i]] 
                                       + reproyesno_TempSummerCoef*TempSummer[rows.w.sz[i]] + reproyesno_TempWinterCoef*TempWinter[rows.w.sz[i]] 
                                       + reproyesno_PptFallCoef*PptFall[rows.w.sz[i]] + reproyesno_PptSummerCoef*PptSummer[rows.w.sz[i]] 
                                       + reproyesno_Transect_randomeffect[TransectNew.num[rows.w.sz[i]]]

    repro_amount[rows.w.sz[i]] <- exp(repro_intercept + repro_RosCoef*log(RosNew[rows.w.sz[i]]) + repro_PptFallCoef*PptFall[rows.w.sz[i]] 
                                  + repro_PptSummerCoef*PptSummer[rows.w.sz[i]] + repro_TempWinterCoef*TempWinter[rows.w.sz[i]] 
                                  + repro_TempFallCoef*TempFall[rows.w.sz[i]] + repro_TempSummerCoef*TempSummer[rows.w.sz[i]] 
                                  + repro_Transect_randomeffect[TransectNew.num[rows.w.sz[i]]])
    
    ## Inflorescences is either the observed or predicted inflorescence number to be used to make a total reproduction estimates to fit predicted inflorescences to new plants
    ## For rows with measured inflorescences, this is just the observed 
    inflors[rows.w.sz[i]] = InflNew[rows.w.sz[i]]
  }
  #####################################################################
  
  #####################################################################
  ## Do the same for reproduction estimates, for years without an observation, so based on inferred size
  for(i in 1:Nindirectszcases){

    logit(repro_prob[rows.wo.sz[i]]) <- reproyesno_intercept + reproyesno_RosCoef*log(regression_mean[rows.wo.sz[i]]) + reproyesno_TempFallCoef*TempFall[rows.wo.sz[i]] 
                                        + reproyesno_TempSummerCoef*TempSummer[rows.wo.sz[i]] + reproyesno_TempWinterCoef*TempWinter[rows.wo.sz[i]] 
                                        + reproyesno_PptFallCoef*PptFall[rows.wo.sz[i]] + reproyesno_PptSummerCoef*PptSummer[rows.wo.sz[i]]
                                        + reproyesno_Transect_randomeffect[TransectNew.num[rows.wo.sz[i]]]

    repro_amount[rows.wo.sz[i]] <- exp(repro_intercept + repro_RosCoef*log(regression_mean[rows.wo.sz[i]]) + repro_PptFallCoef*PptFall[rows.wo.sz[i]]   
                                   + repro_PptSummerCoef*PptSummer[rows.wo.sz[i]] + repro_TempWinterCoef*TempWinter[rows.wo.sz[i]] 
                                   + repro_TempFallCoef*TempFall[rows.wo.sz[i]] + repro_TempSummerCoef*TempSummer[rows.wo.sz[i]]
                                   + repro_Transect_randomeffect[TransectNew.num[rows.wo.sz[i]]])
    ## Inflorescences is either the observed or predicted inflorescence number to be used to make a total reproduction estimates to fit predicted inflorescences to new plants 
    ## For rows without measured inflorescences, this is just the predicted number, accounting for survival & reproduction yes/no
    ## Max statement enforces that if a plant must have been alive, its estimated inforescences should be counted in the total for estimating new plants 
    inflors[rows.wo.sz[i]] = max(Surv_mu[rows.wo.sz[i]], rows.wo.sz.alive[i])*repro_prob[rows.wo.sz[i]]*repro_amount[rows.wo.sz[i]]
  }
  #####################################################################
  
  
#####################################################################
  ## Only loop over rows with >0 inflorescences, to compare estimated inflorescences to observed
  for(i in 1:Nrows.w.inflors) {
 
   p.infls[i] <- r.infls/(r.infls + repro_amount[rows.w.inflors[i]])

   InflNew[rows.w.inflors[i]] ~ dnegbin(p.infls[i], r.infls)
    
  }
  #####################################################################
  
  #####################################################################
  ## Then, a loop that matches the estimated sizes with the observed sizes, just for rows with observed sizes
  for(i in 1:Ngrowcases){
    ## These lines describe the response distribution and linear model terms
    
    regression_residual[i] <- RosNew[goodgrowrows[i]] - regression_mean[goodgrowrows[i]]
    
    p[goodgrowrows[i]] <- r.growth[goodgrowrows[i]]/(r.growth[goodgrowrows[i]]+regression_mean[goodgrowrows[i]])
    
    RosNew[goodgrowrows[i]] ~ dnegbin( p[goodgrowrows[i]], r.growth[goodgrowrows[i]])
    
  } #End of going through cases for sizes
  #####################################################################

  #####################################################################
  ## Now a loop to do survival predictions, for the ending years with observations
  for(i in 1:Ncases){
    ## Fit the survival function:
    Survs[goodrows[i]] ~ dbern(Surv_mu[goodrows[i]])
  } #End loop to do survival predictions
#####################################################################
  
  
#####################################################################
  ## This loop is getting the predicted number of total inflorescences for each transect-year combo, & fitting to number of new plants a year later 
  ## The inprod & the 1st logical are making 0,1 for whether a given line is part of each transect-year combo
  for (i in 1:newPltlines){

       pred.tot.inflors[i] = inprod(yrtranscombo==newplt.yrtranscombo[i], inflors)
       
       mn.new.plts[i] = exp(newplt_intercept + log(pred.tot.inflors[i]+0.1))
       
       p.newplts[i] <- r.newplts/(r.newplts+mn.new.plts[i])
       newplts[i] ~  dnegbin(p.newplts[i], r.newplts)
      }
#####################################################################
  
  
  
  
  
########################################
########################################    
## These lines give the prior distributions for the parameters to be estimated
  #r_disp ~ dunif(0,50) 
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


## PROBABILITY OF REPRODUCTION
reproyesno_intercept ~  dnorm(0, 10^-6)
reproyesno_RosCoef ~  dnorm(0, 10^-6)
reproyesno_PptFallCoef ~  dnorm(0, 10^-6)
reproyesno_PptSummerCoef ~  dnorm(0, 10^-6)
reproyesno_TempFallCoef ~  dnorm(0, 10^-6)
reproyesno_TempSummerCoef ~  dnorm(0, 10^-6)
reproyesno_TempWinterCoef ~  dnorm(0, 10^-6)

## Reproduction probability transect random effects
for(TransectNew.num_iterator in 1:12){
  reproyesno_Transect_randomeffect[TransectNew.num_iterator] ~ dnorm(0, reproyesno_Transect_precision) }
reproyesno_Transect_precision ~ dunif(0,1)


# AMOUNT OF REPRODUCTION
r.infls ~ dunif(0,10) 
repro_intercept ~ dunif(-5,5) 
repro_RosCoef ~ dunif(0.5,1.5) 
repro_PptSummerCoef ~ dnorm(0, 10^-6)
repro_PptFallCoef ~ dnorm(0, 10^-6)
repro_TempWinterCoef ~ dnorm(0, 10^-6)   
repro_TempFallCoef ~ dnorm(0, 10^-6)  
repro_TempSummerCoef ~ dnorm(0, 10^-6)   

## Reproduction amount transect random effects
for(TransectNew.num_iterator in 1:12){
  repro_Transect_randomeffect[TransectNew.num_iterator] ~ dnorm(0, repro_Transect_precision) }
repro_Transect_precision ~ dunif(0,1)


## Negative binomial variance for predicted inflorescence number
inflors_precision ~ dgamma(0.001, 0.001)
var.inflors = 1/inflors_precision


## NEW PLANTS
newplt_intercept ~ dnorm(0, 10^-6)
r.newplts ~ dgamma(0.01,0.01)


#resid.sum.sq <- sum(regression_residual^2)
} # end of model specification 


# These lines are hooks to be read by runjags (they are ignored by JAGS):
#monitor# deviance,grwth_Transect_randomeffect,surv_Transect_randomeffect,reproyesno_Transect_randomeffect,repro_Transect_randomeffect,repro_Transect_precision,reproyesno_Transect_precision,surv_Transect_precision,grwth_Transect_precision,grwth_intercept,grwth_RosCoef,grwth_TempFallCoef,grwth_TempSummerCoef,grwth_TempWinterCoef,grwth_PptFallCoef,grwth_PptSummerCoef,grwth_PptWinterCoef,grwthvar_intercept,grwthvar_RosCoef,surv_intercept,surv_RosCoef,surv_PptWinterCoef,surv_TempFallCoef,surv_TempSummerCoef,surv_TempWinterCoef,reproyesno_intercept,reproyesno_RosCoef,reproyesno_PptFallCoef,reproyesno_PptSummerCoef,reproyesno_TempFallCoef,reproyesno_TempSummerCoef,reproyesno_TempWinterCoef,repro_intercept,repro_RosCoef,repro_PptFallCoef,repro_PptSummerCoef,repro_TempWinterCoef,repro_TempFallCoef,repro_TempSummerCoef,newplt_intercept,r.infls,r.newplts
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

  "reproyesno_intercept" <- -1
  "reproyesno_RosCoef" <- 0.1
  "reproyesno_TempFallCoef" <- 0.01
  "reproyesno_TempSummerCoef" <- 0.01
  "reproyesno_TempWinterCoef" <- 0.01   
  "reproyesno_PptFallCoef" <- 0.01
  "reproyesno_PptSummerCoef" <- 0.01  
  "reproyesno_Transect_precision" <- 0.01
 
 "repro_intercept" <- 2
 "repro_RosCoef" <- 1
 "repro_PptFallCoef" <- 0.01
 "repro_PptSummerCoef" <- 0.01
 "repro_TempWinterCoef" <- 0.01
 "repro_TempFallCoef" <- 0.01
 "repro_TempSummerCoef" <- 0.01
 "repro_Transect_precision" <- 0.01
 
 "newplt_intercept" <- -1
 "r.newplts" <- 1
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

  "reproyesno_intercept" <- -1
  "reproyesno_RosCoef" <- 1
  "reproyesno_TempFallCoef" <- 0.01
  "reproyesno_TempSummerCoef" <- 0.01
  "reproyesno_TempWinterCoef" <- 0.01   
  "reproyesno_PptFallCoef" <- 0.01
  "reproyesno_PptSummerCoef" <- 0.01  
  "reproyesno_Transect_precision" <- 0.01
  
  "repro_precision" <- 0.01
  "repro_intercept" <- 2
  "repro_RosCoef" <- 1
  "repro_PptFallCoef" <- 0.01
  "repro_PptSummerCoef" <- 0.01
  "repro_TempWinterCoef" <- 0.01
  "repro_TempFallCoef" <- 0.01
  "repro_TempSummerCoef" <- 0.01
  "repro_Transect_precision" <- 0.01
  
  "newplt_intercept" <- -0.1
  "r.newplts" <- 1
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

  "reproyesno_intercept" <- -5
  "reproyesno_RosCoef" <- 0.1
  "reproyesno_TempFallCoef" <- 0.01
  "reproyesno_TempSummerCoef" <- 0.01
  "reproyesno_TempWinterCoef" <- 0.01   
  "reproyesno_PptFallCoef" <- 0.01
  "reproyesno_PptSummerCoef" <- 0.01   
  "reproyesno_Transect_precision" <- 0.01
  
  "repro_precision" <- 0.01
  "repro_intercept" <- 2
  "repro_RosCoef" <- 0.5
  "repro_PptFallCoef" <- 0.01
  "repro_PptSummerCoef" <- 0.01
  "repro_TempWinterCoef" <- 0.01
  "repro_TempFallCoef" <- 0.01
  "repro_TempSummerCoef" <- 0.01
  "repro_Transect_precision" <- 0.01
  
  "newplt_intercept" <- 0
  "r.newplts" <- 1
}