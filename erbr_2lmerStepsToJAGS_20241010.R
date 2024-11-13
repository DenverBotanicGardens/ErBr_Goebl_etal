## Code used for analyses in Goebl et al. manuscript on Eriogonum brandegeii modeling
## Modify data and assign variables needed for JAGS model with data lags (missing years)
## Associated JAGS script models growth, survival, reproduction, and recruitment 
## Use run.jags to run associated JAGS script



rm(list=ls())
graphics.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(lme4)
library(ggplot2)
library(rjags)
library(runjags)
library(dplyr)
library(coda)
library(corrplot)
library(robustbase)
library(resample)
library(gplots)
library(matrixStats)
## ------------------------------------------------------------------------------------------------



## SET WD (WHERE JAGS SCRIPT IS LOCATED) ----------------------------------------------------------
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
dats <- read.csv("erbr_TagClust2022_20230408.csv", header = TRUE)
#erbr_TagClust4to13.csv
#erbr_TagClust4to13evn.csv
#erbr_TagClust4to13odd.csv
#erbr_TagClust4to8.csv
#erbr_TagClust9to13.csv
## ------------------------------------------------------------------------------------------------



## RESCALE ALL CLIMATE VARIABLES ------------------------------------------------------------------
clim32yrMAXES <- read.csv("erbr_climData3seas32yr_MAXES.csv", header=TRUE)
dats$Tot_fall_ppt	 =dats$Tot_fall_ppt/clim32yrMAXES$Tot_fall_ppt
dats$Tot_winter_ppt	= dats$Tot_winter_ppt/clim32yrMAXES$Tot_winter_ppt
dats$Tot_summer_ppt =	dats$Tot_summer_ppt/clim32yrMAXES$Tot_summer_ppt
dats$Mean_fall_temp	= dats$Mean_fall_temp/clim32yrMAXES$Mean_fall_temp
dats$Mean_winter_temp	=dats$Mean_winter_temp/clim32yrMAXES$Mean_winter_temp
dats$Mean_summer_temp =dats$Mean_summer_temp/clim32yrMAXES$Mean_summer_temp
## ------------------------------------------------------------------------------------------------





## The following lines assemble and clean the necessary input files for running the JAGS model. 
## These steps are particular to our study, but show the general steps needed to clean a data set 
## and check for errors that can prevent JAGS from running properly. 
## Please skip to line 185 to see a summary of the necessary input data structures and then the running of the JAGS model



## CHANGE NAMES OF CLIMATE COLUMNS ------------------------------------------------------------------
dats <- rename(dats, PptFall=Tot_fall_ppt, PptWinter=Tot_winter_ppt, PptSummer=Tot_summer_ppt,
               TempFall=Mean_fall_temp, TempWinter=Mean_winter_temp, TempSummer=Mean_summer_temp)
## ------------------------------------------------------------------------------------------------



## SETTING UP NECESSARY VARIABLES -----------------------------------------------------------------
## Setting up the jags model with lagged values
Nallrows <- length(dats$Site)

numyears <- length(unique(dats$Year))
numtrans <- length(unique(dats$TransectNew))

## Identify rows that are good dependent values (ending sizes) for survival or growth  
lagforsurv <- dats$lagforsurv       #Full list of lags or of -1 for first observation rows
goodrows <- which(dats$lagforsurv > 0)
goodgrowrows <- which(dats$lagsrtsz > 0)

## Identify the lag for these rows: how far back is the last good size measurement?
lagvals <- dats$lagforsurv[which(dats$lagforsurv > 0)]
Ncases <- length(goodrows)
Ngrowcases <- length(goodgrowrows)
Survs <- dats$surv
RosNew <- dats$RosNew
InflNew <-  dats$InflNew 
InflYesNo <- dats$InflYesNo


## Setting up variables for use in reproduction fitting: 
dats$InflYesNo <- dats$InflNew
dats$InflYesNo[dats$InflNew>1] <- 1
rows.w.sz <- which(is.na(dats$RosNew)==FALSE)
rows.wo.sz <- which(is.na(dats$RosNew)==TRUE)
Ndirectszcases <- length(rows.w.sz)           #Direct measures of size upon which to base reproduction
Nindirectszcases <- length(rows.wo.sz)        #No direct measures of size; reproduction to be inferred from estimated size 
rows.w.inflors <- which(dats$InflNew>0)       #Non-zero estimates of inflorescences so reproduction amount can be estimated, if reproductive
Nrows.w.inflors <- length(rows.w.inflors)


## Add vector to indicate if alive or dead after missing year(s)
dats$RowNum <- 1:nrow(dats)                           #Add a column to indicate row number
rows.wo.sz.alive <- as.data.frame(matrix(NA, nrow=length(rows.wo.sz), ncol=2))
colnames(rows.wo.sz.alive) <- c("Rows", "Alive")
rows.wo.sz.alive$Rows <- rows.wo.sz
 
for (ww in rows.wo.sz) {                                                  #Loop over all tags with 1 or more missing years                      
  tag.val <- dats$TagNew[ww]
  tag.each <- subset(dats, dats$TagNew==tag.val)                          #Process each tag 
  tag.surv <- tag.each$surv[!is.na(tag.each$surv) & tag.each$RowNum>ww]   #Store survival for 1st non-missing year post each missed year
  rows.wo.sz.alive$Alive[rows.wo.sz.alive$Rows==ww] <- tag.surv[1] 
}

rows.wo.sz.alive$Alive[is.na(rows.wo.sz.alive$Alive)] <- 0  #Change NAs to 0, these are lines where missed year was last & recorded as dead
rows.wo.sz.alive <- rows.wo.sz.alive$Alive                  #Change to vector



## Make transect & year values numerical to use in jags as random effects 
dats$TransectNew.num <- as.factor(dats$TransectNew)
dats$TransectNew.num <- as.numeric(dats$TransectNew.num)
TransectNew.num <- dats$TransectNew.num

dats$Year.num <- as.factor(dats$Year)
dats$Year.num <- as.numeric(dats$Year.num)   
Year.num <- dats$Year.num

## Make a linear index of transect-year combos
yrtranscombo=100*dats$TransectNew.num+dats$Year.num

## Set up a logical variable that is whether there is survival or growth data in a year (0,1): 
## this is needed for the summing up of inflorescence numbers to predict new plants
datayesno <- rep(1,length(dats$surv))
datayesno[which(is.na(dats$surv)==TRUE)] <- 0 
## ------------------------------------------------------------------------------------------------



## Make dataframe with new plants (that are likely recent seedlings) for each transect & year ------------
## Make dataframe that will hold data containing new plants
years <- unique(dats$Year.num)
years <- years[order(years)]
dats.newPlts <- as.data.frame(rep(unique(dats$TransectNew.num), each=length(years)))
colnames(dats.newPlts) <- "TransectNew.num"
dats.newPlts$Year.num <- rep(years)

## Identify new plants
newPlts <- dats %>% group_by(TagNew) %>% slice(which.min(Year))   #Identify rows with 1st appearance for each plant
newPlts <- newPlts[newPlts$Year!=2004,]                           #Remove 2004 (first year of data collection)
sz.cutoff <- 5                                                    #Size cutoff, above which plant was likely not recently a seedling 
newPlts <- newPlts[newPlts$RosNew < sz.cutoff,]                   #Remove if >X rosettes (these were likely missed and are not new)
num.newPlts <- newPlts %>% group_by(TransectNew.num, Year.num) %>% summarise(num.newPlts=n())  #Count number new plants per year & transect

## Add number of new plants to dataframe of each transect & year
dats.newPlts <- left_join(dats.newPlts, num.newPlts, by=c("TransectNew.num", "Year.num"))
dats.newPlts$num.newPlts[is.na(dats.newPlts$num.newPlts)] <- 0   #Change NAs (no new plants) to zeros
dats.newPlts$num.newPlts[dats.newPlts$Year.num==1] <- NA         #Change new plts in 2004 (year 1) to NA

## Add column so new plants in t+1 match year t
dats.newPlts <- dats.newPlts %>% mutate(num.newPlts1=lead(num.newPlts))  

## Add climate variables to new plants data 
dats.clim <- dats %>% dplyr::select(c(Year.num, PptFall, PptWinter, PptSummer, TempFall, TempWinter, TempSummer)) 
clim <- unique(dats.clim)
dats.newPlts <- dplyr::left_join(dats.newPlts, clim, by="Year.num")
dats.newPlts <- dats.newPlts %>% dplyr::mutate(PptFall1=lead(PptFall), PptWinter1=lead(PptWinter), PptSummer1=lead(PptSummer),
                                        TempFall1=lead(TempFall), TempWinter1=lead(TempWinter), TempSummer1=lead(TempSummer))  
## Note: variables with '1' at the end (e.g. num.newPlts1, PptFall1) represent t+1 data   


## Remove lines that correspond to transect-year combos that are not in the main data file
## Note: GP East Transect 6 & 7 do not have data from 2004, 2005, or 2006. 
dats.newPlts$yrtranscombo=100*dats.newPlts$TransectNew.num+dats.newPlts$Year.num
yrtrans.unq <- unique(yrtranscombo)
dats.newPlts <- dats.newPlts[dats.newPlts$yrtranscombo %in% yrtrans.unq,]

dats.newPlts <- dats.newPlts[is.na(dats.newPlts$num.newPlts1) == FALSE,]

newplts <- dats.newPlts$num.newPlts1
newplt.trans <- dats.newPlts$TransectNew.num
newplt.yr <- dats.newPlts$Year.num

newPltlines <- length(dats.newPlts$TransectNew.num)

## Make a linear index of transect-year combos for new plts
newplt.yrtranscombo=100*newplt.trans+newplt.yr 
## ------------------------------------------------------------------------------------------------


## In response to error, make sure any rows with data for RosNew don't have NA for InflNew and InflYesNo
dats[1217,]$InflNew <- 0
dats[1217,]$InflYesNo <- 0


## --------------------------------------------------------------------------
## List of required data structures that are used by the JAGS model ======
#'Ncases' and 'Ngrowcases' are the number of rows used for estimates of survival and growth 
#'goodrows' and 'goodgrowrows' are the rows used for estimates of survival and growth
#'lagvals' are the number of years since the last data collection   
#'TransectNew.num' is the transect number designation to be used as random effects
#'RosNew' are the number of rosettes used as a measure of plant size  
#'InflNew' are the number of inflorescences 
#'InflYesNo' is if a plant flowered (yes) or not (no)
#'Survs' is the 0 or 1 variable indicating if a plant is alive (1) or not (0)
#'rows.w.sz', 'rows.wo.sz', and 'rows.wo.sz.alive' are the rows with or without a measure of plant size
#'Ndirectszcases' and 'Nindirectszcases' are the number of rows with or without a measure of plant size
#'numyears' and 'numtrans' are the number of different years and transects in the study. 
#'rows.w.inflors' and 'Nrows.w.inflors' are the rows and number of rows with a non-zero measure of reproduction
#'yrtranscombo' are all the year-transect combinations
#'newplt.yrtranscombo' are all the year-transect combinations of new plants 
#'newplts' and 'newPltlines' are the rows and number of rows with the 1st appearance for each plant
#'PptFall','PptWinter','PptSummer','TempWinter','TempFall', and 'TempSummer' are the seasonal climate variables
###--------------------------------------------------------------------------



## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------
jags.mod <- run.jags('erbr_3JAGSmodBest.R', n.chains=3, data=dats, burnin=10000, thin=10, sample=30000, adapt=500, method='parallel')

## Save output
saveRDS(jags.mod, "erbr_JAGSmodBest_c3t10s30b10.rds")
#chains <- jags.mod$mcmc
#chains <- bind_rows(lapply(chains, as.data.frame))
#saveRDS(chains, file="chains.c3t10s30b10.rds")
## ------------------------------------------------------------------------------------------------



## LOOK AT MODEL OUTPUT ---------------------------------------------------------------------------
#jags.mod <- readRDS("erbr_JAGSmodBest_c3t10s30b10.rds")
summary(jags.mod)
plot(jags.mod)
summ.mod <- summary(jags.mod)
tail(summ.mod[,1:3], n=31)
gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)
## -----------------------------------------------------------------------------------------------

