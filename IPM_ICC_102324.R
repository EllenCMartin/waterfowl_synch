<<<<<<<< HEAD:NOPI.R
# Preliminary Bayesian IPM of waterfowl across Prairie Pothole Region Strata
# from early 2000s - 2020 in North America
library(jagsUI)
library(HDInterval)
library(jagshelper)

setwd("C:/Users/setashc/Documents/Projects/PopSynch_Waterfowl")
setwd("C:/Users/ema/Desktop/Waterfowl/1 NOPI IPM") ### ECM add own working directory



sink("ipm_ppr.jags")
cat("
model {

#########################################################
# 1. Define the priors and constraints for the parameters
#########################################################


## 1.1 Priors for initial abundances ##

for(s in 1:(n.strata)){      
     Nasy[s,1] ~ dgamma(n1[s], 0.001)  ### ECM totally made this up, need to check and change for a logical tau? THIS IS FOR DENSITY.
     #Nasy[s,1] ~ dpois(n1[s]) # This is for abundance.
   }



## 1.2 Priors and model constraints for survival and recovery

# recovery probabilities:
for (t in 1:(n.years)){                                            
for(s in 1:(n.strata)){
    logit(rd[s,t]) <- alpha.r[s] + b1[s]  # Direct recovery of AHY S
    logit(r[s,t]) <- alpha.r[s]  #indirect recovery of AHY birds
  } # t
} #s

# intercept for recovery rates
for(s in 1:(n.strata)){
tauslope1[s] <- pow(1.3, -2) # Northrup and Gerber 2018 prior for logit regression
alpha.r[s] ~ dnorm(0, tauslope1[s]) # beta coefficient for direct recovery of AHY 

# offset for direct recovery
tauslope2[s] <- pow(1.3, -2) # Northrup and Gerber 2018 prior for logit regression
b1[s] ~ dnorm(0, tauslope2[s]) # beta coefficient for direct recovery of AHY 
}




# mortality hazard rates
for(s in 1:(n.strata)){
for (t in 1:(n.years)){
     log(hr[s,t]) <- alpha0 + lin.alpha.s[s] + betat[t] + beta[s,t]     
}
     log(hrmean[s]) <- alpha0 + lin.alpha.s[s]     # mean hazard for stratum s
 }


# use strong prior on survival and transform to hazard scale
mu.band<-0.957
sd.band<-0.1                     
alpha.band<- ((mu.band^2)-(mu.band^3)-(mu.band*sd.band^2))/sd.band^2
beta.band<- (mu.band-(2*mu.band^2)+(mu.band^3)-(sd.band^2)+(mu.band*sd.band^2))/sd.band^2
band_surv ~ dbeta(alpha.band, beta.band)
alpha0 <-log(-log(band_surv))





# weakly informative prior on the hazard scale for the offset for each stratum
for(s in 1:(n.strata-1)){
alpha.s[s] ~ dexp(1) 
lin.alpha.s[s] <- log(alpha.s[s])
}
lin.alpha.s[n.strata]<-0 # last stratum is the baseline



# random time effect - shared across stratum
for(t in 1:(n.years)){
betat[t] ~ dnorm(0, tau.year.betat)
}
sig.year.betat ~ dt(0, pow(2.5, -2), 1)T(0,)  # Half-cauchy - use for random effect
tau.year.betat <- 1/ sig.year.betat^2
}



# random time effect - different for each stratum
for(s in 1:(n.strata)){
for(t in 1:(n.years)){
beta[s,t] ~ dnorm(0, tau.year)
}
sig.year ~ dt(0, pow(2.5, -2), 1)T(0,)  # Half-cauchy - use for random effect
tau.year <- 1/ sig.year^2
# beta[s,n.years]<-0 # use this for fixed effect and make for loop 1:n.years-1


# Intraclass correlation in survival
ICC[1] <- sig.year.betat^2 / ((sig.year.betat^2) + (sig.year^2))







# Derived parameters; this way the likelihoods do not have to be changed
for (s in 1:(n.strata)){
  for (t in 1:(n.years)){
    sm[s,t] <- exp(-(hr[s,t]))                 # monthly survival
    anns[s,t] <- pow(sm[s,t],12)               # AHY annual survival   
    annh[s,t] <- -log(anns[s,t])               # AHY annual hazard rate
    le[s,t] <- 1/(-log(anns[s,t]))             # life expectancy to use as a metric of life history strategy
  }
    smmean[s] <- exp(-(hrmean[s]))             # time-averaged monthly survival
    annsmean[s] <- pow(smmean[s],12)           # time-averaged annual survival   
}





## 1.3 Priors and model constraints for fecundity 
## (fledged young (females) per breeding female of a given age ##

# Priors and constraints
 # mean parameter values (regression coefficients)
 b0.mu ~ dnorm(0, 0.25) # intercepts
 b1.mu ~ dnorm(0, 0.25) # Site Latitude
 b2.mu ~ dnorm(0, 0.25) # Site Longitude

 # Priors for SD to describe stratum-specific variation in regression coefficients
 b0.sigma ~ dunif(0,5) # numerically indexed as above

 # convert SD to precision (tau) by inverting and squaring
 b0.tau <- pow(b0.sigma,-2)

 # priors for random year, year within stratum, and banding location (site) effects
 eta.site.sigma ~ dunif(0,5)
 eta.strat.sigma ~ dunif(0,5)
 eta.yr.sigma ~ dunif(0,5)
 eta.site.tau <- pow(eta.site.sigma,-2)
 eta.strat.tau <- pow(eta.strat.sigma,-2)
 eta.yr.tau <- pow(eta.yr.sigma,-2)

 for (s in 1:n.strata){ # where s indexes to the maximum number of strata
 b0[s] ~ dnorm(b0.mu, b0.tau) # generate stratum-specific regression coefficients

 for (t in 1:n.years){ # where t indexes to number of years surveyed in each stratum
 eta.yr[t] ~ dnorm(0,eta.yr.tau) # generate random year effects 
 
 eta.strat[s,t] ~ dnorm(0,eta.strat.tau) # generate random year effects within each stratum

 for (k in 1:maxsites){ # where k indexes to the maximum number of banding sites
 eta.site[s,t,k] ~ dnorm(0,eta.site.tau) # generate random banding site effects
 } # close sites loop
 } # close years loop
 } # close strata loop

# Intraclass correlation for fecundity
ICC[2] <- eta.yr.sigma^2 / ((eta.yr.sigma^2) + (eta.strat.sigma^2))



############################################
# 2. The likelihoods of the single data sets
############################################

# 2.1 Likelihood of count data and process of population dynamics
# 2.1.1 System process including demographic stochasticity
for (s in 1:(n.strata)){      
### Nfl[s,1] ~ dpois(f[s,1] * Nasy[s,1])    # fledged females at end of summer ## ECM FOR ABUNDANCE
### Nasy[s,2] ~ dbin(anns[s,1], (Nasy[s,1]+Nfl[s,1])) # assumption for first year ## ECM FOR ABUNDANCE


## WHEN USING DENSITY:
Nfl[s,1] <- (f[s,1] * Nasy[s,1])    # ECM CHANGED  fledged females at end of summer ## No stochasticity for density.
Nasy[s,2] <- (anns[s,1] * (Nasy[s,1]+Nfl[s,1]))  ## ECM FOR DENSITY # assumption for first year ## No stochasticity for density.
} # ECM close n.strata



for (s in 1:(n.strata)){   
for (t in 3:n.years){  
  Nfl[s,t-1] <- (f[s,t-1]*Nasy[s,t-1])
  Nasy[s,t] <-  ((pow(sm[s,t-2],3)*pow(sm[s,t-1],9)) * (Nfl[s,t-1] + Nasy[s,t-1]))
  } ### ECM close n.years
  Nfl[s,n.years] <- (f[s,n.years]*Nasy[s,n.years])  # last year since it's not estimated in loop above
} ### ECM close n.strata


for (s in 1:(n.strata)){    
for (t in 1:n.years){      
  N[s,t] <- Nasy[s,t] + Nfl[s,t]      
} ### ECM close n.years
} ### ECM close n.strata


# 2.1.2 Derived population growth rates
for (s in 1:(n.strata)){
for(t in 1:(n.years-1)){
 lambda[s,t] <- N[s,t+1]/N[s,t]
}
}


# 2.1.3 Observation model that incorporates external estimates
# of counts that are treated as data

for (s in 1:(n.strata)){    
for(t in 1:n.years){        
  ## y[s,t] ~ dnorm(N[s,t], pow(y.sd[t],-2)) ## ECM if the z.sd is based on sd across time.
   y[s,t] ~ dnorm(N[s,t], pow(y.sds[s],-2)) ## ECM if the y.sd is based on sd across strata. 

  } ### ECM close n.years
} ### ECM close n.strata



#-------------------------------------------------
# 2.2 The likelihoods of the Seber model m-array data sets
#-------------------------------------------------

### Code for looping over years with no band releases
rel_a_ind = which(relAHYf !=0, arr.ind=T) # indices of rows (column 1) and columns (column 2) in m-array that HAVE data
rel_a_ind_n = nrow(rel_a_ind) # the total number of elements in the m-array with data so you can loop over the ones that don't

# 2.2.1 Define the multinomial likelihood 
for(i in 1:rel_a_ind_n){
  marrAHYf[rel_a_ind[i,1],1:(n.years+1),rel_a_ind[i,2]] ~ dmulti(prAHYf[rel_a_ind[i,1],,rel_a_ind[i,2]],
                                                          relAHYf[rel_a_ind[i,1],rel_a_ind[i,2]])
}


# # 2.2.1 Define the multinomial likelihood 
# for(s in 1:(n.strata)){
# for (t in 1:(n.years)){
#   marrAHYf[t,1:(n.years+1),s] ~ dmulti(prAHYf[t,,s], relAHYf[t,s])
# }}

# 2.2.2 Define the cell probabilities of the m-array
# Main diagonal
for(s in 1:(n.strata)){
for (t in 1:(n.years)){
  prAHYf[t,t,s] <- (1 - anns[s,t]) * rd[s,t]

  # Above main diagonal
  for (j in (t+1):(n.years)){
    prAHYf[t,j,s] <- prod(anns[s,t:(j-1)]) * (1-anns[s,j]) * r[s,j]
  } #j

  # Below main diagonal
  for (j in 1:(t-1)){
    prAHYf[t,j,s] <- 0
  } #j
 } #t
} #s

# Last column: probability of non-recovery
for(s in 1:(n.strata)){
for (t in 1:(n.years)){
  prAHYf[t,n.years+1,s] <- 1-sum(prAHYf[t,1:n.years,s])
} #t
} #s



#-------------------------------------------------
# 2.3 The likelihood for fecundity
#-------------------------------------------------

# 2.3.1 Observation model

 for (i in 1:nrows){ # i indexes each line of data (from a banding site within a year within a stratum)
 # predict proportion of juveniles (Pjuv) based on stratum-specific random intercepts,
 # fixed effects (lat, long), and random effects of year, year/stratum and site/year/stratum
 
 logit(Pjuv[i]) <- b0[stratum[i]] + b1.mu*Lat[i] + b2.mu*Long[i] + eta.yr[yr[i]] +
                   eta.strat[stratum[i],yr[i]] + eta.site[stratum[i],yr[i],site[i]] 


 # 2.3.2 compare observed juveniles to predicted juveniles based on above model
 
 HY_bands[i] ~ dbinom(Pjuv[i],bands[i]) # Juv ~ Binomial(P=Pjuv, N=Juv+Ad)

 }

 
 # 2.3.3 year- and stratum-specific fecundity derived quantity
 
 for (s in 1:n.strata){ # where i indexes to the maximum number of strata
 for (t in 1:n.years){ # where t indexes to number of years surveyed in each stratum
    logit(Pjuv.strat[s,t]) <- b0[s] + eta.yr[s,t] + eta.strat[stratum[i],yr[i]] 

 
 f_init[s,t]<-Pjuv.strat[s,t]/(1-Pjuv.strat[s,t])  
 
 # 2.3.4 ECM: Setting upper limit on fecundity values: Within for loop
 f[s,t] <- ifelse(f_init[s,t]>6, 5.9, f_init[s,t])
 
 }
 }
 } ",fill = TRUE)
sink()








#################################################################################
################ Read in data and organize it to run model ######################
#################################################################################


###############################
# survival/dead recovery data #
###############################

#*********************** BANDING M-ARRAYS **************************************
# band recovery m-arrays with numbers 'not recovered' in the final column,
# arranged by age, sex, and season release classes

## NOPI
#load("megamarraynopi_synch.RData")
load("megamarraynopi_synch_density.RData")## density marray 
dim(megamarray) # rows (cohort), columns (recovery years and final column of never-recovered birds), stratum, age


marrHYf<-array(NA, dim=c(dim(megamarray)[1], dim(megamarray)[2], dim(megamarray)[3]))
marrAHYf<-array(NA, dim=c(dim(megamarray)[1], dim(megamarray)[2], dim(megamarray)[3]))




n.strata<-dim(megamarray)[3]
for(s in 1:n.strata){
  marrHYf[,,s]<- megamarray[,,s,1]
  marrAHYf[,,s]<- megamarray[,,s,2]
  
}

relHYf<-matrix(NA, dim(marrHYf)[1], n.strata)
relAHYf<-matrix(NA, dim(marrAHYf)[1], n.strata)

for(s in 1:n.strata){
  relHYf[,s]<- rowSums(marrAHYf[,,s])
  relAHYf[,s]<- rowSums(marrAHYf[,,s])
  
}





##################
# fecundity data #
##################

load("nopi_f_dat.RData")
dat<-subset(dat, year<2020) # abundance and banding only go up to 2019
strat<-read.csv("nopi.strat.f.csv")
dat<-merge(dat, strat, by="stratno", all.y=TRUE)
for(k in 2:5){
  dat[is.na(dat[,k]), k]<-mean(dat[,k], na.rm=TRUE) #replace missing values with column means (lat, long, site, year)
}

for(k in 7:8){
  dat[is.na(dat[,k]), k]<-mean(dat[,k], na.rm=TRUE) #replace missing values with column means (no. HY bands and no. total bands)
}

# round the categorical/integer column means
dat$year<-round(dat$year)
dat$site<-round(dat$site)
dat$HY_bands<-round(dat$HY_bands)
dat$bands<-round(dat$bands)

dat<-subset(dat, year>1999) # abundance and banding only go up to 2019

library(tidyverse)
dat %>% select('year','STRAT') %>% ftable()


##################
# abundance data #
##################

load("NOPI_DATA_ABUNDANCE.Rda") ### ECM import data. object name "DATA"
load("NOPI_DATA_ABUNDANCE_MCMC.Rda") ### ECM import mcmc gssm data. object name "mcmcdf"
load("NOPI_DATA_DENSITY.Rda") ### ECM import data. object name "DATA"
load("NOPI_DATA_DENSITY_MCMC.Rda") ### ECM import data. object name "DATA"

# y <- mcmcdf ### ECM THIS IS FOR DENSITY. Set one of these as the data of choice for the analysis. Either DATA or mcmcdf
y <- mcmcdf ### ECM THIS IS FOR nonMCMC density or abundance



# abundance/density count data
n1 <- y[,2] # ECM. Old CS comment: use the count data to make this the initial pop size of adult females 
y.sd <- as.vector(apply(y[,2:ncol(y)], 2, sd)) # CS - I changed this to the SD across strata for each year for now, 
y.sds <- as.vector(apply(y[,2:ncol(y)], 1, sd)) # ECM - SD within strata across years. 
# although it's really supposed to represent observation error so it would have to come from some kind of double observer info or something



# Bundle data: n.years = number of years for survival dataset starting in 2000
# df data should be 1 dimension > dim of parameter matrix. R data provide
# priors to process variance and could be multiplied by 2 to give more
# vague prior information.

# remove stratrum 23 (actually 76)
y<-y[-23,]
n1<-n1[-23]
y.sd <- as.vector(apply(y[,2:ncol(y)], 2, sd)) # CS - I changed this to the SD across strata for each year for now, 
dat<-subset(dat, stratno!=23)
marrAHYf<-marrAHYf[,,-23]
relAHYf<-relAHYf[,-23]

###replace all zeros in year they were surveyed with the mean value per row.

y.sd <- as.vector(apply(y[,2:ncol(y)], 2, sd)) # CS - I changed this to the SD across strata for each year for now, 

## remove all years before 2000 in marray and fecundity stuff



jags.data <- list(n1=n1, y=as.matrix(y[2:ncol(y)]), y.sd=y.sd, # abundance
                  
                  y.sds=y.sds,  
                  n.years = ncol(marrAHYf)-1, marrAHYf = marrAHYf, # survival
                  relAHYf = relAHYf, n.strata=ncol(relAHYf),
                  
                  stratum=dat$stratno, 
                  yr=as.numeric(as.factor(dat$year)),  # fecundity
                  
                  site=dat$site, maxsites=max(dat$site, na.rm = TRUE), 
                  Lat=scale(dat$lat)[,], Long=scale(dat$long)[,],
                  HY_bands=dat$HY_bands, bands=dat$bands,
                  nrows=nrow(dat))




# Initial values: 
inits <- function(){list()}

# Parameters monitored
parameters <- c("anns", "f", "sm", "hr", "b0", "b1", "alpha.r", "lin.alpha.s", "alpha0",
                "rd", "r",  "annh", "beta", "annsmean", "N", "Nasy", "Pjuv.strat", "f_init", "Nfl",
                "et.f","et.hr","ets.f","ets.hr","ICC", "var.tt","var.ts","tau.tt", "tau.ts")

# MCMC settings
ni <- 10000
nt <- 1000
nb <- 500
nc <- 3

# Call JAGS from R 
set.seed(runif(1,1,900))  # sometimes this has to be reset


nopi_ipm <- jags(jags.data, inits, parameters,"ipm_ppr.jags",n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

nopi_ipm$Rhat 
traceplot(nopi_ipm, parameters=c("b1"))
traceplot(nopi_ipm, parameters=c("f"))
traceplot(nopi_ipm, parameters=c("anns")) # a lot of times derived quantities say they haven't converged when they have,
# so I never trust traceplots/R-hats for them and try to just look at ACTUAL parameters

save(nopi_ipm, file="nopi_ipm_density_04062024")


plotRhats(nopi_ipm)



################# SOME PLOTS OF SAVED PARAMETERS ##########################
library(RColorBrewer)
nstrata <- 22
n <- nstrata
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
qu <- function(x) quantile(x, c(0.025, 0.975), na.rm=TRUE)
ptc <- 1.2
op <- par(las=1, mfrow=c(2, 1), mar=c(2,4,2,1))
years <- 20



# PLOT 1:  Annual survival "anns"
plot(x=seq(1:years), y=nopi_ipm$mean$anns[1,], type='b', pch=16, ylab='Annual Mean Survival',
     ylim=c(0, 1), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector[1], cex=ptc)
segments(seq(1:years), nopi_ipm$q2.5$anns[1,], seq(1:years), nopi_ipm$q97.5$anns[1,], col=col_vector[1])
for (s in 1:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$anns[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(seq(1:years), nopi_ipm$q2.5$anns[s,], seq(1:years), nopi_ipm$q97.5$anns[s,], col=col_vector[s])
}
axis(2)
axis(1, at=1:30)


# PLOT 2:  Monthly? survival "sm"
plot(x=seq(1:years), y=nopi_ipm$mean$sm[1,], type='b', pch=16, ylab='sm',
     ylim=c(0, 1), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector[1], cex=ptc)
segments(seq(1:years), nopi_ipm$q2.5$sm[1,], seq(1:years), nopi_ipm$q97.5$sm[1,], col=col_vector[1])
for (s in 1:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$sm[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(seq(1:years), nopi_ipm$q2.5$sm[s,], seq(1:years), nopi_ipm$q97.5$sm[s,], col=col_vector[s])
}
axis(2)
axis(1, at=1:30)


# PLOT 3: fecundity by strata 
plot(x=seq(1:years), y=nopi_ipm$mean$Nfl[1,], type='b', pch=16, ylab='f',
     ylim=c(0, 5), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector[1], cex=ptc)
segments(seq(1:years), nopi_ipm$q2.5$Nfl[1,], seq(1:years), nopi_ipm$q97.5$Nfl[1,], col=col_vector[1])
for (s in 1:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$Nfl[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(seq(1:years), nopi_ipm$q2.5$Nfl[s,], seq(1:years), nopi_ipm$q97.5$Nfl[s,], col=col_vector[s])
}
axis(2)
axis(1, at=1:30)



# PLOT 4:  "rd" direct recovery (first year after release) per strata
plot(x=seq(1:nstrata), y=nopi_ipm$mean$rd[,1], type='b', pch=16, ylab='rd',
     ylim=c(0, 1), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector, cex=ptc)
segments(seq(1:nstrata), nopi_ipm$q2.5$rd[,1], seq(1:nstrata), nopi_ipm$q97.5$rd[,1], col=col_vector)
axis(2)
axis(1, at=1:23)


# PLOT 5:  "r" indirect recovery rates (prob of being recovered and reported) after first year after release per strata
plot(x=seq(1:nstrata), y=nopi_ipm$mean$r[,1], type='b', pch=16, ylab='r',
     ylim=c(0, 1), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector, cex=ptc)
segments(seq(1:nstrata), nopi_ipm$q2.5$r[,1], seq(1:nstrata), nopi_ipm$q97.5$r[,1], col=col_vector)
axis(2)
axis(1, at=1:23)


# PLOT 6:  hazard rate "annh"
plot(x=seq(1:years), y=nopi_ipm$mean$annh[1,], type='b', pch=16, ylab='annh',
     ylim=c(0, 5), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector[1], cex=ptc)
segments(seq(1:years), nopi_ipm$q2.5$annh[1,], seq(1:years), nopi_ipm$q97.5$annh[1,], col=col_vector[1])
for (s in 1:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$annh[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(seq(1:years), nopi_ipm$q2.5$annh[s,], seq(1:years), nopi_ipm$q97.5$annh[s,], col=col_vector[s])
}
axis(2)
axis(1, at=1:30)


# Annual survival by strata "annsmean"
plot(x=1, y=nopi_ipm$mean$annsmean[1], type='b', pch=16, ylab='Annual Mean Survival',
     ylim=c(0, 1), xlab=NA, xlim=c(0, nstrata), axes=FALSE, col=col_vector[1], cex=ptc)
for (s in 1:nstrata){
  points(x=s, y=nopi_ipm$mean$annsmean[s], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(s, nopi_ipm$q2.5$annsmean[s], s, nopi_ipm$q97.5$annsmean[s],  col=col_vector[s])
}
axis(2)
axis(1, at=1:30)



# PLOT 1:  Annual all Rhats "anns"
plot(x=seq(1:years), y=nopi_ipm$Rhat$anns[1,], type='b', pch=16, ylab='Rhat Annual Mean Survival',
     ylim=c(0.95, 1.25), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$anns[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)


# PLOT 2:  Annual all Rhats "f"
plot(x=seq(1:years), y=nopi_ipm$Rhat$f[1,], type='b', pch=16, ylab='Rhat Fecundity',
     ylim=c(0.95, 2), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$f[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)


# PLOT 3:  Annual all Rhats "sm"
plot(x=seq(1:years), y=nopi_ipm$Rhat$sm[1,], type='b', pch=16, ylab='Rhat Monthly Survival',
     ylim=c(0.9, 2), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$sm[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)


# PLOT 3:  Annual all Rhats "N"
plot(x=seq(1:years), y=nopi_ipm$Rhat$N[1,], type='b', pch=16, ylab='Rhat N',
     ylim=c(0.95, 1.3), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$N[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)



# PLOT 3:  Annual all Rhats "N"
plot(x=seq(1:years), y=nopi_ipm$mean$N[1,], type='b', pch=16, ylab='Strata Density Estimate',
     ylim=c(0, 7), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$N[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)




# PLOT 3:  Annual all Rhats "N"
plot(x=seq(1:years), y=nopi_ipm$Rhat$N[1,], type='b', pch=16, ylab='Strata Density Estimate',
     ylim=c(0.95, 1.4), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$N[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)



#### CALCULATE CROSS CORRELATIONS

CC_survival <- rcorr(t(nopi_ipm$mean$anns))
CC_survival <- CC_survival$r
TF<- lower.tri(CC_survival, diag = FALSE)
CC_survival <- CC_survival*TF
CC_survival[CC_survival==0] = NA
CC_survival[CC_survival=="NaN"] = NA


CC_fec <-  rcorr(t(nopi_ipm$mean$f))
CC_fec <- CC_fec$r
TF<- lower.tri(CC_fec, diag = FALSE)
CC_fec <- CC_fec*TF
CC_fec[CC_fec==0] = NA
CC_fec[CC_fec=="NaN"] = NA


## calculate distances
library(sf)
library(geosphere)
setwd("C:/Users/ema/Desktop/Waterfowl")
## Define projections for sf objects. Layers are EPSG and WGS84: want to transform to UTM33
lonlatproj <- "+proj=longlat +unit=dd +datum=WGS84"
test <- "EPSG:4326"

BirdData <- read.csv("wbphs_stratum_densities/wbphs_stratum_densities_forDistribution.csv")# Load the bird location data
BirdData[is.na(BirdData)] <- 0 # Include NA as 0, but that is likely not correct

BirdLocs <- read_sf(dsn = "WBPHS_Stratum_Boundaries/WBPHS_Stratum_Boundaries.shp")
BirdLocsLONLAT <- st_transform(BirdLocs, crs=lonlatproj) # Transform Kantons to UTM

StrataLocs <- read.csv("WBPHS_Segment_Counts/wbphs_segment_effort_forDistribution.csv")# Load the bird location data
StrataLocs$midpoint_latitude<- as.numeric(StrataLocs$midpoint_latitude)
StrataLocs$midpoint_longitude<- as.numeric(StrataLocs$midpoint_longitude)
StrataLocs <- na.omit(StrataLocs)
StrataLocs_sf <- st_as_sf(x=StrataLocs, coords=c("midpoint_longitude", "midpoint_latitude"), crs=test) # Make it an sf object
StrataLocs_sf <- st_transform(StrataLocs_sf, crs=lonlatproj) # Transform Kantons to UTM


## Calculate mean survey location coordinates per strata
Coords <- aggregate(st_coordinates(StrataLocs_sf), 
                    by=list(stratumcoords=data.frame(StrataLocs_sf)[,"stratum"]), FUN=mean)
Coordssf <- st_as_sf(x=Coords, coords=c("X", "Y"), crs=lonlatproj) # Make it an sf object
distmat <- distm(Coords[,2:3], fun=distGeo)/1000

disttab<-disttab[-23,]
disttab<-disttab[,-23]
========
# Preliminary Bayesian IPM of waterfowl across Prairie Pothole Region Strata
# from early 2000s - 2020 in North America
library(jagsUI)
library(HDInterval)
library(jagshelper)

setwd("C:/Users/setashc/Documents/Projects/PopSynch_Waterfowl")
setwd("C:/Users/ema/Desktop/Waterfowl/1 NOPI IPM") ### ECM add own working directory



sink("ipm_ppr.jags")
cat("
model {

#########################################################
# 1. Define the priors and constraints for the parameters
#########################################################


## 1.1 Priors for initial abundances ##

for(s in 1:(n.strata)){      
     Nasy[s,1] ~ dgamma(n1[s], 0.001)  ### ECM totally made this up, need to check and change for a logical tau? THIS IS FOR DENSITY.
     #Nasy[s,1] ~ dpois(n1[s]) # This is for abundance.
   }



## 1.2 Priors and model constraints for survival and recovery

# recovery probabilities:
for (t in 1:(n.years)){                                            
for(s in 1:(n.strata)){
    logit(rd[s,t]) <- alpha.r[s] + b1[s]  # Direct recovery of AHY S
    logit(r[s,t]) <- alpha.r[s]  #indirect recovery of AHY birds
  } # t
} #s

# intercept for recovery rates
for(s in 1:(n.strata)){
tauslope1[s] <- pow(1.3, -2) # Northrup and Gerber 2018 prior for logit regression
alpha.r[s] ~ dnorm(0, tauslope1[s]) # beta coefficient for direct recovery of AHY 

# offset for direct recovery
tauslope2[s] <- pow(1.3, -2) # Northrup and Gerber 2018 prior for logit regression
b1[s] ~ dnorm(0, tauslope2[s]) # beta coefficient for direct recovery of AHY 
}

# # Priors for temp rand eff:
# for (t in 1:(n.years)){                                            
#      et.hr[t] ~ dnorm(0, tau.tt[1])
#      et.f[t] ~ dnorm(0, tau.tt[2])
# } # t
# 
# # Priors for within site temp rand eff:
# for(s in 1:(n.strata)){
# for (t in 1:(n.years)){                                            
#       ets.hr[s,t] ~ dnorm(0, tau.ts[1])
#       ets.f[s,t] ~ dnorm(0, tau.ts[2])
# } # t nyears for fecundity
# } # s
# 
# 
# # Priors for precisions
#      for (i in 1:2){
#         tau.tt[i] ~ dgamma(0.001, 0.001)
#         tau.ts[i] ~ dgamma(0.001, 0.001)
#         
#         # Variances
#         var.tt[i] <- 1 / tau.tt[i]
#         var.ts[i] <- 1 / tau.ts[i]
# 
#         # Intraclass correlation 
#         ICC[i] <- var.tt[i] / ( var.tt[i] + var.ts[i])
#         } # i



# mortality hazard rates
for(s in 1:(n.strata)){
for (t in 1:(n.years)){
     log(hr[s,t]) <- alpha0 + lin.alpha.s[s] + betat[t] + beta[s,t]     # with random time effect
  #  log(hr[s,t]) <- meanhr[s] + et.hr[t] + ets.hr[s,t]     # ECM NEW FOR SYNCHRONY CALC: mean + temp rand eff + within site temp ran eff
}
     log(hrmean[s]) <- alpha0 + lin.alpha.s[s]     # mean hazard for stratum s
 }


# use strong prior on survival and transform to hazard scale
mu.band<-0.957
sd.band<-0.1                     ## ECM changed this from 0.005 to 0.1 on May 27 2024
alpha.band<- ((mu.band^2)-(mu.band^3)-(mu.band*sd.band^2))/sd.band^2
beta.band<- (mu.band-(2*mu.band^2)+(mu.band^3)-(sd.band^2)+(mu.band*sd.band^2))/sd.band^2
band_surv ~ dbeta(alpha.band, beta.band)
alpha0 <-log(-log(band_surv))

# # Mean intercept per stratum
# for (s in 1:n.strata){
#   meanhr[s] ~ dnorm(alpha0, 0.001)
# } 


# weakly informative prior on the hazard scale for the offset for each stratum
for(s in 1:(n.strata-1)){
alpha.s[s] ~ dexp(1) 
lin.alpha.s[s] <- log(alpha.s[s])
}
lin.alpha.s[n.strata]<-0 # last stratum is the baseline



# random time effect - shared across stratum
for(t in 1:(n.years)){
betat[t] ~ dnorm(0, tau.year.betat)
}
sig.year.betat ~ dt(0, pow(2.5, -2), 1)T(0,)  # Half-cauchy - use for random effect
tau.year.betat <- 1/ sig.year.betat^2
}



# random time effect - different for each stratum
for(s in 1:(n.strata)){
for(t in 1:(n.years)){
beta[s,t] ~ dnorm(0, tau.year)
}
sig.year ~ dt(0, pow(2.5, -2), 1)T(0,)  # Half-cauchy - use for random effect
tau.year <- 1/ sig.year^2
# beta[s,n.years]<-0 # use this for fixed effect and make for loop 1:n.years-1


ICC[1] <- sig.year.betat^2 / ((sig.year.betat^2) + (sig.year^2))







# Derived parameters; this way the likelihoods do not have to be changed
for (s in 1:(n.strata)){
  for (t in 1:(n.years)){
    sm[s,t] <- exp(-(hr[s,t]))                 # monthly survival
    anns[s,t] <- pow(sm[s,t],12)               # AHY annual survival   
    annh[s,t] <- -log(anns[s,t])               # AHY annual hazard rate
    le[s,t] <- 1/(-log(anns[s,t]))             # life expectancy to use as a metric of life history strategy
  }
    smmean[s] <- exp(-(hrmean[s]))             # time-averaged monthly survival
    annsmean[s] <- pow(smmean[s],12)           # time-averaged annual survival   
}





## 1.3 Priors and model constraints for fecundity 
## (fledged young (females) per breeding female of a given age ##

# Priors and constraints
 # mean parameter values (regression coefficients)
 b0.mu ~ dnorm(0, 0.25) # intercepts
 b1.mu ~ dnorm(0, 0.25) # Site Latitude
 b2.mu ~ dnorm(0, 0.25) # Site Longitude

 # Priors for SD to describe stratum-specific variation in regression coefficients
 b0.sigma ~ dunif(0,5) # numerically indexed as above

 # convert SD to precision (tau) by inverting and squaring
 b0.tau <- pow(b0.sigma,-2)

 # priors for random year and banding location (site) effects
 eta.yr.sigma ~ dunif(0,5)
 eta.site.sigma ~ dunif(0,5)
 eta.yr.tau <- pow(eta.yr.sigma,-2)
 eta.site.tau <- pow(eta.site.sigma,-2)

 for (s in 1:n.strata){ # where s indexes to the maximum number of strata
 b0[s] ~ dnorm(b0.mu, b0.tau) # generate stratum-specific regression coefficients

 for (t in 1:n.years){ # where t indexes to number of years surveyed in each stratum
 eta.yr[s,t] ~ dnorm(0,eta.yr.tau) # generate random year effects within each stratum

 for (k in 1:maxsites){ # where k indexes to the maximum number of banding sites
 eta.site[s,t,k] ~ dnorm(0,eta.site.tau) # generate random banding site effects
 } # close sites loop
 } # close years loop
 } # close strata loop

ICC[2] <- eta.yr.sigma^2 / ((eta.yr.sigma^2) + (eta.site.sigma^2))



############################################
# 2. The likelihoods of the single data sets
############################################

# 2.1 Likelihood of count data and process of population dynamics
# 2.1.1 System process including demographic stochasticity
for (s in 1:(n.strata)){      
### Nfl[s,1] ~ dpois(f[s,1] * Nasy[s,1])    # fledged females at end of summer ## ECM FOR ABUNDANCE
### Nasy[s,2] ~ dbin(anns[s,1], (Nasy[s,1]+Nfl[s,1])) # assumption for first year ## ECM FOR ABUNDANCE


## WHEN USING DENSITY:
Nfl[s,1] <- (f[s,1] * Nasy[s,1])    # ECM CHANGED  fledged females at end of summer ## No stochasticity for density.
Nasy[s,2] <- (anns[s,1] * (Nasy[s,1]+Nfl[s,1]))  ## ECM FOR DENSITY # assumption for first year ## No stochasticity for density.
} # ECM close n.strata



for (s in 1:(n.strata)){   
for (t in 3:n.years){  
  Nfl[s,t-1] <- (f[s,t-1]*Nasy[s,t-1])
  Nasy[s,t] <-  ((pow(sm[s,t-2],3)*pow(sm[s,t-1],9)) * (Nfl[s,t-1] + Nasy[s,t-1]))
  } ### ECM close n.years
  Nfl[s,n.years] <- (f[s,n.years]*Nasy[s,n.years])  # last year since it's not estimated in loop above
} ### ECM close n.strata


for (s in 1:(n.strata)){    
for (t in 1:n.years){      
  N[s,t] <- Nasy[s,t] + Nfl[s,t]      
} ### ECM close n.years
} ### ECM close n.strata


# 2.1.2 Derived population growth rates
for (s in 1:(n.strata)){
for(t in 1:(n.years-1)){
 lambda[s,t] <- N[s,t+1]/N[s,t]
}
}


# 2.1.3 Observation model that incorporates external estimates
# of counts that are treated as data

for (s in 1:(n.strata)){    
for(t in 1:n.years){        
  ## y[s,t] ~ dnorm(N[s,t], pow(y.sd[t],-2)) ## ECM if the z.sd is based on sd across time.
   y[s,t] ~ dnorm(N[s,t], pow(y.sds[s],-2)) ## ECM if the y.sd is based on sd across strata. 

  } ### ECM close n.years
} ### ECM close n.strata



#-------------------------------------------------
# 2.2 The likelihoods of the Seber model m-array data sets
#-------------------------------------------------


# 2.2.1 Define the multinomial likelihood 
for(s in 1:(n.strata)){
for (t in 1:(n.years)){
  marrAHYf[t,1:(n.years+1),s] ~ dmulti(prAHYf[t,,s], relAHYf[t,s]) #band only
}}

# 2.2.2 Define the cell probabilities of the m-array
# Main diagonal
for(s in 1:(n.strata)){
for (t in 1:(n.years)){
  prAHYf[t,t,s] <- (1 - anns[s,t]) * rd[s,t]

  # Above main diagonal
  for (j in (t+1):(n.years)){
    prAHYf[t,j,s] <- prod(anns[s,t:(j-1)]) * (1-anns[s,j]) * r[s,j]
  } #j

  # Below main diagonal
  for (j in 1:(t-1)){
    prAHYf[t,j,s] <- 0
  } #j
 } #t
} #s

# Last column: probability of non-recovery
for(s in 1:(n.strata)){
for (t in 1:(n.years)){
  prAHYf[t,n.years+1,s] <- 1-sum(prAHYf[t,1:n.years,s])
} #t
} #s



#-------------------------------------------------
# 2.3 The likelihood for fecundity
#-------------------------------------------------

# 2.3.1 Observation model
 for (i in 1:nrows){ # i indexes each line of data (from a banding site within a year within a stratum)
 # predict proportion of juveniles (Pjuv) based on stratum-specific random covariates,
 # fixed effects (lat, long), and random effects of year/stratum and site/year/stratum
 
 logit(Pjuv[i]) <- b0[stratum[i]] + b1.mu*Lat[i] + b2.mu*Long[i] + 
                   eta.yr[stratum[i],yr[i]] + eta.site[stratum[i],yr[i],site[i]] 


 # 2.3.2 compare observed juveniles to predicted juveniles based on above model
 
 HY_bands[i] ~ dbinom(Pjuv[i],bands[i]) # Juv ~ Binomial(P=Pjuv, N=Juv+Ad)

 }

 
 # 2.3.3 year- and stratum-specific fecundity derived quantity
 
 for (s in 1:n.strata){ # where i indexes to the maximum number of strata
 for (t in 1:n.years){ # where t indexes to number of years surveyed in each stratum
    logit(Pjuv.strat[s,t]) <- b0[s] + eta.yr[s,t]
  # logit(Pjuv.strat[s,t]) <- b0[s] + et.f[t] + ets.f[s,t] # ECM NEW   
 
 
 f_init[s,t]<-Pjuv.strat[s,t]/(1-Pjuv.strat[s,t])  
 
 # 2.3.4 ECM: Setting upper limit on fecundity values: Within for loop
 f[s,t] <- ifelse(f_init[s,t]>6, 5.9, f_init[s,t])
 
 }
 }
 } ",fill = TRUE)
sink()








#################################################################################
################ Read in data and organize it to run model ######################
#################################################################################


###############################
# survival/dead recovery data #
###############################

#*********************** BANDING M-ARRAYS **************************************
# band recovery m-arrays with numbers 'not recovered' in the final column,
# arranged by age, sex, and season release classes

## NOPI
#load("megamarraynopi_synch.RData")
load("megamarraynopi_synch_density.RData")## density marray 
dim(megamarray) # rows (cohort), columns (recovery years and final column of never-recovered birds), stratum, age


marrHYf<-array(NA, dim=c(dim(megamarray)[1], dim(megamarray)[2], dim(megamarray)[3]))
marrAHYf<-array(NA, dim=c(dim(megamarray)[1], dim(megamarray)[2], dim(megamarray)[3]))




n.strata<-dim(megamarray)[3]
for(s in 1:n.strata){
  marrHYf[,,s]<- megamarray[,,s,1]
  marrAHYf[,,s]<- megamarray[,,s,2]
  
}

relHYf<-matrix(NA, dim(marrHYf)[1], n.strata)
relAHYf<-matrix(NA, dim(marrAHYf)[1], n.strata)

for(s in 1:n.strata){
  relHYf[,s]<- rowSums(marrAHYf[,,s])
  relAHYf[,s]<- rowSums(marrAHYf[,,s])
  
}



##################
# fecundity data #
##################

load("nopi_f_dat.RData")
dat<-subset(dat, year<2020) # abundance and banding only go up to 2019
strat<-read.csv("nopi.strat.f.csv")
dat<-merge(dat, strat, by="stratno", all.y=TRUE)
for(k in 2:5){
  dat[is.na(dat[,k]), k]<-mean(dat[,k], na.rm=TRUE) #replace missing values with column means (lat, long, site, year)
}

for(k in 7:8){
  dat[is.na(dat[,k]), k]<-mean(dat[,k], na.rm=TRUE) #replace missing values with column means (no. HY bands and no. total bands)
}

# round the categorical/integer column means
dat$year<-round(dat$year)
dat$site<-round(dat$site)
dat$HY_bands<-round(dat$HY_bands)
dat$bands<-round(dat$bands)

dat<-subset(dat, year>1999) # abundance and banding only go up to 2019

library(tidyverse)
dat %>% select('year','STRAT') %>% ftable()


##################
# abundance data #
##################

load("NOPI_DATA_ABUNDANCE.Rda") ### ECM import data. object name "DATA"
load("NOPI_DATA_ABUNDANCE_MCMC.Rda") ### ECM import mcmc gssm data. object name "mcmcdf"
load("NOPI_DATA_DENSITY.Rda") ### ECM import data. object name "DATA"
load("NOPI_DATA_DENSITY_MCMC.Rda") ### ECM import data. object name "DATA"

# y <- mcmcdf ### ECM THIS IS FOR DENSITY. Set one of these as the data of choice for the analysis. Either DATA or mcmcdf
y <- mcmcdf ### ECM THIS IS FOR nonMCMC density or abundance



# abundance/density count data
n1 <- y[,2] # ECM. Old CS comment: use the count data to make this the initial pop size of adult females 
y.sd <- as.vector(apply(y[,2:ncol(y)], 2, sd)) # CS - I changed this to the SD across strata for each year for now, 
y.sds <- as.vector(apply(y[,2:ncol(y)], 1, sd)) # ECM - SD within strata across years. 
# although it's really supposed to represent observation error so it would have to come from some kind of double observer info or something



# Bundle data: n.years = number of years for survival dataset starting in 2000
# df data should be 1 dimension > dim of parameter matrix. R data provide
# priors to process variance and could be multiplied by 2 to give more
# vague prior information.

# remove stratrum 23 (actually 76)
y<-y[-23,]
n1<-n1[-23]
y.sd <- as.vector(apply(y[,2:ncol(y)], 2, sd)) # CS - I changed this to the SD across strata for each year for now, 
dat<-subset(dat, stratno!=23)
marrAHYf<-marrAHYf[,,-23]
relAHYf<-relAHYf[,-23]

###replace all zeros in year they were surveyed with the mean value per row.

y.sd <- as.vector(apply(y[,2:ncol(y)], 2, sd)) # CS - I changed this to the SD across strata for each year for now, 

## remove all years before 2000 in marray and fecundity stuff



jags.data <- list(n1=n1, y=as.matrix(y[2:ncol(y)]), y.sd=y.sd, # abundance
                  
                  y.sds=y.sds,  
                  n.years = ncol(marrAHYf)-1, marrAHYf = marrAHYf, # survival
                  relAHYf = relAHYf, n.strata=ncol(relAHYf),
                  
                  stratum=dat$stratno, 
                  yr=as.numeric(as.factor(dat$year)),  # fecundity
                  
                  site=dat$site, maxsites=max(dat$site, na.rm = TRUE), 
                  Lat=scale(dat$lat)[,], Long=scale(dat$long)[,],
                  HY_bands=dat$HY_bands, bands=dat$bands,
                  nrows=nrow(dat))




# Initial values: 
inits <- function(){list()}

# Parameters monitored
parameters <- c("anns", "f", "sm", "hr", "b0", "b1", "alpha.r", "lin.alpha.s", "alpha0",
                "rd", "r",  "annh", "beta", "annsmean", "N", "Nasy", "Pjuv.strat", "f_init", "Nfl",
                "et.f","et.hr","ets.f","ets.hr","ICC", "var.tt","var.ts","tau.tt", "tau.ts")

# MCMC settings
ni <- 10000
nt <- 1000
nb <- 500
nc <- 3

# Call JAGS from R 
set.seed(runif(1,1,900))  # sometimes this has to be reset


nopi_ipm <- jags(jags.data, inits, parameters,"ipm_ppr.jags",n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

nopi_ipm$Rhat 
traceplot(nopi_ipm, parameters=c("b1"))
traceplot(nopi_ipm, parameters=c("f"))
traceplot(nopi_ipm, parameters=c("anns")) # a lot of times derived quantities say they haven't converged when they have,
# so I never trust traceplots/R-hats for them and try to just look at ACTUAL parameters

save(nopi_ipm, file="nopi_ipm_density_04062024")


plotRhats(nopi_ipm)



################# SOME PLOTS OF SAVED PARAMETERS ##########################
library(RColorBrewer)
nstrata <- 22
n <- nstrata
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
qu <- function(x) quantile(x, c(0.025, 0.975), na.rm=TRUE)
ptc <- 1.2
op <- par(las=1, mfrow=c(2, 1), mar=c(2,4,2,1))
years <- 20



# PLOT 1:  Annual survival "anns"
plot(x=seq(1:years), y=nopi_ipm$mean$anns[1,], type='b', pch=16, ylab='Annual Mean Survival',
     ylim=c(0, 1), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector[1], cex=ptc)
segments(seq(1:years), nopi_ipm$q2.5$anns[1,], seq(1:years), nopi_ipm$q97.5$anns[1,], col=col_vector[1])
for (s in 1:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$anns[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(seq(1:years), nopi_ipm$q2.5$anns[s,], seq(1:years), nopi_ipm$q97.5$anns[s,], col=col_vector[s])
}
axis(2)
axis(1, at=1:30)


# PLOT 2:  Monthly? survival "sm"
plot(x=seq(1:years), y=nopi_ipm$mean$sm[1,], type='b', pch=16, ylab='sm',
     ylim=c(0, 1), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector[1], cex=ptc)
segments(seq(1:years), nopi_ipm$q2.5$sm[1,], seq(1:years), nopi_ipm$q97.5$sm[1,], col=col_vector[1])
for (s in 1:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$sm[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(seq(1:years), nopi_ipm$q2.5$sm[s,], seq(1:years), nopi_ipm$q97.5$sm[s,], col=col_vector[s])
}
axis(2)
axis(1, at=1:30)


# PLOT 3: fecundity by strata 
plot(x=seq(1:years), y=nopi_ipm$mean$Nfl[1,], type='b', pch=16, ylab='f',
     ylim=c(0, 5), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector[1], cex=ptc)
segments(seq(1:years), nopi_ipm$q2.5$Nfl[1,], seq(1:years), nopi_ipm$q97.5$Nfl[1,], col=col_vector[1])
for (s in 1:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$Nfl[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(seq(1:years), nopi_ipm$q2.5$Nfl[s,], seq(1:years), nopi_ipm$q97.5$Nfl[s,], col=col_vector[s])
}
axis(2)
axis(1, at=1:30)



# PLOT 4:  "rd" direct recovery (first year after release) per strata
plot(x=seq(1:nstrata), y=nopi_ipm$mean$rd[,1], type='b', pch=16, ylab='rd',
     ylim=c(0, 1), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector, cex=ptc)
segments(seq(1:nstrata), nopi_ipm$q2.5$rd[,1], seq(1:nstrata), nopi_ipm$q97.5$rd[,1], col=col_vector)
axis(2)
axis(1, at=1:23)


# PLOT 5:  "r" indirect recovery rates (prob of being recovered and reported) after first year after release per strata
plot(x=seq(1:nstrata), y=nopi_ipm$mean$r[,1], type='b', pch=16, ylab='r',
     ylim=c(0, 1), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector, cex=ptc)
segments(seq(1:nstrata), nopi_ipm$q2.5$r[,1], seq(1:nstrata), nopi_ipm$q97.5$r[,1], col=col_vector)
axis(2)
axis(1, at=1:23)


# PLOT 6:  hazard rate "annh"
plot(x=seq(1:years), y=nopi_ipm$mean$annh[1,], type='b', pch=16, ylab='annh',
     ylim=c(0, 5), xlab=NA, xlim=c(0,years), axes=FALSE, col=col_vector[1], cex=ptc)
segments(seq(1:years), nopi_ipm$q2.5$annh[1,], seq(1:years), nopi_ipm$q97.5$annh[1,], col=col_vector[1])
for (s in 1:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$annh[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
  segments(seq(1:years), nopi_ipm$q2.5$annh[s,], seq(1:years), nopi_ipm$q97.5$annh[s,], col=col_vector[s])
}
axis(2)
axis(1, at=1:30)


# Annual survival by strata "annsmean"
plot(x=1, y=nopi_ipm$mean$annsmean[1], type='b', pch=16, ylab='Annual Mean Survival',
     ylim=c(0, 1), xlab=NA, xlim=c(0, nstrata), axes=FALSE, col=col_vector[1], cex=ptc)
for (s in 1:nstrata){
points(x=s, y=nopi_ipm$mean$annsmean[s], type='b', pch=16, col=col_vector[s], cex=ptc)
segments(s, nopi_ipm$q2.5$annsmean[s], s, nopi_ipm$q97.5$annsmean[s],  col=col_vector[s])
}
axis(2)
axis(1, at=1:30)



# PLOT 1:  Annual all Rhats "anns"
plot(x=seq(1:years), y=nopi_ipm$Rhat$anns[1,], type='b', pch=16, ylab='Rhat Annual Mean Survival',
     ylim=c(0.95, 1.25), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$anns[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)


# PLOT 2:  Annual all Rhats "f"
plot(x=seq(1:years), y=nopi_ipm$Rhat$f[1,], type='b', pch=16, ylab='Rhat Fecundity',
     ylim=c(0.95, 2), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$f[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)


# PLOT 3:  Annual all Rhats "sm"
plot(x=seq(1:years), y=nopi_ipm$Rhat$sm[1,], type='b', pch=16, ylab='Rhat Monthly Survival',
     ylim=c(0.9, 2), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$sm[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)


# PLOT 3:  Annual all Rhats "N"
plot(x=seq(1:years), y=nopi_ipm$Rhat$N[1,], type='b', pch=16, ylab='Rhat N',
     ylim=c(0.95, 1.3), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$N[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)



# PLOT 3:  Annual all Rhats "N"
plot(x=seq(1:years), y=nopi_ipm$mean$N[1,], type='b', pch=16, ylab='Strata Density Estimate',
     ylim=c(0, 7), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$mean$N[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)




# PLOT 3:  Annual all Rhats "N"
plot(x=seq(1:years), y=nopi_ipm$Rhat$N[1,], type='b', pch=16, ylab='Strata Density Estimate',
     ylim=c(0.95, 1.4), xlab="Year", xlim=c(0,years), axes=TRUE, col=col_vector[1], cex=ptc)
for (s in 2:nstrata){
  points(x=seq(1:years), y=nopi_ipm$Rhat$N[s,], type='b', pch=16, col=col_vector[s], cex=ptc)
}
axis(2)
axis(1, at=1:30)



#### CALCULATE CROSS CORRELATIONS

CC_survival <- rcorr(t(nopi_ipm$mean$anns))
CC_survival <- CC_survival$r
TF<- lower.tri(CC_survival, diag = FALSE)
CC_survival <- CC_survival*TF
CC_survival[CC_survival==0] = NA
CC_survival[CC_survival=="NaN"] = NA


CC_fec <-  rcorr(t(nopi_ipm$mean$f))
CC_fec <- CC_fec$r
TF<- lower.tri(CC_fec, diag = FALSE)
CC_fec <- CC_fec*TF
CC_fec[CC_fec==0] = NA
CC_fec[CC_fec=="NaN"] = NA


## calculate distances
library(sf)
library(geosphere)
setwd("C:/Users/ema/Desktop/Waterfowl")
## Define projections for sf objects. Layers are EPSG and WGS84: want to transform to UTM33
lonlatproj <- "+proj=longlat +unit=dd +datum=WGS84"
test <- "EPSG:4326"

BirdData <- read.csv("wbphs_stratum_densities/wbphs_stratum_densities_forDistribution.csv")# Load the bird location data
BirdData[is.na(BirdData)] <- 0 # Include NA as 0, but that is likely not correct

BirdLocs <- read_sf(dsn = "WBPHS_Stratum_Boundaries/WBPHS_Stratum_Boundaries.shp")
BirdLocsLONLAT <- st_transform(BirdLocs, crs=lonlatproj) # Transform Kantons to UTM

StrataLocs <- read.csv("WBPHS_Segment_Counts/wbphs_segment_effort_forDistribution.csv")# Load the bird location data
StrataLocs$midpoint_latitude<- as.numeric(StrataLocs$midpoint_latitude)
StrataLocs$midpoint_longitude<- as.numeric(StrataLocs$midpoint_longitude)
StrataLocs <- na.omit(StrataLocs)
StrataLocs_sf <- st_as_sf(x=StrataLocs, coords=c("midpoint_longitude", "midpoint_latitude"), crs=test) # Make it an sf object
StrataLocs_sf <- st_transform(StrataLocs_sf, crs=lonlatproj) # Transform Kantons to UTM


## Calculate mean survey location coordinates per strata
Coords <- aggregate(st_coordinates(StrataLocs_sf), 
                    by=list(stratumcoords=data.frame(StrataLocs_sf)[,"stratum"]), FUN=mean)
Coordssf <- st_as_sf(x=Coords, coords=c("X", "Y"), crs=lonlatproj) # Make it an sf object
distmat <- distm(Coords[,2:3], fun=distGeo)/1000

disttab<-disttab[-23,]
disttab<-disttab[,-23]
>>>>>>>> 5e82e6b5c208670730fa0c0c13965aac3a3f6c17:IPM_ICC_102324.R
