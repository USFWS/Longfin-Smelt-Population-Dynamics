#### 1. INTRODUCTION ####
# 2024-03-14
# Vanessa Tobias <vanessa_tobias@fws.gov>

# The purpose of this code is to document the population dynamics model
#    appendix/technical note, as it appears in the Longfin Smelt Species
#    Status Assessment.

# This code is part of LFS Technical Note 4 (v. 3.2).
# It fits the population model and creates figures 5-10, which are the
#    graphs related to the population model and the population viability
#    analysis related to it.

# The code file called "Data_Exploration_Figures.R" creates figures 2, 3, and 4,
#   which are the data exploration graphs.


#### 2. SETUP ####
library(R2jags)
library(dplyr)
library(RColorBrewer)
library(plotrix)
library(lattice)


#### 3. FUNCTIONS ####
timePlot <- function(yVar, time, ylab = "", main = "", altLim = FALSE){
  # yVar = call to a variable (use $ notation) e.g., bhMod.cov$BUGSoutput$sims.list$a
  # time = vector of years; must be the same length as yVar
  # ylab = text for y axis label
  
  if(altLim){
    yVar <- yVar[, 2:36]
  }
  
  ul = apply(yVar, 2, function(x) quantile(x, 0.975))
  ll = apply(yVar, 2, function(x) quantile(x, 0.025))
  ml = apply(yVar, 2, mean)
  
  plot(time,
       ul,
       type = "n",
       ylab = ylab, 
       xaxt = "n",
       xlab = "",
       ylim = if(altLim) c(0, 500) else c(min(ll), max(ul)),
       xlim = c(1980, 2015),
       main = main)
  polygon(c(time, rev(time)),
          c(ul, rev(ll)),
          border = NA,
          col = rgb(0, 0, 0, 1/4))
  lines(time, ml, lwd = 2)
  axis(side = 1,
       at = 1980:2015, labels = NA)
  axis(side = 1,
       at = seq(1980, 2015, by = 5),
       labels = seq(1980, 2015, by = 5))
}

maxDens <- function(x){
  # use this to find the peak in the posterior distributions
  # mean doesn't really work because some of them are very skewed
  # be careful using this on distributions that aren't unimodal
  dens <- density(x)
  maxY <- which.max(dens$y)
  maxX <- dens$x[maxY]
  return(maxX)
}

indPlot <- function(mod, x, main = x,
                    params, axisLabs = params,
                    labAngle = TRUE,
                    abline = TRUE){
  # makes a plot of inclusion probabilities for
  #  covariates in a model, based on the spike and
  #  slab setup
  m <- apply(eval(parse(text = 
                          paste0(mod, 
                                 "$BUGSoutput$sims.list$",
                                 x))),
             2, 
             mean)
  l <- eval(
    parse(text = paste0("dim(",
                        mod, 
                        "$BUGSoutput$sims.list$",
                        x,
                        ")[[1]]")))
  plotCI(x = m,
         uiw = 1.96*sqrt((m*(1-m))/l),
         pch = 16,
         xaxt = "n",
         xlab = "",
         ylab = "Inclusion Probability",
         main = main)
  if(labAngle){
    axis(side = 1, 
         at = 1:length(params), 
         labels = NA)
    text(x = 1:length(params),
         y = par("usr")[3] - abs(0.2*(par("usr")[3]-par("usr")[4])),
         labels = axisLabs,
         srt = 45,
         adj = 0.9,
         cex = 0.75,
         xpd = TRUE)
  } else{
    axis(side = 1, 
         at = 1:length(params), 
         labels = axisLabs,
         cex.axis = 0.75)
  }
  
  if(abline) abline(h = 0.2, lty = 2, col = "grey")
}

inv.logit <- function(x) {
  # if x > 709  exp(x) returns inf
  # this fix isn't strictly mathematically correct because inf/inf is undefined,
  # but it's rare and this can be updated later
  ifelse(x < 710, exp(x)/(1+exp(x)), 1)
  }


#### 4. DATA ####
#### 4.1 Abundance Data ####

# new SFBS data through 2020 from Jillian Burns
# This file goes through 2017, but I'm not using all of the years (because of NAs)
#   so you could use one of the older files that goes to 2015 here.

sfbs <- read.csv("./Data_Original/LFSIndices1980to2020_BayStudy23Sep2021.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE)

#### 4.2 Environmental Covariates ####

yrClass <- read.csv("./Data_Original/yrClass.csv",
                    header = TRUE)
# year starts at 1980
# fix typo in name
names(yrClass)[grep("ngpo", names(yrClass))] <- sub("ngpo", "npgo", names(yrClass)[grep("ngpo", names(yrClass))])

# Data compiled by Joe Miller and formatted by Bryan Matthias:
# Temperature and salinity
tempSal <- read.csv("./Data_Original/Temperature_and_Salinity_Summary.csv")
tempSal2 <- tempSal[, c("WY",
                        "MRZ_days_under_15deg",
                        "MRZ_cond_Spawn_season_mean",
                        "MRZ_cond_mean",
                        "MRZ_cond_Jul_mean",
                        "MRZ_cond_Oct_mean",
                        "MRZ_temp_Spawn_season_mean",
                        "MRZ_temp_mean",
                        "MRZ_temp_Jul_mean",
                        "MRZ_temp_Oct_mean")]

yrClass <- merge(yrClass, tempSal2,
                 by.x = "brood.year",
                 by.y = "WY",
                 all.x = TRUE)

# Dayflow
flow <- read.csv("./Data_Original/Dayflow_Data_Summary.csv")
flow <- flow[, c("WY", "OUT_median", "TOT_median")]
# OUT = North Delta Outflow Index
# TOT = Total Delta Inflow


# Food
# data from zooper package (Bashevkin), compiled by Joe Miller for his DFP project
# station NZ048 is in Suisun Bay, in the channel
# We've been using Martinez for the water quality measurements. This was the
# closest zooplankton station that we had data for
mysids <- read.csv("./Data_Original/Mysida_Summary.csv")
# mysids are all mean CPUE
eury <- read.csv("./Data_Original/Eurytemora_affinis_Summary.csv")
# eury are all mean_Eurytemora_affinis_Carbon_micrograms_per_cubic_meter

names(mysids)
mysids <- mysids[, c("WY", "NZ048_Nov_Mysida_CPUE")]
mysids$year <- mysids$WY #+ 1 #because November data
# november and march were good for having few NAs. I chose november so it wasn't
# as likely to be 100% correlated with march eury data
mysids <- mysids[which(mysids$year %in% 1980:2015),]

names(eury)
eury <- eury[,c("WY", "NZ054_Mar_mean_Eurytemora_affinis_Carbon_micrograms_per_cubic_meter")]
names(eury) <- c("WY", "Eury_C_Mar")
eury$year <- eury$WY #because March data
# march because LFS <25mm like them so we need spring data
# february might also have been good, but there were lots of NAs
eury <- eury[which(eury$year %in% 1980:2015),]

covars3 <- merge(yrClass[, c("brood.year", "upwellAug.mean")],
                 flow[, c("WY", "OUT_median", "TOT_median")],
                 by.x = "brood.year", by.y = "WY",
                 all.x = TRUE)
names(covars3)[1] <- "year"
covars3 <- merge(covars3,
                 tempSal[, c("WY", 
                             "MRZ_days_under_15deg", 
                             "MRZ_days_over_20deg",
                             "MRZ_days_under_2PSU")],
                 by.x = "year", by.y = "WY",
                 all.x = TRUE)
covars3 <- merge(covars3,
                 mysids[, c("year", "NZ048_Nov_Mysida_CPUE")],
                 all.x = TRUE)
covars3 <- merge(covars3,
                 eury[, c("year", "Eury_C_Mar")],
                 all.x = TRUE)

# get the correct years
covars3 <- covars3[which(covars3$year %in% 1980:2015),]

# drop the year column
covars4 <- covars3[, -c(1)]

# get a list of names
covarList <- names(covars4)

# There is one NA in the eurytemora column. Need to fix that.
# scale and center covars
covars4 <- apply(covars4, 2,
                 FUN = function(x) scale(x))
# replace missing values with means
# since these are centered, mean = 0 for all columns
covars4[which(is.na(covars4))] <- 0


#### 4.3 Package Data for JAGS ####

bhDat <- list(Y = as.matrix(log(sfbs[1:36, c(     "mwt.0",
                                                  "mwt.1",
                                                  "mwt.2",
                                                  "ot.0",
                                                  "ot.1",
                                                  "ot.2")]+0.01)),
              #covars = as.matrix(covars),
              covarsa = as.matrix(covars4),
              covarsb = as.matrix(covars4),
              covarsS0 = as.matrix(covars4),
              covarsS1 = as.matrix(covars4),
              nYears = 36,
              #p = ncol(covars),
              pa = length(covarList),
              pb = length(covarList),
              pS0 = length(covarList),
              pS1 = length(covarList)
)
covarLista <- covarList
covarListb <- covarList
covarListS0 <- covarList
covarListS1 <- covarList



### 4.6 Data Exploration ####
png("./Figures/splom.png",
    height = 10, width = 10, units = "in", res = 300)
pairs(covars4, 
      pch = 16, col = rgb(0, 0, 0, 1/4),
      cex.axis = 0.6)
dev.off()

png("./Figures/ab_cov.png",
    height = 11, width = 11, units = "in", res = 300)
par(mfrow = c(9, 6),
    mar= c(4.1, 4.1, 1.1, 0.1))
for(j in dimnames(bhDat$covarsa)[[2]]){
  for(i in dimnames(bhDat$Y)[[2]]){
    plot(bhDat$covarsa[,j],
         bhDat$Y[,i],
         pch = 21,
         bg = "grey",
         xlim = c(-3, 4),
         ylim = c(2, 13),
         xlab = j,
         ylab = i)
  }
}
par(mfrow = c(1, 1),
    mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()
# ylim removes low outliers in mwt.2 & ot.2 (one each) from graphs

# Are days over 10K cfs correlated with mysids?
plot(covars4[,1], covars4[,2])
cor(covars4[,1], covars4[,2])
# 0.287
# sort of, but probably not enough to be concerned about using them both

plot(yrClass$Daysover10k.0, yrClass$totalMys)
plot(log(yrClass$Daysover10k.0), log(yrClass$totalMys))
cor(x = log(yrClass$Daysover10k.0), 
    y = log(yrClass$totalMys),
    use = "pairwise.complete.obs")
# still not super high: 0.3170838

#### 5. MODELS ################################################

#### 5.3 "bh_cov" with logit link for survival ####
# with SSVS for covariates
# THIS IS THE MODEL THAT I USED IN THE TN FOR VARIABLE SELECTION

cat("
    model {

    # PRIORS

    logN.est[1, 1] ~ dnorm(9, 0.01)      # Log initial population size
    logN.est[1, 2] ~ dnorm(8, 0.01)
    logN.est[1, 3] ~ dnorm(2, 0.01)

    mean.a[1] ~ dnorm(0, 0.0001)T(0,) #dunif(0.01, 100000) # 
    mean.b[1] ~ dnorm(0, 0.0001)T(0,) #dunif(0.001, 5)

    sigma.a ~ dunif(0, 10)            
    tau.a <- pow(sigma.a, -2)
    sigma.b ~ dunif(0, 10)            
    tau.b <- pow(sigma.b, -2)

    sigma.S0 ~ dunif(0, 10)            
    tau.S0 <- pow(sigma.S0, -2)
    sigma.S1 ~ dunif(0, 10)            
    tau.S1 <- pow(sigma.S1, -2)

    for(i in 1:6){
    sigma.obs[i] ~ dunif(0, 10)          # SD of observation processes
    tau.obs[i] <- pow(sigma.obs[i], -2)
    }

    for (i in 1:3){
    k[i] ~ dnorm(0, 0.0001)
    }

    # Age-0 a parameter covariates
    int.a ~ dunif(0, 1000) #dnorm(0, 0.001)
    int.b ~ dunif(0, 10) # dnorm(0, 0.001)
    # dbeta(5, 5) = bell-shaped, centered at 0.5
    int.S0 ~ dbeta(5, 5) # dnorm(0, 0.001)
    int.S1 ~ dbeta(5, 5) # dnorm(0, 0.001)

    # LIKELIHOOD
    # -- STATE PROCESS --
# from DSM LCM for estimating reprod. directly
#nPL[i] given nA[i-1]
# logmu.R[i] <- zeta[1]+inprod(X_rho[i,],zeta[2:(1+num_unique_cov_rho)])
# rho[i] ~ dlnorm(logmu.R[i],1/(pn.sig.R^2))
# nPL[i] <- nA[i-1]*rho[i]

    for(t in 2:(nYears)){
# inference is on a and b - used to calculate age-0 abundance
    logN.est[t, 1] <- log(exp(logN.est[t, 3])/(a_prime[t]+b_prime[t]*exp(logN.est[t, 3])))

    a_prime[t] <- 1/a[t]
    b_prime[t] <- b[t]/a[t]

    a[t] ~ dnorm(mean.a[t], tau.a)T(0,)
    b[t] ~ dnorm(mean.b[t], tau.b)T(0,)

    mean.a[t] <- int.a + inprod(covarsa[t,], beta.a)
    
    mean.b[t] <- int.b + inprod(covarsb[t,], beta.b)
    }
    
    # variable selection for linear models
    for(j in 1:pa){ #p = number of variables in covars matrix
    ind.a[j] ~ dbern(pind.a)
    betaT.a[j] ~ dnorm(0, taubeta.a)
    beta.a[j] <- ind.a[j]*betaT.a[j]
    }
    for(j in 1:pb){   
    ind.b[j] ~ dbern(pind.b)
    betaT.b[j] ~ dnorm(0, taubeta.b)
    beta.b[j] <- ind.b[j]*betaT.b[j]
    }
    taubeta.a ~ dgamma(1, 0.001)
    pind.a ~ dbeta(2, 8) #dunif(0, 1) # peak at 20% of vars included with long right tail
    taubeta.b ~ dgamma(1, 0.001)
    pind.b ~ dbeta(2, 8) #dunif(0, 1) # peak at 20% of vars included with long right tail


    for (t in 1:(nYears-1)){

    logit.S0.mean[t] <- int.S0 + inprod(covarsS0[t,], beta.S0)
    logit.S1.mean[t] <- int.S1 + inprod(covarsS1[t,], beta.S1)

    S0.logit[t] ~ dnorm(logit.S0.mean[t], tau.S0)
    S1.logit[t] ~ dnorm(logit.S1.mean[t], tau.S1)

    S0[t] <- ilogit(S0.logit[t])
    S1[t] <- ilogit(S1.logit[t])

    logN.est[t+1, 2] <- S0[t]*logN.est[t, 1]
    logN.est[t+1, 3] <- S1[t]*logN.est[t, 2]
    
    }
    # variable selection for linear models
    for(j in 1:pS0){ #p = number of variables in covars matrix
    ind.S0[j] ~ dbern(pind.S0)
    betaT.S0[j] ~ dnorm(0, taubeta.S0)
    beta.S0[j] <- ind.S0[j]*betaT.S0[j]
    }
    for(j in 1:pS1){
    ind.S1[j] ~ dbern(pind.S1)
    betaT.S1[j] ~ dnorm(0, taubeta.S1)
    beta.S1[j] <- ind.S1[j]*betaT.S1[j]
    }
    taubeta.S0 ~ dgamma(1, 0.001)
    pind.S0 ~ dbeta(2, 8) #dunif(0, 1) # peak at 20% of vars included with long right tail
    taubeta.S1 ~ dgamma(1, 0.001)
    pind.S1 ~ dbeta(2, 8) #dunif(0, 1) # peak at 20% of vars included with long right tail


    # -- OBSERVATION PROCESS --
    for (t in 1:nYears){
    Y[t, 1] ~ dnorm(logN.est[t, 1] + k[1], tau.obs[1])
    Y[t, 2] ~ dnorm(logN.est[t, 2] + k[2], tau.obs[2])
    Y[t, 3] ~ dnorm(logN.est[t, 3] + k[3], tau.obs[3])

    Y[t, 4] ~ dnorm(logN.est[t, 1], tau.obs[4])
    Y[t, 5] ~ dnorm(logN.est[t, 2], tau.obs[5])
    Y[t, 6] ~ dnorm(logN.est[t, 3], tau.obs[6])
    }

    # Population sizes on real scale
    for (t in 1:nYears){
    N.est[t, 1] <- exp(logN.est[t, 1])
    N.est[t, 2] <- exp(logN.est[t, 2])
    N.est[t, 3] <- exp(logN.est[t, 3])

    # calculate recruits per spawner
    #logrps[t] <- logN.est[t, 1]/logN.est[t, 3]
    rps[t] <- N.est[t, 1]/N.est[t, 3]
              #exp(logrps[t])
    }
    }
    ", file = "./bh_cov.txt")




#### run "bh_cov ####
bhParams.cov <- c("a", "b",
              "mean.a", "mean.b",
              "int.a", "int.b",
              "int.S0", "int.S1",
              "S0", "S1",
              "S0.logit", "S1.logit",
              "logit.S0.mean", "logit.S1.mean",
              #"mean.S0", "mean.S1",
              "N.est", "logN.est",
              "beta.a", "beta.b",
              "rps",
              "beta.S0", "beta.S1",
              "ind.a", "ind.b",
              "ind.S0", "ind.S1",
              "k", 
              "tau.a", "tau.b",
              "tau.S0", "tau.S1")

bhInits <- function(){list(
  logN.est = matrix(c(rnorm(1, 9, 0.01),
                      rnorm(1, 8, 0.01),
                      rnorm(1, 2, 0.01),
                      rep(NA, 3*(36-1))),
                    ncol = 3, byrow = TRUE))}


# This model takes a little over an hour to run on my laptop.
timestamp()
bhMod.cov <- jags.parallel(model.file = "./bh_cov.txt",
                  data = bhDat,
                  inits = bhInits,
                  parameters.to.save = bhParams.cov,
                  n.chains = 3,
                  n.iter = 1000000,
                  n.burnin = 10000,
                  n.thin = 100,
                  DIC = TRUE)
timestamp()

#### Graph parameter densities ####
# png("./Figures/IEP Poster/paramsCov_dens.png",
png("./Figures/paramsCov_dens.png",
    height = 4, width = 6, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(2.1, 4.1, 3.1, 2.1))
plot(density(bhMod.cov$BUGSoutput$sims.list$a), 
     main = "Intrinsic Productivity")#,
     #xlim = c(0, 200))
plot(density(bhMod.cov$BUGSoutput$sims.list$b), main = "Density Dependence")#,
#xlim = c(0, 5))
# plot(density(bhMod.cov$BUGSoutput$sims.list$rps), main = "Recruits per Spawner",
#      xlim = c(0, 1000))

plot(density(bhMod.cov$BUGSoutput$sims.list$S0), main = "Age-0 Survival")#,
#xlim = c(0, 1))
plot(density(bhMod.cov$BUGSoutput$sims.list$S1), main = "Age-1 Survival")#,
#xlim = c(0, 1))
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

plot(density(bhMod.cov$BUGSoutput$sims.list$int.S0), main = "SO Intercept")
plot(density(bhMod.cov$BUGSoutput$sims.list$int.S1), main = "S1 Intercept")

#### Graph params over time ####
# png("./Figures/IEP Poster/params_time_cov_all.png",
png("./Figures/params_time_cov_all.png",
    height = 10, width = 7.5, units = "in", res = 300)
par(mfcol = c(5, 1), mar = c(2.1, 4.1, 1.5, 1.1),
    cex = 1)
timePlot(yVar = bhMod.cov$BUGSoutput$sims.list$a, 
         time = 1981:2015, 
         ylab = "a",
         main = "Intrinsic Productivity")
timePlot(yVar = bhMod.cov$BUGSoutput$sims.list$b, 
         time = 1981:2015, 
         ylab = "b",
         "Density Dependence")
timePlot(yVar = bhMod.cov$BUGSoutput$sims.list$rps, #c(NA, bhMod.cov$BUGSoutput$mean$rps[2:36]), 
         time = 1981:2015, 
         ylab = "RPS",
         "Recruits per Spawner",
         altLim = TRUE)
timePlot(yVar = bhMod.cov$BUGSoutput$sims.list$S0, 
         time = 1981:2015,
         ylab = "S0",
         "Age-0 Survival")
timePlot(yVar = bhMod.cov$BUGSoutput$sims.list$S1, 
         time = 1981:2015, 
         ylab = "S1",
         "Age-1 Survival")
par(mfcol = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1),
    cex = 1)
dev.off()

# Table for variable selection metrics ####

write.csv(
  data.frame(time = 1981:2015,
             a = bhMod.cov$BUGSoutput$mean$a,
             b = bhMod.cov$BUGSoutput$mean$b,
             rps = bhMod.cov$BUGSoutput$mean$rps[1:35],
             S0 = bhMod.cov$BUGSoutput$mean$S0,
             S1 = bhMod.cov$BUGSoutput$mean$S1),
  "./Data_Derived/bhMod_paramMeans.csv", row.names = FALSE)

data.frame(
  param = c(
    rep("a", length(bhMod.cov$BUGSoutput$mean$ind.a)),
    rep("b", length(bhMod.cov$BUGSoutput$mean$ind.b)),
    rep("S0", length(bhMod.cov$BUGSoutput$mean$ind.S0)),
    rep("S1", length(bhMod.cov$BUGSoutput$mean$ind.S1))
  ),
  vars = c( #!!!BRITTLE!!!
    covarLista,
    covarListb,
    covarListS0,
    covarListS1
  ),
  ind = c(
    apply(bhMod.cov$BUGSoutput$sims.list$ind.a,
          2, 
          mean),
    apply(bhMod.cov$BUGSoutput$sims.list$ind.b,
          2, 
          mean),
    apply(bhMod.cov$BUGSoutput$sims.list$ind.S0,
          2, 
          mean),
    apply(bhMod.cov$BUGSoutput$sims.list$ind.S1,
          2, 
          mean)),
  beta = c(
    apply(bhMod.cov$BUGSoutput$sims.list$beta.a, 2,
          FUN = function(x) maxDens(x)),
    apply(bhMod.cov$BUGSoutput$sims.list$beta.b, 2,
          FUN = function(x) maxDens(x)),
    apply(bhMod.cov$BUGSoutput$sims.list$beta.S0, 2,
          FUN = function(x) maxDens(x)),
    apply(bhMod.cov$BUGSoutput$sims.list$beta.S1, 2,
          FUN = function(x) maxDens(x))
  ))

#### *** Fig. 5 *** ####
# Graph Variable Inclusion probabilities ####
# inclusion probabilities with simplified variable names

png("./Figures/Fig5_paramsCov_ind.png",
    height = 6, width = 8.5, units = "in", res = 600)
par(mfrow = c(2, 2), mar = c(5.1, 4.1, 1.1, 1.1),
    cex = 1)
# covarList = vars in model bhMod.cov
indPlot("ind.a", 
        mod = "bhMod.cov", 
        params = covarLista,
        axisLabs = c("Upwelling", "Outflow", "Inflow", "Low Temp", "High Temp", 
                     "Salinity", "Mysids", "Copepods"),
        labAngle = TRUE,
        abline = FALSE,
        main = "Intrinsic Productivity")
indPlot("ind.b", 
        mod = "bhMod.cov", 
        params = covarListb,
        axisLabs = c("Upwelling", "Outflow", "Inflow", "Low Temp", "High Temp", 
                     "Salinity", "Mysids", "Copepods"),
        labAngle = TRUE,
        abline = FALSE,
        main = "Density Dependence")
indPlot("ind.S0", 
        mod = "bhMod.cov", 
        params = covarListS0,
        axisLabs = c("Upwelling", "Outflow", "Inflow", "Low Temp", "High Temp", 
                     "Salinity", "Mysids", "Copepods"),
        labAngle = TRUE,
        abline = FALSE,
        main = "Age-0 Survival")
indPlot("ind.S1", 
        mod = "bhMod.cov", 
        params = covarListS1,
        axisLabs = c("Upwelling", "Outflow", "Inflow", "Low Temp", "High Temp", 
                     "Salinity", "Mysids", "Copepods"),
        labAngle = TRUE,
        abline = FALSE,
        main = "Age-1 Survival")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()



#### Graph States ####
# make a graph of states with truncated y axis and *s for high points
# png("./Figures/States_covall.png",
#     height = 4, width = 6, units = "in", res = 300)
plot(bhMod.cov$BUGSoutput$mean$N.est[, 1],
     pch = 16, col = brewer.pal(3, "Dark2")[1],
     ylim = c(0, 35000),
     xlab = "Year",
     ylab = "Abundance",
     xaxt = "n")
points(bhMod.cov$BUGSoutput$mean$N.est[, 2],
       pch = 16, col = brewer.pal(3, "Dark2")[2])
points(bhMod.cov$BUGSoutput$mean$N.est[, 3],
       pch = 16, col = brewer.pal(3, "Dark2")[3])
axis(side = 1, at = 1:36, labels = 1980:2015)
points(c(which(bhMod.cov$BUGSoutput$mean$N.est[, 1]>35000), 
         which(bhMod.cov$BUGSoutput$mean$N.est[, 2]>35000)),
       rep(35000, length(c(which(bhMod.cov$BUGSoutput$mean$N.est[, 1]>35000), 
                           which(bhMod.cov$BUGSoutput$mean$N.est[, 2]>35000)))), 
       pch = "*", cex = 1.25,
       col = c(rep(brewer.pal(3, "Dark2")[1], 
                   length(which(bhMod.cov$BUGSoutput$mean$N.est[, 1]>35000))),
               rep(brewer.pal(3, "Dark2")[2],
                   length(which(bhMod.cov$BUGSoutput$mean$N.est[, 2]>35000)))))
legend("topright", 
       pch = 16, col = brewer.pal(3, "Dark2"),
       legend = c("Age-0", "Age-1", "Age-2"),
       bty = "n")
text(x = 31, y = 20000, pos = 4, offset = 0,
     labels = "* = values \n > 35000", cex = 0.85)
# dev.off()




#### 5.4 "bh_cov2" - selected vars ####

covarList2a <- covarList[c(1, 2, 4:8)] 
covarList2b <- covarList[c(2, 4, 6, 7)]
covarList2S0 <- covarList[c(4, 5, 7, 8)] 
covarList2S1 <- covarList[c(1, 4:8)] 

bhDat2 <- list(Y = as.matrix(log(sfbs[1:36, c(     "mwt.0",
                                                   "mwt.1",
                                                   "mwt.2",
                                                   "ot.0",
                                                   "ot.1",
                                                   "ot.2")]+0.01)),
               covarsa = as.matrix(covars4[, c(1, 2, 4:8)]),
               covarsb = as.matrix(covars4[, c(2, 4, 6, 7)]), 
               covarsS0 = as.matrix(covars4[,c(4, 5, 7, 8)]), 
               covarsS1 = as.matrix(covars4[,c(1, 4:8)]),
               nYears = 36,
               pa = length(covarList2a),
               pb = length(covarList2b),
               pS0 = length(covarList2S0),
               pS1 = length(covarList2S1)
)

# Note: this model took about 2 hours to run on my laptop.
timestamp()
bhMod.cov2 <- jags.parallel(model.file = "./bh_cov.txt",
                   data = bhDat2,
                   inits= bhInits,
                   parameters.to.save= bhParams.cov,
                   n.chains=3,
                   n.iter=1000000,
                   n.burnin=10000,
                   n.thin=100,
                   DIC=TRUE)
timestamp()

# Graph parameter densities ####
par(mfrow = c(2, 2), mar = c(2.1, 4.1, 3.1, 2.1))
plot(density(bhMod.cov2$BUGSoutput$sims.list$a), 
     main = "Intrinsic Productivity")#,
#xlim = c(0, 200))
plot(density(bhMod.cov2$BUGSoutput$sims.list$b), main = "Density Dependence")#,
#xlim = c(0, 5))
# plot(density(bhMod.cov$BUGSoutput$sims.list$rps), main = "Recruits per Spawner",
#      xlim = c(0, 1000))

plot(density(bhMod.cov2$BUGSoutput$sims.list$S0), main = "Age-0 Survival")#,
#xlim = c(0, 1))
plot(density(bhMod.cov2$BUGSoutput$sims.list$S1), main = "Age-1 Survival")#,
#xlim = c(0, 1))
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))


# Graph Variable Inclusion probabilities ####
# png("./Figures/paramsCov_ind.png",
#     height = 6, width = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(7.1, 4.1, 1.1, 1.1),
    cex = 1)
indPlot("ind.a", mod = "bhMod.cov2", params = covarList2a)
indPlot("ind.b", mod = "bhMod.cov2", params = covarList2b)
indPlot("ind.S0", mod = "bhMod.cov2", params = covarList2S0)
indPlot("ind.S1", mod = "bhMod.cov2", params = covarList2S1)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()



#### *** Fig. 6 *** ####
#### Graph params over time ####

png("./Figures/Fig6_params_time_cov2.png",
    height = 8, width = 6, units = "in", res = 300)
par(mfrow = c(5, 1), mar = c(2.1, 4.1, 1.5, 1.1),
    cex = 1)
timePlot(yVar = bhMod.cov2$BUGSoutput$sims.list$a, 
         time = 1981:2015, 
         ylab = "a",
         main = "Intrinsic Productivity")
timePlot(yVar = bhMod.cov2$BUGSoutput$sims.list$b, 
         time = 1981:2015, 
         ylab = "b",
         main = "Density Dependence")
timePlot(yVar = bhMod.cov2$BUGSoutput$sims.list$rps[,2:36], 
         time = 1981:2015, 
         ylab = "RPS",
         "Recruits per Spawner")
timePlot(yVar = bhMod.cov2$BUGSoutput$sims.list$S0, 
         time = 1981:2015, 
         ylab = "S0",
         main = "Age-0 Survival")
timePlot(yVar = bhMod.cov2$BUGSoutput$sims.list$S1, 
         time = 1981:2015, 
         ylab = "S1",
         main = "Age-1 Survival")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1),
    cex = 1)
dev.off()

#### *** Fig. 7 *** ####
#### Graph States ####
# make a graph of states with truncated y axis and *s for high points
png("./Figures/Fig7_States_cov2.png",
    height = 4, width = 6, units = "in", res = 300)
plot(bhMod.cov2$BUGSoutput$mean$N.est[, 1],
     #pch = 16, 
     type = "l", lwd = 3,
     col = brewer.pal(3, "Dark2")[1],
     #ylim = c(0, 35000),
     xlab = "Year",
     ylab = "Abundance",
     xaxt = "n")
points(bhMod.cov2$BUGSoutput$mean$N.est[, 2],
       #pch = 16, 
       type = "l", lwd = 3, 
       col = brewer.pal(3, "Dark2")[2])
points(bhMod.cov2$BUGSoutput$mean$N.est[, 3],
       #pch = 16, 
       type = "l", lwd = 3,
       col = brewer.pal(3, "Dark2")[3])
axis(side = 1, at = 1:36, labels = 1980:2015)
legend("topright", 
       #pch = 16, 
       lty = 1, lwd = 3,
       col = brewer.pal(3, "Dark2"),
       legend = c("Age-0", "Age-1", "Age-2"),
       bty = "n")
dev.off()

plot(bhMod.cov2$BUGSoutput$summary[1:35, 7],
     type = "n",
     ylim = c(0, 250000),
     ylab = "Estmated Abundance",
     xlab = "Year",
     xaxt = "n")
# Age 0
polygon(c(1:36, 36:1),
        c(bhMod.cov2$BUGSoutput$summary[1:36, 7], 
          rev(bhMod.cov2$BUGSoutput$summary[1:36, 3])),
        col = adjustcolor("blue", alpha.f = 0.5),
        border = NA)
lines(1:36,
      bhMod.cov2$BUGSoutput$summary[1:36, 1],
      col = "blue", lwd = 3)
# Age 1
polygon(c(1:36, 36:1),
        c(bhMod.cov2$BUGSoutput$summary[37:72, 7], 
          rev(bhMod.cov2$BUGSoutput$summary[37:72, 3])),
        col = adjustcolor("red", alpha.f = 0.5),
        border = NA)
lines(1:36,
      bhMod.cov2$BUGSoutput$summary[37:72, 1],
      col = "red", lwd = 3)
# Age 2+
polygon(c(1:36, 36:1),
        c(bhMod.cov2$BUGSoutput$summary[73:108, 7], 
          rev(bhMod.cov2$BUGSoutput$summary[73:108, 3])),
        col = adjustcolor("green", alpha.f = 0.5),
        border = NA)
lines(1:36,
      bhMod.cov2$BUGSoutput$summary[73:108, 1],
      col = "green", lwd = 3)
# x axis
axis(side = 1, at = 1:36, 1980:2015)


plot(bhMod.cov2$BUGSoutput$summary[71:71+35, 7],
     type = "n")


#### 6. PVA ####
#### 6.1 Projection Matrix Functions ####
# 1981:2070 is 90 years
# last 10 years of the dataset is 2005:2015 = 25:35
# 2006:2070 is 65 years

#### RUN PVA with selected covars ####
GraphMod <- bhMod.cov2

apply(GraphMod$BUGSoutput$sims.list$mean.a[,27:36],
      2,
      FUN = function(x) sample(x, 5))

summary(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,1]))
summary(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,2]))
summary(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,3]))

# set Quasi-Extinction thresholds based on population estimates
# -- it's a little circular, but we're only using the lowest values
QE <- data.frame("QE0" = quantile(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,1]), 0.1),
                 "QE1" = quantile(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,2]), 0.1),
                 "QE2" = quantile(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,3]), 0.1))

projectMat <- function(years.in, years.out, year.start, 
                       QE, GraphMod){
  # samples the posteriors from the state space model
  # projects them into the future
  # returns the vector of projected total population size = sum(ages)
  # years.in = sequence from:to of index values for the years that you want
  #             to use to inform the projection (i.e., select the past years you want to use)
  # years.out = a number of years to project into the future (a single value)
  # year.start = an index value indicating the starting year for the simulation
  #              gives values for N0 for each age in the matrix. 
  #              Note: the user is not directly supplying a starting value for 
  #              the population.
  proj <- data.frame("a" =  sample(GraphMod$BUGSoutput$sims.list$a[,years.in], years.out),
                     "b" =  sample(GraphMod$BUGSoutput$sims.list$b[,years.in], years.out),
                     "S0" = sample(GraphMod$BUGSoutput$sims.list$S0[,years.in], years.out),
                     "S1" = sample(GraphMod$BUGSoutput$sims.list$S1[,years.in], years.out))
  N0 <- data.frame("N00" =  sample(GraphMod$BUGSoutput$sims.list$N.est[ , year.start, 1], 1),
                   "N01" =  sample(GraphMod$BUGSoutput$sims.list$N.est[ , year.start, 2], 1),
                   "N02" =  sample(GraphMod$BUGSoutput$sims.list$N.est[ , year.start, 3], 1))
  Nt <- data.frame("N0t" = NA,
                   "N1t" = NA,
                   "N2t" = NA)
  
  Nt[1, 1] <- (proj$a[1] * N0$N02)/(1 + proj$b[1] * N0$N02)
  Nt[1, 2] <- N0$N00*proj$S0[1]
  Nt[1, 3] <- N0$N01*proj$S1[1]
  # for(i in 1:nrow(proj)){
  #   Nt[i+1, 1] <- (proj$a[i] * Nt[i, 3])/(1 + proj$b[i] * Nt[i, 3])
  #   Nt[i+1, 2] <- Nt[i, 1]*proj$S0[i]
  #   Nt[i+1, 3] <- Nt[i, 2]*proj$S1[i]
  # }
  # apply quasiextinction - if any drop below the threshold they stick
  for(i in 1:nrow(proj)){
    
    N0 <- (proj$a[i] * Nt[i, 3])/(1 + proj$b[i] * Nt[i, 3])
    Nt[i+1, 1] <- if(N0 < QE$QE0) 0 else N0
    
    N1 <- Nt[i, 1]*proj$S0[i]
    Nt[i+1, 2] <- if(N1 < QE$QE1) 0 else N1
    
    N2 <- Nt[i, 2]*proj$S1[i]
    Nt[i+1, 3] <- if(N2 < QE$QE2) 0 else N2
  }
  # calculate total population
  Nt$Nt <- apply(Nt[,1:3], 1, sum)
  return(Nt$Nt)
}


# popMat <- matrix(NA, nrow = 46, ncol = 100)
# for(i in 1:ncol(popMat)){
#   popMat[,i] <- projectMat()
# }

makePopMat <- function(){
  popMat <- matrix(NA, nrow = 46, ncol = 1000)
  for(i in 1:ncol(popMat)){
    popMat[,i] <- projectMat(years.in = 1:35, 
                             years.out = 45,
                             year.start = 36, #25)
                             QE = QE,
                             GraphMod = GraphMod) 
  }
  return(popMat)
}



setQE <- function(prob){
  QE = data.frame("QE0" = quantile(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,1]), probs = prob),
                   "QE1" = quantile(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,2]), probs = prob),
                   "QE2" = quantile(as.vector(GraphMod$BUGSoutput$sims.list$N.est[,,3]), probs = prob))
  return(QE)
}



#### 6.2 PVA Graphs ####
# GraphMod <- bhMod.cov3
QE <- setQE(prob = 0.05)
QE
popMat <- makePopMat()
options(scipen = 10)
png("./Figures/PVA_ab_time.png",
    height = 4, width = 6, units = "in", res = 300)
plot(0, 0, type = "n",
     xlim = c(1, 46),
     ylim = c(0, 100000),
     xlab = "",
     ylab = "Total Estimated Abundance",
     xaxt = "n")
polygon(c(1:46, 46:1),
        c(apply(popMat, 1, 
                FUN = function(x){summary(x)[5]}),
          rev(apply(popMat, 1, 
                    FUN = function(x){summary(x)[2]}))),
        col = "grey", border = NA)
lines(apply(popMat, 1, 
      FUN = function(x){summary(x)[3]}), 
      lty = 1, lwd = 3)
# # 3rd Q
# lines(apply(popMat, 1, 
#             FUN = function(x){summary(x)[5]}),
#       lty = 2, lwd = 2)
# # 1st Q
# lines(apply(popMat, 1, 
#             FUN = function(x){summary(x)[2]}),
#       lty = 2, lwd = 2)
axis(side = 1, at = seq(1, 56, by = 5),
     labels = seq(2015, 2070, by = 5))
dev.off()

png("./Figures/PVA_ab_time95CI.png",
    height = 4, width = 6, units = "in", res = 300)
# main plot that focuses on median line:
par(cex = 1.25, new = FALSE,
    mar = c(2.1, 4.1, 2.1, 2.1),
    fig = c(0, 1, 0, 1)) #default for fig = graph space goes from zero to one
plot(0, 0, type = "n",
     xlim = c(1, 46),
     ylim = c(1, 1164218),
     # log = "y",
     xlab = "",
     ylab = "Total Estimated Abundance",
     xaxt = "n")
# median
lines(apply(popMat, 1, 
            FUN = function(x){quantile(x, probs = 0.5)}),
      lty = 1, lwd = 3)
# mean
# lines(apply(popMat, 1, 
#             FUN = function(x){mean(x)}), 
#       lty = 1, lwd = 3)
# Upper CL
lines(apply(popMat, 1, 
            FUN = function(x){quantile(x, probs = 0.975)}),
      lty = 2, lwd = 1)
lines(apply(popMat, 1, 
            FUN = function(x){quantile(x, probs = 0.95)}),
      lty = 2, lwd = 2)
# Lower CL
lines(apply(popMat, 1, 
            FUN = function(x){quantile(x, probs = 0.025)}),
      lty = 2, lwd = 1)
lines(apply(popMat, 1, 
            FUN = function(x){quantile(x, probs = 0.05)}),
      lty = 2, lwd = 2)
axis(side = 1, at = seq(1, 56, by = 5),
     labels = seq(2015, 2070, by = 5))
legend("topright",
       lwd = c(2, 2, 1),
       lty = c(1, 2, 2),
       legend = c("Median", "90% CI", "95% CI"))
dev.off()

# inset of whole y range:
par(cex = 1.25,
    fig = c(0.5, 1, 0.5, 1), # dock in upper right
    mar = c(2.1, 2.1, 2.1, 2.1), # keep last two the same to align with docking side
    new = TRUE)
plot(0, 0, type = "n",
     xlim = c(1, 46),
     ylim = c(1,
              max(apply(popMat, 1,
                        FUN = function(x){quantile(x, probs = 0.975)}))),
     # log = "y",
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     new = TRUE)
lines(apply(popMat, 1,
            FUN = function(x){mean(x)}), #{quantile(x, probs = 0.5)}),
      lty = 1, lwd = 3)
# Upper CL
lines(apply(popMat, 1,
            FUN = function(x){quantile(x, probs = 0.975)}),
      lty = 2, lwd = 2)
# Lower CL
lines(apply(popMat, 1,
            FUN = function(x){quantile(x, probs = 0.025)}),
      lty = 2, lwd = 2)
axis(side = 1, at = seq(1, 46, by = 5),
     labels = seq(2005, 2050, by = 5))
axis(side = 2, at = c(0,
                     max(apply(popMat, 1,
                               FUN = function(x){
                                 quantile(x, probs = 0.975)}))),
     labels = c(0, round(max(apply(popMat, 1,
                                   FUN = function(x){
                                     quantile(x, probs = 0.975)}))/100000, 0)))
# dev.off()

# *** Fig. 8 *** ####
png("./Figures/Fig8_PVA_ab_time95CI_inset.png",
    height = 4, width = 6, units = "in", res = 300)
# main plot
par(cex = 1.25, new = FALSE,
    mar = c(2.1, 4.1, 2.1, 2.1),
    fig = c(0, 1, 0, 1)) #default for fig = graph space goes from zero to one
plot(0, 0, type = "n",
     xlim = c(1, 46),
     ylim = c(0, 100000),
     xlab = "",
     ylab = "Total Estimated Abundance",
     xaxt = "n")
# grey polygon for middle 50%
polygon(c(1:46, 46:1),
        c(apply(popMat, 1, 
                FUN = function(x){summary(x)[5]}),
          rev(apply(popMat, 1, 
                    FUN = function(x){summary(x)[2]}))),
        col = "grey", border = NA)
lines(apply(popMat, 1, 
            FUN = function(x){summary(x)[3]}), 
      lty = 1, lwd = 3)
axis(side = 1, at = seq(1, 56, by = 5),
     labels = seq(2015, 2070, by = 5))

# inset plot with CIs
par(cex = 1.25,
    fig = c(0.4, 1, 0.5, 1), # dock in upper right
    mar = c(2.1, 2.1, 2.1, 2.1), # keep last two the same to align with docking side
    new = TRUE)
plot(0, 0, type = "n",
     xlim = c(1, 46),
     ylim = c(1,
              max(apply(popMat, 1,
                        FUN = function(x){quantile(x, probs = 0.975)}))),
     # log = "y",
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     new = TRUE)
lines(apply(popMat, 1, # {quantile(x, probs = 0.5)}), #),{mean(x)}
            FUN = function(x){summary(x)[3]}), 
      lty = 1, lwd = 3)
# Upper CL
lines(apply(popMat, 1,
            FUN = function(x){quantile(x, probs = 0.975)}),
      lty = 2, lwd = 2)
# Lower CL
lines(apply(popMat, 1,
            FUN = function(x){quantile(x, probs = 0.025)}),
      lty = 2, lwd = 2)
axis(side = 1, at = seq(1, 56, by = 5),
     labels = seq(2015, 2070, by = 5))
axis(side = 2, 
     at = c(0, 10e5),
     labels = c(0, 10))
mtext(side = 2, expression("Abundance x10"^5), line = 2)
dev.off()


# *** Fig. 9 *** ####
# This figure takes a little while because making the matrices is slow
dev.off() #reset par after insets above
png("./Figures/Fig9_QE_prob_time.png",
    height = 4, width = 6, units = "in", res = 300)
QE <- setQE(prob = 0.15)
popMat <- makePopMat()
par(mfrow = c(1, 1),
    mar = c(5.1, 4.1, 4.1, 2.1))
plot(1:nrow(popMat),
     rep(0, nrow(popMat)),
     type = "n",
     xaxt = "n",
     xlab = "",
     ylab = "Probability of Quasi-Extinction",
     #xlim = c(1, 36), #2050 = foreseeable future
     ylim = c(0, 1),
     main = "")
lines(1:nrow(popMat),
      apply(popMat, 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
      lwd = 2, col = brewer.pal(5, "YlOrRd")[5])
# QE <- setQE(prob = 0.125)
# lines(1:nrow(popMat),
#       apply(makePopMat(), 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
#       lwd = 2, col = brewer.pal(7, "YlOrRd")[6])
QE <- setQE(prob = 0.1)
popMat <- makePopMat()
lines(1:nrow(popMat),
      apply(popMat, 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
      lwd = 2, col = brewer.pal(5, "YlOrRd")[4])
QE <- setQE(prob = 0.05)
popMat <- makePopMat()
lines(1:nrow(popMat),
      apply(popMat, 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
      lwd = 2, col = brewer.pal(5, "YlOrRd")[3])
# QE <- setQE(prob = 0.025)
# lines(1:nrow(popMat),
#       apply(makePopMat(), 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
#       lwd = 2, col = brewer.pal(7, "YlOrRd")[3])
QE <- setQE(prob = 0.01)
popMat <- makePopMat()
lines(1:nrow(popMat),
      apply(popMat, 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
      lwd = 2, col = brewer.pal(5, "YlOrRd")[2])
axis(side = 1, at = seq(1, 56, by = 5),
     labels = seq(2015, 2070, by = 5), 
     cex.axis = 0.75)
legend("topleft", lwd = 2,
       col = rev(brewer.pal(5, "YlOrRd")[2:5]),
       legend = rev(c(1, 5, 10, 15)),
       title = "Quantile for QE Threshold",
       cex = 0.75,
       bty = "n")
dev.off()

# Running the PVA out longer
# -- 100 years instead of 45
# makePopMat2 <- function(){
#   popMat <- matrix(NA, nrow = 101, ncol = 1000)
#   for(i in 1:ncol(popMat)){
#     popMat[,i] <- projectMat(years.in = 1:35, 
#                              years.out = 100,
#                              year.start = 36, #25)
#                              QE = QE,
#                              GraphMod = GraphMod) 
#   }
#   return(popMat)
# }
# # GraphMod = bhMod.cov3
# QE <- setQE(prob = 0.15)
# popMat <- makePopMat2()
# plot(1:nrow(popMat),
#      rep(0, nrow(popMat)),
#      type = "n",
#      xaxt = "n",
#      xlab = "",
#      ylab = "Probability of Quasi-Extinction",
#      #xlim = c(1, 36), #2050 = foreseeable future
#      ylim = c(0, 1),
#      main = "")
# lines(1:nrow(popMat),
#       apply(popMat, 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
#       lwd = 2, col = brewer.pal(5, "YlOrRd")[5])
# QE <- setQE(prob = 0.1)
# popMat <- makePopMat2()
# lines(1:nrow(popMat),
#       apply(popMat, 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
#       lwd = 2, col = brewer.pal(5, "YlOrRd")[4])
# QE <- setQE(prob = 0.05)
# popMat <- makePopMat2()
# lines(1:nrow(popMat),
#       apply(popMat, 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
#       lwd = 2, col = brewer.pal(5, "YlOrRd")[3])
# QE <- setQE(prob = 0.01)
# popMat <- makePopMat2()
# lines(1:nrow(popMat),
#       apply(popMat, 1, FUN = function(x){sum(x == 0)})/ncol(popMat),
#       lwd = 2, col = brewer.pal(5, "YlOrRd")[2])
# axis(side = 1, at = seq(1, 101, by = 5),
#      labels = seq(2015, 2115, by = 5), 
#      cex.axis = 0.75)
# legend("topleft", lwd = 2,
#        col = rev(brewer.pal(5, "YlOrRd")[2:5]),
#        legend = rev(c(1, 5, 10, 15)),
#        title = "Quantile for QE Threshold",
#        cex = 0.75,
#        bty = "n")


QETable <- data.frame("Percentile" = c(0.15, 0.10, 0.05, 0.01),
                      "QE0" = NA,
                      "QE1" = NA,
                      "QE2" = NA)
for(i in 1:length(QETable$Percentile)){
  QETable[i, 2:4] <- round(setQE(QETable[i, 1]), 0)
}

#### 6.X Pop matrix with environmetal vars ####

#### 6.X Assemble Covariates for Projection ####
# make a data frame of projected covariates
# has to have all the covariates that are used for any vital rate

#  in this case c(1, 2, 8, 14)
#  covarList[c(1, 2, 8, 14)]
# For now, use the actual values to set it up
# later, add climate-chnage specific projections based on what's in the SSA
# covarsa = as.matrix(lengths$Age.24[which(lengths$Year %in% 1980:2015)]), #as.matrix(covars4[, covarList3a]),
# covarsb = as.matrix(covars4[, covarList3b]), 
# covarsS0 = as.matrix(covars4[,covarList3S0]), 
# covarsS1 = as.matrix(covars4[,covarList3S1])

#### Historical Covariates ####
covarProj <- covars4 #data.frame(cbind(covars4, lengths$Age.24[which(lengths$Year %in% 1980:2015)]))
# names(covarProj) <- c(names(covars3)[-1], "Age.24")

# Data for bh_cov2
# covarsa = as.matrix(covars4[, c(1, 2, 4:8)]),
# covarsb = as.matrix(covars4[, c(2, 4, 6, 7)]), 
# covarsS0 = as.matrix(covars4[,c(4, 5, 7, 8)]), 
# covarsS1 = as.matrix(covars4[,c(1, 4:8)]),

covarProjNames <- names(covarProj)
# names(covarProj) <- covarList[c(1, 2, 8, 14)]
covarProj <- as.matrix(covarProj)


# scale based on the scaling frame for the original data
# covarProj <- apply(covarProj, 2,
#       FUN = function(x){
#         (x - mean(x[1:36], na.rm = TRUE))/sd(x[1:36], na.rm = TRUE)
#       })
# replace missing values with mean = 0
# covarProj[is.na(covarProj)] <- 0
# length data weren't standardized in the original model...
# covarProj[,9] <- lengths$Age.24[which(lengths$Year %in% 1980:2015)]






#### Function that calculates projections ####
projectMatCov <- function(covarsa,
                          covarsb,
                          covarsS0,
                          covarsS1,
                          n,
                          year.start){
  # 1. make a way to set specific vars for making projections
  # 2. test function
  
  # n <- length(covarsa) #!!!ONLY WORKS IF COVARSA IS A VECTOR
  
  #vars.a <- vars.a
  # below works for a single beta value:
  int.a <- rep(GraphMod$BUGSoutput$mean$int.a, n)
  beta.a <- GraphMod$BUGSoutput$mean$beta.a
  tau.a <- GraphMod$BUGSoutput$mean$tau.a
  
  # mean.a <- int.a + covarsa*beta.a
  mean.a <- if(length(dim(beta.a)) > 1) {
    int.a + covarsa %*% beta.a} else{
      int.a + covarsa * rep(beta.a, n)}
  
  a <- rep(NA, n)
  for(i in 1:n){
    a[i] <- rnorm(1, mean.a[i], sd = sqrt(1/tau.a))
    }
  a[which(a < 0)] <- 0 #zero truncation part
  a_prime <- 1/a
  
 # covarsb <- SET ABOVE
  # 
  int.b <- rep(GraphMod$BUGSoutput$mean$int.b, n)       # vector (n)
  beta.b <- as.matrix(GraphMod$BUGSoutput$mean$beta.b)            # vector (4)
  #dim(beta.b) <- c(36, 4) #turn beta.b into a matrix
  tau.b <- GraphMod$BUGSoutput$mean$tau.b               # scalar
  
  mean.b <- if(length(dim(beta.b)) > 1) {
    int.b + covarsb %*% beta.b} else{
      int.b + covarsS0 * rep(beta.b, n)}
  # mean.b <- as.vector(mean.b)
  b <- rep(NA, n)
  # b <- matrix(NA, nrow = dim(beta.b)[1], ncol = dim(beta.b)[2])
  for(i in 1:n){
    b[i] <- rnorm(1, mean.b[i], sd = sqrt(1/tau.b))
  }
  b[which(b < 0)] <- 0.01 #zero truncation part
  b_prime <- b/a
  
  int.S0 <- rep(GraphMod$BUGSoutput$mean$int.S0, n)
  beta.S0 <- GraphMod$BUGSoutput$mean$beta.S0
  tau.S0  <- GraphMod$BUGSoutput$mean$tau.S0
  
  int.S1 <- rep(GraphMod$BUGSoutput$mean$int.S1, n)
  beta.S1 <- GraphMod$BUGSoutput$mean$beta.S1
  tau.S1  <- GraphMod$BUGSoutput$mean$tau.S1
  
  
  logit.S0.mean <- if(length(dim(beta.S0)) > 1) {
    int.S0 + covarsS0 %*% beta.S0} else{
      int.S0 + covarsS0 * rep(beta.S0, n)}
  
  # logit.S1.mean <- int.S1 + covarsS1 %*% beta.S1
  logit.S1.mean <- if(length(dim(beta.S1)) > 1) {
    int.S1 + covarsS1 %*% beta.S1} else{
      int.S1 + covarsS1 * rep(beta.S1, n)}
  
  S0.logit <- rep(NA, n)
  S1.logit <- rep(NA, n)
  
  for(i in 1:n){
    S0.logit[i] <- rnorm(1, logit.S0.mean[i], sqrt(1/tau.S0))
    S1.logit[i] <- rnorm(1, logit.S1.mean[i], sqrt(1/tau.S1))
  }
  S0 <- inv.logit(S0.logit)
  S1 <- inv.logit(S1.logit)
  
  proj <- data.frame(a,
                     b,
                     S0,
                     S1)
  
  # !!! Every simulation starts at the N0 estimate
  # !!! need to be able to select which year (or value) to start at
  # N0 <- data.frame("N00" = sample(GraphMod$BUGSoutput$sims.list$N.est[,1,1], 1),
  #                  "N01" = sample(GraphMod$BUGSoutput$sims.list$N.est[,1,2], 1),
  #                  "N02" = sample(GraphMod$BUGSoutput$sims.list$N.est[,1,3], 1))
  N0 <- data.frame("N00" =  sample(GraphMod$BUGSoutput$sims.list$N.est[ , year.start, 1], 1),
                   "N01" =  sample(GraphMod$BUGSoutput$sims.list$N.est[ , year.start, 2], 1),
                   "N02" =  sample(GraphMod$BUGSoutput$sims.list$N.est[ , year.start, 3], 1))
  Nt <- data.frame("N0t" = rep(NA, nrow(proj)), # Age-01
                   "N1t" = rep(NA, nrow(proj)), # Age-1
                   "N2t" = rep(NA, nrow(proj))) # Age-2
  
  #logN.est[t, 1] <- log(exp(logN.est[t, 3])/(a_prime[t]+b_prime[t]*exp(logN.est[t, 3])))
  
  Nt[1, 1] <- (proj$a[1] * N0$N02)/(1 + proj$b[1] * N0$N02) #N0$N02/(a_prime[1]+b_prime[1]*N0$N02) 
  Nt[1, 2] <- N0$N00*proj$S0[1]
  Nt[1, 3] <- N0$N01*proj$S1[1]

  # apply quasiextinction - if any drop below the threshold they stick
  for(i in 1:(nrow(proj)-1)){
    
    Nt$N1t[i+1] <- Nt$N0t[i]*proj$S0[i]
    Nt$N1t[i+1] <- if(Nt$N1t[i+1] < QE$QE1) 0 else Nt$N1t[i+1]
    
    # N2 <- Nt[i, 2]*proj$S1[i]
    # Nt[i+1, 3] <- if(N2 < QE$QE2) 0 else N2
    
    Nt$N2t[i+1] <- Nt$N1t[i]*proj$S1[i]
    Nt$N2t[i+1] <- if(Nt$N2t[i+1] < QE$QE2) 0 else Nt$N2t[i+1]

    # Nt$N0t[i+1] <- Nt$N2t[i]/(a_prime[i]+b_prime[i]*Nt$N2t[i])
    Nt$N0t[i+1] <- (proj$a[i] * Nt$N2t[i+1])/(1 + proj$b[i] * Nt$N2t[i+1])
    # Nt$N0t[i+1] <- (proj$a[i] * Nt$N2t[i])/(1 + proj$b[i] * Nt$N2t[i])
    Nt$N0t[i+1] <- if(Nt$N0t[i+1] < QE$QE0) 0 else Nt$N0t[i+1]
    
    # N1 <- Nt[i, 1]*proj$S0[i]
    # Nt[i+1, 2] <- if(N1 < QE$QE1) 0 else N1

  }
  # calculate total population
  Nt$Nt <- apply(Nt[,1:3], 1, sum)
  return(Nt$Nt)
}

#### Function that makes a matrix of projections ####
# GraphMod <- bhMod.cov3
# makePopMatCov <- function(){
#   popMatCov <- matrix(NA, nrow = 36, ncol = 1000)
#   for(i in 1:ncol(popMatCov)){
#     popMatCov[,i] <- projectMatCov(covarsa = covarProj[,which(covarProjNames %in% covarList3a)],
#                                    covarsb = covarProj[,which(covarProjNames %in% covarList3b)],
#                                    covarsS0 =  covarProj[,which(covarProjNames %in% covarList3S0)],
#                                    covarsS1 =  covarProj[,which(covarProjNames %in% covarList3S1)])
#   }
#   return(popMatCov)
# }

makePopMatCov <- function(){
  popMatCov <- matrix(NA, nrow = 36, ncol = 1000)
  for(i in 1:ncol(popMatCov)){
    popMatCov[,i] <- projectMatCov(covarProj,#[,which(covarProjNames %in% covarLista)],
                                   covarProj,#[,which(covarProjNames %in% covarListb)],
                                   covarProj,#[,which(covarProjNames %in% covarListS0)],
                                   covarProj)#[,which(covarProjNames %in% covarListS1)])
  }
  return(popMatCov)
}
# makePopMatCov()


#### RUN THIS: ####
QE <- setQE(prob = 0.10)
popMatCov <- matrix(NA, nrow = 36, ncol = 1000)

# !!! Double check that these are the correct columns ####
for(i in 1:ncol(popMatCov)){
  # covarsa = as.matrix(covars4[, c(1, 2, 4:8)]),
  # covarsb = as.matrix(covars4[, c(2, 4, 6, 7)]), 
  # covarsS0 = as.matrix(covars4[,c(4, 5, 7, 8)]), 
  # covarsS1 = as.matrix(covars4[,c(1, 4:8)]),
  

  popMatCov[,i] <- projectMatCov(covarsa = covarProj[, c(1, 2, 4:8)], #covarProj,#covarsa = covarProj[,which(covarProjNames %in% covarList3a)],
                                 covarsb = covarProj[, c(2, 4, 6, 7)],#covarProj,#[,which(covarProjNames %in% covarList3b)],
                                 covarsS0 = covarProj[,c(4, 5, 7, 8)],#covarProj,#covarsS0 =  covarProj[,which(covarProjNames %in% covarList3S0)],
                                 covarsS1 = covarProj[,c(1, 4:8)],#covarProj,#covarsS1 =  covarProj[,which(covarProjNames %in% covarList3S1)],
                                 n = 36,
                                 year.start = 1)
}

plot(apply(popMatCov, 1, FUN = function(x){summary(x)[3]}), 
     log = "y",
     xaxt = "n")
axis(side = 1, at = seq(1, 36, by = 5),
     labels = seq(1980, 2015, by = 5))
points(apply(GraphMod$BUGSoutput$mean$N.est, 1, sum),
       pch = "*", cex = 3)

plot(apply(GraphMod$BUGSoutput$mean$N.est, 1, sum),
     apply(popMatCov, 1, FUN = function(x){summary(x)[3]}),
     xlab = "states",
     ylab = "calculations",
     log = c("x", "y"))



png("./Figures/PVA_ab_timeCov_hist.png",
    height = 4, width = 6, units = "in", res = 300)
### !!! DEPENDS ON WHAT QE IS SET TO ABOVE !!!
plot(0, 0, type = "n",
     xlim = c(1, 36),
     ylim = c(0, 250000),
     xlab = "",
     ylab = "Total Estimated Abundance",
     xaxt = "n")
text(paste0("QE = ", rownames(QE)),
     x = 30, y = 200000)
lines(apply(popMatCov, 1, 
            FUN = function(x){summary(x)[3]}), 
      lty = 1, lwd = 3)
# 3rd Q
lines(apply(popMatCov, 1, 
            FUN = function(x){summary(x)[5]}),
      lty = 2, lwd = 2)
# 1st Q
lines(apply(popMatCov, 1, 
            FUN = function(x){summary(x)[2]}),
      lty = 2, lwd = 2)

axis(side = 1, at = seq(1, 36, by = 5),
     labels = seq(1980, 2015, by = 5)) #seq(2005, 2050, by = 5))
dev.off()

#### Investigate QE values? ####

png("./Figures/PVA_ab_timeCov_hist_QE.png",
    height = 4, width = 6, units = "in", res = 300)
plot(0, 0, type = "n",
     xlim = c(1, 36),
     ylim = c(0, 250000),
     xlab = "",
     ylab = "Total Estimated Abundance",
     xaxt = "n")
for(x in c(0.01, 0.10, 0.15, 0.20)){
  QE <- setQE(prob = x)
  popMatCov <- matrix(NA, nrow = 36, ncol = 1000)
  for(i in 1:ncol(popMatCov)){
    popMatCov[,i] <- projectMatCov(covarsa = covarProj[, c(1, 2, 4:8)], #covarProj,#covarsa = covarProj[,which(covarProjNames %in% covarList3a)],
                                   covarsb = covarProj[, c(2, 4, 6, 7)],#covarProj,#[,which(covarProjNames %in% covarList3b)],
                                   covarsS0 = covarProj[,c(4, 5, 7, 8)],#covarProj,#covarsS0 =  covarProj[,which(covarProjNames %in% covarList3S0)],
                                   covarsS1 = covarProj[,c(1, 4:8)],#covarsS1 =  covarProj[,which(covarProjNames %in% covarList3S1)],
                                   
      # covarsa = covarProj[,which(covarProjNames %in% covarList3a)],
      #                              covarsb = covarProj[,which(covarProjNames %in% covarList3b)],
      #                              covarsS0 =  covarProj[,which(covarProjNames %in% covarList3S0)],
      #                              covarsS1 =  covarProj[,which(covarProjNames %in% covarList3S1)],
                                   n = 36,
                                   year.start = 1)
  }
  lines(apply(popMatCov, 1, 
              FUN = function(x){summary(x)[3]}), 
        lty = 1, lwd = 3,
        col = heat.colors(4)[which(c(0.01, 0.10, 0.15, 0.20) == x)])
}
axis(side = 1, at = seq(1, 36, by = 5),
     labels = seq(1980, 2015, by = 5)) #seq(2005, 2050, by = 5))
legend("topright", 
       lty = 1, lwd = 3, col = heat.colors(4),
       legend = c(0.01, 0.10, 0.15, 0.20),
       bty = "n")
dev.off()



#### Plot SFBS MWT.0 against QE-specific preditions ####
par(cex = 1.25, new = FALSE,
    mar = c(4.1, 4.1, 2.1, 2.1),
    fig = c(0, 1, 0, 1))
plot(0, 0, type = "n",
     xlim = c(0, 5000), # !!! cutting off the outliers
     ylim = c(0, 30000),
     xlab = "SFBS Age-0 MWT",
     ylab = "Abundance Estimate")
abline(a = 0, b = 1, lty = 2)
for(x in c(0.1, 0.15, 0.17)){#c(0.01, 0.10, 0.15)){ # take out 0.2 because we know it's too low
  QE <- setQE(prob = x)
  popMatCov <- matrix(NA, nrow = 36, ncol = 1000)
  for(i in 1:ncol(popMatCov)){
    popMatCov[,i] <- projectMatCov(covarsa = covarProj,#[,which(covarProjNames %in% covarList3a)],
                                   covarsb = covarProj,#[,which(covarProjNames %in% covarList3b)],
                                   covarsS0 =  covarProj,#[,which(covarProjNames %in% covarList3S0)],
                                   covarsS1 =  covarProj,#[,which(covarProjNames %in% covarList3S1)],
                                   n = 36,
                                   year.start = 1)
  }
  points(sfbs$mwt.0[1:36],
         apply(popMatCov, 1, 
               FUN = function(x){summary(x)[3]}), 
         lty = 1, lwd = 3,
         col = heat.colors(4)[which(c(0.1, 0.15, 0.17) == x)])
} 
legend("topright", 
       lty = 1, lwd = 3, col = heat.colors(4)[1:3],
       legend = c(0.1, 0.15, 0.17), #c(0.01, 0.10, 0.15),
       bty = "n")


# *** Fig. 10 *** ####
for(x in c(0.01, 0.15, 0.20)){
  QE <- setQE(prob = x)
  popMatCov <- matrix(NA, nrow = 36, ncol = 1000)
  for(i in 1:ncol(popMatCov)){
    popMatCov[,i] <- projectMatCov(covarsa = covarProj[, c(1, 2, 4:8)], 
                                   covarsb = covarProj[, c(2, 4, 6, 7)],
                                   covarsS0 = covarProj[,c(4, 5, 7, 8)],
                                   covarsS1 = covarProj[,c(1, 4:8)],
                                   n = 36,
                                   year.start = 1)
  }
  png(paste0("./Figures/Fig10_PVA_ab_timeCov_hist_QE",
  x, ".png"),
      height = 4, width = 6, units = "in", res = 300)
  plot(0, 0, type = "n",
       xlim = c(1, 36),
       ylim = c(0, 250000),
       xlab = "",
       ylab = "Total Estimated Abundance",
       xaxt = "n")
  lines(apply(popMatCov, 1, 
              FUN = function(x){summary(x)[3]}), 
        lty = 1, lwd = 3,
        col = "black") #heat.colors(4)[which(c(0.01, 0.10, 0.15, 0.20) == x)])
  # 3rd Q
  lines(apply(popMatCov, 1, 
              FUN = function(x){summary(x)[5]}),
        lty = 2, lwd = 2,
        col = "black") #heat.colors(4)[which(c(0.01, 0.10, 0.15, 0.20) == x)])
  # 1st Q
  lines(apply(popMatCov, 1, 
              FUN = function(x){summary(x)[2]}),
        lty = 2, lwd = 2,
        col = "black") #heat.colors(4)[which(c(0.01, 0.10, 0.15, 0.20) == x)])
  
  axis(side = 1, at = seq(1, 36, by = 5),
     labels = seq(1980, 2015, by = 5)) #seq(2005, 2050, by = 5))
  text(x = 15, y = 240000,
       paste("QE =", x),
       cex = 1.5, col = "black")
  points(sfbs$mwt.0[1:36], pch = "*", col = "grey")
  points(sfbs$ot.0[1:36], pch = "*", col = "lightblue")
  legend("topright",
         pch = "*", col = c("grey", "lightblue"),
         legend = c("MWT", "OT"),
         bty = "n",
         title = "SFBS Age-0 Indices")
# legend("topright", 
#        lty = 1, lwd = 3, col = heat.colors(4),
#        legend = c(0.01, 0.10, 0.15, 0.20),
#        bty = "n") 
dev.off()
}


# These QE thresholds don't match the figure in the TN
png("./Figures/QE_prob_timeCov.png",
    height = 4, width = 6, units = "in", res = 300)
plot(1:nrow(popMatCov),
     apply(popMatCov, 1, 
           FUN = function(x){sum(x == 0)})/ncol(popMatCov),
     type = "n",
     xaxt = "n",
     xlab = "",
     ylab = "Probability",
     ylim = c(0, 1),
     main = "Quasi-Extinction Risk")
QE <- setQE(prob = 0.1)
lines(#1:nrow(popMatCov),
      apply(makePopMatCov(), 1, FUN = function(x){sum(x == 0)})/ncol(popMatCov),
      lwd = 2, col = brewer.pal(5, "YlOrRd")[5])
QE <- setQE(prob = 0.05)
lines(#1:nrow(popMatCov),
      apply(makePopMatCov(), 1, FUN = function(x){sum(x == 0)})/ncol(popMatCov),
      lwd = 2, col = brewer.pal(5, "YlOrRd")[4])
QE <- setQE(prob = 0.025)
lines(#1:nrow(popMatCov),
      apply(makePopMatCov(), 1, FUN = function(x){sum(x == 0)})/ncol(popMatCov),
      lwd = 2, col = brewer.pal(5, "YlOrRd")[3])
QE <- setQE(prob = 0.01)
lines(#1:nrow(popMatCov),
      apply(makePopMatCov(), 1, FUN = function(x){sum(x == 0)})/ncol(popMatCov),
      lwd = 2, col = brewer.pal(5, "YlOrRd")[2])
axis(side = 1, at = seq(1, 66, by = 5),
     labels = seq(2005, 2070, by = 5))
legend("topleft", lwd = 2, 
       col = rev(brewer.pal(5, "YlOrRd")[2:5]),
       legend = rev(c(1, 2.5, 5, 10)),
       title = "Percentile",
       cex = 0.75)
dev.off()


#### 6.3 PVA Tables ####
# *** Table 2 *** ####
# These thresholds match Table 2 in the TN
# The values are slightly different, though.
QEtab <- round(rbind(setQE(prob = 0.5),
                     setQE(prob = 0.15),
                     setQE(prob = 0.10),
                     setQE(prob = 0.05),
                     setQE(prob = 0.01)),
               0)

plot(Nt$Nt)

plot(Nt$N0t)
plot(Nt$N1t)
plot(Nt$N2t)


#### 7. Model Checking ####
library(coda)
library(MCMCvis)
# get convergence stats
print(bhMod.cov, digits = 3)

MCMCtrace(bhMod.cov)

# Gelman-Rubin stats
MCMCsummary(bhMod.cov)

coda::gelman.diag(bhMod.cov)
