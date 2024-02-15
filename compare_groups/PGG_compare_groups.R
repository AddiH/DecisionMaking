seed_id = 1983
set.seed(seed_id)

install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, ggplot2, glue, dplyr)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}


##### Preprocessing #####

groupSize <- 4
ntrials <- 10
pi <- 1.6 # used to be 1.4, but the original paper (and Josh' preprint both say 1.6)
ntokens <- 20
vals <- seq(0,ntokens,1)
#vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens

rawDat <- read.csv("data/HerrmannThoeniGaechterDATA.csv", skip = 3) # Public goods game

#- create covariates in raw data matrix
# nation index
rawDat$nation <- c()
rawDat$nation[rawDat$city=="Melbourne"]=1
rawDat$nation[rawDat$city=="Minsk"]=2
rawDat$nation[rawDat$city=="Chengdu"]=3
rawDat$nation[rawDat$city=="Copenhagen"]=4
rawDat$nation[rawDat$city=="Bonn"]=5
rawDat$nation[rawDat$city=="Athens"]=6
rawDat$nation[rawDat$city=="Seoul"]=7
rawDat$nation[rawDat$city=="Samara"]=8
rawDat$nation[rawDat$city=="Zurich"]=9
rawDat$nation[rawDat$city=="St. Gallen"]=9
rawDat$nation[rawDat$city=="Istanbul"]=10
rawDat$nation[rawDat$city=="Nottingham"]=11
rawDat$nation[rawDat$city=="Dnipropetrovs'k"]=12
rawDat$nation[rawDat$city=="Boston"]=13

# create variable for GINI. Old data from 
# http://hdr.undp.org/sites/default/files/reports/269/hdr_2009_en_complete.pdf,
# see the suppl material in the published paper for the updated values:
# https://www.sciencedirect.com/science/article/pii/S2666622723000254#sec0014

# # Old pre-published Gini-values
# rawDat$gini <- c()
# rawDat$gini[rawDat$city=="Melbourne"]=34.3
# rawDat$gini[rawDat$city=="Minsk"]=25.3
# rawDat$gini[rawDat$city=="Chengdu"]=38.5
# rawDat$gini[rawDat$city=="Copenhagen"]=28.7
# rawDat$gini[rawDat$city=="Bonn"]=31.9
# rawDat$gini[rawDat$city=="Athens"]=34.4
# rawDat$gini[rawDat$city=="Seoul"]=31.6
# rawDat$gini[rawDat$city=="Samara"]=37.5
# rawDat$gini[rawDat$city=="Zurich"]=32.7
# rawDat$gini[rawDat$city=="St. Gallen"]=32.7
# rawDat$gini[rawDat$city=="Istanbul"]=41.9
# rawDat$gini[rawDat$city=="Nottingham"]=34.8
# rawDat$gini[rawDat$city=="Dnipropetrovs'k"]=26.1
# rawDat$gini[rawDat$city=="Boston"]=41.1

# Gini-values taken directly from the supplementary materials in the published paper
rawDat$gini <- c()
rawDat$gini[rawDat$city=="Melbourne"]=33.3
rawDat$gini[rawDat$city=="Minsk"]=28.3
rawDat$gini[rawDat$city=="Chengdu"]=41.5
rawDat$gini[rawDat$city=="Copenhagen"]=25.4
rawDat$gini[rawDat$city=="Bonn"]=30.7
rawDat$gini[rawDat$city=="Athens"]=34
rawDat$gini[rawDat$city=="Seoul"]=31.7
rawDat$gini[rawDat$city=="Samara"]=38
rawDat$gini[rawDat$city=="Zurich"]=31.7
rawDat$gini[rawDat$city=="St. Gallen"]=31.7
rawDat$gini[rawDat$city=="Istanbul"]=41.4
rawDat$gini[rawDat$city=="Nottingham"]=35.2
rawDat$gini[rawDat$city=="Dnipropetrovs'k"]=29.1
rawDat$gini[rawDat$city=="Boston"]=40.8

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]

# subsetting to only two nations (here: 8 and 9)
redDat <- redDat %>% filter(nation %in% c(8,9))

group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

ngroups_nation <- redDat %>% group_by(nation) %>% summarize(ngroups = length(unique(groupid)))

# # THIS WILL REMOVE SUBJECTS WITH MISSING DATA IN NO PUNISHMENT CONDITION
# ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
c_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_no_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

missing <- array(0,ngroups)

for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Gga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    Gc_no_punish[s,,g] <- colSums(c_no_punish[-s,,g])
    Ga_no_punish[s,,g] <- colMeans(c_no_punish[-s,,g])
    Ggas_no_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# data for punishment condition #
c_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Gga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    Gc_punish[s,,g] <- colSums(c_punish[-s,,g])
    Ga_punish[s,,g] <- colMeans(c_punish[-s,,g])
    Ggas_punish[s,,g] <- colMeans(c_punish[,,g])
  }
}

# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Gga <- array(0,c(ntrials,ngroups,2))
Gga[,,1] <- Gga_no_punish
Gga[,,2] <- Gga_punish

Ggas <- array(0,c(groupSize,ntrials,ngroups,2))
Ggas[,,,1] <- Ggas_no_punish
Ggas[,,,2] <- Ggas_punish


Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

Ga <- array(0,c(groupSize,ntrials,ngroups,2))
Ga[,,,1] <- Ga_no_punish
Ga[,,,2] <- Ga_punish

c_choice_index <- c

Nation <- array(0, ngroups)
for (g in 1:ngroups) {
  Nation[g] <- mean(redDat$nation[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
}

# ignoring the groups with missing values
c <- c[,,!missing,]
Ga <- Ga[,,!missing,]

c_win <- c_no_punish
c_keep <- rep()-c_win

# calculate the winnings (i.e. apply the multiplication-factor to the sum of each groups contributions)
winnings <-  array(0, ngroups)
for (g in 1:ngroups) {
  winnings[g] <- sum(colSums(c_win[,,g])*pi)
}

# organizing data into two groups
# getting a count of the groups in each nation-group - subtracting those with missing values
ngroups_gr1 <- ngroups_nation$ngroups[1]-sum(missing[1:ngroups_nation$ngroups[1]])
ngroups_gr2 <- ngroups_nation$ngroups[2]-sum(missing[(ngroups_nation$ngroups[1]+1):
                                                       (ngroups_nation$ngroups[1]+ngroups_nation$ngroups[2])])

c_gr1 <- c[,,1:ngroups_gr1,]
c_gr2 <- c[,,(ngroups_gr1+1):
             (ngroups_gr1+ngroups_gr2),]

Ga_gr1 <- Ga[,,1:ngroups_gr1,]
Ga_gr2 <- Ga[,,(ngroups_gr1+1):
               (ngroups_gr1+ngroups_gr2),]



################################################################################
########################### Conditional cooperation model ######################
################################################################################

# JZS priors for partial correlation. Method described here
# https://link.springer.com/article/10.3758/s13423-012-0295-x
# Code available here
# https://github.com/MicheleNuijten/BayesMed/blob/master/R/jzs_corSD.R
# Paper where code is used here (mediation paper)
# https://link.springer.com/article/10.3758/s13428-014-0470-2

#################################################################
#------------------ Winnings analysis ---------------------------
#################################################################

### NB! WINNINGS NOT FIXED AS SUCH YET...
#-------------------  Regress Gini on winnings ---------------

# standardise variables

                     
Y <- (winnings-mean(winnings))/sd(winnings)

data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") 
params <- c("beta0","betaX") 

# - run jags code
win.samples <- jags.parallel(data, inits=NULL, params,
                    model.file ="win_corr.txt",
                    n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=3)

# #----------------------- Control for all four national variables using partial correlation _------------
# fullControl.win <- function (X1,X2,X3,X4,X5,ngroups,nnations,winnings,Nation) {
#   
#   # standardise covariates
#   X1 <- (X1-mean(X1))/sd(X1)
#   X2 <- (X2-mean(X2))/sd(X2)
#   X3 <- (X3-mean(X3))/sd(X3)
#   X4 <- (X4-mean(X4))/sd(X4)
#   X5 <- (X5-mean(X5))/sd(X5)
#   
#   X <- cbind(X1,X2,X3,X4,X5)
#   V <- solve(t(X)%*%X)
#   
#   Y <- (winnings-mean(winnings))/sd(winnings)
#   
#   data <- list("ngroups", "Y", "nnations","X","V","Nation") #data inputted into jags
#   params <- c("beta0","betaX","prior_T") #parameters we'll track in jags
#   
#   # - run jags code
#   win.samples <- jags(data, inits=NULL, params,
#                       model.file ="win_partcor_4control.txt",
#                       n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1)
#   
#   return(list(win.samples))
#   
# }
# 
# Full_control.win <- fullControl.win(GDP,Indiv,Trust,Civic,Gini,ngroups,nnations,winnings,Nation)

#################################################################
#------------------ CC model analysis ---------------------------
#################################################################

#-------------------  group diffs on belief weights and slope of prefs in CC model ---------------

# sees to work best with a 0,1 vector
#X <- c(-0.5, 0,5)
X <- c(0, 1)

X <- (X-mean(X))/sd(X)

invSigma <- solve(t(X)%*%X) # required for JZS priors

ngroups=length(X)

data <- list("groupSize", "ngroups_gr1", "ngroups_gr2", "ntrials", "X",
             "c_gr1", "c_gr2", "Ga_gr1", "Ga_gr2", "invSigma", "ngroups") 
params<-c("beta0_alpha","beta0_rho","beta0_omega",
          "betaX_alpha", "betaX_rho", "betaX_omega",
          "prec_alpha", "prec_rho", "prec_omega") 

# - run jags code
start_time = Sys.time()
CC.samples <- jags.parallel(data, inits=NULL, params,
                   model.file ="reparam_compare_JZS.txt",
                   n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=3)
end_time = Sys.time()
end_time - start_time

##### Plotting needs some work ;)

Y <- CC.samples$BUGSoutput$sims.list

##### Plotting the slopes (and intercepts) against their priors (still in linear space)

par(mfrow=c(2,2))
plot(density(rnorm(10000,0,1/sqrt(MPD(Y$prec_alpha)))),ylim=c(0,.7),main="alpha")
lines(density(Y$betaX_alpha),col="red")

plot(density(rnorm(10000,0,1/sqrt(MPD(Y$prec_rho)))),ylim=c(0,.7),main="rho")
lines(density(Y$betaX_rho),col="red")

plot(density(rnorm(10000,0,1/sqrt(MPD(Y$prec_omega)))),ylim=c(0,.7),main="omega")
lines(density(Y$betaX_omega),col="red")


par(mfrow=c(2,2))
plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="alpha")
lines(density(Y$beta0_alpha),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="rho")
lines(density(Y$beta0_rho),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="omega")
lines(density(Y$beta0_omega),col="red")



##### Plotting credible intervals
par(mfrow=c(2,2))
#### ALPHA slope (betaX) (= initial belief)
# estimating and plotting Bayesian credible interval (95%)
MAP_betaX_alpha <- MPD(Y$betaX_alpha)

betaX_alpha_cred = quantile(Y$betaX_alpha,c(0.025,0.975))

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(Y$betaX_alpha),
     main = "alpha slope - initial belief")
lines(betaX_alpha_cred, c(0,0), col="red")
points(MAP_betaX_alpha, 0, pch=19, col="red")


#### RHO slope (betaX) (= degree of undermatching, aka. conditional cooperation)
# estimating and plotting Bayesian credible interval (95%)
MAP_betaX_rho <- MPD(Y$betaX_rho)

betaX_rho_cred = quantile(Y$betaX_rho,c(0.025,0.975))

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(Y$betaX_rho),
     main = "rho slope - conditional cooperation")
lines(betaX_rho_cred, c(0,0), col="red")
points(MAP_betaX_rho, 0, pch=19, col="red")


#### OMEGA slope (betaX) (= attention to others)
# estimating and plotting Bayesian credible interval (95%)
MAP_betaX_omega <- MPD(Y$betaX_omega)

betaX_omega_cred = quantile(Y$betaX_omega,c(0.025,0.975))

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(Y$betaX_omega),
     main = "omega slope - attention to others")
lines(betaX_omega_cred, c(0,0), col="red")
points(MAP_betaX_omega, 0, pch=19, col="red")


#####
##### simple group comparison

data <- list("groupSize", "ngroups_gr1", "ngroups_gr2", "ntrials", "X",
             "c_gr1", "c_gr2", "Ga_gr1", "Ga_gr2")

params<-c("diff_alpha","diff_rho","diff_omega",
          "mu_alpha", "mu_rho", "mu_omega") 

# - run jags code
start_time = Sys.time()
CC.samples <- jags.parallel(data, inits=NULL, params,
                            model.file ="compare_groups/simple_compare.txt",
                            n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=3)
end_time = Sys.time()
end_time - start_time

