install.packages("pacman")
pacman::p_load(tidyverse, R2jags, extraDistr, parallel, truncnorm)
source('utils/hier_PVL_sim.R')
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}
payoff <- read.csv("data/payoff_100_trials.csv")
payoff <- payoff/100
niterations <- 100 # fewer because it takes too long
n_subs <- 60 # avg subs in dataset
n_trials <- 100
ntrials_all <- rep(n_trials, n_subs) # all 60 subs have 100 trials each

# mu
true_mu_w <- array(NA,c(niterations))
true_mu_A <- array(NA,c(niterations))
true_mu_theta <- array(NA,c(niterations))
true_mu_a <- array(NA,c(niterations))

infer_mu_w <- array(NA,c(niterations))
infer_mu_A <- array(NA,c(niterations))
infer_mu_theta <- array(NA,c(niterations))
infer_mu_a <- array(NA,c(niterations))

# sigma (SD for R) / lambda (precision for JAGS)
true_lambda_w <- array(NA,c(niterations))
true_lambda_A <- array(NA,c(niterations))
true_lambda_theta <- array(NA,c(niterations))
true_lambda_a <- array(NA,c(niterations))

infer_lambda_w <- array(NA,c(niterations))
infer_lambda_A <- array(NA,c(niterations))
infer_lambda_theta <- array(NA,c(niterations))
infer_lambda_a <- array(NA,c(niterations))

start_time = Sys.time()
for (i in 1:niterations) {
  n_trials <- ntrials_all
  
  mu_w <- runif(1,.5,2.5)
  mu_A <- runif(1,0,1)
  mu_theta <- runif(1,0,2)
  mu_a <- runif(1,0,1)
  
  sigma_w <- runif(1,0,0.2)
  sigma_A <- runif(1,0,0.1)
  sigma_theta <- runif(1,0,0.2)
  sigma_a <- runif(1,0,0.1)
  
  PVL_sims <- hier_PVL_sim(payoff,n_subs,n_trials,mu_w,mu_A,mu_a,mu_theta,
                           sigma_w, sigma_A, sigma_a, sigma_theta)
  
  x <- PVL_sims$x
  X <- PVL_sims$X
  
  # set up jags and run jags model
  data <- list("x","X","n_trials","n_subs") 
  params<-c("mu_w","mu_A","mu_theta","mu_a","lambda_w","lambda_A","lambda_theta","lambda_a")
  print(Sys.time())
  print("running jags...")
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="utils/hier_PVL.txt", n.chains=3, 
                           n.iter=3000, n.burnin=1000, n.thin=1, n.cluster=4)
  
  # mu
  true_mu_w[i] <- mu_w
  true_mu_A[i] <- mu_A
  true_mu_theta[i] <- mu_theta
  true_mu_a[i] <- mu_a
  
  # find maximum a posteriori
  Y <- samples$BUGSoutput$sims.list
  infer_mu_w[i] <- MPD(Y$mu_w)
  infer_mu_A[i] <- MPD(Y$mu_A)
  infer_mu_theta[i] <- MPD(Y$mu_theta)
  infer_mu_a[i] <- MPD(Y$mu_a)
  
  # sigma
  true_lambda_w[i] <- sigma_w
  true_lambda_A[i] <- sigma_A
  true_lambda_theta[i] <- sigma_theta
  true_lambda_a[i] <- sigma_a
  
  # find maximum a posteriori
  infer_lambda_w[i] <- MPD(Y$lambda_w)
  infer_lambda_A[i] <- MPD(Y$lambda_A)
  infer_lambda_theta[i] <- MPD(Y$lambda_theta)
  infer_lambda_a[i] <- MPD(Y$lambda_a)
  
  print(i)
  print(Sys.time())
  
}
end_time = Sys.time()
end_time - start_time

save(true_mu_w, true_mu_A, true_mu_theta, true_mu_a, infer_mu_w, infer_mu_A, infer_mu_theta, infer_mu_a, true_lambda_w, true_lambda_A, true_lambda_theta, true_lambda_a, infer_lambda_w, infer_lambda_A, infer_lambda_theta, infer_lambda_a,
     file = paste('data/hier_recover_100trials.csv'))


