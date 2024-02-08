load('data/all_studies.rdata')

data_list <- list(
  Fridberg = Fridberg, 
  Horstmann = Horstmann, 
  Kjome = Kjome, 
  Maia = Maia, 
  SteingroverInPrep = SteingroverInPrep, 
  Premkumar = Premkumar, 
  Wood = Wood, 
  Worthy = Worthy, 
  Steingroever2011 = Steingroever2011, 
  Wetzels = Wetzels
)

for (i in 1:length(data_list)) {
  df <- data_list[[i]]
  
  # identify and count unique subject IDs
  subIDs <- unique(df$Subj)
  n_subs <- length(subIDs)
  n_trials <- length(df$Subj)/n_subs
  
  # empty arrays to fill
  ntrials_all <- array(0,c(n_subs))
  x_all <- array(0,c(n_subs,n_trials))
  X_all <- array(0,c(n_subs,n_trials))
  
    for (s in 1:n_subs) {
      ntrials_all[s] <- length(df$x[df$Subj==subIDs[s]])
      
      x_sub <- df$x[df$Subj==subIDs[s]] 
      X_sub <- df$X[df$Subj==subIDs[s]] 
      
      # assign arrays
      x_all[s,] <- x_sub
      X_all[s,] <- X_sub
      
    }
  
  save_path <- paste('estimations/', names(data_list)[i],'.rdata', sep = '')
  
  if (file.exists(save_path)) {
    # If the path exists, skip to the next iteration of the loop
    next
  }
  x <- x_all
  X <- X_all
  
  n_trials <- ntrials_all # this pisses me off
  
  # set up jags and run jags model
  data <- list("x","X","n_trials","n_subs") 
  params<-c("mu_w","mu_A","mu_theta","mu_a","lambda_w","lambda_A","lambda_theta","lambda_a")
  
  start_time = Sys.time()
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="utils/hier_PVL.txt",
                           n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)
  
  end_time = Sys.time()
  save(samples, file = save_path)
  print(paste(names(data_list)[i],'complete. Elapsed time:', end_time - start_time))
}