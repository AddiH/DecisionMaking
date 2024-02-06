PVL <- function(payoff, ntrials, w, A, a, theta) {
  # arrays to populate for simulation
  x <- array(NA, c(ntrials)) # chosen deck (A, B, C or D)
  X <- array(NA, c(ntrials)) # the payoff for the trial
  u <- array(NA, c(ntrials)) # the utility of the payoff
  Ev <- array(NA, c(ntrials, 4)) # expected value of each deck
  exp_p <- array(NA, c(ntrials, 4)) # exponential of expected values
  p <- array(NA, c(ntrials, 4)) # probability of choosing each option

  x[1] <- rcat(1, c(.25, .25, .25, .25)) # assigning a "flat" probability structure to the first choice (i.e. random choice between the four decks)

  X[1] <- payoff[1, x[1]] # assigning the payoff to first random choice

  Ev[1, ] <- rep(0, 4) # assigning zero as the expected value for all four decks at the first "random" choice

  for (t in 2:ntrials) {
    for (d in 1:4) {
      if (x[t - 1] == d) { # if the current deck was chosen last trial

        u[t] <- ifelse(X[t - 1] < 0, -w * abs(X[t - 1])^A, X[t - 1]^A) # calculating subjective utiliy
        Ev[t, d] <- Ev[t - 1, d] + (a * (u[t] - Ev[t - 1, d])) # set the updated EV
      
      } else {
        
        Ev[t, d] <- Ev[t - 1, d] # copy the previous EV
      }

      exp_p[t, d] <- exp(theta * Ev[t, d]) # probability of choosing deck
    }

    for (d in 1:4) {
      p[t, d] <- exp_p[t, d] / sum(exp_p[t, ]) # softmax
    }

    x[t] <- rcat(1, p[t, ]) # find the chosen deck
    X[t] <- payoff[t, x[t]] # get the corresponding payoff
  }

  result <- list(
    x = x,  # chosen deck
    X = X,  # corresponding payoff
    Ev = Ev # Expected utility
  ) 

  return(result)
}
