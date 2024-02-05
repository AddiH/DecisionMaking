create_deck <- function(n_trials, div, reward, loss, freq_of_loss) {
  wins <- rep(reward, div) # make wins
  loss <- rep(loss, div*freq_of_loss) # make losses
  no_loss <- rep(0, div*(1-freq_of_loss)) # make no losses
  losses <- c(loss, no_loss) # combine the two loss lists
  outcome <- wins + losses # add wins and losses
  
  deck <- c() # empty deck
  
  for (i in 1:(n_trials / div)) {  # Iterate over each subdivision
    deck <- c(deck, sample(outcome))  # Shuffle the outcomes and append to deck
  }
  
  deck <- deck/100 # divide the outcomes
  return(deck)
}