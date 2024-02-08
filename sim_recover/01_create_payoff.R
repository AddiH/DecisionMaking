source("utils/payoff.R")

n_trials <- 100
div <- 10

# Bad deck
# Large reward, frequent small loss
A <- create_deck(
  n_trials = n_trials,
  div = div,
  reward = 100,
  loss = -250,
  freq_of_loss = .5
)

# Bad deck
# Large reward, infrequent large loss
B <- create_deck(
  n_trials = n_trials,
  div = div,
  reward = 100, 
  loss = -1250,
  freq_of_loss = .1
)

# Good deck
# Small reward, frequent tiny loss
C <- create_deck(
  n_trials = n_trials,
  div = div,
  reward = 50, 
  loss = -50,
  freq_of_loss = .5
)

# Good deck
# Small reward, infrequent small loss
D <- create_deck(
  n_trials = n_trials,
  div = div,
  reward = 50,
  loss = -250,
  freq_of_loss = .1
)

payoff <- cbind(A, B, C, D) # combine

write.csv(payoff, "data/payoff.csv", row.names=FALSE)