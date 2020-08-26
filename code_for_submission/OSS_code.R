library(tidyverse)

payoff <- function(s1, s2){
  # payoff to individual s1 against individual s2
  # s1 and s2 are two individuals, each represented as a vector with two entries (x, y)
  if (family == "quadratic") {
    bx <- (s1[1] + s2[1]) - b * (s1[1] + s2[1])^2
    by <- (s1[2] + s2[2]) - b * (s1[2] + s2[2])^2
    cx <- c * (s1[1] - d * (s1[1])^2)
    cy <- c * (s1[2] - d * (s1[2])^2)

    return( (1-alpha)*(bx-cx) + alpha*(by-cy))
  } else if (family == "power") {
    # didn't use
  }
}

equilibrium_frequencies <- function(s){
  P <- sapply(1:nrow(s), function(partner) sapply(1:nrow(s), function(focal) payoff(s[focal,], s[partner,])) ) # Payoff matrix. Pij is the payoff of strain i against strain j.
  P <- matrix(as.numeric(P), ncol=nrow(s))
  f <- solve(P) %*% rep(1, nrow(s)) # Use payoff matrix to calculate equilibrium frequencies: (Inverse of P * column vector of 1s)
  return(f / sum(f))
}

invasion_fitness <- function(invader, population){
  mutant_payoffs <- sapply(1:nrow(population), function(r) payoff(invader, population[r,])) %>% as.numeric # payoff of mutant against each resident strain
  resident_payoffs <- sapply(1:nrow(population), function(r) payoff(population[1,], population[r,])) %>% as.numeric # all residents have the same fitness (by definition) so we pick the 1st resident arbitrarily
  mutant_fitness <- sum(mutant_payoffs * f)
  resident_fitness <- sum(resident_payoffs * f)
  return((mutant_fitness - resident_fitness)/mutant_fitness)  # mutant invasion fitness at n-hat
}

for(theta in c(45, 90, 135, 180, 225)) {
  # Parameters --------------------------------------------------------------
  b <- 1.05; c <- 0.9; d <- 1.65; ess <-  - (c - 1)/(2*(2*b-c*d))
  mu <- 1e-2; # mutation rate.
  sd <- 1.5e-3; # mutation st dev.
  starting_distance <- ess/2
  epsilon <- 1e-4 # extinct strain threshold
  family <- "quadratic"
  tmax <- 100
  population_size <- 10000
  alpha <- 0.5
  theta <- 225

  x0 <- ess + starting_distance * cos(theta/180 * acos(-1)) # initial population
  y0 <- ess + starting_distance * sin(theta/180 * acos(-1))

  nr_replicates <- 9000
  output_list <- vector(mode = "list", length = nr_replicates)
  for(simulation in 1:nr_replicates) {
    fitness_coefficient <- 100 # magnifies small fitness differences


    # Simulation --------------------------------------------------------------
    t <- 0
    s <- tibble(x = x0, y = y0)

    # Calculate initial, equilibrium population sizes
    f <- equilibrium_frequencies(s); n <- f * population_size

    invasion_attempt <- 0
    output <- cbind(s, f, t, invasion_attempt, fitness_coefficient)
    row_output <- nrow(s)+1; output[row_output:row_output+100000,] <- NA

    while(t <= tmax){
      # Mutant emergence
      w <- mu * n; W <- sum(w); # rate of mutant emergence

      invasion_attempt <- 1; inv_attempts_since_change <- 1 # iterators used to optimise the fitness coefficient
      repeat{

        if( inv_attempts_since_change > 10 ) {
          fitness_coefficient <- fitness_coefficient*5  # if we're taking too long to get a successful invasion, decrease the fitness coefficient
          inv_attempts_since_change <- 1
        }

        # Introduce mutant
        i <- sample(1:nrow(s), size = 1, prob = w/W) # parent lineage
        mutant <- rnorm(2, mean = as.numeric(s[i,]), sd = sd) # strategy si'
        mutant[mutant < 0] <- 0
        if(family == "quadratic"){
			mutant[mutant > 1/(4*b)] <- 1/(4*b) # value of x that maximizes benefit(2x)
			mutant[mutant > 1/(2*d)] <- 1/(2*d) # value of x that maximizes cost(x)
		}
        # Update time
        t <- t - 1/W * log(runif(1))

        # Invasion fitness of mutant
        mutant_inv_fitness <- invasion_fitness(mutant, s)
        mutant_inv_fitness <- 1-exp(-mutant_inv_fitness * fitness_coefficient) # transform fitness so that small differences are magnified, to speed up code

        # Update iterators and check if invasion was successful
        invasion_attempt <- invasion_attempt + 1; inv_attempts_since_change <- inv_attempts_since_change + 1;
        if( runif(1) < mutant_inv_fitness ) { # if random_nr < f(si') invasion was successful
          if(invasion_attempt < 5) fitness_coefficient <- fitness_coefficient / 5 # if we reached invasion in few attempts, decrease the fitness coefficient
          break
        }
      }

      # Are the mutant and the resident mutually invadable??
      # Create imaginary population where mutant is the resident
      s_prime <- s; s_prime[i,] <- mutant
      f_prime <- equilibrium_frequencies(s_prime); n_prime <- f_prime * population_size

      # Calculate invasion fitness of the resident (invader) against the mutant (noninvader)
      resident_inv_fitness <- invasion_fitness(s[i,], s_prime) # invasion fitness at n-hat of the resident against the mutant

      if (!any(n_prime < epsilon) & resident_inv_fitness < 0) {
        s <- s_prime; n <- n_prime; f <- f_prime # mutant replaces resident
      } else {
        s_prime <- rbind(s, mutant) # mutant invades and coexists
        n_prime <- rep(epsilon+1, nrow(s_prime)) # temporary n_prime used so that no rows are removed in the first round of the repeat loop below

        # Calculate equilibrium population sizes
        f_prime <- equilibrium_frequencies(s_prime); n_prime <- f_prime * population_size

        s <- s_prime; n <- n_prime; f <- f_prime # mutant coexists with resident
      }

      while(any(n < epsilon)) {
        s <- s[n > epsilon, ]

        # Calculate equilibrium population sizes
        f <- equilibrium_frequencies(s); n <- f * population_size
      }

      output[row_output:(row_output+nrow(s)-1),] <- cbind(s, f, t, invasion_attempt, fitness_coefficient)
      row_output <- row_output + nrow(s)

      if (max(s$x) > ess+starting_distance | min(s$x) < ess-starting_distance | max(s$y) > ess+starting_distance | min(s$y) < ess-starting_distance) t <- tmax+1
    }

    output_list[[simulation]] <- (output %>% drop_na())
    print(paste0("Theta ", theta, ", replicate ", simulation, " complete!"))
  }

  saveRDS(output_list, file = paste0("alpha_", alpha, "_theta_", theta, "_output_list"))
}

