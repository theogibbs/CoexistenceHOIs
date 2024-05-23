########################################################
# Functions for three species simulations
# For the three species case, each species only
# experiences one HOI in our parameterization, so we use
# a slightly simpler notation / coding format to keep
# track of these interactions. This file contains the
# functions for these simpler cases.
########################################################

# Functions

# computes the per capita growth rates of the model
ThreeSpeciesGrowthRates <- function(state, pars) {
  
  growth_rates <- with(pars, {
    outer_state <- outer(state, state)
    outer_state <- outer_state[upper.tri(outer_state)]
    outer_state <- rev(outer_state)
    growth_rates <- (r - A %*% state - b * outer_state)
    return(growth_rates)
  })
  
  return(growth_rates)
}

# returns a list of the derivatives of the model given the state and parameters
# using the growth rates function
ThreeSpeciesDynamics <- function(time, state, pars) {
  
  dNdt <- with(pars, {
    outer_state <- outer(state, state)
    outer_state <- outer_state[upper.tri(outer_state)]
    outer_state <- rev(outer_state)
    dNdt <- state * (r - A %*% state - b * outer_state)
    return(dNdt)
  })
  
  return(list(dNdt))
}

# builds a list of the parameters of the three species model
BuildThreeSpeciesPars <- function(in_pars) {
  
  pars <- with(in_pars, {
    pars <- list()
    
    pars$r <- rnorm(S, mean = MuR, sd = SigmaR)
    pars$A <- BuildA(S = S, mu = MuA, sigma = SigmaA, rho = RhoA)
    diag(pars$A) <- rnorm(S, mean = MuD, sd = SigmaD)
    pars$b <- rep(b, times = S)
    return(pars)
  })
  
  return(pars)
}

# creates the vector of the three higher
# order interactions and constrains them
# to counter the pairwise interactions
GetThreeSpeciesB <- function(abd, pars) {
  const_b <- with(pars, {
    
    if(length(abd) == 1) {
      abd <- rep(abd, times = 3)
    }
    
    outer_abd <- outer(abd, abd)
    outer_abd <- outer_abd[upper.tri(outer_abd)]
    outer_abd <- rev(outer_abd)
    
    const_b <- (r - A %*% abd) / outer_abd
    const_b <- as.vector(const_b)
    
    return(const_b)}
  )
  return(const_b)
}

# compues the three species Jacobian using the
# formula from the main text
BuildThreeSpeciesJacobian <- function(eq_abd, pars) {
  J <- with(pars, {
    Bmat <- matrix(b, 3, 3)
    diag(Bmat) <- 0
    
    if(length(eq_abd) == 1) {
      eq_abd <- rep(eq_abd, times = 3)
    }
    
    Bmat[1,] <- Bmat[1,] * c(0, eq_abd[3], eq_abd[2])
    Bmat[2,] <- Bmat[2,] * c(eq_abd[3], 0, eq_abd[1])
    Bmat[3,] <- Bmat[3,] * c(eq_abd[2], eq_abd[1], 0)

    J <- - eq_abd * (A + Bmat)
    return(J)
  })
  return(J)
}

# integrates the dynamics of the model and returns the ending state
# and the stability of this resulting state
GetThreeSpeciesAbds <- function(pars, ini_state, end_time, zero_cutoff) {
  out_abds <- IntegrateDynamics(inistate = ini_state,
                                pars = pars,
                                endtime = end_time,
                                timestep = end_time,
                                fn = ThreeSpeciesDynamics)
  end_pt <- as.numeric(out_abds[2,2:(length(pars$r)+1)])
  
  end_pt[end_pt < zero_cutoff] <- zero_cutoff^2
  end_pt[end_pt > 1/zero_cutoff] <- 1/zero_cutoff
  end_pt[is.na(end_pt)] <- 1/zero_cutoff
  eq_log <- unlist(ThreeSpeciesDynamics(0, end_pt, pars))
  eq_log <- as.logical(prod(eq_log < zero_cutoff))

  stab_log <- BuildThreeSpeciesJacobian(end_pt, pars)
  stab_log <- GetEig(stab_log)
  stab_log <- as.logical(stab_log)
  
  target_abd <- GetTargetAbd(pars)
  target_log <- all.equal(rep(target_abd, times = 3), end_pt)
  target_log <- isTRUE(target_log)
  
  out_abds <- data.frame(Abd = end_pt,
                         Equilibrium = eq_log,
                         Stability = stab_log,
                         TargetAbd = target_log)
  return(out_abds)
}
