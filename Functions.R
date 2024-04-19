########################################################
# Functions
########################################################

# Loading packages

library(tidyverse)
library(deSolve)
library(reshape2)
library(MASS)
library(nleqslv)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(sads)
library(viridis)

# Functions

# integrates the dynamics using deSolve
IntegrateDynamics <- function(inistate, pars, endtime, timestep, fn){
  times <- seq(0, endtime, by = timestep)
  timeseries <- as.data.frame(ode(inistate, times, fn, pars))  
  return(timeseries)
}

# plots time series
PlotSeries <- function(series, title = "Community Dynamics") {
  meltseries <- melt(series, id.vars = "time")
  pl <- ggplot() +
    geom_line(data = meltseries, aes(x = time, y = value, color = variable),
              size = 3, alpha = 0.5) +
    ggtitle(title) +
    xlab("Time") + ylab("Abundance") +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size=30),
          strip.text.x = element_text(size = 25),
          strip.text.y = element_text(size = 25),
          legend.text=element_text(size = 25))
  return(pl)
}

# Finds equilibrium coexistence solution of the LV model
GetPairwiseAbd <- function(pars) {
  abd <- with(pars, {
    abd <- solve(A) %*% r
    return(abd)
  })
  return(abd)
}

# Computes the leading eigenvalue of a matrix
GetEig <- function(J) {
  eigs <- eigen(J, only.values = TRUE)$values
  re_eigs <- Re(eigs)
  lambda <- max(re_eigs)
  return(lambda)
}

# Samples coefficients from bivariate normal dist
SampleCoefficientsFromBivariate <- function(S, mux, muy, sigmax, sigmay, rho){
  NumPairs <- S * (S - 1) / 2
  mus <- c(mux, muy)
  covariance.matrix <- matrix(c(sigmax^2,
                                rho * sigmax * sigmay,
                                rho * sigmax * sigmay,
                                sigmay^2),
                              2, 2)
  Pairs <- mvrnorm(NumPairs, mus, covariance.matrix)
  return(Pairs)
}

# Builds the pairwise interaction matrix
BuildA <- function(S = 100, mu = 0, sigma = 0.1, rho = 0){
  
  Pairs <- SampleCoefficientsFromBivariate(S, mux = mu, muy = mu,
                                           sigmax = sigma, sigmay = sigma,
                                           rho = rho)
  
  A <- matrix(0, S, S)
  k <- 1
  for (i in 1:(S-1)){
    for (j in (i + 1):S){
      A[j,i] <- Pairs[k,1]
      A[i,j] <- Pairs[k,2]
      k <- k + 1
    }
  }
  return(A)
}

# returns the target abundance vector given the input parameterization
GetTargetAbd <- function(pars, choice = "Carrying Capacities", value = 0, min_val = 1e-2) {
  
  # returns the carrying capacioties of the model
  if(choice == "Carrying Capacities") {
    car_cap <- unique(rep(pars$r / diag(pars$A)))
    target_abd <- car_cap + rnorm(pars$S, mean = 0, sd = value)
    target_abd[target_abd <= 0] <- min_val
    target_abd <- target_abd / mean(target_abd) * car_cap
    # simply returns the input values
  } else if(choice == "Specified") {
    target_abd <- rep(value, times = pars$S)
    # computes the equilibrium of the pairwise only model and returns this vector plus noise
  } else if(choice == "Pairwise + Noise") {
    target_abd <- GetPairwiseAbd(pars)
    if(sum(target_abd < 0)) {
      target_abd <- GetAbds(pars = pars,
                            ini_state = runif(pars$S, min = 0, max = 0.01),
                            end_time = 1e5,
                            zero_cutoff = 1e-7)
      target_abd <- target_abd$Abd
      target_abd[target_abd < 1e-7] <- -1
    }
    target_abd <- target_abd + rnorm(pars$S, mean = 0, sd = value)
    target_abd[target_abd <= 0] <- min_val
    # the following parameterizations imlement different realistitc species abundance distributions
  } else if(choice == "Log Normal") {
    target_abd <- rlnorm(n = pars$S, sdlog = value)
    car_cap <- unique(rep(pars$r / diag(pars$A)))
    target_abd[target_abd <= 0] <- min_val
    target_abd <- target_abd / mean(target_abd) * car_cap
  } else if(choice == "Poisson Log Normal") {
    car_cap <- unique(rep(pars$r / diag(pars$A)))
    target_abd <- rpoilog(n = pars$S, mu = car_cap, sig = value)
    target_abd[target_abd <= 0] <- min_val
    target_abd <- target_abd / mean(target_abd) * car_cap
  } else if(choice == "Log Series") {
    car_cap <- unique(rep(pars$r / diag(pars$A)))
    target_abd <- rls(n = pars$S, N = 10 * pars$S, alpha = 1 / value)
    target_abd[target_abd <= 0] <- min_val
    target_abd <- target_abd / mean(target_abd) * car_cap
  } else if(choice == "Geometric Series") {
    car_cap <- unique(rep(pars$r / diag(pars$A)))
    target_abd <- rgs(n = pars$S, k = 1 - value, S = pars$S)
    target_abd[target_abd <= 0] <- min_val
    target_abd <- target_abd / mean(target_abd) * car_cap
  } else if(choice == "Zipf") {
    car_cap <- unique(rep(pars$r / diag(pars$A)))
    target_abd <- rzipf(n = pars$S, N = 100 * pars$S, s = 5 - 3 * value)
    target_abd[target_abd <= 0] <- min_val
    target_abd <- target_abd / mean(target_abd) * car_cap
  }
  return(target_abd)
}

# plots heatmaps of matrices like those in Fig. 4
GridM <- function(M, Title, in_color, in_max) {
  hm.palette <- colorRampPalette(brewer.pal(9, in_color), space='Lab')
  M.melted <- melt(t(M))
  gridM <- ggplot(M.melted, aes(x = Var1, y = dim(M)[1] - Var2, fill = value)) +
    geom_tile() +
    ggtitle(Title) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid = element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 17.5)) +
    coord_equal() +
    scale_fill_gradientn(
      colours=hm.palette(100), limits = c(-in_max,in_max)
    )
  return(gridM)
}
