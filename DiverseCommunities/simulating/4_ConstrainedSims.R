########################################################
# Diverse community simulations while varying the
# correlations between pairwise and higher-order
# interactions for panel (C) of Fig. 4
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

# choosing parameters to iterate over
S <- 20
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(-0.5, 0, 0.5) / S
sigma_A <- seq(0.5, 1.5, length.out = 3)
mu_B <- 0
sigma_B <- 0
rho_A <- 0
intra <- "None"
self_reg <- "Quadratic"
scaling <- FALSE
dist_B <- "Normal"

# building a data frame of the input parameters
in_pars <- crossing(S = S,
                    MuR = mu_R,
                    SigmaR = sigma_R,
                    MuD = mu_D,
                    SigmaD = sigma_D,
                    MuA = mu_A,
                    SigmaA = sigma_A,
                    RhoA = rho_A,
                    MuB = mu_B,
                    SigmaB = sigma_B,
                    Intra = intra,
                    SelfReg = self_reg,
                    scaling = scaling,
                    DistB = dist_B)

num_repl <- 100 # choosing the number of replicates
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

out_data <- data.frame()

# choosing correlations between the pairwise
# and higher-order interaction matrices
ps <- seq(0.1, 1, length.out = 10)

start_time <- Sys.time()

# looping over all the parameters
for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_stats)
  
  # looping over the correlations
  for(p in ps) {
    
    # setting the target abundances
    target_abd <- GetTargetAbd(pars = cur_pars)
    
    # building the Jacobian from A and Btilde
    J <- BuildAbbJacobian(abd = target_abd, pars = cur_pars, p = p)
    
    # recording stability and interaction statistics
    cur_eig <- GetEig(J)
    cur_stable <- (cur_eig < 0)
    cur_data <- data.frame(Stable = cur_stable, p = p)
    cur_data <- cbind(cur_stats, cur_data)
    
    out_data <- rbind(out_data, cur_data)
  }
}

print(Sys.time() - start_time)

# writing out the data
filename <- "SimCorr"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


