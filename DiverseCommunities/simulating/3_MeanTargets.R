########################################################
# Non-equal abundances mean varying simulations  for
# panel (A) of Fig. 3
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

### Simulation script

# choosing parameters to iterate over
S <- 20
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(-0.5, 0.5) / S
sigma_A <- 0.15
rho_A <- 0
mu_B <- 0
sigma_B <- 0
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

# choosing the mean abundance values
abundance_vals <- seq(0.01, 1.99, length.out = 10)

num_repl <- 100 # choosing the number of replicates
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

out_data <- data.frame()

start_time <- Sys.time()

# looping over all the parameters
for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_stats)
  
  # looping over the different mean abundance values
  for(abd_val in abundance_vals) {
    
    # setting the target abundances and sampling the constrained HOIs
    target_abd <- GetTargetAbd(pars = cur_pars, choice = "Specified", value = abd_val)
    cur_pars$B <- GetFeasibleB(target_abd, cur_pars)
    
    J <- BuildJacobian(eq_abd = target_abd, pars = cur_pars)
    
    # recording the eigenvalues and feasibility
    cur_eig <- GetEig(J)
    cur_stable <- (cur_eig < 0)
    feas_data <- data.frame(Stable = cur_stable, TargetAbd = abd_val, Type = "Carrying Capacity")
    
    # recording the interaction statistics
    cur_data <- cbind(cur_stats, feas_data)
    out_data <- rbind(out_data, cur_data)
  }
}

print(Sys.time() - start_time)

# writing out the data
filename <- "MeanTargets"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)
