########################################################
# Non-equal abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

### Simulation script

S <- 20
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- seq(-12.5, 5, length.out = 10) / S
sigma_A <- c(0, 0.1)
rho_A <- 0
mu_B <- 0
sigma_B <- 0
intra <- "None"
self_reg <- "Quadratic"
scaling <- FALSE
dist_B <- "Normal"

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

num_repl <- 20
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

in_dists <- c("Carrying Capacities", "Log Normal", "Geometric Series", "Zipf")
sd_vals <- c(0.1, 0.75)

out_data <- data.frame()

start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_stats)
  
  for(cur_dist in in_dists) {
    for(cur_val in sd_vals) {
      
      target_abd <- GetTargetAbd(pars = cur_pars, choice = cur_dist,
                                 value = ifelse(cur_dist == "Carrying Capacities",
                                                0.5 * cur_val, cur_val))
      
      cur_pars$B <- GetNonEqualFeasibleB(target_abd, cur_pars)
      
      J <- BuildJacobian(eq_abd = target_abd, pars = cur_pars)
      
      cur_eig <- GetEig(J)
      cur_stable <- (cur_eig < 0)
      
      feas_data <- data.frame(Stable = cur_stable, SAD = cur_dist, SD = cur_val)
      
      cur_data <- cbind(cur_stats, feas_data)
      
      out_data <- rbind(out_data, cur_data)
    }
  }
}

print(Sys.time() - start_time)

# writing out the data
filename <- "RealisticSADs"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)
