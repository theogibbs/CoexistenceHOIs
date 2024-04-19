########################################################
# Non-equal abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

### Simulation script

S <- 40
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 0
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

in_dists <- "Carrying Capacities"
sd_vals <- seq(0.01, 0.95, length.out = 10)

out_data <- data.frame()

start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_stats)
  print(cur_row)
  
  for(cur_dist in in_dists) {
    for(cur_val in sd_vals) {
      target_abd <- GetTargetAbd(pars = cur_pars, choice = cur_dist, value = cur_val)
      
      #cur_pars$B <- GetNonEqualFeasibleB(target_abd, cur_pars)
      cur_pars$B <- GetEqualAcrossRowsFeasibleB(target_abd, cur_pars)
      mean_B <- mean(cur_pars$B[cur_pars$B != 0])
      
      J <- BuildJacobian(eq_abd = target_abd, pars = cur_pars)
      
      cur_eigs <- eigen(J, only.values = TRUE)$values
      re_eigs <- Re(cur_eigs)
      im_eigs <- Im(cur_eigs)
      
      cur_eig <- GetEig(J)
      cur_stable <- (cur_eig < 0)
      
      feas_data <- data.frame(Stable = cur_stable, SAD = cur_dist, SD = cur_val,
                              Real = re_eigs, Imaginary = im_eigs, MeanB = mean_B)
      
      cur_data <- cbind(cur_stats, feas_data)
      
      out_data <- rbind(out_data, cur_data)
    }
  }
}

print(Sys.time() - start_time)

# writing out the data
filename <- "VaryingMeanZeroTargetAbundances"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)
