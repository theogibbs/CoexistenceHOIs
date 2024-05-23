########################################################
# Diversity versus variance simulations for SI Fig. D
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

S <- seq(10, 100, by = 10)
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 0.15
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

comp_pars <- in_pars
comp_pars$SigmaA <- 0.025
comp_pars$MuA <- 0.025

fac_pars <- comp_pars
fac_pars$MuA <- - comp_pars$MuA

in_pars <- rbind(in_pars, comp_pars, fac_pars)

out_data <- data.frame()

num_repl <- 10
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  
  cur_params <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_params)
  
  target_abd <- GetTargetAbd(cur_pars)
  target_abd <- rep(target_abd, times = cur_params$S)
  
  pw_abds <- target_abd
  pw_eig <- GetEig(BuildJacobian(pw_abds, cur_pars))
  
  pw_feas <- prod(pw_abds > 0)
  pw_stable <- ifelse(pw_feas, (pw_eig < 0), 0)
  
  pw_feas_data <- data.frame(Feasible = pw_feas,
                             Eigenvalue = pw_eig,
                             Stable = pw_stable,
                             Type = "Feasible Pairwise")
  
  pw_feas_data <- cbind(cur_params, pw_feas_data)
  
  cur_pars$B <- GetFeasibleB(target_abd, cur_pars)
  
  hoi_abds <- target_abd
  hoi_eig <- GetEig(BuildJacobian(hoi_abds, cur_pars))
  
  hoi_feas <- prod(target_abd > 0)
  hoi_stable <- (hoi_eig < 0)
  
  hoi_data <- data.frame(Feasible = hoi_feas,
                         Eigenvalue = hoi_eig,
                         Stable = hoi_stable,
                         Type = "Constrained HOIs")
  hoi_data <- cbind(cur_params, hoi_data)
  
  cur_data <- rbind(pw_feas_data, hoi_data)
  out_data <- rbind(out_data, cur_data)
  
}

print(Sys.time() - start_time)

# writing out the data
filename <- "DiverseCommunityDiversityVariance"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


