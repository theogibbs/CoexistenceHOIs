########################################################
# Diverse community alternative sampling
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

S <- 20
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(-0.05, 0, 0.05)
sigma_A <- c(0.01, 0.15)
rho_A <- 0
mu_B <- 0
sigma_B <- seq(0, 0.1, length.out = 10)
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

out_data <- data.frame()

num_repl <- 50
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  
  cur_params <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_params)
  
  target_abd <- GetTargetAbd(cur_pars)
  
  cur_pars$B <- GetFeasibleB(target_abd, cur_pars)
  
  cur_pert <- matrix(rnorm(cur_pars$S * choose(cur_pars$S - 1, 2), mean = 0, sd = cur_params$SigmaB),
                     nrow = cur_pars$S, ncol = choose(cur_pars$S - 1, 2))
  row_sum_mat <- matrix(rowSums(cur_pert), nrow = nrow(cur_pert), ncol = ncol(cur_pert))
  row_sum_mat <- row_sum_mat / choose(cur_pars$S - 1, 2)
  
  cur_pert <- cur_pert - row_sum_mat
  
  cur_pars$B[cur_pars$B != 0] <- cur_pars$B[cur_pars$B != 0] + cur_pert
  
  hoi_abds <- target_abd
  hoi_eig <- GetEig(BuildJacobian(hoi_abds, cur_pars))
  
  hoi_feas <- prod(target_abd > 0)
  hoi_stable <- (hoi_eig < 0)
  
  hoi_data <- data.frame(Feasible = hoi_feas,
                         Stable = hoi_stable,
                         Type = "Constrained HOIs")
  hoi_data <- cbind(cur_params, hoi_data)
  
  out_data <- rbind(out_data, hoi_data)
  
}

print(Sys.time() - start_time)

# writing out the data
filename <- "DiverseCommunityAlternativeSampling"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


