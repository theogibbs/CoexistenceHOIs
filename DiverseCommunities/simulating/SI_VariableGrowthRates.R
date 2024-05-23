########################################################
# Diverse community variable growth rates for SI Fig. B
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

S <- 20
mu_R <- 1
sigma_R <- seq(0, 1.5, length.out = 10)
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- c(0, 0.25)
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
filename <- "DiverseCommunityGrowthRates"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


