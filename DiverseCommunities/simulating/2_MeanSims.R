########################################################
# Diverse community mean correlation simulations for
# panel (A) of Fig. 2
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
mu_A <- seq(-0.75, 0.25, length.out = 10)
sigma_A <- seq(0, 0.1, length.out = 2)
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

out_data <- data.frame()

num_repl <- 100 # choosing the number of replicates
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

start_time <- Sys.time()

# looping over all the parameters
for(cur_row in 1:nrow(in_pars)) {
  
  # extracting the current parameters
  cur_params <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_params)
  
  # generating the target abundances and the constrained HOIs
  target_abd <- GetTargetAbd(cur_pars)
  cur_pars$B <- GetFeasibleB(target_abd, cur_pars)
  
  # finding the eigenvalues
  hoi_abds <- target_abd
  hoi_eig <- GetEig(BuildJacobian(hoi_abds, cur_pars))
  
  # recording feasibility (which is always true in these sims)
  # and stability
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
filename <- "DiverseCommunityMean"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


