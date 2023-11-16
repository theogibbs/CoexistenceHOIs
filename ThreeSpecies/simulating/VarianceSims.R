########################################################
# Three species simulations for varying mean
# interaction strength
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Simulating script

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(0.4, 0, -0.4)
sigma_A <- seq(0.25, 2, length.out = 10)
rho_A <- 0
b <- 0

in_pars <- crossing(S = S,
                    MuR = mu_R,
                    SigmaR = sigma_R,
                    MuD = mu_D,
                    SigmaD = sigma_D,
                    MuA = mu_A,
                    SigmaA = sigma_A,
                    RhoA = rho_A,
                    b = b)

num_repl <- 1e3
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

start_time <- Sys.time()

out_data <- data.frame()

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  cur_pars <- BuildThreeSpeciesPars(cur_stats)
  
  target_abd <- GetTargetAbd(cur_pars)
  
  cur_pars$b <- GetThreeSpeciesB(target_abd, cur_pars)
  
  hoi_abds <- target_abd
  hoi_eig <- GetEig(BuildThreeSpeciesJacobian(hoi_abds, cur_pars))
  
  hoi_feas <- prod(target_abd > 0)
  hoi_stable <- (hoi_eig < 0)
  
  hoi_data <- data.frame(Feasible = hoi_feas,
                         Stable = hoi_stable,
                         Type = "Constrained HOIs")
  
  hoi_data <- cbind(cur_stats, hoi_data)
  out_data <- rbind(out_data, hoi_data)
  
}

print(Sys.time() - start_time)

# writing out the data
filename <- "ThreeSpeciesVarianceSims"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)
