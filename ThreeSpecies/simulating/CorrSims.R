########################################################
# Three species simulations for varying correlations in
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
mu_A <- 0
sigma_A <- seq(0.01, 2, length.out = 20)
rho_A <- c(-1, 0, 1)
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

num_repl <- 500
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

out_data <- data.frame()

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  cur_pars <- BuildThreeSpeciesPars(cur_stats)
  
  target_abd <- GetTargetAbd(cur_pars)
  
  pw_abds <- target_abd
  pw_eig <- GetEig(BuildThreeSpeciesJacobian(pw_abds, cur_pars))
  
  pw_feas <- prod(pw_abds > 0)
  pw_stable <- ifelse(pw_feas, (pw_eig < 0), 0)
  
  feas_pw_data <- data.frame(SigmaA = cur_stats$SigmaA,
                             RhoA = cur_stats$RhoA,
                             Feasible = pw_feas,
                             Stable = pw_stable,
                             Type = "Feasible Pairwise")
  
  cur_pars$b <- GetThreeSpeciesB(target_abd, cur_pars)
  
  hoi_abds <- target_abd
  hoi_eig <- GetEig(BuildThreeSpeciesJacobian(hoi_abds, cur_pars))
  
  hoi_feas <- prod(target_abd > 0)
  hoi_stable <- (hoi_eig < 0)
  
  hoi_data <- data.frame(SigmaA = cur_stats$SigmaA,
                         RhoA = cur_stats$RhoA,
                         Feasible = hoi_feas,
                         Stable = hoi_stable,
                         Type = "Constrained HOIs")
  
  cur_data <- rbind(feas_pw_data, hoi_data)
  out_data <- rbind(out_data, cur_data)
  
}

# writing out the data
filename <- "ThreeSpeciesCorrSims"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)
