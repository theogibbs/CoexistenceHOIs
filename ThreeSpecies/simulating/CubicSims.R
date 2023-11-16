########################################################
# Three species simulations cubic dynamics
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Simulating script

S <- 3
mu_R <- c(0.5, 1, 5)
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(-0.5, 0, 0.5)
sigma_A <- seq(0.01, 2, length.out = 20)
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

num_repl <- 50
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

out_data <- data.frame()

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_stats)
  
  target_abd <- GetTargetAbd(pars, choice = "Specified", value = 1)
  
  pw_abds <- target_abd
  pw_eig <- GetEig(BuildCubicJacobian(pw_abds, cur_pars))
  
  pw_feas <- prod(pw_abds > 0)
  pw_stable <- ifelse(pw_feas, (pw_eig < 0), 0)
  
  feas_pw_data <- data.frame(SigmaA = cur_stats$SigmaA,
                             MuA = cur_stats$MuA,
                             MuR = cur_stats$MuR,
                             Feasible = pw_feas,
                             Stable = pw_stable,
                             Type = "Feasible Pairwise")
  
  cur_pars$b <- GetThreeSpeciesB(target_abd, cur_pars)
  
  hoi_abds <- target_abd
  hoi_eig <- GetEig(BuildCubicJacobian(hoi_abds, cur_pars))
  
  hoi_feas <- prod(target_abd > 0)
  hoi_stable <- (hoi_eig < 0)
  
  hoi_data <- data.frame(SigmaA = cur_stats$SigmaA,
                         MuA = cur_stats$MuA,
                         MuR = cur_stats$MuR,
                         Feasible = hoi_feas,
                         Stable = hoi_stable,
                         Type = "Constrained HOIs")
  
  cur_data <- rbind(feas_pw_data, hoi_data)
  out_data <- rbind(out_data, cur_data)
  
}

# writing out the data
filename <- "ThreeSpeciesCubic"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


