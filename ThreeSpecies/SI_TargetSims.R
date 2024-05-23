########################################################
# Three species simulations for varying target
# abundances for SI Fig. E
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Simulating script

set.seed(3)

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(-0.25, 0, 0.25)
sigma_A <- 0.25
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

num_repl <- 1
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

out_data <- data.frame()

target_abds <- seq(0.5, 2, length.out = 1e3)

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 1
rho_A <- 0
b <- 0

in_stat <- crossing(S = S,
                    MuR = mu_R,
                    SigmaR = sigma_R,
                    MuD = mu_D,
                    SigmaD = sigma_D,
                    MuA = mu_A,
                    SigmaA = sigma_A,
                    RhoA = rho_A,
                    b = b)

pars <- BuildThreeSpeciesPars(in_stat)

for(cur_row in 1:nrow(in_pars)) {
  for(cur_target in target_abds) {
    cur_stats <- in_pars[cur_row,]
    
    cur_sigma <- cur_stats$SigmaA
    cur_mu <- cur_stats$MuA
    
    cur_A <- pars$A
    cur_A <- cur_sigma * cur_A
    cur_A <- cur_A + cur_mu
    diag(cur_A) <- diag(pars$A)
    
    cur_pars <- pars
    cur_pars$A <- cur_A
    
    target_abd <- GetTargetAbd(pars, choice = "Specified", value = cur_target)
    
    cur_pars$b <- GetThreeSpeciesB(target_abd, cur_pars)
    
    hoi_abds <- target_abd
    hoi_eig <- GetEig(BuildThreeSpeciesJacobian(hoi_abds, cur_pars))
    
    hoi_feas <- prod(target_abd > 0)
    hoi_stable <- (hoi_eig < 0)
    
    cur_data <- data.frame(TargetAbd = cur_target,
                           Feasible = hoi_feas,
                           Eigenvalue = hoi_eig,
                           Stable = hoi_stable,
                           Type = "Constrained HOIs")
    cur_data <- cbind(cur_stats, cur_data)
    
    out_data <- rbind(out_data, cur_data)
    
  }
}

# writing out the data
filename <- "ThreeSpeciesTargets"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


