########################################################
# Three species simulations for multiple equilibria
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
mu_A <- seq(0.1, 1, length.out = 10)
sigma_A <- 0
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
in_pars$SigmaA <- in_pars$MuA / sqrt(3)

num_repl <- 100
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

out_abds <- data.frame()
ini_state <- runif(S, min = 0, max = 0.1)
end_time <- 1e4
zero_cutoff <- 1e-7

for(cur_row in 1:nrow(in_pars)) {
  
  cur_pars <- in_pars[cur_row,]
  pars <- BuildThreeSpeciesPars(cur_pars)
  cur_pars$Type <- "No HOIs"
  
  pwi_abds <- GetThreeSpeciesAbds(pars = pars,
                                  ini_state = ini_state,
                                  end_time = end_time,
                                  zero_cutoff = zero_cutoff)
  
  pwi_abds <- cbind(cur_pars, pwi_abds)
  
  cur_pars$Type <- "Constrained HOIs"
  target_abd <- GetTargetAbd(pars)
  opt_b <- GetThreeSpeciesB(target_abd, pars)
  pars$b <- opt_b
  
  hoi_abds <- GetThreeSpeciesAbds(pars = pars,
                                  ini_state = ini_state,
                                  end_time = end_time,
                                  zero_cutoff = zero_cutoff)
  
  hoi_abds <- cbind(cur_pars, hoi_abds)
  cur_abds <- rbind(pwi_abds, hoi_abds)
  
  out_abds <- rbind(out_abds, cur_abds)
}

# writing out the data
filename <- "ThreeSpeciesEquilibria"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_abds, file = cur_file)

