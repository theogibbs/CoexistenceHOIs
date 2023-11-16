########################################################
# Three species simulations for multiple equilibria
# with mean and variance independetly changing
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
mu_A <- c(-0.25, 0, 0.25)
sigma_A <- seq(0, 1, length.out = 100)
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

num_repl <- 500
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

target_abd <- GetTargetAbd(BuildThreeSpeciesPars(in_pars[1,]))

out_abds <- data.frame()
ini_state <- runif(S, min = 0, max = 0.01)
ini_state <- rep(target_abd, times = 3) - runif(S, min = 0, max = 0.0000001)
end_time <- 1e4
zero_cutoff <- 1e-7

compareable_params <- FALSE
start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  
  cur_pars <- in_pars[cur_row,]
  pars <- BuildThreeSpeciesPars(cur_pars)
  
  if(compareable_params) {
    cur_pars$Type <- "No HOIs"
    pwi_abds <- GetThreeSpeciesAbds(pars = pars,
                                    ini_state = ini_state,
                                    end_time = end_time,
                                    zero_cutoff = zero_cutoff)
    
    pwi_abds <- cbind(cur_pars, pwi_abds)
  }
  
  cur_pars$Type <- "Constrained HOIs"
  target_abd <- GetTargetAbd(pars)
  const_b <- GetThreeSpeciesB(target_abd, pars)
  pars$b <- const_b
  
  hoi_abds <- GetThreeSpeciesAbds(pars = pars,
                                  ini_state = ini_state,
                                  end_time = end_time,
                                  zero_cutoff = zero_cutoff)
  
  hoi_abds <- cbind(cur_pars, hoi_abds)
  
  if(compareable_params) {
    cur_abds <- rbind(pwi_abds, hoi_abds)
  } else {
    cur_abds <- hoi_abds
  }
  out_abds <- rbind(out_abds, cur_abds)
}

print(Sys.time() - start_time)

# writing out the data
filename <- "ThreeSpeciesEquilibriaMV_test"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_abds, file = cur_file)

