########################################################
# Diverse community simulations for multiple equilibria
# with mean and variance both changing for SI Fig. J
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
mu_A <- seq(19, 21, length.out = 10) / S
sigma_A <- 0
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
num_repl <- 100
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

target_abd <- GetTargetAbd(BuildPars(in_pars[1,]))

out_abds <- data.frame()
ini_state <- runif(S, min = 0, max = 0.01)
#ini_state <- rep(target_abd, times = S) - runif(S, min = 0, max = 0.0000001)
end_time <- 1e4
zero_cutoff <- 1e-7

start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  
  cur_pars <- in_pars[cur_row,]
  pars <- BuildPars(cur_pars)
  target_abd <- GetTargetAbd(pars)
  cur_r <- pars$r
  
  cur_pars$Type <- "Feasible Pairwise"
  pw_r <- pars$A %*% rep(target_abd, times = pars$S)
  pars$r <- pw_r
  pwi_abds <- GetAbds(pars = pars,
                      ini_state = ini_state,
                      end_time = end_time,
                      zero_cutoff = zero_cutoff)
  
  pwi_abds <- cbind(cur_pars, pwi_abds)
  
  pars$r <- cur_r
  
  cur_pars$Type <- "Constrained HOIs"
  const_B <- GetFeasibleB(target_abd, pars)
  pars$B <- const_B
  
  hoi_abds <- GetAbds(pars = pars,
                      ini_state = ini_state,
                      end_time = end_time,
                      zero_cutoff = zero_cutoff)
  
  hoi_abds <- cbind(cur_pars, hoi_abds)
  
  #cur_pars$Type <- "Permuted HOIs"
  #non_zero_vals <- const_B[const_B != 0]
  #perm_vals <- sample(non_zero_vals, size = length(non_zero_vals))
  
  #perm_B <- const_B
  #perm_B[perm_B != 0] <- perm_vals
  #pars$B <- perm_B
  
  #per_abds <- GetAbds(pars = pars,
  #                    ini_state = ini_state,
  #                    end_time = end_time,
  #                    zero_cutoff = zero_cutoff)
  
  #per_abds <- cbind(cur_pars, per_abds)
  
  cur_abds <- rbind(pwi_abds, hoi_abds)#, per_abds)

  out_abds <- rbind(out_abds, cur_abds)
}

print(Sys.time() - start_time)

# writing out the data
filename <- "MultipleEquilibriaMean"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_abds, file = cur_file)

