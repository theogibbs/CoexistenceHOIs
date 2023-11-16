########################################################
# Diverse community simulations minimizing the
# eigenvalues as a function of the target abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

S <- 5
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(-0.5, 0, 0.5) / S
sigma_A <- seq(0, 0.25, length.out = 50)
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

num_repl <- 5
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

out_data <- data.frame()

start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  
  cur_pars <- BuildPars(cur_stats)
  
  cur_min <- optimize(f = GetMinEig, interval = c(0, 10), cur_pars)
  
  hoi_eig <- cur_min$objective
  hoi_stable <- (hoi_eig < 0)
  
  min_data <- data.frame(TargetAbd = cur_min$minimum,
                         Eigenvalue = hoi_eig,
                         Stable = hoi_stable,
                         Type = "Minimized")
  
  cur_data <- cbind(cur_stats, min_data)
  out_data <- rbind(out_data, cur_data)
}

print(Sys.time() - start_time)

# writing out the data
filename <- "DiverseMinimaTargets"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


