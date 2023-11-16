########################################################
# Diverse community simulations for multiple equilibria
# with mean and variance independetly changing
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
mu_A <- 0
sigma_A <- 1
mu_B <- 0
sigma_B <- 0
rho_A <- 0
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

pars <- BuildPars(in_pars[1,])

target_abd <- GetTargetAbd(pars)
pars$B <- GetFeasibleB(target_abd, pars)
BuildJacobian(target_abd, pars)
pars$r - pars$A %*% target_abd - pars$B %*% as.vector(outer(target_abd, target_abd))
Bmat <- pars
Bmat$A <- matrix(0, nrow = Bmat$S, ncol = Bmat$S)
Bmat <- BuildJacobian(target_abd, Bmat)
Bmat
rowSums(Bmat) / 2
rowSums(pars$A) - 1

const_mat <- BuildConstraintMat(pars)

num_repl <- 1
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

out_data <- data.frame()

ps <- seq(0, 0.5, length.out = 1)

start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_stats)
  for(p in ps) {
    
    target_abd <- GetTargetAbd(pars = cur_pars)
    const_vec <- BuildConstraintVector(target_abd, cur_pars)
    const_B <- GetCorrelatedB(pars = cur_pars, const_mat = const_mat,
                              const_vec = const_vec, p = p)
    cur_pars$B <- const_B
    
    J <- BuildJacobian(eq_abd = target_abd, pars = cur_pars)
    
    cur_eig <- GetEig(J)
    
    cur_stable <- (cur_eig < 0)
    
    cur_corr <- GetCorr(cur_pars)
    
    cur_data <- data.frame(Stable = cur_stable, p = p, Corr = cur_corr)
    cur_data <- cbind(cur_stats, cur_data)
    
    out_data <- rbind(out_data, cur_data)
  }
}

print(Sys.time() - start_time)

# writing out the data
filename <- "Correlations"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


