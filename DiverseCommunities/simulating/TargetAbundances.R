########################################################
# Varying the target abundances
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
pars <- BuildPars(in_pars)
#pars
target_abd <- GetTargetAbd(pars = pars)

const_vec <- BuildConstraintVector(target_abd, pars)
const_mat <- BuildConstraintMat(pars)

const_B <- GetConstrainedB(pars = pars, const_mat = const_mat, const_vec = const_vec, p = 0)
const_mat %*% as.vector(t(const_B)) - const_vec

pars$B <- const_B

print("Should be zero:")
print((pars$r - rowSums(pars$A) - rowSums(pars$B)))

J <- BuildJacobian(eq_abd = rep(1, times = pars$S), pars)
J

GetCorr(pars)

#####