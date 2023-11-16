########################################################
# Diverse community B tensor testing
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
mu_A <- 0
sigma_A <- 0.1
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
pars <- BuildPars(in_pars)
#pars
target_abd <- GetTargetAbd(pars = pars)

const_vec <- BuildConstraintVector(target_abd, pars)
const_mat <- BuildConstraintMat(pars, stab_const = FALSE)

const_B <- GetConstrainedB(pars = pars, const_mat = const_mat, const_vec = const_vec, beta = 0.5, p = 0)
const_mat %*% as.vector(t(const_B)) - const_vec

pars$B <- const_B

print("Should be zero:")
print((pars$r - rowSums(pars$A) - rowSums(pars$B)))

J <- BuildJacobian(eq_abd = rep(1, times = pars$S), pars)
GetEig(J)

pars$B <- GetFeasibleB(target_abd, pars)

J <- BuildJacobian(eq_abd = rep(1, times = pars$S), pars)
GetEig(J)
