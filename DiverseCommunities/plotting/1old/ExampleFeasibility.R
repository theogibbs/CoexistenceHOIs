########################################################
# Example diverse community dynamics
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

set.seed(1)

S <- 10
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 0.1
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
target_abd <- GetTargetAbd(pars)
ini_state <- runif(S, min = 0, max = 0.05)

end_time <- 50
time_step <- 0.1

out_fs <- IntegrateDynamics(inistate = ini_state,
                            pars = pars,
                            endtime = end_time,
                            timestep = time_step,
                            fn = Dynamics)
out_fs <- melt(out_fs, id.vars = "time")
out_fs$Type <- "Feasible"

pars$A <- pars$A * 5
diag(pars$A) <- 1

out_ufs <- IntegrateDynamics(inistate = ini_state,
                             pars = pars,
                             endtime = end_time,
                             timestep = time_step,
                             fn = Dynamics)
out_ufs <- melt(out_ufs, id.vars = "time")
out_ufs$Type <- "Not feasible"

series <- rbind(out_fs, out_ufs)

plExample <- ggplot() +
  geom_line(data = series, aes(x = time, y = value, color = variable),
            size = 2, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  facet_wrap(~Type) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = 0, linetype = "dashed")
show(plExample)

jpeg("../CoexistenceHOIs-Paper/figs/FeasibilityExample.jpeg",
     width = 1800, height = 1000, res = 300)
plExample
dev.off()

