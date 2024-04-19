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
mu_A <- 0.05
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

pars <- BuildPars(in_pars)
target_abd <- GetTargetAbd(pars, value = 0.1)
ini_state <- runif(S, min = 0, max = 0.05)
#ini_state <- target_abd + rnorm(pars$S, sd = 0.01)

end_time <- 50
time_step <- 0.1

out_pw <- IntegrateDynamics(inistate = ini_state,
                            pars = pars,
                            endtime = end_time,
                            timestep = time_step,
                            fn = Dynamics)
out_pw <- melt(out_pw, id.vars = "time")
out_pw$Type <- "No HOIs"

const_B <- GetFeasibleB(target_abd, pars)
#const_B <- GetNonEqualFeasibleB(target_abd, pars)
#const_B <- GetAllTermsB(target_abd, pars)
const_B <- GetEqualAcrossRowsFeasibleB(target_abd, pars)
#const_B <- GetIdCorrelatedB(target_abd, pars, p = 1)
const_B <- GetAbdVarianceFeasibleB(target_abd, pars)

pars$B <- const_B

out_hoi <- IntegrateDynamics(inistate = ini_state,
                             pars = pars,
                             endtime = end_time,
                             timestep = time_step,
                             fn = Dynamics)
out_hoi <- melt(out_hoi, id.vars = "time")
out_hoi$Type <- "Constrained HOIs"

series <- rbind(out_pw, out_hoi) %>%
  mutate(Type = factor(Type, levels = c("No HOIs", "Constrained HOIs")))

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
  ggtitle("B")
show(plExample)

jpeg("../CoexistenceHOIs-Paper/figs/Fig1ConstraintExample.jpeg",
     width = 2000, height = 1100, res = 300)
plExample
dev.off()

