########################################################
# Example stability diverse community dynamics
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

set.seed(5)

S <- 10
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 0.25
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

const_B <- GetFeasibleB(target_abd, pars)
pars$B <- const_B

ini_state <- GetTargetAbd(pars)
end_time_1 <- 5
end_time_2 <- 10
time_step <- 0.1

out_eq <- IntegrateDynamics(inistate = ini_state,
                            pars = pars,
                            endtime = end_time_1,
                            timestep = time_step,
                            fn = Dynamics)

end_eq <- as.numeric(out_eq[nrow(out_eq), -1])
ini_pert <- end_eq
ini_pert <- ini_pert + runif(pars$S, min = -0.25, max = 0.25)

out_pert <- IntegrateDynamics(inistate = ini_pert,
                              pars = pars,
                              endtime = end_time_2,
                              timestep = time_step,
                              fn = Dynamics)

out_pert$time <- out_pert$time + max(out_eq$time)

out_stable <- rbind(out_eq, out_pert)
out_stable$Stability <- "Stable"

sigma_A <- 0.75

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

const_B <- GetFeasibleB(target_abd, pars)
pars$B <- const_B



out_eq <- IntegrateDynamics(inistate = ini_state,
                            pars = pars,
                            endtime = end_time_1,
                            timestep = time_step,
                            fn = Dynamics)

end_eq <- as.numeric(out_eq[nrow(out_eq), -1])
ini_pert <- end_eq
ini_pert <- ini_pert + runif(pars$S, min = -0.25, max = 0.25)

out_pert <- IntegrateDynamics(inistate = ini_pert,
                              pars = pars,
                              endtime = end_time_2,
                              timestep = time_step,
                              fn = Dynamics)

out_pert$time <- out_pert$time + max(out_eq$time)

out_unstable <- rbind(out_eq, out_pert)
out_unstable$Stability <- "Unstable"

out <- rbind(out_stable, out_unstable)
melt_data <- melt(out, id.vars = c("time", "Stability")) %>%
  mutate(Stability = factor(Stability, levels = c("Stable", "Unstable")))

plExample <- ggplot() +
  geom_line(data = melt_data, aes(x = time, y = value, color = variable),
            size = 2, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  theme_classic() +
  facet_wrap(~Stability) +
  theme(legend.position = "none",
        text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_vline(xintercept = end_time_1, linetype = "dashed") +
  ggtitle("C")
show(plExample)

jpeg("../CoexistenceHOIs-Paper/figs/Fig1StabilityExample.jpeg",
     width = 2250, height = 1100, res = 300)
plExample
dev.off()

