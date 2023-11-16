########################################################
# Example diverse community dynamics
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

#set.seed(1)

S <- 10
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0.2
sigma_A <- abs(mu_A) / sqrt(3)
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

ini_state <- runif(S, min = 0, max = 0.05)
end_time_1 <- 200
end_time_2 <- 200
end_time_3 <- 200
time_step <- 0.1

out_ta <- IntegrateDynamics(inistate = ini_state,
                             pars = pars,
                             endtime = end_time_1,
                             timestep = time_step,
                             fn = Dynamics)

end_ta <- as.numeric(out_ta[nrow(out_ta), -1])
ini_rm <- end_ta
ini_rm[1] <- 0

out_rm <- IntegrateDynamics(inistate = ini_rm,
                            pars = pars,
                            endtime = end_time_2,
                            timestep = time_step,
                            fn = Dynamics)

out_rm$time <- out_rm$time + max(out_ta$time)

end_rm <- as.numeric(out_rm[nrow(out_rm), -1])
ini_in <- end_rm
ini_in[1] <- 1e-5

out_in <- IntegrateDynamics(inistate = ini_in,
                            pars = pars,
                            endtime = end_time_3,
                            timestep = time_step,
                            fn = Dynamics)

out_in$time <- out_in$time + max(out_rm$time)

out <- rbind(out_ta, out_rm, out_in)

melt_data <- melt(out, id.vars = "time")

plExample <- ggplot() +
  geom_line(data = melt_data, aes(x = time, y = value, color = variable),
            size = 2, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_vline(xintercept = c(end_time_1, end_time_1 + end_time_2), linetype = "dashed")
show(plExample)

jpeg("../CoexistenceHOIs-Paper/figs/DiverseExampleInvasions.jpeg",
     width = 2000, height = 1000, res = 300)
plExample
dev.off()

