########################################################
# Example dynamics solving and plotting
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 2
sigma_D <- 0
mu_A <- c(1.0, 1.05, 1.1, 1.5, 2.1, 2.5)
sigma_A <- 0
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

out_data <- data.frame()
ini_state <- runif(S, min = 0, max = 0.1)

for(cur_row in 1:nrow(in_pars)) {
  
  pars <- BuildThreeSpeciesPars(in_pars[cur_row,])
  pars$A[upper.tri(pars$A)] <- 1 / pars$A[upper.tri(pars$A)]
  
  out_pwi <- IntegrateDynamics(inistate = ini_state,
                               pars = pars,
                               endtime = 5e2,
                               timestep = 0.1,
                               fn = ThreeSpeciesDynamics)
  out_pwi$MuA <- in_pars$MuA[cur_row]
  out_pwi$Type = "No HOIs"
  
  target_abd <- GetTargetAbd(pars)
  opt_b <- GetThreeSpeciesB(target_abd, pars)
  pars$b <- opt_b
  
  out_hoi <- IntegrateDynamics(inistate = ini_state,
                               pars = pars,
                               endtime = 5e2,
                               timestep = 0.1,
                               fn = ThreeSpeciesDynamics)
  out_hoi$MuA <- in_pars$MuA[cur_row]
  out_hoi$Type = "Constrained HOIs"
  
  cur_data <- rbind(out_pwi, out_hoi)
  
  out_data <- rbind(out_data, cur_data)
}

out_data <- out_data %>%
  melt(id.vars = c("time", "MuA", "Type")) %>%
  mutate(Type = factor(Type, levels = c("No HOIs", "Constrained HOIs")))

plDominance <- ggplot() +
  geom_line(data = out_data, aes(x = time, y = value, color = variable),
            size = 2, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  facet_grid(Type~MuA, labeller = label_bquote(cols = a == .(MuA))) +
  theme_classic() +
  geom_hline(yintercept = target_abd, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(legend.position = "none",
        text = element_text(size=10),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15))
show(plDominance)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesCompetitiveDominance.jpeg", width = 2750, height = 1300, res = 300)
plDominance
dev.off()


