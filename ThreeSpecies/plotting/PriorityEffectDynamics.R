########################################################
# Plotting symmetric priority effects dynamics for
# three species
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Scripting

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(0.25, 0.75, 1, 1.5, 2.1)
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
  target_abd <- GetTargetAbd(pars)
  
  opt_b <- GetThreeSpeciesB(target_abd, pars)
  pars$b <- opt_b
  
  out_hoi <- IntegrateDynamics(inistate = ini_state,
                               pars = pars,
                               endtime = 50,
                               timestep = 0.1,
                               fn = ThreeSpeciesDynamics)
  out_hoi <- melt(out_hoi, id.vars = "time")
  out_hoi$MuA <- in_pars$MuA[cur_row]
  
  out_data <- rbind(out_data, out_hoi)
}

plPriority <- ggplot() +
  geom_line(data = out_data, aes(x = time, y = value, color = variable),
            size = 2, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  facet_grid(~MuA, labeller = label_bquote(cols = a == .(MuA))) +
  theme_classic() +
  geom_hline(yintercept = target_abd, linetype = "dashed") +
  theme(legend.position = "none",
        text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15))
show(plPriority)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesPriority.jpeg", width = 2400, height = 700, res = 300)
plPriority
dev.off()


