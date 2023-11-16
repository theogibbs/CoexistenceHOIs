########################################################
# Example dynamics solving and plotting
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

#set.seed(123)

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- -0.25
sigma_A <- 1
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

pars <- BuildThreeSpeciesPars(in_pars)
target_abd <- GetTargetAbd(pars)
ini_state <- target_abd + rnorm(S, sd = 0.5)

out_pw <- IntegrateDynamics(inistate = ini_state,
                            pars = pars,
                            endtime = 50,
                            timestep = 0.1,
                            fn = ThreeSpeciesDynamics)
out_pw <- melt(out_pw, id.vars = "time")
out_pw$Type <- "No HOIs"

const_b <- GetThreeSpeciesB(rep(target_abd, times = length(pars$r)), pars)
pars$b <- const_b

out_hoi <- IntegrateDynamics(inistate = ini_state,
                             pars = pars,
                             endtime = 50,
                             timestep = 0.1,
                             fn = ThreeSpeciesDynamics)
out_hoi <- melt(out_hoi, id.vars = "time")
out_hoi$Type <- "Constrained HOIs"

series <- rbind(out_pw, out_hoi) %>%
  mutate(Type = factor(Type, levels = c("No HOIs", "Constrained HOIs")))

plExample <- ggplot() +
  geom_line(data = series, aes(x = time, y = value, color = variable),
            size = 3, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  facet_wrap(~Type) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")
show(plExample)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesExample.jpeg", width = 2750, height = 1500, res = 300)
plExample
dev.off()