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
mu_D <- 1
sigma_D <- 0
mu_A <- c(0.1, 0.25, 0.75, 1, 1.5)
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
in_pars$SigmaA <- mu_A / sqrt(3)

out_data <- data.frame()
ini_abds <- seq(0.01, 1.5, length.out = 20)

end_time <- 1e3
zero_cutoff <- 1e-7

for(cur_row in 1:nrow(in_pars)) {
  pw_pars <- BuildThreeSpeciesPars(in_pars[cur_row,])
  
  target_abd <- GetTargetAbd(pw_pars)
  const_b <- GetThreeSpeciesB(target_abd, pw_pars)
  hoi_pars <- pw_pars
  hoi_pars$b <- const_b
  
  for(cur_ini in ini_abds) {
    
    ini_state <- rep(cur_ini, times = 3) + runif(3, min = -cur_ini, max = cur_ini) * 0.01
    
    out_pwi <- GetThreeSpeciesAbds(pars = pw_pars,
                                   ini_state = ini_state,
                                   end_time = end_time,
                                   zero_cutoff = zero_cutoff)
    
    out_pwi <- cbind(in_pars[cur_row,], out_pwi)
    out_pwi$Type <- "No HOIs"
    out_pwi$IniAbd <- cur_ini
    out_pwi$SpID <- as.factor(1:3)
    
    out_hoi <- GetThreeSpeciesAbds(pars = hoi_pars,
                                   ini_state = ini_state,
                                   end_time = end_time,
                                   zero_cutoff = zero_cutoff)
    
    out_hoi <- cbind(in_pars[cur_row,], out_hoi)
    out_hoi$Type <- "Constrained HOIs"
    out_hoi$IniAbd <- cur_ini
    out_hoi$SpID <- as.factor(1:3)
    
    cur_data <- rbind(out_pwi, out_hoi)
    
    out_data <- rbind(out_data, cur_data)
  }
}

melt_data <- out_data %>%
  mutate(Abd = ifelse(!Equilibrium, NA, Abd)) %>%
  mutate(Type = factor(Type, levels = c("No HOIs", "Constrained HOIs")))

plIC <- ggplot() +
  geom_line(data = melt_data, aes(x = IniAbd, y = Abd, color = SpID),
            size = 2, alpha = 0.5) +
  xlab("Initial Abundance") + ylab("Equilibrium Abundance") +
  facet_grid(Type~MuA, labeller = label_bquote(cols = a == .(MuA))) +
  theme_classic() +
  geom_hline(yintercept = target_abd, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(legend.position = "none",
        text = element_text(size=10),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15))
show(plIC)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesInitialConditions.jpeg", width = 2750, height = 1300, res = 300)
plIC
dev.off()


