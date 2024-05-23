########################################################
# Finds the outcomes of simulations for three species
# dynamics to create panel (A) of SI Fig. J
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

set.seed(3)

# choosing parameters for the simulations

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 1
rho_A <- 0
b <- 0

# building all parameter combinations
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

# standardizing pairwise interaction statistics
off_diag_entries <- pars$A[upper.tri(pars$A)]
off_diag_entries <- c(off_diag_entries, pars$A[lower.tri(pars$A)])

pars$A <- pars$A / sd(off_diag_entries)
off_diag_entries <- pars$A[upper.tri(pars$A)]
off_diag_entries <- c(off_diag_entries, pars$A[lower.tri(pars$A)])

pars$A <- pars$A - mean(off_diag_entries)
diag(pars$A) <- 1

# choosing pairwise interaction statistics to loop over
in_mus <- seq(0.25, 1.5, length.out = 500)
in_sigmas <- seq(0, 0.75, length.out = 4)

# settings for integrating the dynamics
end_time <- 1e5
zero_cutoff <- 1e-7

# looping over all pairwise interaction mean and variances
out_abds <- data.frame()
for(cur_mu in in_mus) {
  for(cur_sigma in in_sigmas) {
    
    # setting up the current interaction statistics and parameterizing the model
    cur_pars <- pars
    cur_pars$A <- cur_sigma * cur_pars$A
    cur_pars$A <- cur_pars$A + cur_mu
    
    diag(cur_pars$A) <- mu_D
    
    target_abd <- GetTargetAbd(cur_pars)
    cur_pars$b <- GetThreeSpeciesB(target_abd, cur_pars)
    
    # choosing initial conditions
    ini_state <- target_abd - runif(S, 0, 0.01)
    
    # integrating the dynamics
    cur_abds <- GetThreeSpeciesAbds(pars = cur_pars,
                                    ini_state = ini_state,
                                    end_time = end_time,
                                    zero_cutoff = zero_cutoff)
    
    # recording simulation data
    cur_abds$SpID <- as.factor(1:3)
    cur_abds$MuA <- cur_mu
    cur_abds$SigmaA <- cur_sigma
    
    out_abds <- rbind(out_abds, cur_abds)
  }
}

# filtering out unstable or non-equilibrium solutions
test_df <- out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]
print(test_df)

# plotting the abundances as a function of the mean interaction strength
plAbds <- ggplot(out_abds, aes(x = MuA, y = Abd, color = SpID)) +
  geom_line(alpha = 0.75, size = 1) +
  theme_classic() +
  labs(x = expression("Pairwise interaction strength"~ (mu[A])),
       y = "Equilibrium Abundance", color = "") +
  theme(text = element_text(size=15),
        legend.position = "none",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(~SigmaA, labeller = label_bquote(cols = sigma[A] == .(SigmaA))) +
  geom_hline(yintercept = GetTargetAbd(pars), linetype = "dashed")
plAbds

out_abds <- out_abds %>%
  filter(SigmaA == in_sigmas[3])

# same as above but focusing on a specific variance in abundances
plBifurcation <- ggplot(out_abds, aes(x = MuA, y = Abd, color = SpID)) +
  geom_line(alpha = 0.75, size = 1) +
  theme_classic() +
  labs(x = expression("Mean strength of pairwise interactions"~ (mu[A])),
       y = "Equilibrium Abundance", color = "") +
  theme(text = element_text(size=13),
        legend.position = "none",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  geom_hline(yintercept = GetTargetAbd(pars), linetype = "dashed")
plBifurcation

# writing out the data
filename <- "ThreeSpeciesBifurcation"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_abds, file = cur_file)

