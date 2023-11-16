########################################################
# Plotting three species simulations for varying
# target abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

set.seed(4)

S <- 10
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 0.3
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
off_diag_entries <- pars$A[upper.tri(pars$A)]
off_diag_entries <- c(off_diag_entries, pars$A[lower.tri(pars$A)])

pars$A <- pars$A - mean(off_diag_entries)
diag(pars$A) <- 1

target_abd <- GetTargetAbd(pars)
pars$B <- GetFeasibleB(target_abd, pars)

in_mus <- seq(0, 0.5, length.out = 500)

end_time <- 1e5
zero_cutoff <- 1e-7

out_abds <- data.frame()

for(cur_mu in in_mus) {
  cur_pars <- pars
  cur_pars$A <- cur_pars$A + cur_mu
  
  diag(cur_pars$A) <- mu_D
  
  cur_pars$B <- pars$B
  non_zero_B <- (pars$B != 0)
  num_non_zero <- unique(rowSums(non_zero_B))
  
  cur_pars$B[non_zero_B] <- cur_pars$B[non_zero_B] - cur_mu * (cur_pars$S-1) / target_abd / num_non_zero
  # NEED TO ADJUST FOR THE NUMBER OF NON ZERO B ENTRIES!
  ini_state <- rep(target_abd, times = cur_pars$S) - runif(cur_pars$S, min = 0, max = 0.1)
  
  cur_abds <- GetAbds(pars = cur_pars,
                      ini_state = ini_state,
                      end_time = end_time,
                      zero_cutoff = zero_cutoff)
  
  cur_abds$SpID <- as.factor(1:S)
  cur_abds$MuA <- cur_mu

  out_abds <- rbind(out_abds, cur_abds)
}

test_df <- out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]
print(test_df)
filt_abds <- out_abds %>%
  filter(Equilibrium) %>%
  filter(Stability)

plBifurcation <- ggplot(filt_abds, aes(x = MuA, y = Abd, color = SpID)) +
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
  geom_hline(yintercept = GetTargetAbd(pars), linetype = "dashed") +
  ggtitle(paste(S, "species"))
plBifurcation

jpeg("../CoexistenceHOIs-Paper/figs/Fig6DiverseBifurcation.jpeg",
     width = 1300, height = 1000, res = 300)
plBifurcation
dev.off()

