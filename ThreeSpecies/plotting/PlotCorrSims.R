########################################################
# Plotting three species simulations for varying
# correlations in interaction strengths
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

out_data <- read_csv2("simdata/ThreeSpeciesCorrSims.csv")

melt_data <- out_data %>%
  dplyr::group_by(RhoA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("RhoA", "SigmaA", "Type"))

plCorr <- ggplot(melt_data, aes(x = SigmaA, y = value, color = Type)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(alpha = 0.75, size = 2) +
  theme_classic() +
  labs(x = expression("Pairwise interaction heterogeneity"~ (sigma[A])),
       y = "Probability of Stability", color = "") +
  theme(text = element_text(size=15),
        legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  facet_grid(~ RhoA,
             labeller = label_bquote(cols = rho[A] == .(RhoA)))
show(plCorr)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesCorr.jpeg", width = 2750, height = 1200, res = 300)
plCorr
dev.off()

