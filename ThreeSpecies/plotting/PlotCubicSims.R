########################################################
# Plotting three species simulations for varying
# target abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

out_data <- read_csv2("simdata/ThreeSpeciesCubic.csv")

melt_data <- out_data %>%
  dplyr::group_by(MuR, MuA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("MuR", "MuA", "SigmaA", "Type"))

plCubic <- ggplot(melt_data, aes(x = SigmaA, y = value, color = Type)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(alpha = 0.75, size = 2) +
  theme_classic() +
  labs(x = expression("Pairwise interaction heterogeneity"~ (sigma[A])),
       y = "Probability of Stability", color = "") +
  theme(text = element_text(size=15),
        legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(MuR ~ MuA, labeller = label_bquote(rows = mu[R] == .(MuR),
                                                cols = mu[A] == .(MuA)))
show(plCubic)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesCubic.jpeg",
     width = 2800, height = 2200, res = 300)
plCubic
dev.off()

