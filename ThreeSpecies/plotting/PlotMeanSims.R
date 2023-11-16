########################################################
# Plotting three species simulations for varying mean
# interaction strength
########################################################

out_data <- read_csv2("simdata/ThreeSpeciesMeanSims.csv")

melt_data <- out_data %>%
  mutate(Type = factor(Type, levels = c("No HOIs", "Constrained HOIs"))) %>%
  dplyr::group_by(MuA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("MuA", "SigmaA", "Type"))

plMean <- ggplot(melt_data, aes(x = MuA, y = value, color = as.factor(SigmaA))) +
  geom_line(size = 1, alpha = 0.5) +
  #geom_point(alpha = 0.75, size = 3) +
  theme_classic() +
  labs(x = expression("Pairwise interaction strength"~ (mu[A])),
       y = "Probability of Stability", color = expression(sigma[A])) +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5))
show(plMean)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesMean.jpeg", width = 1800, height = 1200, res = 300)
plMean
dev.off()

