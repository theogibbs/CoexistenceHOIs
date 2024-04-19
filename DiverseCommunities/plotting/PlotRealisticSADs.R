########################################################
# Plotting all of the target simulations
########################################################

# Loading dependencies

source("./Functions.R")

# Plotting script

out_data <- read_csv2("simdata/RealisticSADs.csv")

melt_data <- out_data %>%
  group_by(MuA, SigmaA, SAD, SD) %>%
  dplyr::summarize(ProbStable = mean(Stable)) %>%
  dplyr::mutate(SigmaA = ifelse(SigmaA == 0, "All Equal", "Variable")) %>%
  dplyr::mutate(SD = ifelse(SD == min(SD), "Low Variance", "High Variance")) %>%
  dplyr::mutate(SD = factor(SD, levels = c("Low Variance", "High Variance"))) %>%
  dplyr::mutate(SAD = ifelse(SAD == "Carrying Capacities", "Normal", SAD)) %>%
  dplyr::mutate(SAD = factor(SAD, levels = c("Normal",
                                             "Log Normal",
                                             "Geometric Series",
                                             "Zipf")))

plSADMean <- ggplot(melt_data, aes(x = MuA, y = ProbStable, color = SigmaA, shape = SigmaA)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  facet_grid(SAD~SD) +
  labs(x = expression("Mean strength of pairwise interactions"~(mu[A])),
       y = "Probability of Stability",
       color = "Pairwise\nInteractions",
       shape = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  scale_color_manual(values = c("ForestGreen", "Black")) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
show(plSADMean)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigRealisticSADs.jpeg",
     width = 2400, height = 2000, res = 300)
plSADMean
dev.off()

