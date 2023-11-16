########################################################
# Plotting the simulations with correlations between
# pairwise and higher-order interactions
########################################################

# Loading dependencies

source("./Functions.R")

# Plotting script

out_data_diverse <- read_csv2("simdata/Correlations.csv")

summ_data <- out_data %>%
  dplyr::group_by(MuA, SigmaA, p) %>%
  dplyr::summarise(ProbStable = mean(Stable), Corr = mean(Corr))

plConstrainedB <- ggplot(summ_data, aes(x = Corr, y = ProbStable,
                                        shape = as.factor(MuA), color = as.factor(SigmaA))) +
  geom_point(size = 1.5, alpha = 0.75) +
  labs(x = expression("Correlation"~(rho)),
       y = "Probability of Stability",
       color = "Variation",
       shape = "Mean Strength") +
  theme_classic() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15))
show(plConstrainedB)

jpeg("../CoexistenceHOIs-Paper/figs/Fig4Correlations.jpeg",
     width = 1200, height = 800, res = 300)
plConstrainedB
dev.off()