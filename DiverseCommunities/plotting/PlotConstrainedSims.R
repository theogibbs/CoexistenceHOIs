########################################################
# Plotting the simulations with correlations between
# pairwise and higher-order interactions
########################################################

# Loading dependencies

source("./Functions.R")

# Plotting script

out_data <- read_csv2("simdata/SimCorr.csv")

summ_data <- out_data %>%
  dplyr::group_by(MuA, SigmaA, p) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  dplyr::mutate(MuA = case_when(MuA == 0 ~ "Zero Mean",
                                MuA < 0 ~ "Facilitative Mean",
                                MuA > 0 ~ "Competitive Mean")) %>%
  dplyr::mutate(MuA = factor(MuA, levels = c("Facilitative Mean", "Zero Mean", "Competitive Mean")))

plConstrainedB <- ggplot(summ_data, aes(x = p, y = ProbStable,
                                        color = as.factor(SigmaA),
                                        shape = as.factor(SigmaA))) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  facet_wrap(~MuA) +
  labs(x = expression("Correlation Strength"~(p)),
       y = "Probability of Stability",
       color = "Variation in \nPairwise\nInteractions",
       shape = "Variation in \nPairwise\nInteractions") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  scale_color_viridis(discrete = TRUE, option = "H")
show(plConstrainedB)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigCorrelations.jpeg",
     width = 3200, height = 1200, res = 300)
plConstrainedB
dev.off()