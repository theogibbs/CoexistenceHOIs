########################################################
# Plotting all of the target simulations
########################################################

# Loading dependencies

source("./Functions.R")

# Plotting script

out_data <- read_csv2("simdata/VaryingMeanZeroTargetAbundances.csv")

melt_data <- out_data %>%
  group_by(MuA, SigmaA, SAD, SD) %>%
  dplyr::summarize(ProbStable = mean(Stable), MeanHOI = mean(MeanB)) %>%
  dplyr::mutate(SigmaA = ifelse(SigmaA == 0, "All Equal", "Variable")) %>%
  melt(id.vars = c("MuA", "SigmaA", "SAD", "SD"))

plVarTA <- ggplot(melt_data, aes(x = SD, y = value)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  facet_wrap(~variable, scales = "free") +
  theme_classic() +
  labs(x = expression("Variance in the target abundances"),
       y = "Probability of Stability") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank())
show(plVarTA)

plot_data <- out_data %>%
  mutate(SD = paste("SD =", round(SD, 3)))

plEigs <- ggplot(plot_data, aes(x = Real, y = Imaginary, color = SD)) +
  geom_point(alpha = 0.5, size = 3) +
  theme_classic() +
  facet_wrap(~SD, scales = "free", nrow = 2) +
  labs(color = "Target\nAbundance\nVariation") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
show(plEigs)

hist_data <- out_data %>%
  mutate(PredDensity = dnorm(Real, mean = -1, sd = SD)) %>%
  mutate(SD = paste("SD =", round(SD, 3)))

plEigHists <- ggplot(hist_data, aes(x = Real, y = after_stat(density), color = SD)) +
  geom_histogram(fill = "white") +
  geom_line(aes(x = Real, y = PredDensity)) +
  theme_classic() +
  facet_wrap(~SD, scales = "free", nrow = 2) +
  labs(color = "Target\nAbundance\nVariation",
       shape = "Target\nAbundance\nVariation") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
show(plEigHists)

