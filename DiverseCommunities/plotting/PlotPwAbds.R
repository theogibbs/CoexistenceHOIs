########################################################
# Plotting all of the mean variance simulations
########################################################

# Loading dependencies

source("./Functions.R")

# Plotting script

out_data <- read_csv2("simdata/PwAbds.csv")

melt_data <- out_data %>%
  dplyr::group_by(SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable))

plNonEqual <- ggplot(melt_data, aes(x = SigmaA, y = ProbStable, color = Type, shape = Type)) +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  labs(x = expression("Variation in pairwise interactions"),
       y = "Probability of Stability",
       color = "",
       shape = "") +
  theme(text = element_text(size=15),
        legend.position = c(0.2, 0.2),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_brewer(palette = "Dark2")
show(plNonEqual)

jpeg("../CoexistenceHOIs-Paper/figs/FigPwAbds.jpeg",
     width = 1800, height = 1400, res = 300)
plNonEqual
dev.off()

