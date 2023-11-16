########################################################
# Plotting the alternative sampling simulations
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_data <- read_csv2("simdata/DiverseCommunityAlternativeSampling.csv")

out_data <- out_data %>%
  dplyr::select(c("S", "MuA", "SigmaA", "RhoA", "SigmaB", "Type", "Stable"))

melt_data <- out_data %>%
  dplyr::group_by(S, MuA, SigmaA, SigmaB, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("S", "MuA", "SigmaA", "SigmaB", "Type")) %>%
  mutate(variable = ifelse(variable == "ProbFeasible",
                           "Probability of Feasibility",
                           "Probability of Stability")) %>%
  filter(variable == "Probability of Stability") %>%
  mutate(S = paste(S, "Species")) %>%
  mutate(S = factor(S, levels = c("3 Species", "20 Species"))) %>%
  mutate(Type = factor(Type, levels = c("Feasible Pairwise", "Constrained HOIs"))) %>%
  mutate(SigmaA = as.factor(SigmaA))

plAlt <- ggplot(melt_data, aes(x = SigmaB, y = value, color = SigmaA, shape = SigmaA)) +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  facet_grid(~MuA ,labeller = label_bquote(cols = mu[A] == .(MuA))) +
  labs(x = expression("Variation in higher-order interactions"~ (sigma[B])),
       y = "Probability of Stability",
       color = expression(sigma[A]),
       shape = expression(sigma[A])) +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_brewer(palette = "Dark2")
show(plAlt)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigAlternativeSampling.jpeg",
     width = 2200, height = 1000, res = 300)
plAlt
dev.off()
