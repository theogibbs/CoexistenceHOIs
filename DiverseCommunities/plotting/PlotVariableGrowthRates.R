########################################################
# Plotting the variable growth rates simulations in SI Fig. B
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_data <- read_csv2("simdata/DiverseCommunityGrowthRates.csv")

out_data <- out_data %>%
  dplyr::select(c("S", "SigmaR", "MuA", "SigmaA", "RhoA", "Type", "Stable"))

melt_data <- out_data %>%
  dplyr::group_by(S, SigmaR, MuA, SigmaA, RhoA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("S", "SigmaR", "MuA", "SigmaA", "RhoA", "Type")) %>%
  mutate(variable = ifelse(variable == "ProbFeasible",
                           "Probability of Feasibility",
                           "Probability of Stability")) %>%
  filter(variable == "Probability of Stability") %>%
  mutate(S = paste(S, "Species")) %>%
  mutate(S = factor(S, levels = c("3 Species", "20 Species"))) %>%
  mutate(Type = factor(Type, levels = c("Feasible Pairwise", "Constrained HOIs"))) %>%
  mutate(SigmaA = as.factor(SigmaA))

plVar <- ggplot(melt_data, aes(x = SigmaR, y = value, color = SigmaA, shape = SigmaA)) +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  labs(x = expression("Variation in growth rates"~ (sigma[R])),
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
show(plVar)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigVariableGrowthRates.jpeg",
     width = 1500, height = 1000, res = 300)
plVar
dev.off()
