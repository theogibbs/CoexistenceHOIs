########################################################
# Plotting the pairwise correlation simulations to make SI Fig. F
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_data <- read_csv2("simdata/DiverseCommunityPairwiseCorrelations.csv")

out_data <- out_data %>%
  dplyr::select(c("S", "MuA", "SigmaA", "RhoA", "Type", "Stable"))

melt_data <- out_data %>%
  dplyr::group_by(S, MuA, SigmaA, RhoA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("S", "MuA", "SigmaA", "RhoA", "Type")) %>%
  mutate(variable = ifelse(variable == "ProbFeasible",
                           "Probability of Feasibility",
                           "Probability of Stability")) %>%
  filter(variable == "Probability of Stability") %>%
  mutate(S = paste(S, "Species")) %>%
  mutate(S = factor(S, levels = c("3 Species", "20 Species"))) %>%
  mutate(Type = factor(Type, levels = c("Feasible Pairwise", "Constrained HOIs"))) %>%
  mutate(RhoA = as.factor(RhoA))

plCorr <- ggplot(melt_data, aes(x = SigmaA, y = value, color = RhoA, shape = RhoA)) +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  facet_wrap(~Type) +
  labs(x = expression("Variation in pairwise interactions"~ (sigma[A])),
       y = "Probability of Stability",
       color = "Pairwise\nCorrelation",
       shape = "Pairwise\nCorrelation") +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_brewer(palette = "Dark2")
show(plCorr)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigPairwiseCorrelations.jpeg",
     width = 2200, height = 1000, res = 300)
plCorr
dev.off()
