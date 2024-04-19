########################################################
# Plotting diverse mean correlation simulations
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_data <- read_csv2("simdata/DiverseCommunityVariance.csv")

melt_data <- out_data %>%
  dplyr::group_by(MuA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("MuA", "SigmaA", "Type")) %>%
  mutate(variable = ifelse(variable == "ProbFeasible",
                           "Probability of Feasibility",
                           "Probability of Stability")) %>%
  filter(variable == "Probability of Stability")

plVar <- ggplot(melt_data, aes(x = SigmaA, y = value, color = as.factor(MuA))) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(alpha = 0.75, size = 3) +
  theme_classic() +
  labs(x = expression("Pairwise interaction heterogeneity"~ (sigma[A])),
       y = "Probability of Stability", color = expression(mu[A])) +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5))
show(plVar)

jpeg("../CoexistenceHOIs-Paper/figs/DiverseCommunityVariance.jpeg",
     width = 1800, height = 1200, res = 300)
plVar
dev.off()
