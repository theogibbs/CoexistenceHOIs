########################################################
# Plotting simulations where the pairwise diagonal is not equal in SI Fig. M
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_data_diverse <- read_csv2("simdata/DiverseCommunityVariableDiagonal.csv")

out_data_diverse <- out_data_diverse %>%
  dplyr::select(c("S", "MuA", "SigmaA", "MuD", "SigmaD", "Type", "Stable"))

out_data <- out_data_diverse

melt_data <- out_data %>%
  dplyr::group_by(S, MuA, SigmaA, MuD, SigmaD, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("S", "MuA", "SigmaA", "MuD", "SigmaD", "Type")) %>%
  mutate(variable = ifelse(variable == "ProbFeasible",
                           "Probability of Feasibility",
                           "Probability of Stability")) %>%
  filter(variable == "Probability of Stability")

plVariableDiagonal <- ggplot(melt_data,
                             aes(x = MuA,
                                 y = value,
                                 color = as.factor(MuD),
                                 shape = as.factor(MuD))) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  facet_grid( ~ SigmaD, scales = "free",
              labeller = label_bquote(cols = sigma[D] == .(SigmaD))) +
  labs(x = expression("Mean strength of pairwise interactions"~ (mu[A])),
       y = "Probability of Stability",
       color = "Mean\nintra-specific\ninteraction",
       shape = "Mean\nintra-specific\ninteraction") +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_manual(values = c("Purple" ,"Orange", "DarkBlue")) +
  geom_vline(xintercept = 0, linetype = "dashed")
show(plVariableDiagonal)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigVariableDiagonal.jpeg",
     width = 2000, height = 800, res = 300)
plVariableDiagonal
dev.off()

