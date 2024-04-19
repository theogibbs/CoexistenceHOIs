########################################################
# Simulates B tensors with all the interactions to create SI Fig. N
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_data_diverse <- read_csv2("simdata/DiverseCommunityAllTermsB.csv")

out_data_diverse <- out_data_diverse %>%
  dplyr::select(c("S", "MuA", "SigmaA", "Type", "Stable"))

out_data <- out_data_diverse

melt_data <- out_data %>%
  dplyr::group_by(S, MuA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("S", "MuA", "SigmaA", "Type")) %>%
  mutate(variable = ifelse(variable == "ProbFeasible",
                           "Probability of Feasibility",
                           "Probability of Stability")) %>%
  filter(variable == "Probability of Stability")

plAllTerms <- ggplot(melt_data,
                             aes(x = MuA,
                                 y = value,
                                 color = as.factor(SigmaA),
                                 shape = as.factor(SigmaA))) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  labs(x = expression("Mean strength of pairwise interactions"~ (mu[A])),
       y = "Probability of Stability",
       color = expression(sigma[A]),
       shape = expression(sigma[A])) +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_manual(values = c("Purple" ,"Orange")) +
  geom_vline(xintercept = 0, linetype = "dashed")
show(plAllTerms)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigAllTerms.jpeg",
     width = 1500, height = 1000, res = 300)
plAllTerms
dev.off()

