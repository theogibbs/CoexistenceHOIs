########################################################
# Plotting all of the mean variance simulations
########################################################

# Loading dependencies

source("./Functions.R")

# Plotting script

out_data <- read_csv2("simdata/DiverseCommunityDiversityVariance.csv")

summ_data <- out_data %>%
  dplyr::group_by(S, MuA, SigmaA, Type) %>%
  dplyr::summarise(AvgEig = mean(Eigenvalue), SdEig = sd(Eigenvalue)) %>%
  dplyr::mutate(TestStat = case_when(MuA == 0 ~ "Variance",
                                     MuA < 0 ~ "Facilitative Mean",
                                     MuA > 0 ~ "Competitive Mean")) %>%
  mutate(TestStat = factor(TestStat, levels = c("Variance",
                                                "Competitive Mean",
                                                "Facilitative Mean"))) %>%
  mutate(Type = factor(Type, levels = c("Feasible Pairwise", "Constrained HOIs")))

plDivVarMean <- ggplot(summ_data, aes(x = S, y = AvgEig, color = Type, shape = Type)) +
  geom_point(alpha = 0.75, size = 3) +
  geom_linerange(alpha = 0.5, aes(ymin = AvgEig - SdEig, ymax = AvgEig + SdEig)) +
  theme_classic() +
  labs(x = expression("Number of species"~ (S)),
       y = "Leading Eigenvalue", color = " ", shape = " ") +
  facet_wrap(~TestStat, scales = "free") +
  theme(text = element_text(size=15),
        legend.position ="top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")
show(plDivVarMean)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigDisentangling.jpeg",
     width = 2500, height = 1000, res = 300)
plDivVarMean
dev.off()

var_data <- summ_data %>%
  filter(TestStat == "Variance")

plDivVar <- ggplot(var_data, aes(x = S, y = AvgEig, color = Type)) +
  #geom_line(size = 1, alpha = 0.5) +
  geom_point(alpha = 0.75, size = 3, shape = 15) +
  geom_linerange(alpha = 0.5, size = 0.5, aes(ymin = AvgEig - SdEig, ymax = AvgEig + SdEig)) +
  theme_classic() +
  labs(x = expression("Number of species"~ (S)),
       y = "Leading Eigenvalue", color = " ") +
  theme(text = element_text(size=12.5),
        legend.position = c(0.7, 0.2),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")
show(plDivVar)
