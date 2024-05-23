########################################################
# Plotting all of the target abundance simulations to make Fig. 3 and SI Fig. E
########################################################

# Loading dependencies

source("./Functions.R")

# Plotting script

# reading in the data for the non equal target abundance simulations
out_data <- read_csv2("simdata/NonEqual.csv")

# computing the probability of stability, reshaping and labeling the data
melt_data <- out_data %>%
  dplyr::group_by(MuA, SigmaA, Noise, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  dplyr::mutate(SigmaA = ifelse(SigmaA == 0, "All Equal", "Variable")) %>%
  dplyr::mutate(MuA = case_when(MuA < 0 ~ "Facilitative",
                                MuA == 0 ~ "Zero",
                                MuA > 0 ~ "Competitive")) %>%
  dplyr::mutate(MuA = factor(MuA, levels = c("Facilitative", "Zero", "Competitive")))

# plotting the probability of stability as a function of the variance in target abundances
# for panel (B) of Figure 2
plNonEqual <- ggplot(melt_data, aes(x = Noise, y = ProbStable, color = SigmaA, shape = SigmaA)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  labs(x = expression("Variation in the target abundances"),
       y = "Probability of Stability",
       color = "Pairwise\nInteractions",
       shape = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_manual(values = c("ForestGreen", "Black")) +
  ggtitle("B")
show(plNonEqual)

# Plotting script

# reading in data where the mean of the target abundances is changing
out_data <- read_csv2("simdata/MeanTargets.csv")

# computing the probability of stability, reshaping and labeling the data
melt_data <- out_data %>%
  dplyr::group_by(MuA, TargetAbd, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  dplyr::mutate(MuA = case_when(MuA < 0 ~ "Facilitative",
                                MuA == 0 ~ "Zero",
                                MuA > 0 ~ "Competitive")) %>%
  dplyr::mutate(MuA = factor(MuA, levels = c("Facilitative", "Zero", "Competitive")))

# plotting the probability of stability as a function of the mean of the target abundances
# for panel (A) of Figure 2
plMeanTargets <- ggplot(melt_data, aes(x = TargetAbd, y = ProbStable, color = MuA, shape = MuA)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  labs(x = expression("Target Abundance"~(bar(x))),
       y = "Probability of Stability",
       color = "Pairwise\nInteractions",
       shape = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_manual(values = c("Purple" ,"Orange")) +
  ggtitle("A") +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5)
show(plMeanTargets)

# reading in simulation data for realistic SADs
out_data <- read_csv2("simdata/RealisticSADsCarryingCapacitiesHigherVariance.csv")

# computing the probability of stability, reshaping and labeling the data
melt_data <- out_data %>%
  group_by(MuA, SigmaA, SAD, SD) %>%
  dplyr::summarize(ProbStable = mean(Stable)) %>%
  dplyr::mutate(SigmaA = ifelse(SigmaA == 0, "All Equal", "Variable")) %>%
  dplyr::mutate(SD = ifelse(SD == min(SD), "Low Variance", "High Variance")) %>%
  dplyr::mutate(SD = factor(SD, levels = c("Low Variance", "High Variance"))) %>%
  dplyr::mutate(SAD = ifelse(SAD == "Carrying Capacities", "Normal", SAD)) %>%
  dplyr::mutate(SAD = factor(SAD, levels = c("Normal",
                                             "Log Normal",
                                             "Geometric Series",
                                             "Zipf")))

# plotting the probability of stability as a function of the mean of the target abundances
# for different choices of the SAD to create panel (C) of Figure 2
plSADMean <- ggplot(melt_data, aes(x = MuA, y = ProbStable, color = SigmaA, shape = SigmaA)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  facet_grid(SAD~SD) +
  labs(x = expression("Mean strength of pairwise interactions"~(mu[A])),
       y = "Probability of Stability",
       color = "Pairwise\nInteractions",
       shape = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        panel.spacing = unit(2, "lines"),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  scale_color_manual(values = c("ForestGreen", "Black")) +
  ggtitle("C") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
show(plSADMean)

########################################################
# the following script creates SI Fig. E and follows
# the same basic workflow as the above code
########################################################

out_data <- read_csv2("simdata/ThreeSpeciesTargets.csv")

out_data <- out_data %>%
  dplyr::mutate(MuA = case_when(MuA < 0 ~ "Facilitative",
                                MuA == 0 ~ "Zero",
                                MuA > 0 ~ "Competitive")) %>%
  dplyr::mutate(MuA = factor(MuA, levels = c("Facilitative", "Zero", "Competitive"))) %>%
  filter(MuA != "Zero")

plMinima <- ggplot(out_data, aes(x = TargetAbd,
                                 y = Eigenvalue,
                                 color = MuA,
                                 linetype = MuA)) +
  geom_line(size = 1) +
  theme_classic() +
  labs(x = expression("Target Abundance"~(bar(x))),
       y = "Leading Eigenvalue",
       color = "Pairwise\nInteractions",
       linetype = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        legend.spacing.y = unit(0.25, 'cm')) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  guides(linetype=guide_legend(keywidth = 2.5, keyheight = 1),
         color=guide_legend(keywidth = 2.5, keyheight = 1, byrow = TRUE)) +
  ggtitle("A")

show(plMinima)

###################################################################################

out_data <- read_csv2("simdata/DiverseMinimaTargets.csv")

out_data <- out_data %>%
  dplyr::mutate(MuA = case_when(MuA < 0 ~ "Facilitative",
                                MuA == 0 ~ "Zero",
                                MuA > 0 ~ "Competitive")) %>%
  dplyr::mutate(MuA = factor(MuA, levels = c("Facilitative", "Zero", "Competitive")))

min_data <- out_data %>%
  dplyr::filter(Type == "Minimized") #%>%
  #dplyr::group_by(MuA, SigmaA, Type) %>%
  #dplyr::summarise(MeanAbd = mean(TargetAbd),
  #                 SdAbd = sd(TargetAbd))

plTargetAbds <- ggplot(min_data, aes(x = SigmaA, y = TargetAbd,
                                     color = MuA, linetype = MuA)) +
  #geom_point(size = 3, alpha = 1) +
  #geom_line(size = 0.5, alpha = 0.5) +
  stat_smooth(method="loess", se=T) +
  #geom_linerange(aes(min = MeanAbd - SdAbd, max = MeanAbd + SdAbd)) +
  theme_classic() +
  labs(x = expression("Variation in pairwise interactions"~(sigma[A])),
       y = "Most Stable Target Abundance",
       color = "Pairwise\nInteractions",
       linetype = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  guides(linetype=guide_legend(keywidth = 2.5, keyheight = 1),
         color=guide_legend(keywidth = 2.5, keyheight = 1, byrow = TRUE)) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("B")
plTargetAbds


###################################################################################

out_data <- read_csv2("simdata/DiverseMinimaStability.csv")

stab_data <- out_data %>%
  dplyr::group_by(MuA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  dplyr::mutate(Type = ifelse(Type == "Minimized", "Most Stable", "Carrying\nCapacities"))

plStability <- ggplot(stab_data, aes(x = SigmaA, y = ProbStable,
                                     color = Type, shape = Type)) +
  geom_point(alpha = 1, size = 3) +
  geom_line(size = 0.5, alpha = 0.5) +
  theme_classic() +
  labs(x = expression("Variation in pairwise interactions"~(sigma[A])),
       y = "Probability of Stability",
       color = expression("Target\nAbundance"~(bar(x))),
       shape = expression("Target\nAbundance"~(bar(x)))) +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.spacing.y = unit(0.5, 'cm')) +
  guides(color = guide_legend(byrow = TRUE)) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("C")
plStability

jpeg("../CoexistenceHOIs-Paper/figs/Fig3TargetAbundances.jpeg",
     width = 3200, height = 3700, res = 300)
grid.arrange(plMeanTargets, plNonEqual, plSADMean,
             layout_matrix = rbind(c(1, 2), c(3, 3), c(3, 3), c(3, 3)))
dev.off()

jpeg("../CoexistenceHOIs-Paper/figs/SIFigTargetAbundances.jpeg",
     width = 4400, height = 1000, res = 300)
grid.arrange(plMinima, plTargetAbds, plStability, nrow = 1)
dev.off()

