########################################################
# Plotting all of the mean variance simulations to create Fig. 2 and supplementary figures
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_data_diverse <- read_csv2("simdata/DiverseCommunityVariance.csv")
out_data_three <- read_csv2("simdata/ThreeSpeciesVarianceSims.csv")

out_data_diverse <- out_data_diverse %>%
  dplyr::select(c("S", "MuA", "SigmaA", "Type", "Stable"))

out_data_three <- out_data_three %>%
  dplyr::select(c("S", "MuA", "SigmaA", "Type", "Stable"))

out_data <- rbind(out_data_diverse, out_data_three)

melt_data <- out_data %>%
  dplyr::group_by(S, MuA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("S", "MuA", "SigmaA", "Type")) %>%
  mutate(variable = ifelse(variable == "ProbFeasible",
                           "Probability of Feasibility",
                           "Probability of Stability")) %>%
  filter(variable == "Probability of Stability") %>%
  mutate(S = paste(S, "Species")) %>%
  mutate(S = factor(S, levels = c("3 Species", "20 Species"))) %>%
  mutate(MeanIntStr = case_when(MuA < 0 ~ "Facilitative",
                                MuA == 0 ~ "Zero",
                                MuA > 0 ~ "Competitive")) %>%
  mutate(MeanIntStr = factor(MeanIntStr,
                             levels = c("Facilitative", "Zero", "Competitive"))) %>%
  filter(S == "20 Species") %>%
  filter(MeanIntStr != "Zero")

plVar <- ggplot(melt_data, aes(x = SigmaA, y = value, color = MeanIntStr, shape = MeanIntStr)) +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
 # facet_wrap(~S, scales = "free") +
  labs(x = expression("Variation in pairwise interactions"~ (sigma[A])),
       y = "Probability of Stability",
       color = "Pairwise\nInteractions",
       shape = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_manual(values = c("Purple" ,"Orange")) +
  ggtitle("B")
show(plVar)

out_data_diverse <- read_csv2("simdata/DiverseCommunityMean.csv")
out_data_three <- read_csv2("simdata/ThreeSpeciesMeanSims.csv")

out_data_diverse <- out_data_diverse %>%
  dplyr::select(c("S", "MuA", "SigmaA", "Type", "Stable"))

out_data_three <- out_data_three %>%
  dplyr::select(c("S", "MuA", "SigmaA", "Type", "Stable"))

out_data <- rbind(out_data_diverse, out_data_three)

melt_data <- out_data %>%
  dplyr::group_by(S, MuA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("S", "MuA", "SigmaA", "Type")) %>%
  mutate(variable = ifelse(variable == "ProbFeasible",
                           "Probability of Feasibility",
                           "Probability of Stability")) %>%
  filter(variable == "Probability of Stability") %>%
  mutate(S = paste(S, "Species")) %>%
  mutate(S = factor(S, levels = c("3 Species", "20 Species"))) %>%
  mutate(VarIntStr = case_when(SigmaA == 0 ~ "All Equal",
                                SigmaA > 0 ~ "Variable")) %>%
  filter(S == "20 Species")

plMean <- ggplot(melt_data, aes(x = MuA, y = value, color = VarIntStr, shape = VarIntStr)) +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
#  facet_wrap(~S, scales = "free") +
  labs(x = expression("Mean strength of pairwise interactions"~ (mu[A])),
       y = "Probability of Stability",
       color = "Pairwise\nInteractions",
       shape = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_color_manual(values = c("ForestGreen", "Black")) +
  ggtitle("A") +
  ylim(c(-0.125, 1)) +
  xlim(-0.75, 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = -0.44,
                   y = -0.125,
                   xend = -0.75,
                   yend = -0.125),
               size = 0.5,
               show.legend = FALSE,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 0.025,
                   y = -0.125,
                   xend = 0.35,
                   yend = -0.125),
               size = 0.5,
               show.legend = FALSE,
               arrow = arrow(length = unit(0.25, "cm")))

grob <- grobTree(textGrob("More Facilitative", x = 0.0775,  y = 0.078, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
plMean <- plMean + annotation_custom(grob)

grob <- grobTree(textGrob("More Competitive", x = 0.69,  y = 0.078, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
plMean <- plMean + annotation_custom(grob)
show(plMean)

jpeg("../CoexistenceHOIs-Paper/figs/Fig2VarMean.jpeg",
     width = 3600, height = 1200, res = 300)
grid.arrange(plMean, plVar, nrow = 1)
dev.off()

#jpeg("../CoexistenceHOIs-Paper/figs/Fig3Var.jpeg",
#     width = 2400, height = 900, res = 300)
#plVar
#dev.off()
