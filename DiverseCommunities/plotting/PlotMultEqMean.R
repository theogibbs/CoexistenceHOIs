########################################################
# Bifurcation plotting
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_abds <- read_csv2("simdata/MultipleEquilibriaMean.csv")

test_df <- out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]
zero_cutoff <- 1e-7

if(nrow(test_df) > 0) {
  out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]$Abd <- zero_cutoff
  #out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]$Abd <- NA
  #out_abds <- out_abds %>% filter(Equilibrium == TRUE) %>% filter(Stability == TRUE)
}

summary(test_df)
dim(test_df)

summ_abds <- out_abds %>%
  mutate(Type = case_when(Type == "Feasible Pairwise" ~ "Feasible\nPairwise",
                          Type == "Permuted HOIs" ~ "Permuted\nHOIs",
                          Type == "Constrained HOIs" ~ "Constrained\nHOIs")) %>%
  mutate(Type = factor(Type, levels = c("Feasible\nPairwise",
                                        "Permuted\nHOIs",
                                        "Constrained\nHOIs"))) %>%
#  filter((Equilibrium == TRUE) & (Stability == TRUE)) %>%
  dplyr::group_by(MuA, SigmaA, Type, ParsID) %>%
  dplyr::summarise(Diversity = sum(Abd > zero_cutoff) / unique(S),
                   Coexist = prod(Abd > zero_cutoff),
                   TargetAbd = unique(TargetAbd)) %>%
  mutate(TargetAbd = ifelse(Coexist, 1 - TargetAbd, NA)) %>%
  dplyr::group_by(MuA, SigmaA, Type) %>%
  dplyr::summarise(MeanDiversity = mean(Diversity, na.rm = T),
                   ProbCoexist = mean(Coexist, na.rm = T),
                   ProbTargetAbd = mean(TargetAbd, na.rm = T)) %>%
  melt(id.vars = c("MuA", "SigmaA", "Type")) %>%
  mutate(variable = case_when(variable == "MeanDiversity" ~ "Average Diversity",
                              variable == "ProbCoexist" ~ "Probability of Coexistence",
                              variable == "ProbTargetAbd" ~ "Prob. of Non-target Coexistence")) %>%
  mutate(variable = factor(variable, levels = c("Probability of Coexistence",
                                                "Average Diversity",
                                                "Prob. of Non-target Coexistence")))

summ_abds <- summ_abds %>%
  filter(variable == "Probability of Coexistence") %>%
  mutate(SigmaA = ifelse(SigmaA == 0, "All equal pairwise interactions", "Variable pairwise interactions"))

plMultEq <- ggplot(summ_abds, aes(x = MuA, y = value,
                                  color = Type, shape = Type)) +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = expression("Mean strength of pairwise interactions"~ (mu[A])),
       y = "Probability of Coexistence",
       color = "", shape = "") +
  #facet_grid( ~ SigmaA) +
  theme_classic() +
  theme(text = element_text(size=15),
        legend.key.size = unit(1.5, "cm"),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  ggtitle("B")
show(plMultEq)

out_abds <- read_csv2("simdata/ThreeSpeciesBifurcation.csv") %>%
  mutate(SpID = as.factor(SpID))

plBifurcation <- ggplot(out_abds, aes(x = MuA, y = Abd, color = SpID)) +
  geom_line(alpha = 0.75, size = 1) +
  theme_classic() +
  labs(x = expression("Mean strength of pairwise interactions"~ (mu[A])),
       y = "Equilibrium Abundance", color = "") +
  theme(text = element_text(size=13),
        legend.position = "none",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  ggtitle("A")
plBifurcation

jpeg("../CoexistenceHOIs-Paper/figs/SIFigBifurcation.jpeg",
     width = 3000, height = 1100, res = 300)
grid.arrange(plBifurcation, plMultEq, layout_matrix = matrix(c(1, 1, 2, 2, 2), nrow = 1))
dev.off()

