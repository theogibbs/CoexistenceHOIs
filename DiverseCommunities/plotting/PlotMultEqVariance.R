########################################################
# Plotting three species simulations for varying
# target abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Plotting script

out_abds <- read_csv2("simdata/MultipleEquilibriaVariance.csv")

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
  mutate(Type = factor(Type, levels = c("Feasible Pairwise",
                                        "Permuted HOIs",
                                        "Constrained HOIs"))) %>%
  filter((Equilibrium == TRUE) & (Stability == TRUE)) %>%
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
  filter(variable == "Probability of Coexistence")

plMultEq <- ggplot(summ_abds, aes(x = SigmaA, y = value,
                                  color = Type, shape = Type)) +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = expression("Variation in pairwise interactions"~ (mu[A])),
       y = "Probability of Coexistence",
       color = "", shape = "") +
#  facet_grid(~variable) +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) + ggtitle("Variation dependence")
show(plMultEq)

jpeg("../CoexistenceHOIs-Paper/figs/Fig6MultEqVariance.jpeg",
     width = 1800, height = 1000, res = 300)
plMultEq
dev.off()

