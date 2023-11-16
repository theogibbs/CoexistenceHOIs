########################################################
# Plotting three species simulations for varying
# target abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

out_abds <- read_csv2("simdata/ThreeSpeciesEquilibriaMV.csv")

test_df <- out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]
zero_cutoff <- 1e-7

if(nrow(test_df) > 0) {
  out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]$Abd <- zero_cutoff
  #out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]$Abd <- NA
  #out_abds <- out_abds %>% filter(Equilibrium == TRUE) %>% filter(Stability == TRUE)
}

summ_abds <- out_abds %>%
  mutate(Type = factor(Type, levels = c("No HOIs", "Constrained HOIs"))) %>%
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
plMultEq <- ggplot(summ_abds, aes(x = SigmaA, y = value, color = as.factor(MuA))) +
  geom_line(size = 1, alpha = 0.5) +
  #geom_point(size = 2, alpha = 0.75, ) +
  labs(x = "Interaction Heterogeneity", y = "", color = expression(mu[A])) +
  facet_wrap(~variable, scales = "free") +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15))
show(plMultEq)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesEquilibriaMV.jpeg",
     width = 3200, height = 900, res = 300)
plMultEq
dev.off()

stab_data <- read_csv2("simdata/ThreeSpeciesVarianceSims.csv")

melt_stab <- stab_data %>%
  mutate(Type = factor(Type, levels = c("No HOIs", "Constrained HOIs"))) %>%
  dplyr::group_by(MuA, SigmaA, Type) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  melt(id.vars = c("MuA", "SigmaA", "Type")) %>%
  mutate(Source = "Equilibrium")

comp_data <- summ_abds %>%
  filter(variable == "Probability of Coexistence") %>%
  mutate(Source = "Dynamics")

melt_data <- rbind(melt_stab, comp_data)

plComp <- ggplot(melt_data, aes(x = SigmaA, y = value, color = Source)) +
  geom_line(size = 1, alpha = 0.5) +
  #geom_point(alpha = 0.75, size = 3) +
  theme_classic() +
  labs(x = expression("Pairwise interaction heterogeneity"~ (sigma[A])),
       y = "Probability of Coexistence", color = "") +
  theme(text = element_text(size=15),
        #legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(~MuA, labeller = label_bquote(cols = mu[A] == .(MuA)))
show(plComp)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesEquilibriumVsDynamics.jpeg",
     width = 2400, height = 900, res = 300)
plComp
dev.off()


