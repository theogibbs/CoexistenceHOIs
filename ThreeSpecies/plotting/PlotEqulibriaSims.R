########################################################
# Plotting three species simulations for multiple
# equilibria
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

out_abds <- read_csv2("simdata/ThreeSpeciesEquilibria.csv")

test_df <- out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]

if(nrow(test_df) > 0) {
  out_abds[((out_abds$Equilibrium == FALSE) | (out_abds$Stability == FALSE)),]$Abd <- zero_cutoff
}

summ_abds <- out_abds %>%
  mutate(Type = factor(Type, levels = c("No HOIs", "Constrained HOIs"))) %>%
  filter((Equilibrium == TRUE) & (Stability == TRUE)) %>%
  dplyr::group_by(MuA, SigmaA, Type, ParsID) %>%
  dplyr::summarise(Phi = sum(Abd > zero_cutoff) / unique(S),
                   TargetAbd = unique(TargetAbd)) %>%
  dplyr::group_by(MuA, SigmaA, Type) %>%
  dplyr::summarise(MeanPhi = mean(Phi),
                   ProbTargetAbd = mean(TargetAbd)) %>%
  melt(id.vars = c("MuA", "SigmaA", "Type")) %>%
  mutate(variable = ifelse(variable == "MeanPhi",
                           "Coexisting Fraction",
                           "Target Abundances"))

plMultEq <- ggplot(summ_abds, aes(x = MuA, y = value, color = variable)) +
  geom_line(size = 2, alpha = 0.5) +
  geom_point(size = 2, alpha = 0.75, ) +
  labs(x = "Interaction Strength", y = "", color = "") +
  facet_grid(~Type) +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15))
show(plMultEq)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesEquilibria.jpeg",
     width = 2200, height = 800, res = 300)
plMultEq
dev.off()

