########################################################
# Plotting three species simulations for varying
# target abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Plotting script

out_data <- read_csv2("simdata/ThreeSpeciesTargets.csv")

out_data <- out_data %>%
  dplyr::mutate(MuA = case_when(MuA < 0 ~ "Facilitative",
                                MuA == 0 ~ "Zero",
                                MuA > 0 ~ "Competitive")) %>%
  dplyr::mutate(MuA = factor(MuA, levels = c("Facilitative", "Zero", "Competitive")))


plTargets <- ggplot(out_data, aes(x = TargetAbd,
                                  y = Eigenvalue,
                                  color = MuA,
                                  linetype = MuA)) +
  geom_line(size = 1.5) +
  theme_classic() +
  labs(x = expression("Target Abundance"~(bar(x))),
       y = "Leading Eigenvalue",
       color = "Pairwise\nInteractions",
       linetype = "Pairwise\nInteractions") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  guides(linetype=guide_legend(keywidth = 2.5, keyheight = 1),
         color=guide_legend(keywidth = 2.5, keyheight = 1)) +
  ggtitle("B")
  
show(plTargets)

jpeg("../CoexistenceHOIs-Paper/figs/Fig3Eigenvalues.jpeg",
     width = 1600, height = 1000, res = 300)
plTargets
dev.off()

