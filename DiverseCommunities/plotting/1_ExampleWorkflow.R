########################################################
# Creates matrices, plots of the dynamics and spectra
# graphs for panel (D) of Fig. 1
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

set.seed(2)

# setting the parameters
S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- -0.25
sigma_A <- 0.75
rho_A <- 0
mu_B <- 0
sigma_B <- 0
intra <- "None"
self_reg <- "Quadratic"
scaling <- FALSE
dist_B <- "Normal"

in_pars <- crossing(S = S,
                    MuR = mu_R,
                    SigmaR = sigma_R,
                    MuD = mu_D,
                    SigmaD = sigma_D,
                    MuA = mu_A,
                    SigmaA = sigma_A,
                    RhoA = rho_A,
                    MuB = mu_B,
                    SigmaB = sigma_B,
                    Intra = intra,
                    SelfReg = self_reg,
                    scaling = scaling,
                    DistB = dist_B)

# building parameters and choosing target abundances
pars <- BuildPars(in_pars)
target_abd <- GetTargetAbd(pars, choice = "Carrying Capacities", value = 0.5)

# setting teh initial conditions and time to integrate the dynamics for
ini_state <- runif(S, min = 0, max = 0.05)
end_time <- 50
time_step <- 0.1

# integrating the dynaics with only pairwise interactions
out_pw <- IntegrateDynamics(inistate = ini_state,
                            pars = pars,
                            endtime = end_time,
                            timestep = time_step,
                            fn = Dynamics)
out_pw <- melt(out_pw, id.vars = "time")
out_pw$Type <- "No HOIs"

# determining the constrained HOIs
const_B <- GetFeasibleB(as.vector(target_abd), pars)
pars$B <- const_B

# building the Jacobian
J <- BuildJacobian(target_abd, pars)

# plotting the target abundances and pairwise interactions
max_val <- max(abs(cbind(pars$A, target_abd, J)))
plTargets <- GridM(- target_abd, Title = expression("Target Abundances"~(bar(x))), in_color = "PuOr", in_max = max_val) +
  theme(plot.title = element_text(size = 14, vjust = -3))
plA <- GridM(- pars$A, Title = expression("Pairwise Interactions"~(A)), in_color = "PuOr", in_max = max_val) +
  theme(plot.title = element_text(size = 14, vjust = -3))
grid.arrange(plTargets, plA, nrow = 1)

# plotting the eigenvalues of the target abundance equilibrium
plot_eigs <- PrepareEigs(J, label = "J")

plEigs <- ggplot(data = plot_eigs, aes(x = Real, y = Imaginary)) +
  geom_point(color = "black", alpha = 0.5, size = 2) +
  theme_bw() +
  theme(
    panel.spacing = unit(4, "lines"),
    legend.position="none",
    panel.background=element_blank(),
    strip.background = element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    text = element_text(size=15),
    legend.text=element_text(size = 15)
  ) +
  ggtitle("") +
  geom_rect(aes(xmin = 0,
                xmax = Inf,
                ymin = -Inf,
                ymax = Inf), fill = "red", alpha = 0.05) +
  xlim(c(-2.9, 0.75))
plEigs


# adding text to the eigenvalue plot
grob <- grobTree(textGrob("Unstable", x = 0.79,  y = 0.05, hjust = 0,
                          gp = gpar(col = "black", fontsize = 10)))
plEigs <- plEigs + annotation_custom(grob)
grob <- grobTree(textGrob("Stable", x = 0.01,  y = 0.05, hjust = 0,
                          gp = gpar(col = "black", fontsize = 10)))
plEigs <- plEigs + annotation_custom(grob)


# adding lines to the dynamics plot
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
end_pts <- out_pw %>% filter(time == max(time))
end_pts <- end_pts$value

# plotting the dynamics and adding arrows
plDyn <- ggplot() +
  geom_line(data = out_pw, aes(x = time, y = value, color = variable),
            size = 2, alpha = 0.5) +
  xlab("Time") + ylab("Abundance") +
  theme_classic() +
  ggtitle("") +
  theme(legend.position = "none",
        text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = target_abd, color = cbPalette[1:3], linetype = "dashed", size = 1, alpha = 0.5) +
  scale_colour_manual(values=cbPalette) +
  geom_segment(aes(x = 42,
                   y = end_pts[1],
                   xend = 42,
                   yend = target_abd[1]),
               size = 1,
               show.legend = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               color = cbPalette[1],
               alpha = 0.75) +
  geom_segment(aes(x = 44,
                   y = end_pts[2],
                   xend = 44,
                   yend = target_abd[2]),
               size = 1,
               show.legend = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               color = cbPalette[2],
               alpha = 0.75) +
  geom_segment(aes(x = 46,
                   y = end_pts[3],
                   xend = 46,
                   yend = target_abd[3]),
               size = 1,
               show.legend = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               color = cbPalette[3],
               alpha = 0.75)
show(plDyn)

# writing out the figure
layout_mat <- matrix(data = c(1, 2, 2, 3, 3, 4, 4), nrow = 1)
jpeg("../CoexistenceHOIs-Paper/figs/Fig1Workflow.jpeg",
     width = 4200, height = 1100, res = 300)
grid.arrange(plTargets, plA, plDyn, plEigs, layout_matrix = layout_mat)
dev.off()

