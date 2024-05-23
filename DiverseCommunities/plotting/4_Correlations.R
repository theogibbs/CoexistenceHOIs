########################################################
# Plotting correlated matrices, spectra and simulation
# results in Fig. 4
########################################################

# Loading dependencies

library(grid)
source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# plotting the matrix addition

set.seed(1)

# choosing parameters
S <- 5
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
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

pars <- BuildPars(in_pars)
target_abd <- GetTargetAbd(pars)

# choosing a skew symmetric matrix to be our Jacobian
C <- BuildA(S = pars$S, mu = 0, sigma = 1, rho = -1)
diag(C) <- 0
# determining the constraints
row_sum_mat <- matrix(pars$r - rowSums(pars$A) - rowSums(C), nrow = pars$S, ncol = pars$S)
diag(row_sum_mat) <- 0
row_sum_mat <- row_sum_mat / (pars$S - 1)
# enforcing the constraints on the Jacobian
C <- C + row_sum_mat

testA <- pars$A
diag(testA) <- 0
# choosing the Btilde matrix to be the difference of the
# Jacobian and the pairwise interactions
Bmat <- C - testA
# parameterizing the Jacobian
J <- - C
diag(J) <- -1

# checking to see if the constraints are satisfied
2 * (pars$r - rowSums(pars$A)) - rowSums(Bmat)

# creating matrix plots of A, B and J
max_val <- max(abs(rbind(pars$A, Bmat, J)))
plA <- GridM(- pars$A, Title = "A", in_color = "PuOr", in_max = max_val)
plPlus <- GridM(matrix(0, nrow = pars$S, ncol = 1),
                Title = "+", in_color = "PuOr", in_max = 0)
plB <- GridM(- Bmat, Title = expression(tilde(B)), in_color = "PuOr", in_max = max_val)
plJ <- GridM(J, Title = "J", in_color = "PuOr", in_max = max_val)

# creating a plot with just a plus sign
plus_data <- data.frame(x = c(0.45, 0.5, 0.55),
                        y = c(0, 0.45, 1),
                        label = c("", "+", ""))

plPlus <- ggplot(plus_data, aes(x = x, y = y, label = label)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        plot.title = element_text(hjust = 0.5, size = 17.5)) +
  coord_fixed() +
  geom_text(size = 10)

# creating a plot with just an equal sign
equa_data <- data.frame(x = c(0.45, 0.5, 0.55),
                        y = c(0, 0.45, 1),
                        label = c("", "=", ""))

plEqual <- ggplot(equa_data, aes(x = x, y = y, label = label)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        plot.title = element_text(hjust = 0.5, size = 17.5)) +
  coord_fixed() +
  geom_text(size = 10)

# creating an empty plot
empty_data <- equa_data
empty_data$label <- " "

plEmpty <- ggplot(empty_data, aes(x = x, y = y, label = label)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        plot.title = element_text(hjust = 0.5, size = 17.5)) +
  coord_fixed() +
  geom_text(size = 10)

# plotting the eigenvalues of each of the matrices with a custom fucntion
plEigs <- PlotEigs(A = -pars$A, B = - Bmat - diag(1, pars$S, pars$S), J = J, Title = "")
plEigs

layout_mat <- rbind(c(1, rep(2, times = 10), 3, 3, rep(4, times = 10), 5, 5, rep(6, times = 10)),
                   c(rep(7, times = 35)))

# binding all the plots together to create panels (A) and (B)
plExamplePart <- grid.arrange(plEmpty, plA, plPlus, plB, plEqual, plJ, plEigs,
                              layout_matrix = layout_mat)

# Plotting script for panel (C)

# reading in simulation data
out_data <- read_csv2("simdata/SimCorr.csv")

# computing the probability of stability for the different cases
# and labeling the data
summ_data <- out_data %>%
  dplyr::group_by(MuA, SigmaA, p) %>%
  dplyr::summarise(ProbStable = mean(Stable)) %>%
  dplyr::mutate(MuA = case_when(MuA == 0 ~ "Zero Mean",
                                MuA < 0 ~ "Facilitative Mean",
                                MuA > 0 ~ "Competitive Mean")) %>%
  dplyr::mutate(MuA = factor(MuA, levels = c("Facilitative Mean", "Zero Mean", "Competitive Mean")))

# plotting the probability of stability as a function of the correlations
# between pairwise and higher order interactions
plConstrainedB <- ggplot(summ_data, aes(x = p, y = ProbStable,
                                        color = as.factor(SigmaA),
                                        shape = as.factor(SigmaA))) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(alpha = 1, size = 3) +
  theme_classic() +
  facet_wrap(~MuA) +
  labs(x = expression("Correlation Strength"~(p)),
       y = "Probability of Stability",
       color = "Variation in Pairwise Interactions",
       shape = "Variation in Pairwise Interactions") +
  theme(text = element_text(size=15),
        legend.position = "top",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  scale_color_viridis(discrete = TRUE, option = "H")
show(plConstrainedB)

# writes out the plot
jpeg("../CoexistenceHOIs-Paper/figs/Fig4Correlations.jpeg",
     width = 3500, height = 3250, res = 300)
grid.arrange(plExamplePart, plConstrainedB,
             layout_matrix = matrix(c(rep(1, times = 2), 2)), nrow = 3, ncol = 1)
dev.off()
