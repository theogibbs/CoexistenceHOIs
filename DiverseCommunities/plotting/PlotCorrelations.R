########################################################
# Plotting correlated matrix spectra
########################################################

# Loading dependencies

library(grid)
source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# plotting the matrix addition

set.seed(1)

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

C <- BuildA(S = pars$S, mu = 0, sigma = 1, rho = -1)
diag(C) <- 0
row_sum_mat <- matrix(pars$r - rowSums(pars$A) - rowSums(C), nrow = pars$S, ncol = pars$S)
diag(row_sum_mat) <- 0
row_sum_mat <- row_sum_mat / (pars$S - 1)
C <- C + row_sum_mat

testA <- pars$A
diag(testA) <- 0
Bmat <- C - testA
J <- - C
diag(J) <- -1

2 * (pars$r - rowSums(pars$A)) - rowSums(Bmat)

max_val <- max(abs(rbind(pars$A, Bmat, J)))
plA <- GridM(- pars$A, Title = "A", in_color = "PuOr", in_max = max_val)
plPlus <- GridM(matrix(0, nrow = pars$S, ncol = 1),
                Title = "+", in_color = "PuOr", in_max = 0)
plB <- GridM(- Bmat, Title = expression(tilde(B)), in_color = "PuOr", in_max = max_val)
plJ <- GridM(J, Title = "J", in_color = "PuOr", in_max = max_val)

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

plEigs <- PlotEigs(A = -pars$A, B = - Bmat - diag(1, pars$S, pars$S), J = J, Title = "")

plEigs

layout_mat <- rbind(c(1, rep(2, times = 10), 3, 3, rep(4, times = 10), 5, 5, rep(6, times = 10)),
                   c(rep(7, times = 35)))

grid.arrange(plEmpty, plA, plPlus, plB, plEqual, plJ, plEigs,
             layout_matrix = layout_mat)

jpeg("../CoexistenceHOIs-Paper/figs/Fig4Correlations.jpeg",
     width = 3200, height = 2000, res = 300)

grid.arrange(plEmpty, plA, plPlus, plB, plEqual, plJ, plEigs,
             layout_matrix = layout_mat)
dev.off()
