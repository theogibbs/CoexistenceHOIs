########################################################
# Plotting all of the mean variance simulations
########################################################

# Loading dependencies

source("./Functions.R")

set.seed(1)

S <- 50
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(-1, 0, 1)
sigma_A <- c(0.1, 1)
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

out_data <- data.frame()

for(cur_row in 1:nrow(in_pars)) {
  
  cur_params <- in_pars[cur_row,]
  cur_pars <- BuildPars(cur_params)
  
  target_abd <- GetTargetAbd(cur_pars)
  
  pw_abds <- target_abd
  pw_eigs <- eigen(BuildJacobian(pw_abds, cur_pars), only.values = TRUE)$values
  
  pw_feas_data <- data.frame(Real = Re(pw_eigs),
                             Imaginary = Im(pw_eigs),
                             Type = "Feasible Pairwise")
  
  pw_feas_data <- cbind(cur_params, pw_feas_data)
  
  cur_pars$B <- GetFeasibleB(target_abd, cur_pars)
  
  hoi_abds <- target_abd
  hoi_eigs <- eigen(BuildJacobian(pw_abds, cur_pars), only.values = TRUE)$values
  

  hoi_data <- data.frame(Real = Re(hoi_eigs),
                         Imaginary = Im(hoi_eigs),
                         Type = "Constrained HOIs")
  hoi_data <- cbind(cur_params, hoi_data)
  
  cur_data <- rbind(pw_feas_data, hoi_data)
  out_data <- rbind(out_data, cur_data)
  
}

out_data <- out_data %>%
  mutate(Type = factor(Type, levels = c("Feasible Pairwise", "Constrained HOIs")))

pred_pw_mus <- crossing(S = S,
                        MuD = mu_D,
                        MuA = mu_A,
                        SigmaA = sigma_A,
                        Type = "Feasible Pairwise",
                        Imaginary = 0)
pred_pw_mus$Real <- -(pred_pw_mus$MuD + (pred_pw_mus$S - 1) * pred_pw_mus$MuA)

pred_ho_mus <- crossing(S = S,
                        MuD = mu_D,
                        MuA = mu_A,
                        SigmaA = sigma_A,
                        Type = "Constrained HOIs",
                        Imaginary = 0)
pred_ho_mus$Real <- -(pred_ho_mus$MuD - (pred_ho_mus$S - 1) * pred_ho_mus$MuA)

pred_mus <- rbind(pred_pw_mus, pred_ho_mus) %>%
  filter(MuA != 0)

plSpectra <- ggplot(data = out_data, aes(x = Real, y = Imaginary,
                                         color = Type, shape = Type)) +
  geom_point(alpha = 0.5, size = 2) +
  theme_classic() +
  facet_grid(SigmaA~MuA, scales = "free",
             labeller = label_bquote(cols = mu[A] == .(MuA),
                                     rows = sigma[A] ==.(SigmaA))) +
  labs(color = "", shape = "") +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        legend.position = "top") +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(data = pred_mus, shape = 1, size = 4)
plSpectra

jpeg("../CoexistenceHOIs-Paper/figs/SIFigExampleSpectra.jpeg",
     width = 2200, height = 1400, res = 300)
plSpectra
dev.off()
