########################################################
# Minimizing the eigenvalues
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

set.seed(1)

# Simulating script

S <- 5
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 0.5
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

target_abd <- GetTargetAbd(pars = pars)
target_abd <- rep(target_abd, times = pars$S)

B <- GetFeasibleB(target_abd, pars)
pars$B <- B

J <- BuildJacobian(eq_abd = target_abd, pars = pars)
GetEig(J)

eig_data <- eigen(J, only.values = TRUE)$values
eig_data <- data.frame(Real = Re(eig_data), Imaginary = Im(eig_data), Source = "Initial")

Bmat <- pars
Bmat$A <- matrix(0, nrow = Bmat$S, ncol = Bmat$S)
Bmat <- BuildJacobian(eq_abd = rep(1, times = Bmat$S), Bmat)
corr_data <- data.frame(PWI = as.vector(pars$A),
                        HOI = as.vector(Bmat),
                        Source = "Initial")

print("Should be zero:")
print((pars$r - rowSums(pars$A) - rowSums(pars$B)))

SwapHOIs <- function(B, pert_ratio = 0.05) {
  retB <- B
  for(i in 1:nrow(B)) {
    cur_vals <- B[i,B[i,] != 0]
    cur_inds <- sample(1:length(cur_vals), size = 2, replace = FALSE)
    new_vals <- cur_vals
    new_vals[cur_inds[1]] <- cur_vals[cur_inds[2]]
    new_vals[cur_inds[2]] <- cur_vals[cur_inds[1]]
    
    cur_pert <- runif(length(new_vals),
                      min = - abs(new_vals) * pert_ratio,
                      max = abs(new_vals) * pert_ratio)
    cur_pert <- cur_pert - mean(cur_pert)
    new_vals <- new_vals + cur_pert
    
    retB[i,retB[i,] != 0] <- new_vals
  }
  return(retB)
}

pars$B <- SwapHOIs(pars$B)

print("Should be zero:")
print((pars$r - rowSums(pars$A) - rowSums(pars$B)))


OldJ <- BuildJacobian(eq_abd = target_abd, pars = pars)
GetEig(OldJ)

GetMinMat <- function(pars, num_trials = 100, pert_ratio = 0.05) {
  
  target_abd <- GetTargetAbd(pars = pars)
  target_abd <- rep(target_abd, times = pars$S)
  
  cur_pars <- pars
  J <- BuildJacobian(eq_abd = target_abd, pars = pars)
  cur_eig <- GetEig(J)
  cur_corr <- GetCorr(J)
  
  out_df <- data.frame(Step = 1:num_trials, Eigenvalue = cur_eig, Correlation = cur_corr)
  
  for(cur_trial in 2:num_trials) {
    new_pars <- cur_pars
    new_pars$B <- SwapHOIs(new_pars$B, pert_ratio = pert_ratio)
    
    J <- BuildJacobian(eq_abd = target_abd, pars = new_pars)
    new_eig <- GetEig(J)
    
    if(new_eig < cur_eig) {
      cur_eig <- new_eig
      cur_pars <- new_pars
      cur_corr <- GetCorr(J)
    }

    out_df$Correlation[cur_trial] <- cur_corr
    out_df$Eigenvalue[cur_trial] <- cur_eig
    
  }
  
  out_list <- list(pars = cur_pars, steps = out_df)
  return(out_list)
}

start_time <- Sys.time()

out_list <- GetMinMat(pars = pars, num_trials = 5e3, pert_ratio = 0.01)

print(Sys.time() - start_time)

out_pars <- out_list$pars
out_steps <- out_list$steps

J <- BuildJacobian(eq_abd = target_abd, pars = out_pars)

opt_data <- eigen(J, only.values = TRUE)$values

opt_data <- data.frame(Real = Re(opt_data), Imaginary = Im(opt_data), Source = "Optimized")
eig_data <- rbind(eig_data, opt_data)

Bmat <- out_pars
Bmat$A <- matrix(0, nrow = Bmat$S, ncol = Bmat$S)
Bmat <- BuildJacobian(eq_abd = rep(1, times = Bmat$S), Bmat)
new_corr_data <- data.frame(PWI = as.vector(out_pars$A),
                            HOI = as.vector(Bmat),
                            Source = "Optimized")

corr_data <- rbind(corr_data, new_corr_data) %>%
  filter(PWI != 1)

ggplot(corr_data, aes(x = PWI, y = HOI)) +
  geom_point() + theme_classic() +
  facet_wrap(~Source)

plEigs <- ggplot(eig_data, aes(x = Real, y = Imaginary)) +
  geom_point() + theme_classic() +
  facet_wrap(~Source) +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        legend.position = "top") +
  geom_vline(xintercept = 0, color = "red")
plEigs

ggplot(out_steps, aes(x = Step, y = Eigenvalue)) +
  geom_line() + theme_classic()

ggplot(out_steps, aes(x = Step, y = Correlation)) +
  geom_line() + theme_classic()

cor(corr_data$PWI[corr_data$Source == "Random"],corr_data$HOI[corr_data$Source == "Random"])
cor(corr_data$PWI[corr_data$Source == "Optimized"],corr_data$HOI[corr_data$Source == "Optimized"])

first_vals <- J[upper.tri(J)]
second_vals <- t(J)[upper.tri(t(J))]
plot(first_vals, second_vals)
cor(first_vals, second_vals)

first_vals <- OldJ[upper.tri(OldJ)]
second_vals <- t(OldJ)[upper.tri(t(OldJ))]
plot(first_vals, second_vals)
cor(first_vals, second_vals)

out_steps$Replicate <- 1
num_opts <- 20
for(i in 2:num_opts) {
  out_list <- GetMinMat(pars = pars, num_trials = 5e3, pert_ratio = 0.01)
  cur_steps <- out_list$steps
  cur_steps$Replicate <- i
  out_steps <- rbind(out_steps, cur_steps)
}

melt_steps <- out_steps %>%
  mutate(Replicate = factor(Replicate)) %>%
  melt(id.vars = c("Step", "Replicate"))

plOpts <- ggplot(melt_steps, aes(x = Step, y = value, color = Replicate)) +
  facet_wrap(~variable, scales = "free") +
  geom_line(alpha = 0.5) + theme_classic() +
  labs(y = "") +
  theme(text = element_text(size=15),
        legend.position = "none",
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15))
plOpts

jpeg("../CoexistenceHOIs-Paper/figs/SIFigOptimizedSpectra.jpeg",
     width = 2200, height = 1000, res = 300)
plEigs
dev.off()

jpeg("../CoexistenceHOIs-Paper/figs/SIFigEmergentCorrelations.jpeg",
     width = 2200, height = 1000, res = 300)
plOpts
dev.off()



