########################################################
# Three species simulations minimizing the eigenvalues
# as a function of the target abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Simulating script

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(-0.25, 0, 0.25, 0.5, 2)
sigma_A <- seq(0, 0.5, length.out = 10)
rho_A <- 0
b <- 0

in_pars <- crossing(S = S,
                    MuR = mu_R,
                    SigmaR = sigma_R,
                    MuD = mu_D,
                    SigmaD = sigma_D,
                    MuA = mu_A,
                    SigmaA = sigma_A,
                    RhoA = rho_A,
                    b = b)

num_repl <- 100
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

out_data <- data.frame()

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  
  cur_pars <- BuildThreeSpeciesPars(cur_stats)
  
  cur_min <- optimize(f = GetThreeSpeciesMinEig, interval = c(0, 10), cur_pars)
  
  hoi_eig <- cur_min$objective
  hoi_stable <- (hoi_eig < 0)
  
  hoi_eig <- cur_min$objective
  hoi_stable <- (hoi_eig < 0)
  
  min_data <- data.frame(TargetAbd = cur_min$minimum,
                         Eigenvalue = hoi_eig,
                         Stable = hoi_stable,
                         Type = "Minimized")
  
  target_abd <- rep(GetTargetAbd(cur_pars), times = length(cur_pars$r))
  const_b <- GetThreeSpeciesB(target_abd, cur_pars)
  cur_pars$b <- const_b
  
  hoi_eig <- GetEig(BuildThreeSpeciesJacobian(target_abd, cur_pars))
  hoi_stable <- (hoi_eig < 0)
  
  car_data <- data.frame(TargetAbd = unique(target_abd),
                         Eigenvalue = hoi_eig,
                         Stable = hoi_stable,
                         Type = "Carrying Capacity")
  
  cur_data <- rbind(min_data, car_data)
  cur_data <- cbind(cur_stats, cur_data)
  
  out_data <- rbind(out_data, cur_data)
}

melt_data <- out_data %>%
  mutate(Stable = Eigenvalue < 0) %>%
  dplyr::group_by(MuA, SigmaA, Type) %>%
  dplyr::summarise(TargetAbd = mean(TargetAbd),
                   Eigenvalue = mean(Eigenvalue),
                   Stable = mean(Stable)) %>%
  melt(id.vars = c("MuA", "SigmaA", "Type")) %>%
  mutate(variable = case_when(variable == "TargetAbd" ~ "Target Abundance",
                              variable == "Eigenvalue" ~ "Leading Eigenvalue",
                              variable == "Stable" ~ "Stable"))

pred_targets <- unique(in_pars$MuR) / (unique(in_pars$MuD) + unique(in_pars$MuA))

plMinima <- ggplot(melt_data, aes(x = SigmaA, y = value,
                                  color = as.factor(MuA))) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(alpha = 0.75, size = 2) +
  theme_classic() +
  labs(x = expression(sigma[A]),
       y = "", color = expression(mu[A])) +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_grid(Type~variable) +
  geom_hline(yintercept = pred_targets, alpha = 0.5)
plMinima

eig_data <- melt_data %>%
  filter(variable == "Leading Eigenvalue")

plEigs <- ggplot(eig_data, aes(x = SigmaA, y = value, color = Type)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(alpha = 0.75, size = 2) +
  theme_classic() +
  labs(x = expression(sigma[A]),
       y = "Leading Eigenvalue", color = "") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_grid(~MuA)
plEigs

stab_data <- melt_data %>%
  filter(variable == "Stable")

plStability <- ggplot(stab_data, aes(x = SigmaA, y = value,
                                     color = Type)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(alpha = 0.75, size = 2) +
  theme_classic() +
  labs(x = expression(sigma[A]),
       y = "Probability of Stability", color = expression(mu[A])) +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  facet_grid(~MuA, labeller = label_bquote(cols = mu[A] == .(MuA)))
plStability
# writing out the data
filename <- "ThreeSpeciesMinima"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


