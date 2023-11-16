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
sigma_A <- seq(0, 0.5, length.out = 100)
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

out_data <- data.frame()

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0
sigma_A <- 1
rho_A <- 0
b <- 0

in_stat <- crossing(S = S,
                    MuR = mu_R,
                    SigmaR = sigma_R,
                    MuD = mu_D,
                    SigmaD = sigma_D,
                    MuA = mu_A,
                    SigmaA = sigma_A,
                    RhoA = rho_A,
                    b = b)

pars <- BuildThreeSpeciesPars(in_stat)

for(cur_row in 1:nrow(in_pars)) {
  cur_stats <- in_pars[cur_row,]
  
  cur_sigma <- cur_stats$SigmaA
  cur_mu <- cur_stats$MuA
  
  cur_A <- pars$A
  cur_A <- cur_sigma * cur_A
  cur_A <- cur_A + cur_mu
  diag(cur_A) <- diag(pars$A)
  
  cur_pars <- pars
  cur_pars$A <- cur_A
  
  cur_min <- optimize(f = GetThreeSpeciesMinEig, interval = c(0, 10), cur_pars)
  
  hoi_eig <- cur_min$objective
  hoi_stable <- (hoi_eig < 0)
  
  cur_data <- data.frame(TargetAbd = cur_min$minimum,
                         Eigenvalue = hoi_eig,
                         Stable = hoi_stable)
  cur_data <- cbind(cur_stats, cur_data)
  
  out_data <- rbind(out_data, cur_data)
}

melt_data <- out_data %>%
  dplyr::select(MuA, SigmaA, TargetAbd, Eigenvalue) %>%
  melt(id.vars = c("MuA", "SigmaA")) %>%
  mutate(variable = ifelse(variable == "TargetAbd", "Target Abundance", "Leading Eigenvalue"))

pred_targets <- unique(in_pars$MuR) / (unique(in_pars$MuD) + unique(in_pars$MuA))

plMinima <- ggplot(melt_data, aes(x = SigmaA, y = value,
                                   color = as.factor(MuA))) +
  geom_line(size = 1, alpha = 0.5) +
  # geom_point(alpha = 0.75, size = 2) +
  theme_classic() +
  labs(x = expression(sigma[A]),
       y = "", color = expression(mu[A])) +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = pred_targets, linetype = "dashed") +
  facet_grid(~variable)
plMinima

# writing out the data
filename <- "ThreeSpeciesMinima"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)


