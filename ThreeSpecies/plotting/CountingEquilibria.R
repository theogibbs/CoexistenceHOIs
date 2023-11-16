########################################################
# Three species simulations to count
# feasible equilibria
########################################################

# Loading dependencies

source("./Functions.R")
source("./ThreeSpecies/ThreeSpeciesFunctions.R")

# Simulating script

GetEquilibria <- function(num_trials, pars) {
  target_abd <- GetTargetAbd(pars)
  out_sols <- data.frame(Abd1 = -1, Abd2 = -1, Abd3 = -1, Stable = 1, TermCd = -1)
  for(i in 1:num_trials) {
    sol_obj <- nleqslv(x = runif(3, min = 0, max = 2 * target_abd),
                       fn = ThreeSpeciesGrowthRates,
                       pars = pars)
    
    if(sol_obj$termcd == 1) {
      
      cur_sol <- data.frame(Abd1 = round(sol_obj$x[1], 5),
                            Abd2 = round(sol_obj$x[2], 5),
                            Abd3 = round(sol_obj$x[3], 5),
                            TermCd = sol_obj$termcd)
      
      row_log <- !nrow(plyr::match_df(out_sols, cur_sol,
                                      on = c("Abd1", "Abd2", "Abd3", "TermCd")))
      
      if(row_log) {
        
        cur_eq <- sol_obj$x
        cur_eig <- GetEig(BuildThreeSpeciesJacobian(cur_eq, pars))
        cur_stable <- cur_eig < 0
        cur_sol$Stable <- cur_stable
        
        out_sols <- rbind(out_sols, cur_sol)
      }
    }
    
  }
  
  out_sols <- out_sols %>%
    filter(TermCd == 1)
  
  if(nrow(out_sols) > 0) {
    out_sols$SolnID <- 1:nrow(out_sols)
    out_sols <- melt(out_sols[,-5], id.vars = c("SolnID", "Stable"))
    out_sols <- out_sols[order(out_sols$SolnID),]
    out_sols <- out_sols[,c("SolnID", "variable", "value", "Stable")]
  } else {
    out_sols <- data.frame(SolnID = -1,
                           variable = c("Abd1", "Abd2", "Abd3"),
                           value = -1,
                           Stable = 1)
  }
  
  return(out_sols)
}

S <- 3
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- 0.75
sigma_A <- seq(0.1, 1.5, length.out = 10)
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

num_trials <- 100
out_data <- data.frame()

num_repl <- 1
in_pars <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
in_pars$ParsID <- runif(nrow(in_pars))

start_time <- Sys.time()

for(cur_row in 1:nrow(in_pars)) {
  cur_pars <- in_pars[cur_row,]
  
  pars <- BuildThreeSpeciesPars(cur_pars)
  
  cur_pw_abds <- GetPairwiseAbd(pars)
  cur_pw <- data.frame(SolnID = 1, 
                       variable = c("Abd1", "Abd2", "Abd3"),
                       value = cur_pw_abds)
  cur_pw_eig <- GetEig(BuildThreeSpeciesJacobian(eq_abd = as.vector(cur_pw_abds), pars = pars))
  cur_pw$Stable <- (cur_pw_eig < 0)
  cur_pw$Type <- "Pairwise"
  
  target_abd <- GetTargetAbd(pars)
  pars$b <- GetThreeSpeciesB(target_abd, pars)
  
  cur_const <- GetEquilibria(num_trials = num_trials, pars = pars)
  cur_const$Type <- "Constrained HOIs"
  
  new_b <- rnorm(3, mean = mean(pars$b), sd = sd(pars$b))
  pars$b <- new_b
  
  cur_rand <- GetEquilibria(num_trials = num_trials, pars = pars)
  cur_rand$Type <- "Random HOIs"
  
  cur_data <- rbind(cur_pw, cur_const, cur_rand)
  cur_data <- cbind(cur_pars, cur_data)
  
  out_data <- rbind(out_data, cur_data)
}

print(Sys.time() - start_time)

# writing out the data
filename <- "ThreeSpeciesCountingEquilibria"
cur_file <- paste0("./simdata/", filename, ".csv")
write_csv2(out_data, file = cur_file)

summ_data <- out_data %>%
  dplyr::group_by(MuA, SigmaA, SolnID, Type, ParsID) %>%
  dplyr::summarise(Feasible = prod(value > 0), Stable = unique(Stable)) %>%
  dplyr::group_by(MuA, SigmaA, Type, ParsID) %>%
  dplyr::summarise(SumFeas = sum(Feasible), StabFeas = sum(Feasible * Stable)) %>%
  dplyr::group_by(MuA, SigmaA, Type) %>%
  dplyr::summarise(MeanSumFeas = mean(SumFeas),
                   SdSumFeas = sd(SumFeas),
                   MeanStabFeas = mean(StabFeas),
                   SdStabFeas = sd(StabFeas)) %>%
  mutate(Type = factor(Type, levels = c("Pairwise", "Random HOIs", "Constrained HOIs")))

plFeasEq <- ggplot(summ_data, aes(x = SigmaA, y = MeanSumFeas, color = Type)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(size = 3) +
  geom_errorbar(size = 0.5, width = 0.05, alpha = 0.5, aes(ymin = MeanSumFeas - SdSumFeas, ymax = MeanSumFeas + SdSumFeas)) +
  labs(x = "Interaction Heterogeneity", y = "Number of Feasible Equilibria", color = "") +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5)
show(plFeasEq)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesFeasEq.jpeg",
     width = 2200, height = 1300, res = 300)
plFeasEq
dev.off()

plStabFeas <- ggplot(summ_data, aes(x = SigmaA, y = MeanStabFeas, color = Type)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(size = 3) +
  geom_errorbar(size = 0.5, width = 0.05, alpha = 0.5, aes(ymin = MeanStabFeas - SdStabFeas, ymax = MeanStabFeas + SdStabFeas)) +
  labs(x = "Interaction Heterogeneity", y = "Number of Feasible and Stable Equilibria", color = "") +
  theme_classic() +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5)
show(plStabFeas)

jpeg("../CoexistenceHOIs-Paper/figs/ThreeSpeciesStabFeasEq.jpeg",
     width = 2200, height = 1300, res = 300)
plStabFeas
dev.off()

