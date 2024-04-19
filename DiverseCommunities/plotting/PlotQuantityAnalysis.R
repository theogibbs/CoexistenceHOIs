PredMeanJacobian <- function(S, mu_A, R, d, x) {
  
  sum_quant <- 0
  
  for(i in 1:S) {
    j_inds <- 1:S
    j_inds <- j_inds[j_inds != i]
    
    for(j in j_inds) {
      
      no_ij_inds <- 1:S
      no_ij_inds <- no_ij_inds[no_ij_inds != i]
      no_ij_inds <- no_ij_inds[no_ij_inds != j]
      
      cur_quant <- x[i] * (R - (d - mu_A) * x[i] - S * mu_A)
      cur_quant <- cur_quant * sum(x[no_ij_inds])
      
      double_sum <- outer(x[j_inds], x[j_inds])
      double_sum <- double_sum[upper.tri(double_sum)]
      double_sum <- sum(double_sum)
      
      cur_quant <- cur_quant / double_sum
      
      sum_quant <- sum_quant + cur_quant
      
    }
  }
  
  ret_quant <- - S * (d + (S-1) * mu_A) - sum_quant
  
  ret_quant <- ret_quant / S
  return(ret_quant)
}

#target_abd <- GetTargetAbd(pars, value = 1)
#pars$B <- GetEqualAcrossRowsFeasibleB(target_abd, pars)
#PredMeanJacobian(S = S, mu_A = mu_A, R = R, d = d, x = target_abd) - mean(rowSums(BuildJacobian(target_abd, pars)))

S <- 10
mu_A <- 0
R <- 1
d <- 1
sd_val <- 0.25

# errors in your notes formula
# should be a x_i and a (1-delta_ij)

x <- rnorm(n = S, mean = 1, sd = sd_val)
x <- x / mean(x)
x[x <= 0] <- 1e-2

PredMeanJacobian(S = S, mu_A = mu_A, R = R, d = d, x = x)

num_repl <- 20

in_sd_val <- seq(0, 0.5, length.out = 5)
in_mu_A <- seq(-10, 5, length.out = 10) / S

out_means <- data.frame()

for(sd_val in in_sd_val) {
  print(sd_val)
  for(mu_A in in_mu_A) {
    for(k in 1:num_repl) {
      
      x <- rnorm(n = S, mean = 1, sd = sd_val)
      x <- x / mean(x)
      x[x <= 0] <- 1e-2
      
      cur_mean <- PredMeanJacobian(S = S, mu_A = mu_A, R = R, d = d, x = x)
      
      out_means <- rbind(out_means, data.frame(MuA = mu_A, SD = sd_val, PredMean = cur_mean))
    }
  }
}

summ_means <- out_means %>%
  group_by(MuA, SD) %>%
  summarise(AvgPredMean = mean(PredMean))

plot_means <- merge(summ_means, melt_data)

ggplot(plot_means, aes(x = AvgPredMean, y = MeanJacobian, color = as.factor(SD))) +
  facet_wrap(~SD, scales = "free", nrow = 1) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(color = "Target\nAbundance\nVariation") +
  theme(text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background = element_blank()) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)



