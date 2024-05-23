########################################################
# Functions for diverse community simulations
########################################################

# Loading packages and dependencies
library(LaplacesDemon)
library(MASS)
library(hitandrun)
library(ggtern)

# Functions

# builds the pairwise interaction matrix from its statistics
BuildA <- function(S = 100, mu = 0, sigma = 0.1, rho = 0){
  
  Pairs <- SampleCoefficientsFromBivariate(S, mux = mu, muy = mu,
                                           sigmax = sigma, sigmay = sigma,
                                           rho = rho)
  
  A <- matrix(0, S, S)
  k <- 1
  for (i in 1:(S-1)){
    for (j in (i + 1):S){
      A[j,i] <- Pairs[k,1] # fills in the pairs of sampled from a possibly correlated
      A[i,j] <- Pairs[k,2] # bivariate distribution sampled above
      k <- k + 1
    }
  }
  return(A)
}

# builds a higher-order interaction matrix from its statistics without constraints
BuildB <- function(S, MuD, SigmaD, MuB, SigmaB, intra = "None",
                   SelfReg = "Quadratic", scaling = TRUE, DistB = "Normal") {
  
  if(DistB == "Normal") {
    B <- array(rnorm(S*S*S, mean = MuB, sd = SigmaB), c(S, S, S))
  } else if(DistB == "Uniform") {
    B <- array(runif(S*S*S, min = 0, max = 2 * MuB), c(S, S, S))
  }
  
  if(intra == "All Terms") {
    
  } else if(intra == "No Cubic") {
    for(i in 1:S) {
      B[i,i,i] <- 0
    }
    
  } else if(intra == "None") {
    for(i in 1:S) {
      B[i,i,] <- 0
      B[i,,i] <- 0
      B[,i,i] <- 0
    }
  }
  if(scaling) {
    B <- B / S^2
  }
  
  if(SelfReg == "Cubic") {
    diagB <- rnorm(S, mean = MuD, sd = SigmaD)
    for(i in 1:S) {
      B[i,i,i] <- diagB[i]
    }
  }
  
  B <- matrix(as.vector(B), nrow = S, ncol = S^2)
  return(B)
}

# builds each of the ingredients of the model and returns a list with them
BuildPars <- function(in_pars) {
  
  pars <- with(in_pars, {
    pars <- list()
    
    pars$S <- S
    pars$r <- rnorm(S, mean = MuR, sd = SigmaR)
    pars$A <- BuildA(S = S, mu = MuA, sigma = SigmaA, rho = RhoA)
    diag(pars$A) <- rnorm(S, mean = MuD, sd = SigmaD)
    pars$B <- BuildB(S, MuD, SigmaD, MuB, SigmaB, Intra, SelfReg, scaling, DistB)
    return(pars)
  })
  
  return(pars)
}

# computes the growth rates of the current state and parameters for the ODEs
Dynamics <- function(time, state, pars) {
  
  dNdt <- with(pars, {
    state_combs <- as.vector(outer(state, state))
    growth_rates <- r - A %*% state - B %*% state_combs
    dNdt <- state * growth_rates
    return(dNdt)
  })
  
  return(list(dNdt))
}

# returns a higher-order interaction matrix that enforces feasibility when the target abundances are equal
GetFeasibleB <- function(abd, pars) {
  B <- with(pars, {
    total_B <- r / abd^2 - rowSums(A) / abd
    
    num_uniq <- choose(S - 1, 2)
    opt_vals <- rdirichlet(S, alpha = rep(1, times = num_uniq))
    
    B <- array(0, c(S, S, S))
    for(i in 1:S) {
      
      # preventing the HOI terms with repeated indices
      # from being included
      inds <- crossing(j = 1:S, k = 1:S) %>%
        filter(j != k) %>%
        filter(j != i) %>%
        filter(k != i) %>%
        filter(j < k)
      
      # filling in the tensor
      for(cur_row in 1:nrow(inds)) {
        cur_inds <- inds[cur_row,]
        B[i, cur_inds$j, cur_inds$k] <- opt_vals[i, cur_row]
      }
      
    }
    
    B <- matrix(as.vector(B), nrow = S, ncol = S^2)
    B <- B * total_B # multiplying by the appropriate constraint values
    return(B)
  })
  return(B)
}

# returns a higher-order interaction matrix that enforces feasibility when the target abundances are not equal
GetNonEqualFeasibleB <- function(abd, pars) {
  B <- with(pars, {
    total_B <- r - A %*% abd
    
    num_uniq <- choose(S - 1, 2)
    opt_vals <- rdirichlet(S, alpha = rep(1, times = num_uniq))
    
    B <- array(0, c(S, S, S))
    for(i in 1:S) {
      
      # preventing the HOI terms with repeated indices
      # from being included
      inds <- crossing(j = 1:S, k = 1:S) %>%
        filter(j != k) %>%
        filter(j != i) %>%
        filter(k != i) %>%
        filter(j < k)
      
      for(cur_row in 1:nrow(inds)) {
        cur_inds <- inds[cur_row,]
        B[i, cur_inds$j, cur_inds$k] <- opt_vals[i, cur_row]
      }
      
    }
    
    # fills in the tensor using the appropriate multiplied abundances
    outer_abds <- as.vector(outer(abd, abd))
    outer_mat <- diag(1 / outer_abds, nrow = S^2, ncol = S^2)
    B <- matrix(as.vector(B), nrow = S, ncol = S^2) %*% outer_mat
    B <- B * as.vector(total_B)
    return(B)
  })
  return(B)
}

# returns a higher-order interaction matrix that enforces feasibility and ensures the B entries are equal across rows
GetEqualAcrossRowsFeasibleB <- function(abd, pars) {
  B <- with(pars, {
    total_B <- r - A %*% abd
    
    num_uniq <- choose(S - 1, 2)
    
    B <- array(0, c(S, S, S))
    for(i in 1:S) {
      
      # preventing the HOI terms with repeated indices
      # from being included
      inds <- crossing(j = 1:S, k = 1:S) %>%
        filter(j != k) %>%
        filter(j != i) %>%
        filter(k != i) %>%
        filter(j < k)
      
      for(cur_row in 1:nrow(inds)) {
        cur_inds <- inds[cur_row,]
        B[i, cur_inds$j, cur_inds$k] <- 1
      }
      
    }
    
    # fills in the tensor while maintaing equality within each row
    outer_abds <- as.vector(outer(abd, abd))
    outer_mat <- diag(outer_abds, nrow = S^2, ncol = S^2)
    B <- matrix(as.vector(B), nrow = S, ncol = S^2)
    row_sum_abds <- rowSums(B %*% outer_mat)
    B <- B * as.vector(total_B) / row_sum_abds
    })
  return(B)
}

# returns a B matrix with no missing terms
GetAllTermsB <- function(abd, pars) {
  B <- with(pars, {
    total_B <- r / abd^2 - rowSums(A) / abd
    
    num_uniq <- S^2 - 1
    opt_vals <- rdirichlet(S, alpha = rep(1, times = num_uniq))

    B <- array(0, c(S, S, S))
    for(i in 1:S) {
      
      inds <- crossing(j = 1:S, k = 1:S)
      inds <- inds[!((inds$j == inds$k) & (inds$j == i)),]
      for(cur_row in 1:nrow(inds)) {
        cur_inds <- inds[cur_row,]
        B[i, cur_inds$j, cur_inds$k] <- opt_vals[i, cur_row]
      }
      
    }

    B <- matrix(as.vector(B), nrow = S, ncol = S^2)
    B <- B * total_B
    return(B)
  })
  return(B)
}

# returns a B matrix where the variation between entries is driven only by the variance in abundances
GetAbdVarianceFeasibleB <- function(abd, pars) {
  B <- with(pars, {
    total_B <- r - A %*% abd
    
    num_uniq <- choose(S - 1, 2)
    
    B <- array(0, c(S, S, S))
    for(i in 1:S) {
      
      inds <- crossing(j = 1:S, k = 1:S) %>%
        filter(j != k) %>%
        filter(j != i) %>%
        filter(k != i) %>%
        filter(j < k)
      
      for(cur_row in 1:nrow(inds)) {
        cur_inds <- inds[cur_row,]
        B[i, cur_inds$j, cur_inds$k] <- 1 / num_uniq
      }
      
    }
    
    # fills in the B tensor with the only differences coming form the abundance values
    outer_abds <- as.vector(outer(abd, abd))
    outer_mat <- diag(1 / outer_abds, nrow = S^2, ncol = S^2)
    B <- matrix(as.vector(B), nrow = S, ncol = S^2) %*% outer_mat
    B <- B * as.vector(total_B)
  })
  return(B)
}

# returns the Jacobian of the model
BuildJacobian <- function(eq_abd, pars) {
  
  J <- with(pars, {
    
    tensorB <- array(B, c(S, S, S))
    J <- matrix(0, S, S)
    
    for(i in 1:S) {
      for(j in 1:S) {
        hoi_sum <- sum(tensorB[i,j,] * eq_abd + tensorB[i,,j] * eq_abd)
        if(i != j) {
          hoi_sum <- hoi_sum - eq_abd[j] * tensorB[i,j,j]
        }
        
        if(i == j) {
          hoi_sum <- hoi_sum - eq_abd[i] * tensorB[i,i,i]
        }
        
        J[i, j] <- - eq_abd[i] * (A[i,j] + hoi_sum)
      }
    }
    
    return(J)
  })
  return(J)
}

# returns the abundances realized in simulations of the dynamics
# as well as the stability of the resulting state
GetAbds <- function(pars, ini_state, end_time, zero_cutoff) {
  out_abds <- IntegrateDynamics(inistate = ini_state,
                                pars = pars,
                                endtime = end_time,
                                timestep = end_time,
                                fn = Dynamics)
  end_pt <- as.numeric(out_abds[2,2:(pars$S+1)])
  
  end_pt[end_pt < zero_cutoff] <- zero_cutoff^2
  eq_log <- unlist(Dynamics(0, end_pt, pars))
  eq_log <- as.logical(prod(eq_log < zero_cutoff))
  
  stab_log <- BuildJacobian(end_pt, pars)
  stab_log <- GetEig(stab_log)
  stab_log <- as.logical(stab_log)
  
  target_abd <- GetTargetAbd(pars)
  target_log <- all.equal(rep(target_abd, times = pars$S), end_pt)
  target_log <- isTRUE(target_log)
  
  out_abds <- data.frame(Abd = end_pt,
                         Equilibrium = eq_log,
                         Stability = stab_log,
                         TargetAbd = target_log)
  return(out_abds)
}

# gets the real part of the leading eigenvalue of the Jacobian
GetMinEig <- function(target_abd, pars) {
  const_B <- GetFeasibleB(target_abd, pars)
  new_pars <- pars
  pars$B <- const_B
  in_target <- rep(target_abd, times = pars$S)
  J <- BuildJacobian(eq_abd = in_target, pars = pars)
  eig <- GetEig(J)
  return(eig)
}

# translates between vectors and tensors
Vec2Tensor <- function(vec) {
  inS <- (length(vec))^(1/3)
  Bmat <- matrix(vec, nrow = inS, ncol = inS^2, byrow = T)
  B <- array(Bmat, c(inS, inS, inS))
  return(B)
}

# creates a matrix from a vector
Vec2WideMat <- function(num_spec, vec) {
  Bmat <- matrix(vec, nrow = num_spec, ncol = num_spec^2, byrow = T)
  return(Bmat)
}

# returns the index in a vector from tensor indices
GetLongInd <- function(num_spec, i, j, k) {
  
  ret <- num_spec^2 * (i-1) + num_spec * (j-1) + k
  return(ret)
  
}

# implements another method of getting B with variable abundances and additional constraints
BuildConstraintMat <- function(pars) {
  num_spec <- pars$S
  
  const_mat <- matrix(0, nrow = num_spec * (num_spec + 1), ncol = num_spec^3)
  for(i in 1:num_spec) {
    const_mat[i, (1:num_spec^2) + num_spec^2 * (i-1)] <- 1
  }
  
  cur_count <- num_spec + 1
  for(i in 1:num_spec) {
    for(j in 1:num_spec) {
      
      ind_range <- 1:num_spec
      
      first_ind <- GetLongInd(num_spec = num_spec, i = i, j = j, k = ind_range)
      second_ind <- GetLongInd(num_spec = num_spec, i = i, j = ind_range, k = j)
      common_ind <- intersect(first_ind, second_ind)
      
      const_mat[cur_count, first_ind] <- 1
      const_mat[cur_count, second_ind] <- 1
      
      cur_count <- cur_count + 1
    }
  }
  
  return(const_mat)
}

# builds a constraint vector for getting B with the alternative method above
BuildConstraintVector <- function(abd, pars) {
  const_vec <- with(pars, {
    total_B <- r / abd^2 - rowSums(A) / abd
    return(total_B)
  })
  return(const_vec)
}

# gets B using hit and run sampling
GetConstrainedB <- function(target_abd, pars) {
  
  constr_list <- simplexConstraints(n = choose(pars$S - 1, 2))
  
  test_pars <- pars
  test_pars$A <- matrix(0.5, nrow = test_pars$S, ncol = test_pars$S)
  diag(test_pars$A) <- 1
  test_B <- GetFeasibleB(rep(1, times = pars$S), test_pars)
  test_B[test_B != 0] <- 1
  out_B <- matrix(0, nrow = pars$S, ncol = pars$S^2)
  
  total_Bs <- pars$r - pars$A %*% target_abd
  
  state_combs <- as.vector(outer(target_abd, target_abd))
  
  for(i in 1:pars$S) {
    
    cur_constr <- constr_list
    
    cur_outer_abds <- test_B[i,] * state_combs
    cur_outer_abds <- cur_outer_abds[cur_outer_abds != 0]
    cur_constr$constr[1,] <- cur_outer_abds
    
    if(total_Bs[i] < 0) {
      cur_constr$constr[-1,] <- - cur_constr$constr[-1,]
    }
    
    cur_constr$rhs[1] <- total_Bs[i]
    
    cur_B_vals <- hitandrun(constr = cur_constr, n.samples = 1)
    cur_B <- test_B[i,]
    cur_B[cur_B != 0] <- cur_B_vals
    out_B[i,] <- cur_B
    
  }
  
  return(out_B)
}

# returns a B tensor with additional correlations besides the original constraints
GetCorrelatedB <- function(pars, const_mat, const_vec, p) {
  const_inv <- ginv(const_mat)
  const_B <- const_inv %*% const_vec
  
  new_mat <- diag(nrow = nrow(const_inv)) - const_inv %*% const_mat
  rand_sam <- runif(length(const_B))
  const_B <- const_B + new_mat %*% rand_sam
  
  noise_added <- runif(length(const_B),
                       min = - p * max(const_B),
                       max = p * max(const_B))
  
  const_B <- Vec2WideMat(num_spec = pars$S, vec = const_B)
  noise_added <- matrix(noise_added, nrow = nrow(const_B), ncol = ncol(const_B))
  
  noise_added <- noise_added - rowSums(noise_added) / ncol(noise_added)
  const_B <- const_B + noise_added
  
  return(const_B)
}

# computes the Jacobian using Btilde rather than the full B tensor
BuildAbbJacobian <- function(abd, pars, p) {
  
  rand_part <- matrix(runif(pars$S^2, min = 0, max = 1 / pars$S), pars$S, pars$S)
  diag(rand_part) <- 0
  
  row_sum_mat <- matrix(2 * (pars$r - rowSums(pars$A)) - rowSums(rand_part), nrow = pars$S, ncol = pars$S)
  diag(row_sum_mat) <- 0
  row_sum_mat <- row_sum_mat / (pars$S - 1)
  rand_part <- rand_part + row_sum_mat
  
  rand_part <- rand_part - rowSums(rand_part)
  rand_part <- rand_part + 2 * (pars$r / abd^2 - rowSums(pars$A) / abd)
  
  corr_part <- - pars$A
  diag(corr_part) <- 0
  
  Btilde <- p * corr_part + (1-p) * rand_part
  
  row_sum_mat <- matrix(2 * (pars$r - rowSums(pars$A)) - rowSums(Btilde), nrow = pars$S, ncol = pars$S)
  diag(row_sum_mat) <- 0
  row_sum_mat <- row_sum_mat / (pars$S - 1)
  Btilde <- Btilde + row_sum_mat
  
  J <- - abd * (pars$A + Btilde)
  
  return(J)
}

# determines the correlation between the upper and lower triangles of the Jacobian
GetCorr <- function(J) {
  
  first_vals <- J[upper.tri(J)]
  second_vals <- t(J)[upper.tri(t(J))]
  out_cor <- cor(first_vals, second_vals)
  return(out_cor)
}

# reformats eigenvalues into a dataframe
PrepareEigs <- function(M, label) {
  eigs <- eigen(M, only.values = TRUE)$values
  deM <- data.frame(Real = Re(eigs), Imaginary = Im(eigs), Label = label)
  return(deM)
}

# plots eigenvalues like in Fig. 4
PlotEigs <- function(A, B, J, Title = "Spectrum") {
  
  eigs <- rbind(PrepareEigs(A, "Pairwise"),
                PrepareEigs(B, "Higher-Order"),
                PrepareEigs(J, "Pairwise and Higher-Order")) %>%
    mutate(Label = factor(Label, levels = c("Pairwise",
                                            "Higher-Order",
                                            "Pairwise and Higher-Order")))
  
  plJ <- ggplot(data = eigs, aes(x = Real, y = Imaginary)) +
    geom_point(color = "black", alpha = 0.5, size = 2) +
    theme_bw() +
    facet_wrap(~Label) +
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
    ggtitle(Title) +
    geom_rect(aes(xmin = 0,
                  xmax = Inf,
                  ymin = -Inf,
                  ymax = Inf), fill = "red", alpha = 0.05) +
    xlim(c(-2.9, 0.75))
  grob <- grobTree(textGrob("Unstable", x = 0.79,  y = 0.05, hjust = 0,
                            gp = gpar(col = "black", fontsize = 10)))
  # Plot
  plJ <- plJ + annotation_custom(grob)
  grob <- grobTree(textGrob("Stable", x = 0.01,  y = 0.05, hjust = 0,
                            gp = gpar(col = "black", fontsize = 10)))
  plJ <- plJ + annotation_custom(grob)
  return(plJ)
}

