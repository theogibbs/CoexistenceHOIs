########################################################
# Diverse community simulations  with correlated HOIs
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

S <- 5
mu_R <- 1
sigma_R <- 0
mu_D <- 1
sigma_D <- 0
mu_A <- c(0)
sigma_A <- c(1)
mu_B <- 0
sigma_B <- 0
rho_A <- 0
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

pars <- BuildPars(in_pars[1,])
target_abd <- GetTargetAbd(pars)
C <- matrix(runif(pars$S^2), pars$S, pars$S)
C <- C - t(C)
diag(C) <- 0



total_Bs <- pars$r - pars$A %*% target_abd
total_Bs <- total_Bs / target_abd^2

test_pars <- pars
test_pars$A <- matrix(0.5, nrow = test_pars$S, ncol = test_pars$S)
diag(test_pars$A) <- 1
test_B <- GetFeasibleB(rep(1, times = pars$S), test_pars)
test_B[test_B != 0] <- 1

array_B <- array(test_B, c(pars$S, pars$S, pars$S))
long_inds <- which(array_B == 1, arr.ind = TRUE)
long_inds <- as.data.frame(long_inds)
long_inds <- long_inds[order(long_inds[,1]),]
long_inds$IntID <- 1:nrow(long_inds)

num_param <- pars$S * choose(pars$S - 1, 2)
num_const <- pars$S + 2 * num_param + pars$S * (pars$S - 1)

constr_mat <- matrix(0, nrow = num_const, ncol = num_param)
constr_dir <- rep("=", times = num_const)
constr_rhs <- rep(0, times = num_const)

# first the positivity constraints

sign_constr <- as.vector(sign(total_Bs))
sign_constr <- rep(sign_constr, each = num_param / pars$S)

constr_mat[1:num_param,] <- - diag(sign_constr, nrow = num_param, ncol = num_param)
constr_dir[1:num_param] <- "<="

# now the feasibility constraints

feas_constr <- rep(diag(1, pars$S,  pars$S), each = num_param / pars$S)
feas_constr <- matrix(feas_constr, pars$S, num_param, byrow = TRUE)

constr_mat[(num_param + 1:pars$S),] <- feas_constr
constr_rhs[(num_param + 1:pars$S)] <- total_Bs

# last the correlation constraints

cur_count <- num_param + pars$S + 1

for(i in 1:pars$S) {
  js <- (1:pars$S)
  js <- js[js != i]
  
  for(j in js) {
    
    curA <- pars$A[i,j]
    curC <- C[i, j]
    
    cur_inds <- long_inds %>%
      filter(dim1 == i) %>%
      filter((dim2 == j) | (dim3 == j))
    
    cur_inds <- cur_inds$IntID
    constr_mat[cur_count, cur_inds] <- 1

    constr_rhs[cur_count] <- - curA / unique(target_abd) + curC
    
    cur_count <- cur_count + 1
  }
}

# putting upper bounds on the parameter search

constr_mat[cur_count:nrow(constr_mat),] <- - constr_mat[1:num_param,]
constr_dir[cur_count:length(constr_dir)] <- "<="
constr_rhs[cur_count:length(constr_dir)] <- rep(total_Bs, each = num_param / pars$S)

# putting it all together

constr_list <- list(constr = constr_mat, dir = constr_dir, rhs = constr_rhs)

out_B_vals <- hitandrun(constr = constr_list, n.samples = 1)

out_B <- array(0, dim = c(pars$S, pars$S, pars$S))
out_B[as.matrix(long_inds[,-4])] <- out_B_vals
out_B <- matrix(out_B, nrow = pars$S, ncol = pars$S^2)
out_B

pars$B <- out_B

pars$r - pars$A %*% target_abd - pars$B %*% (as.vector(outer(target_abd, target_abd)))
BuildJacobian(target_abd, pars)



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




BuildConstraintVector <- function(pars) {
  testA <- pars$A
  diag(testA) <- 0
  testA <- - testA
  out_sums <- rowSums(testA)
  return(out_sums)
}

GetSkewSymmetry <- function(in_mat) {
  out_mat <- in_mat + t(in_mat)
  out_mat <- abs(out_mat)
  out_mat <- sum(out_mat)
  return(out_mat)
}

RowSumConst <- function(in_mat, const_vec) {
  in_sums <- rowSums(in_mat)
  out_diffs <- in_sums - const_vec
  return(out_diffs)
}

OptSkewSymMat <- function(const_vec,
                          num_trials = 100,
                          pert_ratio = 0.05) {
  
  num_sp <- length(const_vec)
  cur_mat <- matrix(const_vec / (num_sp-1), nrow = num_sp, ncol = num_sp)
  diag(cur_mat) <- 0
  cur_skew <- GetSkewSymmetry(cur_mat)
  
  out_df <- data.frame(Step = 1:num_trials, SkewSymmetry = NA)
  
  for(cur_trial in 1:num_trials) {
    
    new_pert <- matrix(runif(num_sp^2), nrow = num_sp, ncol = num_sp)
    diag(new_pert) <- 0
    row_sum_mat <- matrix(rowSums(new_pert), nrow = num_sp, ncol = num_sp)
    diag(row_sum_mat) <- 0
    row_sum_mat <- row_sum_mat / (num_sp - 1)
    new_pert <- new_pert - row_sum_mat
    
    new_mat <- cur_mat + pert_ratio * new_pert
    
    new_skew <- GetSkewSymmetry(new_mat)
    
    if(new_skew < cur_skew) {
      cur_skew <- new_skew
      cur_mat <- new_mat
    }
    
    out_df$SkewSymmetry[cur_trial] <- cur_skew
    
  }
  
  out_list <- list(mat = cur_mat, steps = out_df)
  return(out_list)
}

const_vec <- BuildConstraintVector(pars)
#opt_out <- OptSkewSymMat(const_vec, num_trials = 1e3, pert_ratio = 0.1)
opt_out$mat
tail(opt_out$steps)

C <- opt_out$mat

C <- BuildA(S = pars$S, mu = 0, sigma = 0.5, rho = -1)
diag(C) <- 0
row_sum_mat <- matrix(const_vec - rowSums(C), nrow = pars$S, ncol = pars$S)
diag(row_sum_mat) <- 0
row_sum_mat <- row_sum_mat / (pars$S - 1)
C <- C - row_sum_mat

ggplot(opt_out$steps, aes(x = Step, y = SkewSymmetry)) +
  geom_line() + theme_classic()

### The commented out code can never work!!! Because skew symmetric matrices
### sum to zero overall, so there are constraints on the row sums they can
### achieve ie. the row sums must sum to zero...

#corr_vec <- BuildCorrelatedVector(pars)
#corr_vec


#local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-15 )
#opts <- list( "algorithm"= "NLOPT_GN_ISRES",
#              "xtol_rel"= 1.0e-15,
#              "maxeval"= 160000,
#              "local_opts" = local_opts,
#              "print_level" = 0 )

#const_vec <- BuildConstraintVector(pars)
#ini_vec <- rep(const_vec, times = pars$S)
#nloptr(x0 = ini_vec,
#       eval_f = GetSkewSymmetry,
#       eval_g_eq = RowSumConst,
#       opts = opts,
#       const_vec = const_vec)



BuildIndMat <- function(pars) {
  ind_mat <- matrix(0, pars$S, pars$S)
  ind_mat[upper.tri(ind_mat)] <- 1
  ind_lookup <- which(ind_mat == 1, arr.ind = TRUE)
  ind_lookup <- as.data.frame(ind_lookup)
  ind_lookup <- ind_lookup[order(ind_lookup$row),]
  ind_lookup$ID <- 1:nrow(ind_lookup)
  return(ind_lookup)
}

BuildCorrelatedMatrix <- function(pars) {
  
  ind_lookup <- BuildIndMat(pars)
  
  num_entries <- pars$S * (pars$S - 1) / 2
  corr_mat <- matrix(0, nrow = pars$S, ncol = num_entries)
  
  for(i in 1:pars$S) {
    
    cur_vals <- rep(0, times = num_entries)
    
    pos_inds <- ind_lookup$row == i
    cur_vals[pos_inds] <- 1
    
    neg_inds <- ind_lookup$col == i
    cur_vals[neg_inds] <- -1
    
    corr_mat[i,] <- cur_vals
  }
  
  return(corr_mat)
}


#corr_vec <- BuildCorrelatedVector(pars)
#corr_vec <- c(1, -1, 1, -1)

#corr_mat <- BuildCorrelatedMatrix(pars)
#corr_mat

#corr_inv <- ginv(corr_mat)
#C_vals <- corr_inv %*% corr_vec

#corr_mat %*% C_vals
#corr_vec

#num_param <- ncol(corr_mat)

#corr_mat <- rbind(corr_mat,
#                  diag(1, nrow = num_param, ncol = num_param),
#                  diag(-1, nrow = num_param, ncol = num_param))

#corr_rhs <- c(corr_vec, rep(10, times = num_param), rep(-10, times = num_param))

#corr_dir <- c(rep("=", times = pars$S), rep("<=", times = 2 * num_param))

#constr_list <- list(constr = corr_mat,  dir = corr_dir, rhs = corr_rhs)

#hitandrun(constr = constr_list, n.samples = 1)



####


#new_mat <- diag(nrow = nrow(corr_inv)) - corr_inv %*% corr_mat
#rand_sam <- runif(length(C_vals))
#C_vals <- C_vals + new_mat %*% rand_sam

#ind_lookup <- BuildIndMat(pars)
#C <- matrix(0, pars$S, pars$S)
#C[as.matrix(ind_lookup[,1:2])] <- C_vals

#C <- C - t(C)

#C
#rowSums(C)
#corr_vec

#GetCorrelatedC <- function(pars, const_mat, const_vec) {
#  corr_inv <- ginv(corr_mat)
#  C <- corr_inv %*% const_vec
#  
#  new_mat <- diag(nrow = nrow(const_inv)) - const_inv %*% const_mat
#  rand_sam <- runif(length(const_B))
#  const_B <- const_B + new_mat %*% rand_sam
#  
#  noise_added <- runif(length(const_B),
#                       min = - p * max(const_B),
#                       max = p * max(const_B))
#  
#  const_B <- Vec2WideMat(num_spec = pars$S, vec = const_B)
#  noise_added <- matrix(noise_added, nrow = nrow(const_B), ncol = ncol(const_B))
#  
#  noise_added <- noise_added - rowSums(noise_added) / ncol(noise_added)
#  const_B <- const_B + noise_added
#  
#  return(const_B)
#}



#const_B <- GetConstrainedB(target_abd, pars)
#pars$B <- const_B

#Bmat <- pars
#Bmat$A <- matrix(0, nrow = Bmat$S, ncol = Bmat$S)
#Bmat <- BuildJacobian(eq_abd = rep(1, times = Bmat$S), Bmat)


GetCorrelatedB <- function(target_abd, pars) {
  
  
  return(out_B)
}
