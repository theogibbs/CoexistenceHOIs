########################################################
# Non-equal abundances
########################################################

# Loading dependencies

source("./Functions.R")
source("./DiverseCommunities/DiverseCommunityFunctions.R")

# Simulating script

S <- 4
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

target_abd <- c(0.1, 0.2, 0.9, 0.8)

B <- GetFeasibleB(target_abd, pars)
pars$B <- B
out_Bs <- matrix(0, nrow = 0, ncol = pars$S^2)

for(k in 1:100) {
  
  out_B <- GetConstrainedB(target_abd, pars)
  # need to merge all these data into one reasonable format over the for loop
  out_Bs <- rbind(out_Bs, out_B)
}

# do a test of the three species simplex to check if this is working to get a uniform sample

tern_data <- matrix(0, nrow = nrow(out_Bs), ncol = choose(pars$S - 1, 2))

for(i in 1:nrow(out_Bs)) {
  curB <- out_Bs[i,]
  tern_data[i,] <- curB[curB != 0]
}

tern_data <- tern_data / (pars$r - as.vector(pars$A %*% target_abd))
tern_data <- as.data.frame(tern_data)
colnames(tern_data) <- c("first_hoi", "second_hoi", "third_hoi")

plSimplex <- ggtern(tern_data, aes(first_hoi, second_hoi, third_hoi)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  theme(text = element_blank(),
        #legend.position = "top",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size = 15),
        axis.line=element_blank(),
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
        plot.margin=unit(c(0,0,0,0),"mm"))

show(plSimplex)

jpeg("../CoexistenceHOIs-Paper/figs/SIFigSimplex.jpeg",
     width = 1000, height = 1000, res = 300)
plSimplex
dev.off()

