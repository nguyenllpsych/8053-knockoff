#########################
## STAT 8053 - Model-X ##
##  Simulation script  ##
##     Linh Nguyen     ##
##     2021-10-20      ##
#########################

# Set up ----
# load libraries
libraries <- c("MASS",       # mvrnorm
               "knockoff",   # model-X knockoff
               "tidyverse")  # general wrangling
lapply(libraries, library, character = TRUE)

# simulation conditions
ns   <- c(100, 200, 500)       # sample size
ps   <- c(200, 500, 1000)      # number of covariates
ks   <- c(0.1)                 # % of covar with non-zero coef/total covar
rhos <- c(0, 0.01, 0.1, 0.5)   # off-diag entries for Sigma
amp  <- 5                      # signal amplitude
fdr  <- 0.1                    # false discovery rate                   

# Simulation 1 ----
# set seed
set.seed(8053)

# simulation conditions
# first conditions as initial test
n   <- ns[1]       # sample size
p   <- ps[1]       # number of covariates
k   <- ks[1]       # % of covar with non-zero coef/total covar
rho <- rhos[1]     # off-diag entries for Sigma

# correlation matrix
Sigma       <- matrix(rho, nrow = p, ncol = p)
diag(Sigma) <- 1

# simulate covariates
X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

# simulate responses from linear model 
hit <- sample(x = p, size = p*k)
beta <- amp * (1:p %in% hit)/sqrt(n)
Y <- X %*% beta + rnorm(n) # linear model + error

# knock-off procedures
# somehow, don't know why, need line 49 for line 50 to run..... 
# even though line 49 produced error
selected <- knockoff.filter(X = X, y = Y, knockoffs = create.gaussian(X = X, mu = rep(0, p), Sigma = Sigma))
selected <- knockoff.filter(X = X, y = Y, fdr = fdr)$selected
selected
hit

# FDR 
sum(beta[selected] == 0) / max(1, length(selected))

# benjamini-hochberg
p.val <- vector(mode = "numeric", length = p)
for (x in seq_len(p)) {
  # fit individual linear models
  mod      <- lm(Y ~ X[, x])
  # store p values
  p.val[x] <- summary(mod)$coefficients[2, 4]
}

# sort from smallest to largest p values
p.val.sort <- data.frame(p.val = sort(p.val),
                         rank = 1:p)
p.val <- data.frame(id = 1:p,
                    p.val = p.val)
p.val <- merge(p.val, p.val.sort)

# calculate B-H critical values
p.val$BH <- (p.val$rank/200)*fdr

# B-H hits
BH <- p.val[which(p.val$BH < 0.05), ]$id

# FDR for B-H
sum(beta[BH] == 0) / max(1, length(BH))
