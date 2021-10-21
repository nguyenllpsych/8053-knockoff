#########################
## STAT 8053 - Model-X ##
##  Simulation script  ##
##     Linh Nguyen     ##
##     2021-10-20      ##
#########################

#### ---- Set up ---- ####

# load libraries
libraries <- c("MASS",       # mvrnorm
               "knockoff",   # model-X knockoff
               "tidyverse")  # general wrangling
lapply(libraries, library, character = TRUE)

#### ---- Function definitions ####

# function for Benjamini-Hochberg
bh <- function(X, Y, fdr, beta) {
  
  # arguments:
  # X = nxp matrix of covariates
  # Y = nx1 vector of responses
  # fdr = numeric false discovery rate
  # beta = px1 vector of true beta coefficients in linear model
  
  # vector to store p.val
  p.val <- vector(mode = "numeric", length = ncol(X))
  
  for (x in seq_len(ncol(X))) {
    
    # fit individual linear models
    mod      <- lm(Y ~ X[, x])
    
    # store p values
    p.val[x] <- summary(mod)$coefficients[2, 4]
  }
  
  # sort from smallest to largest p values
  p.val.sort <- data.frame(p.val = sort(p.val),
                           rank = 1:ncol(X))
  p.val <- data.frame(id = 1:ncol(X),
                      p.val = p.val)
  p.val <- merge(p.val, p.val.sort)
  
  # calculate B-H critical values
  p.val$BH <- (p.val$rank/ncol(X))*fdr
  
  # B-H hits
  bh_selected <- p.val[which(p.val$BH < 0.05), ]$id
  
  # FDR for B-H
  bh_fdr <- sum(beta[bh_selected] == 0) / max(1, length(bh_selected))
  
  # return list of hits and fdr
  return(list(bh_selected = bh_selected,
              bh_fdr      = bh_fdr))
}

# function for knock off procedures
ko_gaussian <- function(ns, ps, rhos, 
                        k, amp, fdr,
                        seed) {
  # arguments:
  # ns = vector of different sample sizes
  # ps = vector of different numbers of covariates
  # rhos = vector of different off-diag entries for Sigma
  # k = numeric proportion of covariates with non-zero coefficient
  # amp = numeric signal amplitude for noise = 1
  # fdr = numeric false discovery rate
  # seed = seed
  
  # returns:
  # list of length (simulation times)
  # results of knock-off procedures and Benjamini-Hochberg procedures
  # including selected variables, fdr, run time
  # using same simulated data set
  
  # set seed
  set.seed(seed)
  
  # initialize return list
  results <- vector(mode = "list", length = length(ns)*length(ps)*length(rhos))
  i <- 0
  
  # create gaussian function
  gaussian_knockoffs <- function(X) create.gaussian(X, mu, Sigma)

  # looping through all simulation conditions
  for (n in ns) {
    for (p in ps) {
      for (rho in rhos) {
        
        ### SIMULATE DATA ###
        
        # simulation count
        i = i + 1
        
        # correlation matrix
        Sigma       <- matrix(rho, nrow = p, ncol = p)
        diag(Sigma) <- 1

        # simulate covariates
        mu <- rep(0, p)
        X <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
        
        # simulate responses from linear model 
        hit <- sample(x = p, size = p*k)
        beta <- amp * (1:p %in% hit)/sqrt(n)
        Y <- X %*% beta + rnorm(n) # linear model + error

        ### MODEL-X KNOCKOFF PROCEDURES ###
        ko_start <- Sys.time()
        ko_selected <- knockoff.filter(X = X, y = Y, fdr = fdr, 
                                       knockoff = gaussian_knockoffs)$selected
        ko_end   <- Sys.time()
        ko_fdr   <- sum(beta[ko_selected] == 0) / max(1, length(ko_selected))
        ko_miss  <- (length(hit) - sum(hit %in% ko_selected))/length(hit)
        ko_time  <- ko_end - ko_start

        ### BENJAMINI-HOCHBERG PROCEDURES ###
        bh_start    <- Sys.time()
        bh_results  <- bh(X = X, Y = Y, fdr = fdr, beta = beta)
        bh_end      <- Sys.time()
        bh_selected <- bh_results$bh_selected
        bh_fdr      <- bh_results$bh_fdr
        bh_miss     <- (length(hit) - sum(hit %in% bh_selected))/length(hit)
        bh_time     <- bh_end - bh_start
        
        # return results
        results[[i]] <- list(n   = n, 
                             p   = p,
                             rho = rho, 
                             k   = k,
                             amp = amp,
                             fdr = fdr,
                             ko_fdr  = ko_fdr,
                             bh_fdr  = bh_fdr,
                             ko_miss = ko_miss,
                             bh_miss = bh_miss,
                             hit     = sort(hit),
                             ko_selected = ko_selected,
                             bh_selected = bh_selected,
                             ko_time     = ko_time,
                             bh_time     = bh_time)
      }
    }
  } # END for ns, ps, rhos LOOPS
  
  # return full simulation results
  return(results)
}

#### ---- Simulation ---- ####

# simulation conditions
ns   <- c(100, 200)            # sample size
ps   <- c(200, 500)            # number of covariates
rhos <- c(0.01, 0.1, 0.5)      # off-diag entries for Sigma
k    <- c(0.1)                 # % of covar with non-zero coef/total covar
amp  <- 5                      # signal amplitude
fdr  <- 0.1                    # false discovery rate       

# simulation results
sim_results <- ko_gaussian(ns = ns, ps = ps, rhos = rhos, 
                           k = k, amp = amp, fdr = fdr,
                           seed = 8053)
sim_results[[1]]
