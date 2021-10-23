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
               "doMC",       # suggest for knockoff
               "Matrix",     # required for knockoff
               "glmnet",     # required for knockoff
               "ridge")      # ridge regression in n > p cases
lapply(libraries, library, character = TRUE)

#### ---- Function definitions ---- ####

# function for Benjamini-Hochberg
bh <- function(X, Y, fdr, beta) {
  
  # arguments:
  # X = nxp matrix of covariates
  # Y = nx1 vector of responses
  # fdr = numeric false discovery rate
  # beta = px1 vector of true beta coefficients in linear model

  # if n > p
  if (ncol(X) < nrow(X)) {
    
    # fit multiple regression and extract p-val
    mod <- lm(Y ~ X + 0)
    p.val <- summary(mod)$coefficient[, 4]
  
  } else {
    
    # fit ridge regression and extract p-val
    mod <- summary(linearRidge(Y ~ X + 0))
    mod.chosen <- mod$chosen.nPCs
    mod.chosen <- mod$summaries[mod.chosen][[1]]
    p.val <- mod.chosen$coefficients[, 5]
  }

  # B-H hits  
  bh_selected <- which(p.adjust(p.val, method = "BH") <= fdr)

  # FDR for B-H
  bh_fdr <- sum(beta[bh_selected] == 0) / max(1, length(bh_selected))

  # return list of hits and fdr
  return(list(bh_selected = bh_selected,
              bh_fdr      = bh_fdr))
}

# function for knock off procedures
ko_gaussian <- function(ns, ps, 
                        k, amp, fdr,
                        reps, seed) {
  # arguments:
  # ns = vector of different sample sizes
  # ps = vector of different numbers of covariates
  # rhos = vector of different off-diag entries for Sigma
  # k = numeric proportion of covariates with non-zero coefficient
  # amp = numeric signal amplitude for noise = 1
  # fdr = numeric false discovery rate
  # reps = number of replications
  # seed = seed
  
  # returns:
  # list of length (simulation times)
  # results of knock-off procedures and Benjamini-Hochberg procedures
  # including selected variables, fdr, run time
  # using same simulated data set
  
  # set seed
  set.seed(seed)
  
  # initialize return list
  results <- vector(mode = "list", length = length(ns)*length(ps))
  
  # create gaussian function
  gaussian_knockoffs <- function(X) create.second_order(X)

  # initialize simulation count and result list
  i <- 0
  
  # looping through all simulation conditions
  for (n in ns) {
    for (p in ps) {
      
      # simulation count
      i = i + 1
      ko_time <- bh_time <- vector(mode = "numeric", length = reps)
      ko_fdr <- bh_fdr <- vector(mode = "numeric", length = reps)
      ko_miss <- bh_miss <- vector(mode = "numeric", length = reps)
      
      # progress bar
      pb <- txtProgressBar(max   = reps,
                           char  = paste0("sim", i),
                           style = 3)
      
      for (rep in seq_len(reps)){
         
        ### SIMULATE DATA ###
        
        # progress bar
        setTxtProgressBar(pb    = pb,
                          value = rep)
        
        # simulate covariates
        X <- matrix(rnorm(n*p), nrow = n)
        Sigma <- cor(X)
        
        # simulate responses from linear model 
        hit <- sample(x = p, size = p*k)
        beta <- amp * (1:p %in% hit)/sqrt(n)
        Y <- X %*% beta + rnorm(n) # linear model + error
        
        ### MODEL-X KNOCKOFF PROCEDURES ###
        ko_start     <- Sys.time()
        ko_selected  <- knockoff.filter(X = X, y = Y, fdr = fdr, 
                                        knockoff = gaussian_knockoffs,
                                        offset = 0)$selected
        ko_end       <- Sys.time()
        ko_fdr[rep]  <- sum(beta[ko_selected]==0)/max(1, length(ko_selected))
        ko_miss[rep] <- (length(hit) - sum(hit %in% ko_selected))/length(hit)
        ko_time[rep] <- ko_end - ko_start
        
        ### BENJAMINI-HOCHBERG PROCEDURES ###
        bh_start     <- Sys.time()
        bh_results   <- bh(X = X, Y = Y, fdr = fdr, beta = beta)
        bh_end       <- Sys.time()
        bh_selected  <- bh_results$bh_selected
        bh_fdr[rep]  <- bh_results$bh_fdr
        bh_miss[rep] <- (length(hit) - sum(hit %in% bh_selected))/length(hit)
        bh_time[rep] <- bh_end - bh_start
      }
      
      # return results -> CHANGE TO AVERAGE ACROSS REPS
      results[[i]] <- list(n   = n, 
                           p   = p,
                           k   = k,
                           amp = amp,
                           fdr = fdr,
                           ko_fdr  = ko_fdr,
                           bh_fdr  = bh_fdr,
                           ko_miss = ko_miss,
                           bh_miss = bh_miss,
                           ko_time = ko_time,
                           bh_time = bh_time)
    }
  } # END for ns, ps, rhos LOOPS
  
  # return full simulation results
  return(results)
}

#### ---- Simulation p > n ---- ####

# simulation conditions
ns   <- c(400)                 # sample size
ps   <- c(1000)                # number of covariates
k    <- 0.1                    # % of covar with non-zero coef/total covar
amp  <- 5                      # signal amplitude
fdr  <- 0.1                    # false discovery rate      

# simulation results
sim_results_hidim <- ko_gaussian(ns = ns, ps = ps, 
                                 k = k, amp = amp, fdr = fdr,
                                 reps = 100, seed = 8053)[[1]]

#### ---- Simulation n > p ---- ####

# simulation conditions
ns   <- c(1000)                # sample size
ps   <- c(400)                 # number of covariates
k    <- 0.1                    # % of covar with non-zero coef/total covar
amp  <- 5                      # signal amplitude
fdr  <- 0.1                    # false discovery rate       

# simulation results
sim_results_fullrank <- ko_gaussian(ns = ns, ps = ps, 
                                    k = k, amp = amp, fdr = fdr,
                                    reps = 100, seed = 8053)[[1]]
