#########################
## STAT 8053 - Model-X ##
##   Analysis script   ##
##     Linh Nguyen     ##
##     2021-11-30      ##
#########################

#### ---- Set up ---- ####

# libraries
libraries <- c("dplyr",      # general wrangling
               "knockoff",   # model-X knockoff
               "doMC",       # suggest for knockoff
               "Matrix",     # required for knockoff
               "glmnet",     # required for knockoff
               "ridge")      # ridge regression 

lapply(libraries, library, character = TRUE)

###########################
## ONLY NON MISSING DATA ##
###########################

# data
data <- readRDS(file = "data_clean.RData")

# controls
fdr <- 0.1 # false discovery rate

# change all genotypes to factor to numeric
char_to_num <- function(x) as.numeric(as.factor(x))
data <- data %>% mutate_at(c(5:6019), char_to_num)

# X and Y
X <- as.matrix(data[, -c(1:4)], nrow = 224)
Y <- data[, 4]

#### ---- MODEL-X KNOCKOFF ---- ####

# set seed 
set.seed(8053)

# create gaussian function
gaussian_knockoffs <- function(X) create.second_order(X)

# knock-off procedures timed
ko_start <- Sys.time()
ko_selected  <- knockoff.filter(X = X, 
                                y = Y, 
                                fdr = fdr, 
                                knockoff = gaussian_knockoffs,
                                offset = 0)$selected
ko_end <- Sys.time()
ko_time <- ko_end - ko_start


#### ---- BENJAMINI-HOCHBERG ---- ####

# fit ridge regression and extract p-val
bh_start <- Sys.time()
mod <- summary(linearRidge(Y ~ X + 0))
mod.chosen <- mod$chosen.nPCs
mod.chosen <- mod$summaries[mod.chosen][[1]]
p.val <- mod.chosen$coefficients[, 5]

# B-H selected
bh_selected <- which(p.adjust(p.val, method = "BH") <= fdr)
bh_end <- Sys.time()
bh_time <- bh_end - bh_start

#### ---- BONFERRONI ---- ####

# alpha = 5e-8
alpha <- 5e-8
bon_selected <- vector(mode = "character")
ps <- colnames(X)

bon_start <- Sys.time()
for (p in 1:ncol(X)) {
  fstat <- summary(lm(Y ~ as.factor(X[, p])))$fstatistic
  pval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

  if(pval < alpha) {bon_selected <- append(bon_selected, ps[p])}
}
bon_end <- Sys.time()
bon_time <- bon_end - bon_start


######################
## 5% MISSING DATA ##
######################

# data
rm(data)
gc()
data_5 <- readRDS(file = "data_5.RData")

# controls
fdr <- 0.1 # false discovery rate

# change all genotypes to factor to numeric
data_5 <- data_5 %>% mutate_at(c(5:ncol(data_5)), char_to_num)

# X and Y
X <- as.matrix(data_5[, -c(1:4)], nrow = 224)
Y <- data_5[, 4]

#### ---- MODEL-X KNOCKOFF ---- ####

# set seed 
set.seed(8053)

# knock-off procedures timed
ko_start_5  <- Sys.time()
ko_selected_5  <- knockoff.filter(X = X, 
                                y = Y, 
                                fdr = fdr, 
                                knockoff = gaussian_knockoffs,
                                offset = 0)$selected
ko_end_5 <- Sys.time()
ko_time_5 <- ko_end_5 - ko_start_5

#### ---- BONFERRONI ---- ####

# alpha = 5e-8
alpha <- 5e-8
bon_selected_5 <- vector(mode = "character")
ps <- colnames(X)

bon_start_5 <- Sys.time()
for (p in 1:ncol(X)) {
  fstat <- summary(lm(Y ~ as.factor(X[, p])))$fstatistic
  pval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

  if(pval < alpha) {bon_selected_5 <- append(bon_selected_5, ps[p])}
}
bon_end_5 <- Sys.time()
bon_time_5 <- bon_end_5 - bon_start_5

saveRDS(object = c(ko_time, ko_selected, bon_time, bon_selected), 
        file = "02-results.RData")

