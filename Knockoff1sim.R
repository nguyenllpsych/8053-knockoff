library(knockoff)

set.seed(8053)
# Problem parameters
n = 500          # number of observations
p = 1000          # number of variables
k = 60            # number of variables with nonzero coefficients
amplitude = 4.5   # signal amplitude (for noise level = 1)

# Generate the variables from a multivariate normal distribution
mu = rep(0,p)
rho = 0.25
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n) %*% chol(Sigma)

# Generate the response from a linear model
nonzero = sample(p, k)
beta = amplitude * (1:p %in% nonzero) / sqrt(n)
y = X %*% beta + rnorm(n)

knockpower = function(selected) sum(beta[selected]>0) / k
fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))

fdr = 0
kp = 0
for (i in 1:100){ #run for 24 hours
  result1 = knockoff.filter(X, y, offset = 0)
  fdr = fdr + fdp(result1$selected)
  kp = kp + knockpower(result1$selected)
  print(paste0("sim ", i," having ",fdr/i, kp/i))
}


gaussian_knockoffs = function(X) create.gaussian(X, mu, Sigma)
result2 = knockoff.filter(X, y, knockoffs=gaussian_knockoffs)
result2

fdp(result2$selected)

result3 = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, statistic = stat.random_forest, fdr=0.2)
print(result3)
fdp(result3$selected)

#High demension BH, not working.
#lm.fit = lm(y ~ X - 1) # 
#p.values = coef(summary(lm.fit))[,4]
#cutoff = max(c(0, which(sort(p.values) <= 0.1 * (1:p) / p)))
#bhq_selected = names(which(p.values <= 0.1 * cutoff / p))
#fdp(bhq_selected)
