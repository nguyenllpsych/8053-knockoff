library(knockoff)
#install.packages("doMC", repos="http://R-Forge.R-project.org")
library(dplyr)
library(stringr)

n = 1000
p = 400
k = 40
rep = 100

norm_fdp = matrix(NA,rep,4)
norm_pow = matrix(NA,rep,4)
exp_fdp = matrix(NA,rep,4)
exp_pow = matrix(NA,rep,4)
pois_fdp = matrix(NA,rep,4)
pois_pow = matrix(NA,rep,4)
logis_fdp = matrix(NA,rep,4)
logis_pow = matrix(NA,rep,4)
cauchy_fdp = matrix(NA,rep,4)
cauchy_pow = matrix(NA,rep,4)

i=1
while(i<=100){

ind = sample(p,k)
beta = rep(0,p)
beta[ind] = 4.5 / sqrt(n) * sample(c(1,-1),k , replace=T)
fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
power = function(x) length(intersect(ind,x))/length(ind)

### normal
X = matrix(rnorm(n*p),n)
y = X%*%beta + rnorm(n)

fit = lm(y~X+0)
p_values = summary(fit)$coefficient[,4]
bh = which(p.adjust(p_values, method="fdr") <= 0.1)
norm_fdp[i,1] = fdp(bh)
norm_pow[i,1] = power(bh)

result = knockoff.filter(X, y, offset=0)
norm_fdp[i,2] = fdp(result$selected)
norm_pow[i,2] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='sdp', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
norm_fdp[i,3] = fdp(result$selected)
norm_pow[i,3] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='equi', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
norm_fdp[i,4] = fdp(result$selected)
norm_pow[i,4] = power(result$selected)

### exponential
X = matrix(rexp(n*p),n)
y = X%*%beta + rnorm(n)

fit = lm(y~X+0)
p_values = summary(fit)$coefficient[,4]
bh = which(p.adjust(p_values, method="fdr") <= 0.1)
exp_fdp[i,1] = fdp(bh)
exp_pow[i,1] = power(bh)

result = knockoff.filter(X, y, offset=0)
exp_fdp[i,2] = fdp(result$selected)
exp_pow[i,2] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='sdp', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
exp_fdp[i,3] = fdp(result$selected)
exp_pow[i,3] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='equi', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
exp_fdp[i,4] = fdp(result$selected)
exp_pow[i,4] = power(result$selected)

### poisson
X = matrix(rexp(n*p),n)
y = X%*%beta + rnorm(n)

fit = lm(y~X+0)
p_values = summary(fit)$coefficient[,4]
bh = which(p.adjust(p_values, method="fdr") <= 0.1)
pois_fdp[i,1] = fdp(bh)
pois_pow[i,1] = power(bh)

result = knockoff.filter(X, y, offset=0)
pois_fdp[i,2] = fdp(result$selected)
pois_pow[i,2] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='sdp', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
pois_fdp[i,3] = fdp(result$selected)
pois_pow[i,3] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='equi', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
pois_fdp[i,4] = fdp(result$selected)
pois_pow[i,4] = power(result$selected)

###logistic
X = matrix(rlogis(n*p,s=sqrt(3)/pi),n)
y = X%*%beta + rnorm(n)

fit = lm(y~X+0)
p_values = summary(fit)$coefficient[,4]
bh = which(p.adjust(p_values, method="fdr") <= 0.1)
logis_fdp[i,1] = fdp(bh)
logis_pow[i,1] = power(bh)

result = knockoff.filter(X, y, offset=0)
logis_fdp[i,2] = fdp(result$selected)
logis_pow[i,2] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='sdp', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
logis_fdp[i,3] = fdp(result$selected)
logis_pow[i,3] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='equi', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
logis_fdp[i,4] = fdp(result$selected)
logis_pow[i,4] = power(result$selected)

###cauchy
X = matrix(rcauchy(n*p),n)
y = X%*%beta + rnorm(n)

fit = lm(y~X+0)
p_values = summary(fit)$coefficient[,4]
bh = which(p.adjust(p_values, method="fdr") <= 0.1)
cauchy_fdp[i,1] = fdp(bh)
cauchy_pow[i,1] = power(bh)

result = knockoff.filter(X, y, offset=0)
cauchy_fdp[i,2] = fdp(result$selected)
cauchy_pow[i,2] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='sdp', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
cauchy_fdp[i,3] = fdp(result$selected)
cauchy_pow[i,3] = power(result$selected)

gaussian_knockoffs = function(X) create.second_order(X, method='equi', shrink=T)
result = knockoff.filter(X, y, knockoffs = gaussian_knockoffs, offset=0)
cauchy_fdp[i,4] = fdp(result$selected)
cauchy_pow[i,4] = power(result$selected)

print(i)
i= i+1
}

fdp = matrix(NA,5,4)
pow = matrix(NA,5,4)
rownames(fdp) = c("normal","exponential",'poisson','logistic','cauchy')
rownames(pow) = c("normal","exponential",'poisson','logistic','cauchy')
colnames(fdp) = c('BH','default','sdp','equi')
colnames(pow) = c('BH','default','sdp','equi')

fdp[1,] = apply(norm_fdp,2,mean)
fdp[2,] = apply(exp_fdp,2,mean)
fdp[3,] = apply(pois_fdp,2,mean)
fdp[4,] = apply(logis_fdp,2,mean)
fdp[5,] = apply(cauchy_fdp,2,mean)

pow[1,] = apply(norm_pow,2,mean)
pow[2,] = apply(exp_pow,2,mean)
pow[3,] = apply(pois_pow,2,mean)
pow[4,] = apply(logis_pow,2,mean)
pow[5,] = apply(cauchy_pow,2,mean)

# (procedure predict nonzero and it is actually non zero coef) / (total number of nonzero )
# p-values are exact (correct)
# p-values are asymptotic (knockoffs)
