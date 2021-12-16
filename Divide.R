libraries <- c("knockoff",   # model-X knockoff
               "doMC",       # suggest for knockoff
               "Matrix",     # required for knockoff
               "glmnet",     # required for knockoff
               "ridge")      # ridge regression in n > p cases
lapply(libraries, library, character = TRUE)

set.seed(8053)
# Problem parameters
n = 500          # number of observations
p = 1000          # number of variables
k = 60            # number of variables with nonzero coefficients
amplitude = 5   # signal amplitude (for noise level = 1)

# Generate the variables from a multivariate normal distribution
mu = rep(0,p)
rho = 0.25
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n) %*% chol(Sigma)

# Generate the response from a linear model
nonzero = sample(p,k)


beta = amplitude * (1:p %in% nonzero) / sqrt(n)
Y = X %*% beta + rnorm(n)

knockpower = function(selected) sum(beta[selected]>0) / k
fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))


#split in advance
kt <- matrix(, nrow = 30, ncol = 6)
k_ds <- rep(list(c()),30)
k_os <- rep(list(c()),30)

for (j in 1:30){
  start_time <- Sys.time()
  k_selected <- knockoff.filter(X,Y,offset = 0)
  end_time <- Sys.time()
  kt[j,6] <- end_time - start_time
  k_os[[j]] <- k_selected$selected
  
  for (i in 1:5){
    Xd <- X[,(200*i-199):(200*i)]
    start_time <- Sys.time()
    k_selected <- knockoff.filter(Xd,Y,offset=0)
    end_time <- Sys.time()
    kt[j,i] <- end_time - start_time
    k_ds[[j]] <- append(k_ds[[j]],k_selected$selected+200*(i-1))
    print(paste0("K Iteration ",j,' ',i," with ", k_ds[[j]]))
    gc()
  }
}

kdspower = matrix(, nrow = 30, ncol = 5)
kdsfdr = matrix(, nrow = 30, ncol = 5)
klist = c(16,12,12,12,8)
for (j in 1:5){
  for (i in 1:30){
    selected = k_ds[[i]][k_ds[[i]]< (200*j+1)&k_ds[[i]]>(200*j-200)]
    kdspower[i,j] = sum(beta[selected]>0) / klist[j]
    kdsfdr[i,j] = sum(beta[selected] == 0) / max(1, length(selected))
  }
}



kospower = 0
for (i in 1:30){
  kospower = kospower + knockpower(k_os[[i]])
}
kospower = kospower/30

kosfdr = 0
for (i in 1:30){
  kosfdr = kosfdr + fdp(k_os[[i]])
}
kosfdr = kosfdr/30

koscount = rep(0,30)
for (i in 1:30){
  koscount[i] = length(k_os[[i]])
}
mean(koscount)

kdscount = matrix(, nrow = 30, ncol = 5)
for (i in 1:30){
  kdscount[i,1] = sum(k_ds[[i]]<201)
  kdscount[i,2] = sum(k_ds[[i]]<401)-kdscount[i,1]
  kdscount[i,3] = sum(k_ds[[i]]<601)-kdscount[i,1]-kdscount[i,2]
  kdscount[i,4] = sum(k_ds[[i]]<801)-kdscount[i,1]-kdscount[i,2]-kdscount[i,3]
  kdscount[i,5] = sum(k_ds[[i]]<1001)-kdscount[i,1]-kdscount[i,2]-kdscount[i,3]-kdscount[i,4]
}
mean(kdscount)

kdslist = 1
for (i in 1:30){
  kdslist = append(kdslist,k_ds[[i]])
}
kdslist3 = sort(unique(kdslist))[which(table(kdslist)>12)]
fdp(kdslist3)
knockpower(kdslist3)

#random split
kt2 <- matrix(, nrow = 30, ncol = 6)
k_ds2 <- rep(list(c()),30)
k_os2 <- rep(list(c()),30)

for (j in 1:30){
  start_time <- Sys.time()
  k_selected <- knockoff.filter(X,Y,offset = 0)
  end_time <- Sys.time()
  kt2[j,6] <- end_time - start_time
  k_os2[[j]] <- k_selected$selected
  gc()
  
  index = sample(1000,1000)
  for (i in 1:5){
    split = index[(200*i-199):(200*i)]
    Xd <- X[,split]
    start_time <- Sys.time()
    k_selected <- knockoff.filter(Xd,Y,offset=0)
    end_time <- Sys.time()
    kt2[j,i] <- end_time - start_time
    realSelected <- split[k_selected$selected]
    k_ds2[[j]] <- append(k_ds2[[j]],realSelected)
    print(paste0("K Iteration ",j,' ',i," with ", k_ds2[[j]]))
    gc()
  }
}

kds2fdr1 <- rep(0,30)
kds2power1 <- rep(0,30)
kds2count <- rep(0,30)
kds2list = c()
kos2fdr1 <- rep(0,30)
kos2power1 <- rep(0,30)
kos2count <- rep(0,30)
kos2list = c()
for (i in 1:30){
  kds2fdr1[i] = fdp(k_ds2[[i]])
  kds2power1[i] = knockpower(k_ds2[[i]])
  kds2count[i] = length(k_ds2[[i]])
  kds2list = append(kds2list,k_ds2[[i]])
  kos2fdr1[i] = fdp(k_os2[[i]])
  kos2power1[i] = knockpower(k_os2[[i]])
  kos2count[i] = length(k_os2[[i]])
  kos2list = append(kos2list,k_os2[[i]])
}


kds2list2 <- rep(list(c()),30)
kds2fdr2 <- rep(0,30)
kds2power2 <- rep(0,30)
kds2count2 <- rep(0,30)
for (i in 1:30){
  kds2list2[[i]] = sort(unique(kds2list))[which(table(kds2list)>i-1)]
  kds2fdr2[i] <- fdp(kds2list2[[i]])
  kds2power2[i] <- knockpower(kds2list2[[i]])
  kds2count2[i] <- length(kds2list2[[i]])
}
