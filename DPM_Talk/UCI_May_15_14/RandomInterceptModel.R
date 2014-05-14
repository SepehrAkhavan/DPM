# Random Intercept Model:
# Author: Sepehr Akhavan 

# Model:
# Y_i = b_0i + Beta_1 * T + e 
    # -- where:
    # 1) m_i: number of measurements for subject i
    # 2) Y_i: a vector of length m_i of albumin measures
    # 3) T: a vector of time 
    # 4) Beta_1: a common covariate across all subjects
    # 5) e ~ N_{m_i}(0, sigma2 = diag(m_i)) --> Inverse gamma on sigma2

# nSub =  30
# b0.unique = c(-5, 0, 5) 
# b0 | G ~ G
# G ~ DP(alpha, G0 = N(0, sd = 10)) ---> we fix alpha to 1.5 for now


#-------------
# Libraries --
#-------------
library(mvtnorm)

#-------------------
# Simulating Data --
#-------------------
set.seed(1234)
nSub <- 30
b0.true <- rep(c(-5, 0, 5), each = nSub/3)
B1.true <- 1
eps.sigma2 <- 0.1
m = sample(5:10, nSub, replace = TRUE)

# data: for the case with different m_i, the dataset SHOULD BE in LONG format -- each row for a subject/time
data <- data.frame(id = rep(NA, sum(m)), b0 = rep(NA, sum(m)), t = rep(NA, sum(m)), eps = rep(NA, sum(m)), albumin = rep(NA, sum(m)))
id.tmp <- b0.tmp <- t.tmp <- albumin.tmp <- eps.tmp <- NULL

for (i in 1:nSub){
  id.tmp <- c(id.tmp, rep(i, m[i]))
  b0.tmp <- c(b0.tmp, rep(b0.true[i], m[i]))
  t.tmp <- c(t.tmp, sort(runif(m[i], 0, 30)))
  eps.tmp <- c(eps.tmp, as.vector(rmvnorm(1, rep(0, m[i]), eps.sigma2*diag(m[i]))))
}
data$id <- id.tmp
data$b0 <- b0.tmp
data$t <- t.tmp
data$eps <- eps.tmp
data$albumin <- data$b0 + B1.true*data$t + data$eps

# remove additional columns (only leave in what we have in reality)
data <- data[,c("id", "t", "albumin")]

#-------------------------------------------------
# Estimating the model via the Bayesian approach -
#-------------------------------------------------
# prior on Beta ~ N_1(mean = 0, sd = 2)
# prior on b0 ~ DP(alpha = 1.5, G0 = N(G0_mu = 0, G0_sd = 15))

# In sampling from the posterior distribution, we use Gibbs sampling:
# getPosterior(): evaluates the posterior distribution up to some constants:

# Given b0_i, we know that: Y_ij ~ N(b0_i + Beta1 T_ij, sigma2 = eps.sigma2)
# Measurements for each subject (given a subject) are all iid
# different subjects are independent
getLogPosterior <- function(Beta1, b0Vec, data, nSub){ # conditional posterior (.|bo)
  b0.tmp <- NULL
  for (i in 1:(length(b0Vec))){
    b0.tmp <- c(b0.tmp, rep(b0Vec[i], m[i])) # m is the num measurements per id !
  }
  logLik <- sum((-1/(2*eps.sigma2)) * (data$albumin - b0.tmp - Beta1*data$t)^2) 
  logPrior.Beta <- (-1/2*(Beta1.sigma2))*(Beta1 - Beta1.mu)^2
  return( logLik + logPrior.Beta )
}

# getBeta(): uses MH algorithm to update Beta's
getBeta <- function(Beta1, b0Vec){ # Beta1 and b0Vec are current values
  newParam <- rnorm(1, mean = Beta1, sd = Beta1.MH.stepSize)
  acc.prob <- exp(getLogPosterior(newParam, b0Vec, data, nSub) - getLogPosterior(Beta1, b0Vec, data, nSub))
  
  if (runif(1) < acc.prob){
    return(list(newParam = newParam, accept = 1))
  }else{
    return(list(newParam = Beta1, accept = 0))
  }
}

# This function uses Neal Alg.1 to sample from the posterior of our random intercept
get.b0 <- function(Beta1, b0Vec, nSub, data){ # Beta1 and b0Vec are current values
  for (i in 1:nSub){
    n.i <- nrow(data[data$id == i,])
    b0.j <- b0Vec[-i]
    q.ij.tmp <- NULL # q.ij = b F(Y_i, beta0_j) note that here Y_i is a vector ! MultiVariateNormal with diagonal variance
    for (k in 1:length(b0.j)){
      mean.tmp <- b0.j[k] + Beta1*data$t[data$id == i]
      dens.tmp <- dmvnorm(data$albumin[data$id == i], mean.tmp, eps.sigma2*diag(n.i), log=FALSE)
      q.ij.tmp <- c(q.ij.tmp, dens.tmp)
    }
    zbar <- sum(data$albumin[data$id == i] - Beta1*data$t[data$id == i])/n.i
    ri.tmp <- DP.alpha*(2*eps.sigma2*G0.sigma2/n.i)*dnorm(zbar, mean = G0.mu, sd = sqrt(G0.sigma2 + eps.sigma2/n.i))
    # ri.tmp <- alpha*dnorm(y[i], mean = Xmat[i,]%*%currentBeta. sd = sqrt(1/sigma2 + 1/sigma_0.2))
    P <- c(ri.tmp, q.ij.tmp)/(sum(c(ri.tmp, q.ij.tmp)))
    # we sample a new b0 here we decide whether to choose it or not later:
    # sample from Hi
    mean_new_b0 <- (n.i/(2*eps.sigma2))/(n.i/(2*eps.sigma2) + 1/(2*G0.sigma2)) * zbar + (eps.sigma2/(n.i * G0.sigma2 + eps.sigma2)) * G0.mu
    sd_new_b0 <- 2*eps.sigma2*G0.sigma2/(n.i * G0.sigma2 + eps.sigma2)
    new_b0 <- rnorm(1, mean = mean_new_b0, sd = sd_new_b0)
    b0.pool <- c(new_b0, b0.j)
    b0Vec[i] <- sample(b0.pool, 1, prob = P) # here we decide which b0 to choose
  }
  return(b0Vec)
}

#------------
# Remixing --
#------------
Remixing.get.b0 <- function(remixingData, b0.star.i, Beta1){
  zbar <- sum(remixingData$albumin - Beta1*remixingData$t)/nrow(remixingData)
  mean_new_b0 <- (nrow(remixingData)/(2*eps.sigma2))/(nrow(remixingData)/(2*eps.sigma2) + 1/(2*G0.sigma2)) * zbar + (eps.sigma2/(nrow(remixingData) * G0.sigma2 + eps.sigma2)) * G0.mu
  sd_new_b0 <- 2*eps.sigma2*G0.sigma2/(nrow(remixingData) * G0.sigma2 + eps.sigma2)
  new_b0 <- rnorm(1, mean = mean_new_b0, sd = sd_new_b0)
  return(new_b0)
}

remixingEngine <- function(Beta1, b0Vec, data){
  b0.star <- unique(b0Vec)
  indSet <- NULL
  for (i in 1:(length(b0Vec))){
    indSet <- c(indSet, which(b0.star == b0Vec[i]))
  }
  for (i in 1:(length(b0.star))){
    remixingData <- data[data$id %in% which(indSet == i),]
    b0.star[i] <- Remixing.get.b0(remixingData, b0.star[i], Beta1)
  }
  b0Vec <- b0.star[indSet]
  return(b0Vec) 
}

#---------------------------------------------------------------------
# MainEngine(): this is our main engine in running our Bayesian code -
#---------------------------------------------------------------------
mainEngine <- function(nIter){
    
  # place holder for future MCMC samples:
  post.Beta1.samp <- rep(NA, nIter)
  post.b0.samp <- matrix(rep(NA, nSub*nIter), ncol = nSub)
  
  # Acceptance Rate - for parameters with MH:
  Beta1.acc <- 0
  
  # initializing parameters:
  Beta1 <- 0
  b0 <- rep(1, nSub)
  
  for(iter in 1:nIter){
    print(iter)
    # First update b0:
    b0 <- get.b0(Beta1, b0, nSub, data)
    
    # remixing:
    b0 <- remixingEngine(Beta1, b0, data)
    
    post.b0.samp[iter,] <- b0
    
    # Then update Beta:
    Beta1.post.tmp <- getBeta(Beta1, b0)
    Beta1.acc <- Beta1.acc + Beta1.post.tmp$accept
    Beta1 <- Beta1.post.tmp$newParam
    post.Beta1.samp[iter] <- Beta1
  }
  
  # returning all samples
  return(list(b0 = post.b0.samp, Beta1 = post.Beta1.samp, Beta1.acc = round(Beta1.acc/nIter, 2)))
}

#---------------------------
# Time to run the code :) --
#---------------------------
Beta1.mu <- 0
Beta1.sigma2 <- 2
Beta1.MH.stepSize <- 0.01
# eps.sigma2: defined earlier when I simulate data
G0.mu <- 0
G0.sigma2 <- 10 # variance of the G0 dist
DP.alpha <- 1.5

start.time <- Sys.time()
rslt <- mainEngine(nIter = 5000)
Sys.time() - start.time
