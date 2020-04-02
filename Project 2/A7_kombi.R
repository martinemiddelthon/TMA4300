# attempting combined block alogrithm, with alternating blocks. 

library(boot)
library(ggplot2)
library(invgamma)
library(tidyverse)
library(ggpubr)

data(coal)

block.1 <- function(theta,i,sigma,t0,t2){
  # block 1 update
  th <- theta
  
  #help vars:
  t1 <- th["t1",i-1]
  y1 <- th["y1",i-1]
  y0 <- th["y0",i-1]
  beta <- th["beta",i-1]
  
  # keep beta constant
  th["beta",i] <- th["beta",i-1]
  
  # sample new t1 from normal distribution
  t1.n <- rnorm(1,t1,sigma)
  if(t1.n > t2 | t1.n < t0){
    t1.n <- t1
  }
  
  #find corresponding y0.n and y1.n
  y0.n <- max(which(coal$date <= t1.n)) - 1
  y1.n <- length(coal$date) - y0.n - 2
  
  # sample lambda0 from Gamma distirbution with new t1
  l0.n <- rgamma(1,shape = y0.n + 2, scale = 1/(t1.n - t0 + 1/beta))
  
  # sample lambda1 from Gamma distribution with new t1
  l1.n <- rgamma(1,shape = y1.n + 2, scale = 1/(t2 - t1.n + 1/beta))
  
  # find alpha
  g.fact <- sum(log(1:y0.n + 1)) + sum(log(1:y1.n + 1)) - sum(log(1:y0+1)) - sum(log(1:y1+1))
  
  t.fact <- (y0 + 2)*log(t1 - t0 + 1/beta) + 
    (y1 + 2)*log(t2 - t1 + 1/beta) - 
    (y0.n + 2)*log(t1.n - t0 + 1/beta) - 
    (y1.n + 2)*log(t2 - t1.n + 1/beta)
  
  alpha <- exp(t.fact + g.fact)
  alpha <- min(1,alpha)

  
  # accept or reject
  u <- runif(1)
  if(u < alpha){
    # accept
    th["t1",i] <- t1.n
    th["l0",i] <- l0.n
    th["l1",i] <- l1.n
    th["y0",i] <- y0.n
    th["y1",i] <- y1.n
    th["accept",i] <- 1
  }
  else{
    # reject
    th["t1",i] <- th["t1",i-1]
    th["l0",i] <- th["l0",i-1]
    th["l1",i] <- th["l1",i-1]
    th["y0",i] <- th["y0",i-1]
    th["y1",i] <- th["y1",i-1]
    th["accept",i] <- 0
  }
  return(th)
}

block.2 <- function(theta,i,sigma, t0, t2){
  # block 2 update
  th <- theta
  
  #help vars:
  t1 <- th["t1",i-1]
  y1 <- th["y1",i-1]
  y0 <- th["y0",i-1]
  beta <- th["beta",i-1]
  
  # keep t1 constant
  th["t1",i] <- th["t1",i-1]
  th["y0",i] <- th["y0",i-1]
  th["y1",i] <- th["y1",i-1]
  
  # sample beta from normal distribution
  beta.n <- rnorm(1,beta,sigma)
  if(beta.n <= 0){
    beta.n <- beta  # if beta out of bounds, reject 
  }
  
  # sample l0 from Gamma distribution with new beta
  l0.n <- rgamma(1,shape = y0 + 2, scale = 1/(t1 - t0 + 1/beta.n))
  
  # sample l1 from Gamma distribution with new beta
  l1.n <- rgamma(1,shape = y1 + 2, scale = 1/(t2 - t1 + 1/beta.n))
  
  # find alpha
  #alpha <- (beta^5*(t1 - t0 + 1/beta)^(y0 + 2)*(t2 - t1 + 1/beta)^(y1 + 2))/(beta.n^5^(t1 - t0 + 1/beta.n)^(y0 + 2)*(t2 - t1 + 1/beta.n)^(y1 + 2))
  
  alpha.lg <- 5*log(beta) - 5*log(beta.n) + 
    (y0 + 2)*log(t1 - t0 + 1/beta) + 
    (y1 + 2)*log(t2 - t1 + 1/beta) - 
    (y0 + 2)*log(t1 - t0 + 1/beta.n) - 
    (y1 + 2)*log(t2 - t1 + 1/beta.n)
  alpha <- exp(alpha.lg)
  alpha <- min(1, alpha)
  
  #accept or reject
  u <- runif(1)
  if(u < alpha){
    # accept
    th["beta",i] <- beta.n
    th["l0",i] <- l0.n
    th["l1",i] <- l1.n
    th["accept",i] <- 1
  }
  else{
    # reject
    th["beta",i] <- th["beta",i-1]
    th["l0",i] <- th["l0",i-1]
    th["l1",i] <- th["l1",i-1]
    th["accept",i] <- 0
  }
  return(th)
}

MH <- function(n,t1.0,l0.0,l1.0,beta.0, obs, sigma1, sigma2){
  # initializing parameters with prior estimates
  
  theta <- matrix(0L, nrow=7, ncol=n)
  rownames(theta) <- c("t1", "l0", "l1", "beta", "y0", "y1", "accept")
  
  t0 <- coal$date[1]
  t2 <- tail(coal$date,1)
  
  theta["t1",1] <- t1.0
  theta["y0",1] <- max(which(coal$date < theta["t1",1])) - 1 # disasters before t1, subtracting the first element which is not a disaster
  theta["y1",1] <- length(coal$date) - theta["y0",1] - 2 # disasters at and after t1, subtracting start time and end time element
  
  theta["l0",1] <- l0.0
  theta["l1",1] <- l1.0
  theta["beta",1] <- beta.0
  
  for(i in 2:n){
    # alternating block updates
    if(i%%2 == 0){
      # block 1 update
      theta <- block.1(theta, i, sigma1, t0, t2)
    }
    else{
      # block 2 update
      theta <- block.2(theta, i, sigma2, t0, t2)
    }
  }
  return(theta)
}


test <- MH(n=100000, t1.0=1900,l0.0 = 3,l1.0 = 1,beta.0  = 1, obs = coal, sigma1 = 10, sigma2 = 1)
test.df <- as.data.frame(t(test))

acceptance <- mean(test.df$accept)
acceptance

relevant.plots(test.df,coal,10000)
