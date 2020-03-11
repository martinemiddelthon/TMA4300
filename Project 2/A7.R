# Metropolis-Hastings block algorithm 1: 

data(coal)

MH.block.1 <- function(n,t1.0,l0.0,l1.0,beta.0, obs, sigma, beta.int){
  # initiate:
  th <- matrix(0L, nrow=6, ncol=n)
  rownames(th) <- c("t1", "l0", "l1", "beta", "y0", "y1")
  
  t0 <- coal$date[1]
  t2 <- tail(coal$date,1)
  
  th["t1",1] <- t1.0
  th["y0",1] <- max(which(coal$date < th["t1",1])) - 1 # disasters before t1, subtracting the first element which is not a disaster
  th["y1",1] <- length(coal$date) - th["y0",1] - 2 # disasters at and after t1, subtracting start time and end time element
  
  th["l0",1] <- l0.0
  th["l1",1] <- l1.0
  th["beta",1] <- beta.0
  
  for(i in 2:n){
    if(i %% beta.int == 0){
      # update beta
      th["beta",i] <- rinvgamma(1,shape=4, scale=1/(th["l0",i-1]+ th["l1",i-1] + 1))
    }
    else{
      # keep beta constant 
      th["beta",i] <- th["beta",i-1]
    }
    # sample t1 by random walk
    y <- rnorm(1,th["t1",i-1],sigma)       # proposal for t1
    if (y > t2 | y < t0){
      th["t1",i] <- th["t1",i-1]
      th["beta",i] <- th["beta",i-1]    # not sure if we should keep the "reset" of beta here, tbc
      th["l0",i] <- th["l0",i-1]
      th["l1",i] <- th["l1",i-1]
    }
    else {
      y0_prop <- max(which(coal$date <= y)) - 1
      y1_prop <- length(coal$date) - y0_prop - 2
      u <- runif(1)        
      #lg.alpha <- min(0, (th["l1",i-1]-th["l0",i-1])*(y - th["t1",i-1]))  # log acceptance probability
      #alpha = exp(lg.alpha)                 # acceptance probability
      alpha <- min(1, th["l0",i-1]^(y0_prop - th["y0",i-1])
                   *th["l1",i-1]^(y1_prop  - th["y1",i-1])
                   *exp((th["l1",i-1]-th["l0",i-1])*(y - th["t1",i-1])))
      
      if (u < alpha){
        th["t1",i] <- y
      }
      else{
        th["t1",i] <- th["t1",i-1]
      }
      
      # update y0 and y1 with current value of t1
      th["y0",i] <- max(which(coal$date <= th["t1",i])) - 1
      th["y1",i] <- length(coal$date) - th["y0",i] - 2
      
      # sample l0 from gamma with current value for t1
      th["l0",i] <- rgamma(1,shape=th["y0",i] + 2, scale=1/(th["t1",i] - t0 + 1/th["beta",i]))
      #th["l0",i] <- rgamma(1,shape=th["y0",i-1] + 2, scale=th["t1",i-1] - t0 + 1/th["beta",i-1]) #DEBUG
      
      # sample l1 from gamma with current value for t1
      th["l1",i] <- rgamma(1,shape = th["y1",i] + 2, scale=1/(t2 - th["t1",i] + 1/th["beta",i]))
      #th["l1",i] <- rgamma(1,shape = th["y1",i-1] + 2, scale=(t2 - th["t1",i-1] + 1/th["beta",i-1])) # DEBUG
      
    }
  }
  return(th)
}

theta <- MH.block.1(n=1000,t1.0 = 1900,l0.0 = 2,l1.0 = 1,beta.0 = 1, obs = coal, sigma = 10, beta.int=100)
theta.df = as.data.frame(t(theta))

relevant.plots(theta.df,coal,50)
