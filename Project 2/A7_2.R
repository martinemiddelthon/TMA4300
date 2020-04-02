# Metropolis-Hastings block algorithm 1: 

data(coal)

sample.beta <- function(beta.0,sigma,l0,l1){
  beta <- rnorm(1,mean=beta.0,sigma)
  if(beta <= 0){
    return(beta.0)
  }
  else{
    # check if we accept new value
    u <- runif(1)
    prop.rat = beta.0^(5)/beta^(5)*exp((l0+l1+1)*(1/beta.0-1/beta))
    alpha = min(1,prop.rat)
    if(u < alpha){
      return(beta)
    }
    else{
      return(beta.0)
    }
  }
}

MH.block.2 <- function(n,t1.0,l0.0,l1.0,beta.0, obs, sigma1, sigma2, int){
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
    if(i %% int == 0){
      y <- rnorm(1,th["t1",i-1],sigma1)                # proposal for t1
      if (y > t2 | y < t0){
        # do not accept step, keep old state
        th["t1",i] <- th["t1",i-1] 
        th["beta",i] <- th["beta",i-1]
        th["l0",i] <- th["l0",i-1]
        th["l1",i] <- th["l1",i-1]
        th["y0",i] <- th["y0",i-1]
        th["y1",i] <- th["y1",i-1]
      }
      else {
        y0_prop <- max(which(coal$date <= y)) - 1
        y1_prop <- length(coal$date) - y0_prop - 2
        u <- runif(1)                  
        # acceptance probability:
        alpha <- min(1, th["l0",i-1]^(y0_prop - th["y0",i-1])
                     *th["l1",i-1]^(y1_prop  - th["y1",i-1])
                     *exp((th["l1",i-1]-th["l0",i-1])*(y - th["t1",i-1])))
        
        if (u < alpha){
          th["t1",i] <- y
        }
        else{
          th["t1",i] <- th["t1",i-1]
        }
        
        # update y0 and y1 based on current t1 and obs
        th["y0",i] <- max(which(coal$date <= th["t1",i])) - 1
        th["y1",i] <- length(coal$date) - th["y0",i] - 2
        
        #sample beta from normal distribution
        th["beta",i] <- sample.beta(th["beta",i-1],sigma2,l0=th["l0",i-1],l1=th["l1",i-1])
        
        # sample l0 from gamma
        th["l0",i] <- rgamma(1,shape=th["y0",i] + 2, scale=1/(th["t1",i] - t0 + 1/th["beta",i]))

        # sample l1 from gamma
        th["l1",i] <- rgamma(1,shape = th["y1",i] + 2, scale=1/(t2 - th["t1",i] + 1/th["beta",i]))

      }
    }
    else{
      # keep t1 constant 
      th["t1",i] <- th["t1",i-1]
      
      th["y0",i] <- max(which(coal$date <= th["t1",i])) - 1
      th["y1",i] <- length(coal$date) - th["y0",i] - 2
      
      #sample beta from normal distribution
      th["beta",i] <- sample.beta(th["beta",i-1],sigma2,l0=th["l0",i-1],l1=th["l1",i-1])
      
      # sample l0 from gamma
      th["l0",i] <- rgamma(1,shape=th["y0",i] + 2, scale=1/(th["t1",i] - t0 + 1/th["beta",i]))
      
      # sample l1 from gamma
      th["l1",i] <- rgamma(1,shape = th["y1",i] + 2, scale=1/(t2 - th["t1",i] + 1/th["beta",i]))
    }
  }
  return(th)
}

theta <- MH.block.2(n=10000,t1.0 = 1900,l0.0 = 2,l1.0 = 1,beta.0 = 1, obs = coal, sigma1 = 10, sigma2 = 2, int=200)
theta.df = as.data.frame(t(theta))

sample.beta(0,1,2,1)

relevant.plots(theta.df,coal,1000)
max(theta.df$beta)
mean(theta.df$beta)
min(theta.df$beta)

