library(boot)
library(ggplot2)

#install.packages("invgamma")
library(invgamma)
  
data(coal)

coal.df <- coal #du kan nok bare hente ting rett ut fra coal 

Alg.4 <- function(n,t1.0,l0.0,l1.0,beta.0, obs, sigma){
  # initializing parameters with prior estimates
  
  theta <- matrix(0L, nrow=6, ncol=n)
  rownames(theta) <- c("t1", "l0", "l1", "beta", "y0", "y1")
  
  t0 <- coal$date[1]
  t2 <- tail(coal$date,1)
  
  theta["t1",1] <- t1.0
  theta["y0",1] <- max(which(coal$date < t1)) - 1 # disasters before t1, subtracting the first element which is not a disaster
  theta["y1",1] <- length(coal$date) - y0 - 2 # disasters at and after t1, subtracting start time and end time element
  
  theta["l0",1] <- l0.0
  theta["l1",1] <- l1.0
  theta["beta",1] <- beta.0
  
  for(i in 2:n){
    # sample beta from invGamma
    theta["beta",i] <- rinvgamma(1,shape=4, scale=theta["l0",i-1]+ theta["l1",i-1] + 1)
    
    # sample l0 from gamma
    theta["l0",i] <- rgamma(1,shape=theta["y0",i-1] + 2, scale=1/(theta["t1",i-1] - t0 + 1/theta["beta",i]))
    
    # sample l1 from gamma
    theta["l1",i] <- rgamma(1,shape = theta["y1",i-1] + 2, scale=1/(t2 - theta["t1",i-1] + 1/theta["beta",i]))
    
    # MH-step; sample t1 (f.ex. from normal with mean t1)
  
    y <- rnorm(1,t1,sigma)                # proposal for t1
    u <- runif(1)        
    lg.alpha <- min(0, (theta["l1",i]-theta["l0",i])*(y - theta["t1",i-1]))  # log acceptance probability
    alpha = exp(lg.alpha)                 # acceptance probability
    if (u < alpha){
      theta["t1",i] <- y
    }
    else{
      theta["t1",i] <- theta["t1",i-1]
    }
    # update y0 based on t1 and obs
    theta["y0",i] <- max(which(coal$date < theta["t1",i])) - 1
    # update y1 based on t1 and obs
    theta["y1",i] <- length(coal$date) - theta["y0",i] - 2
  }
  return(theta)
}

theta.test <- Alg.4(50,1870,5,1,3,coal,15)




















