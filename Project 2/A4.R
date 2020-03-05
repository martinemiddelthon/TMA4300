library(boot)
library(ggplot2)

#install.packages("invgamma")
library(invgamma)

#install.packages("tidyverse")
library(tidyverse)

#install.packages("ggpubr")
library(ggpubr)
  
data(coal)

set.seed(124)

Alg.4 <- function(n,t1.0,l0.0,l1.0,beta.0, obs, sigma){
  # initializing parameters with prior estimates
  
  theta <- matrix(0L, nrow=6, ncol=n)
  rownames(theta) <- c("t1", "l0", "l1", "beta", "y0", "y1")
  
  t0 <- coal$date[1]
  t2 <- tail(coal$date,1)
  
  theta["t1",1] <- t1.0
  theta["y0",1] <- max(which(coal$date < theta["t1",1])) - 1 # disasters before t1, subtracting the first element which is not a disaster
  theta["y1",1] <- length(coal$date) - theta["y0",1] - 2 # disasters at and after t1, subtracting start time and end time element
  
  theta["l0",1] <- l0.0
  theta["l1",1] <- l1.0
  theta["beta",1] <- beta.0
  
  for(i in 2:n){
    
    # MH-step; sample t1 (f.ex. from normal with mean t1)
  
    y <- rnorm(1,theta["t1",i-1],sigma)                # proposal for t1
    if (y > t2 | y < t0){
      #cat("y utenfor")
      theta["t1",i] <- theta["t1",i-1] # etterpå: prøv evt å ha dette først
      theta["beta",i] <- theta["beta",i-1]
      theta["l0",i] <- theta["l0",i-1]
      theta["l1",i] <- theta["l1",i-1]
    }
    else {
      y0_prop <- max(which(coal$date <= y)) - 1
      y1_prop <- length(coal$date) - y0_prop - 2
      u <- runif(1)        
      #lg.alpha <- min(0, (theta["l1",i-1]-theta["l0",i-1])*(y - theta["t1",i-1]))  # log acceptance probability
      #alpha = exp(lg.alpha)                 # acceptance probability
      alpha <- min(1, theta["l0",i-1]^(y0_prop - theta["y0",i-1])
                   *theta["l1",i-1]^(y1_prop  - theta["y1",i-1])
                   *exp((theta["l1",i-1]-theta["l0",i-1])*(y - theta["t1",i-1])))  # DEBUG
          
      if (u < alpha){
        theta["t1",i] <- y
      }
      else{
        theta["t1",i] <- theta["t1",i-1]
      }
      # sample beta from invGamma
      #theta["beta",i] <- rinvgamma(1,shape=4, scale=theta["l0",i-1]+ theta["l1",i-1] + 1)
      theta["beta",i] <- rinvgamma(1,shape=4, scale=1/(theta["l0",i-1]+ theta["l1",i-1] + 1))  # DEBUG
      
      # sample l0 from gamma
      theta["l0",i] <- rgamma(1,shape=theta["y0",i-1] + 2, scale=1/(theta["t1",i-1] - t0 + 1/theta["beta",i-1]))
      #theta["l0",i] <- rgamma(1,shape=theta["y0",i-1] + 2, scale=theta["t1",i-1] - t0 + 1/theta["beta",i-1]) #DEBUG
      
      # sample l1 from gamma
      theta["l1",i] <- rgamma(1,shape = theta["y1",i-1] + 2, scale=1/(t2 - theta["t1",i-1] + 1/theta["beta",i-1]))
      #theta["l1",i] <- rgamma(1,shape = theta["y1",i-1] + 2, scale=(t2 - theta["t1",i-1] + 1/theta["beta",i-1])) # DEBUG
      
    }
    # update y0 based on current t1 and obs
    theta["y0",i] <- max(which(coal$date <= theta["t1",i])) - 1
    # update y1 based on current t1 and obs
    theta["y1",i] <- length(coal$date) - theta["y0",i] - 2
  }
  return(theta)
}

relevant.plots <- function(theta.df,coal,burn.in){

  t1.mean <- mean(theta.df[burn.in:length(theta.df$t1),]$t1)
  cat("t1.mean", t1.mean)
  
  hist <- ggplot(theta.df,aes(x=t1)) + geom_histogram(binwidth = 1)
  hist <- hist + geom_vline(xintercept = t1.mean, col="red")
  #hist
  
  time.proc <- ggplot(theta.df,aes(x=seq(1,length(t1),1),y=t1)) + geom_line()
  #time.proc
  
  l0.proc <- ggplot(theta.df,aes(x=seq(1,length(t1),1),y=l0)) + geom_line()
  #l0.proc
  
  l1.proc <- ggplot(theta.df,aes(x=seq(1,length(t1),1),y=l1)) + geom_line()
  #l1.proc
  
  beta.proc <- ggplot(theta.df,aes(x=seq(1,length(t1),1),y=beta)) + geom_line()
  #beta.proc
  
  l0.mean <- mean(theta.df[burn.in:length(theta.df$t1),]$l0)
  cat("l0.mean",l0.mean)
  
  l1.mean <- mean(theta.df[burn.in:length(theta.df$t1),]$l1)
  cat("l1.mean", l1.mean)
  
  #values in helplines:
  x0 <- coal$date[1]
  yend0 <- l0.mean*(t1.mean - x0)
  xend1 <- tail(coal$date,1)
  yend1 <- l1.mean*(xend1 - t1.mean) + yend0
  
  
  compare <- ggplot() + geom_point(data=coal, aes(x=date,y=seq(1,length(date),1)))
  compare <- compare + geom_segment(aes(x=x0, xend=t1.mean, y = 0, yend = yend0))
  compare <- compare + geom_segment(aes(x = t1.mean, xend=xend1, y = yend0, yend = yend1))
  #compare
  
  ggarrange(ggarrange(hist,compare,ncol=2, labels=c("hist","compare")),
            ggarrange(l0.proc,l1.proc,ncol=2, labels=c("l0.porc", "l1.proc")),
            ggarrange(time.proc,beta.proc,ncol=2, labels=c("time.proc","beta.proc")),nrow=3)
}



test <- Alg.4(100000,1890,2.5,1,0.2,coal,10)
test.df <- as.data.frame(t(test))

theta <- Alg.4(10000,1930,5,1,3,coal,5)
theta.df <- as.data.frame(t(theta))
   
relevant.plots(theta.df,coal)
relevant.plots(test.df,coal, 7400)




t1.mean <- mean(test.df[7400:length(test.df$t1),]$t1)
cat("t1.mean", t1.mean)

hist <- ggplot(test.df,aes(x=t1)) + geom_histogram(binwidth = 1)
hist <- hist + geom_vline(xintercept = t1.mean, col="red")
hist

time.proc <- ggplot(test.df,aes(x=seq(1,length(t1),1),y=t1)) + geom_line()
time.proc

l0.proc <- ggplot(test.df,aes(x=seq(1,length(t1),1),y=l0)) + geom_line()
l0.proc

l1.proc <- ggplot(test.df,aes(x=seq(1,length(t1),1),y=l1)) + geom_line()
l1.proc

beta.proc <- ggplot(test.df,aes(x=seq(1,length(t1),1),y=beta)) + geom_line()
beta.proc

l0.mean <- mean(test.df[7400:length(test.df$t1),]$l0)
cat("l0.mean",l0.mean)

l1.mean <- mean(test.df[7400:length(test.df$t1),]$l1)
cat("l1.mean", l1.mean)

#values in helplines:
x0 <- coal$date[1]
yend0 <- l0.mean*(t1.mean - x0)
xend1 <- tail(coal$date,1)
yend1 <- l1.mean*(xend1 - t1.mean) + yend0


compare <- ggplot() + geom_point(data=coal, aes(x=date,y=seq(1,length(date),1)))
compare <- compare + geom_segment(aes(x=x0, xend=t1.mean, y = 0, yend = yend0))
compare <- compare + geom_segment(aes(x = t1.mean, xend=xend1, y = yend0, yend = yend1))
compare

ggarrange(ggarrange(hist,compare,ncol=2, labels=c("hist","compare")),
          ggarrange(l0.proc,l1.proc,ncol=2, labels=c("l0.porc", "l1.proc")),
          ggarrange(time.proc,beta.proc,ncol=2, labels=c("time.proc","beta.proc")),nrow=3)
                                  



















