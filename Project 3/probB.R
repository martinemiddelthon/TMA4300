

bilirubin <- read.table("https://www.math.ntnu.no/emner/TMA4300/2020v/Exercise/ex3-additionalFiles/bilirubin.txt", header=T)


#  1:

#  boxplot of logarithms of concentration for each indivicual:
boxplot(log(meas)~pers, data=bilirubin,
        main="Boxplot of Logarithm of Concentration", xlab="Individual", ylab="log(meas)")

#  fitted regression model:
reg <- lm(log(meas)~pers, data=bilirubin)
summary.lm <- summary(reg); summary.lm
Fval <- summary.lm$fstatistic; Fval

#  compare with 0.95 quantile of F distribution:
F.q <- qf(.95, df1=2, df2=26); F.q  # reject null hypothesis

#  defining function PermTest():

PermTest <- function(bilirubin){
  perm.df <- data.frame(meas=sample(bilirubin$meas), pers = bilirubin$pers)
  reg <- lm(log(meas)~pers, data=perm.df)
  Fval <- summary(reg)$fstatistic["value"]
  return(Fval)
}

#  generate a sample of 999 values of the F-statistic

F.samp <- rep(0,999)
for(i in 1:999){
  F.samp[i] <- PermTest(bilirubin)
}

#  compute p-value for Fval using F.samp:
p.val <- length(F.samp[F.samp>=Fval["value"]])/999; p.val
