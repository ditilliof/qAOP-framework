library(dplyr)
#DEFINE FUNCTION TO GET BIFURCARION POINTS
determine_BP <- function(df, xvar, yvar, res = "y") {
  sf <- splinefun(df[[yvar]], df[[xvar]], method = "natural", ties = mean)
  roots <- rootSolve::uniroot.all(sf, interval = c(min(df[[yvar]]), max(df[[yvar]])), deriv = 1)
  if(res == "y") return(roots)
  else if(res == "x") return(sf(roots))
  else stop("invalid option for res")
}
#CLEAN THE DATA
posbif2AOP = posbif2AOP[,1:3]
posbif2AOP = na.omit(posbif2AOP)
posbif2AOP = unique(posbif2AOP)

#RENAME COLUMNS
names(posbif2AOP)[1] <- 'MIE'
names(posbif2AOP)[2] <- 'KE1'
names(posbif2AOP)[3] <- 'KE2'

plot(posbif2AOP$MIE, posbif2AOP$KE2, xlab = "MIE", ylab="KE2")


bifpoints <- determine_BP(posbif2AOP, xvar="MIE", yvar="KE2", res = "x")

posbif2AOP_low <- posbif2AOP %>% filter(KE2<=0.345560)
posbif2AOP_up <- posbif2AOP %>% filter(KE2>3.1)
posbif2AOP_mid <- setdiff(setdiff(posbif2AOP, posbif2AOP_low), posbif2AOP_up)
plot(1, xlab="MIE", ylab="KE2", type="n", xlim=c(0,0.30), ylim=c(0, 4.1))
lines(posbif2AOP_low$MIE, posbif2AOP_low$KE2,col="green")
lines(posbif2AOP_up$MIE, posbif2AOP_up$KE2, col="red")
lines(posbif2AOP_mid$MIE, posbif2AOP_mid$KE2, col="black")
abline( v=0.23125, lty=2)
abline( v=0.083, lty=2)
text(0.05, 0.4, 'stable')
text(0.15,4, 'stable')
text(0.15,1.2,'unstable')
