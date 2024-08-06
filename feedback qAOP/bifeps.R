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
bifeps = bifeps[,1:3]
bifeps = na.omit(bifeps)
bifeps = unique(bifeps)

#RENAME COLUMNS
names(bifeps)[1] <- 'MIE'
names(bifeps)[2] <- 'KE1'
names(bifeps)[3] <- 'KE2'

plot(bifeps$MIE, bifeps$KE2, xlab = "MIE", ylab="KE2")


bifpoints <- determine_BP(bifeps, xvar="MIE", yvar="KE2", res = "x")

bifeps_low <- bifeps %>% filter(KE2<=0.345560)
bifeps_up <- bifeps %>% filter(KE2>3.1)
bifeps_mid <- setdiff(setdiff(bifeps, bifeps_low), bifeps_up)
plot(1, xlab="MIE", ylab="KE2", type="n", xlim=c(0,0.30), ylim=c(0, 4.1))
lines(bifeps_low$MIE, bifeps_low$KE2,col="green")
lines(bifeps_up$MIE, bifeps_up$KE2, col="red")
lines(bifeps_mid$MIE, bifeps_mid$KE2, col="black")
abline( v=0.23125, lty=2)
abline( v=0.083, lty=2)
text(0.05, 0.4, 'stable')
text(0.15,4, 'stable')
text(0.15,1.2,'unstable')
