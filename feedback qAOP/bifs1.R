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
bifs1 = bifs1[,1:2]
bifs1 = na.omit(bifs1)
bifs1 = unique(bifs1)

#RENAME COLUMNS
names(bifs1)[1] <- 's1'
names(bifs1)[2] <- 'KE2'

plot(bifs1$s1, bifs1$KE2, xlab = "s1", ylab="KE2")


bifpoints <- determine_BP(bifs1, xvar="s1", yvar="KE2", res = "x")

bifs1_low <- bifs1 %>% filter(KE2<=0.335725)
bifs1_up <- bifs1 %>% filter(KE2>3)
bifs1_mid <- setdiff(setdiff(bifs1, bifs1_low), bifs1_up)
plot(1, xlab="s1", ylab="KE2", type="n", xlim=c(0,6.2), ylim=c(0, 4.5))
lines(bifs1_low$s1, bifs1_low$KE2,col="green")
lines(bifs1_up$s1, bifs1_up$KE2, col="red")
lines(bifs1_mid$s1, bifs1_mid$KE2, col="black")
abline( v=5.16, lty=2)
abline( v=1.85, lty=2)
text(4, 4.37, 'stable')
text(1,0.4, 'stable')
text(4,0.85,'unstable')
