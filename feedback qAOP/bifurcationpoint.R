x = posbif2AOP[,1] 
y = posbif2AOP[,2]
z= posbif2AOP[,3]
colnames(x) = "MIE"
colnames(y)= "KE1"
colnames(z)= "KE2"
names(posbif2AOP)[1] <- 'MIE'
names(posbif2AOP)[2] <- 'KE1'
names(posbif2AOP)[3] <- 'KE2'
plot(x$MIE, y$KE1, xlab = "MIE", ylab="KE1")
plot(x$MIE, z$KE2, xlab = "MIE", ylab="KE2")
abline(v= 0.08139822)
abline(v= 0.23428119)

determine_BP <- function(df, xvar, yvar, res = "y") {
  sf <- splinefun(df[[yvar]], df[[xvar]], method = "natural", ties = mean)
  roots <- rootSolve::uniroot.all(sf, interval = c(min(df[[yvar]]), max(df[[yvar]])), deriv = 1)
  if(res == "y") return(roots)
  else if(res == "x") return(sf(roots))
  else stop("invalid option for res")
}

det_BP <- function(df, xvar, yvar) {
  sf <- splinefun(df[[yvar]], df[[xvar]], method = "natural", ties = mean)
  roots <- rootSolve::uniroot.all(sf, interval = c(min(df[[yvar]]), max(df[[yvar]])), deriv = 1)
 return(roots)
}

posbif2AOP = posbif2AOP[,1:3]
posbif2AOP = na.omit(posbif2AOP)

determine_BP(df=posbif2AOP, xvar="MIE", yvar="KE2", res = "x" )

