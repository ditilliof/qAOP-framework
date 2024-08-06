library("deSolve")
library("phaseR")
AOP = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dKE1 = -d1 * KE1 + s1*MIE*KE2 +eps
    dKE2 = s2*(KE1**5)/(1+k2*(KE1**5)) -d2*KE2 + eps
    list(c(dKE1, dKE2))
  })
}

# set parameter values to be used in numerical simulation
pars = c(d1 = 1, s1 = 4, MIE=0.18, eps = 0.2, s2 = 4, k2 = 1, d2=1)
# set minimum and maximum values of the axes in the plots
xlim = c(0, 50)
ylim = c(0, 0.36)
# set initial state of variables
inival = c(0.2, 0.20127959)
# set time frame of simulation
finish = c(0, 5)
# set axis labels
axlabs = c("KE1", "KE2")
# perform numerical simulation and plot end result
nfkb.numericalSolution = numericalSolution(AOP, y0 = inival, tlim = finish, type = "one", parameters = pars, col = c("green", "orange"), ylim = ylim, state.names = c("KE1","KE2"), add.legend = F)
# make legend for plot
legend("topright", col = c("green", "orange"), legend = axlabs, lty = 1)

# Draw flowfield
AOP.flowfield <- flowField(xlab = "KE1", ylab = "KE2", xlim = c(0, 5), ylim = c(0, 5), AOP, parameters = pars,
                                  add = FALSE, points = 19, state.names = c("KE1", "KE2"))
# Draw nullclines
AOP.nullclines <- nullclines(AOP ,xlim = c(-10,10), ylim = c(-10,10),
                                    state.names = c("KE1", "KE2"), parameters = pars)

AOP.trajectory = trajectory(AOP, y0 = c(0.2, 0.2012), tlim = c(0, 10), parameters = pars,
                             state.names = c("KE1", "KE2"))
AOP.trajectory2 = trajectory(AOP, y0 = c(1, 0.5), tlim = c(0, 10), parameters = pars,
                            state.names = c("KE1", "KE2"))
AOP.trajectory3 = trajectory(AOP, y0 = c(0.5, 1), tlim = c(0, 10), parameters = pars,
                             state.names = c("KE1", "KE2"))

eq = c(tail(AOP.trajectory$x, 1), tail(AOP.trajectory$y, 1))
AOP.stability = stability(AOP, ystar = eq, parameters = pars, state.names = c("KE1","KE2"))
???
