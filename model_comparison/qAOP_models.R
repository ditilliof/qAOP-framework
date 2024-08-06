library(deSolve)

# Define model functions
qAOPA = function(timepoint, state, parameters){
  with(as.list(c(state, parameters)), {
    dMIE = -tau * MIE
    dKE1 = eps + MIE - d1 * KE1
    dKE2 = k_12 * KE1 - d2 * KE2 + k_2 * MIE
    
    list(c(dMIE, dKE1, dKE2))
  })
}

qAOPB = function(timepoint, state, parameters){
  with(as.list(c(state, parameters)), {
    dMIE = -tau * MIE
    dKE1 = eps + MIE - d1 * KE1
    dKE2 = k_12 * KE1 - d2 * KE2
    
    list(c(dMIE, dKE1, dKE2))
  })
}

qAOPC = function(timepoint, state, parameters){
  with(as.list(c(state, parameters)), {
    dMIE = -tau * MIE
    dKE1 = eps + MIE - d1 * KE1
    dKE2 = k_12 * KE1 - d2 * KE2 + k_2 / (h_2 + MIE)
    
    list(c(dMIE, dKE1, dKE2))
  })
}
