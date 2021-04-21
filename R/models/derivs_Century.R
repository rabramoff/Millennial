derivs_Century <- function(step.num,state,parameters) {
  with(as.list(c(state,parameters)), {

  t_scalar <- (parameters$t2 + (parameters$t3 / pi) * atan(pi * parameters$t4 * (forc_st(step.num) - parameters$t1))) /
    (parameters$t2 + (parameters$t3 / pi) * atan(pi * parameters$t4 *(30.0 - parameters$t1)))

  w_scalar <- 1.0 / (1.0 + parameters$w1 * exp(-parameters$w2 * forc_sw(step.num)/0.39))

  f_TEX = parameters$c1 - parameters$c2*parameters$clay_silt*0.01 #converts from % to fraction
  
  f_ACTIVE <- ACTIVE * parameters$k_active * t_scalar * w_scalar * f_TEX
  
  f_SLOW <- SLOW * parameters$k_slow * t_scalar * w_scalar
  
  f_PASSIVE <- PASSIVE * parameters$k_passive * t_scalar * w_scalar
  
  dACTIVE <- parameters$forc_npp*(1-parameters$litter_to_slow) + f_SLOW * parameters$slow_to_active + f_PASSIVE * parameters$passive_to_active - f_ACTIVE
  
  dSLOW <- parameters$forc_npp*parameters$litter_to_slow + f_ACTIVE * (1-f_TEX-parameters$active_to_passive) - f_SLOW
  
  dPASSIVE <- f_ACTIVE * parameters$active_to_passive + f_SLOW * parameters$slow_to_passive - f_PASSIVE
            
    return(list(c(dACTIVE, dSLOW, dPASSIVE)))
  })
}
