derivs_Century <- function(step.num,state,parameters) {
  with(as.list(c(state,parameters)), {

  t_scalar <- (parameters$t2 + (parameters$t3 / pi) * atan(pi * parameters$t4 * (forc_st(step.num) - parameters$t1))) /
    (parameters$t2 + (parameters$t3 / pi) * atan(pi * parameters$t4 *(30.0 - parameters$t1)))

  w_scalar <- 1.0 / (1.0 + parameters$w1 * exp(-parameters$w2 * forc_sw(step.num)))

  f_DOC_leaching <- DOC * parameters$k_leaching * (DOC / (DOC + 100)) * t_scalar * w_scalar

  f_DOC_ATM <- DOC * parameters$k_doc * t_scalar * w_scalar

  f_ACTIVE_ATM <- ACTIVE * parameters$k_active * t_scalar * w_scalar

  f_SLOW_ATM <- SLOW * parameters$k_slow * t_scalar * w_scalar

  f_PASSIVE_ATM <- PASSIVE * parameters$k_passive * t_scalar * w_scalar

  f_ACTIVE_DOC <- ACTIVE * parameters$active_to_doc * t_scalar * w_scalar

  f_ACTIVE_SLOW <- ACTIVE * parameters$active_to_slow * t_scalar * w_scalar

  f_SLOW_PASSIVE <- SLOW * parameters$slow_to_passive * t_scalar

  f_ACTIVE_PASSIVE <- ACTIVE * t_scalar * parameters$active_to_passive

  f_PASSIVE_ACTIVE <- parameters$passive_to_active * t_scalar * w_scalar

  dDOC <- f_ACTIVE_DOC - f_DOC_ATM - f_DOC_leaching
           
  dACTIVE <- forc_npp(step.num) + f_PASSIVE_ACTIVE - f_ACTIVE_DOC - f_ACTIVE_SLOW - f_ACTIVE_PASSIVE - f_ACTIVE_ATM

  dSLOW <- f_ACTIVE_SLOW - f_SLOW_PASSIVE - f_SLOW_ATM

  dPASSIVE <- f_SLOW_PASSIVE + f_ACTIVE_PASSIVE - f_PASSIVE_ACTIVE - f_PASSIVE_ATM
            
    return(list(c(dDOC, dACTIVE, dSLOW, dPASSIVE)))
  })
}
