#Title: derivs_V2_Century.R
#Author: Rose Abramoff, after Parton et al. (1987)
#Date: Sep 11, 2021

#This function contains the system of equations for the Century model, which was developed in Parton et al. (1987). 
#The equation numbers correspond to those in Abramoff et al. (2021) Appendix B, with parameters defined in Table A2.

##This function takes as arguments:
##1 number of time steps
##2 initial states of pools
##3 parameters
##4 forc_st: soil temperature in degrees Celcius
##5 forc_sw: volumetric soil moisture in mm3/mm3
##6 forc_npp: daily plant inputs in gC/m2/day
derivs_Century <- function(step.num,state,parameters,forc_st,forc_sw,forc_npp) {
  with(as.list(c(state,parameters)), {

#Equation B1
  t_scalar <- (parameters$t2 + (parameters$t3 / pi) * atan(pi * parameters$t4 * (forc_st(step.num) - parameters$t1))) /
    (parameters$t2 + (parameters$t3 / pi) * atan(pi * parameters$t4 *(30.0 - parameters$t1)))

#Equation B2
  w_scalar <- 1.0 / (1.0 + parameters$w1 * exp(-parameters$w2 * forc_sw(step.num)/0.39))

#Equation B3
  f_TEX = parameters$c1 - parameters$c2*0.8
  
#Equation B4
  f_StrLitter = StrLitter * parameters$k_strlitter * t_scalar * w_scalar * exp(-3*parameters$LigFrac)
  
#Equation B5
  f_MetLitter = MetLitter * parameters$k_metlitter * t_scalar * w_scalar  
  
#Equation B6
  f_ACTIVE <- ACTIVE * parameters$k_active * t_scalar * w_scalar * f_TEX

#Equation B7 
  f_SLOW <- SLOW * parameters$k_slow * t_scalar * w_scalar

#Equation B8
  f_PASSIVE <- PASSIVE * parameters$k_passive * t_scalar * w_scalar
  
#Equation B9
  dStrLitter = parameters$input_to_strlitter * forc_npp(step.num) - f_StrLitter

#Equation B10
  dMetLitter = (1-parameters$input_to_strlitter) * forc_npp(step.num) - f_MetLitter
  
#Equation B11
  dACTIVE <- (1-parameters$LigFrac) * parameters$strlitter_to_active * f_StrLitter + parameters$metlitter_to_active * f_MetLitter  + f_SLOW * parameters$slow_to_active + f_PASSIVE * parameters$passive_to_active - f_ACTIVE

#Equation B12
  dSLOW <-  parameters$LigFrac * parameters$strlitter_to_slow * f_StrLitter + f_ACTIVE * (1-f_TEX-parameters$active_to_passive) - f_SLOW
  
#Equation B13
  dPASSIVE <- f_ACTIVE * parameters$active_to_passive + f_SLOW * parameters$slow_to_passive - f_PASSIVE
            
    return(list(c(dStrLitter, dMetLitter, dACTIVE, dSLOW, dPASSIVE)))
  })
}
