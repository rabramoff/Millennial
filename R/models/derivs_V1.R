#Title: derivs_V1.R
#Author: Rose Abramoff translated Xiaofeng Xu's original fortran code into R. 
#Ben Bond-Lamberty also has a prototype R translation that is probably more organized than this one: https://github.com/PNNL-TES/millenial/tree/master/R
#Date: Sep 11, 2021

derivs_V1 <- function(step.num,state,parameters,forc_st,forc_sw,forc_npp) {
  with(as.list(c(state,parameters)), {

  # Equation 3 in Abramoff et al. (2017)
  t_scalar <- (parameters$t2 + (parameters$t3 / pi) * atan(pi * parameters$t4 * (forc_st(step.num) - parameters$t1))) /
    (parameters$t2 + (parameters$t3 / pi) * atan(pi * parameters$t4 *(30.0 - parameters$t1)))

  w_scalar <- 1.0 / (1.0 + parameters$w1 * exp(-parameters$w2 * forc_sw(step.num)/0.35))

  # Equation 8 in Abramoff et al. (2017) except Xiaofeng removed water impact after review at GBC, June 2017
  f_LM_leaching <- LMWC * parameters$k_leaching * t_scalar #* w_scalar ! Xiaofeng removed water impact, after review at GBC June,2017         

  # if (MAOM > parameters$M_Lmin) {
  #   f_MI_LM_des <- parameters$Vm_l * (MAOM - parameters$M_Lmin) / (parameters$km_l + MAOM - parameters$M_Lmin) * t_scalar * w_scalar
  # } else {
  #   f_MI_LM_des <- 0
  # }
#not in equations .. removed, but in main.F90 is still in ODEs
  
  # LMWC -> MINERAL: This desorption function is from Mayes 2012, SSAJ
  klmc_min <- (10.0 ^ (-0.186 * parameters$param_pH - 0.216)) / 24.0

  # Equation 11 in Abramoff et al. (2017)
  # "Qmax is the maximum sorption capacity (mg C kg-1 dry soil) and is converted to
  # C density (g C m-2) by multiplying soil bulk density (BD = 1350 kg m-3), assuming
  # a 1 m soil profile. The parameters c1 and c2 are the coefficients for computing
  # Qmax from the clay content in percent, derived from Mayes et al. (2012)."
  # 2019-03-09 note: per Xiaofeng Xu, they ended up using 1.0 for bulk density
  Qmax <- 10.0 ^ (parameters$c1 * log(parameters$param_clay * 100.0) + parameters$c2 - 0.50)

  f_LM_MI_sor <- (((klmc_min * Qmax * LMWC ) / (2 + klmc_min * LMWC) - MAOM) / Qmax + 0.0015) * LMWC / 50 * t_scalar * w_scalar

  if (f_LM_MI_sor < (LMWC * 0.9)) {
    f_LM_MI_sor <- f_LM_MI_sor
  } else {
    f_LM_MI_sor <- LMWC * 0.9
  }

  f_LM_MB_uptake <- LMWC * parameters$klmc * t_scalar * w_scalar * MIC / (MIC + parameters$kes) * LMWC / (20 + LMWC)
  temp2 = f_LM_MB_uptake * (1 - (parameters$CUEref - parameters$CUET * (forc_st(step.num) - parameters$Taeref)))
temp2 <- max(temp2, 0)
f_LM_MB_uptake = f_LM_MB_uptake - temp2

  # POM -> LMWC
  # "Decomposition of POM is governed by a double Michaelis–Menten equation,
  # where Vpl is the maximum rate of POM decomposition, Kpl is the half-
  # saturation constant, B is the microbial biomass carbon, and Kpe is the
  # half-saturation constant of microbial control on POM mineralization.
  # Equation 2 in Abramoff et al. (2017)
  f_PO_LM_dep <- parameters$Vpom_lmc * POM / (POM + parameters$kpom) * t_scalar * w_scalar #* (1. - MB / (MB + k_POMes))

  f_PO_LM_dep <- min(f_PO_LM_dep,  0.9 * POM)
  
  if (MIC > 0 & MAOM < Qmax) {
    f_MB_MI_sor <- MIC * parameters$kmic * 0.15 * t_scalar * w_scalar  #* (MB / 200) * (MB / 200)
  } else {
    f_MB_MI_sor <- 0
  }

  f_MB_MI_sor <- min(f_MB_MI_sor, 0.9 * MIC)
  
  f_MB_MI_sor <- max(f_MB_MI_sor, 0)

  f_MB_atm <- temp2 + MIC * parameters$kmic * t_scalar * w_scalar

  # POM -> SOILAGG
  # "The formation of aggregate C (A) from POM follows Michaelis–Menten dynamics,
  # where Vpa is the maximum rate of aggregate formation, Kpa is the half-
  # saturation constant of aggregate formation, and Amax is the maximum capacity
  # of C in soil aggregates.
  # Eauation 5 in Abramoff et al. (2017)
  f_PO_SO_agg <- parameters$Vpom_agg * POM / (parameters$kpom_agg + POM) * (1 - AGG / parameters$AGGmax) * t_scalar * w_scalar

  f_PO_SO_agg <- min(f_PO_SO_agg, 0.9 * POM)

  f_MI_SO_agg <- parameters$Vmin_agg * MAOM / (parameters$kmin_agg + MAOM) * (1 - AGG / parameters$AGGmax) #* t_scalar * w_scalar
 
  f_MI_SO_agg <- min(f_MI_SO_agg, 0.9 * MAOM)

  # Soil aggregate C breakdown is partitioned to POM and MAOM,
  # where kb is the rate of breakdown.
  # Eauation 6 in Abramoff et al. (2017)
  f_SO_break <- AGG * parameters$kagg * t_scalar * w_scalar
  f_SO_PO_break <- f_SO_break * 1.5 / 3
  f_SO_MI_break <- f_SO_break * 1.5 / 3

  if((f_PO_LM_dep + f_PO_SO_agg) > POM) {
    temp3 <- POM / (f_PO_LM_dep + f_PO_SO_agg)
    f_PO_LM_dep <- f_PO_LM_dep * temp3
    f_PO_SO_agg <- f_PO_SO_agg * temp3
  }
  
# Update state variables
  # Equation 7 in Abramoff et al. (2017)
  dLMWC <- (f_PO_LM_dep - f_LM_leaching - f_LM_MI_sor - f_LM_MB_uptake - temp2) + forc_npp(step.num) / 3

  # Equation 1 in Abramoff et al. (2017)
  dPOM <- (f_SO_PO_break - f_PO_LM_dep - f_PO_SO_agg) + forc_npp(step.num) * 2 / 3

  dMIC <- (f_LM_MB_uptake - f_MB_MI_sor - f_MB_atm)

  dMAOM <- (f_LM_MI_sor + f_MB_MI_sor + f_SO_MI_break - f_MI_SO_agg)

  dAGG <- (f_PO_SO_agg + f_MI_SO_agg - f_SO_PO_break - f_SO_MI_break)

    return(list(c(dPOM, dLMWC, dAGG, dMIC, dMAOM, f_MB_atm)))
  })
}