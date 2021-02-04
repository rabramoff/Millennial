derivs_V2_MM <- function(step.num,state,parameters) {
  with(as.list(c(state,parameters)), {
          
# Soil type properties  
  #Equation 11
  #Mayes 2012, SSAJ
  kaff_lm = exp(-parameters$param_p1 * 7 - parameters$param_p2) * parameters$kaff_des

  #Equation 12
  #Georgiou in review
  param_qmax = parameters$param_bulkd * parameters$param_c1 * parameters$param_claysilt 
              
# Hydrological properties

 #Equation 5
 scalar_wd = (forc_sw(step.num) / parameters$porosity)^0.5

 #Equation 4
 scalar_wb = exp(parameters$lambda * -parameters$matpot) * (parameters$kamin + (1 - parameters$kamin) * ((parameters$porosity - forc_sw(step.num)) / parameters$porosity)^0.5) * scalar_wd

# Decomposition

  gas_const <- 8.31446
   
  #Equation 3
    vmax_pl = parameters$alpha_pl * exp(-parameters$eact_pl / (gas_const * (forc_st(step.num) + 273.15)))
    
  #Equation 2
    # POM -> LMWC
    if(POM>0 && MIC>0){
      f_PO_LM = vmax_pl * scalar_wd * POM * MIC / (parameters$kaff_pl + MIC)
    }else{
      f_PO_LM=0
    }
    
  #Equation 6
    # POM -> AGG
    if(POM>0){
      f_PO_AG = parameters$rate_pa * scalar_wd * POM
    }else{
      f_PO_AG=0
    }
    
  #Equation 7
    # AGG -> MAOM
    if(AGG>0){
      f_AG_break = parameters$rate_break * scalar_wd * AGG
    }else{
      f_AG_break=0
    }
  
  #Equation 9
    # LMWC -> out of system leaching
    if(LMWC>0){
      f_LM_leach = parameters$rate_leach * scalar_wd * LMWC
    }else{
      f_LM_leach=0
    }
  
  #Equation 10
    # LMWC -> MAOM
    if(LMWC>0 && MAOM>0){
      f_LM_MA = scalar_wd * kaff_lm * LMWC * (1 - MAOM / param_qmax)
    }else{
      f_LM_MA=0
    }
  
  #Equation 13
    # MAOM -> LMWC
    if(MAOM>0){
      f_MA_LM = parameters$kaff_des * MAOM / param_qmax
    }else{
      f_MA_LM=0
    }
    
  #Equation 15
    vmax_lb = parameters$alpha_lb * exp(-parameters$eact_lb / (gas_const * (forc_st(step.num) + 273.15)))
  
  #Equation 14
    # LMWC -> MIC
    if(LMWC>0 && MIC>0){
      f_LM_MB = vmax_lb * scalar_wb * MIC * LMWC / (parameters$kaff_lb + LMWC)
    }else{
      f_LM_MB=0
    }
  
  #Equation 16
    # MIC -> MAOM/LMWC
    if(MIC>0){
      f_MB_turn = parameters$rate_bd * MIC^2.0
    }else{
      f_MB_turn=0
    }
    
  #Equation 18
    # MAOM -> AGG
    if(MAOM>0){  
      f_MA_AG = parameters$rate_ma * scalar_wd * MAOM
    }else{
      f_MA_AG=0
    }
  
    #Equation 22
    # microbial growth flux, but is not used in mass balance

  #Equation 21
    # MIC -> atmosphere
    if(MIC>0 && LMWC>0){ 
    f_MB_atm = f_LM_MB * (1 - (parameters$cue_ref - parameters$cue_t * (forc_st(step.num) - parameters$tae_ref) ) )
    }else{
      f_MB_atm=0
    }
    
# Update state variables
    
  #Equation 1
    dPOM = forc_npp(step.num) * parameters$param_pi + f_AG_break * parameters$param_pa - f_PO_AG - f_PO_LM
  
  #Equation 8
    dLMWC = forc_npp(step.num) * (1. - parameters$param_pi) - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + f_MB_turn * (1. - parameters$param_pb) + f_MA_LM
  
  #Equation 17
    dAGG = f_MA_AG + f_PO_AG - f_AG_break
  
  #Equation 19
    dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * parameters$param_pb - f_MA_AG + f_AG_break * (1. - parameters$param_pa)
    
  #Equation 20
    dMIC = f_LM_MB - f_MB_turn - f_MB_atm

      return(list(c(dPOM, dLMWC, dAGG, dMIC, dMAOM, f_MB_atm)))
  })
}