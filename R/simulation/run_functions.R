#Title: run_functions.R
#Author: Rose Abramoff
#Date: Sep 11, 2021

##This function takes as arguments:
##1 input data
##2 model equations
##3 parameters
##4 length of model run in years (optional; default = 100)
##5 initial states of pools (optional; default = 1 for all pools and 0 for CO2 flux)
Run_Model  <- function(inputdata,derivs,parameters,num.years=100,state=c(POM = 1, LMWC = 1, AGG = 1, MIC = 1, MAOM=1,CO2=0)) {

    #Define runtime
    run.steps <- seq(1,num.years*365)
    run.steps.minus.one <- run.steps[1:length(run.steps)-1]

    #Define forcing functions
    forc_st <- approxfun(run.steps, rep(inputdata$forc_st,num.years)) #temperature input function
    forc_sw <- approxfun(run.steps, rep(inputdata$forc_sw,num.years)) #moisture input function
    forc_npp <- approxfun(run.steps,rep(inputdata$forc_npp,num.years)) #NPP input function

    #Run model
    output <- ode(y = state, times=run.steps, func= derivs, parms = parameters, forc_st=forc_st, forc_sw=forc_sw, forc_npp=forc_npp, method="rk4") #solve ode, return output
    return(output)
}

##This function takes as arguments:
##1 input data
##2 model equations
##3 parameters
##4 whether or not to calculate the eigenvalues (optional; default = 0; 0 = no, 1 = yes)
##5 initial states of pools (optional; default = 1 for all pools and 0 for CO2 flux)
Solve_Model  <- function(inputdata,derivs,parameters,calc_eigens=0,state=c(POM = 1, LMWC = 1, AGG = 1, MIC = 1, MAOM=1)) {

    SStime = 1000*365 #time at which steady state solution is evaluated

    #Define constant forcing functions
    forc_st <- approxfun(1:SStime, rep(mean(inputdata$forc_st),SStime)) #temperature input function
    forc_sw <- approxfun(1:SStime, rep(mean(inputdata$forc_sw),SStime)) #moisture input function
    forc_npp <- approxfun(1:SStime, rep(mean(inputdata$forc_npp),SStime)) #NPP input function

    derivs_SS_wrapper <- function(step.num,state,parameters,forc_st,forc_sw,forc_npp) {
      output <- derivs(step.num,state,parameters,forc_st,forc_sw,forc_npp)
      return(list(output[[1]][1:5]))
    }
    
    if(calc_eigens==1){
      #Calculate and save eigenvalues
      eigens <- eigen(jacobian.full(y = state, func = derivs_SS_wrapper, parms = parameters, forc_st=forc_st, forc_sw=forc_sw, forc_npp=forc_npp, time = 1))
      write.csv(eigens$values, file= paste0("eigens_",as.character(substitute(derivs)),".csv"))
    }
    
    #Solve steady state
    SS_output <- stode(y = state, time = SStime, func = derivs_SS_wrapper, parms = parameters, forc_st=forc_st, forc_sw=forc_sw, forc_npp=forc_npp, positive = TRUE)
    return(SS_output)
}
