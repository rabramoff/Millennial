#f(T) and f(W) comparison

#century temperature function
soilTemp <- seq(-20,40,length.out = 100)
teff <- c(15.4, 11.75, 29.7, 0.031)
tfunc <- (teff[2] + (teff[3]/pi)* atan(pi*teff[4]*(soilTemp - teff[1]))) / (teff[2] + (teff[3]/pi)* atan(pi*teff[4]*(30 - teff[1])))

#century water function
relwc <- seq(0,1,length.out = 100)
wfunc <- 1/(1 + 30 * exp(-9*relwc))

#millennial temperature function
Q10 <- 2
Tref <- 15
mill.tfunc <- Q10^((soilTemp - Tref)/10)

#millennial water function
psimin <- -10
psimax <- -0.01
psiopt <- -3.33
psi <- seq(-10,-0.01,length.out = 100)
mill.wfunc <- ( ( (psi - psimin)*(psi- psimax) ) / ( (psi - psimin)*(psi- psimax) - (psi - psiopt)*(psi - psiopt) ) )^0.5

#plots
#all temps in degrees C
pdf(file="/Users/xxuadmin/BUSINESS/PUBLICATIONS/WorkingOn_Abramoff_Perspective/millennial_code_2017July_Century_WT/tandwfuncs.pdf", height = 4.52756*1, width = 4.52756*2)

par(mfrow = c(1,2))
plot(soilTemp, tfunc, ylim = c(0,6), xlab = expression("Soil Temperature (" ~ degree ~ "C)"), ylab = "f(T)", type = "l")
lines(soilTemp, mill.tfunc, lty = 2, col=2)
legend("topleft", c("Century","Millennial"), col=c(1,2), lty=c(1,2))

#relative water content or soil water potential
plot(relwc, wfunc, xlab = "Relative Water Content", ylab = "f(W)", type = "l") 
par(new=T)
plot(psi, mill.wfunc, axes=F, xlab="", ylab="", type= "l", col = 2, lty = 2)
axis(3, ylim=c(-10,0), col="red",col.axis="red")
mtext("Soil Water Potential (MPa)", side=3, col="red", line=2.5)

dev.off()
