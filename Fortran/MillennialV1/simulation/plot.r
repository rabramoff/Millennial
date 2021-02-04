setwd("/Users/xxuadmin/BUSINESS/PUBLICATIONS/WorkingOn_Abramoff_Perspective/millennial_code_2017July_Century_WT")

# output
data <- read.table("outputcontrol.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
outputcontrol <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
#par(mai = c(1,1,1,1,1))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")


pdf(file = "outputcontrol.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("Control run:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("output5cwarm.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output5cwarm <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output5cwarm.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("5 C warming:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("output10clay.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output10clay <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output10clay.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("10% Clay content:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("output20clay.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output20clay <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output20clay.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("20% Clay content:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("output30clay.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output30clay <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output30clay.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("30% clay content:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("output40clay.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output40clay <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output40clay.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("40% clay content:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("output50clay.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output50clay <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output50clay.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("50% clay content:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("output60clay.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output60clay <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output60clay.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("60% clay content:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end
# output
data <- read.table("output70clay.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output70clay <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output70clay.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("70% Clay content:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("output80clay.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
output80clay <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "output80clay.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("80% Clay content:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("outputdoubleCinput.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
outputdoubleCinput <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "outputdoubleCinput.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("Double C input:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("outputhalfwater.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
outputhalfwater <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "outputhalfwater.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("Half water moisture:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("outputTandCinput.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
outputTandCinput <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "outputTandCinput.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("5C Warming and Double C input:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end

# output
data <- read.table("outputTandW.txt")
LMC <- data$V2[729636:730000]
POM <- data$V3[729636:730000]
MB <- data$V4[729636:730000]
MINERAL <- data$V5[729636:730000]
AGG <- data$V6[729636:730000]
outputTandW <- data$V2 + data$V3 + data$V4 + data$V5 + data$V6

summary(LMC)
summary(POM)
summary(MB)
summary(MINERAL)
summary(AGG)
lmcaverage = mean(LMC)
pomaverage = mean(POM)
mbaverage = mean(MB)
mineralaverage = mean(MINERAL)
aggaverage = mean(AGG)

carbon <- c(lmcaverage, pomaverage, mbaverage, mineralaverage, aggaverage)
#lbls <- c("LMWC", "POM", "MB", "MINERAL", "AGGREGATE")
#pct <- round(carbon/sum(carbon)*100)
#lbls <- paste(lbls,pct) # add percents to labels
#lbls <- paste(lbls, "%", sep="")
#pie(carbon, labels=lbls, col=rainbow(length(lbls)),main=round(sum(carbon)))
colors <- c("black","red","green","blue","pink")
numb_labels <- round(carbon/sum(carbon) * 100,1)
numb_labels <- paste(numb_labels, "%", sep="")
xx <- c("LMWC","POM","MB","MINERAL","AGGREGATE")
xx <- paste(round(carbon/sum(carbon) * 100,1),xx,sep=" ")

pdf(file = "outputTandW.pdf")
pie(carbon, labels="", col=colors,clockwise=TRUE,main=paste("5C Warming and Half Water:", round(sum(carbon)),"gC",sep=" "))
legend(-1.25,-0.5, legend = xx, fill=colors,bty="n")
dev.off()

rm(data)
rm(LMC)
rm(POM)
rm(MB)
rm(MINERAL)
rm(AGG)
# end


#plot time series of total carbon for all scenarios
#outputcontrol
#output5cwarm
#output10clay
#output20clay
#output30clay
#output40clay
#output50clay
#output60clay
#output70clay
#output80clay
#outputdoubleCinput
#outputhalfwater
#outputTandCinput
#outputTandW
pdf(file = "timeseries.pdf")
day=c(1:136500)
plot(c(1,175000),c(-2000,2000),type="n",xlab="Day",ylab="changed in total C storage")
lines(day,output5cwarm[1:136500] - outputcontrol[1:136500],col="blue")
lines(day,output10clay[1:136500] - outputcontrol[1:136500],col="brown")
lines(day,output20clay[1:136500] - outputcontrol[1:136500],col="yellow")
lines(day,output30clay[1:136500] - outputcontrol[1:136500],col="cyan")
lines(day,output40clay[1:136500] - outputcontrol[1:136500],col="darkblue")
lines(day,output50clay[1:136500] - outputcontrol[1:136500],col="darkgreen")
lines(day,output60clay[1:136500] - outputcontrol[1:136500],col="darkgray")
lines(day,output70clay[1:136500] - outputcontrol[1:136500],col="darkorange")
lines(day,output80clay[1:136500] - outputcontrol[1:136500],col="darkred")
lines(day,outputdoubleCinput[1:136500] - outputcontrol[1:136500],col="deeppink")
lines(day,outputhalfwater[1:136500] - outputcontrol[1:136500],col="deepskyblue")
lines(day,outputTandCinput[1:136500] - outputcontrol[1:136500],col="darkgreen")
lines(day,outputTandW[1:136500] - outputcontrol[1:136500],col="lightgreen")
legend(137500,1000,c("warming","10% clay","20% clay","30% clay","40% clay","50% clay","60% clay","70% clay","80% clay","double C input","half water","warming and double C input","warming and half water"),lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","brown","yellow","cyan","darkblue","darkgreen","darkgray","darkorange","darkred","deeppink","deepskyblue","darkgreen","lightgreen"))
dev.off()

