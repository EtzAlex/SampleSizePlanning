set.seed(2)
options(scipen=20) #disable scientific notation for numbers

nSim<-10 #numbber of simulated studies

library(pwr)
library(MBESS)
library(gsDesign) # The group sequential design package
library(BayesFactor)

#Settings for sequential analysis
n<-200 #total number of datapoints you are willing to collect
Looks<-4 #set number of looks at the data
#Determine boundaries pocock function. 
seqdesign <- gsDesign(k=Looks, test.type=1, alpha=0.05, sfu=sfLDPocock, beta=.2)
plot(seqdesign)
#get alpha boundaries for each look
alphalook<-(pnorm(-abs(seqdesign$upper$bound)))
#get N at each look
LookN<-(ceiling(n*seqdesign$timing))

#define parameters for N*2.5 calculations
D<-0.5 #True effect size (Keep SD below to 1, otherwise, this is just mean dif, not d)
n_orig<-20 #sample size in original study
SD<-1 #Set True standard deviation.
alphalookn25<-0.05 #Type 1 error rate for single test after n*2.5

#Define scale for Bayes Factor. Adjust depending on expected effect size. Old default was 1 (put more weight on larger effect sizes), changed to 0.707 (Morey et al., 2011). Can be set to 0.5 when expecting small effect sizes (Rouder et al., 2009).  
rscaleBF<- 0.5

#define some variables needed for calculations
mat1<-matrix(NA, nrow=nSim, ncol=Looks) #Matrix for p-values at sequential tests
matdobs<-matrix(NA, nrow=nSim, ncol=Looks) #Matrix for d at sequential tests
matBF<-matrix(NA, nrow=nSim, ncol=Looks) #Matrix for BayesFactors
mat2<-matrix(NA, nrow=nSim, ncol=Looks)
TotalN<-numeric(Looks)
supportH0<-numeric(Looks) #Variable that stores how many studies support H0 (Bayes Factor < 0.3)
supportH1<-numeric(Looks) #Variable that stores how many studies support H0 (Bayes Factor > 3)
lookn25sig<-numeric(nSim)
p_n2.5<-numeric(nSim)
pmat<-matrix(NA, nrow=nSim, ncol=n) #matrix for all p-values (for plot)
dmat<-matrix(NA, nrow=nSim, ncol=n) #matrix for all effect sizes d (for plot) 
p<-numeric(n) 
t<-numeric(n) 
x<-numeric(n)
y<-numeric(n)
dlist<-numeric(n)

#Run simulation - create 1 participant at a time, store all data in a list, so I can do tests after any number of participants I want
for (j in 1:nSim){
  #Calculate initial 10 datapoints before performing first t-test, cannot do test on n = 1 and small samples vary widely.
  for(i in 1:10){ #for each simulated experiment
    x[i]<-rnorm(n = 1, mean = 0, sd = SD)
    y[i]<-rnorm(n = 1, mean = D, sd = SD) #Use D as a difference score, because sd=1 equals Cohen's D
    meanX<-mean(x[1:i])
    meanY<-mean(y[1:i])
    sdX<-sd(x[1:i])
    sdY<-sd(y[1:i])
    observedD<-(meanY-meanX)/(sqrt((((length(x[1:i]) - 1)*((sdX^2))) + (length(y[1:i]) - 1)*((sdY^2)))/((length(x[1:i])+length(y[1:i])-2))))#calculate observed D
    dlist[i]<-observedD
  }  

#After the first 10, I perform one-sided statistical tests after every participant. 
  for(i in 11:n){ #for each simulated participants after the first 10
    x[i]<-rnorm(n = 1, mean = 0, sd = SD)
    y[i]<-rnorm(n = 1, mean = D, sd = SD)
    z<-t.test(x[1:i],y[1:i], alternative = "less", var.equal=TRUE) #perform the t-test
    meanX<-mean(x[1:i])
    meanY<-mean(y[1:i])
    sdX<-sd(x[1:i])
    sdY<-sd(y[1:i])
    observedD<-(meanY-meanX)/(sqrt((((length(x[1:i]) - 1)*((sdX^2))) + (length(y[1:i]) - 1)*((sdY^2)))/((length(x[1:i])+length(y[1:i])-2))))
    dlist[i]<-observedD
    dmat[j,i]<-observedD #store all effect sizes for all simulations
    p[i]<-z$p.value 
    pmat[j,i]<-z$p.value #store all p-values for all simulations
  }
  
  #n*2.5 test
  
  meanX<-mean(x[1:(ceiling(n_orig*2.5))])
  meanY<-mean(y[1:(ceiling(n_orig*2.5))])
  sdX<-sd(x[1:(ceiling(n_orig*2.5))])
  sdY<-sd(y[1:(ceiling(n_orig*2.5))])
  observedD<-(meanY-meanX)/(sqrt((((length(x[1:(ceiling(n_orig*2.5))]) - 1)*((sdX^2))) + (length(y[1:(ceiling(n_orig*2.5))]) - 1)*((sdY^2)))/((length(x[1:(ceiling(n_orig*2.5))])+length(y[1:(ceiling(n_orig*2.5))])-2))))
  z<-t.test(x[1:(ceiling(n_orig*2.5))],y[1:(ceiling(n_orig*2.5))], alternative = c("less"), var.equal=TRUE) #perform the t-test
  p_n2.5[j]<-z$p.value
  if (z$p.value<alphalookn25){
    lookn25sig[j]<-1
  } else {
    lookn25sig[j]<-0
  }
  
  #Perform tests at each look, store p-value
  for (k in 1:Looks){
    meanX<-mean(x[1:LookN[k]])  
    meanY<-mean(y[1:LookN[k]])
    sdX<-sd(x[1:LookN[k]])
    sdY<-sd(y[1:LookN[k]])
    ObsD<-(meanY-meanX)/(sqrt((((LookN[k] - 1)*((sdX^2))) + (LookN[k] - 1)*((sdY^2)))/((LookN[k]+LookN[k]-2))))
    z<-t.test(x[1:LookN[k]],y[1:LookN[k]], alternative = c("less"), var.equal=TRUE) #perform the t-test
    mat1[j,k]<-z$p.value
    matdobs[j,k]<-ObsD
    bf <- exp(ttest.tstat(z$statistic, n1=LookN[k], n2=LookN[k], nullInterval = c(0, -Inf), rscale = rscaleBF)[['bf']])
    matBF[j,k]<-bf
  }
}


#Save only signigicant studies with N 
for (i in 1:Looks){
  mat2[,i] <- ifelse(mat1[,i]<alphalook[i],LookN[i],NA)
}

#count number of significant studies
SigSeq<-numeric(Looks)
for (i in 1:Looks){
  SigSeq[i]<-nSim-sum(is.na(mat2[,i]))
}

#Then replace last columns NA with max sample size to determine overal max sample size
mat2[,Looks] <- LookN[Looks]

#Determine the total sample size for all the sequential analyses (stopping when significant)
MinSample<-numeric(nSim)
for (i in 1:nSim){
  MinSample[i]<-min(mat2[i,], na.rm=TRUE)
}

#Count How Many studies support H0 (with Bayes < 0.3) 
for (i in 1:Looks){
  supportH0[i] <- sum(matBF[,i]<0.3)
}

#Count How Many studies support H1 (with Bayes > 3) 
for (i in 1:Looks){
  supportH1[i] <- sum(matBF[,i]>3)
}

#Determine minimum sample size N*2.5
sum(nSim*n_orig*2.5)
#significant studies using n*2.5
sum(lookn25sig)

#function to plot power (one-sided)
Power <- (function(D, n)
{
  ncp <- D*(n*n/(n+n))^0.5 #formula to calculate t from d from Dunlap, Cortina, Vaslow, & Burke, 1996, Appendix B
  t <- qt(0.95,df=n-2)
  1-(pt(t,df=n-2,ncp=ncp)-pt(-t,df=n-2,ncp=ncp))
}
)

#functions to plot 95% CI around true effect size d   
a<-0.05 #set alpha
CIplotU <- Vectorize(function(Dtrue,n, a)
{
  ci_u_d<-ci.smd(smd = Dtrue, n.1 = n, n.2 = n, conf.level=1-a)$Upper.Conf.Limit.smd
}
)

CIplotL <- Vectorize(function(Dtrue,n, a)
{
  ci_l_d<-ci.smd(smd = Dtrue, n.1 = n, n.2 = n, conf.level=1-a)$Lower.Conf.Limit.smd
}
)

#Create the plot
png(file="mygraphic.png",width=1500,height=1250, res = 150)
plot.new()
par(mfrow=c(2,1))
par(mar=c(4.0, 6.0, 0.5, 1.5)) 
plot(0, col="red", lty=1, lwd=3, ylim=c(0,1), xlim=c(20,n), type="l", xlab='', ylab='observed p-value \n & a-priori power')
for (i in 1:nSim){
  lines(pmat[i,], xlim=c(20,n), lwd=2)
}
abline(h=0.05, col="red", lty=2, lwd=2)
abline(h=0.8, col="forestgreen", lty=2, lwd=2)
for (i in 1:Looks){
  abline(v=LookN[i], col=i)
}
abline(v=ceiling(n_orig*2.5), col="dodgerblue", lwd=2)
par(new=TRUE)
curve(Power(D=D, n=x), 3, n, type="l", lty=1, lwd=3, ylim=c(0,1), xlim=c(20,n), col="forestgreen", xlab='', ylab='', axes = FALSE) 
plot(0, ylim=c(D-1,D+1), lty=1, type="l", xlim=c(20,n), xlab='sample size', ylab='Cohens d (observed, \n 95% CI, and true)')
for (i in 1:nSim){
  lines(dmat[i,], lwd=2)
}
abline(h=D, lty=3, lwd=3, col="gold")
for (i in 1:Looks){
  abline(v=LookN[i], col=i)
}
abline(v=ceiling(n_orig*2.5), col="dodgerblue", lwd=2)
par(new=TRUE)
curve(CIplotU(Dtrue=D, n=x, a=.05), 2, n, type="l", lty=3, lwd=3, ylim=c(D-1,D+1), xlim=c(20,n), col="gold", xlab='', ylab='', , axes = FALSE)
par(new=TRUE)
curve(CIplotL(Dtrue=D, n=x, a=.05), 2, n, type="l", lty=3, lwd=3, ylim=c(D-1,D+1), xlim=c(20,n), col="gold", xlab='', ylab='', axes = FALSE)
D2=0 #to draw 95 CI around 0
par(new=TRUE)
curve(CIplotU(Dtrue=0, n=x, a=0.04), 2, n, type="l", lty=4, lwd=2, ylim=c(D-1,D+1), xlim=c(20,n), col="brown", xlab='', ylab='', , axes = FALSE)
#par(new=TRUE)
#curve(CIplotL(Dtrue=0, n=x), 2, n, type="l", lty=4, lwd=2, ylim=c(D-1,D+1), xlim=c(20,n), col="brown", xlab='', ylab='', axes = FALSE)
abline(h=0, col="brown", ylim=c(D-1,D+1), lty=4, lwd=2)
par(mfrow=c(1,1))
dev.off()

#return significant studies based on N*2.5
sum(nSim*n_orig*2.5)
sum(lookn25sig)

#Return significant studies using sequential analyses
sum(MinSample)
SigSeq

#Return conclusions based on BF
supportH0
supportH1