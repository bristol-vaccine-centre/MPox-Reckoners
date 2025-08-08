library(tidyverse)
library(boot)
library(cowplot)

setwd("~/OneDrive - University of Bristol/projects2022/MPX/poxreckoners")

natsal=read.delim("~/OneDrive - University of Bristol/projects2022/MPX/UKDA-7799-tab/tab/eul_natsal_2010_for_archive.tab")
natsalCOVIDwave2=read.delim("~/OneDrive - University of Bristol/projects2022/MPX/UKDA-8865-tab/tab/natsal_covid_wave_2_archive.tab")


myweightedmean=function(DD,ix)
{
  return(sum(DD[ix,2]*DD[ix,1]^2)/sum(DD[ix,1]*DD[ix,2]))
}


natsal %>%
  select(sin2,total_wt,dage,hhsize,het1yr,sam1yr,sam5yrs) -> RRdata

attackrate = c(0.12,0.12,0.013) #from Packer et al. 

RRdata %>%
  filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
  mutate(vaccineeffect = 1) %>%
  mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> RRdata



RRdata %>%
  filter(sam5yrs>0) -> MSM

####how to calculate Rt: 
rind2=cbind(MSM$indR,MSM$total_wt)
myboot=boot(rind2,statistic=myweightedmean,R=500)
mult1*mean(myboot$t)
mult1*range(myboot$t)




MSM %>% 
  mutate(totalcontacts = het1yr+sam1yr) %>% 
  mutate(propsame = sam1yr/totalcontacts)->MSM



#what about a hard to reach population? Randomly distributed through the population 
#Vaccination
mult1=1
deltaX = seq(100,0,-5)
RpopMSM_target=matrix(NA,nrow=length(deltaX),ncol=3)
RpopMSM_random=matrix(NA,nrow=length(deltaX),ncol=3)
numreps=10
j=1
VE=0.8
numMSM = length(MSM$total_wt)

MSM %>%
  filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
  mutate(vaccineeffect = 1) %>% #this is whether the vaccine is effective in this person
  mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> MSM



rind2=cbind(MSM$indR,MSM$total_wt)
myboot=boot(rind2,statistic=myweightedmean,R=500)
print(paste("R = ",mult1*mean(myboot$t)))
print(mult1*range(myboot$t))

for(i in deltaX)
{
  MSM$vaccineeffect=1
  #print(paste("i=",i))
  #calculate number of individuals to vaccinate 
  num2vac = round(numMSM*i/100)
  #print(paste("num to vaccinated=",num2vac))


  rind_rand=c()
  rind_target=c()
  
  for(rx in 1:numreps)
  {
    MSM$vaccineeffect=1
    MSM %>%
      mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> MSM
    
    MSM2 = MSM
    
    rand1 = sample(1:numMSM, num2vac, replace=F)
    MSM$vaccineeffect[rand1] <- 1 - VE
    
    MSM %>%
      mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> MSM
    
    #rind2=cbind(MSM$indR,MSM$total_wt)
    #myboot=boot(rind2,statistic=myweightedmean,R=500)
    #print(mult1*mean(myboot$t))
    
    #extract individual R numbers
    rind1=cbind(MSM$indR,MSM$total_wt)
    rind_rand=rbind(rind_rand,rind1) 
    
    rand2 = sample(1:numMSM, num2vac, replace=F, prob = (MSM$totalcontacts+1)/sum(MSM$totalcontacts+numMSM))
    MSM2$vaccineeffect[rand2] = 1 - VE

    MSM2 %>%
      mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> MSM2
    
    rind1=cbind(MSM2$indR,MSM2$total_wt)
    rind_target=rbind(rind_target,rind1) 
  }
  
  myboot=boot(rind_rand,statistic=myweightedmean,R=1000)
  bci=boot.ci(myboot,type="perc",conf=myinterval)
  
  RpopMSM_random[j,1] = mult1*mean(myboot$t)
  RpopMSM_random[j,2:3] = mult1*bci$percent[1,4:5]
  
  #print(mult1*mean(myboot$t))
  
  
  myboot=boot(rind_target,statistic=myweightedmean,R=1000)
  bci=boot.ci(myboot,type="perc",conf=myinterval)
  
  RpopMSM_target[j,1] = mult1*mean(myboot$t)
  RpopMSM_target[j,2:3] = mult1*bci$percent[1,4:5]

  j=j+1
}

png('figs/MSM_cladeIIb.png',width=1.3*480,height=480)
mycols2=c(rgb(0,0.5,1,0.5),rgb(0.5,0.8,0.8),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.5,0.5,0.5))
mycols1=c(rgb(0,0.5,1,0.5),rgb(0.2,0.5,0.5),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.5,0.5,0.5))

par(mar=c(5,5,1,1),mfrow=c(1,1))
x=rev(deltaX)
y2=rev(RpopMSM_random[,3]); y1=rev(RpopMSM_random[,2]) 

plot(x,rev(RpopMSM_random[,1]),ylim=range(RpopMSM_random[,],0),pch=19,xlab="",col="white",
     ylab=expression(paste(R[t]," in individuals with same sex partnerships")),cex.lab=1.5,cex.axis=1.3,cex=1.0)

lines(x,rev(RpopMSM_random[,2]))
lines(x,rev(RpopMSM_random[,3]))

mtext("% individuals vaccinated",side=1,line=2.5,cex=1.3)

grid()

abline(h=1,col="grey",lwd=3)

polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)
lines(x,rev(RpopMSM_random[,1]),lty=2,col=mycols1[2],lwd=3)

y2=rev(RpopMSM_target[,3]); y1=rev(RpopMSM_target[,2]) 
polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)
lines(x,rev(RpopMSM_target[,1]),lty=2,col=mycols1[4],lwd=3)

legend('bottomleft',c('Random vaccination',
                   'Targeted vaccination'),
       col=mycols2[c(2,4)],pt.bg=mycols2[c(2,4)],pch=22,pt.cex=2,cex=1.0)
dev.off()

png('figs/MSM_cladeIIb_percent.png',width=1.3*480,height=480)
par(mar=c(5,5,1,1),mfrow=c(1,1))
baseline=rev(RpopMSM_random[,1])[1]

aspercent <- function(R)
{
  P1 = 100*(R-baseline)/baseline
  return(P1)
}

x=rev(deltaX)
y2=aspercent(rev(RpopMSM_random[,3])); y1=aspercent(rev(RpopMSM_random[,2]) )

plot(x,aspercent(rev(RpopMSM_random[,1])),ylim=range(aspercent(RpopMSM_random[,]),0),pch=19,xlab="",col="white",
     ylab=expression(paste("% change in ",R[t])),cex.lab=1.5,cex.axis=1.3,cex=1.0)

lines(x,aspercent(rev(RpopMSM_random[,2])))
lines(x,aspercent(rev(RpopMSM_random[,3])))

mtext("% individuals vaccinated",side=1,line=2.5,cex=1.3)

grid()

polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)

lines(x,aspercent(rev(RpopMSM_random[,1])),lty=2,col=mycols1[2],lwd=3)

y2=aspercent(rev(RpopMSM_target[,3])); y1=aspercent(rev(RpopMSM_target[,2]) )
polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)
lines(x,aspercent(rev(RpopMSM_target[,1])),lty=2,col=mycols1[4],lwd=3)
legend('bottomleft',c('Random vaccination',
                      'Targeted vaccination'),
       col=mycols2[c(2,4)],pt.bg=mycols2[c(2,4)],pch=22,pt.cex=2,cex=1.0)
dev.off()



deltaX = seq(100,0,-5)
RpopMSM_target=matrix(NA,nrow=length(deltaX),ncol=3)
RpopMSM_random=matrix(NA,nrow=length(deltaX),ncol=3)
numreps=10
j=1
VE=0.8

numMSM = length(MSM$total_wt)
attackrate = c(0.12,0.12,0.13) #from Packer et al. 

attackrate = c(0.12,0.12,0.3) #from Packer et al. 

MSM %>%
  filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
  mutate(vaccineeffect = 1) %>%
  mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> MSM

summary(MSM$indR)

rind2=cbind(MSM$indR,MSM$total_wt)
myboot=boot(rind2,statistic=myweightedmean,R=5000)
print(paste("R = ",mult1*mean(myboot$t)))
print(mult1*range(myboot$t))


for(i in deltaX)
{
  MSM$vaccineeffect=1
  #print(paste("i=",i))
  #calculate number of individuals to vaccinate 
  num2vac = round(numMSM*i/100)
  #print(paste("num to vaccinated=",num2vac))
  
  
  rind_rand=c()
  rind_target=c()
  
  for(rx in 1:numreps)
  {
    MSM$vaccineeffect=1
    MSM %>%
      mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> MSM
    
    MSM2 = MSM
    
    rand1 = sample(1:numMSM, num2vac, replace=F)
    MSM$vaccineeffect[rand1] <- 1 - VE
    
    MSM %>%
      mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> MSM
    
    #rind2=cbind(MSM$indR,MSM$total_wt)
    #myboot=boot(rind2,statistic=myweightedmean,R=500)
    #print(mult1*mean(myboot$t))
    
    #extract individual R numbers
    rind1=cbind(MSM$indR,MSM$total_wt)
    rind_rand=rbind(rind_rand,rind1) 
    
    rand2 = sample(1:numMSM, num2vac, replace=F, prob = (MSM$totalcontacts+1)/sum(MSM$totalcontacts+numMSM))
    MSM2$vaccineeffect[rand2] = 1 - VE
    
    MSM2 %>%
      mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> MSM2
    
    rind1=cbind(MSM2$indR,MSM2$total_wt)
    rind_target=rbind(rind_target,rind1) 
  }
  
  myboot=boot(rind_rand,statistic=myweightedmean,R=1000)
  bci=boot.ci(myboot,type="perc",conf=myinterval)
  
  RpopMSM_random[j,1] = mult1*mean(myboot$t)
  RpopMSM_random[j,2:3] = mult1*bci$percent[1,4:5]
  
  #print(mult1*mean(myboot$t))
  
  
  myboot=boot(rind_target,statistic=myweightedmean,R=1000)
  bci=boot.ci(myboot,type="perc",conf=myinterval)
  
  RpopMSM_target[j,1] = mult1*mean(myboot$t)
  RpopMSM_target[j,2:3] = mult1*bci$percent[1,4:5]
  
  j=j+1
}

png('figs/MSM_cladeIb.png',width=1.3*480,height=480)
mycols2=c(rgb(0,0.5,1,0.5),rgb(0.5,0.8,0.8),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.5,0.5,0.5))
mycols1=c(rgb(0,0.5,1,0.5),rgb(0.2,0.5,0.5),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.5,0.5,0.5))

par(mar=c(5,5,1,1),mfrow=c(1,1))

x=rev(deltaX)
y2=rev(RpopMSM_random[,3]); y1=rev(RpopMSM_random[,2]) 

plot(x,rev(RpopMSM_random[,1]),ylim=range(RpopMSM_random[,],0),pch=19,xlab="",col="white",
     ylab=expression(paste(R[t]," in individuals with same sex partnerships")),cex.lab=1.5,cex.axis=1.3,cex=1.0)

lines(x,rev(RpopMSM_random[,2]))
lines(x,rev(RpopMSM_random[,3]))

mtext("% individuals vaccinated",side=1,line=2.5,cex=1.3)

grid()

abline(h=1,col="grey",lwd=3)

polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)
lines(x,rev(RpopMSM_random[,1]),lty=2,col=mycols1[2],lwd=3)

y2=rev(RpopMSM_target[,3]); y1=rev(RpopMSM_target[,2]) 
polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)
lines(x,rev(RpopMSM_target[,1]),lty=2,col=mycols1[4],lwd=3)

legend('bottomleft',c('Random vaccination',
                      'Targeted vaccination'),
       col=mycols2[c(2,4)],pt.bg=mycols2[c(2,4)],pch=22,pt.cex=2,cex=1.0)

dev.off()

png('figs/MSM_cladeIb_percent.png',width=1.3*480,height=480)
par(mar=c(5,5,1,1),mfrow=c(1,1))
baseline=rev(RpopMSM_random[,1])[1]

aspercent <- function(R)
{
  P1 = 100*(R-baseline)/baseline
  return(P1)
}

x=rev(deltaX)
y2=aspercent(rev(RpopMSM_random[,3])); y1=aspercent(rev(RpopMSM_random[,2]) )

plot(x,aspercent(rev(RpopMSM_random[,1])),ylim=range(aspercent(RpopMSM_random[,]),0),pch=19,xlab="",col="white",
     ylab=expression(paste("% change in ",R[t])),cex.lab=1.5,cex.axis=1.3,cex=1.0)

lines(x,aspercent(rev(RpopMSM_random[,2])))
lines(x,aspercent(rev(RpopMSM_random[,3])))

mtext("% individuals vaccinated",side=1,line=2.5,cex=1.3)

grid()

polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)

lines(x,aspercent(rev(RpopMSM_random[,1])),lty=2,col=mycols1[2],lwd=3)

y2=aspercent(rev(RpopMSM_target[,3])); y1=aspercent(rev(RpopMSM_target[,2]) )
polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)
lines(x,aspercent(rev(RpopMSM_target[,1])),lty=2,col=mycols1[4],lwd=3)

legend('bottomleft',c('Random vaccination',
                      'Targeted vaccination'),
       col=mycols2[c(2,4)],pt.bg=mycols2[c(2,4)],pch=22,pt.cex=2,cex=1.0)

dev.off()







attackrate = c(0.12,0.12,0.013) #from Packer et al. 
hhar=seq(0,0.2,0.01)
Rt=matrix(NA,nrow=length(hhar),ncol=3)
k=1
for(ar in hhar)
{
  MSM %>%
    filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
    filter(sam1yr>=0 & het1yr>=0 & hhsize>=0) %>%
    mutate(vaccineeffect = 1) %>%
    mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*ar*hhsize) %>%
    mutate(riskinf = vaccineeffect*attackrate[1]*(sam1yr-1) + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*ar*hhsize)-> MSM
  
  
    Rt[k,1] = sum(MSM$total_wt*MSM$indR*MSM$riskinf)/sum(MSM$indR)
  
  #rind2=cbind(MSM$indR,MSM$total_wt)
  #myboot=boot(rind2,statistic=myweightedmean,R=5000)
  #Rt[k,1] = mean(myboot$t)
  #Rt[k,2:3] = range(myboot$t)
  k=k+1
}



y2=Rt[,3]; y1=Rt[,2]
plot(hhar,Rt[,1],lty=2,col=mycols1[4],lwd=3)



#polygon(c(hhar,rev(hhar)),c(y2,rev(y1)),col=mycols2[4],border=NA)
#lines(x,aspercent(rev(RpopMSM_target[,1])),lty=2,col=mycols1[4],lwd=3)
