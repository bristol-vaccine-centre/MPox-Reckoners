#natsal=read.delim("~/OneDrive - University of Bristol/Shared Documents//MPX/UKDA-7799-tab/tab/eul_natsal_2010_for_archive.tab")
library(tidyverse)

natsal=read.delim("../UKDA-7799-tab/tab/eul_natsal_2010_for_archive.tab")

natsalCOVIDwave2=read.delim("../UKDA-8865-tab/tab/natsal_covid_wave_2_archive.tab")


myweightedmean=function(DD,ix)
{
  return(sum(DD[ix,2]*DD[ix,1]^2)/sum(DD[ix,1]*DD[ix,2]))
}

colnames(natsal)

#het5yrs no. of het. sex partners, last 5 years

head(natsal$het5yrs)

#9999 == not answered
table(natsal$sam5yrs)

#table(natsal$overlp5y)
table(natsal$het1yr)

#9999 == not answered
table(natsal$sam1yr)


table(natsal$hhsize)
#9995: at least one if hhsize is not answered

summary(natsal$dage)

natsal %>%
  filter(het1yr<900) %>%
  group_by(dage) %>%
  summarise(meantot = mean(het1yr)) %>%
  ggplot(aes(x=dage,y=meantot)) +
  geom_point()

names(natsalCOVIDwave2)

natsal %>%
  filter(het1yr<900) %>%
  group_by(dage) %>%
  summarise(meantot = mean(het1yr)) %>%
  ggplot(aes(x=dage,y=meantot)) +
  geom_point()

natsal %>%
  filter(het1yr<900) %>%
  group_by(dage) %>%
  summarise(meantot = mean(total_wt)) %>%
  ggplot(aes(x=dage,y=meantot)) +
  geom_point()

summary(natsal$total_wt)
natsal %>%
  select(sin2,total_wt,dage,hhsize,het1yr,sam1yr) -> RRdata

attackrate = c(0.12,0.12,0.013)

RRdata %>%
  filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
  mutate(indR = attackrate[1]*sam1yr + attackrate[2]*het1yr + attackrate[3]*hhsize) -> RRdata


mult1=1
rind=cbind(RRdata$indR,RRdata$total_wt)
myboot=boot(rind,statistic=myweightedmean,R=500)
mult1*mean(myboot$t)
bci=boot.ci(myboot,type="basic",conf=myinterval)
mult1*bci$basic[1,4:5]
mult1*range(myboot$t)

sum(rind[,2]*rind[,1]^2)/sum(rind[,1])

summary(RRdata$indR)




RRdata %>%
  group_by(dage) %>%
  summarise(meantot = mean(indR)) %>%
  ggplot(aes(x=dage,y=meantot)) +
  geom_point()


RRdata %>%
  group_by(sam1yr) %>%
  summarise(meantot = mean(indR)) %>%
  ggplot(aes(x=sam1yr,y=meantot)) +
  geom_point()

RRdata %>%
  ggplot(aes(x=indR)) + 
  geom_histogram()


RRdata %>%
  filter(sam1yr>0) -> MSM


rind2=cbind(MSM$indR,MSM$total_wt)
myboot=boot(rind2,statistic=myweightedmean,R=500)
mult1*mean(myboot$t)
mult1*range(myboot$t)


head(MSM)

sum(MSM$sam1yr)/sum(MSM$sam1yr+MSM$het1yr)

MSM %>% 
  mutate(totalcontacts = het1yr+sam1yr) %>% 
  mutate(propsame = sam1yr/totalcontacts)->MSM

head(MSM)
MSM %>% 
  ggplot(aes(x=propsame))+
  geom_histogram()

MSM %>% 
  summary(propsame)

MSM %>% 
  ggplot(aes(x=totalcontacts,y=propsame))+
  geom_point()

MSM %>% 
  group_by(sam1yr) %>% 
  summarise(meanR = mean(indR)) %>% 
  ggplot(aes(x=sam1yr,y=meanR))+
  geom_point()



totalhet = sum(RRdata$het1yr)
totalsame = sum(RRdata$sam1yr)
num1=length(RRdata$het1yr)
mult1=1 #use this later if we need to scale contacts
deltaX = 1:100
Rpop=matrix(NA,nrow=length(deltaX),ncol=3)
for(i in deltaX)
{
  #calculate number of contacts to reduce
  contactra = RRdata$het1yr
  totalcontacts=sum(contactra)
  
  #distribute them randomly across all people (not accounting for weights here!)
  r1=runif(length(contactra),min=0,max=contactra)
  r1=(1-i/100)*totalcontacts*r1/sum(r1)
  
  new1=c(round(r1[1]),diff(round(cumsum(r1))))

  RRdata2 = RRdata
  RRdata2$het1yr = new1
  
  #recalculate individual R numbers
  RRdata2 %>%
    filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
    mutate(indR = attackrate[1]*sam1yr + attackrate[2]*het1yr + attackrate[3]*hhsize) -> RRdata2
  
  #extract individual R numbers
  rind=cbind(RRdata2$indR,RRdata2$total_wt)
  myboot=boot(rind,statistic=myweightedmean,R=2000)
  Rpop[i,1] = mult1*mean(myboot$t)
  Rpop[i,2:3] = mult1*range(myboot$t)
}

deltaX = 1:100
RpopSAME=matrix(NA,nrow=length(deltaX),ncol=3)
for(i in deltaX)
{
  #calculate number of contacts to reduce
  contactra = RRdata$sam1yr
  totalcontacts=sum(contactra)
  
  #distribute them randomly across all people (not accounting for weights here!)
  r1=runif(length(contactra),min=0,max=contactra)
  r1=(1-i/100)*totalcontacts*r1/sum(r1)
  
  new1=c(round(r1[1]),diff(round(cumsum(r1))))
  
  RRdata2 = RRdata
  RRdata2$sam1yr = new1
  
  #recalculate individual R numbers
  RRdata2 %>%
    filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
    mutate(indR = attackrate[1]*sam1yr + attackrate[2]*het1yr + attackrate[3]*hhsize) -> RRdata2
  
  #extract individual R numbers
  rind=cbind(RRdata2$indR,RRdata2$total_wt)
  myboot=boot(rind,statistic=myweightedmean,R=500)
  RpopSAME[i,1] = mult1*mean(myboot$t)
  RpopSAME[i,2:3] = mult1*range(myboot$t)
}

#plot
mycols2=c(rgb(0,0.5,1,0.5),rgb(0,0.8,0.8,0.5),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.3,0.3,0.3,0.6))

par(mar=c(5,5,1,1),mfrow=c(1,1))
x=deltaX
y2=rev(Rpop[,3]); y1=rev(Rpop[,2]) 
plot(x,Rpop[,1],ylim=range(Rpop[,],0),pch=19,xlab="",col="white",
     ylab=expression(paste("Reproduction number, ",R[t])),cex.lab=1.5,cex.axis=1.3,cex=1.0)

mtext("% active sexual contacts",side=1,line=2.5,cex=1.3)

grid()

abline(h=1,col="grey",lwd=3)

polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[6],border=NA)


y2=rev(RpopSAME[,3]); y1=rev(RpopSAME[,2])
polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[5],border=NA)

#test
# mycontacts = round(runif(10,min=0,max=100))
# mycontacts
# totalcontacts=sum(mycontacts)
# 
# r1=runif(10,min=0,max=mycontacts)
# r1=0.5*totalcontacts*r1/sum(r1)
# new1=c(round(r1[1]),diff(round(cumsum(r1))))


deltaX = 1:100
RpopMSM=matrix(NA,nrow=length(deltaX),ncol=3)
numreps=10
for(i in deltaX)
{
  #calculate number of contacts to reduce
  contactra = MSM$sam1yr
  totalcontacts=sum(contactra)
  
  contactrahet = MSM$het1yr
  totalcontactshet=sum(contactrahet)
  
  rind=c()
  
  for(rx in 1:numreps)
  {
    #distribute them randomly across all people (not accounting for weights here!)
    r1=runif(length(contactra),min=0,max=contactra)
    r1=(1-i/100)*totalcontacts*r1/sum(r1)
  
    new1=c(round(r1[1]),diff(round(cumsum(r1))))
    #new1 = c(round((1-i/100)*contactra[1]),round(diff(cumsum(((1-i/100)*contactra)))))
  
    MSM2 = MSM
    MSM2$sam1yr = new1
    
    r1=runif(length(contactrahet),min=0,max=contactrahet)
    r1=(1-i/100)*totalcontactshet*r1/sum(r1)
    
    new1=c(round(r1[1]),diff(round(cumsum(r1))))
    
    MSM2$het1yr = new1
  
    #recalculate individual R numbers
    MSM2 %>%
      filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
      mutate(indR = attackrate[1]*sam1yr + attackrate[2]*het1yr + attackrate[3]*hhsize) -> MSM2
  
    #extract individual R numbers
    rind1=cbind(MSM2$indR,MSM2$total_wt)
    
    rind=rbind(rind,rind1)  
  }
  
  myboot=boot(rind,statistic=myweightedmean,R=2000)
  bci=boot.ci(myboot,type="perc",conf=myinterval)
  
  RpopMSM[i,1] = mult1*mean(myboot$t)
  #RpopMSM[i,2:3] = mult1*range(myboot$t)
  RpopMSM[i,2:3] = mult1*bci$percent[1,4:5]
}

mycols2=c(rgb(0,0.5,1,0.5),rgb(0,0.8,0.8,0.5),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.3,0.3,0.3,0.6))

par(mar=c(5,5,1,1),mfrow=c(1,1))
x=deltaX
y2=rev(RpopMSM[,3]); y1=rev(RpopMSM[,2]) 

plot(x,rev(RpopMSM[,1]),ylim=range(RpopMSM[,],0),pch=19,xlab="",col="red",
     ylab=expression(paste("Reproduction number, ",R[t])),cex.lab=1.5,cex.axis=1.3,cex=1.0)

lines(x,rev(RpopMSM[,2]))
lines(x,rev(RpopMSM[,3]))

mtext("% active same sex sexual contacts",side=1,line=2.5,cex=1.3)

grid()

abline(h=1,col="grey",lwd=3)

polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[6],border=NA)


#####eliminating number of contacts over X

deltaX = 1:max(MSM$sam1yr+MSM$het1yr)
RpopMSM=matrix(NA,nrow=length(deltaX),ncol=3)
numreps=1
for(i in deltaX)
{

  rind=c()
  
  for(rx in 1:numreps)
  {
    #MSM %>%
    #  mutate(sam1yr = min(sam1yr,i)) %>%
    #  mutate(het1yr = min(het1yr ,i)) -> MSM2

    
    #recalculate individual R numbers
    MSM2 %>%
      filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
      mutate(indR = attackrate[1]*sam1yr + attackrate[2]*het1yr + attackrate[3]*hhsize) -> MSM2
    
    #extract individual R numbers
    rind1=cbind(MSM2$indR,MSM2$total_wt)
    
    rind=rbind(rind,rind1)  
  }
  
  myboot=boot(rind,statistic=myweightedmean,R=2000)
  bci=boot.ci(myboot,type="perc",conf=myinterval)
  
  RpopMSM[i,1] = mult1*mean(myboot$t)
  #RpopMSM[i,2:3] = mult1*range(myboot$t)
  RpopMSM[i,2:3] = mult1*bci$percent[1,4:5]
}

mycols2=c(rgb(0,0.5,1,0.5),rgb(0,0.8,0.8,0.5),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.3,0.3,0.3,0.6))

par(mar=c(5,5,1,1),mfrow=c(1,1))
x=deltaX
y2=rev(RpopMSM[,3]); y1=rev(RpopMSM[,2]) 

plot(x,rev(RpopMSM[,1]),ylim=range(RpopMSM[,],0),pch=19,xlab="",col="red",
     ylab=expression(paste("Reproduction number, ",R[t])),cex.lab=1.5,cex.axis=1.3,cex=1.0)

lines(x,rev(RpopMSM[,2]))
lines(x,rev(RpopMSM[,3]))

mtext("% active same sex sexual contacts",side=1,line=2.5,cex=1.3)

grid()

abline(h=1,col="grey",lwd=3)

polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[6],border=NA)
