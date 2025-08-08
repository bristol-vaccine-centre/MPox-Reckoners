setwd("~/OneDrive - University of Bristol/projects2022/MPX/poxreckoners")

natsal=read.delim("~/OneDrive - University of Bristol/projects2022/MPX/UKDA-7799-tab/tab/eul_natsal_2010_for_archive.tab")
natsalCOVIDwave2=read.delim("~/OneDrive - University of Bristol/projects2022/MPX/UKDA-8865-tab/tab/natsal_covid_wave_2_archive.tab")

OUTPUTFIGS=0
VI = 0.0
VT=0.0
mult1=1
attackrate = c(0.12,0.12,0.013) #from Packer et al. 
attackrate = mult1*c(0.12,0.12,0.12) #from Packer et al. 
myinterval=0.95
vaccinecoverage=100.0 #percentage vaccinated

myweightedmean=function(DD,ix)
{
  return(sum(DD[ix,1]*DD[ix,2]^2)/sum(DD[ix,1]*DD[ix,2]))
}

myweightedmeanvac=function(DD,ix)
{
  #return(sum(DD[ix,1]*DD[ix,2]*DD[ix,3])/sum(DD[ix,3]))
  #return(sum(DD[ix,1]*DD[ix,2]*DD[ix,3])
  return(sum(DD[ix,1]*DD[ix,2]))
}


resetMSM <- function(MSM,probvac1=NULL,ve_inf=0.0,ve_trans=0.0)
{
  if(!is.null(probvac1)){MSM$probvac=probvac1}
  
  MSM %>%
    filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
    mutate(totalcontacts = het1yr+sam1yr) %>%
    mutate(propsame = sam1yr/totalcontacts)%>%
    #mutate(probvac = probvac1) %>%
    mutate(degree = sam1yr + het1yr + hhsize) %>%
    mutate(degreedenom = sum(total_wt*degree)) %>%
    mutate(ci = total_wt*degree/degreedenom) %>%
    # mutate(wi_vac = probvac*total_wt) %>%
    # mutate(wi_unvac = (1-probvac)*total_wt) %>%
     mutate(wi_vac = probvac) %>%
     mutate(wi_unvac = (1-probvac)) %>%
    mutate(ci_vac = (1-ve_inf)*ci) %>%
    mutate(ci_unvac = ci) %>%
    mutate(totalAR = attackrate[1]*sam1yr + attackrate[2]*het1yr + attackrate[3]*hhsize) %>%
    mutate(Ri_vac = (1-ve_trans)*totalAR) %>%
    mutate(Ri_unvac = totalAR)%>%
    mutate(Rsquared = Ri_unvac*ci_unvac)%>%
    mutate(Rsquared2 = Ri_vac*ci_vac)-> MSM
  return(MSM)
}
calcRvac <- function(MSM)
{
  rindvac=cbind(MSM$Rsquared2,MSM$wi_vac)
  rindunvac = cbind(MSM$Rsquared,MSM$wi_unvac)
  rind = rbind(rindvac,rindunvac)
  myboot=boot(rind,statistic=myweightedmeanvac,R=1000)
  bci=boot.ci(myboot,type="perc",conf=0.95)
  return(c(mean(myboot$t),quantile(myboot$t,c(0.025,0.975))))
}




natsal %>%
  select(sin2,total_wt,dage,hhsize,het1yr,sam1yr,sam5yrs) -> RRdata


RRdata %>%
  filter(sam1yr <100 & het1yr<100 & hhsize<100) %>%
  mutate(vaccineeffect = 1) %>%
  mutate(indR = vaccineeffect*attackrate[1]*sam1yr + vaccineeffect*attackrate[2]*het1yr + vaccineeffect*attackrate[3]*hhsize) -> RRdata



RRdata %>%
  filter(sam5yrs>0) -> MSM

rind2=cbind(MSM$indR,MSM$total_wt)
myboot=boot(rind2,statistic=myweightedmean,R=500)
mult1*mean(myboot$t)
mult1*range(myboot$t)


mult1=1
deltaX = seq(0,100,5)
RpopMSM_target=matrix(NA,nrow=length(deltaX),ncol=3)
RpopMSM_random=matrix(NA,nrow=length(deltaX),ncol=3)



MSM = resetMSM(MSM, probvac1=1)  
mean(MSM$indR)
calcRvac(MSM)

#no vaccination
MSM = resetMSM(MSM, probvac1=0.0,ve_inf = 0.9,ve_trans =0.3)  
MSM %>%
  select(Rsquared,Rsquared2,wi_unvac,wi_vac,total_wt) -> temp1

head(temp1)

sum(temp1$Rsquared)
sum(temp1$wi_vac)
sum(temp1$Rsquared*temp1$wi_unvac) + sum(temp1$Rsquared2*temp1$wi_vac)
calcRvac(MSM)

sum(MSM$total_wt)
sum(MSM$wi_unvac) + sum(MSM$wi_vac)

rindvac=cbind(MSM$Ri_vac,MSM$ci_vac,MSM$wi_vac)
rindunvac = cbind(MSM$Ri_unvac,MSM$ci_unvac,MSM$wi_unvac)
rind = rbind(rindvac,rindunvac)
myboot=boot(rind,statistic=myweightedmeanvac,R=1000)

calcRvac(MSM)

#everyone vaccinated but effectiveness is zero (no effect)
MSM = resetMSM(MSM, probvac1=1,ve_inf=0, ve_trans = 0)  
MSM %>%
  select(Rsquared,Rsquared2,wi_unvac,wi_vac,total_wt) -> temp1

sum(temp1$Rsquared)
sum(temp1$Rsquared2)
sum(temp1$wi_vac)
sum(temp1$Rsquared*temp1$wi_unvac) + sum(temp1$Rsquared2*temp1$wi_vac)
calcRvac(MSM)
calcRvac(MSM)

#everyone vaccinated but effectiveness is zero (no effect)
MSM = resetMSM(MSM, probvac1=0.2,ve_inf=0, ve_trans = 0)  
calcRvac(MSM)

#everyone vaccinated and effectiveness against infection is 50%
MSM = resetMSM(MSM, probvac1=1,ve_inf=0.5, ve_trans = 0)  
calcRvac(MSM)

#50% vaccinated and effectiveness against infection is 50%
MSM = resetMSM(MSM, probvac1=0.5,ve_inf=0.5, ve_trans = 0)  
calcRvac(MSM)

#everyone vaccinated and effectiveness against transmission is 50%
MSM = resetMSM(MSM, probvac1=1,ve_inf=0.5, ve_trans = 0.5)  
calcRvac(MSM)
i=vaccinecoverage
j=1
VI=0.6
VT=0.3
for(i in deltaX)
{
  #random vaccination
  MSM = resetMSM(MSM, probvac1=i/100,ve_inf=VI,ve_trans=VT)  
  RpopMSM_random[j,1:3] =   calcRvac(MSM)
  
  
  #targeted vaccination
  #original version:
  # k = sum(MSM$total_wt)*i/(100*sum(MSM$totalcontacts*MSM$total_wt))
  # f=k*MSM$totalcontacts
  # sum(f*MSM$total_wt)/sum(MSM$total_wt)  #check: this should equal i/100
  # MSM$probvac <- f
  # MSM %>%
  #   mutate(probvac = ifelse(probvac>1,1,probvac)) -> MSM
  # sum(MSM$probvac*MSM$total_wt)/sum(MSM$total_wt) #warning - this is the actual fraction that were vaccinated
  MSM = resetMSM(MSM)  
  
  N=sum(MSM$total_wt) #total number of people (we want i/100 to be vaccinated)
  num2vac = sum(MSM$total_wt)*i/100
  #probvac1=num2vac*MSM$ci/max(MSM$ci) #needs to be strictly <=1
  
  O1=order(MSM$ci,decreasing = TRUE)
  CS1=cumsum(MSM$total_wt[O1])
  Vthreshold = min(which(CS1>num2vac))
  probvac1=rep(0,length(O1))
  probvac1[O1[1:Vthreshold]]=1
  MSM$probvac <- probvac1
  #sum(MSM$probvac*MSM$total_wt)/sum(MSM$total_wt) 
  
  MSM = resetMSM(MSM, ve_inf=VI,ve_trans=VT)  
  RpopMSM_target[j,1:3] =   calcRvac(MSM)
  j=j+1
}

if(OUTPUTFIGS==1)png('figs/MSM_cladeIIb.png',width=1.3*480,height=480)

mycols2=c(rgb(0,0.5,1,0.5),rgb(0.5,0.8,0.8),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.5,0.5,0.5))
mycols1=c(rgb(0,0.5,1,0.5),rgb(0.2,0.5,0.5),rgb(0,0.8,0,0.5),rgb(0.8,0.0,0.0,0.5),rgb(1,1,0,0.5),rgb(0.5,0.5,0.5))

par(mar=c(5,5,1,1),mfrow=c(1,1))
x=rev(deltaX)
y2=rev(RpopMSM_random[,3]); y1=rev(RpopMSM_random[,2]) 

plot(x,rev(RpopMSM_random[,1]),ylim=range(RpopMSM_random[,],0,RpopMSM_target[,]),pch=19,xlab="",col="white",
     ylab=expression(paste(R[t]," in individuals with same sex partnerships")),
     cex.lab=1.5,cex.axis=1.3,cex=1.0)

lines(x,rev(RpopMSM_random[,2]))
lines(x,rev(RpopMSM_random[,3]))

mtext("VE against transmission",side=1,line=2.5,cex=1.3)

grid()

abline(h=1,col="grey",lwd=3)

polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[2],border=NA)
lines(x,rev(RpopMSM_random[,1]),lty=2,col=mycols1[2],lwd=3)

y2=rev(RpopMSM_target[,3]); y1=rev(RpopMSM_target[,2]) 
polygon(c(x,rev(x)),c(y2,rev(y1)),col=mycols2[4],border=NA)
lines(x,rev(RpopMSM_target[,1]),lty=2,col=mycols1[4],lwd=3)

text(80,0.2,paste("VT =",VT,",VI =",VI))

legend('bottomleft',c('Random vaccination',
                      'Targeted vaccination'),
       col=mycols2[c(2,4)],pt.bg=mycols2[c(2,4)],pch=22,pt.cex=2,cex=1.0)
if(OUTPUTFIGS==1)dev.off()