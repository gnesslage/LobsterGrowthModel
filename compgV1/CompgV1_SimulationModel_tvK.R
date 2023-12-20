#Basic simulation model to test performance of Comprehensive Growth Model Version 1
#Genny Nesslage & Mike Wilberg Dec 20, 2023

#Setup
setwd("C:/compgV1")
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)

graphics.off()    
rm(list=ls(all=TRUE))
source("CompgV1_functions.R")

#### SIMULATION SET UP ####
#Time step setup
nyrs=6 #Number of years
nseas=4 #Number of seasons in a year
remainder=1 #modulo remainder = argument setting season with R pulse
ts=nseas*nyrs #Number of time steps

#Specify population info and mortality
fage=1 #First age
lage=6 #Last age
nages=lage-fage+1 #Number of ages
M=as.vector(rep(0.3/nages,nages)/nseas) #0.2/nages #natural mortality rate
#M=as.vector(c(0.0125,0.0125,0.0125,0.0125,0.025,0.025)/nseas) #0.2/nages #natural mortality rate
F=0.0/nseas #Seasonal fishing mortality rate
nsex=2 #Number of sexes
fsexratio=0.5 #Proportion female
Rmean=1000 #Mean recruitment
Rsd=0 #SD of recruitment
Ninit=2000 #Initial pop size

#Specify growth
startsize=12.5 #Mean length for R pulse in sd simulation
#temps=rep(20,ts) #Average daily temp in each time step constant
temps=c(20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,21,21)
GDD<-ifelse(temps>5.0,(temps-5.0)*90,0) #GDD per season given 5deg min
fkp=0.0002222/4 #Female proportionality constant to convert GDD to K
mkp=0.0002222/4 #Male proportionality constant to convert GDD to K
fK=(fkp*GDD) #Female vonB K
mK=(mkp*GDD) #Male vonB K
fkCV=0.2 #Female CV of vonB K
mkCV=0.2 #Male CV of vonB K
fKSd=fK[1]*fkCV #Female SD of vonB K
mKSd=mK[1]*mkCV #Male SD of vonB K
fLinf=250 #Female vonB Linf
mLinf=250 #Male vonB Linf
fLinfCV=0.1 #Female CV of vonB Linf
mLinfCV=0.1 #Male CV of vonB Linf
fLinfSd=fLinf*fLinfCV #Female SD of vonB Linf
mLinfSd=mLinf*mLinfCV #Male SD of vonB Linf
fRho=0.01 #Female rho (cov bet vonB Linf and K)
mRho=0.01 #Male rho (cov bet vonB Linf and K)

#Set up length bins
fbin=10 #First bin
lbin=1.01*max(fLinf,mLinf) #Last bin (make sure last bin close to Linf)
binwidth=5 #Binwidth
midbins=seq(fbin+binwidth/2,lbin+binwidth/2,binwidth) #Vector of middle of each bin
nbins=length(midbins) #Number of bins

#Set up observation model
LFAsamps=0.5 #Proportion of population sampled by LFA in obs model
TAGsamps=0.1 #Proportion of population tagged in obs model
fEffN = 20 #effective sample size for females in LFA
mEffN = 20 #effective sample size for males in LFA
fa=fage #first age observed in LFA sampling
la=lage #last age observed in LFA sampling
tfa=fage #first  age observed in tagging study
tla=lage #last age observed in tagging study

#### GENERATE SEASONAL GROWTH MATRICES ####
fincsd=0.25*midbins[1] #sd of initial length distribution for females
mincsd=0.25*midbins[1] #sd of initial length distribution for males

##grow out one pulse of recruits to max age
ages=seq(fage,lage) #Vector of ages
Z=F+M #Total mortality rate
N=Rmean #Set Ninit equal to recruitment pulse
fN=N*fsexratio #Divide population into females and males
mN=N-fN
fmL=dnorm(midbins,startsize,fincsd) #initial numbers at length 
mmL=dnorm(midbins,startsize,mincsd) #initial numbers at length 
fpopLF=(fmL/sum(fmL)) #normalized female initial LF
mpopLF=(mmL/sum(mmL)) #normalized male initial LF
fNAL=round(fpopLF*fN,0) #Female initial numbers at length
mNAL=round(mpopLF*mN,0) #Male initial numbers at length

#Calculate seasonal growth matrices
fG<-array(0,c(ts,nbins,nbins)) #Female GTMs
for(t in 1:(ts)){
  if(fK[t]<=0){
    fG[t,,]<-diag(nbins)
  }
  else{
    fG[t,,]<-growth_matrix(Linf=fLinf,K=fK[t],Linf_Var=fLinfSd^2,K_Var=fKSd^2,Linf_K_Cov=fRho*fKSd^2*fLinfSd^2,bin_mids=midbins)
  }#close else
} #close for

mG<-array(0,c(ts,nbins,nbins)) #Male GTMs
for(t in 1:(ts)){
  if(mK[t]<=0){
    mG[t,,]<-diag(nbins)
  }
  else{
    mG[t,,]<-growth_matrix(Linf=mLinf,K=mK[t],Linf_Var=mLinfSd^2,K_Var=mKSd^2,Linf_K_Cov=mRho*mKSd^2*mLinfSd^2,bin_mids=midbins)
  }#close else
} #close m

#use GTM to grow initial pulse of recruits out to max age 
fN<-matrix(0,nrow=ts,ncol=nbins) #Female N matrix
fN[1,]=fNAL #Populate N matrix with initial # at length
mN<-matrix(0,nrow=ts,ncol=nbins) #Male N matrix
mN[1,]=mNAL #Populate N matrix with initial # at length

#Populate remaining N matrix
for(t in 1:(ts-1)){
  #m=t%%nseas+1  #determine month using the modulus function
  fN[t+1,]=fG[t,,]%*%(fN[t,])
  mN[t+1,]=mG[t,,]%*%(mN[t,])
}
colnames(fN)=midbins;colnames(mN)=midbins

#Calculate mean and sd for all time steps
timeseries=seq(1,ts,1)
falldat=matrix(0,ts,2)
malldat=matrix(0,ts,2)
for (i in 1:nrow(fN)){
  ftmp2=as.matrix(cbind(midbins,fN[i,]))
  mtmp2=as.matrix(cbind(midbins,mN[i,]))
  ftmp=matrix(0,1,1)
  mtmp=matrix(0,1,1)
  for(j in 1:nrow(ftmp2)){
    if(ftmp2[j,2]>0) {ftmp3=as.vector(rep(ftmp2[j,1],ftmp2[j,2]))}
    else {ftmp3=0}
    if(mtmp2[j,2]>0) {mtmp3=as.vector(rep(mtmp2[j,1],mtmp2[j,2]))}
    else {mtmp3=0}
    ftmp=c(ftmp,ftmp3)
    mtmp=c(mtmp,mtmp3)
  }
  ftmp=ftmp[-1]
  mtmp=mtmp[-1]
  ftmp=subset(ftmp,ftmp>0)
  mtmp=subset(mtmp,mtmp>0)
  falldat[i,1]=mean(ftmp)
  malldat[i,1]=mean(mtmp)
  falldat[i,2]=sd(ftmp)
  malldat[i,2]=sd(mtmp)
}

#Summarize mean and sd for each length bin by season
fmvec=matrix(0,nseas,nages)
mmvec=matrix(0,nseas,nages)
fsdvec=matrix(0,nseas,nages)
msdvec=matrix(0,nseas,nages)
for(s in 1:nseas){
  getrows=timeseries[seq(s, length(timeseries), nseas)]
  fsm=falldat[getrows,1]
  msm=malldat[getrows,1]
  fsd=falldat[getrows,2]
  msd=malldat[getrows,2]
  fmvec[s,]=fsm
  mmvec[s,]=msm
  fsdvec[s,]=fsd
  msdvec[s,]=msd
}
colnames(fmvec)=ages;colnames(mmvec)=ages
colnames(fsdvec)=ages;colnames(msdvec)=ages

#Organize and plot mean length and sd/CV of length at age
fvissmLAA=data.frame(cbind(ages,t(fmvec)));colnames(fvissmLAA)=c("Age","S1","S2","S3","S4")
fvissmLAA_long <- gather(fvissmLAA, Season, mLAA, S1:S4,factor_key=TRUE)
ggplot(data=fvissmLAA_long,aes(x=Age,y=mLAA,colour=Season))+geom_line()+theme_bw(base_size = 20)+xlab("Age")+ylab("Female Length")
mvissmLAA=data.frame(cbind(ages,t(mmvec)));colnames(mvissmLAA)=c("Age","S1","S2","S3","S4")
mvissmLAA_long <- gather(mvissmLAA, Season, mLAA, S1:S4,factor_key=TRUE)
ggplot(data=mvissmLAA_long,aes(x=Age,y=mLAA,colour=Season))+geom_line()+theme_bw(base_size = 20)+xlab("Age")+ylab("Male Length")

vissfsdMLvec=data.frame(cbind(ages,t(fsdvec)));colnames(vissfsdMLvec)=c("Age","S1","S2","S3","S4")
vissfsdMLvec_long <- gather(vissfsdMLvec, Season, sdLAA, S1:S4,factor_key=TRUE)
vissfsdMLvec_long$CV=vissfsdMLvec_long$sdLAA/fvissmLAA_long$mLAA
vissmsdMLvec=data.frame(cbind(ages,t(msdvec)));colnames(vissmsdMLvec)=c("Age","S1","S2","S3","S4")
vissmsdMLvec_long <- gather(vissmsdMLvec, Season, sdLAA, S1:S4,factor_key=TRUE)
vissmsdMLvec_long$CV=vissmsdMLvec_long$sdLAA/mvissmLAA_long$mLAA
ggplot(data=vissfsdMLvec_long,aes(x=Age,y=sdLAA,colour=Season))+geom_line()+theme_bw(base_size = 20)+xlab("Age")+ylab("SD Female Length")
ggplot(data=vissmsdMLvec_long,aes(x=Age,y=sdLAA,colour=Season))+geom_line()+theme_bw(base_size = 20)+xlab("Age")+ylab("SD Male Length")
ggplot(data=vissfsdMLvec_long,aes(x=Age,y=CV,colour=Season))+geom_line()+theme_bw(base_size = 20)+xlab("Age")+ylab("CV Female Length")
ggplot(data=vissmsdMLvec_long,aes(x=Age,y=CV,colour=Season))+geom_line()+theme_bw(base_size = 20)+xlab("Age")+ylab("CV Male Length")

#### POPULATION MODEL ####
#calculate initial numbers at length (same for M&F)
popN=seasgrowpop(ts,fage,nages,Rmean,Rsd,Z,nseas,remainder)
popN=cbind(rep(seq(1,4,1),ts/nseas),popN)
fipopLF=as.matrix(genLF(nages,midbins,fmvec,fsdvec,popN,ts));colnames(fipopLF)=c(midbins)
mipopLF=as.matrix(genLF(nages,midbins,mmvec,msdvec,popN,ts));colnames(mipopLF)=c(midbins)
fnewpopN=fipopLF*(Ninit/2)
mnewpopN=fipopLF*(Ninit/2)

#NAL
fLFfin=fnewpopN/rowSums(fnewpopN)
mLFfin=mnewpopN/rowSums(mnewpopN)
fpopLFplot=data.frame(t(rbind(fLFfin,midbins,rep(1,length(midbins)))));colnames(fpopLFplot)=c(seq(1,ts),"Midbin","Model")#1=pop model
mpopLFplot=data.frame(t(rbind(mLFfin,midbins,rep(1,length(midbins)))));colnames(mpopLFplot)=c(seq(1,ts),"Midbin","Model")#1=pop model

#Visualize simulated population LFs
fpopLFplot$Bin=as.factor(fpopLFplot$Midbin)
mpopLFplot$Bin=as.factor(mpopLFplot$Midbin)
popLFplot=fpopLFplot
for(t in 1:24){
  allLFsub=popLFplot[,t]
  allLFsub=data.frame(cbind(allLFsub,popLFplot[,ts+1]));colnames(allLFsub)=c("Proportion","Midbin")
  g=ggplot(allLFsub,aes(x=Midbin,y=Proportion))+
    geom_line()+
    geom_point(size=4)+
    scale_shape_manual(values=c(16))+
    scale_color_manual(values=c("#000000"))+
    theme_bw(base_size = 20)+
    xlab("Length(midbin)")+
    theme(
      legend.position = c(0.15,.95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  print(g)
}

#### OBSERVATION MODEL ####
#### Set up observation model for LFA submodel ####
fobsLF=data.frame(samplepopLFA(ts,fLFfin,LFAsamps,fnewpopN))
mobsLF=data.frame(samplepopLFA(ts,mLFfin,LFAsamps,mnewpopN))
fsamps=round(LFAsamps*sum(fnewpopN,0))/ts
msamps=round(LFAsamps*sum(mnewpopN,0))/ts
colnames(fobsLF)=midbins
colnames(mobsLF)=midbins
fmLFobs=data.frame(rbind(fobsLF,mobsLF))

#visualize pop vs obs LFA
fobsLFplot=data.frame(t(rbind(fobsLF,midbins,rep(2,length(midbins)))));colnames(fobsLFplot)=c(seq(1,ts),"Midbin","Model")#2=obs model
mobsLFplot=data.frame(t(rbind(mobsLF,midbins,rep(2,length(midbins)))));colnames(mobsLFplot)=c(seq(1,ts),"Midbin","Model")#2=obs model
obsLFplot=fobsLFplot
LFApopLF=fLFfin
LFApopLF=LFApopLF/rowSums(LFApopLF)
LFApopLFplot=data.frame(t(rbind(LFApopLF,midbins,rep(1,length(midbins)))));colnames(LFApopLFplot)=c(seq(1,ts),"Midbin","Model")#1=pop model
allLF=rbind(LFApopLFplot,obsLFplot)
allLF$Bin=as.factor(allLF$Midbin);allLF$Model=as.factor(allLF$Model)
for(t in 1:ts){
  allLFsub=allLF[,t]
  allLFsub=data.frame(cbind(allLFsub,allLF[,ts+1],allLF[,ts+2]));colnames(allLFsub)=c("Proportion","Midbin","Model");allLFsub$Model=as.factor(allLFsub$Model);levels(allLFsub$Model) = c("Pop","Obs")
  g=ggplot(allLFsub,aes(x=Midbin,y=Proportion,color=Model))+
    geom_line()+
    geom_point(size=4,aes(shape=Model))+
    scale_shape_manual(values=c(16,15))+
    scale_color_manual(values=c("#000000", "#FF33B5"))+
    theme_bw(base_size = 20)+xlab("Length(midbin)")+
    theme(
      legend.position = c(0.2,.95),#c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  print(g)
}

##Tagging submodel
#Simulate multinomial dist sampling the population N@L 
fobstag=data.frame(samplepopTAG(nbins,ts,fLFfin,TAGsamps,fnewpopN))
mobstag=data.frame(samplepopTAG(nbins,ts,mLFfin,TAGsamps,mnewpopN))

#Expand out the sampled data to create tagging data sets (initial CL)
ftagdat=tagTAL(ts,fobstag,nbins,midbins)
mtagdat=tagTAL(ts,mobstag,nbins,midbins)

#Grow tagged individuals 
ftagsimdat=tagsimTVK(ftagdat,fLinf,fK,fKSd,fLinfSd,fRho,ts)
mtagsimdat=tagsimTVK(mtagdat,mLinf,mK,mKSd,mLinfSd,mRho,ts)
ftagsimdat$tRecap=ftagsimdat$tCap+ftagsimdat$DeltaT
mtagsimdat$tRecap=mtagsimdat$tCap+mtagsimdat$DeltaT
ftagsimdat=filter(ftagsimdat,tRecap<=24&DeltaT>0)
mtagsimdat=filter(mtagsimdat,tRecap<=24&DeltaT>0)

#Generate tagging dataset for analysis 
nrecs=as.vector(c(nrow(ftagsimdat),nrow(mtagsimdat)))
simtagdat=rbind(ftagsimdat[,c(1,8,2,4,5)],mtagsimdat[,c(1,8,2,4,5)])

############## SEND TO ADMB ####################
#est model settings
guessnages=6

#create .dat file
shell("del compgV1.dat")
write(c(nseas,
        nsex,
        nage=guessnages,
        ts,
        fage=fa,
        nbins,
        minbin=fbin,
        binw=binwidth,
        midbins,
        GDD,
        fEffN, 
        mEffN,
        flfdata_mlfdata=t(as.matrix(fmLFobs)),
        TaggingDim=nrecs,
        TaggingData=t(as.matrix(simtagdat)),
        999),
      file = "compgV1.dat", sep = " ")

#parameter vectors - 6 cols: phase, logscale lower_boundary, logscale upper_boundary, type_of_distribution, mean, std_deviation, arith scale starting_value
fkpvec=as.vector(c(1, -9.9, -9.7, 0, 0.0, 0.0, fkp)) #0.00005555
mkpvec=as.vector(c(1, -9.9, -9.7, 0, 0.0, 0.0, mkp))
fLinfvec=as.vector(c(-1, 5.0, 6.0, 0, 0.0, 0.0, fLinf))
mLinfvec=as.vector(c(-1, 5.0, 6.0, 0, 0.0, 0.0, mLinf))
fkSevec=as.vector(c(-1, -15, 0, 0, 0, 0.0, fKSd))
mkSevec=as.vector(c(-1, -15, 0, 0, 0, 0.0, mKSd))
fLinfSevec=as.vector(c(-1, -15, 10, 0, 0, 0.0, fLinfSd))
mLinfSevec=as.vector(c(-1, -15, 10, 0, 0, 0.0, mLinfSd))
fRhovec=as.vector(c(-1, -1, 1, 0, 0, 0.0, fRho)) 
mRhovec=as.vector(c(-1, -1, 1, 0, 0, 0.0, mRho))

fsdmat=matrix(rep(t(fsdvec),ts/nseas),ncol=ncol(fsdvec),byrow=TRUE)
msdmat=matrix(rep(t(msdvec),ts/nseas),ncol=ncol(msdvec),byrow=TRUE)

#Use to modify mean length at age and sd of length at age if guessed nages is different from true
# fsdmat=data.frame(fsdmat)
# fsdmat=cbind(fsdmat,fsdmat[,6])
# msdmat=data.frame(msdmat)
# msdmat=cbind(msdmat,msdmat[,6])
#fsdmat=data.frame(fsdmat)
#fsdmat=fsdmat[,1:5]
#msdmat=data.frame(msdmat)
#msdmat=msdmat[,1:5]

#create .ctl file 
shell("del compgV1.ctl")
write(c(fkpvec,
        mkpvec,
        fLinfvec,
        mLinfvec,
        lfirstphase=-1,
        iflfirst=fmvec[1,1],
        imlfirst=mmvec[1,1],
        t(as.matrix(fsdmat)),
        t(as.matrix(msdmat)),
        fLinfSevec,
        mLinfSevec,
        fkSevec,
        mkSevec,
        fRhovec,
        mRhovec,
        999),
      file = "compgV1.ctl", sep = " ")

shell("del compgV1.std")
shell("del LFA_est.rep")
shell("del LFA_obs.rep")
shell("del compgV1.rep")
shell("del GI_est.rep")
shell("del GI_obs.rep")
#shell("compgV1.exe -nohess")
shell("compgV1.exe")

############## VIEW Tagging model RESULTS ####################
#Grab and plot obs vs est tagging DeltaLs
estdeltaLs=data.frame((read.delim("GI_est.rep",sep="",header=FALSE)))
festdeltaLs=t(estdeltaLs[1,]);festdeltaLs=festdeltaLs[!is.na(festdeltaLs)]
mestdeltaLs=t(estdeltaLs[2,]);mestdeltaLs=mestdeltaLs[!is.na(mestdeltaLs)]
obsdeltaLs=data.frame((read.delim("GI_obs.rep",sep="",header=FALSE)))
fobsdeltaLs=t(obsdeltaLs[1,]);fobsdeltaLs=fobsdeltaLs[!is.na(fobsdeltaLs)]
mobsdeltaLs=t(obsdeltaLs[2,]);mobsdeltaLs=mobsdeltaLs[!is.na(mobsdeltaLs)]

fobsDL=data.frame(cbind(as.vector(fobsdeltaLs),rep(2,length(fobsdeltaLs))));colnames(fobsDL)=c("DeltaL","Results")
mobsDL=data.frame(cbind(as.vector(mobsdeltaLs),rep(2,length(mobsdeltaLs))));colnames(mobsDL)=c("DeltaL","Results")
estfDL=data.frame(cbind(as.vector(festdeltaLs),rep(1,length(festdeltaLs))));colnames(estfDL)=c("DeltaL","Results")
estmDL=data.frame(cbind(as.vector(mestdeltaLs),rep(1,length(mestdeltaLs))));colnames(estmDL)=c("DeltaL","Results")

fDLres=data.frame(rbind(fobsDL,estfDL))
mDLres=data.frame(rbind(mobsDL,estmDL))
fDLres$Results=as.factor(fDLres$Results);levels(fDLres$Results) = c("Est","Obs")
mDLres$Results=as.factor(mDLres$Results);levels(mDLres$Results) = c("Est","Obs")

#Obs vs est DeltaLs
test=subset(fDLres,fDLres$Results=="Est")
test2=subset(fDLres,fDLres$Results=="Obs")
plot(test2$DeltaL,test$DeltaL,xlab="Obs",ylab="Est")
abline(1,1)
testy=subset(mDLres,mDLres$Results=="Est")
testy2=subset(mDLres,mDLres$Results=="Obs")
plot(testy2$DeltaL,testy$DeltaL,xlab="Obs",ylab="Est")
abline(1,1)

ggplot(fDLres,aes(x=DeltaL ,fill = Results))+geom_histogram(color="#D5DBDB", position = "identity", alpha = 0.4,binwidth = 0.5)+xlab("DeltaL")+ylab("Count")+theme(legend.position="right")+scale_fill_manual(values=c("#FF33B5","#660882"))+theme_bw()+theme(text = element_text(size = 20))+
  theme(
    legend.position =  c(0.3,.95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

ggplot(mDLres,aes(x=DeltaL ,fill = Results))+geom_histogram(color="#D5DBDB", position = "identity", alpha = 0.4,binwidth = 0.5)+xlab("DeltaL")+ylab("Count")+theme(legend.position="right")+scale_fill_manual(values=c("#FF33B5","#660882"))+theme_bw()+theme(text = element_text(size = 20))+
  theme(
    legend.position =  c(0.3,.95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

############## VIEW LFA RESULTS ####################
estLFs=data.frame(read.delim("LFA_est.rep",sep="",header=FALSE))
fMFLFsplot=t(rbind(estLFs[1:ts,],midbins,rep(3,length(midbins))));colnames(fMFLFsplot)=c(seq(1,ts),"Midbin","Model")
mMFLFsplot=t(rbind(estLFs[(ts+1):nrow(estLFs),],midbins,rep(3,length(midbins))));colnames(mMFLFsplot)=c(seq(1,ts),"Midbin","Model")
obsLFs=data.frame(read.delim("LFA_obs.rep",sep="",header=FALSE))
nfobsLF=t(rbind(obsLFs[1:ts,],midbins,rep(2,length(midbins))));colnames(nfobsLF)=c(seq(1,ts),"Midbin","Model")
nmobsLF=t(rbind(obsLFs[(ts+1):nrow(obsLFs),],midbins,rep(2,length(midbins))));colnames(nmobsLF)=c(seq(1,ts),"Midbin","Model")
fresLF=data.frame(rbind(fMFLFsplot,nfobsLF))
mresLF=data.frame(rbind(mMFLFsplot,nmobsLF))
fresLF$Model=as.factor(fresLF$Model);levels(fresLF$Model)= c("Obs","Est")
mresLF$Model=as.factor(mresLF$Model);levels(mresLF$Model)= c("Obs","Est")

#Plot tvK
estks=as.matrix(read.delim("tvK.rep",sep="",header=FALSE))
festks=estks[1,];mestks=estks[2,]
allfks=cbind(fK,as.numeric(estks[1,]),seq(1:ts))
plot(allfks[,3],allfks[,1],col="red",pch=20,ylab="Simulated k",xlab="Time step")
points(allfks[,2],col="black",pch=1)

#Obs vs est LFs all time steps
resLF=fresLF
for(t in 1:ts){
  resLFsub=resLF[,t]
  resLFsub=data.frame(cbind(resLFsub,resLF[,ts+1],resLF[,ts+2]));colnames(resLFsub)=c("Proportion","Midbin","Model");resLFsub$Model=as.factor(resLFsub$Model);levels(resLFsub$Model) = c("Obs","Est")
  g=ggplot(resLFsub,aes(x=Midbin,y=Proportion,color=Model))+
    geom_line()+geom_point(size=4,aes(shape=Model))+
    scale_shape_manual(values=c(16,15))+
    scale_color_manual(values=c("#FF33B5","#660882"))+
    theme_bw(base_size = 20)+
    xlab("Length(midbin)")+
    theme(
      legend.position = c(.95, .95),# c(0.2,.95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  print(g)
}