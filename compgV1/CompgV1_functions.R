#Functions for compgV1 growth simulation model

#Generate a growth transition matrix
growth_matrix<-function(Linf,K,Linf_Var,K_Var,Linf_K_Cov,bin_mids)
{
  nb<-length(bin_mids)
  bw<-bin_mids[2]-bin_mids[1] #calculate width of bins
  mat<-matrix(0,nrow=nb, ncol=nb) #set up empty matrix
  mean_g<-(Linf-bin_mids)*(1-exp(-K)) #calculate mean growth as a function of starting length
  dmidsize<-Linf-bin_mids
  
  #calculate variance for starting length bin i
  L_SD<-Linf_Var*(1-exp(-K))^2
  L_SD<-L_SD+K_Var*(dmidsize*exp(-K))^2
  L_SD<-L_SD+2.0*Linf_K_Cov*Linf_Var*K_Var*(1-exp(-K))*dmidsize*exp(-K);
  L_SD<-sqrt(L_SD) #convert from var to SD
  
  for(i in 1:(nb-1))  #loop over starting length bins
  {
    up_lim=bw/2
    mat[i,i]=pnorm(up_lim,mean_g[i],sd=L_SD[i])
    
    for(j in (i+1):(nb-1)) #loop over ending length bins
    {
      if((i+1)>(nb-1)) {
        break
      }
      #calculate increment
      inc=bin_mids[j]-bin_mids[i]
      #calculate lower bound of growth increment
      lo_lim=up_lim
      #calculate upper bound of growth increment
      up_lim=lo_lim+bw
      
      #Calculate proportion that grow to bin j
      mat[j,i]=pnorm(up_lim,mean_g[i],sd=L_SD[i])-pnorm(lo_lim,mean_g[i],sd=L_SD[i])
      #print(c(i,j))
    }#close j
    #last length bin a plus group
    lo_lim=up_lim
    mat[nb,i]=1-pnorm(lo_lim,mean_g[i],sd=L_SD[i]) 
    #print(i)
    #normalize proportions so they sum to one
     #sp<-sum(mat[,i]) 
     #if(sp>0)
     #{
    #   mat[,i]=mat[,i]/sp
     #}#close if
    #print(i)
  }#close i
  mat[nb,nb]=1
  return(mat)
}

#growth increment sd
incsd=function(Linf,K,bins,KSe,LinfSe,Rho)
{
  incsd=as.vector(bins)
  for(i in 1:length(bins)){
    incsd[i]=LinfSe^2*(1-exp(-K))^2+(Linf-bins[i])^2*KSe^2*exp(-2*(K))+2*Rho*LinfSe*KSe*(1-exp(-(K)))*(Linf-bins[i])*exp(-K)
  }
  return(incsd)
}

#seasonal population dynamics
seasgrowpop=function(ts,fage,nages,Rmean,Rsd,Z,nseas,remainder)
  {
    N=matrix(0,ts,nages)
    #ts=1
    N[1,1]=rnorm(1,mean=Rmean,sd=Rsd)
    for(a in 2:nages)
      {
        N[1,a]=N[1,a-1]*exp(-Z[a])
    }
    #ts=2+
    for(t in 2:ts){
      #if seas=1
      if(t%%nseas==remainder){
        N[t,1]=rnorm(1,mean=Rmean,sd=Rsd)
        for(a in 2:nages){
          #if(a<nages){
            N[t,a]=N[t-1,a-1]*exp(-Z[a])
          #}
          #else{N[t,a]=N[t-1,a-1]*exp(-Z[a])+N[t-1,a]*exp(-Z[a])}
        }
      } else{
        for(a in 1:nages){
            N[t,a]=N[t-1,a]*exp(-Z[a])
          }
        }
      }
    N=round(N,0)
    return(N)
  }

#determined population LF
genLF=function(nages,midbins,mvec,sdvec,N,ts)
  {
    popLF=matrix(0,ts,length(midbins))
    for(t in 1:ts){
        s=N[t,1]
        LF=matrix(0,nages,length(midbins))
        for(a in 1:nages){
          LF[a,]=dnorm(midbins,mvec[s,a],sdvec[s,a])*N[t,a+1]
          }
        popLF[t,]=colSums(LF)/sum(colSums(LF))
      }
    return(popLF)
}

#LFA of observed samples
samplepopLFA=function(ts,popLF,LFAsamps,N)
  {
    #sampLs=matrix(0,ts,highbin-lowbin+1) 
    sampLs=matrix(0,ts,ncol(popLF)) 
    for (t in 1:ts){
      nsamps=round(LFAsamps*sum(N[t,]),0)
      #nal=rmultinom(1,nsamps,prob=c(popLF[t,lowbin:highbin]))
      nal=rmultinom(1,nsamps,prob=c(popLF[t,]))
      sampLs[t,1:ncol(popLF)]=nal/nsamps
  }
  return(sampLs)
  }

#number of lobsters tagged at length
#samplepopTAG=function(tfa,tla,tlowbin,thighbin,ts,popLF,TAGsamps,N)
samplepopTAG=function(nbins,ts,popLF,TAGsamps,N)
  {
    sampLs=matrix(0,ts,nbins) 
    for (t in 1:ts){
      #nsamps=round(TAGsamps*sum(N[t,(tfa+1):(1+tla)]),0)
      nsamps=round(TAGsamps*sum(N[t,1:nbins]),0)
      #nal=rmultinom(1,nsamps,prob=c(popLF[t,tlowbin:thighbin]))
      nal=rmultinom(1,nsamps,prob=c(popLF[t,1:nbins]))
      sampLs[t,]=nal
  }
  return(sampLs)
}

#Expand simulated data to begin creating a tagging dataset
tagTAL=function(ts,obstag,nbins,midbins)
{
  obstag$seas=rep(seq(1:4),ts/nseas)
  tmp=matrix(0,1,7)
  for (t in 1:ts){
    for (i in 1:nbins){
      if(obstag[t,i]>0) {tmp2=cbind(t,rep(midbins[i],obstag[t,i]),obstag[t,nbins+1],0,0,0,0)}
      else {tmp2=cbind(0,0,0,0,0,0,0)}
      tmp=rbind(tmp,tmp2)
    }
  }
  tmp=tmp[-1,]
  tmp=data.frame(subset(tmp,tmp[,1]>0))
  colnames(tmp)=c("tCap","CL","NewSeas","DeltaT","DeltaL","sd","CL2")
  return(tmp)
}

#Simulate growth of tagged animals using GImodel and GIsd and using TVK
tagsimTVK=function(tagdat,Linf,K,KSd,LinfSd,Rho,ts)
{
  tagdat[,4]=round(rlnorm(nrow(tagdat),mean=1,sd=.5),0) 
  #hist(tagdat[,4]);min(tagdat[,4]);max(tagdat[,4])
  for (i in 1:nrow(tagdat)){
    if(tagdat[i,4]>ts){tagdat[i,4]=ts}
  }
  #Calculate GI and sd of GI using random unit normal dist with mean of 0 and GIsd
  tagdat[,5]=(Linf-tagdat[,2])*(1-exp((-(K[tagdat[i,3]])*tagdat[,4])))  
  #tagdat[,5]=incsd(Linf,fK,initCL,KSd,LinfSd,Rho)
  tagdat[,6]=rnorm(nrow(tagdat),0,incsd(Linf,K[tagdat[i,3]],tagdat[,2],KSd,LinfSd,Rho))
  #Calculate CL2 while preventing significant shrinkage
  for (i in 1:nrow(tagdat)){
    if(tagdat[i,5]+tagdat[i,6]>=0){tagdat[i,7]=tagdat[i,2]+tagdat[i,5]+tagdat[i,6]} else {tagdat[i,7]=tagdat[i,2]}
  }
  return(tagdat)
}