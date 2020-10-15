
rm(list=ls())

Rcpp::sourceCpp('/home/rlin/AMS-design/AMS-cal2.cpp')
get.oc <- function(target, pE.true,pT.true,  ncohort, cohortsize, startdose=1, cutoff.eli=0.95, 
                   lfT, lfE, lfC, pEE.sample, ntrial=10)
{ 
  set.seed(6);
  ndose=length(pE.true)	
  npts = ncohort*cohortsize;
  YT=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  YE=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
  dselect = rep(0, ntrial); # store the selected dose level
  sel=rep(0,ndose);
  pts=rep(0,ndose);
  dlt=rep(0,ndose);
  eff=rep(0,ndose);
  tox=0;
  acr=0;
  exc=0;
  
  ################## simulate trials ###################
  for(trial in 1:ntrial)
  {
    ptm <- proc.time()
    yT<-yE<-rep(0, ndose);    ## the number of DLT at each dose level
    n<-rep(0, ndose);    ## the number of patients treated at each dose level
    earlystop=0;         ## indiate whether the trial terminates early
    d=startdose;         ## starting dose level
    elimi = rep(0, ndose);  ## indicate whether doses are eliminated
    likeliT<-matrix(rep(1,NN*(ndose+1)),ncol=ndose+1)
    likeliE<-matrix(rep(1,NN*ndose),ncol=ndose)
    likeliC<-matrix(rep(1,NN*ndose),ncol=ndose)
    likT<-rep(1,ndose+1)
    likE<-rep(1,ndose)
    UT<-rep(0,ndose)
    UE<-rep(0,ndose) 
    for(i in 1:ncohort)  
    { 
      ### generate toxicity outcome
      wT = sum(runif(cohortsize)<pT.true[d])
      yT[d] = yT[d] + wT;
      wE = sum(runif(cohortsize)<pE.true[d])
      yE[d] = yE[d] + wE;
      n[d] = n[d] + cohortsize;
      nc = n[d]/cohortsize;
      temp = post_cal2(likeliT,likeliE,likeliC,lfT[,(d-1)*(cohortsize+1)+wT+1,],
                      lfE[,(d-1)*(cohortsize+1)+wE+1,],lfC[,(d-1)*(cohortsize+1)+wE+1,],ndose);
      likeliT = temp$likeliT;
      likeliE = temp$likeliE;
      likeliC = temp$likeliC;
      posT = temp$posT;
      if(n[d]>=12 & max(yE[d]/(n[d]))>=0.25){
        posE = temp$posC;
      } else {posE = temp$posE;}
      # likeliT<-likeliT*lfT[,(d-1)*(cohortsize+1)+wT+1,]
      # likeliE<-likeliE*lfE[,(d-1)*(cohortsize+1)+wE+1,]
      # likT<-apply(likeliT,2,mean)
      # likE<-apply(likeliE,2,mean)
      # posT<-likT/sum(likT)
      # posE<-likE/sum(likE)
      cumposT<-cumsum(posT)
      dT_opt<-max(min(which(cumposT>0.85))-1,1)
      #dT_opt<-max(which.max(posT)-1,1)
      dE_opt<-which.max(posE)
      d_opt<-min(dE_opt,dT_opt)
      
      if ((1-pbeta(target,1+yT[d],n[d]-yT[d]+1))>0.95) {elimi[d:ndose]<-1
      if (elimi[1]==1) {earlystop=1; break;}}
      if (pbeta(0.25,1+yE[d],1+n[d]-yE[d])>0.9){elimi[d]<-1}
      if (sum(elimi==1)==ndose){earlystop=1;break;}
      
      if (sum(yT+yE)==0 & d<ndose){
        d<-d+1
      } else{
        
        if (d_opt>d & sum(elimi[(d+1):ndose]==0)>0){
          d<-d+which(elimi[(d+1):ndose]==0)[1]
        } else if (d_opt<d & sum(elimi[1:(d-1)]==0)>0){
          d<-max(which(elimi[1:(d-1)]==0))
        }
        
      }            
      if (sum(elimi==1)==ndose){earlystop=1;break;}
      if(elimi[d]==1){earlystop=1;break}           
      
    }
    
    YT[trial,]=yT;
    YE[trial,]=yE;
    N[trial,]=n;
    if (earlystop==0){
      posE<-posteriorHE(yE,n,pEE.sample) 
      dE_opt<-which.max(posE)
      dE_opt<-min(which(posE>=(max(posE)-0.1)))
      dE_opt<-min(which(posE>=(max(posE)-0.05)))		
      d_opt<-min(dE_opt,dT_opt)  		
      dselect[trial]=d_opt
      sel[dselect[trial]]=sel[dselect[trial]]+1/ntrial*100
    } else {dselect[trial]<-99}	
    pts<-pts+n/ntrial
    dlt<-dlt+yT/ntrial
    eff<-eff+yE/ntrial
    #print(c(trial,proc.time()-ptm))
  }	
  pts<-round(pts,1)
  dlt<-round(dlt,1)
  #print(sel)
  #print(pts) 
  results=list(sel=sel,pts=pts,dlt=dlt,eff=eff)
  return(results)	
}



target = 0.3
cohortsize = 3
ncohort = 10

ndose = 5;
NN<-50000
pT.sample<-matrix(0,NN,ndose*(ndose+1));

pT.sample[,1]<-runif(NN,target,1)
for (i in 2:ndose){
  pT.sample[,i]<-runif(NN,pT.sample[,i-1],1)
}

for (i in 1:ndose){
  pT.sample[,i+i*ndose]<-runif(NN,0,target)
  if (i>1){
    for (j in (i-1):1){
      pT.sample[,j+i*ndose]<-runif(NN,0,pT.sample[,j+1+i*ndose])
    }
  }
  if (i<ndose){
    for (j in (i+1):ndose){
      pT.sample[,j+i*ndose]<-runif(NN,sapply(pT.sample[,j-1+i*ndose],function(x) max(x,target)),1)
    }
  }
  
}


pE.sample<-matrix(0,NN,ndose*ndose);

for (i in 1:ndose){
  pE.sample[,i+(i-1)*ndose]<-runif(NN,0.25,1)
  if (i>1){
    for (j in (i-1):1){
      pE.sample[,j+(i-1)*ndose]<-runif(NN,0,pE.sample[,j+1+(i-1)*ndose])
    }
    #pE.sample[,1:i+(i-1)*ndose]<-t(apply(pE.sample[,1:i+(i-1)*ndose],1,sort))
  }
  if (i<ndose){
    for (j in (i+1):ndose){
      pE.sample[,j+(i-1)*ndose]<-runif(NN,0,pE.sample[,j-1+(i-1)*ndose])
    }
    #pE.sample[,i:ndose+(i-1)*ndose]<-t(apply(pE.sample[,i:ndose+(i-1)*ndose],1,sort,decreasing=T))
  }
  
}


pC.sample<-matrix(0,NN,ndose*ndose);

for (i in 1:ndose){
  pC.sample[,i+(i-1)*ndose]<-runif(NN,0.5,1)
  if (i>1){
    for (j in (i-1):1){
      pC.sample[,j+(i-1)*ndose]<-runif(NN,0,pC.sample[,j+1+(i-1)*ndose])
    }
    #pE.sample[,1:i+(i-1)*ndose]<-t(apply(pE.sample[,1:i+(i-1)*ndose],1,sort))
  }
  if (i<ndose){
    for (j in (i+1):ndose){
      pC.sample[,j+(i-1)*ndose]<-runif(NN,0,pC.sample[,j-1+(i-1)*ndose])
    }
    #pE.sample[,i:ndose+(i-1)*ndose]<-t(apply(pE.sample[,i:ndose+(i-1)*ndose],1,sort,decreasing=T))
  }
  
}

pEE.sample<-matrix(0,NN,ndose*ndose);

for (i in 1:ndose){
  pEE.sample[,i+(i-1)*ndose]<-runif(NN,0,1)
  if (i>1){
	for(j in (i-1):1){
	  temp<-rbinom(NN,1,0.7)
      pEE.sample[,j+(i-1)*ndose]<-temp*runif(NN,0,pEE.sample[,j+1+(i-1)*ndose])+(1-temp)*pEE.sample[,j+1+(i-1)*ndose]
	  }
  }
  if (i<ndose){
    for (j in (i+1):ndose){
	  temp<-rbinom(NN,1,0.7)
      pEE.sample[,j+(i-1)*ndose]<-temp*runif(NN,0,pEE.sample[,j-1+(i-1)*ndose])+(1-temp)*pEE.sample[,j-1+(i-1)*ndose]
    }
    #pEE.sample[,i:ndose+(i-1)*ndose]<-t(apply(pEE.sample[,i:ndose+(i-1)*ndose],1,sort,decreasing=T))
  }
  
}

lfT<-array(dim=c(NN,ndose*(cohortsize+1),ndose+1))

for (i in 1:(ndose+1)){
  for (j in 1:ndose){
    for (k in 0:cohortsize){
      lfT[,(j-1)*(cohortsize+1)+k+1,i]<-pT.sample[,j+(i-1)*ndose]^k*(1-pT.sample[,j+(i-1)*ndose])^(cohortsize-k)
    }
  }
}

lfE<-array(dim=c(NN,ndose*(cohortsize+1),ndose))

for (i in 1:ndose){
  for (j in 1:ndose){
    for (k in 0:cohortsize){
      lfE[,(j-1)*(cohortsize+1)+k+1,i]<-pE.sample[,j+(i-1)*ndose]^k*(1-pE.sample[,j+(i-1)*ndose])^(cohortsize-k)
    }
  }
}

lfC<-array(dim=c(NN,ndose*(cohortsize+1),ndose))

for (i in 1:ndose){
  for (j in 1:ndose){
    for (k in 0:cohortsize){
      lfC[,(j-1)*(cohortsize+1)+k+1,i]<-pC.sample[,j+(i-1)*ndose]^k*(1-pC.sample[,j+(i-1)*ndose])^(cohortsize-k)
    }
  }
}


pE.true<-c(0.10,0.27,0.44,0.58,0.69)
pT.true<-c(0.04,0.18,0.37,0.54,0.67)
#pE.true<-c(0.10,0.20,0.25,0.50,0.54)
#pT.true<-c(0.05,0.07,0.10,0.15,0.35)
pE.true<-c(0.6,0.8,0.5,0.4,0.2)
pT.true<-c(0.01,0.05,0.10,0.15,0.30)
#pE.true<-c(0.05,0.08,0.15,0.28,0.43)
#pT.true<-c(0.02,0.05,0.07,0.10,0.12)
#pE.true<-c(0.2,0.3,0.6,0.8,0.55)
#pT.true<-c(0.08,0.12,0.2,0.30,0.4)
#pE.true<-c(0.2,0.4,0.6,0.8,0.55)
#pT.true<-c(0.06,0.08,0.14,0.2,0.3)
#pE.true<-c(0.2,0.4,0.6,0.8,0.55)
#pT.true<-c(0.05,0.1,0.25,0.5,0.6)
#pE.true<-c(0.05,0.25,0.45,0.68,0.80)
#pT.true<-c(0.05,0.10,0.15,0.20,0.50)
#pE.true<-c(0.1,0.3,0.5,0.5,0.5)
#pT.true<-c(0.1,0.2,0.4,0.5,0.6)


p<-matrix(ncol=5,nrow=12)
p[1,]<-c(0.08,0.12,0.20,0.30,0.40)
p[2,]<-c(0.01,0.05,0.10,0.15,0.30)
p[3,]<-c(0.06,0.08,0.14,0.20,0.30)
p[4,]<-c(0.05,0.10,0.25,0.50,0.60)
p[5,]<-c(0.05,0.10,0.15,0.20,0.50)
p[6,]<-c(0.10,0.20,0.40,0.50,0.60)
p[7,]<-c(0.15,0.32,0.45,0.55,0.62)
p[8,]<-c(0.04,0.18,0.37,0.54,0.67)
p[9,]<-c(0.02,0.05,0.07,0.10,0.12)
p[10,]<-c(0.10,0.12,0.15,0.36,0.65)
p[11,]<-c(0.05,0.07,0.10,0.15,0.35)
p[12,]<-c(0.10,0.12,0.15,0.30,0.60)

q<-matrix(ncol=5,nrow=12)
q[1,]<-c(0.20,0.40,0.60,0.80,0.55)
q[2,]<-c(0.60,0.80,0.50,0.40,0.20)
q[3,]<-c(0.20,0.40,0.60,0.80,0.55)
q[4,]<-c(0.20,0.40,0.60,0.80,0.55)
q[5,]<-c(0.05,0.25,0.45,0.68,0.80)
q[6,]<-c(0.10,0.30,0.50,0.50,0.50)
q[7,]<-c(0.28,0.30,0.44,0.60,0.74)
q[8,]<-c(0.10,0.27,0.44,0.58,0.69)
q[9,]<-c(0.05,0.08,0.15,0.28,0.43)
q[10,]<-c(0.15,0.18,0.38,0.40,0.60)
q[11,]<-c(0.10,0.20,0.25,0.50,0.54)
q[12,]<-c(0.02,0.10,0.42,0.45,0.50)




p<-matrix(ncol=5,nrow=20)
p[1,]<-c(0.10,0.35,0.40,0.45,0.50)
p[2,]<-c(0.01,0.05,0.10,0.15,0.30)
p[3,]<-c(0.10,0.25,0.38,0.50,0.60)
p[4,]<-c(0.05,0.10,0.25,0.50,0.60)
p[5,]<-c(0.10,0.12,0.15,0.25,0.35)
p[6,]<-c(0.05,0.10,0.15,0.20,0.35)
p[7,]<-c(0.02,0.05,0.07,0.10,0.15)
p[8,]<-c(0.05,0.10,0.15,0.20,0.25)
p[9,]<-c(0.01,0.02,0.05,0.25,0.35)
p[10,]<-c(0.01,0.02,0.05,0.10,0.25)
p[11,]<-c(0.01,0.02,0.15,0.38,0.45)
p[12,]<-c(0.10,0.12,0.15,0.25,0.35)
p[13,]<-c(0.01,0.02,0.03,0.04,0.05)
p[14,]<-c(0.20,0.30,0.40,0.50,0.60)
p[15,]<-c(0.05,0.17,0.25,0.38,0.52)
p[16,]<-c(0.05,0.10,0.15,0.20,0.25)
p[17,]<-c(0.08,0.10,0.15,0.32,0.40)
p[18,]<-c(0.10,0.23,0.36,0.44,0.51)
p[19,]<-c(0.50,0.55,0.60,0.65,0.70)
p[20,]<-c(0.25,0.55,0.67,0.87,0.95)

q<-matrix(ncol=5,nrow=20)
q[1,]<-c(0.30,0.40,0.55,0.65,0.75)
q[2,]<-c(0.25,0.40,0.40,0.40,0.40)
q[3,]<-c(0.15,0.35,0.40,0.45,0.50)
q[4,]<-c(0.10,0.25,0.45,0.60,0.50)
q[5,]<-c(0.15,0.18,0.38,0.38,0.38)
q[6,]<-c(0.10,0.20,0.30,0.40,0.50)
q[7,]<-c(0.05,0.08,0.15,0.30,0.45)
q[8,]<-c(0.30,0.45,0.40,0.35,0.30)
q[9,]<-c(0.25,0.45,0.65,0.65,0.65)
q[10,]<-c(0.10,0.20,0.35,0.55,0.55)
q[11,]<-c(0.25,0.45,0.45,0.45,0.45)
q[12,]<-c(0.15,0.18,0.50,0.70,0.70)
q[13,]<-c(0.50,0.40,0.30,0.20,0.10)
q[14,]<-c(0.25,0.40,0.55,0.70,0.85)
q[15,]<-c(0.20,0.30,0.45,0.30,0.20)
q[16,]<-c(0.15,0.30,0.45,0.40,0.35)
q[17,]<-c(0.10,0.20,0.55,0.55,0.65)
q[18,]<-c(0.40,0.60,0.70,0.75,0.80)
q[19,]<-c(0.30,0.35,0.45,0.66,0.70)
q[20,]<-c(0.01,0.55,0.55,0.55,0.55)

# 
# p<-matrix(ncol=5,nrow=16)
# p[1,]<-c(0.05,0.12,0.27,0.35,0.50)
# p[2,]<-c(0.05,0.12,0.27,0.35,0.50)
# p[3,]<-c(0.03,0.05,0.20,0.22,0.45)
# p[4,]<-c(0.05,0.10,0.15,0.20,0.27)
# p[5,]<-c(0.03,0.06,0.10,0.30,0.45)
# p[6,]<-c(0.02,0.05,0.10,0.20,0.30)
# p[7,]<-c(0.01,0.15,0.20,0.40,0.50)
# p[8,]<-c(0.15,0.20,0.25,0.35,0.45)
# p[9,]<-c(0.03,0.05,0.07,0.09,0.11)
# p[10,]<-c(0.01,0.08,0.10,0.20,0.25)
# p[11,]<-c(0.08,0.10,0.15,0.32,0.40)
# p[12,]<-c(0.05,0.07,0.10,0.20,0.35)
# p[13,]<-c(0.01,0.05,0.10,0.12,0.27)
# p[14,]<-c(0.02,0.05,0.07,0.10,0.15)
# p[15,]<-c(0.50,0.55,0.60,0.65,0.70)
# p[16,]<-c(0.25,0.55,0.67,0.87,0.95)
# 
# 
# 
# 
# 
# 
# q<-matrix(ncol=5,nrow=16)
# q[1,]<-c(0.20,0.35,0.36,0.37,0.38)
# q[2,]<-c(0.01,0.05,0.30,0.60,0.60)
# q[3,]<-c(0.05,0.10,0.50,0.68,0.70)
# q[4,]<-c(0.01,0.05,0.10,0.20,0.40)
# q[5,]<-c(0.10,0.20,0.40,0.45,0.50)
# q[6,]<-c(0.05,0.15,0.40,0.40,0.40)
# q[7,]<-c(0.45,0.45,0.45,0.65,0.80)
# q[8,]<-c(0.25,0.55,0.40,0.30,0.20)
# q[9,]<-c(0.45,0.30,0.25,0.20,0.10)
# q[10,]<-c(0.07,0.10,0.35,0.20,0.20)
# q[11,]<-c(0.10,0.20,0.70,0.70,0.75)
# q[12,]<-c(0.05,0.15,0.30,0.55,0.55)
# q[13,]<-c(0.10,0.15,0.25,0.50,0.50)
# q[14,]<-c(0.05,0.08,0.15,0.30,0.45)
# q[15,]<-c(0.30,0.35,0.45,0.66,0.70)
# q[16,]<-c(0.01,0.55,0.55,0.55,0.55)
# 



res<-NULL;
for (i in 1:20){
  pT.true<-p[i,]
  pE.true<-q[i,]
  oc=get.oc(target, pE.true,pT.true,  ncohort, cohortsize, startdose=1, lfT=lfT, lfE=lfE, lfC=lfC, pEE.sample=pEE.sample, cutoff.eli=0.95, ntrial=2000)
  print(i)
  print(pT.true)
  print(pE.true)
  print(oc)
  res<-rbind(res,oc)
}

