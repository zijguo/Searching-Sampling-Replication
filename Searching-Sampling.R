########## Start date: Jan 2, 2021; Most updated date: March 30, 2021
########## Author: Zijian Guo
########## Email: zijguo@stat.rutgers.edu
########## Goal: uniform inference for the treatment effect under majority and plurality rule
########## Method: searching and sampling
###### bootstrap threshold
norm.diff<-function(vec,beta.grid,SigmaSqD,SigmaSqY,SigmaYD,InitiSet,SE.norm){
  gamma.s<-vec[1:pz]
  Gamma.s<-vec[-(1:pz)]
  norm.max<-0
  n.beta<-length(beta.grid)
  for(j in 1:n.beta){
    b<-beta.grid[j]
    se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
    temp.diff<-abs(Gamma.s[InitiSet]-b*gamma.s[InitiSet])/(se.b*SE.norm[InitiSet])
    norm.max<-max(norm.max,max(temp.diff))
  }
  return(norm.max)
}
cut.off<-function(SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,pz,beta.grid,N=1000){
  #N=1000   
  n=nrow(WUMat)
  unit.matrix<-t(WUMat)%*%WUMat/n^2
  #(solve(covW)/n)[1:pz,1:pz]
  Cov1<-cbind(SigmaSqD*unit.matrix,SigmaYD*unit.matrix)
  Cov2<-cbind(SigmaYD*unit.matrix,SigmaSqY*unit.matrix)
  Cov.total<-rbind(Cov1,Cov2)
  Gen.mat<-mvrnorm(N, rep(0,2*pz), Cov.total)
  #se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
  SE.norm<-diag(unit.matrix)^{1/2}
  sample.sim<-rep(0,N)
  for(j in 1:N){
    sample.sim[j]<-norm.diff(Gen.mat[j,],beta.grid,SigmaSqD,SigmaSqY,SigmaYD,InitiSet,SE.norm)
  }
  critical.val<-quantile(sample.sim,probs=0.95)
  return(critical.val)
}

cut.off.IVStr<-function(SigmaSqD,WUMat,pz,N=1000,cut.prob=0.99){
   n=nrow(WUMat)
   unit.matrix<-t(WUMat)%*%WUMat/n^2
   Cov.D<-SigmaSqD*unit.matrix
   max.vec<-rep(NA,N)
   for(j in 1:N){
     gamma.s<-mvrnorm(1, rep(0,pz), Cov.D)
     SE.norm<-diag(unit.matrix)^{1/2}
     max.vec[j]<-max(abs(gamma.s/SE.norm))
   }
   critical.val<-quantile(max.vec,probs=0.99)
   return(critical.val)
 }


  
##### Functions: handling the (possible) union of CIs
##### CI.matrix is a matrix, where each row represents a CI
analysis.CI<-function(CI.matrix,true.val=1,grid.size=n^{-0.8}){
  #### number of rows in the CI.matrix
  d<-dim(CI.matrix)[1]
  CI.coverage<-0
  CI.len<-0
  grid.seq<-NULL
  for (l in 1: d){
    CI.len<-CI.len+CI.matrix[l,2]-CI.matrix[l,1]
    if((CI.matrix[l,2]>true.val)*(CI.matrix[l,1]<true.val)==1){
      CI.coverage<-1
    }
    grid.seq<-c(grid.seq,seq(CI.matrix[l,1],CI.matrix[l,2],by=grid.size))
  }
  return(list(CI.coverage=CI.coverage,CI.len=CI.len,grid.seq=grid.seq))
}

##### Functions: CI construction by searching method 
### INPUT: ITT_Y: a numeric, non-missing vector estimating each instrument's effect on the outcome
###        ITT_D: a numeric, non-missing vector estimating each instrument's effect on the treatment
###        SigmaSqY: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of Y in the reduced-form model of Y
###        SigmaSqD: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of D in the reduced-form model of D
###        SigmaYD: a numeric, non-missing, scalar value estimating the covariance
###                 of the error terms in the reduced-form models of Y and D
###        InitiSet: a set of pre-selected IVs (for majority rule: it is the set of relevant IVs; 
###        for plurality rule: it is a set of weakly violated IVs)
Searching.CI<-function(ITT_Y,ITT_D,SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,beta.grid=NULL,bootstrap=TRUE){
  pz<-ncol(WUMat)
  if(length(InitiSet)==0){
    warning("No valid IV! OLS is used.")
    CI.search=t(as.matrix(confint(lm(Y~D+W))[2,]))
    rule<-FALSE
    return(list(CI.search=CI.search,rule=rule))
  }else{
    threshold.size<-(length(InitiSet)/2)
    if(is.null(beta.grid)){
      beta.grid<-seq(-5,5,by=max(n,500)^{-1})
    }
    n.beta<-length(beta.grid)
    valid.grid<-rep(NA,n.beta)
    SE.norm<-sqrt(SigmaSqD * colSums(WUMat^2) /n^2)
        #(diag(solve(covW)/n)^{1/2})[1:pz]
    if(bootstrap==TRUE){
      Tn<-cut.off(SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,pz,beta.grid,N=1000)
    }else{
      Tn<-sqrt(2.005*log(n.beta))
    }
    for(j in 1:n.beta){
      b<-beta.grid[j]
      se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
      valid.grid[j]<-sum(abs(ITT_Y[InitiSet]-b*ITT_D[InitiSet])<Tn*se.b*SE.norm[InitiSet])
    }
    #CI.search<-c(min(beta.grid[which(valid.grid>threshold.size)]),max(beta.grid[which(valid.grid>threshold.size)]))
    if(length(beta.grid[which(valid.grid>threshold.size)])==0){
      ####### rule=FALSE indicates the corresponding rule fails
      rule=FALSE
      warning("Rule Fails. SS will give misleading CIs, SEs, and p-values.")
      sel.index<-which(valid.grid==max(valid.grid))
    }else{
      rule=TRUE
      sel.index<-which(valid.grid>threshold.size)
    }
    CI<-matrix(NA,nrow=length(sel.index),ncol=2)
    CI[,1]<-beta.grid[sel.index]
    upper.index<-sel.index+1
    upper.index[length(upper.index)]<-min(upper.index[length(upper.index)],n.beta)
    CI[,2]<-beta.grid[upper.index]
    if(dim(as.matrix(CI))[1]==1)
    {
      CI.search<-as.matrix(CI)
    }else{
      uni<- Intervals(CI)
      ###### construct the confidence interval by taking a union
      CI.search<-as.matrix(interval_union(uni))
    }
    return(list(CI.search=CI.search,rule=rule,valid.grid=valid.grid))
  }
}

###### Sampling and Searching 
##### Functions: CI construction by sampling and searching method 
### INPUT: ITT_Y: a numeric, non-missing vector estimating each instrument's effect on the outcome
###        ITT_D: a numeric, non-missing vector estimating each instrument's effect on the treatment
###        SigmaSqY: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of Y in the reduced-form model of Y
###        SigmaSqD: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of D in the reduced-form model of D
###        SigmaYD: a numeric, non-missing, scalar value estimating the covariance
###                 of the error terms in the reduced-form models of Y and D
###        InitiSet: a set of pre-selected IVs (for majority rule: it is the set of relevant IVs; 
###        for plurality rule: it is a set of weakly violated IVs)
Searching.CI.Sampling<-function(ITT_Y,ITT_D,SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,beta.grid=NULL,rho=NULL,M=1000,bootstrap=TRUE){
  pz<-ncol(WUMat)
  if(length(InitiSet)==0){
    warning("No valid IV! OLS is used.")
    CI.union=t(as.matrix(confint(lm(Y~D+W))[2,]))
    rule<-FALSE
    CI=CI.union
    return(list(CI.union=CI.union,rule=rule,CI=CI))
  }else{
    threshold.size<-(length(InitiSet)/2)
    if(is.null(rho)){
      rho<-(log(n)/M)^{1/(2*length(InitiSet))}/6
    }
    if(is.null(beta.grid)){
      beta.grid<-seq(-5,5,by=n^{-0.8})
    }
    n.beta<-length(beta.grid)
    valid.grid.sample<-matrix(NA,nrow=M,ncol=n.beta)
    SE.norm<-sqrt(SigmaSqD * colSums(WUMat^2) /n^2)
      #(diag(solve(covW)/n)^{1/2})[1:pz]
    #Cov.Y<-SigmaSqD*solve(covW)[1:pz,1:pz]/n
    #Cov.D<-SigmaSqY*solve(covW)[1:pz,1:pz]/n
    if(bootstrap==TRUE){
      Tn<-cut.off(SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,pz,beta.grid,N=1000)
    }else{
      Tn<-sqrt(2.005*log(n.beta))
    }
    unit.matrix<-t(WUMat)%*%WUMat/n^2
      #(solve(covW)/n)[1:pz,1:pz]
    Cov1<-cbind(SigmaSqD*unit.matrix,SigmaYD*unit.matrix)
    Cov2<-cbind(SigmaYD*unit.matrix,SigmaSqY*unit.matrix)
    Cov.total<-rbind(Cov1,Cov2)
    for(m in 1:M){
      Gen.mat<-mvrnorm(1, rep(0,2*pz), Cov.total)
      ITT_Y.sample<-ITT_Y-Gen.mat[(pz+1):(2*pz)]
      ITT_D.sample<-ITT_D-Gen.mat[1:pz]
      ###### generating indepedent copies
      #ITT_Y.sample<-ITT_Y-mvrnorm(1,rep(0,pz),Cov.Y)
      #ITT_D.sample<-ITT_D-mvrnorm(1,rep(0,pz),Cov.D)
      for(j in 1:n.beta){
        b<-beta.grid[j]
        se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
        valid.grid.sample[m,j]<-sum(abs(ITT_Y.sample[InitiSet]-b*ITT_D.sample[InitiSet])<rho*Tn*se.b*SE.norm[InitiSet])
      }
    }
    CI<-matrix(NA,nrow=M,ncol=2)
    for(m in 1:M){
      if(length(which(valid.grid.sample[m,]>threshold.size))>0){
        CI[m,1]<-min(beta.grid[which(valid.grid.sample[m,]>threshold.size)])
        CI[m,2]<-max(beta.grid[which(valid.grid.sample[m,]>threshold.size)])
      }
    }
    CI<-na.omit(CI)
    while((dim(as.matrix(CI))[1]<min(0.05*M,50)) && (rho<0.5)){
      #print(as.matrix(CI)[1])
      #print(rho)
      rho<-1.25*rho
      for(m in 1:M){
        Gen.mat<-mvrnorm(1, rep(0,2*pz), Cov.total)
        ITT_Y.sample<-ITT_Y-Gen.mat[(pz+1):(2*pz)]
        ITT_D.sample<-ITT_D-Gen.mat[1:pz]
        #ITT_Y.sample<-ITT_Y-mvrnorm(1,rep(0,pz),Cov.Y)
        #ITT_D.sample<-ITT_D-mvrnorm(1,rep(0,pz),Cov.D)
        for(j in 1:n.beta){
          b<-beta.grid[j]
          se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
          valid.grid.sample[m,j]<-sum(abs(ITT_Y.sample[InitiSet]-b*ITT_D.sample[InitiSet])<rho*Tn*se.b*SE.norm[InitiSet])
        }
      }
      CI<-matrix(NA,nrow=M,ncol=2)
      for(m in 1:M){
        if(length(which(valid.grid.sample[m,]>threshold.size))>0){
          CI[m,1]<-min(beta.grid[which(valid.grid.sample[m,]>threshold.size)])
          CI[m,2]<-max(beta.grid[which(valid.grid.sample[m,]>threshold.size)])
        }
      }
      CI<-na.omit(CI)
    }
    rule<-TRUE
    if(dim(as.matrix(CI))[1]==0){
      rule<-FALSE
      CI.union<-t(as.matrix(c(min(beta.grid),max(beta.grid))))
    }else if(dim(as.matrix(CI))[1]==1)
    {
      CI.union<-as.matrix(CI)
    }else{
      uni<- Intervals(CI)
      ###### construct the confidence interval by taking a union
      CI.union<-as.matrix(interval_union(uni))
    }
    #CI.upper<-quantile(CI[,2],0.975)
    #CI.lower<-quantile(CI[,1],0.025)
    #CI.quantile<-c(CI.lower,CI.upper)
    return(list(CI.union=CI.union,rule=rule,CI=CI))
  }
}

### TSHT.Initial (This modifies the original TSHT function and produces an initial esitmator for sampling and searching)
### FUNCTION: Estimates the set of valid and relevant IVs
### INPUT: ITT_Y: a numeric, non-missing vector estimating each instrument's effect on the outcome
###        ITT_D: a numeric, non-missing vector estimating each instrument's effect on the treatment
###        WUMat: an n by pz numeric, non-missing matrix that computes the precision of W 
###               (ex1 t(WUMat) %*% WUMat / n = (W^T W/n)^(-1))
###               (ex2 U = (W^T W/n)^(-1))
###        SigmaSqY: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of Y in the reduced-form model of Y
###        SigmaSqD: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of D in the reduced-form model of D
###        SigmaYD: a numeric, non-missing, scalar value estimating the covariance
###                 of the error terms in the reduced-form models of Y and D
###        tuning, a numeric scalar value tuning parameter for TSHT greater 
###                than 2 (default = 2.01)
### OUTPUT: a list (a) VHat (numeric vector denoting the set of valid and relevant IVs)
###                (b) SHat (numeric vector denoting the set of relevant IVs)
TSHT.Initial <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,bootstrap=FALSE,tuning = 2.01) {
  # Check ITT_Y and ITT_D
  stopifnot(!missing(ITT_Y),!missing(ITT_D),length(ITT_Y) == length(ITT_D))
  stopifnot(all(!is.na(ITT_Y)),all(!is.na(ITT_D)))
  ITT_Y = as.numeric(ITT_Y); ITT_D = as.numeric(ITT_D)
  
  # Check WUMat
  stopifnot(!missing(WUMat),is.matrix(WUMat), nrow(WUMat) > 1, ncol(WUMat) == length(ITT_Y)) 
  stopifnot(all(!is.na(WUMat)))
  
  # Check Sigmas 
  stopifnot(!missing(SigmaSqY), is.numeric(SigmaSqY), length(SigmaSqY) == 1, !is.na(SigmaSqY), SigmaSqY > 0)
  stopifnot(!missing(SigmaSqD), is.numeric(SigmaSqD), length(SigmaSqD) == 1, !is.na(SigmaSqD),SigmaSqD > 0)
  stopifnot(!missing(SigmaYD),is.numeric(SigmaYD), length(SigmaYD) == 1,!is.na(SigmaYD))
  
  # Other Input check
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)
  
  # Constants
  n = nrow(WUMat); 
  pz = length(ITT_Y)
  # First Stage
  if(bootstrap==TRUE){
  Tn<-min(cut.off.IVStr(SigmaSqD,WUMat,pz,cut.prob = 0.95),sqrt(log(n))) ### this can be modified by the user
  SE.norm<-(diag(solve(covW)/n)^{1/2})[1:pz]
  SHat<-(1:pz)[abs(ITT_D)>Tn*sqrt(SigmaSqD)*SE.norm]
  }else{
  SHat = (1:pz)[(abs(ITT_D) >= (sqrt(SigmaSqD * colSums(WUMat^2) /n) * sqrt(tuning*log(pz)/n)))]
  }
  if(length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. TSHT with these IVs will give misleading CIs, SEs, and p-values. Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE
  
  # Second Stage
  # pi.candidate is the estimated value of pi across different candidates
  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat
  for(j in SHat) {
    beta.j = ITT_Y[j] / ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    sigmasq.j = SigmaSqY + beta.j^2 * SigmaSqD - 2* beta.j * SigmaYD
    PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j * colSums( (WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) * 
      sqrt(tuning^2 * log(pz)/n)
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }
  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }
  # Voting
  #VM.1 = apply(VHats.bool,1,sum)
  #VM.2 = apply(VHats.bool,2,sum)
  #VM<-VM.1
  #for(l in 1:length(VM.1)){
  # VM[l]=min(VM.1[l],VM.2[l])
  #}
  VM= apply(VHats.boot.sym,1,sum)
  VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
  VM.p = rownames(VHats.boot.sym)[max(VM) == VM] #Plurality winners
  V.set<-NULL
  for(index in union(VM.m,VM.p)){
    V.set<-union(V.set,names(which(VHats.boot.sym[index,]==1)))
  }
  VHat<-NULL
  for(index in V.set){
    VHat<-union(VHat,names(which(VHats.boot.sym[,index]==1)))
  }
  VHat=as.numeric(VHat)
  return(list(SHat=SHat,VHat=VHat,V.set=V.set))
}