### Load R Packages 
```R
library(AER)
library(MASS)
library(sandwich)
library(intervals);
```
### Load source codes
```R
source('Searching-Sampling.R', encoding = 'UTF-8')
source('TSHT-ldim.R', encoding = 'UTF-8')
```

### Generate the design covariance matrix
```R
# simulate A1 of the matrix (rho)^|i-j|
A1gen<-function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    } 
  }
  A1
}
```
### Set the model parameter
```R
# n: sample size
# IV.str: individual IV strength
# VIO.str: violation strength
# beta:   true treatment effect (set at 1)
# px:number of covariates
# L: number of candiate IVs
n = 500;
IV.str=0.5;
VIO.str=0.4;
pi.value<-IV.str*VIO.str
beta = 1; 
px <- 10;
L = 10; 
p=L+px
phi<-rep(0,px)
psi<-rep(0,px)
phi[1:px]<-(1/px)*seq(1,px)+0.5
psi[1:px]<-(1/px)*seq(1,px)+1
rho=0.5
Cov<-(A1gen(rho,p))
### Generate a setting where only the plurality rule holds
# 6 invalid IVs and 4 valid IVs 
s1 = 3;
s2 = 3;
s=s1+s2
alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(1,s2)/2);
gamma=rep(IV.str,L)
```


### set up the sampling time
```R
M0=1000
```
### set up whether bootstrap is applied to choose the threshold
```R
boot.value=TRUE
```



### Generate the data
```R
# Y: n by 1 vector of outcomes (must be continuous)
# D: n by 1 vector of treatments (continuous or discrete)
# Z: n by L matrix of instruments (continuous or discrete)
# X: n by px matrix of baseline covariates (continuous or discrete)
W<-mvrnorm(n, rep(0, p), Cov)
Z=W[,1:L]
X=W[,(L+1):p]
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D = 0.5 + Z %*% gamma+ X%*% psi + epsilon[,1]
Y = -0.5 + Z %*% alpha + D * beta + X%*%phi+ epsilon[,2]
```

### Compute the reduced form estimaors with the corresponding covariance 
```R
pz<-ncol(Z)
W = cbind(W,1)
# Compute covariance of W and W %*% U
covW = t(W) %*% W /n #this should automatically turn covW into a matrix
WUMat = W %*% (solve(covW))[,1:pz]
qrW = qr(W)
ITT_Y = qr.coef(qrW,Y)[1:pz]
ITT_D = qr.coef(qrW,D)[1:pz]
SigmaSqY = sum(qr.resid(qrW,Y)^2)/(n -p-1)
SigmaSqD = sum(qr.resid(qrW,D)^2)/(n -p-1)
SigmaYD = sum(qr.resid(qrW,Y) * qr.resid(qrW,D)) / (n - p-1)
######## select strongly associated IV
#Tn<-sqrt(log(n)) ### this can be modified by the user
#SE.norm<-(diag(solve(covW)/n)^{1/2})[1:pz]
#Shat<-(abs(ITT_D)>Tn*sqrt(SigmaSqD)*SE.norm)
Tn<-min(cut.off.IVStr(SigmaSqD,WUMat,pz),sqrt(log(n))) ### this can be modified by the user
SE.norm<-(diag(solve(covW)/n)^{1/2})[1:pz]
Shat<-(abs(ITT_D)>Tn*sqrt(SigmaSqD)*SE.norm)
```


### screen out strongly invalid IVs and retain a set of valid and weakly invalid IVs 
```R
V0.hat<-TSHT.Initial(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD)$VHat
V0.hat<-sort(V0.hat)
```

### construct the initial range [L,U]  
```R
temp<-(SigmaSqY)/(ITT_D[V0.hat]^2)+SigmaSqD*(ITT_Y[V0.hat]^2)/(ITT_D[V0.hat]^4)-2*SigmaYD*(ITT_Y[V0.hat])/(ITT_D[V0.hat]^3)
var.beta<-(diag(solve(covW)/n)[1:pz])[V0.hat]*temp
CI.initial<-matrix(NA,nrow=length(V0.hat),ncol=2)
CI.initial[,1]<-(ITT_Y/ITT_D)[V0.hat]-sqrt(log(n)*var.beta)
CI.initial[,2]<-(ITT_Y/ITT_D)[V0.hat]+sqrt(log(n)*var.beta)
uni<- Intervals(CI.initial)
CI.initial.union<-as.matrix(interval_union(uni))
beta.grid.seq<-analysis.CI(CI.initial.union,grid.size=n^{-1})$grid.seq
```

### conduct the initial searching and output a refined range [L,U]

```R
CI.sea<-Searching.CI(ITT_Y,ITT_D,SigmaSqD,SigmaSqY,SigmaYD,V0.hat,WUMat,beta.grid=beta.grid.seq,bootstrap=FALSE)
CI.temp<-CI.sea$CI.search
beta.grid<-analysis.CI(as.matrix(CI.sea$CI.search),beta,n^{-0.6})$grid.seq
```

### conduct the refined searching 
```R
CI.sea.refined<-Searching.CI(ITT_Y,ITT_D,SigmaSqD,SigmaSqY,SigmaYD,V0.hat,WUMat,beta.grid,bootstrap=boot.value)
CI.temp<-CI.sea.refined$CI.search
``` 
### conduct the refined sampling

```R
CI.sampling<-Searching.CI.Sampling(ITT_Y,ITT_D,SigmaSqD,SigmaSqY,SigmaYD,V0.hat,WUMat,beta.grid,M=M0,bootstrap=boot.value)
CI.temp<-c(min(CI.sampling$CI.union[,1]),max(CI.sampling$CI.union[,2]))
```
