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
