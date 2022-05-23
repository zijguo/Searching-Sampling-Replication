########################################
### Date: 05/22/2022
### Author: Zhenyu WANG
### Aim: compute threshold max{T_jk} for S2

library(MASS)
library(intervals)
library(AER)
library(sandwich)
source("src/helpers.R")
library(ggplot2)

### Compute V.Gamma, V.gamma, C
n=30000  # {500, 1000, 2000, 5000}
VIO.str = 0.1
IV.str = 0.5
setting = "S2"

pi.value = IV.str*VIO.str
beta = 1

case = "homo"
px = 10
if(setting=="S2"){
  L = 10; 
  s1 = 2; s2 = 4; s=s1+s2
  alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(1,s2)/3)
}

gamma=rep(IV.str,L)
p=L+px # p stands for the number of total exogeneous variables 
phi<-rep(0,px)
psi<-rep(0,px)
phi[1:px]<-(1/px)*seq(1,px)+0.5
psi[1:px]<-(1/px)*seq(1,px)+1
rho=0.5
A1gen <- function(rho, p){
  A1 = matrix(0, nrow=p, ncol=p)
  for(i in 1:p) for(j in 1:p) A1[i, j] = rho^(abs(i-j))
  return(A1)
}
Cov<-(A1gen(rho,p))

set.seed(1)
W = mvrnorm(n, rep(0, p), Cov)
Z = W[, 1:L]
X = W[, (L+1):p]
if(case=="hetero"){
  epsilon1 = rnorm(n)
  tao1 = rep(NA, n); for(i.n in 1:n) tao1[i.n] = rnorm(n=1, mean=0, sd=0.25+0.5*(Z[i.n, 1])^2)
  tao2 = rnorm(n)
  epsilon2 = 0.3*epsilon1 + sqrt((1-0.3^2)/(0.86^4+1.38072^2))*(1.38072*tao1+0.86^2*tao2)
}else if(case=="homo"){
  epsilonSigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)
  epsilon = mvrnorm(n, rep(0, 2), epsilonSigma)
  epsilon1 = epsilon[,1]
  epsilon2 = epsilon[,2]
}
D = 0.5 + Z %*% gamma+ X%*% psi + epsilon1
Y = -0.5 + Z %*% alpha + D * beta + X%*%phi+ epsilon2
if(is.null(X)) W = Z else W = cbind(Z, X)
n = length(Y); pz = ncol(Z); p = ncol(W)
intercept = TRUE
if(intercept) W = cbind(W, 1)
covW = t(W)%*%W/n
U = solve(covW) # precision matrix
WUMat = (W%*%U)[,1:pz]
## OLS estimators
qrW = qr(W)
ITT_Y = qr.coef(qrW, Y)[1:pz]
ITT_D = qr.coef(qrW, D)[1:pz]
resid_Y = as.vector(qr.resid(qrW, Y))
resid_D = as.vector(qr.resid(qrW, D))
V.Gamma = (t(WUMat)%*%diag(resid_Y^2)%*%WUMat)/n
V.gamma = (t(WUMat)%*%diag(resid_D^2)%*%WUMat)/n
C = (t(WUMat)%*%diag(resid_Y * resid_D)%*%WUMat)/n
out = list()
out$V.Gamma = V.Gamma; out$V.gamma = V.gamma; out$C = C

########## Compute Threshold ###########
VIO.str.options = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
n.options = c(500,2000)
Thres.mat.list=list()
table.out = matrix(NA, nrow=length(n.options)*length(VIO.str.options), ncol=6)
for(i.VIO.str in 1:length(VIO.str.options)){
  VIO.str = VIO.str.options[i.VIO.str]
  for(i.n in 1:length(n.options)){
    n = n.options[i.n]
    ind = (i.VIO.str-1)*length(n.options)+i.n
    
    IV.str = 0.5
    setting = "S2"
    
    pi.value = IV.str*VIO.str
    beta = 1
    
    case = "homo"
    px = 10
    if(setting=="S2"){
      L = 10; 
      s1 = 2; s2 = 4; s=s1+s2
      alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(1,s2)/3)
    }
    
    gamma=rep(IV.str,L)
    p=L+px # p stands for the number of total exogeneous variables 
    phi<-rep(0,px)
    psi<-rep(0,px)
    phi[1:px]<-(1/px)*seq(1,px)+0.5
    psi[1:px]<-(1/px)*seq(1,px)+1
    
    Gamma_true = beta * alpha + gamma
    beta_thres = Gamma_true / gamma
    # print(beta_thres)
    Thres.mat = Thres.mat.0 = matrix(NA, nrow=L, ncol=L)
    for(j in 1:L){
      R = out$V.Gamma + beta_thres[j]^2*out$V.gamma - 2*beta_thres[j]*out$C
      # R = inputs[[i.VIO.str]]$V.Gamma + beta_thres[j]^2*inputs[[i.VIO.str]]$V.gamma - 2*beta_thres[j]*inputs[[i.VIO.str]]$C
      for(k in 1:L){
        Thres.mat.0[j, k] = sqrt(1/n*(R[k, k]/gamma[k]^2 + R[j, j]/gamma[j]^2 - 2*R[j,k]/(gamma[k]*gamma[j])))
      }
    }
    for(j in 1:L) for(k in 1:L) Thres.mat[j,k] = min(Thres.mat.0[j, k], Thres.mat.0[k, j])
    Thres.mat.list[[ind]] = Thres.mat
    num1=max(Thres.mat[1:4,1:10])
    num2=max(Thres.mat)
    num3=max(Thres.mat[1:4,1:4])
    num4=num1*2*sqrt(log(n))
    num5=num2*2*sqrt(log(n))
    num6=num3*2*sqrt(log(n))
    table.out[ind,] = c(num1, num2, num3, num4, num5, num6)
  }
}
col1=rep(VIO.str.options, rep(length(n.options),length(VIO.str.options)))
col2=rep(n.options, length(VIO.str.options))
table.out = cbind(col1, col2, table.out)
colnames(table.out) = c("VIO.str","n","maxV","maxAll","maxVV","n-maxV","n-maxAll","n-maxVV")
table.out

thres_table = table.out
save.image("RDatas/S2_thres.RData")
#save.image("S2_thres_0517.RData") #save.image("S2_thres.RData")
