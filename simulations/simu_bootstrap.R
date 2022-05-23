###########################
### Author: Zhenyu Wang
### Date: 05/13
### Todo: 
### 1.bootstrap vs. not bootstrap
library(MASS)
library(intervals)
source("src/helpers.R")
source("src/TSHT-ldim.R")
source("src/invalidIV.R")
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")
########### Part-1: Specify Settings and Parameters ##########
# parameters changing
n=2000  # {500, 1000, 2000, 5000}
VIO.str = 0.4 # {0.2, 0.4}
setting = "S4" # {"S1", "S2", "S3", "S4", "S5"}
sim.round = 1 # {1}

nsim = 500
IV.str = 0.5
####################
pi.value = IV.str*VIO.str
beta = 1

case = "homo"
px = 10
if(setting=="S1"){
  L = 10; 
  s1 = s2 = 2; s = s1+s2
  alpha = c(rep(0,L-s1-s2),rep(pi.value,s1),-seq(1,s2)/2)
}
if(setting=="S2"){
  L = 10; 
  s1 = 2; s2 = 4; s=s1+s2
  alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(1,s2)/3)
}
if(setting=="S3"){
  L = 10; 
  s1 = 2; s2 = 4; s=s1+s2
  alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(1,s2)/6)
}
if(setting=="S4"){
  L = 6; 
  s1 = s2 = s3 = s4 = 1; s=s1+s2+s3+s4
  alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(0.8,s2),-seq(0.4,s3),seq(0.6,s4))
}
if(setting=="S5"){
  L = 6; 
  s1 = s2 = s3= s4 = 1; s=s1+s2+s3+s4
  alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(0.8,s2),-seq(0.4,s3),seq(pi.value+0.1,s4));
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

########### Part-2: Simulations ##########
# matrix storing info
Cov.mat = Leng.mat = matrix(NA, nrow=nsim, ncol=2)
colnames(Cov.mat) = colnames(Leng.mat) = c("B-Y", "B-N")
CI.mat = matrix(NA, nrow=nsim, ncol=4)
colnames(CI.mat) = c("B-Y-L","B-Y-U","B-N-L","B-N-U")
Rule.mat = matrix(NA, nrow=nsim, ncol=2)
colnames(Rule.mat) = c("B-Y", "B-N")

for(i.sim in 1:nsim){
  set.seed(i.sim+(sim.round-1)*nsim)
  print(i.sim)
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
  
  ######### Preparation ###########
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

  ######## Select VHat ########
  ## Use TSHT method to select
  TSHT.out <- TSHT.Init(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma)
  VHat.TSHT = sort(TSHT.out$VHat)
  SHat.TSHT = sort(TSHT.out$SHat)

  ####### beta.grid  ########
  # TSHT method
  var.beta = 1/n * (diag(V.Gamma)/ITT_D^2 + diag(V.gamma)*ITT_Y^2/ITT_D^4 - 2*diag(C)*ITT_Y/ITT_D^3)
  var.beta = var.beta[VHat.TSHT]
  CI.init = matrix(NA, nrow=length(VHat.TSHT), ncol=2)
  CI.init[,1] = (ITT_Y/ITT_D)[VHat.TSHT] - sqrt(log(n)*var.beta)
  CI.init[,2] = (ITT_Y/ITT_D)[VHat.TSHT] + sqrt(log(n)*var.beta)
  uni = Intervals(CI.init)
  CI.init.union = as.matrix(interval_union(uni))
  beta.grid.TSHT = grid.CI(CI.init.union, grid.size=n^(-0.6))

  ####### Column Group-2 #########
  # TSHT Searching - bootstrap TRUE
  out = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = VHat.TSHT,
                     beta.grid = beta.grid.TSHT, bootstrap = TRUE)
  CI1 = out$CI; rule1 = out$rule
  CI.mat[i.sim, 1:2] = CI1
  Cov.mat[i.sim, 1] = sum((CI1[,1]<beta)*(CI1[,2]>beta))
  Leng.mat[i.sim, 1] = sum(CI1[,2] - CI1[,1])
  Rule.mat[i.sim, 1] = rule1
  
  # TSHT Searching - bootstrap FALSE
  out = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = VHat.TSHT,
                     beta.grid = beta.grid.TSHT, bootstrap = FALSE)
  CI1 = out$CI; rule1 = out$rule
  CI.mat[i.sim, 3:4] = CI1
  Cov.mat[i.sim, 2] = sum((CI1[,1]<beta)*(CI1[,2]>beta))
  Leng.mat[i.sim, 2] = sum(CI1[,2] - CI1[,1])
  Rule.mat[i.sim, 2] = rule1
}

rm(list=setdiff(ls(), c("setting", "IV.str", "VIO.str", "n", "Cov.mat", "Leng.mat", "CI.mat", "Rule.mat")))
filename = paste("Bootstrap-Homo-Setting", setting, "-Strength", IV.str, "-Violation", VIO.str, "-n", n,"-SimRound",sim.round, ".RData", sep="")
save.image(filename)
