###########################
### Author: Zhenyu Wang
### Date: 05/22/2022
### Aim: investigate varying violation strengths for S2 when hetero errors.

library(MASS)
library(intervals)
source("src/helpers.R")
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")
########### Part-1: Specify Settings and Parameters ##########
# parameters changing
n=500  # {500, 1000, 2000, 5000}
VIO.str = 0.025 # {0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0}
sim.round = 1
IV.str = 0.5 
setting = "S2"
# the number of simulations
nsim=500
####################
pi.value = IV.str*VIO.str
beta = 1
case = "hetero"
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
Cov.mat = Leng.mat = matrix(NA, nrow=nsim, ncol=12)
colnames(Cov.mat) = colnames(Leng.mat) = c("oracle-1", "oracle-2" ,"TSHT-1", "TSHT-2", "CIIV-1", "CIIV-2",
                                           "Sear-TSHT", "Sear-CIIV", "Sear-Comb",
                                           "Samp-TSHT", "Samp-CIIV", "Samp-Comb")
Rule.mat = matrix(NA, nrow=nsim, ncol=4)
colnames(Rule.mat) = c("Sear-TSHT", "Sear-CIIV", "Samp-TSHT", "Samp-CIIV")
CI.mat = matrix(NA, nrow=nsim, ncol=12)
colnames(CI.mat) = c("Oracle1-L","Oracle1-U","Oracle2-L","Oracle2-U",
                     "TSHT1-L","TSHT1-U","TSHT2-L","TSHT2-U",
                     "CIIV1-L","CIIV1-U","CIIV2-L","CIIV2-U")

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
  TSHT.out.small = TSHT.Init.small(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma)
  VHat.TSHT.small = sort(TSHT.out.small$VHat)
  ## Use TSHT method to select
  TSHT.out <- TSHT.Init(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma)
  VHat.TSHT = sort(TSHT.out$VHat)
  ## Use CIIV method to select
  CIIV.temp = tryCatch_E(CIIV(Y,D,Z,X,robust=TRUE), ret.obj = list(NULL))
  if (is.null(CIIV.temp$error)) CIIV.result <- CIIV.temp$value
  VHat.CIIV = sapply(CIIV.result$Valid_Instruments, function(x) as.numeric(substring(x, 2)))
  if(!all(VHat.CIIV %in% (1:pz))){
    # Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)/2))
    # SHat = (1:pz)[abs(ITT_D) > (Tn * sqrt(diag(V.gamma)/n))]
    VHat.CIIV = VHat.TSHT
  }
  ####### Column Group-1 ########
  # Oracle Method
  # Oracle Method
  CI = robust_method(Y, regressor=cbind(D, X, Z[,-(1:(L-s))]), Exogenous = cbind(Z, X))
  CI.mat[i.sim, 1:2] = CI
  Cov.mat[i.sim, 1] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 1] = CI[2] - CI[1]
  CI = robust_method_new(WUMat, seq(L-s), ITT_Y, ITT_D, V.Gamma, V.gamma, C)
  CI.mat[i.sim, 3:4] = CI
  Cov.mat[i.sim, 2] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 2] = CI[2] - CI[1]
  # out = ivmodel(Y, D, Z[,1:(L-s)], cbind(Z[,-(1:(L-s))], X), intercept=TRUE, heteroSE=TRUE)
  # CI = confint(out)["TSLS",]
  # Cov.mat[i.sim, 2] = (CI[1]<beta)*(CI[2]>beta)
  # Leng.mat[i.sim, 2] = CI[2] - CI[1]
  
  # TSHT Method
  CI = robust_method(Y, regressor=cbind(D, X, Z[,-VHat.TSHT.small]), Exogenous = cbind(Z, X))
  CI.mat[i.sim, 5:6] = CI
  Cov.mat[i.sim, 3] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 3] = CI[2] - CI[1]
  CI = robust_method_new(WUMat, VHat.TSHT.small, ITT_Y, ITT_D, V.Gamma, V.gamma, C)
  CI.mat[i.sim, 7:8] = CI
  # out = ivmodel(Y, D, Z[,VHat.TSHT], cbind(Z[,-VHat.TSHT], X), intercept=TRUE, heteroSE=TRUE)
  # CI = confint(out)["TSLS",]
  Cov.mat[i.sim, 4] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 4] = CI[2] - CI[1]
  
  # CIIV Method
  CI.mat[i.sim, 9:10] = CIIV.result$ci_CIM
  Cov.mat[i.sim, 5] <- ifelse(beta<CIIV.result$ci_CIM[2]&beta>CIIV.result$ci_CIM[1], 1, 0)
  Leng.mat[i.sim, 5] <- CIIV.result$ci_CIM[2] - CIIV.result$ci_CIM[1]
  CI = robust_method_new(WUMat, VHat.CIIV, ITT_Y, ITT_D, V.Gamma, V.gamma, C)
  CI.mat[i.sim, 11:12] = CI
  Cov.mat[i.sim, 6] <- (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 6] <- CI[2] - CI[1]
  
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
  
  # CIIV method
  var.beta = 1/n * (diag(V.Gamma)/ITT_D^2 + diag(V.gamma)*ITT_Y^2/ITT_D^4 - 2*diag(C)*ITT_Y/ITT_D^3)
  var.beta = var.beta[VHat.CIIV]
  CI.init = matrix(NA, nrow=length(VHat.CIIV), ncol=2)
  CI.init[,1] = (ITT_Y/ITT_D)[VHat.CIIV] - sqrt(log(n)*var.beta)
  CI.init[,2] = (ITT_Y/ITT_D)[VHat.CIIV] + sqrt(log(n)*var.beta)
  uni = Intervals(CI.init)
  CI.init.union = as.matrix(interval_union(uni))
  beta.grid.CIIV = grid.CI(CI.init.union, grid.size=n^(-0.6))
  
  ####### Column Group-2 #########
  # TSHT Searching
  out = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = VHat.TSHT,
                     beta.grid = beta.grid.TSHT)
  CI1 = out$CI; rule1 = out$rule
  Cov.mat[i.sim, 7] = sum((CI1[,1]<beta)*(CI1[,2]>beta))
  Leng.mat[i.sim, 7] = sum(CI1[,2] - CI1[,1])
  # CIIV Searching
  out = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = VHat.CIIV,
                     beta.grid = beta.grid.CIIV)
  CI2 = out$CI; rule2 = out$rule
  Cov.mat[i.sim, 8] = sum((CI2[,1]<beta)*(CI2[,2]>beta))
  Leng.mat[i.sim, 8] = sum(CI2[,2] - CI2[,1])
  # Comb Searching
  CI3 = Intervals(rbind(CI1, CI2))
  CI3 = as.matrix(interval_union(CI3))
  Cov.mat[i.sim, 9] = sum((CI3[,1]<beta)*(CI3[,2]>beta))
  Leng.mat[i.sim, 9] = sum(CI3[,2] - CI3[,1])
  # Rule
  Rule.mat[i.sim, 1:2] = c(rule1, rule2)
  
  ####### Column Group-3 #########
  # TSHT Sampling
  out = Searching.CI.sampling(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet=VHat.TSHT,
                              beta.grid = beta.grid.TSHT)
  CI1 = out$CI; rule1 = out$rule
  Cov.mat[i.sim, 10] = sum((CI1[,1]<beta)*(CI1[,2]>beta))
  Leng.mat[i.sim, 10] = sum(CI1[,2] - CI1[,1])
  # CIIV Sampling
  out = Searching.CI.sampling(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet=VHat.CIIV,
                              beta.grid = beta.grid.CIIV)
  CI2 = out$CI; rule2 = out$rule
  Cov.mat[i.sim, 11] = sum((CI2[,1]<beta)*(CI2[,2]>beta))
  Leng.mat[i.sim, 11] = sum(CI2[,2] - CI2[,1])
  # Comb Sampling
  CI3 = Intervals(rbind(CI1, CI2))
  CI3 = as.matrix(interval_union(CI3))
  Cov.mat[i.sim, 12] = sum((CI3[,1]<beta)*(CI3[,2]>beta))
  Leng.mat[i.sim, 12] = sum(CI3[,2] - CI3[,1])
  # Rule
  Rule.mat[i.sim, 3:4] = c(rule1, rule2)
}

rm(list=setdiff(ls(), c("setting", "IV.str", "VIO.str", "n", "Cov.mat", "Leng.mat", "Rule.mat", "CI.mat", "sim.round")))
filename = paste("Hetero-Setting", setting, "-Strength", IV.str, "-Violation", VIO.str, "-n", n,"-SimRound",sim.round, ".RData", sep="")
save.image(filename)
