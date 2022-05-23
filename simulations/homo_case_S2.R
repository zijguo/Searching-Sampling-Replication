###########################
### Author: Zhenyu Wang
### Date: 05/22/2022
### Aim: investigate varying violation strengths for S2 when homo errors.

library(MASS)
library(intervals)
source("src/helpers.R")
source("src/TSHT-ldim.R")
source("src/invalidIV.R")
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")
########### Part-1: Specify Settings and Parameters ##########
# parameters changing
n=1000  # {500, 1000, 2000, 5000}
VIO.str = 1 # {0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5}
sim.round = 1 # seq(1,5)

nsim = 100
IV.str = 0.5
setting = "S2"
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
Cov.mat = Leng.mat = matrix(NA, nrow=nsim, ncol=9)
colnames(Cov.mat) = colnames(Leng.mat) = c("oracle", "TSHT", "CIIV", "Sear-TSHT", "Sear-CIIV", "Sear-Comb",
                                           "Samp-TSHT", "Samp-CIIV", "Samp-Comb")
CI.mat = matrix(NA, nrow=nsim, ncol=4)
colnames(CI.mat) = c("TSHT-L","TSHT-U","CIIV-L","CIIV-U",)
Rule.mat = matrix(NA, nrow=nsim, ncol=4)
colnames(Rule.mat) = c("Sear-TSHT", "Sear-CIIV", "Samp-TSHT", "Samp-CIIV")
TSHT.VHat.mat = matrix(NA, nrow=nsim, ncol=L)

comp.union=TRUE
comp.union.median=TRUE
method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR")
length.union05<-matrix(NA, nrow=nsim, ncol=10)
cover.union05<-matrix(NA, nrow=nsim, ncol=10)
colnames(length.union05)=c(method,method)
colnames(cover.union05)=c(method,method)

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
  
  ######## Column Group-1 #########
  # Oracle Method
  CI = confint(ivreg(Y ~ D + Z[,-(1:(L-s))]+ X  | Z + X ))[2,]
  Cov.mat[i.sim, 1] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 1] = CI[2] - CI[1]

  # TSHT Method
  TSHT.temp = TSHT(Y, D, Z, X)
  CI = TSHT.temp$ci
  CI.mat[i.sim, 1:2] = CI
  Cov.mat[i.sim, 2] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 2] = CI[2] - CI[1]
  # print(paste("simulation-", i.sim))
  TSHT.VHat.mat[i.sim, 1:length(TSHT.temp$VHat)] = TSHT.temp$VHat

  # CIIV Method
  CIIV.temp = tryCatch_E(CIIV(Y,D,Z,X,robust=FALSE), ret.obj = list(NULL))
  if (is.null(CIIV.temp$error)) CIIV.result <- CIIV.temp$value
  CI.mat[i.sim, 3:4] = CIIV.result$ci_CIM
  Cov.mat[i.sim,3] <- ifelse(beta<CIIV.result$ci_CIM[2]&beta>CIIV.result$ci_CIM[1], 1, 0)
  Leng.mat[i.sim, 3] <- CIIV.result$ci_CIM[2] - CIIV.result$ci_CIM[1]

  ######## Select VHat ########
  ## Use TSHT method to select
  TSHT.out <- TSHT.Init(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma)
  VHat.TSHT = sort(TSHT.out$VHat)

  ## Use CIIV method to select
  VHat.CIIV = sapply(CIIV.result$Valid_Instruments, function(x) as.numeric(substring(x, 2)))
  # CIIV.rule[i.sim] = all(VHat.CIIV %in% (1:pz))
  if(!all(VHat.CIIV %in% (1:pz))){
    VHat.CIIV = VHat.TSHT
  }
  # valid.index = (1:ncol(Z))[colnames(Z) %in% CIIV.result$Valid_Instruments]

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
  Cov.mat[i.sim, 4] = sum((CI1[,1]<beta)*(CI1[,2]>beta))
  Leng.mat[i.sim, 4] = sum(CI1[,2] - CI1[,1])
  # CIIV Searching
  out = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = VHat.CIIV,
                     beta.grid = beta.grid.CIIV)
  CI2 = out$CI; rule2 = out$rule
  Cov.mat[i.sim, 5] = sum((CI2[,1]<beta)*(CI2[,2]>beta))
  Leng.mat[i.sim, 5] = sum(CI2[,2] - CI2[,1])
  # Comb Searching
  CI3 = Intervals(rbind(CI1, CI2))
  CI3 = as.matrix(interval_union(CI3))
  Cov.mat[i.sim, 6] = sum((CI3[,1]<beta)*(CI3[,2]>beta))
  Leng.mat[i.sim, 6] = sum(CI3[,2] - CI3[,1])
  # Rule
  Rule.mat[i.sim, 1:2] = c(rule1, rule2)

  ####### Column Group-3 #########
  # TSHT Sampling
  out = Searching.CI.sampling(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet=VHat.TSHT,
                              beta.grid = beta.grid.TSHT)
  CI1 = out$CI; rule1 = out$rule
  Cov.mat[i.sim, 7] = sum((CI1[,1]<beta)*(CI1[,2]>beta))
  Leng.mat[i.sim, 7] = sum(CI1[,2] - CI1[,1])
  # CIIV Sampling
  out = Searching.CI.sampling(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet=VHat.CIIV,
                              beta.grid = beta.grid.CIIV)
  CI2 = out$CI; rule2 = out$rule
  Cov.mat[i.sim, 8] = sum((CI2[,1]<beta)*(CI2[,2]>beta))
  Leng.mat[i.sim, 8] = sum(CI2[,2] - CI2[,1])
  # Comb Sampling
  CI3 = Intervals(rbind(CI1, CI2))
  CI3 = as.matrix(interval_union(CI3))
  Cov.mat[i.sim, 9] = sum((CI3[,1]<beta)*(CI3[,2]>beta))
  Leng.mat[i.sim, 9] = sum(CI3[,2] - CI3[,1])
  # Rule
  Rule.mat[i.sim, 3:4] = c(rule1, rule2)

  ####### Column 4 #########
  # comp.union
  if(comp.union){
    timestart = Sys.time()
    s_bar=L-1
    print(paste("s_bar is",s_bar))
    resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                                ret.obj = list(NULL))
    if (is.null(resultOut.temp$error)) {
      resultOut <- resultOut.temp$value
    } else {
      next
    }
    resultOutStats = statsForInvalidIVCI(resultOut,trueBetaGrid = 0, method=method, trueBeta = beta)
    length.union05[i.sim,1:5] = resultOutStats$medianLength
    cover.union05[i.sim,1:5] = resultOutStats$coverageProportion
  }
  # comp.union.median
  if(comp.union.median){
    timestart = Sys.time()
    s_bar=ceiling(L/2)
    print(paste("s_bar is",s_bar))
    resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, X, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
                                ret.obj = list(NULL))
    if (is.null(resultOut.temp$error)) {
      resultOut <- resultOut.temp$value
    } else {
      next
    }
    resultOutStats = statsForInvalidIVCI(resultOut,trueBetaGrid = 0, method=method, trueBeta = beta)
    length.union05[i.sim,6:10] = resultOutStats$medianLength
    cover.union05[i.sim,6:10] = resultOutStats$coverageProportion
  }
}

rm(list=setdiff(ls(), c("setting", "IV.str", "VIO.str", "n", "Cov.mat", "Leng.mat", "CI.mat", "TSHT.VHat.mat", "Rule.mat" ,"length.union05", "cover.union05", "sim.round")))
filename = paste("Homo-Setting", setting, "-Strength", IV.str, "-Violation", VIO.str, "-n", n,"-SimRound",sim.round, ".RData", sep="")
save.image(filename)
