#### Date: 4/12/2022
#### Author: Zhenyu WANG
#### Aim: CIIV plurality settings in the paper
#### Date: 05/18/2022
#### Aim: add CI output
library(MASS)
library(intervals)
source("helpers.R")
source("TSHT-ldim.R")
source("invalidIV.R")
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")

########### Part-1: Specify Settings and Parameters ##########
# parameters changing
n=2000  # {500, 1000, 2000, 5000}
VIO.str = 0.2 # {0.2, 0.4}
setting = "CIIV1" # {"CIIV1", "CIIV2"}
sim.round = 1 # {1,2,3,4,5,6,7,8,9,10}
# the number of simulations
nsim=50

IV.str = 0.4
if(setting=="CIIV1"){
  L = 21
  s1 = s2 = 6
  s = s1+s2
  ZSigma = matrix(0,L,L)
  for(j in 1:L){
    for(k in 1:L){
      ZSigma[j, k] = 0.5^abs(j-k)
      ZSigma[k, j] = 0.5^abs(k-j)
    }
  }
  alpha = VIO.str*c(rep(0, L-s1-s2), rep(1, s1), rep(0.5, s2))
}else if(setting=="CIIV2"){
  L = 21
  s1 = s2 = 6
  s = s1+s2
  ZSigma = matrix(0,L,L)
  for(j in 1:L){
    for(k in 1:L){
      ZSigma[j, k] = 0.5^abs(j-k)
      ZSigma[k, j] = 0.5^abs(k-j)
    }
  }
  alpha = VIO.str*c(rep(0, L-s1-s2), rep(1, s1/2), rep(-1, s1/2), rep(0.5, s2/2), rep(-0.5, s2/2))
}
beta = 1
gamma = rep(IV.str, L)

########### Part-2: Simulations ##########
# matrix storing info
CI.mat = matrix(NA, nrow=nsim, ncol=6)
colnames(CI.mat) = c("Oracle-L","Oracle-U","TSHT-L","TSHT-U","CIIV-L","CIIV-U")
Cov.mat = Leng.mat = matrix(NA, nrow=nsim, ncol=9)
colnames(Cov.mat) = colnames(Leng.mat) = c("oracle", "TSHT", "CIIV", "Sear-TSHT", "Sear-CIIV", "Sear-Comb",
                                           "Samp-TSHT", "Samp-CIIV", "Samp-Comb")
CIIV.valid_rule = rep(NA, nsim)
Rule.mat = matrix(NA, nrow=nsim, ncol=4)
colnames(Rule.mat) = c("Sear-TSHT", "Sear-CIIV", "Samp-TSHT", "Samp-CIIV")
comp.union=TRUE
comp.union.median=FALSE
method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR")
length.union05<-matrix(NA, nrow=nsim, ncol=10)
cover.union05<-matrix(NA, nrow=nsim, ncol=10)
colnames(length.union05)=c(method,method)
colnames(cover.union05)=c(method,method)

for(i.sim in 1:nsim){
  set.seed(i.sim+(sim.round-1)*nsim)
  print(i.sim)
  
  epsilonSigma = matrix(c(1,0.25,0.25,1),2,2)
  Z = mvrnorm(n, rep(0,L), ZSigma)
  epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
  D = 0.5 + Z %*% gamma + epsilon[,1]
  Y = -0.5 + Z %*% alpha + D * beta + epsilon[,2]
  
  W = Z
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
  CI = confint(ivreg(Y ~ D + Z[,-(1:(L-s))] | Z))[2,]
  CI.mat[i.sim, 1:2] = CI
  Cov.mat[i.sim, 1] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 1] = CI[2] - CI[1]

  # TSHT Method
  TSHT.temp = TSHT(Y, D, Z)
  CI = TSHT.temp$ci
  CI.mat[i.sim, 3:4] = CI
  Cov.mat[i.sim, 2] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 2] = CI[2] - CI[1]

  # CIIV Method
  CIIV.temp = tryCatch_E(CIIV(Y,D,Z,robust=FALSE), ret.obj = list(NULL))
  if (is.null(CIIV.temp$error)) CIIV.result <- CIIV.temp$value
  CI.mat[i.sim, 5:6] = CIIV.result$ci_CIM
  Cov.mat[i.sim,3] <- ifelse(beta<CIIV.result$ci_CIM[2]&beta>CIIV.result$ci_CIM[1], 1, 0)
  Leng.mat[i.sim, 3] <- CIIV.result$ci_CIM[2] - CIIV.result$ci_CIM[1]

  ######## Select VHat ########
  ## Use TSHT method to select
  TSHT.out <- TSHT.Init(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma)
  VHat.TSHT = sort(TSHT.out$VHat)
  
  ## Use CIIV method to select
  VHat.CIIV = sapply(CIIV.result$Valid_Instruments, function(x) as.numeric(substring(x, 2)))
  CIIV.valid_rule[i.sim] = all(VHat.CIIV %in% (1:pz))
  if(!all(VHat.CIIV %in% (1:pz))){
    VHat.CIIV = VHat.TSHT
  }
  
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
    resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
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
    resultOut.temp = tryCatch_E(invalidIVCI(Y, D, Z, U = s_bar, alpha = 0.05, alpha.overid = 0.05/2, intercept=TRUE,convexHull=FALSE,method=method),
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

rm(list=setdiff(ls(), c("setting", "IV.str", "VIO.str", "n", "Cov.mat", "Leng.mat", "Rule.mat" ,"length.union05", "cover.union05", "CIIV.valid_rule", "CI.mat", "sim.round")))
filename = paste("CIIV-Setting", setting, "-Strength", IV.str, "-Violation", VIO.str, "-n", n,"-SimRound",sim.round, ".RData", sep="")
save.image(filename)