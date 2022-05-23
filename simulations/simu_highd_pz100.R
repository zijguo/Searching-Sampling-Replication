###########################
### Date: 05/02/2022
### Author: Zhenyu Wang
### Aim: High Dimension Simulations
library(MASS)
library(intervals)
library(glmnet)
library(flare)
library(CVXR)
source("src/helpers.R")
source("src/SIHR-hdim.R")
source("src/TSHT-ldim.R")
source("src/invalidIV.R")
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")
########### Part-1: Specify Settings and Parameters ##########
# parameters changing
n=1000  # {200, 300, 1000, 2500}
VIO.str =0.5 # {0.5, 1, 2}
sim.round = 20 # seq(1,25)
nsim = 25 # the number of simulations

####################
case = "homo"
IV.str = 0.5
pi.value = IV.str*VIO.str
beta = 1
L = 100
# alpha = c(rep(0,4), rep(pi.value, 2), -seq(1,4)/3, rep(0, L-10))
# gamma = c(rep(IV.str, 10), rep(0, L-10))
alpha = c(rep(0, 5), rep(pi.value, 2), rep(0, L-7))
gamma = c(rep(IV.str, 7), rep(0, L-7))

px = 150
p=L+px # p stands for the number of total exogeneous variables 
phi<-rep(0,px)
psi<-rep(0,px)
phi[1:10] <- 1/10*seq(1, 10)+0.5
psi[1:10] <- 1/10*seq(1, 10)+1

rho=0.5
A1gen <- function(rho, p){
  A1 = matrix(0, nrow=p, ncol=p)
  for(i in 1:p) for(j in 1:p) A1[i, j] = rho^(abs(i-j))
  return(A1)
}
Cov<-(A1gen(rho,p))

########### Part-2: Simulations ##########

# matrix storing info
Cov.mat = Leng.mat = CI.mat = matrix(NA, nrow=nsim, ncol=4)
colnames(Cov.mat) = colnames(Leng.mat) = c("oracle", "TSHT", "Sear", "Samp")
colnames(CI.mat) = c("oracle-L","oracle-U","TSHT-L","TSHT-U")
Rule.mat = matrix(NA, nrow=nsim, ncol=2)
colnames(Rule.mat) = c("Sear","Samp")

start_time = Sys.time()
for(i.sim in 1:nsim){
  set.seed(i.sim+(sim.round-1)*nsim)
  print(i.sim)
  W = mvrnorm(n, rep(0, p), Cov)
  Z = W[, 1:L]
  X = W[, (L+1):p]
  if(case=="homo"){
    # epsilonSigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)
    epsilonSigma = matrix(c(1.5, 0.75, 0.75, 1.5), 2, 2)
    epsilon = mvrnorm(n, rep(0, 2), epsilonSigma)
    epsilon1 = epsilon[,1]
    epsilon2 = epsilon[,2]
  }
  D = 0.5 + Z %*% gamma+ X%*% psi + epsilon1
  Y = -0.5 + Z %*% alpha + D * beta + X%*%phi+ epsilon2
  
  # ######### Preparation ###########
  if(is.null(X)) W = Z else W = cbind(Z, X)
  n = length(Y); pz = ncol(Z); p = ncol(W)
  intercept = TRUE
  if(intercept) W = cbind(W, 1)

  ## LF estimators
  init_Y = Lasso.init(W, Y)
  init_D = Lasso.init(W, D)

  ## residual
  resid_Y = as.vector(Y - W%*%init_Y)
  resid_D = as.vector(D - W%*%init_D)

  ## Debias
  ITT_Y = rep(NA, pz); ITT_D = rep(NA, pz)
  U = matrix(NA, nrow=ncol(W), ncol=pz)
  W_no_intercept = W[,-ncol(W)]
  for(i in 1:pz){
    print(paste("Debias--",i,sep=""))
    loading = rep(0, p)
    loading[i] = 1
    ITT_Y[i] = GLM_LF(W_no_intercept, Y, loading=loading, intercept.loading=FALSE, intercept=TRUE, init.coef = NULL, verbose=FALSE)$prop.est
    model_D= GLM_LF(W_no_intercept, D, loading=loading, intercept.loading=FALSE, intercept=TRUE, init.coef = NULL, verbose=FALSE)
    ITT_D[i] = model_D$prop.est
    U[-nrow(U),i] = (model_D$proj)[-1]
    U[nrow(U),i] = (model_D$proj)[1]
  }
  WUMat = W%*%U

  SigmaSqY = sum(resid_Y^2)/(n-1)
  SigmaSqD = sum(resid_D^2)/(n-1)
  SigmaYD = sum(resid_Y * resid_D)/(n-1)

  Temp = t(WUMat)%*%WUMat / n
  V.Gamma = SigmaSqY * Temp
  V.gamma = SigmaSqD * Temp
  C = SigmaYD * Temp

  ######## Column Group-1 #########
  # Oracle Method
  CI = confint(ivreg(Y ~ D + Z[,5:10] + X[,1:10] | Z[,1:10] + X[,1:10]))[2,]
  CI.mat[i.sim, 1:2] = CI
  Cov.mat[i.sim, 1] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 1] = CI[2] - CI[1]

  # TSHT Method
  TSHT.temp = TSHT(Y, D, Z, X, method="DeLasso")
  CI = TSHT.temp$ci
  CI.mat[i.sim, 3:4] = CI
  Cov.mat[i.sim, 2] = (CI[1]<beta)*(CI[2]>beta)
  Leng.mat[i.sim, 2] = CI[2] - CI[1]
  
  ######## Select VHat ########
  ## Use TSHT method to select
  TSHT.out <- TSHT.Init(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma)
  VHat.TSHT = sort(TSHT.out$VHat)

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
  # TSHT Searching
  out = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = VHat.TSHT,
                     beta.grid = beta.grid.TSHT)
  CI1 = out$CI; rule1 = out$rule
  Cov.mat[i.sim, 3] = sum((CI1[,1]<beta)*(CI1[,2]>beta))
  Leng.mat[i.sim, 3] = sum(CI1[,2] - CI1[,1])
  # Rule
  Rule.mat[i.sim, 1] = rule1

  ####### Column Group-3 #########
  # TSHT Sampling
  out = Searching.CI.sampling(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet=VHat.TSHT,
                              beta.grid = beta.grid.TSHT)
  CI1 = out$CI; rule1 = out$rule
  Cov.mat[i.sim, 4] = sum((CI1[,1]<beta)*(CI1[,2]>beta))
  Leng.mat[i.sim, 4] = sum(CI1[,2] - CI1[,1])
  # Rule
  Rule.mat[i.sim, 2] = rule1
}
end_time = Sys.time()
time_diff = end_time - start_time

rm(list=setdiff(ls(), c("n", "IV.str", "VIO.str", "Cov.mat", "Leng.mat", "CI.mat", "time_diff", "sim.round")))
filename <- paste("Highd-n",n,"-Violation",VIO.str,"-SimRound",sim.round,".RData",sep="")
save.image(filename)