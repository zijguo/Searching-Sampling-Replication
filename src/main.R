##################################
### Date: 4/21/2022
### Author: Zhenyu Wang
### Aim: 
### 1.revise functions for real datasets (No MaxClique, + intermediate outputs)

## This version only Plurality Algorithm without CIIV option, but integrates MaxClique option
SearchingSampling <- function(Y, D, Z, X, intercept=TRUE, 
                              method=c("OLS","LF"), 
                              noise.mode=c("hetero","homo"),
                              CI.init = NULL,
                              a=0.6,
                              Sampling=TRUE, 
                              rho=NULL, M=1000, prop=0.1){
  ###############################################
  ## Arguments:
  ## CI.init   initial range [L, U], if not specified we provided a default method to construct
  ## Sampling  Do Sampling method or not?
  ## rho       For sampling method, the initial rho (corresponding to lambda in paper)
  ## M         For sampling method, sampling times
  method = match.arg(method)
  noise.mode = match.arg(noise.mode)

  if(is.null(X)) W = Z else W = cbind(Z, X)
  n = length(Y); pz = ncol(Z); p = ncol(W)
  
  ## Preparation
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
  
  ## Computing (co)variance matrices
  if(noise.mode=="homo"){
    SigmaSqY = sum(resid_Y^2)/(n-1)
    SigmaSqD = sum(resid_D^2)/(n-1)
    SigmaYD = sum(resid_Y * resid_D)/(n-1)
    V.Gamma = SigmaSqY * U[1:pz, 1:pz]
    V.gamma = SigmaSqD * U[1:pz, 1:pz]
    C = SigmaYD * U[1:pz, 1:pz]
  }else{
    V.Gamma = (t(WUMat)%*%diag(resid_Y^2)%*%WUMat)/n
    V.gamma = (t(WUMat)%*%diag(resid_D^2)%*%WUMat)/n
    C = (t(WUMat)%*%diag(resid_Y * resid_D)%*%WUMat)/n
  }
  
  TSHT.out <- TSHT.Init(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma, noise.mode=noise.mode)
  V0.hat = sort(TSHT.out$VHat)
  ## Construct range [L, U]
  if(is.vector(CI.init)){
    CI.init.union = t(as.matrix(sort(CI.init)))
  }else{
    ## current method to select initial [L, U]
    var.beta = 1/n * (diag(V.Gamma)/ITT_D^2 + diag(V.gamma)*ITT_Y^2/ITT_D^4 - 2*diag(C)*ITT_Y/ITT_D^3)
    var.beta = var.beta[V0.hat]
    CI.init = matrix(NA, nrow=length(V0.hat), ncol=2)
    CI.init[,1] = (ITT_Y/ITT_D)[V0.hat] - sqrt(log(n)*var.beta)
    CI.init[,2] = (ITT_Y/ITT_D)[V0.hat] + sqrt(log(n)*var.beta)
    uni = Intervals(CI.init)
    CI.init.union = as.matrix(interval_union(uni))
  }
  
  # Construct beta.grid
  beta.grid = grid.CI(CI.init.union, grid.size=n^{-a})
  
  if(Sampling){
    ## Sampling Method
    CI.sampling = Searching.CI.sampling(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet=V0.hat,
                                        beta.grid = beta.grid, rho=rho, M=M, prop=prop, filtering=TRUE)
    CI=CI.sampling$CI
    rule=CI.sampling$rule
  }else{
    ## Searching Method
    CI.searching = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = V0.hat,
                                beta.grid = beta.grid)
    CI=CI.searching$CI
    rule=CI.searching$rule
  }
  returnList <- list(CI=CI, rule=rule, VHat=V0.hat, CI.init=CI.init.union, TSHT.out = TSHT.out)
  
  return(returnList)
}

