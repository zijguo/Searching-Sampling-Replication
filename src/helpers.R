grid.CI <- function(CI.matrix, grid.size){
  d = dim(CI.matrix)[1]
  grid.seq = NULL
  for(l in 1:d) grid.seq = c(grid.seq, seq(CI.matrix[l, 1], CI.matrix[l, 2], by=grid.size))
  return(grid.seq)
}

cutoff <- function(n, V.Gamma, V.gamma, C, InitiSet, beta.grid, N){
  pz = nrow(V.Gamma)
  n.beta = length(beta.grid)
  ## Covariance Matrix
  Cov1 = cbind(V.Gamma/n, C/n)
  Cov2 = cbind(t(C/n), V.gamma/n)
  Cov.total = rbind(Cov1, Cov2)
  
  Gen.mat = mvrnorm(N, rep(0, 2*pz), Cov.total)
  sample.sim = rep(0, N)
  for(i in 1:N){
    Gamma.i = Gen.mat[i, 1:pz]
    gamma.i = Gen.mat[i, -(1:pz)]
    temp.i = 0
    for(j in 1:n.beta){
      b = beta.grid[j]
      temp.diff = abs(Gamma.i[InitiSet]-b*gamma.i[InitiSet])/(sqrt(diag(V.Gamma + b^2*V.gamma - 2*b*C)/n)[InitiSet])
      temp.i = max(temp.i, temp.diff)
    }
    sample.sim[i] = temp.i
  }
  critical.val = quantile(sample.sim, probs=0.95)
  return(critical.val)
}

Searching.CI <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet, beta.grid, bootstrap=FALSE){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V.Gamma)[1]
  
  ## new rho method
  if(bootstrap){
    Tn = cutoff(n, V.Gamma, V.gamma, C, InitiSet, beta.grid, N=1000)
  }else{
    Tn = qnorm(1-0.05/(2*pz))
  }
  
  ## valid grid
  valid.grid = rep(NA, n.beta)
  for(j in 1:n.beta){
    b = beta.grid[j]
    temp = sqrt(diag(V.Gamma + b^2*V.gamma - 2*b*C)/n)
    valid.grid[j] = sum(abs(ITT_Y[InitiSet] - b*ITT_D[InitiSet]) < Tn*temp[InitiSet])
  }
  
  ## select beta
  if(length(beta.grid[which(valid.grid > threshold.size)])==0){
    rule=FALSE
    warning("Rule Fails. SS will give misleading CIs, SEs, and p-values.")
    sel.index = which(valid.grid==max(valid.grid))
  }else{
    rule=TRUE
    sel.index = which(valid.grid>threshold.size)
  }
  CI = t(as.matrix(c(min(beta.grid[sel.index]), max(beta.grid[sel.index]))))
  
  return(list(CI=CI, rule=rule))
}

Searching.CI.sampling <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet,
                                  beta.grid, rho=NULL, M=1000, prop=0.1, filtering=TRUE, bootstrap=FALSE){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V.Gamma)[1]
  
  ## new rho method
  if(bootstrap){
    Tn = cutoff(n, V.Gamma, V.gamma, C, InitiSet, beta.grid, N=1000)
  }else{
    Tn = qnorm(1-0.05/(2*pz))
  }
  
  ## Covariance Matrix
  Cov1 = cbind(V.Gamma/n, C/n)
  Cov2 = cbind(t(C/n), V.gamma/n)
  Cov.total = rbind(Cov1, Cov2)
  
  valid.grid.sample = matrix(NA, nrow=M, ncol=n.beta)
  Gen.mat = MASS::mvrnorm(M, rep(0, 2*pz), Cov.total)
  if(is.null(rho)) rho = (log(n)/M)^(1/(2*length(InitiSet)))/6 # initial rho if not specified
  
  if(filtering){
    temp1 = abs(t(t(Gen.mat[,InitiSet])/sqrt(diag(V.Gamma)[InitiSet]/n)))
    temp2 = abs(t(t(Gen.mat[,pz+InitiSet])/sqrt(diag(V.gamma)[InitiSet]/n)))
    temp = cbind(temp1, temp2)
    temp = apply(temp, MARGIN=1, FUN=max)
    temp = temp <= qnorm(1-0.05/(4*length(InitiSet)))
    Gen.mat = Gen.mat[temp,]
    M = sum(temp)
  }
  
  while(rho < 0.5){
    for(m in 1:M){
      ITT_Y.sample = ITT_Y - Gen.mat[m, 1:pz]
      ITT_D.sample = ITT_D - Gen.mat[m, (pz+1):(2*pz)]
      for(j in 1:n.beta){
        b = beta.grid[j]
        temp = sqrt(diag(V.Gamma + b^2*V.gamma - 2*b*C)/n)
        valid.grid.sample[m, j] = sum(abs(ITT_Y.sample[InitiSet] - b*ITT_D.sample[InitiSet])<
                                        rho*Tn*temp[InitiSet])
      }
    }
    CI = matrix(NA, nrow=M, ncol=2)
    for(m in 1:M){
      if(length(which(valid.grid.sample[m, ] > threshold.size))>0){
        CI[m, 1] = min(beta.grid[which(valid.grid.sample[m, ] > threshold.size)])
        CI[m, 2] = max(beta.grid[which(valid.grid.sample[m, ] > threshold.size)])
      }
    }
    CI = CI[!rowSums(is.na(CI)), , drop=FALSE] # CI = na.omit(CI)
    
    ## check CI dims, stop iterations if CI dim is big enough
    if(dim(as.matrix(CI))[1] >= prop*M) break
    rho = 1.25 * rho # increase rho with iterations
  }
  
  rule = TRUE
  if(dim(as.matrix(CI))[1] < prop*M){
    warning("Sampling Criterion not met, trasfer to Searching Method.")
    CI.searching = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet, beta.grid, bootstrap)
    rule = CI.searching$rule
    CI = CI.searching$CI
  }else{
    uni = Intervals(CI)
    CI = as.matrix(interval_union(uni))
    CI = t(as.matrix(c(min(CI[,1]), max(CI[,2]))))
  }
  
  return(list(CI=CI, rule=rule))
}

TSHT.Init <- function(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma, noise.mode="hetero"){
  n = nrow(WUMat); pz = ncol(WUMat)
  
  if(noise.mode=="homo"){
    ## compute Sigmas
    SigmaSqY = sum(resid_Y^2)/(n-1)
    SigmaSqD = sum(resid_D^2)/(n-1)
    SigmaYD = sum(resid_Y * resid_D)/(n-1)
  }
  
  ## First stage
  Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)/2))
  SHat = (1:pz)[abs(ITT_D) > (Tn * sqrt(diag(V.gamma)/n))]
  if(length(SHat)==0){
    warning("First Thresholding Warning: IVs individually weak. 
            TSHT with these IVs will give misleading CIs, SEs, and p-values. 
            Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE, pz); SHat.bool[SHat] = TRUE
  
  ## Second stage
  nCand = length(SHat)
  VHats.bool = matrix(FALSE, nCand, nCand)
  colnames(VHats.bool) = rownames(VHats.bool) = SHat
  
  for(j in SHat){
    beta.j = ITT_Y[j]/ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    if(noise.mode=="homo"){
      sigmasq.j = SigmaSqY + beta.j^2 * SigmaSqD - 2* beta.j * SigmaYD
      PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j * colSums((WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) * 
        sqrt(tuning^2 * log(pz)/n)  ### tuning
    }else if(noise.mode=="hetero"){
      mid.j = resid_Y - beta.j*resid_D
      side.j = WUMat - outer(WUMat[, j]/ITT_D[j], ITT_D)
      PHat.bool.j = abs(pi.j) <= sqrt(colSums((mid.j * side.j)^2)/n^2) * sqrt(log(n))
    }
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }
  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }
  diag(VHats.boot.sym) = 1
  
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
  VHat=sort(as.numeric(VHat))
  
  out <- list(SHat=SHat, VHat=VHat, voting.mat=VHats.boot.sym)
  return(out)
}

TSHT.Init.small <- function(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma, noise.mode="hetero"){
  n = nrow(WUMat); pz = ncol(WUMat)
  
  if(noise.mode=="homo"){
    ## compute Sigmas
    SigmaSqY = sum(resid_Y^2)/(n-1)
    SigmaSqD = sum(resid_D^2)/(n-1)
    SigmaYD = sum(resid_Y * resid_D)/(n-1)
  }
  
  ## First stage
  Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)/2))
  SHat = (1:pz)[abs(ITT_D) > (Tn * sqrt(diag(V.gamma)/n))]
  if(length(SHat)==0){
    warning("First Thresholding Warning: IVs individually weak. 
            TSHT with these IVs will give misleading CIs, SEs, and p-values. 
            Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE, pz); SHat.bool[SHat] = TRUE
  
  ## Second stage
  nCand = length(SHat)
  VHats.bool = matrix(FALSE, nCand, nCand)
  colnames(VHats.bool) = rownames(VHats.bool) = SHat
  
  for(j in SHat){
    beta.j = ITT_Y[j]/ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    if(noise.mode=="homo"){
      sigmasq.j = SigmaSqY + beta.j^2 * SigmaSqD - 2* beta.j * SigmaYD
      PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j * colSums((WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) * 
        sqrt(tuning^2 * log(pz)/n)  ### tuning
    }else if(noise.mode=="hetero"){
      mid.j = resid_Y - beta.j*resid_D
      side.j = WUMat - outer(WUMat[, j]/ITT_D[j], ITT_D)
      PHat.bool.j = abs(pi.j) <= sqrt(colSums((mid.j * side.j)^2)/n^2) * sqrt(log(n))
    }
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }
  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }
  diag(VHats.boot.sym) = 1
  
  # Voting
  VM= apply(VHats.boot.sym,1,sum)
  VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
  VM.p = rownames(VHats.boot.sym)[max(VM) == VM] #Plurality winners
  VHat = as.numeric(union(VM.m,VM.p))
  
  # Error check
  if(length(VHat) == 0){
    warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being valid")
    VHat = 1:pz
  }
  
  return(list(VHat = VHat,SHat=SHat))
}

## Change threshold for TSHT.Init
TSHT.Init_1 <- function(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma, noise.mode="hetero"){
  n = nrow(WUMat); pz = ncol(WUMat)
  
  if(noise.mode=="homo"){
    ## compute Sigmas
    SigmaSqY = sum(resid_Y^2)/(n-1)
    SigmaSqD = sum(resid_D^2)/(n-1)
    SigmaYD = sum(resid_Y * resid_D)/(n-1)
  }
  
  ## First stage
  Tn = sqrt(log(n))
  SHat = (1:pz)[abs(ITT_D) > (Tn * sqrt(diag(V.gamma)/n))]
  if(length(SHat)==0){
    warning("First Thresholding Warning: IVs individually weak. 
            TSHT with these IVs will give misleading CIs, SEs, and p-values. 
            Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE, pz); SHat.bool[SHat] = TRUE
  
  ## Second stage
  nCand = length(SHat)
  VHats.bool = matrix(FALSE, nCand, nCand)
  colnames(VHats.bool) = rownames(VHats.bool) = SHat
  
  for(j in SHat){
    beta.j = ITT_Y[j]/ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    if(noise.mode=="homo"){
      sigmasq.j = SigmaSqY + beta.j^2 * SigmaSqD - 2* beta.j * SigmaYD
      PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j * colSums((WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) * 
        sqrt(tuning^2 * log(pz)/n)  ### tuning
    }else if(noise.mode=="hetero"){
      mid.j = resid_Y - beta.j*resid_D
      side.j = WUMat - outer(WUMat[, j]/ITT_D[j], ITT_D)
      PHat.bool.j = abs(pi.j) <= sqrt(colSums((mid.j * side.j)^2)/n^2) * sqrt(log(n))
    }
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }
  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }
  diag(VHats.boot.sym) = 1
  
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
  VHat=sort(as.numeric(VHat))
  
  out <- list(SHat=SHat, VHat=VHat, voting.mat=VHats.boot.sym)
  return(out)
}
TSHT.Init_2 <- function(ITT_Y, ITT_D, resid_Y, resid_D, WUMat, V.gamma, noise.mode="hetero"){
  n = nrow(WUMat); pz = ncol(WUMat)
  
  if(noise.mode=="homo"){
    ## compute Sigmas
    SigmaSqY = sum(resid_Y^2)/(n-1)
    SigmaSqD = sum(resid_D^2)/(n-1)
    SigmaYD = sum(resid_Y * resid_D)/(n-1)
  }
  
  ## First stage
  Tn = sqrt(2*log(n))
  SHat = (1:pz)[abs(ITT_D) > (Tn * sqrt(diag(V.gamma)/n))]
  if(length(SHat)==0){
    warning("First Thresholding Warning: IVs individually weak. 
            TSHT with these IVs will give misleading CIs, SEs, and p-values. 
            Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE, pz); SHat.bool[SHat] = TRUE
  
  ## Second stage
  nCand = length(SHat)
  VHats.bool = matrix(FALSE, nCand, nCand)
  colnames(VHats.bool) = rownames(VHats.bool) = SHat
  
  for(j in SHat){
    beta.j = ITT_Y[j]/ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    if(noise.mode=="homo"){
      sigmasq.j = SigmaSqY + beta.j^2 * SigmaSqD - 2* beta.j * SigmaYD
      PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j * colSums((WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) * 
        sqrt(tuning^2 * log(pz)/n)  ### tuning
    }else if(noise.mode=="hetero"){
      mid.j = resid_Y - beta.j*resid_D
      side.j = WUMat - outer(WUMat[, j]/ITT_D[j], ITT_D)
      PHat.bool.j = abs(pi.j) <= sqrt(colSums((mid.j * side.j)^2)/n^2) * sqrt(log(n))
    }
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }
  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }
  diag(VHats.boot.sym) = 1
  
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
  VHat=sort(as.numeric(VHat))
  
  out <- list(SHat=SHat, VHat=VHat, voting.mat=VHats.boot.sym)
  return(out)
}

Lasso.init <- function(X, y, lambda = "CV.min", intercept = FALSE) {
  p <- ncol(X)
  n <- nrow(X)
  
  htheta <- if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  }
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}

tryCatch_E <- function(expr, ret.obj){
  # error handler
  errhandler <- function(e) {
    err <<- conditionMessage(e)
    ret.obj # Return argument ret.obj if an error occcurs.
  }
  
  # evaluate the expression
  err <- NULL
  value <- withCallingHandlers(tryCatch(expr, error = errhandler))
  
  # return a list with value, error and warning
  list(value = value, error = err)
}

## Rewritten functions from the source code of CIIV
robust_method <- function(Y, regressor, Exogenous){
  length_regressor = ncol(regressor)
  Coefficients_CIM_GMM <- (CIM.HansenJTest(Y, regressor, Exogenous)[[2]])[1:length_regressor]
  sd_CIM_GMM <- (CIM.HansenJTest(Y, regressor, Exogenous)[[3]])[1:length_regressor];
  alpha = 0.05
  ci_CIM_GMM <- c(
    Coefficients_CIM_GMM[1] - qnorm(1-alpha/2) * sd_CIM_GMM[1],
    Coefficients_CIM_GMM[1] + qnorm(1-alpha/2) * sd_CIM_GMM[1]
  )
  return(ci_CIM_GMM)
}

## Functions from the source code of CIIV
# Function for Two Step GMM and Hansen-J Test
CIM.HansenJTest <- function(Y,X,Z){
  res_FirstStep <- residuals(AER::ivreg(Y ~ X - 1 | Z));
  
  Weight_SecondStep <- crossprod(res_FirstStep * Z);
  
  Coef_SecondStep <- solve(
    t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
  ) %*% t(X) %*% Z %*% solve(Weight_SecondStep) %*% t(Z) %*% Y;
  
  res_SecondStep <- as.vector(Y - X %*% Coef_SecondStep);
  
  sd_SecondStep <- sqrt(diag(solve(
    t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
  ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)%*%crossprod(res_SecondStep * Z)%*%t(
    solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)
  )));
  
  HansenJ_Stat <- t(res_SecondStep) %*% Z %*% solve(Weight_SecondStep) %*%
    t(Z) %*% res_SecondStep;
  
  list(HansenJ_Stat, Coef_SecondStep,sd_SecondStep)
}
