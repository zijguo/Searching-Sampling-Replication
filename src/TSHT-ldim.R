### TSHT.R
### Function: Implements the two-stage hard thresholing
###           procedure for estimation and inference of 
###           treatment effects in the presence of invalid IVs
###           as developed in Guo, Kang, Cai, and Small (2016)             
### Maintainer: Hyunseung Kang
### E-mail: hyunseung@stat.wisc.edu

### TSHT
### FUNCTION: Point estimate, SE, and CI for treatment effect with invalid IVs 
###           using a single-sample, individual-level data via TSHT. 
### INPUT: Y, continuous, non-missing, numeric outcome vector (u by 1 vector)
###        D, continuous or discrete, non-missing, numeric treatment vector (n by 1 vector)
###        Z, continuous or discrete, non-missing, numeric instrument matrix containing p_z 
###           instruments (n by p_z matrix)
###        X, optional continuous or discrete, non-missing, numeric matrix containing p_x 
###           covariates (n by p_x matrix)
###        intercept, a boolean scalar asking "should the intercept term be included?" 
###                   with TRUE/FALSE (default = TRUE)
###        alpha, a numeric scalar value between 0 and 1 indicating the significance level for 
###               the confidence interval (default = 0.05)
###        tuning, a numeric scalar value tuning parameter for TSHT greater 
###                than 2 (default = 2.01)
###        method, a character scalar declaring the method used to estimate the inputs in TSHT
###                (default = "OLS")
### OUTPUT: a list (a) VHat (numeric vector denoting the set of valid and relevant IVs)
###                (b) SHat (numeric vector denoting the set of relevant IVs)
###                (c) betaHat (scalar numeric value:
###                             the estimate of the treatment effect)
###                (d) varHat (scalar numeric value:
###                            estimated variance of betaHat)
###                (e) ci (two dimensional vector, 
###                           the 1-alpha confidence interval for beta
###                           with the lower and upper endpoints)
TSHT <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05,tuning=2.01,method=c("OLS","DeLasso")) {
  method = match.arg(method)
  # Check and Clean Input Type #
  # Check Y
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  Y = as.numeric(Y)
  
  # Check D
  stopifnot(!missing(D),(is.numeric(D) || is.logical(D)),(is.matrix(D) || is.data.frame(D)) && ncol(D) == 1)
  stopifnot(all(!is.na(D)))
  D = as.numeric(D)
  
  # Check Z
  stopifnot(!missing(Z),(is.numeric(Z) || is.logical(Z)),is.matrix(Z))
  stopifnot(all(!is.na(Z)))
  
  # Check dimesions
  stopifnot(length(Y) == length(D), length(Y) == nrow(Z))
  
  # Check X, if present
  if(!missing(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),is.matrix(X) && nrow(X) == nrow(Z))
    stopifnot(all(!is.na(X)))
    
    W = cbind(Z,X)
  } else {
    W = Z
  }
  
  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.numeric(alpha),length(alpha) == 1,alpha <= 1,alpha >= 0)
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)
  stopifnot(is.character(method))
  
  # Derive Inputs for TSHT
  n = length(Y); pz=ncol(Z)
  if(method == "OLS") {
    inputs = TSHT.OLS(Y,D,W,pz,intercept)
    A = t(inputs$WUMat) %*% inputs$WUMat / n
    # A = diag(pz)
  } else{
    inputs = TSHT.DeLasso(Y,D,W,pz,intercept)
    A = diag(pz)
  }
  
  Temp = t(inputs$WUMat)%*%inputs$WUMat / n
  V.Gamma = inputs$SigmaSqY * Temp
  V.gamma = inputs$SigmaSqD * Temp
  C = inputs$SigmaYD * Temp
  
  # Estimate Valid IVs
  SetHats = TSHT.VHat(ITT_Y = inputs$ITT_Y,ITT_D = inputs$ITT_D,WUMat = inputs$WUMat,
                      SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,SigmaYD=inputs$SigmaYD,tuning=tuning)
  VHat = SetHats$VHat; SHat = SetHats$SHat
  
  # Obtain point est, se, and ci
  AVHat = solve(A[VHat,VHat])
  betaHat = (t(inputs$ITT_Y[VHat]) %*% AVHat %*% inputs$ITT_D[VHat]) / (t(inputs$ITT_D[VHat]) %*% AVHat %*% inputs$ITT_D[VHat])
  SigmaSq = inputs$SigmaSqY + betaHat^2 * inputs$SigmaSqD - 2*betaHat * inputs$SigmaYD
  betaVarHat = SigmaSq * (t(inputs$ITT_D[VHat]) %*% AVHat %*% (t(inputs$WUMat) %*% inputs$WUMat/ n)[VHat,VHat] %*% AVHat %*% inputs$ITT_D[VHat]) / (t(inputs$ITT_D[VHat]) %*% AVHat %*% inputs$ITT_D[VHat])^2
  ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVarHat / n),betaHat + qnorm(1-alpha/2) * sqrt(betaVarHat/n))
  
  return(list(VHat = VHat,SHat=SHat,betaHat=betaHat,betaVarHat = betaVarHat,ci=ci,
              V.Gamma=V.Gamma, V.gamma=V.gamma, C=C))
}

TSHT.OLS <- function(Y,D,W,pz,intercept=TRUE) {
  n = nrow(W)
  # Include intercept
  if(intercept) {
    W = cbind(W,1)
  }
  p = ncol(W) 
  
  # Compute covariance of W and W %*% U
  covW = t(W) %*% W /n #this should automatically turn covW into a matrix
  WUMat = W %*% (solve(covW))[,1:pz]
  
  # First Part (OLS Estimation)
  qrW = qr(W)
  ITT_Y = qr.coef(qrW,Y)[1:pz]
  ITT_D = qr.coef(qrW,D)[1:pz]
  SigmaSqY = sum(qr.resid(qrW,Y)^2)/(n -p)
  SigmaSqD = sum(qr.resid(qrW,D)^2)/(n -p)
  SigmaYD = sum(qr.resid(qrW,Y) * qr.resid(qrW,D)) / (n - p)
  
  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD))
}

TSHT.DeLasso <- function(Y,D,W,pz,intercept=TRUE) {
  # Fit Reduced-Form Model for Y and D
  model_Y <- SSLasso(X=W,y=Y,intercept=intercept,verbose=FALSE)
  model_D = SSLasso(X=W,y=D,intercept=intercept,verbose=FALSE) 
  ITT_Y = model_Y$unb.coef[1:(pz)]
  ITT_D = model_D$unb.coef[1:(pz)]
  resid_Y = model_Y$resid.lasso; resid_D = model_D$resid.lasso
  SigmaSqY=sum(resid_Y^2)/n
  SigmaSqD=sum(resid_D^2)/n
  SigmaYD =sum(resid_Y * resid_D)/n
  WUMat = model_D$WUMat[,1:pz]
  
  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD))
}


### TSHT.VHat
### FUNCTION: Estimates the set of valid and relevant IVs
### INPUT: ITT_Y: a numeric, non-missing vector estimating each instrument's effect on the outcome
###        ITT_D: a numeric, non-missing vector estimating each instrument's effect on the treatment
###        WUMat: an n by pz numeric, non-missing matrix that computes the precision of W 
###               (ex1 t(WUMat) %*% WUMat / n = (W^T W/n)^(-1))
###               (ex2 U = (W^T W/n)^(-1))
###        SigmaSqY: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of Y in the reduced-form model of Y
###        SigmaSqD: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of D in the reduced-form model of D
###        SigmaYD: a numeric, non-missing, scalar value estimating the covariance
###                 of the error terms in the reduced-form models of Y and D
###        tuning, a numeric scalar value tuning parameter for TSHT greater 
###                than 2 (default = 2.01)
### OUTPUT: a list (a) VHat (numeric vector denoting the set of valid and relevant IVs)
###                (b) SHat (numeric vector denoting the set of relevant IVs)
TSHT.VHat <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,tuning = 2.01) {
  # Check ITT_Y and ITT_D
  stopifnot(!missing(ITT_Y),!missing(ITT_D),length(ITT_Y) == length(ITT_D))
  stopifnot(all(!is.na(ITT_Y)),all(!is.na(ITT_D)))
  ITT_Y = as.numeric(ITT_Y); ITT_D = as.numeric(ITT_D)
  
  # Check WUMat
  stopifnot(!missing(WUMat),is.matrix(WUMat), nrow(WUMat) > 1, ncol(WUMat) == length(ITT_Y)) 
  stopifnot(all(!is.na(WUMat)))
  
  # Check Sigmas 
  stopifnot(!missing(SigmaSqY), is.numeric(SigmaSqY), length(SigmaSqY) == 1, !is.na(SigmaSqY), SigmaSqY > 0)
  stopifnot(!missing(SigmaSqD), is.numeric(SigmaSqD), length(SigmaSqD) == 1, !is.na(SigmaSqD),SigmaSqD > 0)
  stopifnot(!missing(SigmaYD),is.numeric(SigmaYD), length(SigmaYD) == 1,!is.na(SigmaYD))
  
  # Other Input check
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)
  
  # Constants
  n = nrow(WUMat); pz = length(ITT_Y)
  
  # First Stage
  SHat = (1:pz)[(abs(ITT_D) >= (sqrt(SigmaSqD * colSums(WUMat^2) /n) * sqrt(tuning*log(pz)/n)))]
  if(length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. TSHT with these IVs will give misleading CIs, SEs, and p-values. Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE
  
  # Second Stage
  # pi.candidate is the estimated value of pi across different candidates
  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat
  for(j in SHat) {
    beta.j = ITT_Y[j] / ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    sigmasq.j = SigmaSqY + beta.j^2 * SigmaSqD - 2* beta.j * SigmaYD
    PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j * colSums( (WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) * 
      sqrt(tuning^2 * log(pz)/n)
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
SSLasso <- function (X, y, lambda = NULL, mu = NULL, intercept = TRUE,
                     resol=1.3, maxiter=50, threshold=1e-2, verbose = TRUE) {
  #
  # Compute confidence intervals and p-values.
  #
  # Args:
  #   X     :  design matrix
  #   y     :  response
  #   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
  #   mu    :  Linfty constraint on U (if null, searches)
  #   intercept: Should the intercept term be included?
  #   resol :  step parameter for the function that computes U
  #   maxiter: iteration parameter for computing U
  #   threshold : tolerance criterion for computing U
  #   verbose : verbose?
  #
  # Returns:
  #   coef    : Lasso estimated coefficients
  #   unb.coef: Unbiased coefficient estimates
  #   WUMat: projection of the inverse covariance matrix.
  #   resid.lasso: residual based on Lasso
  #
  p <- ncol(X);
  n <- nrow(X);
  pp <- p;
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
  X <- X %*% diag(col.norm);
  
  # Solve Lasso problem using FLARE package
  htheta <- Lasso (X,y,lambda=lambda,intercept=intercept);
  
  # Format design matrix to include intercept and standardize.
  if (intercept==TRUE){
    Xb <- cbind(rep(1,n),X);
    col.norm <- c(1,col.norm);
    pp <- (p+1);
  } else {
    Xb <- X;
  }
  resid.lasso = (y - Xb %*% htheta)
  sigma.hat <- (1/n)*(t(Xb)%*%Xb);
  
  # Estimation of U (or M in Javanard and Montanari's original paper)
  # Check to see if this is a low dimensional problem
  if ((n>=2*p)){
    tmp <- eigen(sigma.hat)
    tmp <- min(tmp$values)/max(tmp$values)
  }else{
    tmp <- 0
  }
  
  # If low-dimensional problem, use inverse of covariance as an estiamte of precision matrix
  # Otherwise, solve the convex optimizatio problem for U
  if ((n>=2*p)&&(tmp>=1e-4)){
    U <- solve(sigma.hat)
  }else{
    U <- InverseLinfty(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose);
  }
  
  # Debias Lasso
  unbiased.Lasso <- as.numeric(htheta + (U%*%t(Xb)%*%(y - Xb %*% htheta))/n);
  
  # Scale them back to the original scaling.
  htheta <- htheta*col.norm;
  unbiased.Lasso <- unbiased.Lasso*col.norm;
  WUMat = Xb %*% t(U) %*% diag(col.norm)
  
  if (intercept==TRUE){
    htheta <- htheta[2:pp];
    unbiased.Lasso <- unbiased.Lasso[2:pp];
    WUMat <- WUMat[,2:pp]
  }
  
  returnList <- list("coef" = htheta,
                     "unb.coef" = unbiased.Lasso,
                     "WUMat" = WUMat,
                     "resid.lasso" = resid.lasso)
  return(returnList)
}

Lasso <- function( X, y, lambda = NULL, intercept = TRUE){
  #
  # Compute the Lasso estimator:
  # - If lambda is given, use glmnet and standard Lasso
  # - If lambda is not given, use square root Lasso
  #
  p <- ncol(X);
  n <- nrow(X);
  
  if  (is.null(lambda)){
    lambda <- sqrt(qnorm(1-(0.1/p))/n);
    outLas <- slim(X,y,lambda=c(lambda),method="lq",q=2,verbose=FALSE);
    # Objective : sqrt(RSS/n) +lambda *penalty
    if (intercept==TRUE) {
      return (c(as.vector(outLas$intercept),as.vector(outLas$beta)))
    }  else {
      return (as.vector(outLas$beta));
    }
  } else {
    outLas <- glmnet(X, y, family = c("gaussian"), alpha =1, intercept = intercept );
    # Objective :1/2 RSS/n +lambda *penalty
    if (intercept==TRUE){
      return (as.vector(coef(Las,s=lambda)));
    } else {
      return (as.vector(coef(Las,s=lambda))[2:(p+1)]);
    }
  }
}

SoftThreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}