#' Inference for linear combination of the regression vector in high dimensional generalized linear regression
#'
#' @description Computes the bias-corrected estimator of thelinearcombination of the regression vector for the high dimensional generalized linear regression and the corresponding standard error.
#'
#' @param X Design matrix, of dimension \eqn{n} x \eqn{p}
#' @param y Outcome vector, of length \eqn{n}
#' @param model The high dimensional regression model, either \code{linear} or \code{logistic} or \code{logistic alternative} or \code{probit} or \code{inverse t1}
#' @param loading Loading, of length \eqn{p}
#' @param intercept.loading Should intercept be included for the \code{loading} (default = \code{TRUE})
#' @param intercept Should intercept be fitted for the initial estimator (default = \code{TRUE})
#' @param init.coef Initial estimator of the regression vector (default = \code{NULL})
#' @param lambda The tuning parameter in the construction of \code{init.coef} (default = \code{NULL})
#' @param mu The dual tuning parameter used in the construction of the projection direction (default = \code{NULL})
#' @param step The step size used to compute \code{mu}; if set to \code{NULL} it is
#' computed to be the number of steps (< \code{maxiter}) to obtain the smallest \code{mu}
#' such that the dual optimization problem for constructing the projection direction converges (default = \code{NULL})
#' @param resol The factor by which \code{mu} is increased/decreased to obtain the smallest \code{mu}
#' such that the dual optimization problem for constructing the projection direction converges (default = 1.5)
#' @param maxiter Maximum number of steps along which \code{mu} is increased/decreased to obtain the smallest \code{mu}
#' such that the dual optimization problem for constructing the projection direction converges (default = 6)
#' @param b0 The null value to be tested against                                                              #removed
#' @param alpha Level of significance to test the null hypothesis which claims that the linear combination of the regression coefficients
#' is less than or equal to \code{b0} (default = 0.05)                                                        #removed
#' @param verbose Should inetrmediate message(s) be printed (default = \code{TRUE})
#'
#' @return
#' \item{prop.est}{The bias-corrected estimator for the linear combination of regression coefficients}
#' \item{se}{The standard error of the bias-corrected estimator}
#' \item{CI}{The confidence interval for the linear combination}                                                     #removed
#' \item{decision}{\code{decision}\eqn{=1} implies the linear combination is above \code{b0} \eqn{\newline}          #removed
#' \code{decision}\eqn{=0} implies the linear combination is not above \code{b0}}
#' \item{proj}{The projection direction, of length \eqn{p}}
#' \item{plug.in}{The plug-in estimator for the linear combination}
#' @export
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats coef na.omit dnorm pnorm dt pt
#' @import CVXR Matrix glmnet
#'
#' @examples
#'
#' \donttest{
#' sp = 20
#' n = 400
#' p = 800
#' f = function(x){
#' pnorm(x)
#' }
#' sig1 = toeplitz(seq(0.6, 0,length.out = p/10))
#' Sig = Matrix::bdiag(rep(list(sig1),10))+diag(rep(0.4,p))
#' set.seed(1203)
#' X = MASS::mvrnorm(n, mu=rep(0,p), Sigma=Sig)
#' b = rep(0,p)
#' b[1:sp] = rep(c(0.4,-0.4), sp/2)
#' prob = f(X %*% b)
#' y = array(dim = 1)
#' for(i in 1:n){
#' set.seed(i)
#' y[i] = rbinom(1,1,prob[i])
#' }
#' Est = SIHR::GLM_LF(X = X, y = y, loading = c(1,rep(0,799)), model = "probit", intercept = FALSE)
#' }
GLM_LF<-function(X, y, model = "linear", loading, intercept.loading = TRUE, intercept = TRUE, 
                 init.coef = NULL, lambda = NULL, mu = NULL, step = NULL, resol = 1.5, maxiter = 6, b0 = 0, alpha = 0.05, verbose = TRUE){
  xnew <- loading
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  n_y <- length(y)
  
  if(n_y!=n)
  {
    stop("Error: Check dimensions of X and y")
  } else {
    data <- na.omit(data.frame(y,X))
    X <- as.matrix(data[,-1])
    y <- as.vector(data[,1])
    p <- ncol(X)
    n <- nrow(X)
    mean = colMeans(X)
    M = matrix(rep(mean,nrow(X)),byrow = T, nrow = nrow(X), ncol = ncol(X))
    X = X - M
    
    if(intercept.loading == TRUE && intercept == TRUE){
      xnew = xnew - mean
    }
    if(model == "logistic"){
      f = function(x){
        exp(x)/(1+exp(x))
      }
      fp = function(x){
        exp(x)/(1+exp(x))^2
      }
      g = function(x){
        fp(x)/f(x)/(1-f(x))
      }
    }else if(model =="probit"){
      f = function(x){
        pnorm(x)
      }
      fp = function(x){
        dnorm(x)
      }
      g = function(x){
        out=c()
        out[abs(x)<5]=(fp(x)/f(x)/(1-f(x)))[abs(x)<5]
        out[abs(x)>5]=((abs(x)+sqrt(x^2+8/pi))/2)[abs(x)>5]
        return(out)
      }
    } else if(model == "inverse t1"){
      f = function(x){
        pt(x,1)
      }
      fp = function(x){
        dt(x,1)
      }
      g = function(x){
        fp(x)/f(x)/(1-f(x))
      }
    }
    
    ############ TALK ############
    
    if(is.null(init.coef)){
      if(model == "probit"){
        init.coef <-  Initialization.step(X*1.75, y, model = "glm", lambda = 0.06*sqrt(log(p)/n), intercept)
        htheta <- init.coef$lasso.est
      } else if(model == "linear"){
        init.coef <-  Initialization.step(X, y, model = "linear", lambda, intercept)
        htheta <- init.coef$lasso.est
      } else {
        init.coef <-  Initialization.step(X, y, model = "glm", lambda = lambda, intercept)
        htheta <- init.coef$lasso.est
      }
    } else {
      htheta = init.coef
    }
    col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
    Xnor <- X %*% diag(col.norm)
    if (intercept == TRUE){
      support<-(abs(htheta)>0.001)
      Xc <- cbind(rep(1,n),X);
      Xb <- cbind(rep(1,n),Xnor)
      pp <- (p+1);
    } else {
      support<-(abs(htheta)>0.001)
      Xc <- X
      Xb <- Xnor
      pp <- p
    }
    htheta <- as.vector(htheta)
    if(model == "linear"){
      sparsity <- sum(abs(htheta) > 0.001)
      sd.est <- sqrt(sum((y - Xc %*% htheta)^2) / max(0.9*n, n - sparsity))
    }
    
    if(intercept==TRUE){
      loading <- rep(0,pp)
      if(intercept.loading == TRUE){
        loading[1] <- 1
      }
      if(intercept.loading == FALSE){
        loading[1] <- 0
      }
      loading[-1] <- xnew
    } else {
      if(intercept.loading == TRUE){
        cat(paste("Setting intercept = FALSE and intercept.loading = FALSE"))
      }
      loading <- xnew
    }
    loading.norm = sqrt(sum(loading^2))
    lasso.plugin <- sum(loading*htheta)
    exp_pred = Xc%*%(htheta)
    if(model == "linear"){
      X.weight = Xc
    }else if(model == "logistic alternative"){
      deriv.vec <- exp(exp_pred)/(1+exp(exp_pred))^2
      weight <- 1/deriv.vec
      X.weight = diag(c(sqrt(deriv.vec*weight))) %*% Xc
    } else {
      X.weight = diag(c(sqrt(fp(exp_pred)*g(exp_pred)))) %*% Xc
    }
    
    count=0
    for(i in 1:ncol(X)){
      if(length(unique(X[,i])) == 1){
        count = count+1
      }
    }
    if(count!=0 && intercept==TRUE)
    {
      stop("Data is singular")
    } else {
      if ((n>=6*p)){
        sigma.hat <- (1/n)*(t(X.weight)%*%X.weight);
        tmp <- eigen(sigma.hat)
        tmp <- min(tmp$values)/max(tmp$values)
      } else {
        tmp <- 0
      }
      
      if ((n >= 6*p) && (tmp >= 1e-4)){
        direction <- solve(sigma.hat)%*%loading/loading.norm
      } else {
        if(is.null(step)){
          step.vec <- rep(NA,3)
          for(t in 1:3){
            index.sel <- sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
            if(model == "linear"){
              Direction.Est.temp <-  Direction_searchtuning(Xc[index.sel,], loading, model = "linear", resol, maxiter)
            } else if(model == "logistic alternative"){
              Direction.Est.temp <-  Direction_searchtuning(Xc[index.sel,], loading, model = "glm", weight = weight[index.sel], deriv.vec = deriv.vec[index.sel], resol, maxiter)
            } else{
              Direction.Est.temp <-  Direction_searchtuning(Xc[index.sel,], loading, model = "glm", weight = g(exp_pred)[index.sel], deriv.vec = fp(exp_pred)[index.sel], resol, maxiter)
            }
            step.vec[t] <- Direction.Est.temp$step
          }
          step<- getmode(step.vec)
        }
        if(model == "linear"){
          Direction.Est <-  Direction_fixedtuning(Xc, loading, model = "linear", mu = sqrt(2.01*log(pp)/n)*resol^{-(step-1)})
        } else if(model == "logistic alternative"){
          Direction.Est <-  Direction_fixedtuning(Xc, loading, model = "glm", mu = sqrt(2.01*log(pp)/n)*resol^{-(step-1)}, weight = weight, deriv.vec = deriv.vec)
        } else{
          Direction.Est <-  Direction_fixedtuning(Xc, loading, model = "glm", mu = sqrt(2.01*log(pp)/n)*resol^{-(step-1)}, weight = g(exp_pred), deriv.vec = fp(exp_pred))
        }
        
        while(is.na(Direction.Est) || length(Direction.Est$proj)==0){
          step <- step-1
          if(model == "linear"){
            Direction.Est <-  Direction_fixedtuning(Xc, loading, model = "linear", mu = sqrt(2.01*log(pp)/n)*resol^{-(step-1)})
          } else if(model == "logistic alternative"){
            Direction.Est <-  Direction_fixedtuning(Xc, loading, model = "glm", mu = sqrt(2.01*log(pp)/n)*resol^{-(step-1)}, weight = weight, deriv.vec = deriv.vec)
          } else{
            Direction.Est <-  Direction_fixedtuning(Xc, loading, model = "glm", mu = sqrt(2.01*log(pp)/n)*resol^{-(step-1)}, weight = g(exp_pred), deriv.vec = fp(exp_pred))
          }
        }
        if(verbose == TRUE){
          cat(paste("step is", step))
        }
        direction <- Direction.Est$proj
      }
      if(model == "linear"){
        weighed.residual = (y - exp_pred)
      } else if(model == "logistic alternative"){
        weighed.residual <- (y - exp(exp_pred)/(1+ exp(exp_pred)))*weight
      } else {
        weighed.residual=(y - f(exp_pred))*g(exp_pred)
      }
      
      correction = sum((Xc%*%direction)*weighed.residual)/n;
      debias.est=lasso.plugin+correction*loading.norm
      
      if(model == "linear"){
        se <- sd.est*sqrt(sum((Xc%*%direction)^2)/(n)^2)*loading.norm
      } else if(model == "logistic alternative"){
        se <- sqrt(mean((Xc%*%direction)^2*weight^2*deriv.vec))*loading.norm/sqrt(n)
      } else {
        X.weight2 = diag(c(g(exp_pred)*sqrt(f(exp_pred)*(1-f(exp_pred))))) %*% Xc
        se<-sqrt(mean((X.weight2%*%direction)^2))/sqrt(n)
      }
      CI <- c(debias.est - qnorm(1-alpha/2)*se, debias.est + qnorm(1-alpha/2)*se)
      if(model == "linear" || model == "logistic alternative"){
        b0 = 0
      }
      if(debias.est - qnorm(1-alpha)*se > b0){
        dec <- 1
      } else {
        dec <- 0
      }
      #returnList <- list("prop.est" = debias.est,
      #                   "se" = se,
      #                   "CI" = CI,
      #                   "decision" = dec,
      #                   "proj" = direction,
      #                   "step" = step,
      #                   "plug.in" = lasso.plugin
      #)
      #return(returnList)
      out <- list(prop.est = debias.est,
                  se = se,
                  proj = direction,
                  step = step,
                  plug.in = lasso.plugin)
      structure(out, class = "LF")
    }
  }
}


getmode <- function(v) {
  tbl <- table(v)
  if (all(tbl == 1)) {
    median(v)
  } else {
    as.numeric(names(which.max(tbl)))
  }
}
#getmode_log <- function(v) {
#  uniqv <- unique(v)
#  uniqv[which.max(tabulate(match(v, uniqv)))]
#}
Initialization.step <- function(X, y, model = "linear", lambda = NULL, intercept = FALSE) {
  n <- nrow(X)
  col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
  Xnor <- X %*% diag(col.norm)
  
  htheta <- Lasso.highd(Xnor, y, model = model, lambda = lambda, intercept = intercept)
  
  if (intercept == TRUE) {
    Xb <- cbind(rep(1, n), Xnor)
    col.norm <- c(1, col.norm)
  } else {
    Xb <- Xnor
  }
  htheta <- htheta * col.norm
  returnList <- list("lasso.est" = htheta)
  return(returnList)
}

Lasso.highd<- function(X, y, model = "linear", lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)
  if(model == "linear"){
    family = "gaussian"
    standardize = TRUE
  } else if(model == "glm"){
    family = "binomial"
    standardize = FALSE
  }
  
  htheta <- if (is.null(lambda)) {
    outLas <- cv.glmnet(X, y, family = family, alpha = 1,
                        intercept = intercept, standardize = standardize)
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = family, alpha = 1,
                        intercept = intercept, standardize = standardize)
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "scalreg" && model == "linear") {
    Xc <- if (intercept) {
      cbind(rep(1, n), X)
    } else {
      X
    }
    outLas <- scalreg(Xc, y)
    if (intercept) {
      outLas$coefficients
    } else {
      c(0, outLas$coefficients)
    }
  }else {
    outLas <- glmnet(X, y, family = family, alpha = 1,
                     intercept = intercept, standardize = standardize)
    as.vector(coef(outLas, s = lambda))
  }
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}

#' Constructs the projection direction with a fixed tuning parameter
#'
#' @description Constructs the projection direction, used for bias-correction, with a fixed tuning parameter
#'
#' @param X Design matrix, of dimension \eqn{n} x \eqn{p}
#' @param loading Loading, of length \eqn{p}
#' @param mu The dual tuning parameter used in the construction of the projection direction
#' @param model The high dimensional regression model, either \code{linear} or \code{glm} (default = \code{linear})
#' @param weight The weight vector of length \eqn{n}; to be supplied if \code{model="glm"} (default=\code{NULL} when \code{model=linear})
#' @param deriv.vec The first derivative vector of the logit function at \eqn{X\%*\%}(\code{init.coef}), of length \eqn{n} ; to be supplied if \code{model="glm"}. Here \code{init.coef} is the initial estimate of the regression vector. (default = \code{NULL} when \code{model=linear})
#' @return
#' \item{proj}{The projection direction, of length \eqn{p}}
#' @export
#'
#' @examples
#' n <- 100
#' p <- 400
#' set.seed(1203)
#' X <- matrix(sample(-2:2,n*p,replace = TRUE),nrow = n,ncol = p)
#' resol <- 1.5
#' step <- 3
#' Est <- Direction_fixedtuning(X,loading=c(1,rep(0,(p-1))),mu=sqrt(2.01*log(p)/n)*resol^{-(step-1)})

Direction_fixedtuning <- function(X, loading, mu = NULL, model = "linear", weight = NULL, deriv.vec = NULL){
  pp <- ncol(X)
  n <- nrow(X)
  if(is.null(mu)){
    mu <- sqrt(2.01*log(pp)/n)
  }
  loading.norm <- sqrt(sum(loading^2))
  
  if (loading.norm == 0){
    H <- cbind(loading, diag(1, pp))
  }else{
    H <- cbind(loading / loading.norm, diag(1, pp))
  }
  
  v <- Variable(pp+1)
  if(model == "linear")
  {
    obj <- 1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  } else if(model == "glm") {
    obj <- 1/4*sum(((X%*%H%*%v)^2)*weight*deriv.vec)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  } else {
    stop("Method not yet developed")
  }
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  if(result$status=="optimal" || result$status == "unbounded"){
    opt.sol<-result$getValue(v)
    cvxr_status<-result$status
    direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  }else{
    direction <- numeric(0)
  }
  returnList <- list("proj"=direction)
  return(returnList)
}

#' Searches for the best step size and computes the projection direction with the searched best step size
#'
#' @description Searches for the best step size and computes the projection direction with the searched best step size
#'
#' @param X Design matrix, of dimension \eqn{n} x \eqn{p}
#' @param loading Loading, of length \eqn{p}
#' @param model The high dimensional regression model, either \code{linear} or \code{glm} (default = \code{linear})
#' @param weight The weight vector of length \eqn{n}; to be supplied if \code{model="glm"} (default=\code{NULL} when \code{model=linear})
#' @param deriv.vec The first derivative vector of the logit function at \eqn{X\%*\%}(\code{init.coef}), of length \eqn{n} ; to be supplied if \code{model="glm"}. Here \code{init.coef} is the initial estimate of the regression vector. (default = \code{NULL} when \code{model=linear})
#' @param resol The factor by which \code{mu} is increased/decreased to obtain the smallest \code{mu}
#' such that the dual optimization problem for constructing the projection direction converges (default = 1.5)
#' @param maxiter Maximum number of steps along which \code{mu} is increased/decreased to obtain the smallest \code{mu}
#' such that the dual optimization problem for constructing the projection direction converges (default = 6)
#'
#' @return
#' \item{proj}{The projection direction, of length \eqn{p}}
#' \item{step}{The best step size}
#' @export
#'
#' @examples
#' n <- 100
#' p <- 400
#' set.seed(1203)
#' X <- matrix(sample(-2:2,n*p,replace = TRUE),nrow = n,ncol = p)
#' Est <- Direction_searchtuning(X,loading=c(1,rep(0,(p-1))))
Direction_searchtuning <- function(X, loading, model = "linear", weight = NULL, deriv.vec = NULL, resol = 1.5, maxiter = 6){
  pp <- ncol(X)
  n <- nrow(X)
  tryno <- 1
  opt.sol <- rep(0,pp+1)
  lamstop <- 0
  cvxr_status <- "optimal"
  mu <- sqrt(2.01*log(pp)/n)
  while (lamstop == 0 && tryno < maxiter){
    lastv <- opt.sol;
    lastresp <- cvxr_status;
    loading.norm <- sqrt(sum(loading^2))
    if (loading.norm == 0){
      H <- cbind(loading, diag(1, pp))
    }else{
      H <- cbind(loading / loading.norm, diag(1, pp))
    }
    v <- Variable(pp+1)
    if(model == "linear")
    {
      obj <- 1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
    } else if(model=="glm") {
      obj <- 1/4*sum(((X%*%H%*%v)^2)*weight*deriv.vec)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
    } else {
      stop("Method not yet developed")
    }
    prob <- Problem(Minimize(obj))
    result <- solve(prob)
    cvxr_status <- result$status
    if(tryno == 1){
      if(cvxr_status == "optimal"){
        incr = 0
        mu=mu/resol
        opt.sol <- result$getValue(v)
        temp.vec <- (-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
        if(model == "linear")
        {
          initial.sd <- sqrt(sum(((X%*% temp.vec)^2))/(n)^2)*loading.norm
        } else {
          initial.sd <- sqrt(sum(((X%*% temp.vec)^2)*weight*deriv.vec)/(n)^2)*loading.norm
        }
        temp.sd <- initial.sd
      } else {
        incr <- 1
        mu <- mu*resol
      }
    } else {
      if(incr == 1){
        if(cvxr_status == "optimal"){
          opt.sol <- result$getValue(v)
          lamstop <- 1
        } else {
          mu <- mu*resol
        }
      } else {
        if(cvxr_status == "optimal" && temp.sd < 3*initial.sd){
          mu <- mu/resol
          opt.sol <- result$getValue(v)
          temp.vec <- (-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
          if(model == "linear")
          {
            temp.sd <- sqrt(sum(((X%*% temp.vec)^2))/(n)^2)*loading.norm
          } else {
            temp.sd <- sqrt(sum(((X%*% temp.vec)^2)*weight*deriv.vec)/(n)^2)*loading.norm
          }
        } else {
          mu <- mu*resol
          opt.sol <- lastv
          lamstop <- 1
          tryno <- tryno-1
        }
      }
    }
    tryno = tryno + 1
  }
  direction <- (-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  step <- tryno-1
  returnList <- list("proj"=direction,
                     "step"=step)
  return(returnList)
}

