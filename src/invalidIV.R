library(MASS)
library(Matrix)
library(ivmodel)

ConcentrationParameterToZtoD <- function(concenParameter,sigma_firstStage, sigma_z, n) {
	return(sqrt(concenParameter * sigma_firstStage^2 / (n*sigma_z^2) ) )
}

# unionCIs: takes unions of CIs 
#   CIMatrix: n by 2 matrix (or a two-dimensional vector) where each row represents (Lower 95%, Upper 95%)
#   convexHull: if TRUE, takes the convex hull of the intervals
#           e.g. (1,3) Union (5, 7) --> (1,7)
#           if FALSE, produces a matrix of combined intervals that do not overlap
# OUTPUT: a q by 2 matrix (regardless of whether a vector was supplied) of q disjoint intervals
unionCIs <- function(CIMatrix,convexHull=TRUE)  {
  if(missing(CIMatrix)) stop("Supply CIMatrix")
  
  if(!is.matrix(CIMatrix)) CIMatrix = matrix(CIMatrix,1,length(CIMatrix))
  if(ncol(CIMatrix) != 2) stop("CIMatrix must be an n by 2 matrix!")
  
  # Remove empty intervals
  CIMatrix.emptyIndex = is.na(CIMatrix[,1])
  if(sum(CIMatrix.emptyIndex) == nrow(CIMatrix)) return(matrix(c(NA,NA),1,2))
  CIMatrix.nonempty = CIMatrix[!CIMatrix.emptyIndex,,drop=FALSE]
  
  # Check whether the first column of CIMatrix is always smaller than the second column of CIMatrix
  wrongOrderedIntervals = CIMatrix.nonempty[,1] > CIMatrix.nonempty[,2]
  if(any(wrongOrderedIntervals)) {
    print("Some of your intervals in CIMatrix have wrong ordering (i.e. first column is bigger than second column)")
	print("We'll automatically fix that for you :)")
	CIMatrix.nonempty[wrongOrderedIntervals,] = CIMatrix.nonempty[wrongOrderedIntervals,c(2,1)]
  }
  
  # If there is only one interval
  if(nrow(CIMatrix.nonempty) == 1) return(CIMatrix.nonempty)
  
  # Finally do computation
  minCIMatrix = min(CIMatrix.nonempty) 
  maxCIMatrix = max(CIMatrix.nonempty)
  if(convexHull) {
    return(matrix(c(minCIMatrix,maxCIMatrix),1,2))
  } else {
	#Sort the confidence intervals based on the lower end of the intervals
	CIMatrix.nonempty.sorted = CIMatrix.nonempty[sort(CIMatrix.nonempty[,1],index.return=TRUE)$ix,,drop=FALSE]
	  
	# declare the grand intervals for unioning
	CIMatrix.final = CIMatrix.nonempty.sorted
	prev = 1
	for(k in 2:(nrow(CIMatrix.nonempty.sorted))) {
	  # check to see whether the kth interval's left-end is less than the grand interval's right-end
      if(CIMatrix.nonempty.sorted[k,1] <= CIMatrix.final[prev,2]) {
        # If we have the case where (l1,u1) and (l2,u2) where l2 <= u1
        CIMatrix.final[prev,2] = max(CIMatrix.nonempty.sorted[k,2],CIMatrix.final[prev,2])
      } else {
		# If we have the case where(l1,u1) and (l2,u2) such that u1 < l2
		# we increment the prev, copy current kth CImatrix element, and loop again.
		prev = prev + 1
		CIMatrix.final[prev,] = CIMatrix.nonempty.sorted[k,]
	  }
    }
	return(CIMatrix.final[1:prev,,drop=FALSE])
  }
}

invalidIVCI <- function(Y,D,Z,X,U=ceiling(ncol(Z)/2),alpha=0.05,alpha.overid = 0.01, method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR"), intercept=TRUE,convexHull = FALSE) {
  if(missing(Y) || missing(D) || missing(Z)) stop("Missing outcome, Y, and exposure, D.")
  if(alpha < alpha.overid) stop("Overidentifying test cannot have alpha greater than the overall Type I error level")
  if(U >= ncol(Z) + 1 || U < 1) stop("U needs to be between 1,2,3,...,L")
 
  L = ncol(Z)
  nUnions = choose(L,L-U+1)
  
  # Preallocate output vector
  resultOut = c()
  for(j in 1:length(method)) {
  	if(method[j] == "AR") {
      resultOut = c(resultOut,list(AR = matrix(NA,nUnions*2,2)))
      curIndexAR = 1
  	}
  	if(method[j] == "CLR") {
  	  resultOut = c(resultOut,list(CLR = matrix(NA,nUnions*2,2)))
  	  curIndexCLR = 1
  	}
  	if(method[j] == "TSLS") {
  	  resultOut = c(resultOut,list(TSLS = matrix(NA,nUnions,2)))
  	  curIndexTSLS = 1
  	}
  	if(method[j] == "SarganTSLS") {
  	  resultOut = c(resultOut,list(SarganTSLS = matrix(NA,nUnions,2)))
  	  curIndexSarganTSLS = 1
  	}
  	if(method[j] == "SarganCLR") {
      resultOut = c(resultOut,list(SarganCLR = matrix(NA,nUnions*2,2)))
      curIndexSarganCLR = 1
    }
  }
  if(length(resultOut) == 0) stop("Specify a correct method(s)")
  
  # Go through each subset
  idSubsetSize = combn(1:L,L-U+1)
  
  for(i in 1:ncol(idSubsetSize)) {
    valid = idSubsetSize[,i]; invalid = setdiff(1:L,valid) # valid subset size: L - U + 1; invalid subset size: U - 1
	Zin = Z[,valid,drop=FALSE]
	if(missing(X) || length(X) == 0) {
	  Xin = Z[,invalid,drop=FALSE]
    } else {
      Xin = cbind(Z[,invalid,drop=FALSE],X)
    }
    if(ncol(Xin) == 0) {
      modelResult = ivmodel(Y,D,Z=Zin,intercept=intercept,alpha = alpha,k=1)
    } else {
      modelResult = ivmodel(Y,D,Z=Zin,X = Xin, intercept=intercept,alpha = alpha,k=1)
    }	
	if(length(intersect(c("SarganTSLS","SarganCLR"),method)) > 0) {
	  if(ncol(Xin) == 0) {
        modelResultSargan = ivmodel(Y,D,Z=Zin,intercept=intercept,alpha = alpha - alpha.overid,k=1)
      } else {
        modelResultSargan = ivmodel(Y,D,Z=Zin,X = Xin, intercept=intercept,alpha = alpha - alpha.overid,k=1)
      }	
	  # Sargan test #
	  TSLS=sum(modelResult$Dadj*qr.fitted(modelResult$ZadjQR, modelResult$Yadj))/sum(modelResult$Dadj*qr.fitted(modelResult$ZadjQR, modelResult$Dadj))
      e=modelResult$Yadj-modelResult$Dadj*TSLS
      Sargan=sum(qr.fitted(modelResult$ZadjQR, e)^2)/(sum(e^2)/length(e))
      if(modelResult$L - 1 == 0){
        pval = NA
        }else{
          pval=1-pchisq(Sargan, df=modelResult$L-1) 
        }
	}   
	for(j in 1:length(method)) {
	 if(method[j] == "AR") {
	   resultOut$AR[curIndexAR:(curIndexAR+nrow(modelResult$AR$ci) -1),] = modelResult$AR$ci
	   curIndexAR = curIndexAR + nrow(modelResult$AR$ci)
	 }
  	 if(method[j] == "CLR") {
  	   resultOut$CLR[curIndexCLR:(curIndexCLR+nrow(modelResult$CLR$ci) -1),] = modelResult$CLR$ci
  	   curIndexCLR = curIndexCLR + nrow(modelResult$CLR$ci)
  	 }
  	 if(method[j] == "TSLS") {
  	   resultOut$TSLS[curIndexTSLS:(curIndexTSLS+nrow(modelResult$kClass$ci) -1),] = modelResult$kClass$ci
  	   curIndexTSLS = curIndexTSLS + nrow(modelResult$kClass$ci)
  	 }	
  	 if(method[j] == "SarganTSLS") {
  	   if(!is.na(pval) & pval > alpha.overid)  {
  	     resultOut$SarganTSLS[curIndexSarganTSLS:(curIndexSarganTSLS+nrow(modelResultSargan$kClass$ci) - 1),] = modelResultSargan$kClass$ci
  	     curIndexSarganTSLS = curIndexSarganTSLS + nrow(modelResultSargan$kClass$ci)
  	   } else {
  	   	 resultOut$SarganTSLS[curIndexSarganTSLS,] = c(NA,NA)
  	     curIndexSarganTSLS = curIndexSarganTSLS + 1
  	   }
  	} 
  	if(method[j] == "SarganCLR") {
  	  if(!is.na(pval) & pval > alpha.overid)  {
  	     resultOut$SarganCLR[curIndexSarganCLR:(curIndexSarganCLR+nrow(modelResultSargan$CLR$ci) - 1),] = modelResultSargan$CLR$ci
  	     curIndexSarganCLR = curIndexSarganCLR + nrow(modelResultSargan$CLR$ci)
  	   } else {
  	   	 resultOut$SarganCLR[curIndexSarganCLR,] = c(NA,NA)
  	     curIndexSarganCLR = curIndexSarganCLR + 1
  	   }
  	}
  }
  }
  for(j in 1:length(method)) {
    if(method[j] == "AR") resultOut$AR = unionCIs(resultOut$AR, convexHull = convexHull)
    if(method[j] == "CLR") resultOut$CLR = unionCIs(resultOut$CLR, convexHull = convexHull)
    if(method[j] == "TSLS") resultOut$TSLS = unionCIs(resultOut$TSLS, convexHull = convexHull)
    if(method[j] == "SarganTSLS") resultOut$SarganTSLS = unionCIs(resultOut$SarganTSLS, convexHull = convexHull)
    if(method[j] == "SarganCLR") resultOut$SarganCLR = unionCIs(resultOut$SarganCLR, convexHull = convexHull)
  }
  return(resultOut) 
}

  
invalidIVCIParallel <- function(Y,D,Z,X,U=ceiling(ncol(Z)/2),alpha=0.05,alpha.overid = 0.01, method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR"), intercept=TRUE,convexHull = FALSE) {
  if(missing(Y) || missing(D) || missing(Z)) stop("Missing outcome, Y, and exposure, D.")
  if(alpha < alpha.overid) stop("Overidentifying test cannot have alpha greater than the overall Type I error level")
  if(U >= ncol(Z) + 1 || U < 1) stop("U needs to be between 1,2,3,...,L")
  if(missing(X)) X = c()
  L = ncol(Z)
  
  # Go through each subset
  idSubsetSize = combn(1:L,L-U+1)
  
  # Paralllel version of this
  resultOutParallel = foreach(i=1:ncol(idSubsetSize),.combine="rbind") %dopar% {
    resultOut = matrix(NA,2,length(method)*2)
    colnames(resultOut) = c("AR.L","AR.U","CLR.L","CLR.U","TSLS.L","TSLS.U","SarganTSLS.L","SarganTSLS.U","SarganCLR.L","SarganCLR.U")
  	valid = idSubsetSize[,i]; invalid = setdiff(1:L,valid) # valid subset size: L - U + 1; invalid subset size: U - 1
	Zin = Z[,valid,drop=FALSE]
	if(length(X) == 0) {
	  Xin = Z[,invalid,drop=FALSE]
    } else {
      Xin = cbind(Z[,invalid,drop=FALSE],X)
    }
    if(ncol(Xin) == 0) {
      modelResult = ivmodel(Y,D,Z=Zin,intercept=intercept,alpha = alpha,k=1)
    } else {
      modelResult = ivmodel(Y,D,Z=Zin,X = Xin, intercept=intercept,alpha = alpha,k=1)
    }	
	if(length(intersect(c("SarganTSLS","SarganCLR"),method)) > 0) {
	  if(ncol(Xin) == 0) {
        modelResultSargan = ivmodel(Y,D,Z=Zin,intercept=intercept,alpha = alpha - alpha.overid,k=1)
      } else {
        modelResultSargan = ivmodel(Y,D,Z=Zin,X = Xin, intercept=intercept,alpha = alpha - alpha.overid,k=1)
      }	
	  # Sargan test #
	  TSLS=sum(modelResult$Dadj*qr.fitted(modelResult$ZadjQR, modelResult$Yadj))/sum(modelResult$Dadj*qr.fitted(modelResult$ZadjQR, modelResult$Dadj))
      e=modelResult$Yadj-modelResult$Dadj*TSLS
      Sargan=sum(qr.fitted(modelResult$ZadjQR, e)^2)/(sum(e^2)/length(e))
      pval=1-pchisq(Sargan, df=modelResult$L-1)
	}   
	for(j in 1:length(method)) {
	 if(method[j] == "AR") resultOut[1:nrow(modelResult$AR$ci),c("AR.L","AR.U")] = modelResult$AR$ci
  	 if(method[j] == "CLR") resultOut[1:nrow(modelResult$CLR$ci),c("CLR.L","CLR.U")] = modelResult$CLR$ci
  	 if(method[j] == "TSLS") resultOut[1,c("TSLS.L","TSLS.U")] =  modelResult$kClass$ci
  	 if(method[j] == "SarganTSLS") {
  	   if(pval > alpha.overid)  {
  	     resultOut[1,c("SarganTSLS.L","SarganTSLS.U")] = modelResultSargan$kClass$ci
  	   }
  	} 
  	if(method[j] == "SarganCLR") {
  	  if(pval > alpha.overid)  {
  	     resultOut[1:nrow(modelResultSargan$CLR$ci),c("SarganCLR.L","SarganCLR.U")] = modelResultSargan$CLR$ci
  	   } 
  	}
    }
    resultOut
  }
  result = list()
  for(j in 1:length(method)) {
    if(method[j] == "AR") result$AR = unionCIs(resultOutParallel[,c("AR.L","AR.U")], convexHull = convexHull)
    if(method[j] == "CLR") result$CLR = unionCIs(resultOutParallel[,c("CLR.L","CLR.U")], convexHull = convexHull)
    if(method[j] == "TSLS") result$TSLS = unionCIs(resultOutParallel[,c("TSLS.L","TSLS.U")], convexHull = convexHull)
    if(method[j] == "SarganTSLS") result$SarganTSLS = unionCIs(resultOutParallel[,c("SarganTSLS.L","SarganTSLS.U")], convexHull = convexHull)
    if(method[j] == "SarganCLR") result$SarganCLR = unionCIs(resultOutParallel[,c("SarganCLR.L","SarganCLR.U")], convexHull = convexHull)
  }
  return(result)
}

statsForInvalidIVCI <- function(resultOut,trueBetaGrid = seq(-10,10,0.1),method=c("AR","CLR","TSLS","SarganTSLS","SarganCLR"),trueBeta) {
  rejectNull = matrix(0,length(trueBetaGrid),length(method))
  medianLength = matrix(0,1,length(method))
  coverageProportion = matrix(0,1,length(method))
  colnames(rejectNull) = colnames(medianLength) = colnames(coverageProportion) = method
  for(i in 1:length(method)) {
  	if(method[i] == "AR") {
  		if(any(is.na(resultOut$AR[,1]))) {
  			rejectNull[,"AR"] = 1
  			medianLength[,"AR"] = 0
  			coverageProportion[,"AR"] = 0
  		} else {
  			rejectNull[,"AR"] = (trueBetaGrid <= min(resultOut$AR[,1]) | trueBetaGrid >= max(resultOut$AR[,2]))
  			medianLength[,"AR"] = sum( abs(resultOut$AR[,2] - resultOut$AR[,1]))
  			coverageProportion[,"AR"] = any((resultOut$AR[,1] <= trueBeta) & (trueBeta <= resultOut$AR[,2]))
  		}
  	}
  	if(method[i] == "CLR") {
  		rejectNull[,"CLR"] = (trueBetaGrid <= min(resultOut$CLR[,1]) | trueBetaGrid >= max(resultOut$CLR[,2]))
  		medianLength[,"CLR"] = sum( abs(resultOut$CLR[,2] - resultOut$CLR[,1]))
  		coverageProportion[,"CLR"] = any((resultOut$CLR[,1] <= trueBeta) & (trueBeta <= resultOut$CLR[,2]))
  	}
  	if(method[i] == "TSLS") {
  		rejectNull[,"TSLS"] = (trueBetaGrid <= min(resultOut$TSLS[,1]) | trueBetaGrid >= max(resultOut$TSLS[,2]))
  		medianLength[,"TSLS"] = sum( abs(resultOut$TSLS[,2] - resultOut$TSLS[,1]))
  		coverageProportion[,"TSLS"] = any((resultOut$TSLS[,1] <= trueBeta) & (trueBeta <= resultOut$TSLS[,2]))
  	}
  	if(method[i] == "SarganTSLS") {
  		if(any(is.na(resultOut$SarganTSLS[,1]))) {
  			rejectNull[,"SarganTSLS"] = 1
  			medianLength[,"SarganTSLS"] = 0
  			coverageProportion[,"SarganTSLS"] = 0
  		} else {
  		  rejectNull[,"SarganTSLS"] = (trueBetaGrid <= min(resultOut$SarganTSLS[,1]) | trueBetaGrid >= max(resultOut$SarganTSLS[,2]))
  		  medianLength[,"SarganTSLS"] = sum( abs(resultOut$SarganTSLS[,2] - resultOut$SarganTSLS[,1]))
  		  coverageProportion[,"SarganTSLS"] = any((resultOut$SarganTSLS[,1] <= trueBeta) & (trueBeta <= resultOut$SarganTSLS[,2]))
  		}
  	}
  	if(method[i] == "SarganCLR") {
  		if(any(is.na(resultOut$SarganCLR[,1]))) {
  			rejectNull[,"SarganCLR"] = 1
  			medianLength[,"SarganCLR"] = 0
  			coverageProportion[,"SarganCLR"] = 0
  		} else {
  		  rejectNull[,"SarganCLR"] = (trueBetaGrid <= min(resultOut$SarganCLR[,1]) | trueBetaGrid >= max(resultOut$SarganCLR[,2]))
  		  medianLength[,"SarganCLR"] = sum( abs(resultOut$SarganCLR[,2] - resultOut$SarganCLR[,1]))
  		  coverageProportion[,"SarganCLR"] = any((resultOut$SarganCLR[,1] <= trueBeta) & (trueBeta <= resultOut$SarganCLR[,2]))
  		}
  	}
  }
  return(list(rejectNull = rejectNull,medianLength = medianLength,coverageProportion = coverageProportion))	
}