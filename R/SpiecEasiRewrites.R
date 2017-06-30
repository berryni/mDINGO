## SPiecEasi rewritten or copied functions because they werent working for me, or I needed a workaround.

igraph2prec = function (Graph, posThetaLims = c(2, 3), negThetaLims = -posThetaLims, 
                        targetCondition = 100, epsBin = 0.01, numBinSearch = 100) 
{
  if (class(Graph) != "igraph") 
    stop("input is not an igraph")
  Graph = as.matrix(Graph[])
  n <- ncol(Graph)
  posThetaLims <- sort(posThetaLims)
  negThetaLims <- sort(negThetaLims)
  Theta <- Graph + diag(n)
  degVec <- colSums(Graph)
  utri <- SpiecEasi::triu(Theta)
  nzind <- which(utri != 0)
  if (length(posThetaLims) > 2 || length(negThetaLims) > 2) 
    stop("theta_max and theta_min should be a numeric vector of length 1 or 2")
  rands <- runif(length(nzind))
  mapToRange <- function(x, lim) {
    span <- diff(sort(lim))
    min(lim) + (x * span)
  }
  boolind <- sample(c(TRUE, FALSE), length(nzind), replace = TRUE)
  rands[boolind] <- mapToRange(rands[boolind], posThetaLims)
  rands[!boolind] <- mapToRange(rands[!boolind], negThetaLims)
  utri[nzind] <- rands
  #####Theta <- triu2diag(utri, 1) <- this wasn't working. The next 2 lines do what the function did I think.
  Theta[lower.tri(Theta)] = t(Theta)[lower.tri(Theta)]
  diag(Theta) = 1
  
  eigVals <- eigen(Theta)$values
  minEig <- min(eigVals)
  maxEig <- max(eigVals)
  if (minEig < 0.01) 
    Theta <- Theta + abs(minEig) * diag(n)
  diagConst <- .binSearchCond(Theta, targetCondition, numBinSearch, 
                              epsBin)
  Theta <- Theta + diagConst * diag(n)
  return(Theta)
}

.binSearchCond <- function(Theta, condTheta, numBinSearch, epsBin) {
  # Internal function that determines the constant in the diagonal to satisfy the
  # condition constraint on the Precision/Covariance matrix
  
  n <- nrow(Theta)
  currCondTheta <- kappa(Theta)
  if (currCondTheta < condTheta) {
    # Max entry in the diagonal (lower bound)
    currLB   <- -max(diag(Theta))
    stepSize <- currLB+.Machine$double.eps
    
    while (currCondTheta < condTheta) {
      currCondTheta <- kappa(Theta+stepSize*diag(n))
      stepSize      <- stepSize/2
    }
    currUB <- stepSize
  } else {
    currLB <- 0
    stepSize = 0.1
    
    while (currCondTheta > condTheta) {
      currCondTheta <- kappa(Theta + stepSize*diag(n))
      stepSize      <- 2*stepSize
    }
    currUB <- stepSize
  }
  
  for (i in 1:numBinSearch) {
    diagConst <- (currUB+currLB)/2
    currCondTheta <- kappa(Theta+diagConst*diag(n))
    
    if (currCondTheta < condTheta) currUB <- diagConst
    else currLB <- diagConst
    
    if (abs(currCondTheta-condTheta)<epsBin) break
  }
  diagConst
}


