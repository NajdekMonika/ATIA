makeSquare <- function(A,b){
  AtA  <- t(A)%*%A
  Atb <- t(A)%*%b
  return(list(A=AtA,b=Atb))
}


makeMatrix <- function(m,n){
  return(matrix(runif(m*n), m))
}

forwardElimination <- function(A, b) {
  n <- nrow(A)
  aMod <- A
  bMod <- b
  for (k in 1:(n-1)) {
    for (i in (k+1):n) {
      if (aMod[i, k] != 0) {
        factor <- aMod[i, k] / aMod[k, k]
        aMod[i, k:n] <- aMod[i, k:n] - factor * aMod[k, k:n]
        bMod[i] <- bMod[i] - factor * bMod[k]
      }
    }
  }
  
  return(list(A = aMod, b = bMod))
}


backSubstitution <- function(A, b) {
  n <- length(b)
  x <- numeric(n)
  
  x[n] <- b[n] / A[n, n]
  
  for (i in (n - 1):1) {
    x[i] <- b[i]
    for (j in (i + 1):n) {
      x[i] <- x[i] - A[i, j] * x[j]
    }
    x[i] <- x[i] / A[i, i]
  }
  
  return(x)
}

solveGauss <- function(A, b, isSquare=TRUE) {
  if (!isSquare) {
    squared <- makeSquare(A,b) 
    A <- squared$A 
    b <- squared$b}
  elim <- forwardElimination(A, b)
  aTriangular <- elim$A
  bModified <- elim$b
  x <- backSubstitution(aTriangular, bModified)
  
  return(x)
}


solveLu <- function(A, b, isSquare=TRUE) {
  if (!isSquare) {
    squared <- makeSquare(A,b) 
    A <- squared$A 
    b <- squared$b}
  luDecomp <- expand(lu(A))
  L <- luDecomp$L
  U <- luDecomp$U
  P <- luDecomp$P
  PB <- solve(P) %*% b
  Y <- forwardsolve(L, PB, k = ncol(L), upper.tri = FALSE,
                    transpose = FALSE) # LY = Pb
  X <- backsolve(U, Y, k = ncol(U), upper.tri = TRUE,
                 transpose = FALSE) # UX = Y
  
  return(X)
}


solveQr <- function(A, b) {
  X <- qr.solve(A,b)
  
  return(X)
}


solveSvd <- function(A, b) {
  svdResult <- svd(A)
  U <- svdResult$u
  S <- svdResult$d
  V <- svdResult$v
  
  tol <- 1e-6
  r <- sum(S > tol)
  U <- U[, 1:r]
  S <- S[1:r]
  V <- V[, 1:r]
  
  
  Si <- diag(1 / S)
  x <- V %*% Si %*% t(U) %*% b
  
  return(x)
}

solveChol <- function(A,b, isSquare = TRUE){
  if (!isSquare) {
    squared <- makeSquare(A,b) 
    A <- squared$A 
    b <- squared$b}
  L <- chol(A)
  y <- forwardsolve(L, b, k = ncol(L), upper.tri = FALSE,
                    transpose = FALSE) #Ly = b
  LT <- t(L)
  x <- backsolve(LT, y, k = ncol(LT), upper.tri = TRUE,
                 transpose = FALSE) #L^T x = y
  return(x)
}


solveChol2Inv <- function(A,b, isSquare=TRUE){
  if (!isSquare) {
    squared <- makeSquare(A,b) 
    A <- squared$A 
    b <- squared$b}
  L <- chol(A)
  Ai <- chol2inv(L)
  x <- Ai %*% b
  return(x)
}


checkTime <- function(matrixList, vectorList, testedFun, isSquare=TRUE) {
  timeDf <- data.frame(N = numeric(), `t[s]` = numeric())
  for (i in seq_along(matrixList)) {
    startTime <- Sys.time()
    if(!isSquare){
      testedFun(matrixList[[i]], vectorList[[i]], isSquare)
    }else{
      testedFun(matrixList[[i]], vectorList[[i]])
    }
    endTime <- Sys.time()
    time <- as.numeric(endTime - startTime)
    timeDf <-
      rbind(timeDf, data.frame(N = dim(matrixList[[i]])[1], t = time))
  }
  timeDf
}


drawChartSquare <- function() {
  myMatrices <- list()
  symmetricMatrices <- list()
  myVectors <- list()
  for (num in seq(5, 1000, by = 100)) {
    newMatrix <- makeMatrix(num,num)
    newVector <- makeMatrix(num,1)
    myMatrices[[length(myMatrices) + 1]] <- newMatrix
    symmetricMatrices[[length(symmetricMatrices) + 1]] <- t(newMatrix)%*%newMatrix
    myVectors[[length(myVectors) + 1]] <- newVector
  }
  gauss <- checkTime(myMatrices, myVectors, solveGauss)
  lu <- checkTime(myMatrices, myVectors, solveLu)
  chol <- checkTime(symmetricMatrices, myVectors, solveChol)
  chol2inv <- checkTime(symmetricMatrices, myVectors, solveChol2Inv)
  qr <- checkTime(myMatrices, myVectors, solveQr)
  svd <- checkTime(myMatrices, myVectors, solveSvd)
  ggplot() +
    geom_line(data = gauss, aes(N, t, color = "gauss")) +
    geom_line(data = lu, aes(N, t, color = "lu")) +
    geom_line(data = chol, aes(N, t, color = "chol")) +
    geom_line(data = chol2inv, aes(N, t, color = "chol2inv")) +
    geom_line(data = qr, aes(N, t, color = "qr")) +
    geom_line(data = svd, aes(N, t, color = "svd")) +
    labs(color = "Method")
}

drawChartRectangular <- function() {
  myMatrices <- list()
  myVectors <- list()
  for (num in seq(5, 1000, by = 100)) {
    newMatrix <- makeMatrix(num,num-1)
    newVector <- makeMatrix(num,1)
    myMatrices[[length(myMatrices) + 1]] <- newMatrix
    myVectors[[length(myVectors) + 1]] <- newVector
  }
  gauss <- checkTime(myMatrices, myVectors, solveGauss, FALSE)
  lu <- checkTime(myMatrices, myVectors, solveLu, FALSE)
  chol <- checkTime(myMatrices, myVectors, solveChol, FALSE)
  chol2inv <- checkTime(myMatrices, myVectors, solveChol2Inv, FALSE)
  qr <- checkTime(myMatrices, myVectors, solveQr)
  svd <- checkTime(myMatrices, myVectors, solveSvd)
  
  ggplot() +
    geom_line(data = gauss, aes(N, t, color = "gauss")) +
    geom_line(data = lu, aes(N, t, color = "lu")) +
    geom_line(data = chol, aes(N, t, color = "chol")) +
    geom_line(data = chol2inv, aes(N, t, color = "chol2inv")) +
    geom_line(data = qr, aes(N, t, color = "qr")) +
    geom_line(data = svd, aes(N, t, color = "svd")) +
    labs(color = "Method")
}

