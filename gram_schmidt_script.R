tolerance <- 1e-4


checkOrthogonality <- function(myMatrix) {
  crossprod <- crossprod(myMatrix)
  all(crossprod[!diag(nrow(crossprod))] < tolerance)
}


checkOrthonormality <- function(myMatrix) {
  checkOrthogonality(myMatrix) &
    all(abs(Matrix::colSums(myMatrix ** 2) - 1.0) <= tolerance)
}


createMatrix <- function(num) {
  return(matrix(sample(0:100, num * num, replace = TRUE), nrow = num))
}


checkIndependence <- function(vMatrix) {
  dim <- dim(vMatrix)
  if (dim[1] == dim[2]) {
    return(det(vMatrix) != 0)
  } else if (dim[2] > dim[1]) {
    return(FALSE)
  } else{
    return(qr(vMatrix)$rank == dim[2])
  }
}


gramSchmidt <- function(vMatrix) {
  if (is.matrix(vMatrix) &
      checkIndependence(vMatrix) &
      !checkOrthogonality(vMatrix) & !checkOrthonormality(vMatrix)) {
    uMatrix <- cbind(vMatrix[, 1] * (1 / sqrt(sum(vMatrix[, 1] ^ 2))))
    if (ncol(vMatrix) > 1) {
      for (vNum in 2:ncol(vMatrix)) {
        u <- vMatrix[, vNum]
        for (uNum in 1:ncol(uMatrix)) {
          scalarQuot <-
            (vMatrix[, vNum] %*% uMatrix[, uNum]) / crossprod(uMatrix[, uNum])
          u <- u - uMatrix[, uNum] %*% scalarQuot
        }
        e <- u * (1 / sqrt(sum(u ^ 2)))
        uMatrix <- cbind(uMatrix, e)
      }
    }
    uMatrix
  }
}



gramSchmidtNoLoop <- function(vMatrix) {
  if (is.matrix(vMatrix) &
      checkIndependence(vMatrix) &
      !checkOrthogonality(vMatrix) & !checkOrthonormality(vMatrix)) {
    QR <- qr(vMatrix)
    qr.Q(QR)
  }
}



checkTime <- function(matrixList, testedFun) {
  timeDf <- data.frame(N = numeric(), `t[s]` = numeric())
  for (matrix in matrixList) {
    startTime <- Sys.time()
    testedFun(matrix)
    endTime <- Sys.time()
    time <- as.numeric(endTime - startTime)
    timeDf <-
      rbind(timeDf, data.frame(N = dim(matrix)[1], t = time))
  }
  timeDf
}


drawChart <- function() {
  myMatrices <- list()
  nList <- list()
  for (num in seq(1, 1000, by = 100)) {
    newMatrix = createMatrix(num)
    myMatrices[[length(myMatrices) + 1]] <- newMatrix
    nList <- c(nList, num)
  }
  loopDf <- checkTime(myMatrices, gramSchmidt)
  noLoopDf <- checkTime(myMatrices, gramSchmidtNoLoop)
  
  ggplot() +
    geom_line(data = loopDf, aes(N, t, color = "Loop")) +
    geom_line(data = noLoopDf, aes(N, t, color = "No Loop")) +
    labs(color = "Method")
}
