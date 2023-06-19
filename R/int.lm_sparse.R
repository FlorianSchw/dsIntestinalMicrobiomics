#### Needs proper Documentation as an internal function




int.lm_sparse <- function(x, y,intercept=FALSE, tol = 1e-07){

  nobs   <- length(y)
  nvar   <- ncol(x) + intercept

  if(!methods::is(x,"dgCMatrix")) x <- MatrixExtra::as.csc.matrix(x)
  if (intercept) {
    x <- MatrixExtra::cbind2(1,x)
    qrX <- base::qr(x, tol = tol)
    dfr  <- nobs - qrX$rank
    keep <- qrX$pivot[seq_len(qrX$rank)]
    rm(qrX)
    newX   <- x[ , keep]
  }else{
    qrX1 <- Matrix::qr(x)
    qrX <- base::qr( as.matrix(Matrix::qrR(qrX1)), tol = tol )
    rm(qrX1)
    dfr  <- nobs - qrX$rank
    keep <- qrX$pivot[seq_len(qrX$rank)]
    rm(qrX)
    newX   <- x[ , keep]
  }

  xnames <- colnames(x)
  rm(x)

  XTX    <- MatrixExtra::crossprod(newX)
  Xy     <- MatrixExtra::crossprod(newX, y)
  yy     <- MatrixExtra::crossprod(y)

  #### output needs to be XTX, XY,yy,nvar, dfr


  outcome <- list(XTX, Xy, yy, nvar, dfr, keep, xnames, intercept)


  return(outcome)
}
