#'
#' @title Computes the association of microbiome data with covariates
#' @description This function calls a custom version of the native R function IFAA from the IFAA package.
#' @details internal function for microbiomeIFAAPooledDS and microbiomeMZILNPooled functions.
#' @param x Needs update
#' @param y Needs update
#' @param intercept Needs update
#' @param tol Needs update
#' @importFrom MatrixExtra as.csc.matrix
#' @importFrom MatrixExtra cbind2
#' @importFrom Matrix qrR
#' @importFrom Matrix qr
#' @importFrom MatrixExtra crossprod
#'

int.lm_sparse <- function(x, y,intercept=FALSE, tol = 1e-07){

  nobs   <- length(y)
  nvar   <- ncol(x) + intercept

  if(!methods::is(x,"dgCMatrix")) x <- MatrixExtra::as.csc.matrix(x)
  if (intercept) {
    x <- MatrixExtra::cbind2(1,x)
    qrX <- base::qr(x, tol = tol)
    dfr  <- nobs - qrX$rank
    keep <- qrX$pivot[seq_len(qrX$rank)]
    newX   <- x[ , keep]
  }else{
    qrX1 <- Matrix::qr(x)
    qrX <- base::qr(as.matrix(Matrix::qrR(qrX1)), tol = tol )
    dfr  <- nobs - qrX$rank
    keep <- qrX$pivot[seq_len(qrX$rank)]
    newX   <- x[ , keep]
  }

  xnames <- colnames(x)

  XTX    <- MatrixExtra::crossprod(newX)
  Xy     <- MatrixExtra::crossprod(newX, y)
  yy     <- MatrixExtra::crossprod(y)

  outcome <- list(XTX, Xy, yy, nvar, dfr, keep, xnames, intercept)


  return(outcome)
}
