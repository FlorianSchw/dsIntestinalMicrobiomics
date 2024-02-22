#'
#' @title Computes the association of microbiome data with covariates
#' @description This function calls a custom version of the native R function IFAA from the IFAA package.
#' @details internal function for microbiomeIFAAPooledDS and microbiomeMZILNPooled functions.








int.allUserFunc <- function() {
  return(
    c(
      "int.dataRecovTrans",
      "int.AIcalcu",
      "int.lm_sparse"
    )
  )
}


