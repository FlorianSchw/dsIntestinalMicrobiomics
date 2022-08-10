#'
#' @title Computes the multivariate zero-inflated logistic normal model
#' @description This function calls the native R function from the IFAA package
#' @details The function computes a model from a SummerizedExperiment object with a given set of
#' microbiome taxa and covariates.
#' @param SumExp is a string character of the data.frame
#' @param taxa is a string character for the microbiome variable denominator (can also be a vector of microbiome variables)
#' @param covariates is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @return XXXXX
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import IFAA
#' @export
#'


microbiomeMZILNDS <- function(SumExp, taxa, covariates){

  SumExp <- eval(parse(text=SumExp), envir = parent.frame())


  # Computes the model

  outome <- IFAA::MZILN(experiment_dat = SumExp,
                      refTaxa = taxa,
                      allCov = covariates)


  # the resulting model will be returned to the analyst
  return(outcome)



}

# aggregate function
