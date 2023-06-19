#'
#' @title Computes the SummarizedExperiment object
#' @description This function calls the native R function from Bioconductor
#' @details The function computes a SummarizedExperiment object from a data.frame with a given set of
#' microbiome and covariate data.
#' @param df is a string character of the data.frame
#' @param microbiomeData is a string character of microbiome variables of interest (can be a vector of names)
#' @param covariateData is a string character of covariate variables of interest (can be a vector of names)
#' @return the object specified by the \code{newobj} argument of \code{ds.summarizedExperiment} or default name \code{sumExp.newobj}
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import SummarizedExperiment
#' @export
#'


summarizedExperimentDS <- function(microbiomeData, covariateData){

  microbiomeData <- eval(parse(text=microbiomeData), envir = parent.frame())
  covariateData <- eval(parse(text=covariateData), envir = parent.frame())


  # Computes the summarizedExperiment object

  outcome <- SummarizedExperiment::SummarizedExperiment(assays = list(MicrobiomeData = t(microbiomeData)), colData =covariateData)


  # the SummarizedExperiment object is assigned to the data servers
  return(outcome)



}


