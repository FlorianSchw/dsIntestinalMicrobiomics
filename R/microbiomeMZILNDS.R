#'
#' @title Computes the multivariate zero-inflated logistic normal model
#' @description This function calls the native R function from the IFAA package
#' @details The function computes a model from a SummarizedExperiment object with a given set of
#' microbiome taxa and covariates.
#' @param SumExp is a string character of the SummarizedExperiment object.
#' @param microbVar_ds This takes a single or vector of microbiome variable names (e.g., taxa, OTU and ASV names) of interest. Default is "all" meaning all microbiome variables will be analyzed. If a subset of microbiome variables is specified, the output will only contain the specified variables, and p-value adjustment for multiple testing will only be applied to the subset.
#' @param refTaxa_ds is a string character for the microbiome variable denominator (can also be a vector of microbiome variables).
#' @param allCov_ds is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param sampleIDname is a string character for the sample ID variable.
#' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for the p.adjust function in R.
#' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with 'covariates'.
#' @param paraJobs If 'sequentialRun' is FALSE, this specifies the number of parallel jobs that will be registered to run the algoithm. If specified as NULL, it will automatically detect the cores to decide the number of parallel jobs.
#' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. Default is 0 which means that taxon without any sequencing reads will be dropped from the analysis.
#' @param standardize is a logical. If TRUE, the design matrix for X will be standardized in the analyses and the results. Default is FALSE.
#' @param verbose Whether the process message is printed out to the console. Default is TRUE.
#' @return {microbiomeMZILNDS} returns the outcome of the specified multivariate zero-inflated logistic normal model
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import IFAA
#' @export
#'


microbiomeMZILNDS <- function(SumExp,
                              microbVar_ds,
                              refTaxa_ds,
                              allCov_ds,
                              sampleIDname,
                              adjust_method,
                              fdrRate,
                              paraJobs,
                              taxDropThresh,
                              standardize,
                              verbose){

  SumExp <- eval(parse(text=SumExp), envir = parent.frame())
  sequentialRun = TRUE



  if(!(is.null(microbVar_ds))){
    microbVar <- unlist(strsplit(microbVar_ds, split=","))
  } else {
    microbVar <- microbVar_ds
  }

  if(!(is.null(refTaxa_ds))){
    refTaxa <- unlist(strsplit(refTaxa_ds, split=","))
  } else {
    refTaxa <- refTaxa_ds
  }

  if(!(is.null(allCov_ds))){
    allCov <- unlist(strsplit(allCov_ds, split=","))
  } else {
    allCov <- allCov_ds
  }

  # Computes the model

  outcome <- IFAA::MZILN(experiment_dat = SumExp,
                        microbVar = microbVar,
                        refTaxa = refTaxa,
                        allCov = allCov,
                        sampleIDname = sampleIDname,
                        adjust_method = adjust_method,
                        fdrRate = fdrRate,
                        paraJobs = paraJobs,
                        taxDropThresh = taxDropThresh,
                        standardize = standardize,
                        sequentialRun = sequentialRun,
                        verbose = verbose)



  RefTaxa <- unlist(outcome[[1]][1])
  Taxa <- unlist(outcome[[1]][2])
  Covariate <- unlist(outcome[[1]][3])
  Estimate <- unlist(outcome[[1]][4])
  Std.Error <- unlist(outcome[[1]][5])
  CI_lower <- unlist(outcome[[1]][6])
  CI_upper <- unlist(outcome[[1]][7])
  p.value_adj <- unlist(outcome[[1]][8])
  p.value_unadj <- unlist(outcome[[1]][9])
  Significance <- unlist(outcome[[1]][10])

  FDR <- fdrRate
  Adjustment_Method <- adjust_method



  output <- data.frame(Taxa,
                       RefTaxa,
                       Covariate,
                       Estimate,
                       Std.Error,
                       CI_lower,
                       CI_upper,
                       p.value_unadj,
                       p.value_adj,
                       FDR,
                       Significance,
                       Adjustment_Method)

  # the resulting model will be returned to the analyst
  return(output)



}

# aggregate function
