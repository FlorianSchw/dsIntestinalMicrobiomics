#'
#' @title Computes the multivariate zero-inflated logistic normal model
#' @description This function calls the native R function from the IFAA package
#' @details The function computes a model from a SummerizedExperiment object with a given set of
#' microbiome taxa and covariates.
#' @param SumExp is a string character of the data.frame
#' @param taxa is a string character for the microbiome variable denominator (can also be a vector of microbiome variables)
#' @param covariates is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param sampleIDname is a string character for the sample ID variable.
#' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for the p.adjust function in R.
#' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with 'covariates'.
#' @param paraJobs If 'sequentialRun' is FALSE, this specifies the number of parallel jobs that will be registered to run the algoithm. If specified as NULL, it will automatically detect the cores to decide the number of parallel jobs.
#' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates for the high dimensional regression.Default is 500.
#' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. Default is 0 which means that taxon without any sequencing reads will be dropped from the analysis.
#' @param standardize is a logical. If TRUE, the design matrix for X will be standardized in the analyses and the results. Default is FALSE.
#' @param sequentialRun is a logical. Default is TRUE. It can be set to be FALSE to increase speed if there are multiple taxa in the argument 'taxa'.
#' @param verbose Whether the process message is printed out to the console. Default is TRUE.
#' @param seed Random seed for reproducibility. Can be set to NULL to remove seeding.
#' @return {microbiomeMZILNDS} returns the outcome of the specified multivariate zero-inflated logistic normal model
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import IFAA
#' @export
#'


microbiomeMZILNDS <- function(SumExp, taxa, covariates, sampleIDname, adjust_method, fdrRate, paraJobs, bootB, taxDropThresh, standardize, sequentialRun, verbose, seed){

  SumExp <- eval(parse(text=SumExp), envir = parent.frame())


  # Computes the model

  outome <- IFAA::MZILN(experiment_dat = SumExp,
                        refTaxa = taxa,
                        allCov = covariates,
                        sampleIDname = sampleIDname,
                        adjust_method = adjust_method,
                        fdrRate = fdrRate,
                        paraJobs = paraJobs,
                        bootB = bootB,
                        taxDropThresh = taxDropThresh,
                        standardize = standardize,
                        sequentialRun = sequentialRun,
                        verbose = verbose,
                        seed = seed)



  Results.Ref_Tax <- unlist(results_1[[1]][1])
  Results.Taxon <- unlist(results_1[[1]][2])
  Results.Covariate <- unlist(results_1[[1]][3])
  Results.Estimate <- unlist(results_1[[1]][4])
  Results.SE.est <- unlist(results_1[[1]][5])
  Results.CI.low <- unlist(results_1[[1]][6])
  Results.CI.high <- unlist(results_1[[1]][7])
  Results.adj.p.value <- unlist(results_1[[1]][8])
  Results.unadjust.p.value <- unlist(results_1[[1]][9])
  Results.sig_ind <- unlist(results_1[[1]][10])

  Results.FDR <- fdrRate
  Results.Adjust_Method <- adjust_method
  Results.Seed <- seed
  Results.Boots <- bootB



  output <- data.frame(Results.Ref_Tax, Results.Taxon, Results.Covariate, Results.Estimate, Results.SE.est, Results.CI.low, Results.CI.high, Results.adj.p.value, Results.unadjust.p.value, Results.sig_ind, Results.FDR, Results.Adjust_Method, Results.Seed, Results.Boots)

  # the resulting model will be returned to the analyst
  return(output)



}

# aggregate function
