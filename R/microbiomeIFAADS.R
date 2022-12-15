#'
#' @title Computes the association of microbiome data with covariates
#' @description This function calls the native R function from the IFAA package
#' @details The function computes an association from a SummerizedExperiment object with a given set of
#' confounders and covariates.
#' @param SumExp is a string character of the data.frame
#' @param microbVar This takes a single or vector of microbiome variable names (e.g., taxa, OTU and ASV names) of interest. Default is "all" meaning all microbiome variables will be analyzed. If a subset of microbiome variables is specified, the output will only contain the specified variables, and p-value adjustment for multiple testing will only be applied to the subset.
#' @param covariates is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param confounders is a string character for the covariates that will be adjusted in the model (can also be a vector of confounders)
#' @param sampleIDname is a string character for the sample ID variable.
#' @param covariatesMany is a logical. If 'TRUE' and 'covariates' are set to NULL, then all variables in the 'covariates' will be used.
#' @param confoundersMany is a logical. If 'TRUE' and 'confounders' are set to NULL, then all variables except the 'coviariates' will be used as confounders.
#' @param nRef number of randomly picked reference taxa used in phase 1.
#' @param nRefMaxForEsti maximum number of final reference taxa used in phase 2.
#' @param taxa vector of taxa or OTU or ASV names. Theses are reference taxa specified by the user to be used in phase 1.
#' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for the p.adjust function in R.
#' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with 'covariates'.
#' @param paraJobs If 'sequentialRun' is FALSE, this specifies the number of parallel jobs that will be registered to run the algoithm. If specified as NULL, it will automatically detect the cores to decide the number of parallel jobs.
#' @param bootB number of bootstrap samples for obtaining confidence interval of estimates in phase 2 for the high dimensional regression. Default is 500.
#' @param standardize is a logical. If 'TRUE', the design matrix for X will be standardized in the analyses and the results. Default is FALSE.
#' @param sequentialRun is a logical. Defines whether there are parallel jobs or not.
#' @param refReadsThresh The threshold of proportion of non-zero sequencing reads for choosing the reference taxon in phase 2. Default is 0.2 meaning that at least 20% non-zero sequencing reads are necessary.
#' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. Default is 0 which means that taxon without any sequencing reads will be dropped from the analysis.
#' @param SDTresh The threshold of standard deviations of sequencing reads for being chosen as the reference taxon in phase 2. The default is 0.05 which means the standard deviation of sequencing reads should be at least 0.05 in order to be chosen as a reference taxon.
#' @param SDquantilThresh The threshold of the quantile of standard deviation of sequencing reads, above which could be selected as a reference taxon. Default is 0.
#' @param balanceCut The threshold of the proportion of non-zero sequencing reads in each group of a binary variable for choosing the final reference taxa in phase 2. The default is 0.2 meaning at least 20% non-zero sequencing reads in each group are needed to be eligible for being chosen as a final reference taxon.
#' @param verbose Whether the process message is printed out to the console. Default is TRUE.
#' @param seed Random seed for reproducibility. Can be set to NULL to remove seeding.
#' @return \code{microbiomeIFAADS} returns the association of the microbiome data with the covariates
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import IFAA
#' @export
#'


microbiomeIFAADS <- function(SumExp, microbVar, covariates, confounders, sampleIDname, covariatesMany, confoundersMany, nRef, nRefMaxforEsti, taxa, adjust_method,
                             fdrRate, paraJobs, bootB, standardize, sequentialRun, refReadsThresh, taxDropThresh, SDThresh, SDquantilThresh, balanceCut, verbose){

  SumExp <- eval(parse(text=SumExp), envir = parent.frame())


  # Computes the model

  outome <- IFAA::IFAA(experiment_dat = SumExp,
                       microbVar = microbVar,
                       testCov = covariates,
                       ctrlCov = confounders,
                       sampleIDname = sampleIDname,
                       testMany = covariatesMany,
                       ctrlMany = confoundersMany,
                       nRef = nRef,
                       nRefMaxForEsti = nRefMaxForEsti,
                       refTaxa = taxa,
                       adjust_method = adjust_method,
                       fdrRate = fdrRate,
                       paraJobs = paraJobs,
                       bootB = bootB,
                       standardize = standardize,
                       sequentialRun = sequentialRun,
                       refReadsThresh = refReadsThresh,
                       taxDropThresh = taxDropThresh,
                       SDThresh = SDThresh,
                       SDquantilThresh = SDquantilThresh,
                       balanceCut = balanceCut,
                       verbose = verbose)


  # the resulting model will be returned to the analyst
  return(outcome)


  Results.Taxon <- unlist(outcome[[1]][1])
  Results.Covariate <- unlist(outcome[[1]][2])
  Results.Estimate <- unlist(outcome[[1]][3])
  Results.SE.est <- unlist(outcome[[1]][4])
  Results.CI.low <- unlist(outcome[[1]][5])
  Results.CI.high <- unlist(outcome[[1]][6])
  Results.adj.p.value <- unlist(outcome[[1]][7])
  Results.unadjust.p.value <- unlist(outcome[[1]][8])
  Results.sig_ind <- unlist(outcome[[1]][9])
  Results.Final_Ref_Taxon <- unlist(outcome[[2]][2])


  Results.FDR <- fdrRate
  Results.Adjust_Method <- adjust_method
  Results.Boots <- bootB


  output <- data.frame(Results.Taxon,
                       Results.Covariate,
                       Results.Estimate,
                       Results.SE.est,
                       Results.CI.low,
                       Results.CI.high,
                       Results.adj.p.value,
                       Results.unadjust.p.value,
                       Results.sig_ind,
                       Results.FDR,
                       Results.Adjust_Method,
                       Results.Boots,
                       as.list(Results.Final_Ref_Taxon))


}

# aggregate function


