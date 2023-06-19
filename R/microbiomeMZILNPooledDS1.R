#'
#' @title Computes the multivariate zero-inflated logistic normal model
#' @description This function calls the native R function from the IFAA package
#' @details The function computes a model from a SummerizedExperiment object with a given set of
#' microbiome taxa and covariates.
#' @param SumExp is a string character of the data.frame
#' @param microbVar This takes a single or vector of microbiome variable names (e.g., taxa, OTU and ASV names) of interest. Default is "all" meaning all microbiome variables will be analyzed. If a subset of microbiome variables is specified, the output will only contain the specified variables, and p-value adjustment for multiple testing will only be applied to the subset.
#' @param refTaxa is a string character for the microbiome variable denominator (can also be a vector of microbiome variables)
#' @param allCov is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param sampleIDname is a string character for the sample ID variable.
#' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for the p.adjust function in R.
#' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with 'covariates'.
#' @param paraJobs If 'sequentialRun' is FALSE, this specifies the number of parallel jobs that will be registered to run the algoithm. If specified as NULL, it will automatically detect the cores to decide the number of parallel jobs.
#' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates for the high dimensional regression.Default is 500.
#' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. Default is 0 which means that taxon without any sequencing reads will be dropped from the analysis.
#' @param standardize is a logical. If TRUE, the design matrix for X will be standardized in the analyses and the results. Default is FALSE.
#' @param sequentialRun is a logical. Default is TRUE. It can be set to be FALSE to increase speed if there are multiple taxa in the argument 'taxa'.
#' @param verbose Whether the process message is printed out to the console. Default is TRUE.
#' @return {microbiomeMZILNPooledDS1} returns the outcome of the specified multivariate zero-inflated logistic normal model
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @export
#' @import IFAA
#' @import doRNG
#'


microbiomeMZILNPooledDS1 <- function(SumExp, microbVar, refTaxa, allCov, sampleIDname, adjust_method, fdrRate, paraJobs, bootB, taxDropThresh, standardize, sequentialRun, verbose){

  experiment_dat <- eval(parse(text=SumExp), envir = parent.frame())



  allFunc <- int.allUserFunc()

  results <- list()

  start.time <- proc.time()[3]

  assay_name <- names(SummarizedExperiment::assays(experiment_dat))
  MicrobData <- data.frame(t(SummarizedExperiment::assays(experiment_dat)[[assay_name]]))

  MicrobData$ID_char_1234 <- rownames(MicrobData)
  CovData <- data.frame(SummarizedExperiment::colData(experiment_dat))
  CovData$ID_char_1234 <- rownames(CovData)
  linkIDname <- "ID_char_1234"

  if (verbose) {
    runMeta <- int.metaData(
      MicrobData = MicrobData,
      CovData = CovData,
      linkIDname = linkIDname,
      sampleIDname = sampleIDname,
      testCov = allCov,
      taxDropThresh = taxDropThresh,
      standardize = standardize
    )
  } else {
    runMeta <- suppressMessages(
      int.metaData(
        MicrobData = MicrobData,
        CovData = CovData,
        linkIDname = linkIDname,
        sampleIDname = sampleIDname,
        testCov = allCov,
        taxDropThresh = taxDropThresh,
        standardize = standardize
      )
    )
  }


  data <- runMeta$data
  covariatesData <- runMeta$covariatesData
  binaryInd <- runMeta$binaryInd
  covsPrefix <- runMeta$covsPrefix
  Mprefix <- runMeta$Mprefix
  testCovInd <- runMeta$testCovInd
  testCovInOrder <- runMeta$testCovInOrder
  testCovInNewNam <- runMeta$testCovInNewNam
  ctrlCov <- runMeta$ctrlCov
  microbName <- runMeta$microbName
  newMicrobNames <- runMeta$newMicrobNames
  results$covriateNames <- runMeta$xNames
  rm(runMeta)

  binaryInd_test <- testCovInd[testCovInd %in% binaryInd]

  nRef <- length(refTaxa)
  refTaxa_newNam <- newMicrobNames[microbName %in% refTaxa]


  if (length(refTaxa) > 0) {
    if (sum(refTaxa %in% microbName) != length(refTaxa)) {
      stop(
        "One or more of the specified reference taxa have no sequencing reads
      or are not in the data set. Double check the names of the reference taxa and their
           sparsity levels."
      )
    }
  }

  if (any(microbVar!="all")) {
    if (any(!(microbVar %in% microbName))) {
      stop("One or more taxon in microbVar is not available. Please double check")
    }
  }




  if (verbose) {
    results$analysisResults <- int.Regulariz_MZILN(
      data = data,
      nRef = nRef,
      sub_taxa = microbVar,
      testCovInd = testCovInd,
      testCovInOrder = testCovInOrder,
      testCovInNewNam = testCovInNewNam,
      microbName = microbName,
      binaryInd = binaryInd,
      binaryInd_test = binaryInd_test,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      refTaxa = refTaxa_newNam,
      adjust_method = adjust_method,
      paraJobs = paraJobs,
      fdrRate = fdrRate,
      bootB = bootB,
      sequentialRun = sequentialRun,
      allFunc = allFunc
    )
  } else {
    results$analysisResults <- suppressMessages(
      Regulariz_MZILN(
        data = data,
        nRef = nRef,
        sub_taxa = microbVar,
        testCovInd = testCovInd,
        testCovInOrder = testCovInOrder,
        testCovInNewNam = testCovInNewNam,
        microbName = microbName,
        binaryInd = binaryInd,
        binaryInd_test = binaryInd_test,
        covsPrefix = covsPrefix,
        Mprefix = Mprefix,
        refTaxa = refTaxa_newNam,
        adjust_method = adjust_method,
        paraJobs = paraJobs,
        fdrRate = fdrRate,
        bootB = bootB,
        sequentialRun = sequentialRun,
        allFunc = allFunc
      )
    )
  }



  results$linkIDname <- linkIDname


  #### the return object needs to be specified better

  return(results)



}

# aggregate function
