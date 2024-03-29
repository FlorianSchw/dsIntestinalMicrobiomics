#'
#' @title Computes the association of microbiome data with covariates
#' @description This function calls a custom version of the native R function IFAA from the IFAA package.
#' @details The function computes an association from a SummarizedExperiment object with a given set of
#' microbiome taxa, covariates and confounders. In contrast to the native R function, there is no phase1 for random refTaxa selection.
#' @param SumExp is a string character of the SummarizedExperiment object.
#' @param microbVar_ds is a single microbiome name or vector of microbiome variable names (e.g., taxa, OTU and ASV names) of interest. Default is "all" meaning all microbiome variables will be analyzed. If a subset of microbiome variables is specified, the output will only contain the specified variables, and p-value adjustment for multiple testing will only be applied to the subset.
#' @param refTaxa_ds vector of microbiome taxa or OTU or ASV names. Will be denominator(s) to the microbVar_ds.
#' @param testCov_ds is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param ctrlCov_ds is a string character for the covariates that will be adjusted in the model (can also be a vector of confounders).
#' @param sampleIDname is a string character for the sample ID variable.
#' @param testMany is a logical. If 'TRUE' and 'covariates' are set to NULL, then all variables in the 'covariates' will be used.
#' @param ctrlMany is a logical. If 'TRUE' and 'confounders' are set to NULL, then all variables except the 'coviariates' will be used as confounders.
#' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for the p.adjust function in R.
#' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with 'covariates'.
#' @param paraJobs If 'sequentialRun' is FALSE, this specifies the number of parallel jobs that will be registered to run the algoithm. If specified as NULL, it will automatically detect the cores to decide the number of parallel jobs.
#' @param standardize is a logical. If 'TRUE', the design matrix for X will be standardized in the analyses and the results. Default is FALSE.
#' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. Default is 0 which means that taxon without any sequencing reads will be dropped from the analysis.
#' @param verbose Whether the process message is printed out to the console. Default is TRUE.
#' @return \code{microbiomeIFAAPooledDS} returns a list consisting of intermediary results, mostly matrix crossproducts, from which estimates will be calculated on the client-side, and some additional information.
#' @importFrom dplyr %>%
#' @import dplyr
#' @import tidyr
#' @import dsBase
#' @import SummarizedExperiment
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @export
#'


microbiomeIFAAPooledDS <- function(SumExp,
                                   microbVar_ds,
                                   testCov_ds,
                                   ctrlCov_ds,
                                   sampleIDname,
                                   testMany,
                                   ctrlMany,
                                   refTaxa_ds,
                                   adjust_method,
                                   fdrRate,
                                   paraJobs,
                                   standardize,
                                   taxDropThresh,
                                   verbose){

  experiment_dat <- eval(parse(text=SumExp), envir = parent.frame())
  sequentialRun <- TRUE


  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)

  if(!(is.null(microbVar_ds))){
    microbVar <- unlist(strsplit(microbVar_ds, split=","))
  } else {
    microbVar <- microbVar_ds
  }

  if(!(is.null(testCov_ds))){
    testCov <- unlist(strsplit(testCov_ds, split=","))
  } else {
    testCov <- testCov_ds
  }

  if(!(is.null(ctrlCov_ds))){
    ctrlCov <- unlist(strsplit(ctrlCov_ds, split=","))
  } else {
    ctrlCov <- ctrlCov_ds
  }

  if(!(is.null(refTaxa_ds))){
    refTaxa <- unlist(strsplit(refTaxa_ds, split=","))
  } else {
    refTaxa <- refTaxa_ds
  }


  datashield_cov <- names(experiment_dat@colData@listData)
  datashield_microb <- experiment_dat@NAMES


  test_microbdata <- as.data.frame(t(experiment_dat@assays@data@listData[["MicrobiomeData"]])) %>%
    mutate(across(all_of(datashield_microb), ~replace_na(., 0))) %>%
    summarise(across(all_of(datashield_microb), ~sum(. == 0))) %>%
    pivot_longer(data = ., cols = everything(), names_to = "Variable", values_to = "Count")


  test_covariates <- as.data.frame(experiment_dat@colData@listData) %>%
    summarise(across(all_of(datashield_cov), ~sum(is.na(.)))) %>%
    pivot_longer(data = ., cols = everything(), names_to = "Variable", values_to = "Count")

  threshold_ds <- (dim(experiment_dat)[2] - nfilter.tab)


  test_combined <- rbind(test_microbdata,
                         test_covariates)

  test_combined_log <- c()
  for (rr in 1:length(test_combined$Variable)){

    test_combined_log[rr] <- test_combined$Count[rr] <= threshold_ds

  }

  unsuitable_varialbes <- paste0(test_combined$Variable[!test_combined_log], collapse = ", ")

  if(!(all(test_combined_log))){

    stop(paste0("Not all variables of interest pass the disclosure check: ", unsuitable_varialbes),  call.=FALSE)

  }




  allFunc <- int.allUserFunc()

  results <- list()
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
      testCov = testCov,
      ctrlCov = ctrlCov,
      testMany = testMany,
      ctrlMany = ctrlMany,
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
        testCov = testCov,
        ctrlCov = ctrlCov,
        testMany = testMany,
        ctrlMany = ctrlMany,
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

  binaryInd_test <- testCovInd[testCovInd %in% binaryInd]


  rm(runMeta)

  if (length(refTaxa) > 0) {
    if (length(unique(refTaxa)) != length(refTaxa)) {
      message("Duplicated names in refTaxa are removed")
      refTaxa <- unique(refTaxa)
    }

    if (sum(refTaxa %in% microbName) != length(refTaxa)) {
      refTaxa <- refTaxa[refTaxa %in% microbName]
      message(
        "One or more of the specified reference taxa in phase 1 have no sequencing reads
      or are not in the data set."
      )
    }
  }

  if (any(microbVar!="all")) {
    if (any(!(microbVar %in% microbName))) {
      stop("One or more taxon in microbVar is not available. Please double check")
    }
  }


  refTaxa_newNam <- newMicrobNames[microbName %in% refTaxa]

  if (verbose) {
    results$analysisResults <- int.Regulariz_IFAA(
      data = data,
      testCovInd = testCovInd,
      testCovInOrder = testCovInOrder,
      testCovInNewNam = testCovInNewNam,
      microbName = microbName,
      sub_taxa = microbVar,
      binaryInd = binaryInd,
      binaryInd_test = binaryInd_test,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      refTaxa_newNam = refTaxa_newNam,
      refTaxa = refTaxa,
      paraJobs = paraJobs,
      adjust_method = adjust_method,
      fwerRate = fdrRate,
      allFunc = allFunc,
      sequentialRun = sequentialRun)
  } else {
    results$analysisResults <- suppressMessages(
      int.Regulariz_IFAA(
        data = data,
        testCovInd = testCovInd,
        testCovInOrder = testCovInOrder,
        testCovInNewNam = testCovInNewNam,
        microbName = microbName,
        sub_taxa = microbVar,
        binaryInd = binaryInd,
        binaryInd_test = binaryInd_test,
        covsPrefix = covsPrefix,
        Mprefix = Mprefix,
        refTaxa_newNam = refTaxa_newNam,
        refTaxa = refTaxa,
        paraJobs = paraJobs,
        adjust_method = adjust_method,
        fwerRate = fdrRate,
        allFunc = allFunc,
        sequentialRun = sequentialRun)
    )
  }



  return(results)

}

# aggregate function


