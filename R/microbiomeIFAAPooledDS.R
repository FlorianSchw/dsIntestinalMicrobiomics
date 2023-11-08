#'
#' @title Computes the association of microbiome data with covariates
#' @description This function calls the native R function from the IFAA package
#' @details The function computes an association from a SummarizedExperiment object with a given set of
#' confounders and covariates.
#' @param SumExp is a string character of the data.frame
#' @param microbVar_ds This takes a single or vector of microbiome variable names (e.g., taxa, OTU and ASV names) of interest. Default is "all" meaning all microbiome variables will be analyzed. If a subset of microbiome variables is specified, the output will only contain the specified variables, and p-value adjustment for multiple testing will only be applied to the subset.
#' @param testCov_ds is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param ctrlCov_ds is a string character for the covariates that will be adjusted in the model (can also be a vector of confounders)
#' @param sampleIDname is a string character for the sample ID variable.
#' @param testMany is a logical. If 'TRUE' and 'covariates' are set to NULL, then all variables in the 'covariates' will be used.
#' @param ctrlMany is a logical. If 'TRUE' and 'confounders' are set to NULL, then all variables except the 'coviariates' will be used as confounders.
#' @param refTaxa_ds vector of taxa or OTU or ASV names. Theses are reference taxa specified by the user to be used in phase 1.
#' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for the p.adjust function in R.
#' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with 'covariates'.
#' @param paraJobs If 'sequentialRun' is FALSE, this specifies the number of parallel jobs that will be registered to run the algoithm. If specified as NULL, it will automatically detect the cores to decide the number of parallel jobs.
#' @param bootB number of bootstrap samples for obtaining confidence interval of estimates in phase 2 for the high dimensional regression. Default is 500.
#' @param standardize is a logical. If 'TRUE', the design matrix for X will be standardized in the analyses and the results. Default is FALSE.
#' @param sequentialRun is a logical. Defines whether there are parallel jobs or not.
#' @param refReadsThresh The threshold of proportion of non-zero sequencing reads for choosing the reference taxon in phase 2. Default is 0.2 meaning that at least 20\% non-zero sequencing reads are necessary.
#' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. Default is 0 which means that taxon without any sequencing reads will be dropped from the analysis.
#' @param SDTresh The threshold of standard deviations of sequencing reads for being chosen as the reference taxon in phase 2. The default is 0.05 which means the standard deviation of sequencing reads should be at least 0.05 in order to be chosen as a reference taxon.
#' @param SDquantilThresh The threshold of the quantile of standard deviation of sequencing reads, above which could be selected as a reference taxon. Default is 0.
#' @param balanceCut The threshold of the proportion of non-zero sequencing reads in each group of a binary variable for choosing the final reference taxa in phase 2. The default is 0.2 meaning at least 20\% non-zero sequencing reads in each group are needed to be eligible for being chosen as a final reference taxon.
#' @param verbose Whether the process message is printed out to the console. Default is TRUE.
#' @return \code{microbiomeIFAADS} returns the association of the microbiome data with the covariates
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
                                   bootB,
                                   standardize,
                                   sequentialRun,
                                   refReadsThresh,
                                   taxDropThresh,
                                   SDThresh,
                                   SDquantilThresh,
                                   balanceCut,
                                   verbose){

  experiment_dat <- eval(parse(text=SumExp), envir = parent.frame())

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
      bootB = bootB,
      sequentialRun = sequentialRun,
      allFunc = allFunc,
      refReadsThresh = refReadsThresh,
      SDThresh = SDThresh,
      SDquantilThresh = SDquantilThresh,
      balanceCut = balanceCut
    )
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
        bootB = bootB,
        sequentialRun = sequentialRun,
        allFunc = allFunc,
        refReadsThresh = refReadsThresh,
        SDThresh = SDThresh,
        SDquantilThresh = SDquantilThresh,
        balanceCut = balanceCut
      )
    )
  }



  return(results)

}

# aggregate function


