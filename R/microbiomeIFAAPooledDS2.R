#'
#' @export
#'


microbiomeIFAAPooledDS2 <- function(SumExp,
                                    microbVar,
                                    testCov,
                                    ctrlCov,
                                    sampleIDname,
                                    testMany,
                                    ctrlMany,
                                    nRef,
                                    nRefMaxForEsti,
                                    refTaxa,
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
                                    verbose,
                                    nRef_smaller,
                                    refTaxa_smaller){




  experiment_dat <- eval(parse(text=SumExp), envir = parent.frame())



  nRef_smaller <- as.numeric(nRef_smaller)
  refTaxa_smaller <- unlist(strsplit(refTaxa_smaller, split = ","))




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
  results$covriateNames <- runMeta$xNames

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


  if (nRefMaxForEsti < 2) {
    nRefMaxForEsti <- 2
    warning("Needs at least two final reference taxon for estimation.")
  }

  refTaxa_newNam <- newMicrobNames[microbName %in% refTaxa]

  if (verbose) {
    results$analysisResults <- int.Regulariz.IFAA2(
      data = data,
      testCovInd = testCovInd,
      testCovInOrder = testCovInOrder,
      testCovInNewNam = testCovInNewNam,
      microbName = microbName,
      sub_taxa = microbVar,
      nRef = nRef,
      nRefMaxForEsti = nRefMaxForEsti,
      binaryInd = binaryInd,
      binaryInd_test = binaryInd_test,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      refTaxa = refTaxa_newNam,
      paraJobs = paraJobs,
      adjust_method = adjust_method,
      fwerRate = fdrRate,
      bootB = bootB,
      sequentialRun = sequentialRun,
      allFunc = allFunc,
      refReadsThresh = refReadsThresh,
      SDThresh = SDThresh,
      SDquantilThresh = SDquantilThresh,
      balanceCut = balanceCut,
      nRef_smaller = nRef_smaller,
      refTaxa_smaller = refTaxa_smaller
    )
  } else {
    results$analysisResults <- suppressMessages(
      int.Regulariz.IFAA2(
        data = data,
        testCovInd = testCovInd,
        testCovInOrder = testCovInOrder,
        testCovInNewNam = testCovInNewNam,
        microbName = microbName,
        sub_taxa = microbVar,
        nRef = nRef,
        nRefMaxForEsti = nRefMaxForEsti,
        binaryInd = binaryInd,
        binaryInd_test = binaryInd_test,
        covsPrefix = covsPrefix,
        Mprefix = Mprefix,
        refTaxa = refTaxa_newNam,
        paraJobs = paraJobs,
        adjust_method = adjust_method,
        fwerRate = fdrRate,
        bootB = bootB,
        sequentialRun = sequentialRun,
        allFunc = allFunc,
        refReadsThresh = refReadsThresh,
        SDThresh = SDThresh,
        SDquantilThresh = SDquantilThresh,
        balanceCut = balanceCut,
        nRef_smaller = nRef_smaller,
        refTaxa_smaller = refTaxa_smaller
      )
    )
  }





  return(results)







}

# aggregate function


