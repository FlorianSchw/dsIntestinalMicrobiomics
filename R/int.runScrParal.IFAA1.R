#### Needs proper Documentation as an internal function



int.runScrParal.IFAA1 <- function(data,
                        testCovInd,
                        testCovInOrder,
                        testCovInNewNam,
                        nRef,
                        paraJobs,
                        refTaxa,
                        maxDimensionScr = 0.8 * 434 * 10 * 10^4,
                        sequentialRun,
                        allFunc = allFunc,
                        refReadsThresh,
                        SDThresh,
                        SDquantilThresh,
                        balanceCut,
                        Mprefix,
                        covsPrefix,
                        binPredInd,
                        adjust_method) {
  results <- list()
  # load data info
  basicInfo <- int.dataInfo(
    data = data,
    qualifyRefTax = TRUE,
    refReadsThresh = refReadsThresh,
    SDThresh = SDThresh,
    SDquantilThresh = SDquantilThresh,
    balanceCut = balanceCut,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binPredInd
  )

  taxaNames <- basicInfo$taxaNames
  nTaxa <- basicInfo$nTaxa
  nPredics <- basicInfo$nPredics
  nSub <- basicInfo$nSub
  predNames <- basicInfo$predNames

  results$goodRefTaxaCandi <- basicInfo$goodRefTaxaCandi
  rm(basicInfo)
  gc()

  nNorm <- nTaxa - 1
  nAlphaNoInt <- nPredics * nNorm
  nAlphaSelec <- nPredics * nTaxa


  if (nRef > (length(taxaNames))) {
    nRef <- length(taxaNames)
  }

  # make reference taxa list
  if (length(refTaxa) < nRef) {
    taxon_to_be_sample <-
      (results$goodRefTaxaCandi[!(results$goodRefTaxaCandi %in% refTaxa)])
    num_to_be_sample <- (nRef - length(refTaxa))
    if (num_to_be_sample >= length(taxon_to_be_sample)) {
      num_to_be_sample <- length(taxon_to_be_sample)
    }

    refTaxa_extra <- sample(taxon_to_be_sample, num_to_be_sample)
    refTaxa <- c(refTaxa, refTaxa_extra)

    if (length(refTaxa) < nRef) {
      taxa_to_be_sample2 <- (taxaNames[!(taxaNames %in% refTaxa)])
      num_to_be_sample2 <- (nRef - length(refTaxa))
      refTaxa_extra2 <-
        sample(taxa_to_be_sample2, num_to_be_sample2)
      refTaxa <- c(refTaxa, refTaxa_extra2)
      message("Reference taxa are randomly picked in Phase 1.")
    }
  } else if (length(refTaxa) > nRef) {
    refTaxa <- sample(refTaxa, nRef)
  }

  results$refTaxa <- refTaxa





  #
  ## run original data screen
  #
  gc(reset = TRUE)

  if (suppressMessages(int.dataSparsCheck(data, Mprefix))[[1]] == 0) {

    screen1 <- int.originDataScreen_imputed(
      data = data,
      testCovInd = testCovInd,
      nRef = nRef,
      refTaxa = refTaxa,
      paraJobs = paraJobs,
      allFunc = allFunc,
      Mprefix = Mprefix,
      covsPrefix = covsPrefix,
      binPredInd = binPredInd,
      sequentialRun = sequentialRun,
      information_goodRef = results$goodRefTaxaCandi,
      information_refTaxa = results$refTaxa
    )
  } else {

    screen1 <- int.originDataScreen.IFAA1(
      data = data,
      testCovInd = testCovInd,
      nRef = nRef,
      refTaxa = refTaxa,
      paraJobs = paraJobs,
      allFunc = allFunc,
      Mprefix = Mprefix,
      covsPrefix = covsPrefix,
      binPredInd = binPredInd,
      sequentialRun = sequentialRun,
      adjust_method = adjust_method,
      information_goodRef = results$goodRefTaxaCandi,
      information_refTaxa = results$refTaxa
    )
  }




  return(screen1)
}


