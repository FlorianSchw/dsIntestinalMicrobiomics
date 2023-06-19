#### Needs proper Documentation as an internal function






int.getScrResu.IFAA1 <- function(data,
                       testCovInd,
                       testCovInOrder,
                       testCovInNewNam,
                       nRef,
                       paraJobs,
                       refTaxa,
                       goodIndeCutPerc = 0.33,
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
  # run permutation


  message("Directly after int.getScrResu")


  scrParal <- int.runScrParal.IFAA1(
    data = data,
    testCovInd = testCovInd,
    testCovInOrder = testCovInOrder,
    testCovInNewNam = testCovInNewNam,
    nRef = nRef,
    paraJobs = paraJobs,
    refTaxa = refTaxa,
    sequentialRun = sequentialRun,
    allFunc = allFunc,
    refReadsThresh = refReadsThresh,
    SDThresh = SDThresh,
    SDquantilThresh = SDquantilThresh,
    balanceCut = balanceCut,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binPredInd,
    adjust_method = adjust_method
  )




    return(scrParal)
}
