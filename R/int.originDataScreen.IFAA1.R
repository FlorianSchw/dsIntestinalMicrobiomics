#### Needs proper Documentation as an internal function




int.originDataScreen.IFAA1 <- function(data,
                             testCovInd,
                             nRef,
                             paraJobs,
                             refTaxa,
                             maxDimensionScr = 434 * 5 * 10^5,
                             run_linear_thresh = 21121540,
                             sequentialRun,
                             allFunc,
                             Mprefix,
                             covsPrefix,
                             binPredInd,
                             adjust_method,
                             subsamp_cut = 3,
                             dfmax_div = 3,
                             nlambda_num = 10,
                             information_goodRef,
                             information_refTaxa) {
  results <- list()

  # load data info
  basicInfo <- int.dataInfo(
    data = data,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binPredInd
  )

  taxaNames <- basicInfo$taxaNames
  nTaxa <- basicInfo$nTaxa
  nPredics <- basicInfo$nPredics
  rm(basicInfo)
  gc()

  nNorm <- nTaxa - 1
  nAlphaNoInt <- nPredics * nNorm
  nAlphaSelec <- nPredics * nTaxa

  countOfSelec <- rep(0, nAlphaSelec)






  # overwrite nRef if the reference taxon is specified
  nRef <- length(refTaxa)

  startT1 <- proc.time()[3]
  if (length(paraJobs) == 0) {
    availCores <- parallelly::availableCores()
    if (is.numeric(availCores)) {
      paraJobs <- max(1, parallelly::availableCores() - 2)
    }
    if (!is.numeric(availCores)) {
      paraJobs <- 1
    }
  }




  if (!sequentialRun) {
    message(
      paraJobs,
      " parallel jobs are registered for the analysis."
    )
  }

  cl <- parallel::makeCluster(paraJobs)

  parallel::clusterExport(
    cl = cl,
    varlist = allFunc[c(1, 2, 3, 4)],
    envir = environment()
  )




  doParallel::registerDoParallel(cl)

  if (sequentialRun) {
    foreach::registerDoSEQ()
  }

  # start parallel computing


  scr1Resu <- foreach(
    i = seq_len(nRef),
    .multicombine = TRUE,
    .errorhandling = "pass"
  ) %dorng% {

    # for(i in seq_len(nRef)){

    ii <- which(taxaNames == refTaxa[i])

    dataForEst <- int.dataRecovTrans(
      data = data,
      ref = refTaxa[i],
      Mprefix = Mprefix,
      covsPrefix = covsPrefix,
      binPredInd = binPredInd,
      contCovStd = TRUE
    )

    xTildLongTild.i <- dataForEst$xTildalong
    yTildLongTild.i <- dataForEst$UtildaLong
    rm(dataForEst)
    gc()

    maxSubSamplSiz <- min(50000, floor(maxDimensionScr /
                                         ncol(xTildLongTild.i)))
    nToSamplFrom <- nrow(xTildLongTild.i)

    subSamplK <- ceiling(nToSamplFrom / maxSubSamplSiz)
    if (subSamplK == 1) {
      maxSubSamplSiz <- nToSamplFrom
    }

    nRuns <- ceiling(subSamplK / subsamp_cut)

    Penal_list <- list()
    for (k in seq_len(nRuns)){
      rowToKeep <- sample(nToSamplFrom, maxSubSamplSiz)
      x <- xTildLongTild.i[rowToKeep, ]
      y <- yTildLongTild.i[rowToKeep]

      if (nrow(x) * ncol(x) < run_linear_thresh) {
        Penal.i=int.runlinear.IFAA1(x=x,y=y, nPredics=nPredics)
        message("Testing here now.")

      }

      Penal_list[[k]] <- Penal.i

      Penal_list[[k]][[11]] <- nRuns
      Penal_list[[k]][[12]] <- nAlphaSelec
      Penal_list[[k]][[13]] <- nAlphaNoInt
      Penal_list[[k]][[14]] <- nTaxa
      Penal_list[[k]][[15]] <- ii
      Penal_list[[k]][[16]] <- testCovInd
      Penal_list[[k]][[17]] <- taxaNames
      Penal_list[[k]][[18]] <- information_goodRef
      Penal_list[[k]][[19]] <- information_refTaxa
    }
    return(Penal_list)

  }

}
