#### Needs proper Documentation as an internal function

#' @import doRNG
#' @import foreach


int.bootResuHDCI <- function(data,
                             refTaxa,
                             originRefTaxNam,
                             maxDimension = 434 * 5 * 10^5,
                             bootLassoAlpha = 0.05,
                             binPredInd,
                             covsPrefix,
                             Mprefix,
                             allFunc,
                             unbalanceTaxa_ori_name,
                             unbalancePred_ori_name,
                             testCovInOrder,
                             adjust_method,
                             microbName,
                             fwerRate,
                             paraJobs,
                             taxa_sepname_list_arg,
                             subsamp_cut = 3){
  results <- list()
  basicInfo <- int.dataInfo(
    data = data,
    binPredInd = binPredInd,
    covsPrefix = covsPrefix,
    Mprefix = Mprefix)

  taxaNames <- basicInfo$taxaNames
  ii <- which(basicInfo$taxaNames %in% refTaxa)
  predNames <- basicInfo$predNames
  nTaxa <- basicInfo$nTaxa
  nPredics <- basicInfo$nPredics

  nNorm <- nTaxa - 1
  nAlphaNoInt <- nPredics * nNorm
  nAlphaSelec <- nPredics * nTaxa

  countOfSelec <- rep(0, nAlphaSelec)
  resultsByRefTaxon <- list()

  taxa_sepname_list_noref <- lapply(taxa_sepname_list_arg, function(x) x[x!=refTaxa])
  taxa_sepname_list<-lapply(taxa_sepname_list_noref, function(x) c(x,refTaxa))

  taxa_sepname_list_noref_vec <- unlist(taxa_sepname_list_noref)
  taxa_sepname_numeric <- as.numeric(sub(pattern = Mprefix,
                                         replacement = "",
                                         x = taxa_sepname_list_noref_vec))


  if (length(paraJobs) == 0) {
    availCores <- parallelly::availableCores()
    if (is.numeric(availCores)) {
      paraJobs <- max(1, parallelly::availableCores() - 2)
    }
    if (!is.numeric(availCores)) {
      paraJobs <- 1
    }
  }

  if (!sequentialRun){
    message(
      paraJobs,
      " parallel jobs are registered for analyzing reference taxa in Phase 2")
  }

  cl <- parallel::makeCluster(paraJobs)


  parallel::clusterExport(
    cl = cl,
    varlist = allFunc,
    envir = environment()
  )

  doParallel::registerDoParallel(cl)

  if (sequentialRun) {
    foreach::registerDoSEQ()
  }

  i <- numeric(0)

  phase2res <- foreach(
    i = seq_len(length(taxa_sepname_list)),
    .multicombine = TRUE,
    .errorhandling = "pass"
  ) %dorng% {
    data_sub_taxa <- data[, c(taxa_sepname_list[[i]], predNames)]
    dataForEst <- int.dataRecovTrans(
      data = data_sub_taxa,
      ref = refTaxa,
      Mprefix = Mprefix,
      covsPrefix = covsPrefix,
      binPredInd = binPredInd
    )
    x <- dataForEst$xTildalong
    y <- dataForEst$UtildaLong

    xCol <- ncol(x)

    maxSubSamplSiz <- min(100000, floor(maxDimension / xCol))
    nToSamplFrom <- length(y)
    subSamplK <- ceiling(nToSamplFrom / maxSubSamplSiz)
    if (subSamplK == 1) {
      maxSubSamplSiz <- nToSamplFrom
    }

    nRuns <- (ceiling(subSamplK / subsamp_cut))

    bootResu_k <- list()
      for (k in seq_len(nRuns)) {
        rowToKeep <- sample(nToSamplFrom, maxSubSamplSiz)
        xSub <- x[rowToKeep, ]
        ySub <- y[rowToKeep]
        bootResu <- int.lm_sparse(x = xSub, y = ySub)
        bootResu_k[[k]] <- bootResu
        bootResu_k[[k]][[9]] <- originRefTaxNam
        bootResu_k[[k]][[10]] <- testCovInOrder
        bootResu_k[[k]][[11]] <- nRuns
        bootResu_k[[k]][[12]] <- taxa_sepname_list
        bootResu_k[[k]][[13]] <- nPredics
        bootResu_k[[k]][[14]] <- i
        bootResu_k[[k]][[15]] <- microbName
        bootResu_k[[k]][[16]] <- unbalanceTaxa_ori_name
        bootResu_k[[k]][[17]] <- unbalancePred_ori_name
        bootResu_k[[k]][[18]] <- fwerRate
      }

    return(bootResu_k)

  }

}
