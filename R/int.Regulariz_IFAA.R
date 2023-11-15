#' @details internal function for the microbiomeIFAAPooledDS function.






int.Regulariz_IFAA <- function(data,
                      testCovInd,
                      testCovInOrder,
                      testCovInNewNam,
                      microbName,
                      sub_taxa,
                      refTaxa,
                      refTaxa_newNam,
                      paraJobs,
                      binaryInd,
                      binaryInd_test,
                      covsPrefix,
                      Mprefix,
                      fwerRate,
                      allFunc,
                      adjust_method,
                      sequentialRun,
                      phase1_taxon_num = 200,
                      trans_x_col = 200,
                      spar_cutoff = 10) {
  results <- list()

  nTestCov <- length(testCovInd)
  int.dataSparsCheck(data = data, Mprefix = Mprefix)


  data.info <- int.dataInfo(
    data = data,
    qualifyRefTax = FALSE,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binaryInd)

  nSub <- data.info$nSub
  taxaNames <- data.info$taxaNames
  nPredics <- data.info$nPredics
  nTaxa <- data.info$nTaxa
  predNames <- data.info$predNames


  qualifyData <- data


  if (length(binaryInd_test) > 0) {
    qualifyData <- data[rowSums(data[, taxaNames] > 0) >= 2, , drop = FALSE]
    allBinPred <- paste0(covsPrefix, binaryInd_test)
    nBinPred <- length(allBinPred)

    # find the pairs of binary preds and taxa for which the associaton is not identifiable
    AllTaxaNamesNoRefTax <- taxaNames[!microbName %in% refTaxa]
    unbalanceTaxa <- c()
    unbalancePred <- c()
    for (i in AllTaxaNamesNoRefTax) {
      for (j in allBinPred) {
        twoColumns.ij <- qualifyData[, c(i, j)]
        nNonZero <- sum(twoColumns.ij[, 1] > 0)
        sumOfBin <- sum(twoColumns.ij[(twoColumns.ij[, 1] > 0), 2])
        if (sumOfBin %in% c(0, 1, (nNonZero - 1), nNonZero)) {
          unbalanceTaxa <- c(unbalanceTaxa, i)
          unbalancePred <- c(unbalancePred, j)
        }
      }
    }
    if (length(unbalanceTaxa) > 0) {
      unbalanceTaxa_ori_name <- microbName[unlist(lapply(unbalanceTaxa, function(x) {
        which(taxaNames %in% x)
      }))]
      unbalancePred_ori_name <- testCovInOrder[unlist(lapply(unbalancePred, function(x) {
        which(testCovInNewNam %in% x)
      }))]
    } else {
      unbalanceTaxa_ori_name <- NULL
      unbalancePred_ori_name <- NULL
    }
    rm(allBinPred)
  } else {
    unbalanceTaxa_ori_name <- NULL
    unbalancePred_ori_name <- NULL
  }

  ##### Partition in phase 2 ####
  num_taxa_each <- ceiling(trans_x_col / (nPredics + 1))

  rowSpars <- apply(data[, taxaNames], 1, function(x) {
    sum(x == 0) / length(x)
  })

  meadianRowSpars <- min(median(rowSpars), 0.999)

  num_taxa_each <- max(ceiling(spar_cutoff / (1 - meadianRowSpars)), num_taxa_each)
  num_taxa_each <- min(num_taxa_each, nTaxa)

  spar_each_taxon <- apply(data[, taxaNames], 2, function(x) {
    sum(x == 0) / length(x)
  })

  high_spar_taxon <- spar_each_taxon[spar_each_taxon > 0.7]
  low_spar_taxon <- spar_each_taxon[spar_each_taxon <= 0.7]

  num_cut <- floor((nTaxa) / num_taxa_each)

  high_spar_cut_taxon <- suppressWarnings(split(high_spar_taxon, seq_len(num_cut)))
  low_spar_cut_taxon <- suppressWarnings(split(low_spar_taxon, seq_len(num_cut)))

  taxa_sepname_list <- list()
  for (i in seq_len(num_cut)) {
    taxa_sepname_list[[i]] <- c(names(low_spar_cut_taxon[[i]]), names(high_spar_cut_taxon[[num_cut + 1 - i]]))
  }

  results$estiList <- list()

  for (iii in seq_len(length(refTaxa))){
    originTaxNam <- refTaxa[iii]
    newRefTaxNam <- taxaNames[microbName %in% originTaxNam]
    results$estiList[[originTaxNam]] <- int.bootResuHDCI(
      data = data,
      refTaxa = newRefTaxNam,
      originRefTaxNam = originTaxNam,
      binPredInd = binaryInd,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      allFunc = allFunc,
      unbalanceTaxa_ori_name = unbalanceTaxa_ori_name,
      unbalancePred_ori_name = unbalancePred_ori_name,
      testCovInOrder = testCovInOrder,
      adjust_method = adjust_method,
      microbName = microbName,
      fwerRate = fwerRate,
      paraJobs = paraJobs,
      taxa_sepname_list_arg = taxa_sepname_list,
      sequentialRun = sequentialRun)
  }


  results$nSub <- nSub
  results$nTaxa <- nTaxa
  results$qualifyData <- dim(qualifyData)[1]

  return(results)

}
