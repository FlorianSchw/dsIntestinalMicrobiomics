# Needs documentation





int.Regulariz.IFAA2 <- function(data,
                                testCovInd,
                                testCovInOrder,
                                testCovInNewNam,
                                microbName,
                                sub_taxa,
                                nRef,
                                nRefMaxForEsti,
                                refTaxa,
                                paraJobs,
                                binaryInd,
                                binaryInd_test,
                                covsPrefix,
                                Mprefix,
                                fwerRate,
                                bootB,
                                sequentialRun,
                                allFunc,
                                refReadsThresh,
                                SDThresh,
                                SDquantilThresh,
                                balanceCut,
                                adjust_method,
                                phase1_taxon_num = 200,
                                trans_x_col = 200,
                                spar_cutoff = 10,
                                nRef_smaller,
                                refTaxa_smaller) {

  message("nRef here?")
  results <- list()
  regul.start.time <- proc.time()[3]

  nTestCov <- length(testCovInd)

  int.dataSparsCheck(data = data, Mprefix = Mprefix)


  data.info <- int.dataInfo(
    data = data,
    qualifyRefTax = TRUE,
    refReadsThresh = refReadsThresh,
    SDThresh = SDThresh,
    SDquantilThresh = SDquantilThresh,
    balanceCut = balanceCut,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binaryInd
  )

  nSub <- data.info$nSub
  taxaNames <- data.info$taxaNames
  nPredics <- data.info$nPredics
  nTaxa <- data.info$nTaxa
  predNames <- data.info$predNames

  results$goodRefTaxaCandi <- data.info$goodRefTaxaCandi

  num_taxon_phase1 <- max(length(refTaxa), phase1_taxon_num)
  num_taxon_phase1 <- min(num_taxon_phase1, nTaxa)
  phase1_taxon_sample_pool <- results$goodRefTaxaCandi[!(results$goodRefTaxaCandi %in% refTaxa)]
  if (length(phase1_taxon_sample_pool) > (num_taxon_phase1 - length(refTaxa))) {
    phase1_taxon_sample1 <- c(sample(phase1_taxon_sample_pool,
                                     num_taxon_phase1 - length(refTaxa)),
                              refTaxa)
  } else {
    phase1_taxon_sample1 <- c(phase1_taxon_sample_pool, refTaxa)
  }
  phase1_taxon_sample <- phase1_taxon_sample1


  if (length(phase1_taxon_sample1) < num_taxon_phase1) {
    not_good_candi <- taxaNames[!((taxaNames %in% results$goodRefTaxaCandi) | (taxaNames %in% refTaxa))]
    non_zero_per_taxon <- colSums(data[, not_good_candi, drop = FALSE] > 0)
    phase1_taxon_sample2 <-
      names(sort(non_zero_per_taxon, decreasing = TRUE)[seq_len(num_taxon_phase1 -
                                                                  length(phase1_taxon_sample1))])
    phase1_taxon_sample <-
      c(phase1_taxon_sample1, phase1_taxon_sample2)
  }

  data_sub_phase1 <- data[, c(phase1_taxon_sample, predNames)]



  #### here the large nested int. functions from getScrResu start

  message("Directly before int.getScrResu")

  selectRegroup <- int.getScrResu.IFAA1(
    data = data_sub_phase1,
    testCovInd = testCovInd,
    testCovInOrder = testCovInOrder,
    testCovInNewNam = testCovInNewNam,
    nRef = nRef_smaller,
    paraJobs = paraJobs,
    refTaxa = refTaxa_smaller,
    sequentialRun = sequentialRun,
    allFunc = allFunc,
    refReadsThresh = refReadsThresh,
    SDThresh = SDThresh,
    SDquantilThresh = SDquantilThresh,
    balanceCut = balanceCut,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binaryInd,
    adjust_method = adjust_method
  )



  #### there is a while loop that was deleted and could be more tricky...



  return(selectRegroup)






}
