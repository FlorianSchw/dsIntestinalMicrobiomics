############### Documentation needs to be adjusted towards betareg


#'
#' @title Fits beta regression models for rates and proportions via maximum likelihood
#' @description This function is similar to the native R function betareg from the betareg package
#' @details The function calls the server-side function \code{microbiomeBetaregDS} that computes the beta regression model
#' from data.frame. The function uses parametrization with mean and precision parameter.
#' @param df is a string character for the data.frame object
#' @param formula symbolic description of the model.
#' @param na.action is a string character specifying how NA data shall be treated.
#' @param weights optional numeric vector of case weights.
#' @param offset optional numeric vector with an a priori known component to be included in the linear predictor for the mean. In betareg.fit, offset may also be a list of two offsets for the mean and precision equation, respectively.
#' @param link character specification of the link function in the mean model (mu). Currently, "logit", "probit", "cloglog", "cauchit", "log", "loglog" are supported. Alternatively, an object of class "link-glm" can be supplied.
#' @param link.phi character specification of the link function in the precision model (phi). Currently, "identity", "log", "sqrt" are supported. The default is "log" unless formula is of type y ~ x where the default is "identity" (for backward compatibility). Alternatively, an object of class "link-glm" can be supplied.
#' @param type character specification of the type of estimator. Currently, maximum likelihood ("ML"), ML with bias correction ("BC"), and ML with bias reduction ("BR") are supported.
#' @param control a list of control arguments specified via betareg.control.
#' @param model is a logical. If TRUE, the model frame from the model fit will be returned.
#' @param x is a logical. If TRUE, the response from the model fit will be returned.
#' @param y is a logical. If TRUE, the model matrix from the model fit will be returned.
#' @return \code{ds.microbiomeBetareg} returns XXXXXXXXXXXXXXXXXXXXXX
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import betareg
#' @export
#'


microbiomeBetaregDS <- function(df, formula, na.action, weights, offset, link, link.phi, type, control, model, x, y){

  df <- eval(parse(text=df), envir = parent.frame())


  # Computes the model

  outome <- IFAA::MZILN(experiment_dat = SumExp,
                        refTaxa = taxa,
                        allCov = covariates,
                        sampleIDname = sampleIDname,
                        adjust_method = adjust_method,
                        fdrRate = fdrRate,
                        paraJobs = paraJobs,
                        bootB = bootB,
                        taxDropThresh = taxDropThresh,
                        standardize = standardize,
                        sequentialRun = sequentialRun,
                        verbose = verbose,
                        seed = seed)





  Results.Ref_Tax <- unlist(results_1[[1]][1])
  Results.Taxon <- unlist(results_1[[1]][2])
  Results.Covariate <- unlist(results_1[[1]][3])
  Results.Estimate <- unlist(results_1[[1]][4])
  Results.SE.est <- unlist(results_1[[1]][5])
  Results.CI.low <- unlist(results_1[[1]][6])
  Results.CI.high <- unlist(results_1[[1]][7])
  Results.adj.p.value <- unlist(results_1[[1]][8])
  Results.unadjust.p.value <- unlist(results_1[[1]][9])
  Results.sig_ind <- unlist(results_1[[1]][10])

  Results.FDR <- fdrRate
  Results.Adjust_Method <- adjust_method
  Results.Seed <- seed
  Results.Boots <- bootB



  output <- data.frame(Results.Ref_Tax, Results.Taxon, Results.Covariate, Results.Estimate, Results.SE.est, Results.CI.low, Results.CI.high, Results.adj.p.value, Results.unadjust.p.value, Results.sig_ind, Results.FDR, Results.Adjust_Method, Results.Seed, Results.Boots)

  # the resulting model will be returned to the analyst
  return(output)



}

# aggregate function
