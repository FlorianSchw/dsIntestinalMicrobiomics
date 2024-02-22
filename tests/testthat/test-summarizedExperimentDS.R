test_that("summarizedExperimentDS errors", {

  # Creating differing test data.frames -  preparation
  Sex <- as.factor(c(0,1,1,1,0,0,0,0,0,1,1,1,0,1,0))
  Weight <- c(56,67,65,45,78,87,88,66,65,77,61,49,52,56,54)
  Height <- c(168,176,164,189,201,167,172,184,178,192,174,176,182,183,157)
  Microb1 <- c(0.87,0.83,0.19,0.34,0.56,0.44,0.89,0.87,0.92,0.56,0.22,0.45,0.68,0.91,0.14)
  Microb2 <- c(0.13,0.17,0.81,0.66,0.44,0.56,0.11,0.13,0.08,0.44,0.78,0.55,0.32,0.09,0.86)

  test_df1 <- data.frame(Sex,
                         Weight,
                         Height,
                         Microb1,
                         Microb2)


  microbiomeData <- subset(test_df1, select = c(4,5))
  covariateData <- subset(test_df1, select = c(1,2,3))

  # Creating results of distDS to compare to
  res1 <- summarizedExperimentDS(microbiomeData = "microbiomeData",
                                 covariateData = "covariateData")

  # Actual Test Start

  expect_equal("SummarizedExperiment", class(res1)[1])
  expect_equal(15, res1@colData@nrows)
  expect_silent(res1)



})
