










test_that("summarizedExperimentDS errors", {

  # Creating differing test data.frames -  preparation
  id <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
  age <- c(71,45,44,68,54,31,58,65,32,41,31,50,52,32,40,35,64,50,59,50,65,64)
  sex <- c(1,1,1,1,0,0,0,0,1,1,0,1,0,1,1,0,1,0,0,0,1,0)
  microbiome1 <- c(36.67, 21.95, 33.64, 22.93, 12.07, 28.07, 37.7, 41.94, 14.23, 33.33, 55.4, 35.21, 25.54, 14.59, 5.16, 30.57, 38.1, 36.55, 49.34, 16.67, 27.02, 11.11)
  microbiome2 <- c(33.33, 15.24, 10.91, 36.1, 24.14, 32.89, 19.67, 9.68, 28.45, 17.02, 14.39, 16.43, 37.5, 14.59, 49.03, 31.61, 22.22, 13.2, 33.55, 3.13, 12.71, 26.8)
  microbiome3 <- c(10.0, 34.76, 29.09, 15.12, 60.34, 30.26, 31.97, 8.6, 30.96, 12.77, 7.19, 29.11, 12.5, 38.92, 32.26, 30.05, 1.59, 26.4, 0.0, 40.1, 19.34, 41.83)
  microbiome4 <- c(20.0, 28.05, 26.36, 25.85, 3.45, 8.77, 10.66, 39.78, 26.36, 36.88, 23.02, 19.25, 24.46, 31.89, 13.55, 7.77, 38.1, 23.86, 17.11, 40.1, 30.94, 20.26)

  microbiome_df1 <- data.frame(id, age, sex, microbiome1, microbiome2, microbiome3, microbiome4)

  # Creating results of distDS to compare to
  res1 <- summarizedExperimentDS(df = "microbiome_df1", microbiomeData = c("microbiome1", "microbiome2"), covariateData = c("age", "sex"))
  res2 <- summarizedExperimentDS(df = "microbiome_df1", microbiomeData = "microbiome1", covariateData = c("age", "sex"))


  ################################################################################
  ################################################################################
  ########################### Needs to be changed from here ######################
  ################################################################################
  ################################################################################



  outcome_correct1 <- matrix(c(0.00000, 32.18695, 53.88877, 44.05678, 22.13594,
                               32.18695, 0.00000, 32.80244, 27.29469, 37.60319,
                               53.88877, 32.80244, 0.00000, 17.57840, 47.81213,
                               44.05678, 27.29469, 17.57840, 0.00000, 34.04409,
                               22.13594, 37.60319, 47.81213, 34.04409, 0.00000),
                             nrow = 5, ncol = 5, byrow = TRUE,
                             dimnames = list(c("1", "2", "3", "4", "5"),
                                             c("1", "2", "3", "4", "5")))

  outcome_correct2 <- matrix(c(0.00000, 33.42155, 62.22540, 45.92385, 25.56039,
                               33.42155, 0.00000, 37.87699, 27.56810, 43.42043,
                               62.22540, 37.87699, 0.00000, 20.26491, 55.20869,
                               45.92385, 27.56810, 20.26491, 0.00000, 39.29377,
                               25.56039, 43.42043, 55.20869, 39.29377, 0.00000),
                             nrow = 5, ncol = 5, byrow = TRUE,
                             dimnames = list(c("1", "2", "3", "4", "5"),
                                             c("1", "2", "3", "4", "5")))


  # Actual Test Start
  expect_equal(attr(res1, "Size"), 5)
  expect_equal(attr(res1, "Diag"), FALSE)
  expect_equal(attr(res1, "Upper"), FALSE)
  expect_equal(attr(res1, "method"), "euclidean")
  expect_equal(class(res1), "dist")

  expect_equal(attr(res2, "Size"), 5)
  expect_equal(attr(res2, "Diag"), FALSE)
  expect_equal(attr(res2, "Upper"), FALSE)
  expect_equal(attr(res2, "method"), "euclidean")
  expect_equal(class(res2), "dist")

  expect_equal(m1, outcome_correct1, tolerance = 1e-2)
  expect_equal(m2, outcome_correct2, tolerance = 1e-2)

})



