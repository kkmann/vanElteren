context("Sample Size van Elteren")
testthat::test_that("First example from Zhao, Table 3.", {

alpha       <- 0.05
beta        <- 0.2
strata_frac <- c(.78, .22)
n_strata    <- length(strata_frac)
treatment_frac   <- matrix(.5, nrow = 2, ncol = n_strata)
n_cat            <- 3
strata_frac_alt_BIFl28 <- matrix(c(
   0.62, 0.31, 0.07, # Duloxetine
   0.51, 0.39, 0.1   # Placebo
), nrow = 2, ncol = n_cat, byrow = TRUE)
strata_frac_alt_BIFgeq28 <- matrix(c(
   0.67, 0.33, 0.00, # Duloxetine
   0.27, 0.55, 0.18  # Placebo
), nrow = 2, ncol = n_cat, byrow = TRUE)
strata_frac_alt <- array(c(
    strata_frac_alt_BIFl28, strata_frac_alt_BIFgeq28
), dim = c(2, n_cat, n_strata))
n <- samplesize.vanElteren(strata_frac, treatment_frac, strata_frac_alt, alpha, beta)

testthat::expect_equal(ceiling(n), 228)

})