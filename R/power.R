#' Compute required overall sample size for van Elteren test
#'
#' \code{power.vanElteren} uses the procedure suggested by Zhao et al. [1] to
#' compute the required overall sample size for a van Elteren test
#' comparing two groups A and B within strata [2].
#'
#' [1] Zhao et al., (2008), "Sample Size Calculation for the van Elteren Test
#'    Adjusting for Ties", Journal of Biopharmaceutical Statistics, 18:1112–1119.
#'
#' [2] van Elteren, (1960). "On the combination of independent two sample tests
#'    of Wilcoxon", Bulletin of the Institute of International Statistics, 37:351–361.
#'
#' @param strata_fractions numeric vector of length n_strata with fraction of
#'   subjects enrolled to the respective strata
#' @param treatment_fractions 2-by-n_strata matrix where the first row
#'   corresponds to group A and the second to group B. Each column must sum to 1
#'   and the k-th column holds the treatment fractions for the k-th stratum.
#' @param strata_fractions_alternative 2-by-n_categories-by-n_strata array. I.e.
#'   strata_fractions_alternative[g, , k] holds the expected category fractions
#'   in group g (= 1 or 2) for stratum k (= 1 ... n_strata);
#'   strata_fractions_alternative[g, , k] must sum to 1.
#' @param alpha maximal type one error rate
#' @param beta maximal type two error rate on specified point alternative
#'
#' @return \code{N}, single number giving the required overall sample size for the
#'   van Elteren test.
#'
#' @examples
#' # First example from [1], Table 3.
#' alpha       <- 0.05
#' beta        <- 0.2
#' strata_frac <- c(.78, .22)
#' n_strata    <- length(strata_frac)
#' treatment_frac   <- matrix(.5, nrow = 2, ncol = n_strata)
#' n_cat            <- 3
#' strata_frac_alt_BIFl28 <- matrix(c(
#'    0.62, 0.31, 0.07, # Duloxetine
#'    0.51, 0.39, 0.1   # Placebo
#' ), nrow = 2, ncol = n_cat, byrow = T)
#' strata_frac_alt_BIFgeq28 <- matrix(c(
#'    0.67, 0.33, 0.00, # Duloxetine
#'    0.27, 0.55, 0.18  # Placebo
#' ), nrow = 2, ncol = n_cat, byrow = T)
#' strata_frac_alt <- array(c(
#'     strata_frac_alt_BIFl28, strata_frac_alt_BIFgeq28
#' ), dim = c(2, n_cat, n_strata))
#' power.vanElteren(strata_frac, treatment_frac, strata_frac_alt, alpha, beta) # should be 228
#'
#' @export
power.vanElteren <- function(
    strata_fractions,
    treatment_fractions,
    strata_fractions_alternative,
    alpha,
    beta
) {
    # notes:

    n_strata <- length(strata_fractions)
    n_groups <- 2
    n_cat    <- dim(strata_fractions_alternative)[2]

    # test inputs
    assertthat::assert_that(all(strata_fractions >= 0))
    assertthat::assert_that(all(strata_fractions <= 1))
    assertthat::assert_that(sum(strata_fractions) == 1)

    assertthat::assert_that(all(dim(treatment_fractions) == c(n_groups, n_strata)))
    assertthat::assert_that(all(treatment_fractions >= 0))
    assertthat::assert_that(all(treatment_fractions <= 1))
    assertthat::assert_that(all(colSums(treatment_fractions) == 1))

    assertthat::assert_that(all(dim(strata_fractions_alternative)[c(1, 3)] == c(n_groups, n_strata)))
    assertthat::assert_that(all(strata_fractions_alternative >= 0))
    assertthat::assert_that(all(strata_fractions_alternative <= 1))
    for (k in 1:n_strata) {
        assertthat::assert_that(all(rowSums(strata_fractions_alternative[ , , k]) == 1))
    }

    mu0       <- .5
    w         <- as.numeric(rep(NA, n_strata)) # init as NA for easier debugging
    vsq0      <- as.numeric(rep(NA, n_strata))
    pi        <- as.numeric(rep(NA, n_strata))
    for (k in 1:n_strata) {
        s_k  <- strata_fractions[k]      # just renaming for easier
                                         # comparison with Zhao et al.
        t_k  <- treatment_fractions[2, k]    # fraction for group 2 in stratum k
        p_k  <- strata_fractions_alternative[1, , k] # category fractions for
                                                     # group 1 in stratum k
        q_k  <- strata_fractions_alternative[2, , k] # category fractions for
                                                     # group 2 in stratum k
        w[k] <-  s_k * t_k * (1 - t_k)      # still needs to be normalized after
                                            # loop!
        vsq0[k] <- ( 1 - sum(((1 - t_k)*p_k + t_k*q_k)^3) ) / (12 * s_k * t_k * (1 - t_k))
        pi[k]   <- .5 * p_k[1]*q_k[1] # only term for c == 1 (equation (15))
        for (c in 2:n_cat) {
            pi[k] <- pi[k] + p_k[c]*sum(q_k[1:(c - 1)]) + .5 * p_k[c]*q_k[c]
        }
        assertthat::assert_that(pi[k] >= 0)
        assertthat::assert_that(pi[k] <= 1)
    }
    w     <- w/sum(w)     # normalize weights!
    vsq0  <- sum(w^2*vsq0)
    mu1   <- sum(w*pi)
    N     <- (qnorm(alpha/2) + qnorm(beta))^2/(mu1 - mu0)^2 * vsq0
    return(N)
}