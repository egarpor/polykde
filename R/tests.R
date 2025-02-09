
#' @title Homogeneity test for several polyspherical samples
#'
#' @description Permutation tests for the equality of distributions of two or
#' \eqn{k} samples of data on \eqn{\mathcal{S}^{d_1} \times \cdots \times
#' \mathcal{S}^{d_k}}. The Jensen--Shannon distance is used to construct a test
#' statistic measuring the discrepancy between the \eqn{k} kernel density
#' estimators. Tests based on the mean and scatter matrices are also available,
#' but for only two samples (\eqn{k=2}).
#'
#' @inheritParams kde_polysph
#' @param labels vector with \code{k} different levels indicating the group.
#' @param type kind of test to be performed: \code{"mean"}, a simple test
#' for the equality of both means (non-omnibus for testing homogeneity);
#' \code{"scatter"}, a simple test for the equality of both scatter matrices
#' (non-omnibus for testing homogeneity); \code{"jsd"} (default), a test
#' comparing the kernel density estimators for \eqn{k} groups using the
#' Jensen--Shannon distance.
#' @param B number of permutations to use. Defaults to \code{1e3}.
#' @param M number of Monte Carlo replicates to use when approximating the
#' Hellinger/Jensen--Shannon distance. Defaults to \code{1e4}.
#' @param cv_jsd use cross-validation to approximate the Jensen--Shannon
#' distance? Does not require Monte Carlo. Defaults to \code{TRUE}.
#' @param plot_boot flag to display a graphical output of the test decision.
#' Defaults to \code{FALSE}.
#' @param seed_jsd seed for the Monte Carlo simulations used to estimate the
#' integrals in the Jensen--Shannon distance.
#' @details Only \code{type = "jsd"} is able to deal with \eqn{k > 2}.
#'
#' The \code{"mean"} statistic measures the maximum (chordal) distance
#' between the estimated group means. This statistic is bounded in \eqn{[0, 1]}.
#' The \code{"var"} statistic measures the maximum affine invariant Riemannian
#' metric between the estimated scatter matrices. The \code{"jsd"} statistic is
#' the Jensen--Shannon divergence. This statistic is bounded in \eqn{[0, 1]}.
#' The \code{"hd"} statistic computes a monotonic transformation of the
#' Hellinger distance, which is the Bhattacharyya divergence (or coefficient).
#' @return An object of class \code{"htest"}.
#' @examples
#' ## Two-sample case
#' \donttest{
#' # H0 holds
#' n <- c(50, 100)
#' X1 <- rotasym::r_vMF(n = n[1], mu = c(0, 0, 1), kappa = 1)
#' X2 <- rotasym::r_vMF(n = n[2], mu = c(0, 0, 1), kappa = 1)
#' hom_test_poly(X = rbind(X1, X2), labels = rep(1:2, times = n),
#'               d = 2, type = "jsd", h = 0.5)
#'
#' # H0 does not hold
#' X2 <- rotasym::r_vMF(n = n[2], mu = c(0, 1, 0), kappa = 2)
#' hom_test_poly(X = rbind(X1, X2), labels = rep(1:2, times = n),
#'               d = 2, type = "jsd", h = 0.5)
#'
#' ## k-sample case
#'
#' # H0 holds
#' n <- c(50, 100, 50)
#' X1 <- rotasym::r_vMF(n = n[1], mu = c(0, 0, 1), kappa = 1)
#' X2 <- rotasym::r_vMF(n = n[2], mu = c(0, 0, 1), kappa = 1)
#' X3 <- rotasym::r_vMF(n = n[3], mu = c(0, 0, 1), kappa = 1)
#' hom_test_poly(X = rbind(X1, X2, X3), labels = rep(1:3, times = n),
#'               d = 2, type = "jsd", h = 0.5)
#'
#' # H0 does not hold
#' X3 <- rotasym::r_vMF(n = n[3], mu = c(0, 1, 0), kappa = 2)
#' hom_test_poly(X = rbind(X1, X2, X3), labels = rep(1:3, times = n),
#'               d = 2, type = "jsd", h = 0.5)
#' }
#' @export
hom_test_poly <- function(X, d, labels,
                          type = c("mean", "scatter", "jsd", "hd")[3],
                          h = NULL, kernel = 1, kernel_type = 1, k = 10,
                          B = 1e3, M = 1e4, plot_boot = FALSE, seed_jsd = NULL,
                          cv_jsd = TRUE) {

  # Dimensions and sample sizes
  r <- length(d)
  N <- nrow(X)

  # Split X into kg groups
  n <- table(labels)
  stopifnot(all(n > 1))
  labels_levels <- names(n)
  n <- unname(n)
  kg <- length(labels_levels)

  # Checks
  stopifnot(N == length(labels))
  stopifnot(ncol(X) == sum(d + 1))

  # Type of test
  if (type == "mean") {

    # Only two groups
    stopifnot(kg == 2)

    # Fill bandwidths
    h <- NULL

    # Index for accessing each S^dj with ind[j]:(ind[j + 1] - 1)
    ind <- cumsum(c(1, d + 1))

    # Test information
    method <- "Permutation-based test of means equality"
    alternative <- ifelse(r == 1, "both group means are different",
                          "at least two group means are different")

    # Test statistic function
    Tn <- function(perm_index) {

      # Permute labels by perm_index
      perm_labels <- labels[perm_index]

      # Two samples
      X1 <- X[perm_labels == labels_levels[1], , drop = FALSE]
      X2 <- X[perm_labels == labels_levels[2], , drop = FALSE]

      # Discrepancies
      discs <- sapply(seq_len(r), function(j) {

        # Access to S^dj
        ind_j <- ind[j]:(ind[j + 1] - 1)

        # Estimated means
        mu_hat_1 <- DirStats::mu_ml(data = X1[, ind_j])
        mu_hat_2 <- DirStats::mu_ml(data = X2[, ind_j])

        # Chordal distance (monotonic relation to its intrinsic counterpart)
        return(0.5 * (1 - sum(mu_hat_1 * mu_hat_2)))

      })

      # Maximal discrepancy
      max(discs)

    }

  } else if (type == "scatter") {

    # Only two groups
    stopifnot(kg == 2)

    # Fill bandwidths
    h <- NULL

    # Index for accessing each S^dj with ind[j]:(ind[j + 1] - 1)
    ind <- cumsum(c(1, d + 1))

    # Test information
    method <- "Permutation-based test of scatter matrices equality"
    alternative <- ifelse(r == 1, "both group scatter matrices are different",
                          "at least two group scatter matrices are different")

    # Test statistic function
    Tn <- function(perm_index) {

      # Permute labels by perm_index
      perm_labels <- labels[perm_index]

      # Two samples
      X1 <- X[perm_labels == labels_levels[1], , drop = FALSE]
      X2 <- X[perm_labels == labels_levels[2], , drop = FALSE]

      # Discrepancies
      discs <- sapply(seq_len(r), function(j) {

        # Access to S^dj
        ind_j <- ind[j]:(ind[j + 1] - 1)

        # Estimated scatters
        S_1 <- crossprod(X1[, ind_j]) / n[1]
        S_2 <- crossprod(X2[, ind_j]) / n[2]

        # Affine invariant Riemannian metric between SPD matrices
        return(sqrt(sum(log(eigen(x = solve(S_1, S_2), symmetric = FALSE,
                                  only.values = TRUE)$values)^2)))

      })

      # Maximal discrepancy
      max(discs)

    }

  } else {

    # Set bandwidths
    if (is.null(h)) {

      h <- bw_rot_polysph(X = X, d = d, kernel = kernel,
                          kernel_type = kernel_type, k = k)$bw

    } else {

      stopifnot(length(h) == r)

    }

    # Kernel-based tests
    if (type == "jsd") {

      # Test information
      method <- "Permutation-based Jensen--Shannon distance test of homogeneity"
      alternative <- "any alternative to homogeneity"

      # Estimated probabilities
      probs <- n / N

      # Test statistic function
      Tn <- function(perm_index) {

        # The Jensen--Shannon Divergence can be conveniently written as
        # JSD(f_1, ..., f_k) = H(f_0) - \sum_j \pi_j * H(f_j) where
        # H(f_j) = -\int \log(f_j(x)) f_j(x) is the Shannon entropy of f_j
        # and f_0(x) = \sum_j \pi_j * f_j(x)
        # Therefore, to approximate JSD(f_1, ..., f_k) we can:
        # (1) Obtain Monte Carlo samples from f_j
        # (2) Obtain a Monte Carlo sample from f_0 by recycling those extracted
        #     in 1) and respecting the proportion pi_j of samples from f_j in
        #     the final sample of f_0.
        # (3) Approximate H(f_j) by Monte Carlo, evaluating \log(f_j).
        # (4) Approximate H(f_0) by Monte Carlo, evaluating
        #     \log(f_0) = log_f0_max + \sum(\exp(log_f0 - log_f0_max)).

        # Alternatively to (1)-(4), replace the Monte Carlo samples by the
        # observed samples, using log_cv_kde_polysph(). Therefore:
        # (1) Approximate H(f_j) = mean(log_cv_kde_polysph(X = X_j, h = h_j))
        # (2) Approximate H(f_0)

        # Permute labels by perm_index
        perm_labels <- labels[perm_index]

        # Access groups
        ind_j <- lapply(1:kg, function(j)
          which(perm_labels == labels_levels[j]))

        # CV approximation of the H(f_j)'s and H(f_0)?
        if (cv_jsd == 1) {

          # # Estimate H(f_j)'s
          # mc_samp <- list()
          # H_fj <- numeric(kg)
          # for (j in seq_len(kg)) {
          #
          #   # H(f_j)
          #   log_fj_cv <- log_cv_kde_polysph(X = X[ind_j[[j]], ], d = d,
          #                                   h = h, kernel = kernel,
          #                                   kernel_type = kernel_type, k = k,
          #                                   wrt_unif = TRUE)
          #   H_fj[j] <- -mean(log_fj_cv)
          #
          # }
          #
          # # H(f_0)
          # log_f0_cv <- log_cv_kde_polysph(X = X, d = d, h = h, kernel = kernel,
          #                                 kernel_type = kernel_type, k = k,
          #                                 weights = rep(probs / n, times = n),
          #                                 wrt_unif = TRUE)
          # H_f0 <- -mean(log_f0_cv)
          #
          # # Finally, compute Jensen--Shannon divergence
          # jsd_old <- H_f0 - sum(probs * H_fj)

          # Compute pi_j * log_fj_cv and concatenate them
          mc_samp <- list()
          log_fj_cv <- numeric()
          for (j in seq_len(kg)) {

            log_fj_cv <- c(log_fj_cv,
                           log_cv_kde_polysph(X = X[ind_j[[j]], ], d = d,
                                              h = h, kernel = kernel,
                                              kernel_type = kernel_type,
                                              k = k, wrt_unif = TRUE))

          }
          log_fj_cv <- log_fj_cv / N

          # Compute log_f0_cv
          log_f0_cv <- log_cv_kde_polysph(X = X, d = d, h = h, kernel = kernel,
                                          kernel_type = kernel_type, k = k,
                                          weights = rep(probs / n, times = n),
                                          wrt_unif = TRUE)
          log_f0_cv <- log_f0_cv / N

          # H_f0 - sum(probs * H_fj), removing Infs and NaNs from the sum
          H_f0_fj_dif <- -log_f0_cv + log_fj_cv
          H_f0_fj_dif[!is.finite(H_f0_fj_dif)] <- NA
          jsd <- sum(H_f0_fj_dif, na.rm = TRUE)

        } else if (cv_jsd == 2) {

          # Estimate H(f_j)'s
          H_fj <- numeric(kg)
          for (j in seq_len(kg)) {

            # H(f_j)
            log_fj <- kde_polysph(x = X[ind_j[[j]], ], X = X[ind_j[[j]], ],
                                  d = d, h = h, kernel = kernel,
                                  kernel_type = kernel_type, k = k, log = TRUE,
                                  wrt_unif = TRUE)
            H_fj[j] <- -mean(log_fj)

          }

          # # log_f0 for H(f_0)
          # log_f0 <- matrix(nrow = N, ncol = kg)
          # for (j in seq_len(kg)) {
          #
          #   log_f0[, j] <- kde_polysph(x = X, X = X[ind_j[[j]], ],
          #                              d = d, h = h, kernel = kernel,
          #                              kernel_type = kernel_type, k = k,
          #                              log = TRUE, wrt_unif = TRUE) +
          #     log(probs[j])
          #
          # }
          # log_max <- apply(log_f0, 1, max)
          # log_f0 <- log_max + log(rowSums(exp(log_f0 - log_max)))
          # H_f0 <- -mean(log_f0)

          # H(f_0)
          log_f0 <- kde_polysph(x = X, X = X, d = d, h = h, kernel = kernel,
                                kernel_type = kernel_type, k = k,
                                weights = rep(probs / n, times = n),
                                log = TRUE, wrt_unif = TRUE)
          H_f0 <- -mean(log_f0)

          # Finally, compute Jensen--Shannon divergence
          jsd <- H_f0 - sum(probs * H_fj)

        } else {

          # Estimate H(f_j)'s
          mc_samp <- list()
          H_fj <- numeric(kg)
          for (j in seq_len(kg)) {

            # j-th group Monte Carlo sample
            mc_samp[[j]] <- r_kde_polysph(n = M, X = X[ind_j[[j]], ], d = d,
                                          h = h, kernel = kernel,
                                          kernel_type = kernel_type, k = k)

            # H(f_j)
            log_fj <- kde_polysph(x = mc_samp[[j]], X = X[ind_j[[j]], ], d = d,
                                  h = h, kernel = kernel,
                                  kernel_type = kernel_type, k = k, log = TRUE,
                                  wrt_unif = TRUE)
            H_fj[j] <- -mean(log_fj)

          }

          # Pooled Monte Carlo sample
          M_j <- round(M * probs)
          mc_samp_0 <- matrix(nrow = sum(M_j), ncol = ncol(mc_samp[[1]]))
          ind_0 <- c(0, cumsum(M_j))
          for (j in seq_len(kg)) {

            mc_samp_0[(ind_0[j] + 1):ind_0[j + 1], ] <-
              mc_samp[[j]][1:M_j[j], , drop = FALSE]

          }

          # log_f0 for H(f_0)
          log_f0 <- matrix(nrow = nrow(mc_samp_0), ncol = kg)
          for (j in seq_len(kg)) {

            log_f0[, j] <- kde_polysph(x = mc_samp_0, X = X[ind_j[[j]], ],
                                       d = d, h = h, kernel = kernel,
                                       kernel_type = kernel_type, k = k,
                                       log = TRUE, wrt_unif = TRUE) +
              log(probs[j])

          }
          log_max <- apply(log_f0, 1, max)
          log_f0 <- log_max + log(rowSums(exp(log_f0 - log_max)))
          H_f0 <- -mean(log_f0)

          # Finally, compute Jensen--Shannon divergence
          jsd <- H_f0 - sum(probs * H_fj)

        }

        return(jsd)

      }

    } else if (type == "hd") {

      # Only two groups
      stopifnot(kg == 2)

      # Test information
      method <- "Permutation-based Hellinger distance test of homogeneity"
      alternative <- "any alternative to homogeneity"

      # Common Monte Carlo sample for integrals
      mc_samp <- r_unif_polysph(n = M, d = d)

      # Test statistic function
      Tn <- function(perm_index) {

        # Permute labels by perm_index
        perm_labels <- labels[perm_index]

        # Log kdes
        log_kde_1 <- kde_polysph(x = mc_samp, X = X[perm_labels, ], d = d,
                                 h = h, kernel = kernel,
                                 kernel_type = kernel_type, k = k, log = TRUE,
                                 wrt_unif = TRUE)
        log_kde_2 <- kde_polysph(x = mc_samp, X = X[!perm_labels, ], d = d,
                                 h = h, kernel = kernel,
                                 kernel_type = kernel_type, k = k, log = TRUE,
                                 wrt_unif = TRUE)

        # Distance
        hd_mc(log_f = log_kde_1, log_g = log_kde_2, bhatta = TRUE)

      }

    } else {

      stop("\"type\" must be \"means\", \"hd\" or \"jsd\".")

    }

  }

  # Original statistic
  if (!is.null(seed_jsd) && type == "jsd") {

    set.seed(seed_jsd, kind = "Mersenne-Twister")

  }
  Tn_orig <- Tn(perm_index = 1:N)
  if (!is.finite(Tn_orig) || Tn_orig == 0) {

    stop("The test statistic is not finite or is zero. Bandwidth too small?")

  }

  # Permutations
  perms <- t(replicate(B, sample(N)))

  # Set seeds for the Monte Carlos inside the JSD and restore at the exit
  if (!is.null(seed_jsd) && type == "jsd") {

    old_seed <- .Random.seed
    on.exit({.Random.seed <<- old_seed})

  }

  # Perform permutations
  pb <- txtProgressBar(style = 3)
  Tn_star <- rep(NA, B)
  for (b in seq_len(B)) {

    # Bootstrapped statistic
    if (!is.null(seed_jsd) && type == "jsd") {

      set.seed(seed_jsd, kind = "Mersenne-Twister")

    }
    Tn_star[b] <- Tn(perm_index = perms[b, ])
    setTxtProgressBar(pb = pb, value = b / B)

    # Current p-value
    pvalue <- mean(Tn_star > Tn_orig, na.rm = TRUE)

    # Plot the position of the original statistic with respect to the
    # permutation replicates? Do it one out of ten replicates
    if (plot_boot && (b %% 10 == 0 || b == B)) {

      hist(Tn_star, probability = TRUE, main = paste("p-value:",
                                                     sprintf("%.4f", pvalue)),
           xlab = expression(T[n]^"*"),
           xlim = range(c(Tn_star, Tn_orig), na.rm = TRUE))
      rug(Tn_star)
      abline(v = Tn_orig, col = 2)

    }

  }

  # Construct an "htest" result
  Tn_orig <- c("Tn" = Tn_orig)
  result <- list(statistic = Tn_orig, p.value = pvalue,
                 statistic_perm = drop(Tn_star), n = n, h = h, B = B,
                 alternative = alternative, method = method,
                 data.name = deparse(substitute(X)))
  class(result) <- "htest"

  # Return "htest"
  return(result)

}


#' @title Hellinger distance between two densities via Monte Carlo
#'
#' @description Computes the Hellinger distance
#' \deqn{H(f, g) = \sqrt(1 - \int_{\mathcal{S}^{d_1} \times \ldots \times
#' \mathcal{S}^{d_k}} \sqrt(f(\boldsymbol{x}) g(\boldsymbol{x}))
#' d\boldsymbol{x})} between two densities \eqn{f} and \eqn{g} on
#' \eqn{\mathcal{S}^{d_1} \times \ldots \times \mathcal{S}^{d_k}} via
#' Monte Carlo.
#'
#' @param log_f,log_g logarithms of \eqn{f} and \eqn{g} evaluated in a Monte
#' Carlo sample.
#' @inheritParams kde_polysph
#' @param bhatta compute the Bhattacharyya divergence (or coefficient) instead?
#' @return A scalar with the estimated distance.
#' @examples
#' # Example with von Mises--Fisher distributions
#' M <- 1e3
#' d <- c(1, 3)
#' mu <- r_unif_polysph(n = 1, d = d)
#' kappa <- c(1, 5)
#' x_mc <- r_unif_polysph(n = M, d = d)
#' log_f <- d_vmf_polysph(x = x_mc, d = d, mu = mu, kappa = kappa, log = TRUE)
#' log_g <- d_vmf_polysph(x = x_mc, d = d, mu = -mu, kappa = kappa, log = TRUE)
#' polykde:::hd_mc(log_f = log_f, log_g = log_f, d = d)
#' polykde:::hd_mc(log_f = log_f, log_g = log_g, d = d)
#' polykde:::hd_mc(log_f = log_f, log_g = log_f, d = d, bhatta = TRUE)
#' polykde:::hd_mc(log_f = log_f, log_g = log_g, d = d, bhatta = TRUE)
#' @keywords internal
hd_mc <- function(log_f, log_g, d, bhatta = FALSE) {

  # Hellinger distance: H(f, g) = \sqrt(1 - \int \sqrt(f(x) * g(x)) dx)
  # Bhattacharyya divergence: DB(f; g) = -\log(\int \sqrt(f(x) * g(x)) dx)
  # H(f, g) = \sqrt(1 - \int \sqrt(f(x) * g(x)) dx)
  #         = \sqrt(1 - \exp(\log(\int \sqrt(f(x) * g(x)) dx)))
  #         = \sqrt(1 - \exp(-DB(f; g)))
  # DB(f; g) = -\log(\int \sqrt(f(x) * g(x)) dx)
  #          = -\log(\int \exp((\log(f(x)) + \log(g(x))) / 2) dx)
  #          ≈ -\log(w / M * \sum_i \exp((\log(f(X_i)) + \log(g(X_i))) / 2))
  #          = \log(M / w) -
  #            \log(\sum_i \exp((\log(f(X_i)) + \log(g(X_i))) / 2))

  # Bhattacharyya divergence
  M <- length(log_f)
  stopifnot(M == length(log_g))
  logs <- 0.5 * (log_f + log_g)
  max_log <- max(logs)
  bhatta_div <- log(M) - sum(rotasym::w_p(p = d + 1, log = TRUE)) -
    (max_log + log(sum(exp(logs - max_log))))

  # Hellinger distance
  hellin_dis <- sqrt(1 - exp(-bhatta_div))
  return(ifelse(bhatta, bhatta_div, hellin_dis))

}
