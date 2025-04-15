
#' @title Cross-validation bandwidth selection for polyspherical kernel
#' density estimator
#'
#' @description Likelihood Cross-Validation (LCV) and Least Squares
#' Cross-Validation (LSCV) bandwidth selection for the polyspherical kernel
#' density estimator.
#'
#' @inheritParams kde_polysph
#' @param type cross-validation type, either \code{"LCV"} (default) or
#' \code{"LSCV"}.
#' @param M Monte Carlo samples to use for approximating the integral in
#' the LSCV loss.
#' @param bw0 initial bandwidth vector for minimizing the CV loss. If
#' \code{NULL}, it is computed internally by magnifying the
#' \code{\link{bw_rot_polysph}} bandwidths by 50\%. Can be also a matrix of
#' initial bandwidth vectors.
#' @param na.rm remove \code{NA}s in the objective function? Defaults to
#' \code{FALSE}.
#' @param h_min minimum h enforced (componentwise). Defaults to \code{0}.
#' @inheritParams bw_rot_polysph
#' @param upscale rescale the resulting bandwidths to work for derivative
#' estimation? Defaults to \code{FALSE}.
#' @param imp_mc use importance sampling in the Monte Carlo approximation of
#' the integral in the LSCV loss? It is more accurate but also more time
#' consuming. Defaults to \code{TRUE}.
#' @param seed_mc seed for the Monte Carlo simulations used to estimate the
#' integral in the LSCV loss. Defaults to \code{NULL} (no seed is fixed for
#' different bandwidths).
#' @param exact_vmf use the closed-form for the LSCV loss with the
#' von Mises--Fisher kernel? Defaults to \code{FALSE}.
#' @param common_h use the same bandwidth for all dimensions? Defaults to
#' \code{FALSE}.
#' @param spline use a faster spline approximation to compute Bessel functions?
#' Defaults to \code{FALSE}.
#' @param opt optimizer to use; either \code{"\link{optim}"} (default) or
#' \code{"\link{nlm}"}.
#' @param ncores number of cores used during the optimization. Defaults to
#' \code{1}.
#' @param ... further arguments passed to \code{\link{optim}} or
#' \code{\link{nlm}} (if \code{ncores = 1}) or
#' \code{\link[optimParallel]{optimParallel}} (if \code{ncores > 1}).
#' @details If \code{bw0} is a matrix, then the optimization is started at that
#' row of bandwidths that is most promising for the optimization, i.e., the
#' bandwidths that minimized the CV loss.
#' @return A list with entries \code{bw} (optimal bandwidth) and \code{opt},
#' the latter containing the output of \code{\link[stats]{nlm}},
#' \code{\link[stats]{optim}}, or \code{\link[optimParallel]{optimParallel}}.
#' @examples
#' n <- 20
#' d <- 1:2
#' kappa <- rep(10, 2)
#' X <- r_vmf_polysph(n = n, d = d, mu = r_unif_polysph(n = 1, d = d),
#'                    kappa = kappa)
#' bw_cv_polysph(X = X, d = d, type = "LCV")$bw
#' bw_cv_polysph(X = X, d = d, type = "LSCV", exact_vmf = TRUE)$bw
#' @export
bw_cv_polysph <- function(X, d, kernel = 1, kernel_type = 1, k = 10,
                          intrinsic = FALSE, type = c("LCV", "LSCV")[1],
                          M = 1e4, bw0 = NULL, na.rm = FALSE, h_min = 0,
                          upscale = FALSE, deriv = 0, imp_mc = TRUE,
                          seed_mc = NULL, exact_vmf = FALSE, common_h = FALSE,
                          spline = FALSE, opt = c("optim", "nlm")[1],
                          ncores = 1, ...) {

  # Check dimensions
  if (ncol(X) != sum(d + 1)) {

    stop("X and d are incompatible.")

  }
  n <- nrow(X)
  r <- length(d)

  # Objective function
  if (type == "LCV") {

    # LCV loss
    obj <- function(log_h) {

      # Replicate bandwidths if common_h
      if (common_h) {

        log_h <- rep(log_h, r)

      }

      # Ensure positivity and lower bound (if not enforced by optim())
      h_pos <- exp(log_h)
      penalty <- 0
      if (any(h_pos < h_min)) {

        message("h left-truncated to h_min")
        penalty <- 1e6 * sum((h_pos - h_min)^2)
        h_pos <- pmax(h_pos, h_min)

      }

      # Wrap in a tryCatch() to avoid too small bandwidth errors
      loss <- tryCatch({

        # Log-cv kernels
        log_cv <- log_cv_kde_polysph(X = X, d = d, h = h_pos, wrt_unif = TRUE,
                                     kernel = kernel, kernel_type = kernel_type,
                                     k = k, intrinsic = intrinsic)

        # -LCV
        loss <- -sum(log_cv, na.rm = na.rm) + penalty
        ifelse(!is.finite(loss), 1e6, loss)

      }, error = function(e) 1e6)
      return(loss)

    }

  } else if (type == "LSCV") {

    if (kernel == "1" && exact_vmf) {

      if (intrinsic) {

        stop("Intrinsic LSCV loss is not available with exact_vmf = TRUE.")

      }

      # Precompute matrix with the lower triangular parts of the matrices
      # (X_{il}'X_{jl})_{ij}, l = 1, ..., r.
      ind_dj <- comp_ind_dj(d = d)
      Xi_Xj_l <- diamond_rcrossprod(X = X, ind_dj = ind_dj)
      Xi_Xj_l <- sapply(seq_len(r), function(j) {

        Xi_Xj_l[, , j][lower.tri(Xi_Xj_l[, , j], diag = FALSE)]

      })
      Xi_Xj_l <- t(Xi_Xj_l) # For better column recycling later
      norm_Xi_Xj_l <- sqrt(2 * (1 + Xi_Xj_l))

      # Precompute other fixed objects in the LSCV loss
      log_n <- log(n)
      log_2_n1 <- log(2 / (n - 1))
      n2 <- n * (n - 1) / 2

      # LSCV loss
      obj <- function(log_h) {

        # Replicate bandwidths if common_h
        if (common_h) {

          log_h <- rep(log_h, r)

        }

        # Ensure positivity and lower bound (if not enforced by optim())
        h_pos <- exp(log_h)
        penalty <- 0
        if (any(h_pos < h_min)) {

          message("h left-truncated to h_min")
          penalty <- 1e6 * sum((h_pos - h_min)^2)
          h_pos <- pmax(h_pos, h_min)

        }

        # Wrap in a tryCatch() to avoid too small bandwidth errors
        exact_loss <- tryCatch({

          # Log-constants
          h_pos2 <- 1 / h_pos^2
          log_c_h2 <- sum(fast_log_c_vMF(p = d + 1, kappa = h_pos2,
                                         spline = FALSE))
          log_c_2h2 <- sum(fast_log_c_vMF(p = d + 1, kappa = 2 * h_pos2,
                                          spline = FALSE))

          # Compute X_{il}'X_{jl} / h_l^2 and
          # \sum_l \log(c_vMF(||X_{il}'X_{jl}|| / h_l^2))
          Xi_Xj_l_h <- .colSums(Xi_Xj_l * h_pos2, m = r, n = n2)
          log_c_norm_Xi_Xj_l_h <- numeric(n2)
          for (l in seq_len(r)) {

            log_c_norm_Xi_Xj_l_h <- log_c_norm_Xi_Xj_l_h +
              fast_log_c_vMF(p = d[l] + 1, kappa = norm_Xi_Xj_l[l, ] *
                               h_pos2[l], spline = spline)

          }

          # CV loss
          cv_1 <- exp(2 * log_c_h2 - log_c_2h2 - log_n)
          cv_2 <- 2 * sum(exp(Xi_Xj_l_h + log_c_h2 + (log_2_n1 - log_n)) -
                            exp(2 * (log_c_h2 - log_n) - log_c_norm_Xi_Xj_l_h))
          cv <- cv_1 - cv_2
          loss <- cv + penalty
          ifelse(!is.finite(loss), 1e6, loss)

        }, error = function(e) 1e6)
        return(exact_loss)

      }

    } else {

      # Set seeds for the Monte Carlo
      if (!is.null(seed_mc)) {

        # old_seed <- .Random.seed
        # on.exit({.Random.seed <<- old_seed})
        set.seed(seed_mc, kind = "Mersenne-Twister")

      }

      # Common uniform Monte Carlo sample
      if (!imp_mc) {

        mc_samp <- r_unif_polysph(n = M, d = d)

      }

      # LSCV loss
      obj <- function(log_h) {

        # Replicate bandwidths if common_h
        if (common_h) {

          log_h <- rep(log_h, r)

        }

        # Ensure positivity and lower bound (if not enforced by optim())
        h_pos <- exp(log_h)
        penalty <- 0
        if (any(h_pos < h_min)) {

          message("h left-truncated to h_min")
          penalty <- 1e6 * sum((h_pos - h_min)^2)
          h_pos <- pmax(h_pos, h_min)

        }

        # Wrap in a tryCatch() to avoid too small bandwidth errors
        loss <- tryCatch({

          # Integral part with LogSumExp trick

          # Use importance-sampling Monte Carlo?
          if (imp_mc) {

            # Set seeds for the Monte Carlo
            if (!is.null(seed_mc)) {

              set.seed(seed_mc, kind = "Mersenne-Twister")

            }

            # h-dependent Monte Carlo
            mc_kde_samp <- r_kde_polysph(n = M, X = X, d = d,
                                         h = h_pos, kernel = kernel,
                                         kernel_type = kernel_type, k = k,
                                         intrinsic = intrinsic)
            log_kde2_mc <- kde_polysph(x = mc_kde_samp, X = X, d = d,
                                       h = h_pos, wrt_unif = FALSE,
                                       kernel = kernel, kernel_type =
                                         kernel_type, k = k, log = TRUE,
                                       intrinsic = intrinsic) - log(M)

          } else {

            # Uniform Monte Carlo
            log_kde2_mc <- 2 * kde_polysph(x = mc_samp, X = X, d = d,
                                           h = h_pos, wrt_unif = FALSE,
                                           kernel = kernel, kernel_type =
                                             kernel_type, k = k,
                                           log = TRUE, intrinsic = intrinsic) -
              log(M) + sum(rotasym::w_p(p = d + 1, log = TRUE))

          }
          max_log_kde2_mc <- max(log_kde2_mc)
          log_int_kde2 <- ifelse(is.finite(max_log_kde2_mc),
                                 max_log_kde2_mc +
                                   log(sum(exp(log_kde2_mc - max_log_kde2_mc))),
                                 -Inf)

          # Sum part with LogSumExp trick
          log_cv_kde <- log_cv_kde_polysph(X = X, d = d, h = h_pos,
                                           wrt_unif = FALSE, kernel = kernel,
                                           kernel_type = kernel_type, k = k,
                                           intrinsic = intrinsic) -
            log(n) + log(2)
          max_log_cv_kde <- max(log_cv_kde)
          log_sum_cv_kde <- ifelse(is.finite(max_log_cv_kde),
                                   max_log_cv_kde +
                                     log(sum(exp(log_cv_kde - max_log_cv_kde))),
                                   -Inf)

          # CV loss
          cv <- exp(log_int_kde2) - exp(log_sum_cv_kde)
          loss <- cv + penalty
          ifelse(!is.finite(loss), 1e6, loss)

        }, error = function(e) 1e6)
        return(loss)

      }

    }

  } else {

    stop("\"type\" must be \"LCV\" or \"LSCV\".")

  }

  # Set initial bandwidths
  if (is.null(bw0)) {

    # 50% larger ROT bandwidths
    bw0 <- 1.5 * bw_rot_polysph(X = X, d = d, kernel = kernel,
                                kernel_type = kernel_type, k = k,
                                upscale = FALSE, deriv = 0)$bw
    bw0 <- rbind(bw0)

  } else {

    bw0 <- rbind(bw0)
    if (r != ncol(bw0)) {

      stop("bw0 and d are incompatible.")

    }

  }

  # Find the best starting bandwidths to run the optimization on it
  if (nrow(bw0) > 1) {

    if (common_h) {

      bw0 <- cbind(rowMeans(bw0))

    }
    ind_best <- which.min(apply(log(bw0), 1, obj))
    bw0 <- bw0[ind_best, ]

  } else {

    # Replicate bandwidths if common_h
    if (common_h) {

      bw0 <- rowMeans(bw0)

    } else {

      bw0 <- drop(bw0)

    }

  }
  if ((!is.null(list(...)$control$trace) && list(...)$control$trace > 0) ||
      (!is.null(list(...)$print.level) && list(...)$print.level > 0)) {

    message("bw0 = ", bw0)

  }

  # Optimization
  if (opt == "nlm") {

    if (ncores == 1) {

      opt <- nlm(f = obj, p = log(bw0), ...)
      bw <- exp(opt$estimate)

    } else {

      stop("nlm() with ncores > 1 not available.")

    }

  } else {

    if (ncores == 1) {

      opt <- optim(par = log(bw0), fn = obj, ...)
      bw <- exp(opt$par)

    } else {

      cl <- parallel::makeCluster(spec = ncores)
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterExport(cl = cl, varlist = ls(), envir = environment())
      parallel::clusterCall(cl, function() {
        library("polykde")
      })
      opt <- optimParallel::optimParallel(par = log(bw0), fn = obj,
                                          parallel = list(cl = cl,
                                                          forward = FALSE,
                                                          loginfo = TRUE),
                                          ...)
      parallel::stopCluster(cl)
      bw <- exp(opt$par)

    }

  }

  # Upscale?
  if (upscale > 0) {

    n_up <- n^(1 / (d * r + 4)) * n^(-1 / (d * r + 2 * deriv + 4))
    bw <- bw * n_up

  }
  return(list("bw" = unname(bw), "opt" = opt))

}


#' @title Minimum bandwidth allowed in likelihood cross-validation for
#' Epanechnikov kernels
#'
#' @description This function computes the minimum bandwidth allowed in
#' likelihood cross-validation with Epanechnikov kernels, for a given dataset
#' and dimension.
#'
#' @inheritParams bw_cv_polysph
#' @return The minimum bandwidth allowed.
#' @examples
#' n <- 5
#' d <- 1:3
#' X <- r_unif_polysph(n = n, d = d)
#' h_min <- rep(bw_lcv_min_epa(X = X, d = d), length(d))
#' log_cv_kde_polysph(X = X, d = d, h = h_min - 1e-4, kernel = 2) # Problem
#' log_cv_kde_polysph(X = X, d = d, h = h_min + 1e-4, kernel = 2) # OK
#' @export
bw_lcv_min_epa <- function(X, d, kernel_type = c("prod", "sph")[1]) {

  # Make kernel_type character
  if (is.numeric(kernel_type)) {

    kernel_type <- switch(kernel_type, "1" = "prod", "2" = "sph",
                          stop("kernel_type must be 1 or 2."))

  }

  # Check dimensions
  if (ncol(X) != sum(d + 1)) {

    stop("X and d are incompatible.")

  }
  n <- nrow(X)
  r <- length(d)

  # Index for accessing each S^dj with ind[j]:(ind[j + 1] - 1)
  ind <- cumsum(c(1, d + 1))

  if (kernel_type == "prod") {

    # Compute max_k X_i,k' X_j,k
    prods <- matrix(0, nrow = n, ncol = n)
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        prods[i, j] <- max(sapply(1:r, function(k) {
          ind_k <- ind[k]:(ind[k + 1] - 1)
          1 - sum(X[i, ind_k] * X[j, ind_k])
        }))
      }
    }

    # Symmetrize
    prods <- prods + t(prods)
    diag(prods) <- Inf

    # Apply max_i min_{j\neq i}
    return(sqrt(max(apply(prods, 1, min))))

  } else if (kernel_type == "sph") {

    stop("Not implemented yet.")

  } else {

    stop("kernel_type must be either \"prod\" or \"sph\".")

  }

}


#' @title Curvature of a polyspherical von Mises--Fisher density
#'
#' @description Computes the curvature matrix
#' \eqn{\boldsymbol{R}(\boldsymbol{\kappa})} of a product of von Mises--Fisher
#' densities on the polysphere. This curvature is used in the rule-of-thumb
#' selector \code{\link{bw_rot_polysph}}.
#'
#' @inheritParams r_vmf_polysph
#' @param log compute the (entrywise) logarithm of the curvature matrix?
#' Defaults to \code{FALSE}.
#' @return A matrix of size \code{c(length(r), length(r))}.
#' @examples
#' # Curvature matrix
#' d <- 2:4
#' kappa <- 1:3
#' curv_vmf_polysph(kappa = kappa, d = d)
#' curv_vmf_polysph(kappa = kappa, d = d, log = TRUE)
#'
#' # Equivalence on the sphere with DirStats::R_Psi_mixvmf
#' drop(curv_vmf_polysph(kappa = kappa[1], d = d[1]))
#' d[1]^2 * DirStats::R_Psi_mixvmf(q = d[1], mu = rbind(c(rep(0, d[1]), 1)),
#'                                 kappa = kappa[1], p = 1)
#' @export
curv_vmf_polysph <- function(kappa, d, log = FALSE) {

  # Dimension check
  stopifnot(length(d) %in% c(1, length(kappa)))

  # Log-Bessels
  log_I_1 <- log(besselI(x = 2 * kappa, nu = (d + 1) / 2, expon.scaled = TRUE))
  log_I_2 <- log(besselI(x = 2 * kappa, nu = (d - 1) / 2, expon.scaled = TRUE))
  log_I_3 <- log(besselI(x = kappa, nu = (d - 1) / 2, expon.scaled = TRUE))
  I <- exp(log_I_1 - log_I_2)

  # Auxiliary functions
  v_kappa <- d * kappa * (2 * (2 + d) * kappa - (d * d - d + 2) * I)
  u_kappa <- d * kappa * I
  R0_kappa <- sum(((d - 1) / 2) * log(kappa) + log_I_2 -
                    d * log(2) - ((d + 1) / 2) * log(pi) - 2 * log_I_3)

  # log(R(kappa))
  R_kappa <- tcrossprod(u_kappa)
  diag(R_kappa) <- 0.5 * v_kappa
  R_kappa <- log(R_kappa) + log(0.25) + R0_kappa

  # Logarithm or not?
  if (!log) {

    R_kappa <- exp(R_kappa)

  }
  return(R_kappa)

}


#' @title Rule-of-thumb bandwidth selection for polyspherical kernel
#' density estimator
#'
#' @description Computes the rule-of-thumb bandwidth for the polyspherical
#' kernel density estimator using a product of von Mises--Fisher distributions
#' as reference in the Asymptotic Mean Integrated Squared Error (AMISE).
#'
#' @inheritParams kde_polysph
#' @param bw0 initial bandwidth for minimizing the CV loss. If \code{NULL}, it
#' is computed internally by magnifying the \code{\link{bw_mrot_polysph}}
#' bandwidths by 50\%. Can be also a matrix of initial bandwidth vectors.
#' @param upscale rescale bandwidths to work on
#' \eqn{\mathcal{S}^{d_1}\times\cdots\times \mathcal{S}^{d_r}} and for
#' derivative estimation?
#' Defaults to \code{FALSE}. If \code{upscale = 1}, the order \code{n} is
#' upscaled. If \code{upscale = 2}, then also the kernel constant is upscaled.
#' @param deriv derivative order to perform the upscaling. Defaults to \code{0}.
#' @param kappa estimate of the concentration parameters. Computed if not
#' provided (default).
#' @param ... further arguments passed to \code{\link[stats]{nlm}}.
#' @details The selector assumes that the density curvature matrix
#' \eqn{\boldsymbol{R}} of the unknown density is approximable by that of a
#' product of von Mises--Fisher densities,
#' \eqn{\boldsymbol{R}(\boldsymbol{\kappa})}. The estimation of the
#' concentration parameters \eqn{\boldsymbol{\kappa}} is done by maximum
#' likelihood.
#'
#' If \code{bw0} is a matrix, then the optimization is started at that row of
#' bandwidths that is most promising for the optimization, i.e., the bandwidths
#' that minimized the CV loss.
#' @return A list with entries \code{bw} (optimal bandwidth) and \code{opt},
#' the latter containing the output of \code{\link[stats]{nlm}}.
#' @examples
#' n <- 100
#' d <- 1:2
#' kappa <- rep(10, 2)
#' X <- r_vmf_polysph(n = n, d = d, mu = r_unif_polysph(n = 1, d = d),
#'                    kappa = kappa)
#' bw_rot_polysph(X = X, d = d)$bw
#' @export
bw_rot_polysph <- function(X, d, kernel = 1, kernel_type = c("prod", "sph")[1],
                           bw0 = NULL, upscale = FALSE, deriv = 0, k = 10,
                           kappa = NULL, ...) {

  # Make kernel_type character
  if (is.numeric(kernel_type)) {

    kernel_type <- switch(kernel_type, "1" = "prod", "2" = "sph",
                          stop("kernel_type must be 1 or 2."))

  }

  # Check dimensions
  if (ncol(X) != sum(d + 1)) {

    stop("X and d are incompatible.")

  }
  n <- nrow(X)
  r <- length(d)

  # Estimate kappa
  if (is.null(kappa)) {

    ind <- cumsum(c(1, d + 1))
    kappa <- sapply(seq_along(d), function(j) {

      # Prepare data + fit vMF
      data <- X[, ind[j]:(ind[j + 1] - 1)]
      min(DirStats::norm2(movMF::movMF(x = data, k = 1, type = "S",
                                       maxit = 300)$theta),
          5e4) # Breaking point for later Bessels

    })

  } else {

    stopifnot(length(kappa) == r)

  }

  # We use that
  # (h^2)' R h^2 = tr[(h^2)' R h^2]
  #              = tr[(h^2 (h^2)') R] # tr(AB') = 1' (A o B) 1 = sum(A o B)
  #              = sum[(h^2 (h^2)') o R]
  # And then we compute the logarithms inside the sum:
  # log[(h^2 (h^2)') o R] = log[h^2 (h^2)'] + log(R)
  # This gives a matrix of logarithms, which is ready to log-sum-exp it.
  # Therefore,
  # bias2 = (b o h^2)' R (h^2 o b)
  #       = (h^2)' [R o (bb')] h^2
  # log(R o (bb')) = log(R) + log(bb')

  # Common objects
  log_R_kappa <- curv_vmf_polysph(kappa = kappa, d = d, log = TRUE)
  b <- b_d(kernel = kernel, d = d, k = k, kernel_type = kernel_type)
  v <- v_d(kernel = kernel, d = d, k = k, kernel_type = kernel_type)
  log_bias2 <- log_R_kappa +
    ifelse(kernel_type == "prod", log(tcrossprod(b)), 2 * log(b[1]))
  log_var <- sum(log(v)) - log(n)

  # AMISE and gradient functions
  f_log_amise_stable_log_h <- function(log_h) {

    h <- exp(log_h)
    log_var2 <- log_var - sum(d * log_h)
    logs <- log(tcrossprod(h^2)) + log_bias2 - log_var2
    log_obj <- log_var2 + log_sum_exp(logs = c(logs, 0))
    attr(log_obj, "gradient") <-
      (exp(log(4) + log_bias2 - log_obj + log_h) %*% h^2 -
         exp(log(d) + log_var2 - log_obj - log_h)) * h
    return(log_obj)

  }
  # fn_log_amise_stable <- function(log_h) {
  #
  #   log_var2 <- log_var - sum(d * log_h)
  #   logs <- log(tcrossprod(exp(2 * log_h))) + log_bias2 - log_var2
  #   log_obj <- log_var2 + log_sum_exp(logs = c(logs, 0))
  #   return(log_obj)
  #
  # }
  # gr_log_amise_stable <- function(log_h) {
  #
  #   h <- exp(log_h)
  #   log_var2 <- log_var - sum(d * log_h)
  #   logs <- log(tcrossprod(h^2)) + log_bias2 - log_var2
  #   log_obj <- log_var2 + log_sum_exp(logs = c(logs, 0))
  #   gr <- (exp(log(4) + log_bias2 - log_obj + log_h) %*% h^2 -
  #     exp(log(d) + log_var2 - log_obj - log_h)) * h
  #   return(gr)
  #
  # }

  # Set initial bandwidths
  if (is.null(bw0)) {

    # 50% larger ROT bandwidths
    bw0 <- 1.5 * bw_mrot_polysph(X = X, d = d, kernel = kernel,
                                 kappa = kappa, upscale = FALSE, deriv = 0)
    bw0 <- rbind(bw0)

  } else {

    bw0 <- rbind(bw0)
    if (r != ncol(bw0)) {

      stop("bw0 and d are incompatible.")

    }

  }

  # Find the best starting bandwidths to run the optimization on it
  if (nrow(bw0) > 1) {

    ind_best <- which.min(apply(log(bw0), 1, f_log_amise_stable_log_h))
    bw0 <- bw0[ind_best, ]

  }
  if (!is.null(list(...)$print.level) && list(...)$print.level > 0) {

    message("bw0 = ", bw0)

  }

  # Optimization
  opt <- nlm(p = log(bw0), f = f_log_amise_stable_log_h, ...)
  # opt <- optim(par = log(bw0), fn = fn_log_amise_stable,
  #              gr = gr_log_amise_stable, method = "L-BFGS-B", ...)
  bw <- exp(opt$estimate)

  # Upscale?
  if (upscale > 0) {

    n_up <- n^(1 / (d + 4)) * n^(-1 / (d * r + 2 * deriv + 4))
    if (upscale > 1) {

      stop("Not supported!")

    }
    bw <- bw * n_up

  }
  return(list("bw" = unname(bw), "opt" = opt))

}


#' @title Marginal rule-of-thumb bandwidth selection for polyspherical kernel
#' density estimator
#'
#' @description Computes marginal (sphere by sphere) rule-of-thumb bandwidths
#' for the polyspherical kernel density estimator using a von Mises--Fisher
#' distribution as reference.
#'
#' @inheritParams kde_polysph
#' @inheritParams bw_rot_polysph
#' @return A vector of size \code{r} with the marginal optimal bandwidths.
#' @examples
#' n <- 100
#' d <- 1:2
#' kappa <- rep(10, 2)
#' X <- r_vmf_polysph(n = n, d = d, mu = r_unif_polysph(n = 1, d = d),
#'                    kappa = kappa)
#' bw_rot_polysph(X = X, d = d)$bw
#' bw_mrot_polysph(X = X, d = d)
#' @export
bw_mrot_polysph <- function(X, d, kernel = 1, k = 10, upscale = FALSE,
                            deriv = 0, kappa = NULL) {

  # Check dimensions
  if (ncol(X) != sum(d + 1)) {

    stop("X and d are incompatible.")

  }
  n <- nrow(X)
  r <- length(d)

  # Index for accessing each S^dj with ind[j]:(ind[j + 1] - 1)
  ind <- cumsum(c(1, d + 1))

  # Marginal ROTs bandwidths
  bw <- sapply(seq_along(d), function(j) {

    # Prepare data + fit vMF
    data <- X[, ind[j]:(ind[j + 1] - 1)]
    if (is.null(kappa)) {

      kappa_j <- min(DirStats::norm2(movMF::movMF(x = data, k = 1,
                                                  type = "S",
                                                  maxit = 300)$theta),
                     5e4) # Breaking point for later Bessels

    } else {

      kappa_j <- kappa[j]

    }
    q <- ncol(data) - 1

    # Kernel-specific constants
    if (kernel == 1) {

      log_kernel_num <- log(4 * sqrt(pi))
      log_kernel_den <- 0 # Already canceled with constants in numerator

    } else if (kernel == 2 || kernel == 3) {

      log_kernel_num <- log(v_d(kernel = kernel, d = q, k = k)) +
        (q + 2) * log(2) + ((q + 1) / 2) * log(pi)
      log_kernel_den <- log(4) + 2 * log(b_d(kernel = kernel, d = q, k = k))

    } else {

      stop("kernel must be 1, 2, or 3.")

    }

    # Safe computation of ROT (Proposition 2 in
    # https://arxiv.org/pdf/1306.0517.pdf)
    log_num <- log_kernel_num +
      2 * log(besselI(nu = (q - 1) / 2, x = kappa_j, expon.scaled = TRUE))
    log_den <- log_kernel_den +
      (q + 1) / 2 * log(kappa_j) +
      log(2 * q * besselI(nu = (q + 1) / 2, x = 2 * kappa_j,
                          expon.scaled = TRUE) +
            (2 + q) * kappa_j * besselI(nu = (q + 3) / 2, x = 2 * kappa_j,
                                        expon.scaled = TRUE)) + log(n)
    log_bwd <- (1 / (q + 4)) * (log_num - log_den)
    bwd <- ifelse(is.finite(log_bwd), exp(log_bwd), 0.025)
    return(bwd)

  })

  # Upscale?
  if (upscale > 0) {

    n_up <- n^(1 / (d + 4)) * n^(-1 / (d * r + 2 * deriv + 4))
    if (upscale > 1) {

      cte_d <- (
        (d * v_d(d = d, kernel = kernel, k = k)) /
          (4 * b_d(d = d, kernel = kernel, k = k)^2)
        )^(1 / (d + 4))
      cte_dr <- (
        (d * r * prod(v_d(d = d, kernel = kernel, k = k))) /
          (4 * sum(b_d(d = d, kernel = kernel, k = k))^2)
        )^(1 / (d * r + 2 * deriv + 4))
      cte_up <- cte_dr / cte_d

    } else {

      cte_up <- 1

    }
    bw <- bw * n_up * cte_up

  }
  return(bw)

}
