
#' @title Local polynomial estimator for polyspherical-on-scalar regression
#'
#' @description Computes a local constant or local linear estimator with
#' polyspherical response and scalar predictor.
#'
#' @param x a vector of size \code{nx} with the evaluation points.
#' @param X a vector of size \code{n} with the predictor sample.
#' @param Y a matrix of size \code{c(n, sum(d) + r)} with the response sample
#' on the polysphere.
#' @inheritParams kde_polysph
#' @param h a positive scalar giving the bandwidth.
#' @param p degree of local fit, either \code{0} or \code{1}. Defaults to
#' \code{0}.
#' @return A vector of size \code{nx} with the estimated regression curve
#' evaluated at \code{x}.
#' @examples
#' x_grid <- seq(-0.25, 1.25, l = 200)
#' n <- 50
#' X <- seq(0, 1, l = n)
#' Y <- r_path_s2r(n = n, r = 1, sigma = 0.1, spiral = TRUE)[, , 1]
#' h0 <- bw_cv_kre_polysph(X = X, Y = Y, d = 2, p = 0, plot_cv = FALSE)$h_1se
#' sc3 <- scatterplot3d::scatterplot3d(Y, pch = 16, xlim = c(-1, 1),
#'                                     ylim = c(-1, 1), zlim = c(-1, 1),
#'                                     xlab = "", ylab = "", zlab = "")
#' sc3$points3d(kre_polysph(x = x_grid, X = X, Y = Y, d = 2, h = h0, p = 0),
#'              pch = 16, type = "l", col = 2, lwd = 2)
#' sc3$points3d(kre_polysph(x = x_grid, X = X, Y = Y, d = 2, h = h0, p = 1),
#'              pch = 16, type = "l", col = 3, lwd = 2)
#' @export
kre_polysph <- function(x, X, Y, d, h, p = 0) {

  # Quick checks
  stopifnot(is.matrix(Y))
  stopifnot(length(X) == nrow(Y))
  stopifnot(sum(d + 1) == ncol(Y))
  stopifnot(p == 0 || p == 1)

  # Log-kernel weights
  nx <- length(x)
  mh2 <- -0.5 / h^2
  Kx <- matrix(nrow = nx, ncol = length(X))
  for (i in seq_along(X)) {

    Kx[, i] <- mh2 * (x - X[i])^2

  }

  # Safe computation of the Nadaraya--Watson weights
  W <- sdetorus::safeSoftMax(logs = Kx, expTrc = Inf)

  # Local constant or local linear?
  if (p == 0) {

    # Weighted means
    Y_hat <- W %*% Y

  } else {

    Y_hat <- matrix(nrow = nx, ncol = ncol(Y))
    for (j in seq_len(nx)) {

      Y_hat[j, ] <- lm.wfit(x = cbind(1, x[j] - X), y = Y,
                            w = W[j, ])$coefficients[1, ]

    }

  }

  # Project to polysphere
  Y_hat <- proj_polysph(x = Y_hat, ind_dj = comp_ind_dj(d = d))
  return(Y_hat)

}


#' @title Cross-validation bandwidth selection for polyspherical-on-scalar
#' regression
#'
#' @description Computes least squares cross-validation bandwidths for kernel
#' polyspherical-on-scalar regression. It computes both the bandwidth that
#' minimizes the cross-validation loss and its "one standard error" variant.
#'
#' @inheritParams kre_polysph
#' @param h_grid bandwidth grid where to optimize the cross-validation loss.
#' Defaults to \code{bw.nrd(X) * 10^seq(-1, 1, l = 100)}.
#' @param plot_cv plot the cross-validation loss curve? Defaults to \code{TRUE}.
#' @param fast use the faster and equivalent version of the cross-validation
#' loss? Defaults to \code{TRUE}.
#' @details A similar output to \code{glmnet}'s \code{\link[glmnet]{cv.glmnet}}
#' is returned.
#' @return A list with the following fields:
#' \item{h_min}{the bandwidth that minimizes the cross-validation loss.}
#' \item{h_1se}{the largest bandwidth within one standard error of the
#' minimal cross-validation loss.}
#' \item{cvm}{the mean of the cross-validation loss curve.}
#' \item{cvse}{the standard error of the cross-validation loss curve.}
#' @examples
#' n <- 50
#' X <- seq(0, 1, l = n)
#' Y <- r_path_s2r(n = n, r = 1, sigma = 0.1, spiral = TRUE)[, , 1]
#' bw_cv_kre_polysph(X = X, Y = Y, d = 2, p = 0)
#' bw_cv_kre_polysph(X = X, Y = Y, d = 2, p = 1, fast = FALSE)
#' @export
bw_cv_kre_polysph <- function(X, Y, d, p = 0, h_grid = bw.nrd(X) *
                                10^seq(-2, 2, l = 100), plot_cv = TRUE,
                              fast = TRUE) {

  # Addends of the CV loss
  n <- length(X)
  ind_dj <- comp_ind_dj(d = d)
  if (fast) {

    # Local constant or local linear?
    if (p == 0) {

      cv_kre <- function(h) {

        # Log-kernel ij-matrix log(K_h(X_i - X_j))
        K <- (-0.5 / h^2) * dist(X)^2
        K <- as.matrix(K)

        # Weights ij-matrix W_i(X_j): \sum_i W_i(X_j) = 1 -- Skip and compute
        # directly W_minus_i for a better numerical stability
        # diag(K) <- 1
        # W <- sdetorus::safeSoftMax(logs = K, expTrc = 100)

        # Safe computation of the weights ij-matrix without i-th datum.
        # Normalize by row once ii-entry is removed:
        # W_{-i,j}(X_i) = W_j(X_i) / (1 - W_i(X_i)), j = 1, ..., n.
        diag(K) <- -Inf
        W_minus_i <- sdetorus::safeSoftMax(logs = K, expTrc = 100)

        # Fits without the i-th datum
        Y_hat_i <- W_minus_i %*% Y
        Y_hat_i <- proj_polysph(x = Y_hat_i, ind_dj = ind_dj)

        # Distances
        dist_polysph(x = Y, y = Y_hat_i, ind_dj = ind_dj, norm_x = TRUE,
                     norm_y = TRUE, std = FALSE)^2

      }

    } else {

      stop("Not implemented yet")

    }

  } else {

    cv_kre <- function(h) {

      Y_hat_i <- t(sapply(1:n, function(i) {
        kre_polysph(x = X[i], X = X[-i], Y = Y[-i, ], d = d, h = h, p = p)
      }))
      dist_polysph(x = Y, y = Y_hat_i, ind_dj = ind_dj, norm_x = TRUE,
                   norm_y = TRUE, std = FALSE)^2

    }

  }

  # Mean and standard deviation of the CV loss
  cv <- t(sapply(h_grid, function(h) cv_kre(h = h)))
  cvm <- rowMeans(cv)
  cvse <- apply(cv, 1, sd)

  # CV bandwidth
  ind_min <- which.min(cvm)
  h_min <- h_grid[ind_min]

  # CV-1SE bandwidth
  cvmse_min <- (cvm + cvse)[ind_min]
  ind_1se <- which(cvmse_min > cvm)
  ind_1se <- ind_1se[length(ind_1se)]
  h_1se <- h_grid[ind_1se]

  # Add plot of the CV loss
  if (plot_cv) {

    plot(h_grid, cvm, type = "o", log = "xy",
         ylim = c(min(cvm), max(cvm + cvse)))
    lines(h_grid, cvm + cvse, col = "gray")
    lines(h_grid, pmax(cvm - cvse, .Machine$double.eps), col = "gray")
    abline(v = c(h_min, h_1se), col = 2, lwd = 2, lty = c(1, 3))
    abline(h = cvmse_min, lty = 2, col = "gray")
    rug(h_grid)

  }

  # CV bandwidth
  return(list("h_min" = h_min, "h_1se" = h_1se, "cvm" = cvm, "cvse" = cvse))

}
