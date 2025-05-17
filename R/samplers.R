
#' @title Sample uniform polyspherical data
#'
#' @description Simulates from a uniform distribution on the polysphere
#' \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}}.
#'
#' @param n sample size.
#' @inheritParams kde_polysph
#' @return A matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @examples
#' # Simulate uniform data on (S^1)^2
#' r_unif_polysph(n = 10, d = c(1, 1))
#' @export
r_unif_polysph <- function(n, d) {

  do.call(cbind, as.list(sapply(d, function(dj)
    rbind(sphunif::r_unif_sph(n = n, p = dj + 1, M = 1)[, , 1]),
    simplify = FALSE)))

}


#' @title Sample von Mises--Fisher distributed polyspherical data
#'
#' @description Simulates from a product of von Mises--Fisher distributions on
#' the polysphere \eqn{\mathcal{S}^{d_1} \times \cdots \times
#' \mathcal{S}^{d_r}}.
#'
#' @param mu a vector of size \code{sum(d) + r} with the concatenated von
#' Mises--Fisher means.
#' @param kappa a vector of size \code{r} with the von Mises--Fisher
#' concentrations.
#' @inheritParams r_kern_polysph
#' @inheritParams r_unif_polysph
#' @inheritParams kde_polysph
#' @return A matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @examples
#' # Simulate vMF data on (S^1)^2
#' r_vmf_polysph(n = 10, d = c(1, 1), mu = c(1, 0, 0, 1), kappa = c(1, 1))
#' @export
r_vmf_polysph <- function(n, d, mu, kappa, norm_mu = FALSE) {

  # Check dimensions
  r <- length(kappa)
  if (r != length(d)) {

    stop("kappa and d are incompatible.")

  }
  if (length(mu) != sum(d + 1)) {

    stop("mu and d are incompatible.")

  }

  # Index for accessing each S^dj with ind[j]:(ind[j + 1] - 1)
  ind <- cumsum(c(1, d + 1))

  # Normalize mu?
  if (norm_mu) {

    for (j in seq_len(r)) {

      ind_j <- ind[j]:(ind[j + 1] - 1)
      mu[ind_j] <- mu[ind_j] / sqrt(sum(mu[ind_j]^2))

    }

  }

  # Sample
  do.call(cbind, as.list(sapply(1:r, function(j) {
    ind_j <- ind[j]:(ind[j + 1] - 1)
    rotasym::r_vMF(n = n, mu = mu[ind_j], kappa = kappa[j])
  }, simplify = FALSE)))

}


#' @title Sample kernel-distributed polyspherical data
#'
#' @description Simulates from the distribution defined by a kernel on the
#' polysphere \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}}.
#'
#' @param mu a vector of size \code{sum(d) + r} with the concatenated means
#' that define the center of the kernel.
#' @param norm_mu ensure a normalization of \code{mu}? Defaults to \code{FALSE}.
#' @inheritParams r_unif_polysph
#' @inheritParams kde_polysph
#' @details Simulation for non-von Mises--Fisher spherically symmetric kernels
#' is done by acceptance-rejection from a von Mises--Fisher proposal
#' distribution.
#' @return A matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @examples
#' # Simulate kernels in (S^1)^2
#' n <- 1e3
#' h <- c(1, 1)
#' d <- c(1, 1)
#' mu <- rep(DirStats::to_cir(pi), 2)
#' samp_ker <- function(kernel, kernel_type, col, main) {
#'   data <- r_kern_polysph(n = n, d = d, mu = mu, h = h, kernel = kernel,
#'                          kernel_type = kernel_type)
#'   ang <- cbind(DirStats::to_rad(data[, 1:2]),
#'                DirStats::to_rad(data[, 3:4]))
#'   plot(ang, xlim = c(0, 2 * pi), ylim = c(0, 2 * pi), pch = 16, cex = 0.25,
#'        col = col, xlab = expression(Theta[1]), ylab = expression(Theta[2]),
#'        main = main)
#' }
#' old_par <- par(mfcol = c(2, 3))
#' samp_ker(kernel = 2, kernel_type = 2, col = 1, main = "Epa sph. symmetric")
#' samp_ker(kernel = 2, kernel_type = 1, col = 2, main = "Epa product")
#' samp_ker(kernel = 3, kernel_type = 2, col = 1, main = "Sfp sph. symmetric")
#' samp_ker(kernel = 3, kernel_type = 1, col = 2, main = "Sfp product")
#' samp_ker(kernel = 1, kernel_type = 2, col = 1, main = "vMF sph. symmetric")
#' samp_ker(kernel = 1, kernel_type = 1, col = 2, main = "vMF product")
#' par(old_par)
#' \donttest{
#' # Simulate kernels in (S^2)^2
#' n <- 1e3
#' h <- c(0.2, 0.6)
#' d <- c(2, 2)
#' mu <- c(c(0, 0, 1), c(0, -1, 0))
#' samp_ker <- function(kernel, kernel_type, main) {
#'   data <- r_kern_polysph(n = n, d = d, mu = mu, h = h, kernel = kernel,
#'                          kernel_type = kernel_type)
#'   scatterplot3d::scatterplot3d(rbind(data[, 1:3], data[, 4:6]),
#'                                xlim = c(-1, 1), ylim = c(-1, 1),
#'                                zlim = c(-1, 1), pch = 16, xlab = "",
#'                                ylab = "", zlab = "", cex.symbols = 0.5,
#'        color = rep(viridis::viridis(n)[rank(data[, 3])], 2), main = main)
#' }
#' old_par <- par(mfcol = c(2, 3))
#' samp_ker(kernel = 2, kernel_type = 2, main = "Epa sph. symmetric")
#' samp_ker(kernel = 2, kernel_type = 1, main = "Epa product")
#' samp_ker(kernel = 3, kernel_type = 2, main = "Sfp sph. symmetric")
#' samp_ker(kernel = 3, kernel_type = 1, main = "Sfp product")
#' samp_ker(kernel = 1, kernel_type = 2, main = "vMF sph. symmetric")
#' samp_ker(kernel = 1, kernel_type = 1, main = "vMF product")
#' par(old_par)
#'
#' # Plot simulated data
#' n <- 1e3
#' h <- c(1, 1)
#' d <- c(2, 2)
#' samp_ker <- function(kernel, kernel_type, col, main) {
#'   X <- r_kern_polysph(n = n, d = d, mu = mu, h = h, kernel = kernel,
#'                       kernel_type = kernel_type)
#'   S <- cbind((1 - X[, 1:3] %*% mu[1:3]) / h[1]^2,
#'              (1 - X[, 4:6] %*% mu[4:6]) / h[2]^2)
#'   plot(S, xlim = c(0, 2 / h[1]^2), ylim = c(0, 2 / h[2]^2), pch = 16,
#'        cex = 0.25, col = col, xlab = expression(t[1]),
#'        ylab = expression(t[2]), main = main)
#'   t_grid <- seq(0, 2 / min(h)^2, l = 100)
#'   gr <- as.matrix(expand.grid(t_grid, t_grid))
#'   if (kernel_type == "1") {
#'
#'     dens <- prod(c_kern(h = h, d = d, kernel = kernel, kernel_type = 1)) *
#'       L(gr[, 1], kernel = kernel) * L(gr[, 2], kernel = kernel)
#'
#'   } else if (kernel_type == "2") {
#'
#'     dens <- c_kern(h = h, d = d, kernel = kernel, kernel_type = 2) *
#'       L(gr[, 1] + gr[, 2], kernel = kernel)
#'
#'   }
#'   dens <- matrix(dens, nrow = length(t_grid), ncol = length(t_grid))
#'   contour(t_grid, t_grid, dens, add = TRUE, col = col,
#'           levels = seq(0, 0.2, l = 41))
#'
#' }
#' old_par <- par(mfcol = c(2, 3))
#' samp_ker(kernel = 2, kernel_type = 2, col = 1, main = "Epa sph. symmetric")
#' samp_ker(kernel = 2, kernel_type = 1, col = 2, main = "Epa product")
#' samp_ker(kernel = 3, kernel_type = 2, col = 1, main = "Sfp sph. symmetric")
#' samp_ker(kernel = 3, kernel_type = 1, col = 2, main = "Sfp product")
#' samp_ker(kernel = 1, kernel_type = 2, col = 1, main = "vMF sph. symmetric")
#' samp_ker(kernel = 1, kernel_type = 1, col = 2, main = "vMF product")
#' par(old_par)
#' }
#' @export
r_kern_polysph <- function(n, d, mu, h, kernel = 1, kernel_type = 1, k = 10,
                           intrinsic = FALSE, norm_mu = FALSE) {

  # Check dimensions
  r <- length(h)
  if (r != length(d)) {

    stop("h and d are incompatible.")

  }
  if (length(mu) != sum(d + 1)) {

    stop("mu and d are incompatible.")

  }

  # Index for accessing each S^dj with ind[j]:(ind[j + 1] - 1)
  ind <- cumsum(c(1, d + 1))

  # Normalize mu?
  if (norm_mu) {

    for (j in seq_len(r)) {

      ind_j <- ind[j]:(ind[j + 1] - 1)
      mu[ind_j] <- mu[ind_j] / sqrt(sum(mu[ind_j]^2))

    }

  }

  # Stop if intrinsic
  if (intrinsic) {

   stop("Sampling for intrisic kernels is not implemented yet.")

  }

  # General algorithm:
  # 1. Sample V ~ [-1, 1]^r according to the kernel
  # 2. For each j = 1, ..., r, set
  #  X_j = V_j * \mu_j + sqrt(1 - V_j^2) * \Gamma_{\mu_j} * U_j,
  #  where U_j ~ Unif(S^{d_j}).

  # Product kernel or vMF
  if (kernel_type == "1" || kernel == "1") {

    # Exploit rotasym::r_tang_norm()
    samp <- do.call(cbind, as.list(sapply(seq_len(r), function(j) {
      ind_j <- ind[j]:(ind[j + 1] - 1)
      rotasym::r_tang_norm(n = n, theta = mu[ind_j],
                           r_U = function(n) rotasym::r_unif_sphere(n = n,
                                                                    p = d[j]),
                           r_V = function(n) r_g_kern(n = n, d = d[j], h = h[j],
                                                      kernel = kernel, k = k))
    }, simplify = FALSE)))

  # Spherically symmetric kernels
  } else if (kernel_type == "2") {

    # Acceptance-rejection from vMF kernel

    # Get M
    t_grid <- seq(-1, 1, l = 100)
    t_grid <- rowSums((1 - t_grid) %o% h^(-2))
    const_f <- c_kern(h = h, d = d, kernel = kernel, kernel_type = "2", k = k)
    const_g <- c_kern(h = h, d = d, kernel = "1", kernel_type = "2")
    dens_f <- L(t = t_grid, kernel = kernel, k = k)
    dens_g <- L(t = t_grid, kernel = "1")
    dens_ratio <- dens_g / dens_f
    M <- 1 / min(dens_ratio[is.finite(dens_ratio)])
    ratio <- const_f / (const_g * M)

    # Acceptance-rejection
    V <- matrix(NA, nrow = n, ncol = r)
    i <- 1
    while (i <= n) {

      # vMF sample
      y <- sapply(seq_len(r), function(j)
        r_g_kern(n = 1, d = d[j], h = h[j], kernel = "1"))
      s <- sum((1 - y) / h^2)

      # Rejection step
      u <- runif(1)
      if (u < (ratio * L(t = s, kernel = kernel, k = k) /
               L(t = s, kernel = "1"))) {

        V[i, ] <- y
        i <- i + 1

      }

    }

    # Simulate data on the polysphere
    samp <- do.call(cbind, as.list(sapply(seq_len(r), function(j) {
      ind_j <- ind[j]:(ind[j + 1] - 1)
      Vj <- V[, j]
      U <- rotasym::r_unif_sphere(n = n, p = d[j]) %*%
        t(rotasym::Gamma_theta(theta = mu[ind_j]))
      Vj * matrix(mu[ind_j], nrow = n, ncol = d[j] + 1, byrow = TRUE) +
        sqrt(1 - Vj * Vj) * U
    }, simplify = FALSE)))

  } else {

    stop("kernel_type must be either 1 or 2.")

  }
  return(samp)

}


#' @title Sample from polyspherical kernel density estimator
#'
#' @description Simulates from the distribution defined by a polyspherical
#' kernel density estimator on \eqn{\mathcal{S}^{d_1} \times \ldots \times
#' \mathcal{S}^{d_r}}.
#'
#' @inheritParams kde_polysph
#' @param norm_X ensure a normalization of the data?
#' @inheritParams r_unif_polysph
#' @details The function uses \code{\link{r_kern_polysph}} to sample from the
#' considered kernel.
#' @return A matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @examples
#' # Simulated data on (S^1)^2
#' n <- 50
#' samp <- r_path_s1r(n = n, r = 2, k = c(1, 2), angles = TRUE)
#' plot(samp, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n),
#'      axes = FALSE, xlab = "", ylab = "", pch = 16, cex = 0.75)
#' points(torus_to_angles(r_kde_polysph(n = 10 * n, X = angles_to_torus(samp),
#'                                      d = c(1, 1), h = c(0.1, 0.1))),
#'        col = "black", pch = 16, cex = 0.2)
#' sdetorus::torusAxis()
#'
#' # Simulated data on S^2
#' n <- 50
#' samp <- r_path_s2r(n = n, r = 1, sigma = 0.1, kappa = 5,
#'                    spiral = TRUE)[, , 1]
#' sc3d <- scatterplot3d::scatterplot3d(
#'   samp, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' xyz <- r_kde_polysph(n = 10 * n, X = samp, d = 2, h = 0.1)
#' sc3d$points3d(xyz[, 1], xyz[, 2], xyz[, 3], col = "black", pch = 16,
#'               cex = 0.2)
#' @export
r_kde_polysph <- function(n, X, d, h, kernel = 1, kernel_type = 1, k = 10,
                          intrinsic = FALSE, norm_X = FALSE) {

  # Check dimensions
  r <- length(h)
  if (r != length(d)) {

    stop("h and d are incompatible.")

  }
  if (ncol(X) != sum(d + 1)) {

    stop("X and d are incompatible.")

  }

  # Index for accessing each S^dj with ind[j]:(ind[j + 1] - 1)
  ind <- cumsum(c(1, d + 1))

  # Normalize X?
  if (norm_X) {

    for (j in seq_len(r)) {

      ind_j <- ind[j]:(ind[j + 1] - 1)
      X[, ind_j] <- X[, ind_j] / sqrt(rowSums(X[, ind_j, drop = FALSE]^2))

    }

  }

  # Choose centers at random
  sample_ids <- sort(sample(x = nrow(X), size = n, replace = TRUE))
  ns <- rle(sample_ids)

  # Sample
  samp <- do.call(rbind, as.list(sapply(seq_along(ns$values), function(j)
    r_kern_polysph(n = ns$lengths[j], d = d, mu = X[ns$values[j], ,
                                                    drop = TRUE],
                   h = h, kernel = kernel, kernel_type = kernel_type, k = k,
                   intrinsic = intrinsic, norm_mu = FALSE),
    simplify = FALSE)))

  # Randomly permute
  return(samp[sample(nrow(samp)), , drop = FALSE])

}


#' @title Samplers of one-dimensional modes of variation for polyspherical data
#'
#' @description Functions for sampling data on \eqn{(\mathcal{S}^d)^r}, for
#' \eqn{d=1,2}, using one-dimensional modes of variation.
#'
#' @param n sample size.
#' @param r number of spheres in the polysphere \eqn{(\mathcal{S}^d)^r}.
#' @param alpha a vector of size \code{r} valued in \eqn{[-\pi,\pi)} with the
#' initial angles for the linear trend. Chosen at random by default.
#' @param k a vector of size \code{r} with the \bold{integer} slopes defining
#' the angular linear trend. Chosen at random by default.
#' @param sigma standard deviation of the noise about the one-dimensional mode
#' of variation. Defaults to \code{0.25}.
#' @param angles return angles in \eqn{[-\pi, \pi)}? Defaults to \code{FALSE}.
#' @param c \href{https://en.wikipedia.org/wiki/Clélie}{Clélie curve}
#' parameter, changing the spiral wrappings. Defaults to \code{1}.
#' @param t latitude, with respect to \code{Theta}, of the small circle.
#' Defaults to \code{0} (equator).
#' @param Theta a matrix of size \code{c(3, r)} giving the north poles for
#' \eqn{\mathcal{S}^2}. Useful for rotating the sample. Chosen at random by
#' default.
#' @param kappa concentration von Mises--Fisher parameter for longitudes in
#' small circles. Defaults to \code{0} (uniform).
#' @param spiral consider a spiral (or, more precisely, a
#' \href{https://en.wikipedia.org/wiki/Clélie}{Clélie curve}) instead of
#' a small circle? Defaults to \code{FALSE}.
#' @return
#' An array of size \code{c(n, d, r)} with samples on \eqn{(\mathcal{S}^d)^r}.
#' If \code{angles = TRUE} for \code{r_path_s1r}, then a matrix of size
#' \code{c(n ,r)} with angles is returned.
#' @examples
#' # Straight trends on (S^1)^2
#' n <- 100
#' samp_1 <- r_path_s1r(n = n, r = 2, k = c(1, 2), angles = TRUE)
#' plot(samp_1, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n),
#'      axes = FALSE, xlab = "", ylab = "", pch = 16)
#' sdetorus::torusAxis()
#'
#' # Straight trends on (S^1)^3
#' n <- 100
#' samp_2 <- r_path_s1r(n = n, r = 3, angles = TRUE)
#' pairs(samp_2, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n),
#'       pch = 16)
#' sdetorus::torusAxis()
#' scatterplot3d::scatterplot3d(
#'   samp_2, xlim = c(-pi, pi), ylim = c(-pi, pi), zlim = c(-pi, pi),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#'
#' # Small-circle trends on (S^2)^2
#' n <- 100
#' samp_3 <- r_path_s2r(n = n, r = 2, sigma = 0.1, kappa = 5)
#' old_par <- par(mfrow = c(1, 2))
#' scatterplot3d::scatterplot3d(
#'   samp_3[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' scatterplot3d::scatterplot3d(
#'   samp_3[, , 2], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' par(old_par)
#'
#' # Spiral trends on (S^2)^2
#' n <- 100
#' samp_4 <- r_path_s2r(n = n, r = 2, c = 3, spiral = TRUE, sigma = 0.01)
#' old_par <- par(mfrow = c(1, 2))
#' scatterplot3d::scatterplot3d(
#'   samp_4[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' scatterplot3d::scatterplot3d(
#'   samp_4[, , 2], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
#' )
#' par(old_par)
#' @export
r_path_s1r <- function(n, r, alpha = runif(r, -pi, pi),
                       k = sample(-2:2, size = r, replace = TRUE),
                       sigma = 0.25, angles = FALSE) {

  # Angular trend
  t <- sort(runif(n))
  theta <- t(alpha + (2 * pi * k) %o% t)

  # Add noise
  eps <- matrix(rnorm(r * n, sd = sigma), nrow = n, ncol = r)
  theta <- sdetorus::toPiInt(theta + eps)

  # Cartesianize
  if (angles) {

    return(theta)

  } else {

    return(sphunif::Theta_to_X(theta))

  }

}


#' @rdname r_path_s1r
#' @export
r_path_s2r <- function(n, r, t = 0, c = 1,
                       Theta = t(rotasym::r_unif_sphere(n = r, p = 3)),
                       kappa = 0, sigma = 0.25, spiral = FALSE) {

  samp <- array(dim = c(n, 3, r))
  if (spiral) {

    # Sample common latitude
    ph <- sort(runif(n, min = 0, max = pi))

    # Loop on the S^2's
    for (j in 1:r) {

      # Unrotated sample
      eps <- rnorm(n, sd = sigma)
      samp[, , j] <- DirStats::to_sph(th = c * ph + eps, ph = ph)

      # Rotated sample
      samp[, , j] <- samp[, , j] %*%
        cbind(Theta[, j], rotasym::Gamma_theta(theta = Theta[, j]))

    }

  } else {

    # Sample common and unrotated longitudes
    U <- rotasym::r_vMF(n = n, mu = c(1, 0), kappa = kappa)
    th <- sort(sdetorus::toPiInt(DirStats::to_rad(U)))
    U <- DirStats::to_cir(th = th)

    # Loop on the S^2's
    for (j in 1:r) {

      # Use the moment estimators of a Beta(shape1, shape2) to get the
      # parameters yielding
      # m = E[Beta(shape1, shape2)] = (t + 1) / 2 and
      # v = Var[Beta(shape1, shape2)] = sigma^2 / 4
      # (https://en.wikipedia.org/wiki/Beta_distribution#Two_unknown_parameters)
      # In this way,
      # E[2 * Beta(shape1, shape2) - 1] = t and
      # Var[2 * Beta(shape1, shape2) - 1] = sigma^2.
      m <- (t + 1) / 2
      v <- 0.25 * sigma^2
      stopifnot(v < m * (1 - m))
      shape1 <- m * (m * (1 - m) / v - 1)
      shape2 <- (1 - m) * (m * (1 - m) / v - 1)

      # Sample latitudes centered at t and with standard deviation sigma
      V <- 2 * rbeta(n = n, shape1 = shape1, shape2 = shape2) - 1

      # Tangent-normal decompositon about Theta[, j]
      samp[, , j] <- V * matrix(Theta[, j], nrow = n, ncol = 3, byrow = TRUE) +
        sqrt(1 - V^2) * U %*% t(rotasym::Gamma_theta(theta = Theta[, j]))

    }

  }

  return(samp)

}
