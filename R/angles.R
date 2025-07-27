
#' @title Conversion between the angular and Cartesian coordinates of the torus
#'
#' @description Transforms the angles \eqn{(\theta_1,\ldots,\theta_d)} in
#' \eqn{[-\pi,\pi)^d} into the Cartesian coordinates
#' \deqn{(\cos(x_1), \sin(x_1),\ldots,\cos(x_d), \sin(x_d))}
#' of the torus \eqn{(\mathcal{S}^1)^d}, and vice versa.
#'
#' @param theta matrix of size \code{c(n, d)} with the angles.
#' @param x matrix of size \code{c(n, 2 * d)} with the Cartesian coordinates on
#' \eqn{(\mathcal{S}^1)^d}. Assumed to be of unit norm by pairs of coordinates
#' in the rows.
#' @return
#' \itemize{
#' \item{\code{angles_to_torus}: the matrix \code{x}.}
#' \item{\code{torus_to_angles}: the matrix \code{theta}.}
#' }
#' @examples
#' # Check changes of coordinates
#' torus_to_angles(angles_to_torus(c(0, pi / 3, pi / 2)))
#' torus_to_angles(angles_to_torus(rbind(c(0, pi / 3, pi / 2), c(0, 1, -2))))
#' angles_to_torus(torus_to_angles(c(0, 1, 1, 0)))
#' angles_to_torus(torus_to_angles(rbind(c(0, 1, 1, 0), c(0, 1, 0, 1))))
#' @export
angles_to_torus <- function(theta) {

  # Convert to matrix
  if (!is.matrix(theta)) {

    theta <- matrix(theta, nrow = 1)

  }

  # Apply transformation
  d <- ncol(theta)
  ind <- c(rbind(seq_len(d), seq(d + 1, 2 * d, by = 1)))
  return(cbind(cos(theta), sin(theta))[, ind, drop = FALSE])

}


#' @rdname angles_to_torus
#' @export
torus_to_angles <- function(x) {

  # Convert to matrix
  if (!is.matrix(x)) {

    x <- matrix(x, nrow = 1)

  }

  # Check p
  p <- ncol(x)
  if (p %% 2) {

    stop("The dimension is not even.")

  }

  # Apply transformation
  n <- nrow(x)
  ind_cos <- seq(1, p, by = 2)
  matrix(atan2(y = x[, ind_cos + 1], x = x[, ind_cos]), nrow = n,
         ncol = p / 2L)[, , drop = FALSE]

}


#' @title Conversion between the angular and Cartesian coordinates of the
#' (hyper)sphere
#'
#' @description Transforms the angles \eqn{(\theta_1,\ldots,\theta_d)} in
#' \eqn{[0,\pi)^{d-1}\times[-\pi,\pi)} into the Cartesian coordinates
#' \deqn{(\cos(x_1),\sin(x_1)\cos(x_2),\ldots,
#' \sin(x_1)\cdots\sin(x_{d-1})\cos(x_d),
#' \sin(x_1)\cdots\sin(x_{d-1})\sin(x_d))}
#' of the sphere \eqn{\mathcal{S}^{d}}, and vice versa.
#'
#' @param theta matrix of size \code{c(n, d)} with the angles.
#' @param x matrix of size \code{c(n, d + 1)} with the Cartesian coordinates
#' on \eqn{\mathcal{S}^{d}}. Assumed to be of unit norm by rows.
#' @return
#' \itemize{
#' \item{\code{angles_to_sph}: the matrix \code{x}.}
#' \item{\code{sph_to_angles}: the matrix \code{theta}.}
#' }
#' @examples
#' # Check changes of coordinates
#' sph_to_angles(angles_to_sph(c(pi / 2, 0, pi)))
#' sph_to_angles(angles_to_sph(rbind(c(pi / 2, 0, pi), c(pi, pi / 2, 0))))
#' angles_to_sph(sph_to_angles(c(0, sqrt(0.5), sqrt(0.1), sqrt(0.4))))
#' angles_to_sph(sph_to_angles(rbind(c(0, sqrt(0.5), sqrt(0.1), sqrt(0.4)),
#'                                   c(0, sqrt(0.5), sqrt(0.5), 0),
#'                                   c(0, 1, 0, 0),
#'                                   c(0, 0, 0, -1),
#'                                   c(0, 0, 1, 0))))
#' @export
angles_to_sph <- function(theta) {

  # Convert to matrix
  if (!is.matrix(theta)) {

    theta <- matrix(theta, nrow = 1)

  }

  # Apply transformation
  d <- ncol(theta)
  x <- matrix(0, nrow = nrow(theta), ncol = d + 1)
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)
  if (d == 1) {

    x <- cbind(cos_theta, sin_theta, deparse.level = 0)

  } else {

    prod_sin <- t(apply(sin_theta, 1, cumprod))
    x <- cbind(cos_theta[, 1, drop = FALSE],
               prod_sin[, -d, drop = FALSE] * cos_theta[, -1, drop = FALSE],
               prod_sin[, d, drop = FALSE], deparse.level = 0)[, , drop = FALSE]

  }
  return(x)

}


#' @rdname angles_to_sph
#' @export
sph_to_angles <- function(x) {

  # Convert to matrix
  if (!is.matrix(x)) {

    x <- matrix(x, nrow = 1)

  }

  # Apply transformation
  d <- ncol(x) - 1
  i_norm <- t(apply(x^2, 1, function(x) rev(sqrt(cumsum(rev(x))))))[, -(d + 1)]
  theta <- x[, -c(d + 1), drop = FALSE] / i_norm
  theta[is.nan(theta)] <- 1 # Avoid NaNs
  theta <- acos(theta)
  xd <- (x[, d + 1] < 0)
  theta[, d] <- (2 * pi) * xd + (1 - 2 * xd) * theta[, d]
  return(unname(theta[, , drop = FALSE]))

}


#' @title Conversion between the angular and Cartesian coordinates of the
#' polysphere
#'
#' @description Obtain the angular coordinates of points on a polysphere
#' \eqn{\mathcal{S}^{d_1}\times\cdots\times\mathcal{S}^{d_r}}, and vice versa.
#'
#' @param x matrix of size \code{c(n, sum(d + 1))} with the Cartesian
#' coordinates on \eqn{\mathcal{S}^{d_1}\times\cdots\times\mathcal{S}^{d_r}}.
#' Assumed to be of unit norm by blocks of coordinates in the rows.
#' @param theta matrix of size \code{c(n, sum(d))} with the angles.
#' @param d vector with the dimensions of the polysphere.
#' @return
#' \itemize{
#' \item{\code{angles_to_polysph}: the matrix \code{x}.}
#' \item{\code{polysph_to_angles}: the matrix \code{theta}.}
#' }
#' @examples
#' # Check changes of coordinates
#' polysph_to_angles(angles_to_polysph(rep(pi / 2, 3), d = 2:1), d = 2:1)
#' angles_to_polysph(polysph_to_angles(x = c(0, 0, 1, 0, 1), d = 2:1), d = 2:1)
#' @export
angles_to_polysph <- function(theta, d) {

  # Convert to matrix
  if (!is.matrix(theta)) {

    theta <- matrix(theta, nrow = 1)

  }

  # Call angles_to_sph() for each sphere
  x <- matrix(0, nrow = nrow(theta), ncol = length(d) + sum(d))
  ind_p <- 0
  ind_r <- 0
  for (di in d) {

    input <- as.matrix(theta[, (ind_p + 1):(ind_p + di)])
    if (nrow(theta) == 1) {

      input <- t(input)

    }

    x[, (ind_r + 1):(ind_r + di + 1)] <- angles_to_sph(input)
    ind_p <- ind_p + di
    ind_r <- ind_r + di + 1

  }
  return(x)

}


#' @rdname angles_to_polysph
#' @export
polysph_to_angles <- function(x, d) {

  # Convert to matrix
  if (!is.matrix(x)) {

    x <- matrix(x, nrow = 1)

  }

  # Call sph_to_angles() for each sphere
  theta <- matrix(0, nrow = nrow(x), ncol = sum(d))
  ind_p <- 0
  ind_r <- 0
  for (di in d) {

    input <- as.matrix(x[, (ind_p + 1):(ind_p + di + 1)])
    if (nrow(x) == 1) {

      input <- t(input)
    }

    theta[, (ind_r + 1):(ind_r + di)] <- sph_to_angles(input)
    ind_p <- ind_p + di + 1
    ind_r <- ind_r + di

  }
  return(theta)

}


#' @title Fibonacci lattice on the sphere
#'
#' @description Computes the Fibonacci lattice on the sphere \eqn{\mathcal{S}^2}
#' to produce pseudo-equispaced points.
#' @param n number of points to be produced.
#' @return A matrix of size \code{c(n, 3)} with the spherical coordinates.
#' @examples
#' scatterplot3d::scatterplot3d(
#'   fib_latt(n = 200), pch = 19, color = viridis::viridis(200),
#'   xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
#'   xlab = "", ylab = "", zlab = "")
#' @export
fib_latt <- function(n) {

  # Check https://extremelearning.com.au/how-to-evenly-distribute-points-on-a-
  # sphere-more-effectively-than-the-canonical-fibonacci-lattice/
  phi <- (1 + sqrt(5)) / 2
  i <- 0:(n - 1)
  theta <- 2 * pi * i / phi
  phi <- acos(1 - 2 * (i + 0.5) / n)
  return(angles_to_sph(cbind(phi, theta, deparse.level = 0)))

}


#' @title Map Cartesian coordinates into Hammer projection
#'
#' @description Computes the Hammer projection of points on the sphere.
#' @param x matrix of size \code{c(n, 3)} with the Cartesian coordinates on
#' \eqn{\mathcal{S}^2}. Assumed to be of unit norm in the rows.
#' @param y matrix of size \code{c(n, 2)} with the Hammer coordinates.
#' @return
#' \itemize{
#' \item{\code{sph_to_hammer}: the matrix \code{y}.}
#' \item{\code{hammer_to_sph}: the matrix \code{x}.}
#' }
#' @examples
#' # Plot Fibonacci lattice
#' plot(sph_to_hammer(fib_latt(n = 1000)))
#' points(sph_to_hammer(rbind(c(0, 0, 1), c(0, 0, -1))), col = 2, pch = 19)
#' points(sph_to_hammer(rbind(c(1, 0, 0), c(-1, 0, 0))), col = 3, pch = 19)
#' points(sph_to_hammer(rbind(c(0, 1, 0), c(0, -1, 0))), col = 4, pch = 19)
#'
#' # Check changes of coordinates
#' hammer_to_sph(sph_to_hammer(rbind(c(1, 0, 0), c(0, 1, 0))))
#' sph_to_hammer(hammer_to_sph(rbind(c(0, 0), c(0.5, 0.5))))
#' @name hammer_to_sph


#' @rdname hammer_to_sph
#' @export
sph_to_hammer <- function(x) {

  # Convert to matrix
  if (!is.matrix(x)) {

    x <- matrix(x, nrow = 1)

  }
  stopifnot(ncol(x) == 3)

  # Reverse with 3:1 and shift (lambda, phi) because sph_to_angles() has the
  # spherical coordinates in reverse order and the transformation is for
  # standard spherical coordinates with \lambda \in [-pi, pi] and
  # \phi \in [-pi / 2, pi / 2].
  phi_theta <- sph_to_angles(x[, 3:1])
  lambda <- phi_theta[, 2] - pi # Longitude, in [-pi, pi]
  phi <- pi / 2 - phi_theta[, 1] # Latitude, in [-pi / 2, pi / 2]
  x <- 2 * sqrt(2) * cos(phi) * sin(lambda / 2)
  y <- sqrt(2) * sin(phi)
  cbind(x, y, deparse.level = 0) / sqrt(1 + cos(phi) * cos(lambda / 2))

}


#' @rdname hammer_to_sph
#' @export
hammer_to_sph <- function(y) {

  # Convert to matrix
  if (!is.matrix(y)) {

    y <- matrix(y, nrow = 1)

  }
  stopifnot(ncol(y) == 2)

  # Reverse with 3:1 and shift (lambda, phi) because sph_to_angles() has the
  # spherical coordinates in reverse order and the inverse transformation is for
  # standard spherical coordinates with \lambda \in [-pi, pi] and
  # \phi \in [-pi / 2, pi / 2].
  x <- y[, 1]
  y <- y[, 2]
  z <- sqrt(1 - (0.25 * x)^2 - (0.5 * y)^2)
  lambda <- 2 * atan(z * x / (2 * (2 * z^2 - 1)))
  phi <- asin(z * y)
  angles_to_sph(cbind(pi / 2 - phi, lambda + pi,
                      deparse.level = 0))[, 3:1, drop = FALSE]

}
