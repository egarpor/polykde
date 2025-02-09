
r <- rpois(1, lambda = 15) + 1
d <- rep(sample(1:10, size = 1), r)
h <- runif(r, min = 0.2)
X <- r_unif_polysph(n = 3, d = d)
x <- r_unif_polysph(n = 2, d = d)
ind_blocks <- sample(1:3, size = r, replace = TRUE)
ind_blocks <- as.numeric(as.factor(ind_blocks))

test_that("Reordering of blocks works in block_euler_ridge()", {
  res <- block_euler_ridge(x = x, X = X, d = d, h = h,
                           h_euler = h, ind_blocks = ind_blocks,
                           N = 1, keep_paths = TRUE,
                           show_prog = FALSE)
  expect_equal(h, res$h)
  expect_equal(d, res$d)
  expect_equal(x, res$start_x)
  expect_equal(x, res$paths[, , 1])
})

# Transformation functions
s1r_to_angles <- function(x) {

  stopifnot(is.matrix(x))
  r <- ncol(x) / 2
  d <- rep(1, r)
  ind <- cumsum(c(1, d + 1))
  th <- matrix(nrow = nrow(x), ncol = r)
  for (j in seq_along(d)) {

    ind_j <- ind[j]:(ind[j + 1] - 1)
    th[, j] <- DirStats::to_rad(x[, ind_j])

  }
  return(sdetorus::toPiInt(th))

}
angles_to_s1r <- function(th) {

  stopifnot(is.matrix(th))
  x <- sphunif::Theta_to_X(th)
  x <- do.call(cbind, lapply(seq_len(ncol(th)), function(i) rbind(x[, , i])))
  return(x)

}

# Visualize (S^1)^r Euler
viz_euler_s1r <- function(y_old, X, col_X, k = 1, ind_dj, d, h, h_euler,
                          N = 100, eps = 1e-6, frame = TRUE, grad = FALSE,
                          rays = FALSE) {

  # Kde
  stopifnot(all(d == 1))
  r <- length(d)
  # kde_X <- kde_polysph(x = X, X = X, d = d, h = h)
  # cut_X <- cut(kde_X, breaks = seq(min(kde_X), max(kde_X), l = 21),
  #              include.lowest = TRUE)
  # col_X <- viridis::viridis(20, alpha = 0.15)[cut_X]

  # Set index and plot sample
  ind_k <- 1:2 + 2 * (k - 1)
  plot(X[, ind_k], col = col_X, pch = 19, xlim = c(-1.25, 1.25),
       ylim = c(-1.25, 1.25), xlab = "x", ylab = "y")
  lines(DirStats::to_cir(seq(0, 2 * pi, l = 200)), col = "lightblue")

  # Begin
  y_old_k <- y_old[, ind_k, drop = FALSE]
  points(y_old_k, col = 2, pch = 16)

  # Euler
  path <- rbind(y_old)
  delta <- Inf
  j <- 1
  while (j <= N) {

    # Plot y
    y_old_k <- y_old[, ind_k, drop = FALSE]
    points(y_old_k, col = 2, pch = 16)
    if (rays) {

      lines(rbind(0, y_old_k), col = 2)

    }

    # Projected normalized gradient
    proj <- proj_grad_kde_polysph(x = y_old, X = X, d = d, h = h,
                                  wrt_unif = FALSE, norm_x = FALSE,
                                  norm_X = FALSE, kernel = 1, k = 1,
                                  sparse = FALSE)
    gh <- grad_hess_kde_polysph(x = y_old, X = X, d = d, h = h,
                                projected = TRUE, norm_grad_hess = FALSE)

    # Eigenvectors Hessian
    H <- gh$hess[1, , ]
    eig <- eigen(H)
    where_null <- abs(eig$values) < 1e-10
    stopifnot(sum(where_null) == r)
    eig$vectors <- eig$vectors[, c(which(!where_null), which(where_null))]
    eig$values <- eig$values[c(which(!where_null), which(where_null))]
    cat("j =", j, "\n")
    cat("y_k =", y_old_k, "\n")
    cat("vectors = \n")
    print(eig$vectors)
    cat("values = \n")
    print(eig$values)
    cat("y =", y_old, "\n")
    cat("u1_proj =", proj$u1, "\n")
    stopifnot(any(drop(1 - abs(proj$u1 %*% eig$vectors)) < 1e-10))

    # Draw frame of eigenvectors?
    if (frame) {

      for (i in seq_len(ncol(X))) {

        eig_k <- t(eig$vectors[ind_k, i])
        eig_k <- eig_k / sqrt(sum(eig_k^2))
        arrows(x0 = y_old_k[1], y0 = y_old_k[2],
               x1 = y_old_k[1] + 0.25 * eig_k[1],
               y1 = y_old_k[2] + 0.25 * eig_k[2], col = (i == 1) + 1,
               length = 0.2, angle = 10)
        text(x = y_old_k + 0.3 * eig_k, labels =
               substitute(u[ii], list(ii = i)), col = (i == 1) + 1)

      }

    }

    # Unprojected gradient
    if (grad) {

      gk <- gh$grad[ind_k]
      gk <- gk / sqrt(sum(gk^2))
      arrows(x0 = y_old_k[1], y0 = y_old_k[2], x1 = y_old_k[1] + 0.25 * gk[1],
             y1 = y_old_k[2] + 0.25 * gk[2], col = 3, length = 0.2, angle = 10)
      text(x = y_old_k + 0.3 * gk, labels = expression(nabla), col = 3)

    }

    # Projected gradient eta
    eta <- proj$eta
    cat("eta_proj =", eta, "\n")
    eta_k <- eta[ind_k]
    h_euler_k <- h_euler[k]
    arrows(x0 = y_old_k[1], y0 = y_old_k[2],
           x1 = y_old_k[1] + h_euler_k[1] * eta_k[1],
           y1 = y_old_k[2] + h_euler_k[2] * eta_k[2], col = 4,
           length = 0.2, angle = 10)
    text(x = y_old_k + 1.1 * h_euler_k * eta_k, labels = expression(nabla[2]),
         col = 4)

    # Advance
    y_new <- y_old + matrix(t(h_euler * matrix(eta, ncol = 2 * r,
                                               byrow = TRUE)), nrow = 1)
    y_new_k <- y_new[, ind_k, drop = FALSE]
    points(y_new_k, col = 4, pch = 16)

    # Project to polysphere
    y_new2 <- proj_polysph(x = y_new, ind_dj = ind_dj)
    y_new_k2 <- y_new2[, ind_k, drop = FALSE]
    lines(rbind(y_new_k, y_new_k2), col = 4)
    points(y_new_k2, col = 4, pch = 16)
    y_new <- y_new2
    y_new_k <- y_new[, ind_k, drop = FALSE]

    # Save path
    path <- rbind(path, y_new)

    # Convergence?
    delta <- c(delta, dist_polysph(x = y_old, y = y_new, ind_dj = ind_dj))
    cat("delta =", delta[j], "\n")
    if (delta[j] < eps) {

      # Stop
      cat("Converged!\n")
      break

    } else {

      # Update
      y_old <- y_new

    }

    # Increase j
    j <- j + 1

  }

  # End
  points(y_new_k, col = 6, pch = 16)
  # path_k <- path[, ind_k]
  # lines(path_k, col = 6, lwd = 5)

}

# Visualize T^2 Euler
viz_euler_t2 <- function(y_old, X, col_X, d, h, h_euler, N = 100, eps = 1e-6,
                         frame = TRUE, grad = FALSE) {

  # Auxiliary function to project on T^2
  v_to_ang <- function(v) {

    sig <- c(prod(sign(v[1:2])), prod(sign(v[3:4])))
    sig * c(sum(v[1:2]^2), sum(v[3:4]^2))

  }

  # Kde
  stopifnot(all(d == 1))
  stopifnot(ncol(X) == 4)
  # kde_X <- kde_polysph(x = X, X = X, d = d, h = h)
  # cut_X <- cut(kde_X, breaks = seq(min(kde_X), max(kde_X), l = 21),
  #              include.lowest = TRUE)
  # col_X <- viridis::viridis(20, alpha = 0.15)[cut_X]

  # Plot sample
  Th <- s1r_to_angles(x = X)
  plot(Th, col = col_X, pch = 19, xlim = c(-pi, pi), ylim = c(-pi, pi),
       xlab = "x", ylab = "y", axes = FALSE)
  sdetorus::torusAxis()

  # Begin
  th_old <- s1r_to_angles(y_old)
  points(th_old, col = 2, pch = 16)

  # Euler
  path <- rbind(th_old)
  delta <- Inf
  j <- 1
  while (j <= N) {

    # Plot y
    th_old <- s1r_to_angles(y_old)
    points(th_old, col = 2, pch = 16)

    # Projected normalized gradient
    proj <- proj_grad_kde_polysph(x = y_old, X = X, d = d, h = h,
                                  wrt_unif = FALSE, norm_x = FALSE,
                                  norm_X = FALSE, kernel = 1, k = 1,
                                  sparse = FALSE)
    gh <- grad_hess_kde_polysph(x = y_old, X = X, d = d, h = h,
                                projected = TRUE,
                                norm_grad_hess = FALSE)

    # Eigenvectors Hessian
    H <- gh$hess[1, , ]
    eig <- eigen(H)
    where_null <- abs(eig$values) < 1e-10
    stopifnot(sum(where_null) == r)
    eig$vectors <- eig$vectors[, c(which(!where_null), which(where_null))]
    eig$values <- eig$values[c(which(!where_null), which(where_null))]
    cat("j =", j, "\n")
    cat("th =", th_old, "\n")
    cat("vectors = \n")
    print(eig$vectors)
    cat("values = \n")
    print(eig$values)
    cat("y =", y_old, "\n")
    cat("u1_proj =", proj$u1, "\n")
    stopifnot(any(drop(1 - abs(proj$u1 %*% eig$vectors)) < 1e-10))

    # Draw frame of eigenvectors?
    if (frame) {

      for (i in 1:2) {

        eig_k <- t(eig$vectors[, i])
        eig_k <- v_to_ang(eig_k)
        arrows(x0 = th_old[1], y0 = th_old[2],
               x1 = th_old[1] + 0.25 * eig_k[1],
               y1 = th_old[2] + 0.25 * eig_k[2], col = (i == 1) + 1,
               length = 0.2, angle = 10)
        text(x = th_old + 0.3 * eig_k, labels =
               substitute(u[ii], list(ii = i)), col = (i == 1) + 1)

      }

    }

    # Unprojected gradient
    if (grad) {

      gk <- gh$grad
      gk <- v_to_ang(gk)
      arrows(x0 = th_old[1], y0 = th_old[2], x1 = th_old[1] + 0.25 * gk[1],
             y1 = th_old[2] + 0.25 * gk[2], col = 3, length = 0.2, angle = 10)
      text(x = th_old + 0.3 * gk, labels = expression(nabla), col = 3)

    }

    # Projected gradient eta
    eta <- proj$eta
    cat("eta_proj =", eta, "\n")
    eta_k <- v_to_ang(eta)
    arrows(x0 = th_old[1], y0 = th_old[2],
           x1 = th_old[1] + h_euler[1] * eta_k[1],
           y1 = th_old[2] + h_euler[2] * eta_k[2], col = 4,
           length = 0.2, angle = 10)
    text(x = th_old + 1.1 * h_euler * eta_k, labels = expression(nabla[2]),
         col = 4)

    # Advance
    y_new <- y_old + matrix(t(h_euler * matrix(eta, ncol = 4, byrow = TRUE)),
                            nrow = 1)
    y_new <- proj_polysph(x = y_new, ind_dj = ind_dj)
    th_new <- s1r_to_angles(y_new)
    points(th_new, col = 4, pch = 16)

    # Save path
    path <- rbind(path, th_new)

    # Convergence?
    delta <- c(delta, dist_polysph(x = y_old, y = y_new, ind_dj = ind_dj))
    cat("delta =", delta[j], "\n")
    if (delta[j] < eps) {

      # Stop
      cat("Converged!\n")
      break

    } else {

      # Update
      y_old <- y_new

    }

    # Increase j
    j <- j + 1

  }

  # End
  points(th_new, col = 6, pch = 16)
  # path_k <- path[, ind_k]
  # lines(path_k, col = 6, lwd = 5)

}

# Visualize (S^2)^r Euler
viz_euler_s2r <- function(y_old, X, col_X, k = 1, ind_dj, d, h, h_euler,
                          N = 100, eps = 1e-6, frame = TRUE, grad = FALSE,
                          rays = FALSE, ...) {

  # Kde
  stopifnot(all(d == 2))
  r <- length(d)
  # kde_X <- kde_polysph(x = X, X = X, d = d, h = h)
  # cut_X <- cut(kde_X, breaks = seq(min(kde_X), max(kde_X), l = 21),
  #              include.lowest = TRUE)
  # col_X <- viridis::viridis(20, alpha = 0.15)[cut_X]

  # Set index and plot sample
  ind_k <- 1:3 + 3 * (k - 1)
  rgl::plot3d(X[, ind_k], col = col_X, size = 5, xlim = c(-1, 1),
              ylim = c(-1, 1), zlim = c(-1, 1), xlab = "", ylab = "",
              zlab = "", ...)
  rgl::spheres3d(0, 0, 0, radius = 1, col = "lightblue", lit = FALSE,
                 alpha = 0.5)

  # Begin
  y_old_k <- y_old[, ind_k, drop = FALSE]
  rgl::points3d(y_old_k, col = 2, size = 10)

  # Euler
  path <- rbind(y_old)
  delta <- Inf
  j <- 1
  while (j <= N) {

    # Plot y
    y_old_k <- y_old[, ind_k, drop = FALSE]
    rgl::points3d(y_old_k, col = 2, size = 5)
    if (rays) {

      rgl::lines3d(rbind(0, y_old_k), col = 2)

    }

    # Projected normalized gradient
    proj <- proj_grad_kde_polysph(x = y_old, X = X, d = d, h = h,
                                  wrt_unif = FALSE, norm_x = FALSE,
                                  norm_X = FALSE, kernel = 1, k = 1,
                                  sparse = FALSE)
    gh <- grad_hess_kde_polysph(x = y_old, X = X, d = d, h = h,
                                projected = TRUE, norm_grad_hess = FALSE)

    # Eigenvectors Hessian
    H <- gh$hess[1, , ]
    eig <- eigen(H)
    where_null <- abs(eig$values) < 1e-10
    eig$vectors <- eig$vectors[, c(which(!where_null), which(where_null))]
    eig$values <- eig$values[c(which(!where_null), which(where_null))]
    cat("j =", j, "\n")
    cat("y_k =", y_old_k, "\n")
    cat("vectors = \n")
    print(eig$vectors)
    cat("values = \n")
    print(eig$values)
    cat("y =", y_old, "\n")
    cat("u1_proj =", proj$u1, "\n")
    stopifnot(any(drop(1 - abs(proj$u1 %*% eig$vectors)) < 1e-10))

    # Draw frame of eigenvectors?
    if (frame) {

      for (i in seq_len(ncol(X))) {

        eig_k <- t(eig$vectors[ind_k, i])
        eig_k <- eig_k / sqrt(sum(eig_k^2))
        rgl::arrow3d(p0 = y_old_k, p1 = y_old_k + 0.25 * eig_k,
                     col = (i == 1) + 1, width = 0.2, type = "lines",
                     theta = 0.1)
        rgl::text3d(x = y_old_k + 0.3 * eig_k, texts =
                      substitute(u[ii], list(ii = i)), col = (i == 1) + 1,
                    usePlotmath = TRUE)

      }

    }

    # Unprojected gradient
    if (grad) {

      gk <- gh$grad[ind_k]
      gk <- gk / sqrt(sum(gk^2))
      rgl::arrow3d(p0 = y_old_k, p1 = y_old_k + 0.25 * gk, col = 3, width = 0.2,
                   type = "lines", theta = 0.05)
      rgl::text3d(x = y_old_k + 0.3 * gk,
                  texts = expression(nabla), usePlotmath = TRUE, col = 3)

    }

    # Projected gradient eta
    eta <- proj$eta
    cat("eta_proj =", eta, "\n")
    eta_k <- eta[ind_k]
    h_euler_k <- h_euler[k]
    rgl::arrow3d(p0 = y_old_k, p1 = y_old_k + h_euler_k * eta_k,
                 col = 4, width = 0.2, type = "lines", theta = 0.05, lwd = 2)
    rgl::text3d(x = y_old_k + 1.1 * h_euler_k * eta_k,
                texts = expression(nabla[2]), usePlotmath = TRUE, col = 4)

    # Advance
    y_new <- y_old + matrix(t(h_euler * matrix(eta, ncol = 3 * r,
                                               byrow = TRUE)), nrow = 1)
    y_new_k <- y_new[, ind_k, drop = FALSE]
    rgl::points3d(y_new_k, col = 4, size = 5)

    # Project to polysphere
    y_new2 <- proj_polysph(x = y_new, ind_dj = ind_dj)
    y_new_k2 <- y_new2[, ind_k, drop = FALSE]
    rgl::lines3d(rbind(y_new_k, y_new_k2), col = 4)
    rgl::points3d(y_new_k2, col = 4, size = 5)
    y_new <- y_new2
    y_new_k <- y_new[, ind_k, drop = FALSE]

    # Save path
    path <- rbind(path, y_new)

    # Convergence?
    delta <- c(delta, dist_polysph(x = y_old, y = y_new, ind_dj = ind_dj))
    cat("delta =", delta[j], "\n")
    if (delta[j] < eps) {

      # Stop
      cat("Converged!\n")
      break

    } else {

      # Update
      y_old <- y_new

    }

    # Increase j
    j <- j + 1

  }

  # End
  rgl::points3d(y_new_k, col = 6, size = 10)
  # path_k <- path[, ind_k]
  # rgl::lines3d(path_k, col = 6, lwd = 5)

}

skip("No tests for euler ridges, just visualizations")

## (S^1)^2 test

# Sample
# set.seed(132121)
r <- 2
d <- rep(1, r)
n <- 200
ind_dj <- comp_ind_dj(d = d)
Th <- r_path_s1r(n = n, r = r, angles = TRUE, k = c(1, 2), sigma = 0.5)
X <- angles_to_s1r(Th)
col_X_alp <- viridis::viridis(n, alpha = 0.25)
col_X <- viridis::viridis(n)

# Single starting point
i <- 50
y <- X[i, , drop = FALSE]
# y <- angles_to_s1r(th = rbind(c(1, 2)))
col_X_alp <- c(col_X_alp, 1)
col_X <- c(col_X, 1)

# Euler
h_rid <- rep(0.75, r)
h_eu <- h_rid^2
N <- 200
eps <- 1e-6
Xy <- rbind(X, y)
Y <- euler_ridge(x = Xy, X = X, d = d, h = h_rid, h_euler = h_eu,
                 N = N, eps = eps, keep_paths = TRUE)
# kde_X <- kde_polysph(x = Xy, X = X, d = d, h = h_rid)
# cut_X <- cut(kde_X, breaks = seq(min(kde_X), max(kde_X), l = 21),
#              include.lowest = TRUE)
# col_X_alp <- viridis::viridis(20, alpha = 0.25)[cut_X]
# col_X <- viridis::viridis(20)[cut_X]

# Kde in grid
L <- 50
tth <- seq(-pi, pi, l = L + 1)[-(L + 1)]
th <- as.matrix(expand.grid(tth, tth))
x <- angles_to_s1r(th)
kde_tth <- matrix(drop(kde_polysph(x = x, X = X, d = d, h = h_rid)),
                  nrow = L, ncol = L)

# Dynamic visualization
manipulate::manipulate({
# for (i in seq_len(dim(Y$paths)[3])) {
# pdf(paste0("S12_euler_", i, ".pdf"), width = 7, height = 7)

  contour(tth, tth, kde_tth, axes = FALSE)
  sdetorus::torusAxis(1:2)
  points(rbind(s1r_to_angles(Y$paths[, , 1])), col = col_X_alp, pch = 19)
  points(rbind(s1r_to_angles(Y$paths[, , i])), col = col_X, pch = 16,
         cex = 0.75)

  for (k in seq_len(nrow(Y$paths))) {

    xy <- s1r_to_angles(t(Y$paths[k, , ]))
    sdetorus::linesTorus(x = xy[, 1], y = xy[, 2], col = col_X_alp[k])

  }

# dev.off()
# }
}, i = manipulate::slider(1, dim(Y$paths)[3]))

# Visualization on T^2
viz_euler_t2(y_old = y, X = X, col_X = col_X_alp, d = d, h = h_rid,
             h_euler = h_eu, N = N, eps = eps, frame = TRUE, grad = FALSE)
lines(s1r_to_angles(t(Y$paths[n + 1, , ])), col = 6, lwd = 5)

# Visualization on (S^1)^2
viz_euler_s1r(y_old = y, X = X, col_X = col_X_alp, k = 1, ind_dj = ind_dj,
              d = d, h = h_rid, h_euler = h_eu, N = N, eps = eps, frame = TRUE,
              grad = FALSE, rays = FALSE)
lines(t(Y$paths[n + 1, 1:2, ]), col = 6, lwd = 5)
viz_euler_s1r(y_old = y, X = X, col_X = col_X_alp, k = 2, ind_dj = ind_dj,
              d = d, h = h_rid, h_euler = h_eu, N = N, eps = eps, frame = TRUE,
              grad = FALSE, rays = FALSE)
lines(t(Y$paths[n + 1, 3:4, ]), col = 6, lwd = 5)

## (S^1)^3 test

# Sample
# set.seed(132121)
r <- 3
d <- rep(1, r)
n <- 200
ind_dj <- comp_ind_dj(d = d)
Th <- r_path_s1r(n = n, r = r, angles = TRUE, k = c(2, 1, 1), sigma = 0.5)
X <- angles_to_s1r(Th)
col_X_alp <- viridis::viridis(n, alpha = 0.25)
col_X <- viridis::viridis(n)

# Single starting point
i <- 10
y <- X[i, , drop = FALSE]
# y <- angles_to_s1r(th = rbind(c(1, 2)))
col_X_alp <- c(col_X_alp, 1)
col_X <- c(col_X, 1)

# Euler
h_rid <- rep(0.65, r)
h_eu <- h_rid^2
N <- 200
eps <- 1e-6
Xy <- rbind(X, y)
Y <- euler_ridge(x = Xy, X = X, d = d, h = h_rid, h_euler = h_eu,
                 N = N, eps = eps, keep_paths = TRUE)
# kde_X <- kde_polysph(x = Xy, X = X, d = d, h = h_rid)
# cut_X <- cut(kde_X, breaks = seq(min(kde_X), max(kde_X), l = 21),
#              include.lowest = TRUE)
# col_X_alp <- viridis::viridis(20, alpha = 0.25)[cut_X]
# col_X <- viridis::viridis(20)[cut_X]

# Dynamic visualization
manipulate::manipulate({

  sc3d <- scatterplot3d::scatterplot3d(rbind(s1r_to_angles(Y$paths[, , 1])),
                               xlim = c(-pi, pi), ylim = c(-pi, pi),
                               zlim = c(-pi, pi), color = col_X_alp, pch = 19,
                               angle = 20, xlab = expression(theta[1]),
                               ylab = expression(theta[2]),
                               zlab = expression(theta[3]))
  sc3d$points3d(rbind(s1r_to_angles(Y$paths[, , i])), col = col_X, pch = 16,
                cex = 0.75)

  # for (k in seq_len(nrow(Y$paths))) {
  #
  #   xy <- s1r_to_angles(t(Y$paths[k, , ]))
  #   sc3d$points3d(x = xy[, 1], y = xy[, 2], z = xy[, 3],
  #                 col = col_X_alp[k], type = "l")
  #
  # }

# dev.off()
# }
}, i = manipulate::slider(1, dim(Y$paths)[3]))

## S^2 test

# Sample
# set.seed(3242)
r <- 1
d <- 2
n <- 200
ind_dj <- comp_ind_dj(d = d)
X <- r_path_s2r(n = n, r = r, spiral = FALSE, Theta = cbind(c(1, 0, 0)),
                sigma = 0.2)[, , 1]
col_X_alp <- viridis::viridis(n, alpha = 0.25)
col_X <- viridis::viridis(n)

# Single starting point
i <- 10
y <- X[i, , drop = FALSE]
y <- rbind(c(1, 0, 0))
col_X_alp <- c(col_X_alp, 1)
col_X <- c(col_X, 1)

# Euler and kde
h_rid <- 2
h_eu <- h_rid^2
N <- 200
eps <- 1e-6
Xy <- rbind(X, y)
Y <- euler_ridge(x = Xy, X = X, d = d, h = h_rid, h_euler = h_eu,
                 N = N, eps = eps, keep_paths = TRUE)
# kde_X <- kde_polysph(x = Xy, X = X, d = d, h = h_rid)
# cut_X <- cut(kde_X, breaks = seq(min(kde_X), max(kde_X), l = 21),
#              include.lowest = TRUE)
# col_X_alp <- viridis::viridis(20, alpha = 0.25)[cut_X]
# col_X <- viridis::viridis(20)[cut_X]


# # Clean and index ridge
# Y <- clean_euler_ridge(e = Y, X = X, p_out = 0.01)
# ridge <- index_ridge(endpoints = Y$ridge_y, X = X, d = d, verbose = TRUE,
#                      type_bwd = "1se", probs_scores = seq(0, 1, l = 101))
#

# Dynamic visualization
manipulate::manipulate({
# for (i in seq_len(dim(Y$paths)[3])) {
# pdf(paste0("S2_euler_", i, ".pdf"), width = 7, height = 7)

  sc3 <- scatterplot3d::scatterplot3d(Y$paths[, , 1], color = col_X_alp,
                                      pch = 19, xlim = c(-1, 1),
                                      ylim = c(-1, 1), zlim = c(-1, 1),
                                      xlab = "x", ylab = "y", zlab = "z")
  sc3$points3d(rbind(Y$paths[, , i]), col = col_X, pch = 16, cex = 0.75)

  for (k in seq_len(nrow(Y$paths))) {

    sc3$points3d(t(Y$paths[k, , ]), col = col_X_alp[k], type = "l")

  }
  # Lines of the smoothed ridge
  # sc3$points3d(ridge$ridge_grid, type = "l", lwd = 3)

  # # Points on the ridge joined by MDS
  # sc3$points3d(Y$ridge_y, col = 1, pch = 17,
  #              cex = 0.6)

# dev.off()
# }
}, i = manipulate::slider(1, dim(Y$paths)[3]))

# Visualization on S^2
rgl::open3d()
rgl::par3d(windowRect = c(80, 125, 1280, 826), zoom = 0.78)
viz_euler_s2r(y_old = y, X = X, col_X = col_X_alp, k = 1, ind_dj = ind_dj,
              d = d, h = h_rid, h_euler = h_eu, N = N, eps = eps, frame = TRUE,
              grad = FALSE, rays = FALSE, box = FALSE, axes = FALSE)
rgl::lines3d(t(Y$paths[n + 1, , ]), col = 6, lwd = 5)

## (S^2)^2 test

# Sample
set.seed(3242)
r <- 2
d <- rep(2, r)
n <- 200
ind_dj <- comp_ind_dj(d)
X <- r_path_s2r(n = n, r = r[1], spiral = FALSE,
                Theta = cbind(c(1, 0, 0), c(0, 1, 0)), sigma = 0.2)
X <- cbind(X[, , 1], X[, , 2])
col_X_alp <- viridis::viridis(n, alpha = 0.25)
col_X <- viridis::viridis(n)
stopifnot(r == 2)

# Single starting point
i <- 10
y <- X[i, , drop = FALSE]
# y <- rbind(c(-1, 0, 0, 1, 0, 0))
col_X_alp <- c(col_X_alp, 1)
col_X <- c(col_X, 1)

# Euler and kde
h_rid <- rep(0.5, r)
h_eu <- h_rid^2
N <- 200
eps <- 1e-6
Xy <- rbind(X, y)
Y <- euler_ridge(x = Xy, X = X, d = d, h = h_rid, h_euler = h_eu,
                 N = N, eps = eps, keep_paths = TRUE)
# kde_X <- kde_polysph(x = Xy, X = X, d = d, h = h_rid)
# cut_X <- cut(kde_X, breaks = seq(min(kde_X), max(kde_X), l = 21),
#              include.lowest = TRUE)
# col_X_alp <- viridis::viridis(20, alpha = 0.25)[cut_X]
# col_X <- viridis::viridis(20)[cut_X]

# Dynamic visualization
manipulate::manipulate({

  old_par <- par(mfrow = c(1, r))
  for (k in seq_len(r)) {

    ind_k <- 1:3 + 3 * (k - 1)
    sc3 <- scatterplot3d::scatterplot3d(rbind(Y$paths[, ind_k, 1]),
                                        color = col_X_alp, pch = 19,
                                        xlim = c(-1, 1), ylim = c(-1, 1),
                                        zlim = c(-1, 1))
    sc3$points3d(rbind(Y$paths[, ind_k, i]), col = col_X, pch = 16, cex = 0.75)

    for (k in seq_len(nrow(Y$paths))) {

      sc3$points3d(t(Y$paths[k, ind_k, ]), col = col_X_alp[k], type = "l")

    }

  }
  par(old_par)

}, i = manipulate::slider(1, dim(Y$paths)[3]))

# Visualization on (S^2)^2
rgl::open3d()
rgl::par3d(windowRect = c(80, 125, 1280, 826), zoom = 0.78)
viz_euler_s2r(y_old = y, X = X, col_X = col_X_alp, k = 1, ind_dj = ind_dj,
              d = d, h = h_rid, h_euler = h_eu, N = N, eps = eps, frame = FALSE)
rgl::lines3d(t(Y$paths[n + 1, 1:3, ]), col = 6, lwd = 5)
rgl::open3d()
rgl::par3d(windowRect = c(80, 125, 1280, 826), zoom = 0.78)
viz_euler_s2r(y_old = y, X = X, col_X = col_X_alp, k = 2, ind_dj = ind_dj,
              d = d, h = h_rid, h_euler = h_eu, N = N, eps = eps, frame = FALSE)
rgl::lines3d(t(Y$paths[n + 1, 4:6, ]), col = 6, lwd = 5)
