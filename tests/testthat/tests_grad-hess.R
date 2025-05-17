
# Randomize testing

r <- sample(2:3, size = 1)
d <- rpois(r, lambda = 2) + 1
n <- 100
ind_dj <- comp_ind_dj(d = d)
X <- r_unif_polysph(n = n, d = d)
x <- X[sample(n, size = 1), , drop = FALSE]
h <- 2 / (1:r)

# Analytical functions

grad_bar <- function(f, x, d) {

  # Gradient
  grad <- numDeriv::grad(func = f, x = x, method.args = list(eps = 1e-12))

  # Projections
  ind <- cumsum(c(1, d + 1))
  for (j in seq_along(d)) {

    ind_j <- ind[j]:(ind[j + 1] - 1)
    I_j <- diag(1, nrow = d[j] + 1, ncol = d[j] + 1)
    grad[ind_j] <- grad[ind_j] %*% (I_j - tcrossprod(x[ind_j]))

  }
  return(grad)

}

grad_kde_vmf <- function(x, X, d, h) {

  h_tilde <- rep(h, times = d + 1)
  gr <- numeric(ncol(X))
  for (i in seq_len(nrow(X))) {

    gr <- gr + drop(kde_polysph(x = x, X = X[i, , drop = FALSE],
                                d = d, h = h)) * X[i, ]

  }
  gr <- gr / (n * h_tilde^2)
  return(gr)

}

hess_bar <- function(f, x, d) {

  # Full gradient and Hessian
  grad <- numDeriv::grad(func = f, x = x, method.args = list(eps = 1e-12))
  hess <- numDeriv::hessian(func = f, x = x, method.args = list(eps = 1e-12))
  hess_bar <- matrix(NA, nrow = nrow(hess), ncol = ncol(hess))

  # Projections
  ind <- cumsum(c(1, d + 1))
  for (j in seq_along(d)) {

    ind_j <- ind[j]:(ind[j + 1] - 1)
    I_jj <- diag(1, nrow = d[j] + 1, ncol = d[j] + 1)

    for (k in seq_along(d)) {

      ind_k <- ind[k]:(ind[k + 1] - 1)
      I_kk <- diag(1, nrow = d[k] + 1, ncol = d[k] + 1)

      if (j == k) {

        hess_bar[ind_j, ind_j] <- (I_jj - tcrossprod(x[ind_j])) %*%
          hess[ind_j, ind_j] %*% (I_jj - tcrossprod(x[ind_j])) -
          s(x[ind_j] %*% t(grad[ind_j]), add = TRUE) -
          (I_jj - 3 * tcrossprod(x[ind_j])) * drop(grad[ind_j] %*% x[ind_j])

      } else {

        hess_bar[ind_k, ind_j] <- (I_kk - x[ind_k] %*% t(x[ind_k])) %*%
          hess[ind_k, ind_j] %*% (I_jj - x[ind_j] %*% t(x[ind_j]))

      }

    }
  }

  return(hess_bar)

}

hess_kde_vmf <- function(x, X, d, h) {

  ind_dj <- comp_ind_dj(d)
  X_diam <- diamond_crossprod(X, ind_dj)

  h_tilde <- rep(h, times = d + 1)
  h_tilde <- tcrossprod(h_tilde)
  he <- numeric(ncol(X))
  for (i in seq_len(nrow(X))) {

    he <- he + drop(kde_polysph(x = x, X = X[i, , drop = FALSE],
                                d = d, h = h)) * X_diam[i, , ]

  }
  he <- he / (n * h_tilde^2)
  return(he)

}

hezz <- function(f, x, d) {

  # Full gradient and Hessian
  grad <- numDeriv::grad(func = f, x = x, method.args = list(eps = 1e-12))
  hess <- numDeriv::hessian(func = f, x = x, method.args = list(eps = 1e-12))

  # Projections
  ind <- cumsum(c(1, d + 1))
  A <- matrix(0, nrow = ind[length(d) + 1] - 1, ncol = ind[length(d) + 1] - 1)
  for (j in seq_along(d)) {

    ind_j <- ind[j]:(ind[j + 1] - 1)
    I_jj <- diag(1, nrow = d[j] + 1, ncol = d[j] + 1)
    A[ind_j, ind_j] <- drop(x[ind_j] %*% grad[ind_j]) * I_jj

  }

  # Hessian in (10) in https://arxiv.org/pdf/2110.08505.pdf
  P <- proj_P(x = x, d = d)
  P %*% (hess - A) %*% P

}

proj_grad_kde_vmf <- function(x, X, d, h) {

  kde <- kde_polysph(x = x, X = X, d = d, h = h)
  grad <- grad_bar(f = function(y) kde_polysph(x = y, X = X, d = d, h = h,
                                               norm_x = FALSE),
                   x = x, d = d)
  hess <- hess_bar(f = function(y) kde_polysph(x = y, X = X, d = d, h = h,
                                               norm_x = FALSE),
                   x = x, d = d)
  eig <- eigen(hess, symmetric = TRUE)
  eta <- drop(grad %*% tcrossprod(eig$vectors[, -1])) / drop(kde)
  return(eta)

}

proj_P <- function(x, d) {

  # Projections
  ind <- cumsum(c(1, d + 1))
  P <- matrix(0, nrow = ind[length(d) + 1] - 1, ncol = ind[length(d) + 1] - 1)
  for (j in seq_along(d)) {

    ind_j <- ind[j]:(ind[j + 1] - 1)
    P[ind_j, ind_j] <- -tcrossprod(x[ind_j])

  }
  diag(P) <- diag(P) + 1
  return(P)

}

## Gradient

grad_ana <- grad_bar(f = function(y) kde_polysph(x = y, X = X, d = d, h = h,
                                                 norm_x = FALSE),
                     x = x, d = d)

grad_num <- numDeriv::grad(func = function(x) {
  kde_polysph(x = x, X = X, d = d, h = h, norm_x = TRUE)
}, x = x)

grad_num_unproj <- numDeriv::grad(func = function(y)
  kde_polysph(x = y, X = X, d = d, h = h, norm_x = FALSE), x = x)

test_that("Numerical and analytical projected gradients coincide", {
  expect_equal(grad_ana, grad_num, tolerance = 1e-9)
})

test_that("Unprojected analytical and vMF-specific R gradient coincide", {
  expect_equal(grad_kde_vmf(x = x, X = X, d = d, h = h),
               grad_num_unproj, tolerance = 1e-9)
})

test_that("Unprojected analytical and vMF-specific Rcpp gradient coincide", {
  gh_1 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = FALSE,
                                norm_grad_hess = FALSE)
  gh_2 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = FALSE,
                                norm_grad_hess = TRUE)
  expect_equal(drop(gh_1$grad), grad_num_unproj, tolerance = 1e-9)
  expect_equal(drop(gh_2$grad), drop(gh_1$grad) / drop(gh_1$dens),
               tolerance = 1e-9)
})

test_that("Analytical and vMF-specific Rcpp gradient coincide", {
  gh_1 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE,
                                norm_grad_hess = FALSE)
  gh_2 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE,
                                norm_grad_hess = TRUE)
  expect_equal(drop(gh_1$grad), grad_ana, tolerance = 1e-9)
  expect_equal(drop(gh_2$grad), drop(gh_1$grad) / drop(gh_1$dens),
               tolerance = 1e-9)
})

test_that("Gradients vanish in the x direction", {
  expect_equal(drop(x %*% grad_ana), 0)
  expect_equal(drop(x %*% grad_num), 0)
  expect_equal(drop(x %*% drop(
  grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE)$grad)), 0)
})

## Hessian

hess_ana <- hess_bar(f = function(y) kde_polysph(x = y, X = X, d = d, h = h,
                                                 norm_x = FALSE),
                     x = x, d = d)

hess_zz <- hezz(f = function(y) kde_polysph(x = y, X = X, d = d, h = h,
                                            norm_x = FALSE),
                x = x, d = d)

P <- proj_P(x = x, d = d)
test_that("Hezzian is the projection of Hessian", {
  expect_equal(P %*% hess_ana %*% P, hess_zz, tolerance = 1e-6)
  expect_equal(grad_hess_kde_polysph(x = x, X = X, d = d, h = h,
                                     projected = TRUE,
                                     proj_alt = TRUE)$hess[1, , ],
               hess_zz, tolerance = 1e-6)
})

hess_num <- numDeriv::hessian(func = function(x) {
  kde_polysph(x = x, X = X, d = d, h = h, norm_x = TRUE)
}, x = x, method.args = list(eps = 1e-12))

hess_num_unproj <- numDeriv::hessian(func = function(y)
  kde_polysph(x = y, X = X, d = d, h = h, norm_x = FALSE), x = x,
  method.args = list(eps = 1e-12))

test_that("Numerical and analytical projected Hessian coincide", {
  expect_equal(hess_ana, hess_num, tolerance = 1e-6)
})

test_that("Unprojected analytical and vMF-specific R Hessian coincide", {
  expect_equal(hess_kde_vmf(x = x, X = X, d = d, h = h),
                      hess_num_unproj, tolerance = 1e-7)
})

test_that("Unprojected analytical and vMF-specific Rcpp Hessian coincide", {
  gh_1 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = FALSE,
                                norm_grad_hess = FALSE)
  gh_2 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = FALSE,
                                norm_grad_hess = TRUE)
  expect_equal(drop(gh_1$hess[1, , ]), hess_num_unproj, tolerance = 1e-7)
  expect_equal(drop(gh_2$hess[1, , ]), drop(gh_1$hess[1, , ]) / drop(gh_1$dens),
               tolerance = 1e-7)
})

test_that("Analytical and vMF-specific Rcpp Hessian coincide", {
  gh_1 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE,
                                norm_grad_hess = FALSE, proj_alt = FALSE)
  gh_2 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE,
                                norm_grad_hess = TRUE, proj_alt = FALSE)
  expect_equal(drop(gh_1$hess[1, , ]), hess_ana, tolerance = 1e-7)
  expect_equal(drop(gh_2$hess[1, , ]), drop(gh_1$hess[1, , ]) / drop(gh_1$dens),
               tolerance = 1e-7)
})

test_that("Hessian vanish in the x direction", {
  expect_equal(drop(x %*% hess_ana %*% t(x)), 0)
  expect_equal(drop(x %*% hess_num %*% t(x)), 0)
  expect_equal(drop(x %*% hess_zz %*% t(x)), 0)
  expect_equal(drop(x %*% grad_hess_kde_polysph(x = x, X = X, d = d, h = h,
                                                projected = TRUE)$hess[1, , ]
                    %*% t(x)), 0)
})

## Projected gradient

test_that("Analytical and vMF-specific Rcpp projected gradient coincide", {
  expect_equal(proj_grad_kde_vmf(x = x, X = X, d = d, h = h),
               expect_warning(drop(proj_grad_kde_polysph(
                 x = x, X = X, d = d, h = h, proj_alt = FALSE,
                 fix_u1 = FALSE)$eta)),
               tolerance = 1e-6)
})

test_that("vMF-specific Rcpp projected gradient with/without sparsity", {
  expect_equal(proj_grad_kde_polysph(x = x, X = X, d = d, h = h,
                                     sparse = TRUE, fix_u1 = FALSE)$eta,
               proj_grad_kde_polysph(x = x, X = X, d = d, h = h,
                                     sparse = FALSE, fix_u1 = FALSE)$eta,
               tolerance = 1e-6)
  expect_equal(proj_grad_kde_polysph(x = x, X = X, d = d, h = h,
                                     sparse = TRUE)$eta,
               proj_grad_kde_polysph(x = x, X = X, d = d, h = h,
                                     sparse = FALSE)$eta,
               tolerance = 1e-6)
})

test_that("Projected gradient is orthogonal to x", {
  eta_1 <- drop(proj_grad_kde_polysph(x = x, X = X, d = d, h = h,
                                      proj_alt = TRUE)$eta)
  eta_2 <- suppressWarnings(drop(
    proj_grad_kde_polysph(x = x, X = X, d = d, h = h, proj_alt = FALSE)$eta))
  P <- AP(x = x, v = x, ind_dj = ind_dj)$P
  expect_equal(drop(x %*% eta_1), 0, tolerance = 1e-9)
  expect_equal(drop(P %*% eta_1), eta_1, tolerance = 1e-9)
  expect_false(drop(x %*% eta_2) < 1e-9)
  expect_false(max(abs(drop(P %*% eta_2) - eta_2)) < 1e-9)
})

## Normalizing constants

test_that("Normalizing constants in unprojected gradient and Hessians", {
  gh_1 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = FALSE,
                                norm_grad_hess = TRUE, normalized = TRUE)
  gh_2 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = FALSE,
                                norm_grad_hess = TRUE, normalized = FALSE)
  gh_3 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = FALSE,
                                norm_grad_hess = FALSE, normalized = TRUE)
  gh_4 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = FALSE,
                                norm_grad_hess = FALSE, normalized = FALSE)
  cte <- prod(c_kern(d = d, h = h))
  gh_2$dens <- cte * gh_2$dens
  gh_4$dens <- cte * gh_4$dens
  gh_4$grad <- cte * gh_4$grad
  gh_4$hess <- cte * gh_4$hess
  expect_equal(gh_1, gh_2)
  expect_equal(gh_3, gh_4)
})

test_that("Normalizing constants in projected gradient and Hessians", {
  gh_1 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE,
                                norm_grad_hess = TRUE, normalized = TRUE)
  gh_2 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE,
                                norm_grad_hess = TRUE, normalized = FALSE)
  gh_3 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE,
                                norm_grad_hess = FALSE, normalized = TRUE)
  gh_4 <- grad_hess_kde_polysph(x = x, X = X, d = d, h = h, projected = TRUE,
                                norm_grad_hess = FALSE, normalized = FALSE)
  cte <- prod(c_kern(d = d, h = h))
  gh_2$dens <- cte * gh_2$dens
  gh_4$dens <- cte * gh_4$dens
  gh_4$grad <- cte * gh_4$grad
  gh_4$hess <- cte * gh_4$hess
  expect_equal(gh_1, gh_2)
  expect_equal(gh_3, gh_4)
})
