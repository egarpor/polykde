
## Constants kernel

d <- sample(1:4, size = 1, replace = TRUE)

test_that("b_d(L) equals definition for kernel = 1, 2", {
  expect_equal(b_d(d = d, kernel = 1),
               DirStats::b_L(L = function(r) exp(-r), q = d) / d,
               tolerance = 1e-6)
  expect_equal(b_d(d = d, kernel = 2),
               DirStats::b_L(L = function(r) (1 - r) * (r <= 1), q = d) / d,
               tolerance = 1e-6)
})

test_that("v_d(L) equals definition for kernel = 1, 2", {
  expect_equal(v_d(d = d, kernel = 1),
               DirStats::d_L(L = function(r) exp(-r), q = d) /
                 DirStats::lambda_L(L = function(r) exp(-r), q = d),
               tolerance = 1e-6)
  expect_equal(v_d(d = d, kernel = 2),
               DirStats::d_L(L = function(r) (1 - r) * (r <= 1), q = d) /
                 DirStats::lambda_L(L = function(r) (1 - r) * (r <= 1), q = d),
               tolerance = 1e-6)
})

## CV bandwidth selectors

r <- 3
d <- sample(1:4, size = r, replace = TRUE)
h <- rep(0.5, r)
n <- 50
mu <- r_unif_polysph(n = 5, d = d)
X <- r_kde_polysph(n = n, X = mu, d = d, h = h)

test_that("bw_cv_polysph(type = \"LCV\") in sequential and parallel mode", {
  skip("Unstable")
  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                  control = list(maxit = 1e3))$par,
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                  control = list(maxit = 1e3), ncores = 2)$par,
               tolerance = 5e-2)
})

test_that("bw_cv_polysph(type = \"LSCV\") in sequential and parallel mode", {
  skip("Unstable")
  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                  control = list(maxit = 1e3))$par,
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                  control = list(maxit = 1e3), ncores = 2)$par,
    tolerance = 5e-2)
})

## CV equivalence with DirStats::bw_dir_*cv()

r <- 1
d <- 3
h <- 0.5
n <- 50
mu <- r_unif_polysph(n = 5, d = d)
X <- r_kde_polysph(n = n, X = mu, d = d, h = h)

test_that("bw_cv_polysph(type = \"LCV\") equals DirStats::bw_dir_lcv()", {
  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                  method = "L-BFGS-B", bw0 = h)$par,
    DirStats::bw_dir_lcv(data = X, optim = TRUE, optim_par = h)$h_opt,
    tolerance = 1e-3)
})

test_that("bw_cv_polysph(type = \"LSCV\") equals DirStats::bw_dir_lscv()", {
  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                  method = "L-BFGS-B", bw0 = 0.25, M = 1e4)$par,
    DirStats::bw_dir_lscv(data = X, optim = TRUE, optim_par = 0.25)$h_opt,
    tolerance = 5e-2)
})

## Curvature matrices for von Mises--Fisher

r <- rpois(1, lambda = 3) + 1
d <- rpois(r, lambda = 2) + 1
mu <- r_unif_polysph(n = 1, d = d)
kappa <- r * runif(r)
ind_dj <- comp_ind_dj(d = d)

# Approximate curvature using weighted inverse sampling
# \int_{S^{d_1, ..., d_r}} t(x)t(x)' dx
#   = \int \int_{S^{d_1, ..., d_r}} [t(x) / f(x)][t(x) / f(x)]' f(x)^2 dx
#   ≈ 1 / M \sum_{i = 1}^M [t(X_i) / f(X_i)][t(X_i) / f(X_i)]' * f(X_i)
# where t(x) = (tr(H_11 f(x)), ..., tr(H_rr f(x)))'.

# t() for vMF
hess_vmf <- function(x) grad_hess_kde_polysph(x = rbind(x), X = mu, d = d,
                                              h = 1 / sqrt(kappa),
                                              norm_grad_hess = TRUE)$hess[1, , ]
t_vmf <- function(x) {
  H <- hess_vmf(x)
  sapply(1:r, function(j) {
    ind_j <- (ind_dj[j] + 1):ind_dj[j + 1]
    sum(diag(H[ind_j, ind_j]))
    })
}

# Integral estimation
N <- 1e4
X_vmf <- r_vmf_polysph(n = N, d = d, mu = mu, kappa = kappa)
f_X_vmf <- drop(kde_polysph(x = X_vmf, X = mu, d = d, h = 1 / sqrt(kappa)))
tt_X_vmf <- t(rbind(apply(X_vmf, 1, function(xi) tcrossprod(t_vmf(x = xi)))))
R_X_vmf <- matrix(colMeans(tt_X_vmf * f_X_vmf), nrow = r, ncol = r)

test_that("curv_vmf_polysph() curvature matrix against its numerical version", {
  expect_equal(R_X_vmf,
               curv_vmf_polysph(kappa = kappa, d = d),
               tolerance = 1e-2)
})

test_that("sum(curv_vmf_polysph()) curvature term against its numerical
          version", {
  expect_equal(sum(R_X_vmf),
               sum(curv_vmf_polysph(kappa = kappa, d = d)),
               tolerance = 5e-2)
})

test_that("curv_vmf_polysph() is coherent with scalar/vector d", {
  expect_equal(curv_vmf_polysph(kappa = kappa, d = rep(d[1], r)),
               curv_vmf_polysph(kappa = kappa, d = d[1]))
})

test_that("curv_vmf_polysph() equivalence with DirStats::R_Psi_mixvmf()", {
  expect_equal(drop(curv_vmf_polysph(kappa = kappa[1], d = d[1])),
               d[1]^2 * DirStats::R_Psi_mixvmf(q = d[1],
                                               mu = rbind(c(rep(0, d[1]), 1)),
                                               kappa = kappa[1], p = 1))
})

## Objective functions for plug-in bandwidth selectors

# Setup
r <- rpois(1, lambda = 5) + 1
d <- rpois(r, lambda = 2) + 1
h <- runif(r)
n <- 50
kappa <- 1:r
mu <- r_unif_polysph(n = 5, d = d)
X <- r_kde_polysph(n = n, X = mu, d = d, h = h)
ind <- cumsum(c(1, d + 1))

# Common objects
R_kappa <- curv_vmf_polysph(kappa = kappa, d = d)
bias2 <- tcrossprod(b_d(kernel = 1, d = d)) * R_kappa
var <- prod(v_d(kernel = 1, d = d)) / n

# AMISE and derivative
log_amise <- function(h) {

  var2 <- var / abs(prod(h^d))
  obj <- sum(tcrossprod(h^2) * bias2) + var2
  attr(obj, "gradient") <- (bias2 %*% (4 * h^2) * h - d * var2 / abs(h)) / obj
  return(log(obj))

}

# Common objects
log_R_kappa <- curv_vmf_polysph(kappa = kappa, d = d, log = TRUE)
log_bias2 <- log(tcrossprod(b_d(kernel = 1, d = d))) + log_R_kappa
log_var <- sum(log(v_d(kernel = 1, d = d))) - log(n)

# log(exp(log_x) + exp(y))
log_sum_exp <- function(x) {

  M <- max(x)
  M + log(sum(exp(x - M)))

}

# AMISE and derivative
log_amise_stable <- function(h) {

  log_h <- log(abs(h))
  log_var2 <- log_var - sum(d * log_h)
  logs <- log(tcrossprod(h^2)) + log_bias2 - log_var2
  log_obj <- log_var2 + log_sum_exp(x = c(logs, 0))
  attr(log_obj, "gradient") <-
    exp(log(4) + log_bias2 - log_obj + log_h) %*% h^2 -
    exp(log(d) + log_var2 - log_obj - log_h)
  return(log_obj)

}
log_amise_stable_log_h <- function(log_h) {

  h <- exp(log_h)
  log_var2 <- log_var - sum(d * log_h)
  logs <- log(tcrossprod(h^2)) + log_bias2 - log_var2
  log_obj <- log_var2 + log_sum_exp(x = c(logs, 0))
  attr(log_obj, "gradient") <-
    (exp(log(4) + log_bias2 - log_obj + log_h) %*% h^2 -
       exp(log(d) + log_var2 - log_obj - log_h)) * h
  return(log_obj)

}

test_that("Coherence between bw_rot_polysph() and bw_mrot_polysph()", {
  expect_equal(
    sapply(seq_along(d), function(j) {
      data <- X[, ind[j]:(ind[j + 1] - 1)]
      bw_rot_polysph(X = data, d = d[j])$bw
      }),
    bw_mrot_polysph(X = X, d = d),
    tolerance = 1e-6)
})

test_that("Coherence between bw_rot_polysph() and DirStats::bw_dir_rot()", {
  expect_equal(
    sapply(seq_along(d), function(j) {
      data <- X[, ind[j]:(ind[j + 1] - 1)]
      bw_rot_polysph(X = data, d = d[j])$bw
    }),
    sapply(seq_along(d), function(j) {
      data <- X[, ind[j]:(ind[j + 1] - 1)]
      DirStats::bw_dir_rot(data = data)
    }), tolerance = 1e-4)
})

test_that("Derivatives of log_amise_*", {
  for (i in 1:5) {
    h <- 2 * runif(r)
    expect_equal(
      drop(attr(log_amise(h), which = "gradient")),
      numDeriv::grad(func = log_amise, x = h, method = "Richardson",
                     method.args = list(eps = 1e-10)),
      tolerance = 1e-7
      )
    expect_equal(
      drop(attr(log_amise_stable(h), which = "gradient")),
      numDeriv::grad(func = log_amise_stable, x = h, method = "Richardson",
                     method.args = list(eps = 1e-10)),
      tolerance = 1e-7
    )
    expect_equal(
      drop(attr(log_amise_stable_log_h(log(h)), which = "gradient")),
      numDeriv::grad(func = log_amise_stable_log_h, x = log(h),
                     method = "Richardson", method.args = list(eps = 1e-10)),
      tolerance = 1e-7
    )
  }
})

test_that("log_amise vs. log_amise_stable vs. log_amise_stable_log_h", {
  for (i in 1:5) {
    h <- 2 * runif(r)
    expect_equal(
      log_amise(h),
      log_amise_stable(h)
    )
    expect_equal(
      log_amise(h),
      {lh <- log_amise_stable_log_h(log(h))
      attr(lh, "gradient") <- attr(lh, "gradient") / h
      lh}
    )
  }
})

test_that("log_amise vs. log_amise_stable", {
  for (i in 1:5) {
    h <- 2 * runif(r)
    expect_equal(
      log_amise(h),
      log_amise_stable(h)
    )
  }
})

test_that("Same optimization with log_amise_stable or log_amise_stable_log_h", {

  h0 <- 2 * runif(r)
  opt1 <- nlm(f = log_amise_stable_log_h, p = log(h0))
  opt2 <- nlm(f = log_amise_stable, p = exp(opt1$estimate))
  opt3 <- nlm(f = log_amise_stable, p = h0)
  expect_true(
    (abs(opt1$minimum - opt2$minimum) < 1e-10) ||
    (max(abs(exp(opt1$est) - abs(opt2$est))) < 1e-10)
  )
  expect_lt(opt1$minimum - 1e-10, opt3$minimum)

})

## Plug-in bandwidth selectors

test_that("Same result for vMF kernel with kernel_type = 1,2", {

  # Setup
  r <- rpois(1, lambda = 5) + 1
  d <- rpois(r, lambda = 2) + 1
  h <- runif(r)
  n <- 50
  kappa <- 1:r
  mu <- r_unif_polysph(n = 5, d = d)
  X <- r_kde_polysph(n = n, X = mu, d = d, h = h)
  expect_equal(
    bw_rot_polysph(X = X, d = d, kernel = 1, kernel_type = 1)$bw,
    bw_rot_polysph(X = X, d = d, kernel = 1, kernel_type = 2)$bw
  )

})

test_that("Same result for Epa and sfp kernel with kernel_type = 1,2
          and r = 1", {

  # Setup
  d <- rpois(1, lambda = 2) + 1
  mu <- r_unif_polysph(n = 1, d = d)
  X <- r_kde_polysph(n = 50, X = mu, d = d, h = runif(1))
  expect_equal(
    bw_rot_polysph(X = X, d = d, kernel = 2, kernel_type = 1)$bw,
    bw_rot_polysph(X = X, d = d, kernel = 2, kernel_type = 2)$bw
  )
  expect_equal(
    bw_rot_polysph(X = X, d = d, kernel = 3, kernel_type = 1)$bw,
    bw_rot_polysph(X = X, d = d, kernel = 3, kernel_type = 2)$bw
  )

})

test_that("Same result with kappa precomputed", {

    # Setup
    r <- 2
    d <- rpois(r, lambda = 2) + 1
    mu <- r_unif_polysph(n = 1, d = d)
    X <- r_kde_polysph(n = 50, X = mu, d = d, h = runif(r))
    ind <- cumsum(c(1, d + 1))
    kappa <- sapply(seq_along(d), function(j) {

      # Prepare data + fit vMF
      data <- X[, ind[j]:(ind[j + 1] - 1)]
      min(DirStats::norm2(movMF::movMF(x = data, k = 1, type = "S",
                                       maxit = 300)$theta),
          5e4) # Breaking point for later Bessels

    })
    expect_equal(
      bw_rot_polysph(X = X + 1, d = d, kappa = kappa)$bw,
      bw_rot_polysph(X = X, d = d, kappa = NULL)$bw
    )
    expect_equal(
      bw_mrot_polysph(X = X + 1, d = d, kappa = kappa),
      bw_mrot_polysph(X = X, d = d, kappa = NULL)
    )

})
