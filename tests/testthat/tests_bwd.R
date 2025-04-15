
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

seed <- 30
set.seed(seed, kind = "Mersenne-Twister")
r <- 2
d <- sample(1:3, size = r, replace = TRUE)
h <- sample(c(0.25, 0.5, 0.75), size = r, replace = TRUE)
n <- 20
mu <- r_unif_polysph(n = 5, d = d)
X <- r_kde_polysph(n = n, X = mu, d = d, h = h)

# CV helper functions
M <- 1e4
seed <- 30
set.seed(seed, kind = "Mersenne-Twister")
mc_samp <- r_unif_polysph(n = M, d = d)
cv_naive <- function(h, X, d, mc_samp, kde_samp = FALSE) {

  if (kde_samp) {

    set.seed(seed, kind = "Mersenne-Twister")
    mc_kde_samp <- r_kde_polysph(n = M, X = X, d = d, h = h, kernel = 1)
    cv_1 <- mean(kde_polysph(x = mc_kde_samp, X = X, d = d, kernel = 1, h = h,
                             wrt_unif = FALSE))
    cv_2 <- 2 * mean(exp(log_cv_kde_polysph(X = X, d = d, kernel = 1, h = h,
                                            wrt_unif = FALSE)))

  } else {

    cv_1 <- prod(rotasym::w_p(p = d + 1)) *
      mean(kde_polysph(x = mc_samp, X = X, d = d, kernel = 1, h = h,
                       wrt_unif = FALSE)^2)
    cv_2 <- 2 * mean(exp(log_cv_kde_polysph(X = X, d = d, kernel = 1, h = h,
                                            wrt_unif = FALSE)))

  }

  return(cv_1 - cv_2)

}
# cv_naive_curve <- function(h1, kde_samp = FALSE) sapply(h1, function(hh1) {
#   cv_naive(h = rep(hh1, r), X = X, d = d, mc_samp = mc_samp,
#            kde_samp = kde_samp)
# })
# cv_bw_cv_curve <- function(h1) sapply(h1, function(hh1) {
#   bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
#                 bw0 = rep(hh1, r), M = M, control = list(maxit = 0),
#                 method = "BFGS", exact_vmf = FALSE, seed_mc = seed)$opt$value
# })
# cv_bw_cv_vmf_curve <- function(h1) sapply(h1, function(hh1) {
#   bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
#                 bw0 = rep(hh1, r), control = list(maxit = 0),
#                 method = "BFGS", exact_vmf = TRUE)$opt$value
# })
#
# # Visualization of LSCV functions
# curve(cv_naive_curve(x, kde_samp = FALSE), from = 0.2, to = 1, n = 100,
#       ylab = "CV loss")
# curve(cv_naive_curve(x, kde_samp = TRUE), from = 0.2, to = 1, n = 100,
#       add = TRUE, lty = 2)
# curve(cv_bw_cv_curve, from = 0.2, to = 1, n = 100, add = TRUE, col = 2)
# curve(cv_bw_cv_vmf_curve, from = 0.2, to = 1, n = 100, add = TRUE, col = 3)

test_that("bw_cv_polysph(type = \"LCV\") in sequential and parallel mode", {
  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                  control = list(maxit = 1e3))$opt$value,
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                  control = list(maxit = 1e3), ncores = 2)$opt$value,
               tolerance = 1e-2)
})

test_that("bw_cv_polysph(type = \"LSCV\") in sequential and parallel mode", {
  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                  control = list(maxit = 1e3), seed_mc = 1)$opt$value,
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                  control = list(maxit = 1e3), seed_mc = 1,
                  ncores = 2)$opt$value,
    tolerance = 1e-2)
})

test_that("bw_cv_polysph(type = \"LSCV\", imp_mc = TRUE) loss", {

  for (f in c(0.25, 0.5, 1, 2)) {

    expect_equal(
      cv_naive(h = f * h, X = X, d = d, kde_samp = TRUE),
      bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                    bw0 = f * h, M = M, control = list(maxit = 0),
                    method = "BFGS", exact_vmf = FALSE, imp_mc = TRUE,
                    seed_mc = seed)$opt$value)

  }

})

test_that("bw_cv_polysph(type = \"LSCV\", imp_mc = FALSE) loss", {

  for (f in c(0.25, 0.5, 1, 2)) {

    expect_equal(
      cv_naive(h = f * h, X = X, d = d, mc_samp = mc_samp, kde_samp = FALSE),
      bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                    bw0 = f * h, M = M, control = list(maxit = 0),
                    method = "BFGS", exact_vmf = FALSE, imp_mc = FALSE,
                    seed_mc = seed)$opt$value,
      tolerance = 5e-2)

  }

})

test_that("bw_cv_polysph(type = \"LSCV\", exact_vmf = TRUE) loss", {

  for (f in c(0.25, 0.5, 1, 2)) {

    expect_equal(
      bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                    bw0 = f * h, M = M, control = list(maxit = 0),
                    method = "BFGS", exact_vmf = FALSE,
                    seed_mc = seed)$opt$value,
      bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                    bw0 = f * h, M = M, control = list(maxit = 0),
                    method = "BFGS", exact_vmf = TRUE)$opt$value,
      tolerance = 5e-2)

  }

})

test_that("bw_cv_polysph(type = \"LCV\", common_h = TRUE)", {

  for (f in c(0.25, 0.5, 1, 2)) {

    expect_equal(
      bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                    bw0 = f * h, M = M, control = list(maxit = 0),
                    method = "BFGS", common_h = TRUE)$opt$value,
      bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                    bw0 = f * rep(mean(h), r), M = M, control = list(maxit = 0),
                    method = "BFGS", common_h = FALSE)$opt$value,
      tolerance = 5e-2)

  }

})

test_that("bw_cv_polysph(type = \"LSCV\", common_h = TRUE)", {

  for (f in c(0.25, 0.5, 1, 2)) {

    expect_equal(
      bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                    bw0 = f * h, M = M, control = list(maxit = 0),
                    method = "BFGS", exact_vmf = TRUE, common_h = TRUE,
                    seed_mc = seed)$opt$value,
      bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                    bw0 = f * rep(mean(h), r), M = M, control = list(maxit = 0),
                    method = "BFGS", exact_vmf = TRUE,
                    common_h = FALSE)$opt$value,
      tolerance = 5e-2)

  }

})

test_that("bw_cv_polysph() with bw0 vector and bw0 matrix", {

  bw0_vec_1 <- h
  bw0_vec_2 <- 0.1 * h
  bw0_vec_3 <- 5 * h
  bw0_mat <- rbind(bw0_vec_1, bw0_vec_2, bw0_vec_3)
  expect_true(
    bw_cv_polysph(X = X, d = d, bw0 = bw0_vec_1,
                  control = list(maxit = 0))$opt$value >=
    bw_cv_polysph(X = X, d = d, bw0 = bw0_mat,
                  control = list(maxit = 0))$opt$value)
  expect_true(
    bw_cv_polysph(X = X, d = d, bw0 = bw0_vec_2,
                  control = list(maxit = 0))$opt$value >=
    bw_cv_polysph(X = X, d = d, bw0 = bw0_mat,
                  control = list(maxit = 0))$opt$value)
  expect_true(
    bw_cv_polysph(X = X, d = d, bw0 = bw0_vec_3,
                  control = list(maxit = 0))$opt$value >=
    bw_cv_polysph(X = X, d = d, bw0 = bw0_mat,
                  control = list(maxit = 0))$opt$value)

})

test_that("bw_cv_polysph() with bw0 vector and bw0 matrix with common_h", {

  bw0_vec_1 <- h
  bw0_vec_2 <- 0.1 * h
  bw0_vec_3 <- 5 * h
  bw0_mat <- rbind(bw0_vec_1, bw0_vec_2, bw0_vec_3)
  suppressWarnings({
  expect_true(
    bw_cv_polysph(X = X, d = d, bw0 = bw0_vec_1, control = list(maxit = 0),
                  common_h = TRUE)$opt$value >=
    bw_cv_polysph(X = X, d = d, bw0 = bw0_mat, control = list(maxit = 0),
                  common_h = TRUE)$opt$value)
  expect_true(
    bw_cv_polysph(X = X, d = d, bw0 = bw0_vec_2, control = list(maxit = 0),
                  common_h = TRUE)$opt$value >=
    bw_cv_polysph(X = X, d = d, bw0 = bw0_mat, control = list(maxit = 0),
                  common_h = TRUE)$opt$value)
  expect_true(
    bw_cv_polysph(X = X, d = d, bw0 = bw0_vec_3, control = list(maxit = 0),
                  common_h = TRUE)$opt$value >=
    bw_cv_polysph(X = X, d = d, bw0 = bw0_mat, control = list(maxit = 0),
                  common_h = TRUE)$opt$value)
  })
})

test_that("bw_cv_polysph(type = \"LCV\") with optim/nlm", {

  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV", opt = "nlm")$bw,
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV", opt = "optim")$bw,
    tolerance = 1e-2)
  expect_error(bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                             opt = "nlm", ncores = 2))

})

## CV equivalence with DirStats::bw_dir_*cv()

r <- 1
d <- 2
h <- 0.25
n <- 20
mu <- r_unif_polysph(n = 5, d = d)
X <- r_kde_polysph(n = n, X = mu, d = d, h = h)

test_that("bw_cv_polysph(type = \"LCV\") equals DirStats::bw_dir_lcv()", {
  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LCV",
                  method = "L-BFGS-B", bw0 = h)$bw,
    DirStats::bw_dir_lcv(data = X, optim = TRUE, optim_par = h)$h_opt,
    tolerance = 1e-3)
})

test_that("bw_cv_polysph(type = \"LSCV\") equals DirStats::bw_dir_lscv()", {
  expect_equal(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV",
                  method = "L-BFGS-B", bw0 = 0.5, exact_vmf = TRUE)$bw,
    DirStats::bw_dir_lscv(data = X, optim = TRUE, optim_par = 0.5)$h_opt,
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
  log_obj <- log_var2 + polykde:::log_sum_exp(logs = c(logs, 0))
  attr(log_obj, "gradient") <-
    exp(log(4) + log_bias2 - log_obj + log_h) %*% h^2 -
    exp(log(d) + log_var2 - log_obj - log_h)
  return(log_obj)

}
log_amise_stable_log_h <- function(log_h) {

  h <- exp(log_h)
  log_var2 <- log_var - sum(d * log_h)
  logs <- log(tcrossprod(h^2)) + log_bias2 - log_var2
  log_obj <- log_var2 + polykde:::log_sum_exp(logs = c(logs, 0))
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
      tolerance = 5e-6
      )
    expect_equal(
      drop(attr(log_amise_stable(h), which = "gradient")),
      numDeriv::grad(func = log_amise_stable, x = h, method = "Richardson",
                     method.args = list(eps = 1e-10)),
      tolerance = 5e-6
    )
    expect_equal(
      drop(attr(log_amise_stable_log_h(log(h)), which = "gradient")),
      numDeriv::grad(func = log_amise_stable_log_h, x = log(h),
                     method = "Richardson", method.args = list(eps = 1e-10)),
      tolerance = 5e-6
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
    (abs(opt1$minimum - opt2$minimum) < 1e-9) ||
    (max(abs(exp(opt1$est) - abs(opt2$est))) < 1e-9)
  )
  expect_lt(opt1$minimum - 1e-9, opt3$minimum)

})

test_that("bw_rot_polysph() with bw0 vector and bw0 matrix", {

  bw0_vec_1 <- h
  bw0_vec_2 <- 0.1 * h
  bw0_vec_3 <- 5 * h
  bw0_mat <- rbind(bw0_vec_1, bw0_vec_2, bw0_vec_3)
  expect_true(bw_rot_polysph(X = X, d = d, bw0 = bw0_vec_1,
                             iterlim = 1)$opt$minimum + 0.1 >=
                bw_rot_polysph(X = X, d = d, bw0 = bw0_mat,
                               iterlim = 1)$opt$minimum)
  expect_true(bw_rot_polysph(X = X, d = d, bw0 = bw0_vec_2,
                             iterlim = 1)$opt$minimum + 0.1 >=
                bw_rot_polysph(X = X, d = d, bw0 = bw0_mat,
                               iterlim = 1)$opt$minimum)
  expect_true(bw_rot_polysph(X = X, d = d, bw0 = bw0_vec_3,
                             iterlim = 1)$opt$minimum + 0.1 >=
                bw_rot_polysph(X = X, d = d, bw0 = bw0_mat,
                               iterlim = 1)$opt$minimum)

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
