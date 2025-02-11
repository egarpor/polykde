
# Randomize testing
r <- 3
d <- c(1, 2, 3)
h <- runif(r, 0.2, 1.5)
n <- 50
p1 <- runif(1)
p2 <- 1 - p1
n1 <- round(0.3 * n)
n2 <- n - n1
nx <- 1e2
x <- r_unif_polysph(n = nx, d = d)
X <- r_unif_polysph(n = n, d = d)
X1 <- X[1:n1, ]
X2 <- X[(n1 + 1):n, ]
x_int <- rbind(X, r_unif_polysph(n = 1e5, d = d))
# Fine for integrating up to r = 3

## Basics

test_that("Normalization", {
  expect_equal(
    max(abs(kde_polysph(x = x, X = 2 * X, d = d, h = h,
                        norm_X = TRUE, norm_x = FALSE) -
              kde_polysph(x = x, X = X, d = d, h = h))), 0)
  expect_equal(
    max(abs(kde_polysph(x = 2 * x, X = X, d = d, h = h,
                        norm_X = FALSE, norm_x = TRUE) -
              kde_polysph(x = x, X = X, d = d, h = h))), 0)
  expect_equal(
    max(abs(kde_polysph(x = 2 * x, X = 3 * X, d = d, h = h,
                        norm_X = TRUE, norm_x = TRUE) -
              kde_polysph(x, X, d, h))), 0)
})

test_that("Log-kde", {
  expect_equal(
    max(abs(exp(kde_polysph(x = x, X = X, d = d, h = h, log = TRUE)) -
              kde_polysph(x = x, X = X, d = d, h = h))), 0)
})

test_that("Kernel normalization", {
  for (ki in 1:3) {
    expect_equal(max(abs(
      prod(c_kern(h = h, d = d, kernel = ki, inc_sfp = FALSE)) *
        kde_polysph(x = x, X = X, d = d, h = h, kernel = ki,
                    normalized = FALSE) -
        kde_polysph(x = x, X = X, d = d, h = h, kernel = ki,
                    normalized = TRUE))), 0)
  }
})

## Integration

test_that("Integration for vMF kernel", {
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h, kernel = 1)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h,
                                wrt_unif = TRUE, kernel = 1)),
               1, tolerance = 1e-1)
})

test_that("Integration for Epa kernel", {
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h, kernel = 2)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h, wrt_unif = TRUE,
                                kernel = 2)),
               1, tolerance = 1e-1)
})

test_that("Integration for sfp kernel", {
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h, kernel = 3,
                                  k = 5)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h, wrt_unif = TRUE,
                                kernel = 3, k = 5)),
               1, tolerance = 1e-1)
})

test_that("Integration for intrinsic vMF kernel", {
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h, kernel = 1,
                                  intrinsic = TRUE)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h, wrt_unif = TRUE,
                                kernel = 1, intrinsic = TRUE)),
               1, tolerance = 1e-1)
})

test_that("Integration for intrinsic Epa kernel", {
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h, kernel = 2,
                                  intrinsic = TRUE)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h,
                                wrt_unif = TRUE, kernel = 2, intrinsic = TRUE)),
               1, tolerance = 1e-1)
})

test_that("Integration for intrinsic sfp kernel", {
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h, kernel = 3,
                                  k = 10, intrinsic = TRUE)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h, wrt_unif = TRUE,
                                kernel = 3, k = 10, intrinsic = TRUE)),
               1, tolerance = 1e-1)
})

## Softplus gives Epa as limiting case

test_that("Convergence of sfp to Epa", {
  epa <- kde_polysph(x = x, X = X, d = d, h = h, kernel = 2)
  sfp_1 <- kde_polysph(x = x, X = X, d = d, h = h, kernel = 3, k = 1)
  sfp_5 <- kde_polysph(x = x, X = X, d = d, h = h, kernel = 3, k = 5)
  sfp_10 <- kde_polysph(x = x, X = X, d = d, h = h, kernel = 3, k = 10)
  sfp_100 <- kde_polysph(x = x, X = X, d = d, h = h, kernel = 3, k = 100)
  expect_gt(max(abs(epa - sfp_1)), max(abs(epa - sfp_5)))
  expect_gt(max(abs(epa - sfp_5)), max(abs(epa - sfp_10)))
  expect_gt(max(abs(epa - sfp_10)), max(abs(epa - sfp_100)))
  expect_equal(epa, sfp_100, tolerance = 1e-3)
})

## Product and spherically symmetric kernels

test_that("Product and spherically symmetric vMF kernels coincide", {
  expect_equal(kde_polysph(x = x, X = X, d = d, h = h, kernel = 1,
                           kernel_type = 1),
               kde_polysph(x = x, X = X, d = d, h = h, kernel = 1,
                           kernel_type = 2))
})

test_that("Unnormalized product and spherically symmetric kernels coincide
          for r = 1", {
  for (kernel in 1:3) {
    r <- 1
    d <- sample(1:4, size = 1, replace = TRUE)
    h <- runif(r, 0.2, 1.5)
    n <- 5e1
    nx <- 1e2
    x <- r_unif_polysph(n = nx, d = d)
    X <- r_unif_polysph(n = n, d = d)
    expect_equal(kde_polysph(x = x, X = X, d = d, h = h, kernel = kernel,
                             kernel_type = 1, normalized = FALSE),
                 kde_polysph(x = x, X = X, d = d, h = h, kernel = kernel,
                             kernel_type = 2, normalized = FALSE))
  }
})

test_that("Integration for spherically symmetric Epa kernel", {
  skip("Unstable")
  h_small <- rep(0.5, r)
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h_small,
                                  kernel = 2, kernel_type = 2)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h_small,
                                wrt_unif = TRUE, kernel = 2,
                                kernel_type = 2)),
               1, tolerance = 1e-1)
})

test_that("Integration for spherically symmetric sfp kernel", {
  h_small <- rep(0.5, r)
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h_small,
                                  kernel = 3, k = 10, kernel_type = 2)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h_small,
                                wrt_unif = TRUE, kernel = 3, k = 10,
                                kernel_type = 2)),
               1, tolerance = 1e-1)
})

test_that("Integration for intrinsic spherically symmetric Epa
          kernel", {
  h_small <- rep(0.5, r)
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h_small,
                                  intrinsic = TRUE, kernel = 2,
                                  kernel_type = 2)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h_small,
                                wrt_unif = TRUE, intrinsic = TRUE, kernel = 2,
                                kernel_type = 2)),
               1, tolerance = 1e-1)
})

test_that("Integration for intrinsic spherically symmetric sfp kernel", {
  h_small <- rep(0.5, r)
  expect_equal(prod(rotasym::w_p(p = d + 1)) *
                 mean(kde_polysph(x = x_int, X = X, d = d, h = h_small,
                                  intrinsic = TRUE, kernel = 3,
                                  kernel_type = 2)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x_int, X = X, d = d, h = h_small,
                                intrinsic = TRUE, wrt_unif = TRUE, kernel = 3,
                                kernel_type = 2)),
               1, tolerance = 1e-1)
})

## Log-cv functions

test_that("Log-cv kde vMF kernel", {
  for (kernel in 1:3) {
    expect_equal(
      log_cv_kde_polysph(X = X, d = d, h = h, wrt_unif = TRUE,
                         kernel = kernel)[1],
      drop(log(kde_polysph(x = X[1, , drop = FALSE], X = X[-1, , drop = FALSE],
                           d = d, h = h, wrt_unif = TRUE, kernel = kernel))))
  }
})

test_that("Log-cv kde Epa kernel", {
  for (kernel in 1:3) {
    expect_equal(
      log_cv_kde_polysph(X = X, d = d, h = h, wrt_unif = TRUE,
                         kernel = kernel)[1],
      drop(log(kde_polysph(x = X[1, , drop = FALSE], X = X[-1, , drop = FALSE],
                           d = d, h = h, wrt_unif = TRUE, kernel = kernel))))
  }
})

test_that("Log-cv kde with weights", {
  expect_equal(
    log_cv_kde_polysph(X = X, d = d, h = h),
    log_cv_kde_polysph(X = X, d = d, h = h, weights = rep(1, n)))
})

test_that("Log-cv kde with weights is compatible with kde", {
  ws <- n:1
  for (kernel in 1:3) {
    expect_equal(
      drop(log_cv_kde_polysph(X = X, d = d, h = h, wrt_unif = TRUE,
                              weights = ws, kernel = kernel)),
      sapply(1:n, function(i) {
        drop(kde_polysph(x = X[i, , drop = FALSE], X = X[-i, , drop = FALSE],
                         d = d, h = h, wrt_unif = TRUE, weights = ws[-i],
                         log = TRUE, kernel = kernel))
        }))
  }
})

## Weights

test_that("Internal normalization of weights, for both kernels", {
  for (kernel in 1:3) {
    expect_equal(kde_polysph(x = x, X = X, d = d, h = h, wrt_unif = TRUE,
                             kernel = kernel),
                 kde_polysph(x = x, X = X, d = d, h = h, wrt_unif = TRUE,
                             weights = rep(2, n), kernel = kernel))
  }
})

test_that("Mixtures computation with weights, for both kernels", {
  ws <- c(rep(p1 / n1, n1), rep(p2 / n2, n2))
  for (kernel in 1:3) {
    expect_equal(p1 * kde_polysph(x = x, X = X1, d = d, h = h,
                                 wrt_unif = TRUE, kernel = kernel) +
                 p2 * kde_polysph(x, X2, d, h, wrt_unif = TRUE,
                                  kernel = kernel),
               kde_polysph(x = x, X = X, d = d, h = h, weights = ws,
                           wrt_unif = TRUE, kernel = kernel))
  }
})

test_that("Mixtures computation with weights in log-scale, for both kernels", {
  ws <- c(rep(p1 / n1, n1), rep(p2 / n2, n2))
  for (kernel in 1:3) {
    expect_equal(log(p1 * kde_polysph(x = x, X = X1, d = d, h = h,
                                    wrt_unif = TRUE, kernel = kernel) +
                     p2 * kde_polysph(x, X2, d, h, wrt_unif = TRUE,
                                      kernel = kernel)),
                 kde_polysph(x = x, X = X, d = d, h = h, weights = ws,
                             wrt_unif = TRUE, kernel = kernel, log = TRUE))
  }
})

## Bandwidths

h_min <- rep(bw_lcv_min_epa(X = X, d = d), r) + 1e-10
test_that("bw_lcv_min_epa is good enough", {
  pert_1 <- runif(r, 0, 1e-2)
  pert_2 <- runif(r, 0, 1e-3)
  h_good_1 <- h_min + pert_1
  h_good_2 <- h_min + pert_2
  expect_true(all(is.finite(log_cv_kde_polysph(X = X, d = d, h = h_min,
                                               wrt_unif = TRUE, kernel = 2))))
  expect_true(all(is.finite(log_cv_kde_polysph(X = X, d = d, h = h_good_1,
                                               wrt_unif = TRUE, kernel = 2))))
  expect_true(all(is.finite(log_cv_kde_polysph(X = X, d = d, h = h_good_2,
                                               wrt_unif = TRUE, kernel = 2))))
})
test_that("bw_lcv_min_epa is the critical point", {
  pert_1 <- runif(r, 0, 1e-2)
  pert_2 <- runif(r, 0, 1e-3)
  h_bad_1 <- h_min - pert_1
  h_bad_2 <- h_min - pert_2
  expect_true(any(!is.finite(log_cv_kde_polysph(X = X, d = d, h = h_bad_1,
                                                wrt_unif = TRUE, kernel = 2))))
  expect_true(any(!is.finite(log_cv_kde_polysph(X = X, d = d, h = h_bad_2,
                                                wrt_unif = TRUE, kernel = 2))))
})

## Others

test_that("Equivalence with DirStats::kde_dir", {
  r <- 1
  d <- rpois(r, 2) + 1
  h <- runif(r, 0.5, 1.5)
  x <- r_unif_polysph(n = 1e3, d = d)
  X <- r_unif_polysph(n = 1e2, d = d)
  expect_equal(
    drop(kde_polysph(x = x, X = X, d = d, h = h)),
    DirStats::kde_dir(x = x, data = X, h = h)
  )
})

## d = 2 for vMF kernel

test_that("Integration for vMF kernel in d = 2", {
  x <- r_unif_polysph(n = 1e4, d = 2)
  X <- r_unif_polysph(n = 10, d = 2)
  h <- runif(1, 0.1, 1)
  expect_equal(prod(rotasym::w_p(p = 2 + 1)) *
                 mean(kde_polysph(x = x, X = X, d = 2, h = h, kernel = 1)),
               1, tolerance = 1e-1)
  expect_equal(mean(kde_polysph(x = x, X = X, d = 2, h = h,
                                wrt_unif = TRUE, kernel = 1)), 1,
               tolerance = 1e-1)
})
