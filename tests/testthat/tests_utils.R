
## proj_polysph()

r <- 3
d <- c(2, 3, 1)
n <- 2
X <- r_unif_polysph(n = n, d = d)
ind_dj <- comp_ind_dj(d)

test_that("Normalization in dist_polysph()", {
  expect_equal(proj_polysph(x = X, ind_dj = ind_dj), X)
  expect_equal(proj_polysph(x = 1:n * X, ind_dj = ind_dj), X)
})

test_that("proj_polysph() rejects zero-norm blocks", {
  expect_error(
    proj_polysph(x = rbind(c(0, 0, 1, 1)), ind_dj = comp_ind_dj(c(1, 1))),
    "Zero-norm"
  )
  out <- proj_polysph(x = rbind(c(2, 0, 1, 1)), ind_dj = comp_ind_dj(c(1, 1)))
  expect_equal(unname(out[1, 1:2]), c(1, 0))
})

## dist_polysph()

r <- 2
d <- c(2, 3)
n <- 10
X1 <- r_unif_polysph(n = n, d = d)
X2 <- r_unif_polysph(n = n, d = d)
ind_dj <- comp_ind_dj(d)

test_that("Normalization in dist_polysph()", {
  expect_equal(
    dist_polysph(x = X1, y = X2, ind_dj = ind_dj,
                 norm_x = FALSE, norm_y = FALSE),
    dist_polysph(x = X1, y = X2, ind_dj = ind_dj,
                 norm_x = TRUE, norm_y = TRUE))
})

test_that("Simple case in dist_polysph()", {
  expect_equal(
    drop(dist_polysph(x = X1, y = X2, ind_dj = ind_dj, std = FALSE)),
    sqrt(
      acos(rowSums(X1[, 1:3] * X2[, 1:3]))^2 +
        acos(rowSums(X1[, 4:7] * X2[, 4:7]))^2
    ))
})

test_that("Distance between the same points is zero in dist_polysph()", {
  expect_equal(drop(dist_polysph(x = X1, y = X1, ind_dj = ind_dj,
                                 norm_x = TRUE, norm_y = TRUE)),
               rep(0, n))
  expect_equal(drop(dist_polysph(x = X2, y = X2, ind_dj = ind_dj,
                                 norm_x = TRUE, norm_y = TRUE)),
               rep(0, n))
  expect_equal(drop(dist_polysph(x = X1, y = X1, ind_dj = ind_dj,
                                 norm_x = TRUE, norm_y = TRUE)),
               drop(dist_polysph(x = X1, y = X1, ind_dj = ind_dj,
                                 norm_x = FALSE, norm_y = FALSE)))
  expect_equal(drop(dist_polysph(x = X2, y = X2, ind_dj = ind_dj,
                                 norm_x = TRUE, norm_y = TRUE)),
               drop(dist_polysph(x = X2, y = X2, ind_dj = ind_dj,
                                 norm_x = FALSE, norm_y = FALSE)))
})

## dist_polysph_cross()

test_that("dist_polysph_cross() equals dist_polysph_matrix() for x = y", {
  expect_equal(
    c(dist_polysph_matrix(x = X1, ind_dj = ind_dj, std = FALSE)),
    c(dist_polysph_cross(x = X1, y = X1, ind_dj = ind_dj, std = FALSE)),
    tolerance = 1e-6)
})

test_that("dist_polysph_matrix() handles a single row", {
  dm <- dist_polysph_matrix(x = rbind(c(1, 0)), ind_dj = c(0, 2))
  expect_equal(dim(dm), c(1L, 1L))
  expect_equal(unname(dm[1, 1]), 0)
})

test_that("dist_polysph_cross() equals dist_polysph()", {
  expect_equal(
    cbind(dist_polysph(x = X1[1:3, ], y = X2[1, , drop = FALSE],
                       ind_dj = ind_dj),
          dist_polysph(x = X1[1:3, ], y = X2[2, , drop = FALSE],
                       ind_dj = ind_dj)),
    dist_polysph_cross(x = X1[1:3, ], y = X2[1:2, ], ind_dj = ind_dj),
    tolerance = 1e-6)
  expect_equal(
    dist_polysph(x = X1[1, , drop = FALSE], y = X2[1, , drop = FALSE],
                 ind_dj = ind_dj),
    dist_polysph_cross(x = X1[1, , drop = FALSE], y = X2[1, , drop = FALSE],
                       ind_dj = ind_dj),
    tolerance = 1e-6)
})

## diamond_crossprod() and diamond_rcrossprod()

# Randomize testing
r <- 2
d <- c(2, 3)
n <- 3
X <- r_unif_polysph(n = n, d = d)
ind_dj <- comp_ind_dj(d)

test_that("Simple case in diamond_crossprod()", {
  for (i in seq_len(n)) {
    expect_equal(
      diamond_crossprod(X = X, ind_dj = ind_dj)[i, , ],
      rbind(cbind(tcrossprod(X[i, 1:3]), tcrossprod(X[i, 1:3], X[i, 4:7])),
            cbind(tcrossprod(X[i, 4:7], X[i, 1:3]), tcrossprod(X[i, 4:7])))
      )
  }
})

test_that("Simple case in diamond_rcrossprod()", {
  expect_equal(
    diamond_rcrossprod(X = X, ind_dj = ind_dj)[, , 1],
    tcrossprod(X[, 1:3]))
  expect_equal(
    diamond_rcrossprod(X = X, ind_dj = ind_dj)[, , 2],
    tcrossprod(X[, 4:7]))
})

## s()

test_that("Simmetrization using s()", {
  expect_true(isSymmetric(s(r_unif_polysph(n = 4, d = 3))))
})

## proj_P()

# Randomize testing
r <- rpois(1, lambda = 3) + 1
d <- rpois(r, lambda = 2) + 1
x <- r_unif_polysph(n = 1, d = d)
ind_dj <- comp_ind_dj(d)

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

test_that("AP does the job", {
  expect_equal(proj_P(x = x, d = d), AP(x = x, v = x, ind_dj = ind_dj)$P)
  expect_equal(diag(1, nrow = ncol(x), ncol = ncol(x)) - proj_P(x = x, d = d),
               AP(x = x, v = x, ind_dj = ind_dj, orth = TRUE)$P)
})

## polylog_minus_exp_mu()

k <- c(0.1, 0.5, 1, 10, 100)
test_that("polylog_minus_exp_mu() is smooth on integer indexes", {
  for (s in 1:5) {
    expect_equal(polylog_minus_exp_mu(s = s, mu = k),
                 polylog_minus_exp_mu(s = s + 1e-6, mu = k),
                 tolerance = 1e-5)
  }
})
test_that("polylog_minus_exp_mu() is smooth on half arguments", {
  for (s in c(0.5, 1.5, 2.5)) {
    expect_equal(polylog_minus_exp_mu(s = s, mu = k),
                 polylog_minus_exp_mu(s = s + 1e-6, mu = k),
                 tolerance = 1e-5)
  }
})
test_that("polylog_minus_exp_mu() edge cases", {
  expect_equal(polylog_minus_exp_mu(s = c(0.5, 1.5), mu = 1),
               c(polylog_minus_exp_mu(s = 0.5, mu = 1),
                 polylog_minus_exp_mu(s = 1.5, mu = 1)))
  expect_error(polylog_minus_exp_mu(s = 1:3, mu = 1:2))
})

## log_besselI_scaled()

# Exact / asymptotic references for the three computation regimes
besselI_ref <- function(nu, x) log(besselI(x, nu, expon.scaled = TRUE))
xasym_ref <- function(nu, x) {
  Bessel::besselIasym(x = x, nu = nu, k.max = 10, expon.scaled = TRUE,
                      log = TRUE)
}
nuasym_ref <- function(nu, x) {
  Bessel::besselI.nuAsym(x = x, nu = nu, k.max = 5, expon.scaled = TRUE,
                         log = TRUE)
}

# The vectorized call must equal the element-wise evaluation (with recycling)
expect_vectorized <- function(nu, x, spline = FALSE) {
  nx <- max(length(nu), length(x))
  nu_i <- rep_len(nu, nx)
  x_i <- rep_len(x, nx)
  expect_equal(log_besselI_scaled(nu = nu, x = x, spline = spline),
               sapply(seq_len(nx), function(i)
                 log_besselI_scaled(nu = nu_i[i], x = x_i[i], spline = spline)))
}

# spline = TRUE must agree with the direct spline = FALSE evaluation
expect_spline_agrees <- function(nu, x, tolerance = 1e-9) {
  expect_equal(log_besselI_scaled(nu = nu, x = x, spline = TRUE),
               log_besselI_scaled(nu = nu, x = x, spline = FALSE),
               tolerance = tolerance)
}

test_that("Correct vectorizations on nu and x for spline = FALSE", {
  xs <- 1:10
  nus <- c(1:9, 101)
  nus_a <- 101:102
  expect_vectorized(nus, xs)
  expect_vectorized(nus, xs[4])
  expect_vectorized(nus[4], xs)
  expect_vectorized(nus_a, xs[4])
})

test_that("Correct vectorizations on nu and x for spline = FALSE and NA's", {
  xs <- c(1, NA, 3:4, NA)
  nus <- 1:5
  expect_vectorized(nus, xs)
  for (j in 1:2) {
    expect_vectorized(nus, xs[j])
    expect_vectorized(nus[j], xs)
  }
})

test_that("log_besselI_scaled() with NAs", {
  x <- c(1, 2, NA, 1e5)
  expect_spline_agrees(1, x)
  expect_spline_agrees(5.5, x)
})

test_that("Accuracy of log_besselI_scaled(nu = seq(0, 24.5, by = 0.5)) over
          spline and asymptotic x grids", {
  for (x in list(10^seq(-5, 4, length.out = 5000), seq(1e4, 1e5, l = 100))) {
    for (nu in seq(0, 24.5, by = 0.5)) {
      expect_spline_agrees(nu, x)
    }
  }
})

test_that("Small-argument of log_besselI_scaled(spline = TRUE)", {
  x <- c(1e-8, 1e-7, 1e-6)
  for (nu in c(0, 0.5, 1, 2)) {
    ref <- besselI_ref(nu, x)
    expect_true(all(is.finite(ref)))
    expect_equal(log_besselI_scaled(nu = nu, x = x, spline = TRUE), ref,
                 tolerance = 1e-8)
  }
  expect_equal(log_besselI_scaled(nu = 0, x = 0, spline = TRUE), 0)
  expect_identical(log_besselI_scaled(nu = 5, x = 0, spline = TRUE), -Inf)
})

test_that("Asymptotic-kappa Bessel approximation", {
  paper_asymp <- function(x, d) {
    log1p(-d * (d - 2) / (8 * x)) - log(2 * pi * x) / 2
  }
  Bessel_asymp <- function(x, d) {
    Bessel::besselIasym(x = x, nu = (d - 1) / 2, expon.scaled = TRUE,
                        log = TRUE, k.max = 1)
  }
  for (d in 1:10) {
    expect_equal(paper_asymp(x = c(50:100, 1e4, 1e5), d = d),
                 Bessel_asymp(x = c(50:100, 1e4, 1e5), d = d))
  }
})

test_that("Asymptotic-d Bessel approximation", {
  for (x in c(0.1, 1, 10, 100)) {
    expect_lt(
      max(abs(diff(log_besselI_scaled(nu = 95:105, x = x, spline = FALSE),
                   differences = 2))),
      max(abs(diff(log_besselI_scaled(nu = 85:95, x = x, spline = FALSE),
                   differences = 2)))
    )
  }
})

test_that("Mix of asymptotic-d and regular Bessel computations", {
  nu <- seq(1, 200, by = 1)
  expect_no_error(log_besselI_scaled(nu = nu, x = nu))
})

test_that("Edge cases of log_besselI_scaled()", {
  expect_error(log_besselI_scaled(nu = NA, x = 1, spline = FALSE))
  expect_no_error(log_besselI_scaled(nu = 1:3, x = NA, spline = FALSE))
  expect_error(log_besselI_scaled(nu = 1:3, x = 1:2, spline = FALSE))
  expect_error(log_besselI_scaled(nu = c(1, NA), x = 1:2, spline = FALSE))
  expect_no_error(log_besselI_scaled(nu = 25, x = 1, spline = TRUE))
  expect_no_error(log_besselI_scaled(nu = 100, x = 0, spline = TRUE))
  expect_no_error(log_besselI_scaled(nu = 1:3, x = 0, spline = TRUE))
})

test_that("spline = FALSE handles large x (nu <= 100) without underflow", {
  nu <- 5
  x <- c(1e5, 5e5, 1e6)
  res <- log_besselI_scaled(nu = nu, x = x, spline = FALSE)
  expect_true(all(is.finite(res)))
  expect_spline_agrees(nu, x)
  expect_equal(res, xasym_ref(nu, x), tolerance = 1e-9)
})

test_that("nu-asymptotics take priority over x-asymptotics when both large", {
  nu <- 300
  x <- 1e4
  res <- log_besselI_scaled(nu = nu, x = x, spline = FALSE)
  expect_equal(res, nuasym_ref(nu, x), tolerance = 1e-9)
  expect_false(isTRUE(all.equal(res, xasym_ref(nu, x), tolerance = 1e-6)))
})

test_that("Vectorized spline = FALSE spans standard, x-asymp and nu-asymp", {
  nu <- c(5, 300, 5)
  x <- c(2, 1e4, 1e5)
  expect_vectorized(nu, x)
  ref <- c(besselI_ref(5, 2), nuasym_ref(300, 1e4), xasym_ref(5, 1e5))
  expect_equal(log_besselI_scaled(nu = nu, x = x, spline = FALSE), ref,
               tolerance = 1e-9)
})

test_that("Standard spline = FALSE path matches base besselI (scaled)", {
  nu <- 3
  x <- c(1, 5, 50, 500)
  expect_equal(log_besselI_scaled(nu = nu, x = x, spline = FALSE),
               besselI_ref(nu, x))
})

test_that("No discontinuity in spline = FALSE across the x = 1e4 boundary", {
  nu <- 5
  x <- seq(1e4 - 5, 1e4 + 5, length.out = 41)
  res <- log_besselI_scaled(nu = nu, x = x, spline = FALSE)
  expect_equal(res, xasym_ref(nu, x), tolerance = 1e-9)
  expect_true(all(diff(res) < 0))
  expect_lt(max(abs(diff(res, differences = 2))), 1e-6)
})

test_that("spline = TRUE silently falls back for unsupported nu", {
  x <- c(1e-6, 1e-2, 1, 50, 1e3, 1e5)
  for (nu in c(25, 30.5, 200)) {
    expect_spline_agrees(nu, x)
  }
  expect_equal(log_besselI_scaled(nu = 200, x = 1e4, spline = TRUE),
               nuasym_ref(200, 1e4), tolerance = 1e-9)
})

test_that("spline = TRUE accepts vectorized nu (silent fallback)", {
  nu <- c(1, 5.5, 25, 300)
  x <- c(2, 1e3, 1e5, 1e4)
  expect_spline_agrees(nu, x)
  expect_vectorized(nu, x)
})

## log_sum_exp()

test_that("log_sum_exp() is correct", {
  x <- c(1, 2, 3)
  expect_equal(log_sum_exp(logs = x), log(sum(exp(x))))
  expect_equal(log_sum_exp(logs = x, avg = TRUE), log(mean(exp(x))))
})

## log1m_exp()

test_that("log1m_exp() is correct", {
  x <- c(1e-5, 10, NA)
  expect_equal(log1m_exp(x), log(1 - exp(-x)))
  expect_equal(log1m_exp(x[1]), log(1 - exp(-x[1])))
  expect_equal(log1m_exp(x[2]), log(1 - exp(-x[2])))
})

## log_diff_exp()

test_that("log_diff_exp() is correct", {
  log_p <- c(-5, 1, 5, NA)
  log_n <- rev(log_p)
  log_diff <- log_diff_exp(log_p = log_p, log_n = log_n)
    expect_equal(log_diff$sgn * exp(log_diff$log_abs),
                 exp(log_p) - exp(log_n))
  for (i in 1:3) {
    log_diff <- log_diff_exp(log_p = log_p[i], log_n = log_n[i])
    expect_equal(log_diff$sgn * exp(log_diff$log_abs),
                 exp(log_p[i]) - exp(log_n[i]))
  }
})

test_that("log_diff_exp() respects the near-zero tolerance", {
  ld <- log_diff_exp(log_p = c(1, 1 + 1e-16, 2), log_n = c(1, 1, 1))
  expect_equal(ld$sgn, c(0, 0, 1))
})

test_that("log_diff_exp() returns a null difference on both sides of zero", {

  log_p <- c(1, 1, 1 + 3e-16)
  log_n <- c(1, 1 + 3e-16, 1)
  ld <- expect_no_warning(log_diff_exp(log_p = log_p, log_n = log_n))
  expect_equal(ld$sgn, c(0, 0, 0))
  expect_equal(ld$log_abs, rep(-Inf, 3))
  expect_equal(ld$sgn * exp(ld$log_abs), rep(0, 3))

})

## asinh_log()

test_that("asinh_log() works fine for large positive log_abs arguments", {
  log_abs <- c(-20, -10, 0, 50, 500, NA)
  sgn <- c(1, -1, 1, -1, 1, 1)
  x <- -3:3
  expect_equal(asinh(sgn * exp(log_abs)),
               asinh_log(log_abs = log_abs, sgn = sgn))
  expect_equal(x,
               sinh(polykde:::asinh_log(log_abs = log(abs(x)), sgn = sign(x))))
  for (i in 1:3) {
    expect_equal(asinh(sgn[i] * exp(log_abs[i])),
                 asinh_log(log_abs = log_abs[i], sgn = sgn[i]))
  }
})

test_that("asinh_log() works fine for large negative log_abs arguments", {
  log_abs <- c(-50, -500, -1000, NA)
  sgn <- c(1, -1, 1, NA)
  x <- c(1e-22, -1e-22, 1e-43)
  expect_equal(asinh(sgn * exp(log_abs)),
               asinh_log(log_abs = log_abs, sgn = sgn))
})
