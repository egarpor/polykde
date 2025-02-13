
r <- 2
d <- 3
n <- 10
h <- runif(1, min = 0.2, max = 1)
x <- drop(r_unif_polysph(n = 1, d = d))
y <- drop(r_unif_polysph(n = 1, d = d))
X <- r_unif_polysph(n = n, d = d)
X_MC <- rbind(x, r_vmf_polysph(n = 1e4, d = d, mu = x, kappa = 1 / h^2))
d_MC <- kde_polysph(x = X_MC, X = rbind(x), d = d, h = h)

d_r <- c(1, 2)
h_r <- runif(r, min = 0.2, max = 1)
x_r <- drop(r_unif_polysph(n = 1, d = d_r))
X_r <- r_unif_polysph(n = n, d = d_r)
X_MC_r <- rbind(x_r, r_vmf_polysph(n = 1e4, d = d_r, mu = x_r,
                                   kappa = 1 / h_r^2))
d_MC_r <- kde_polysph(x = X_MC_r, X = rbind(x_r), d = d_r, h = h_r)

## Product kernel constants

test_that("Kernel constants for vMF", {
  expect_equal(c_kern(h = h, d = d, kernel = 1) *
                 mean(L(t = drop(1 - X_MC %*% x) / h^2, kernel = 1) / d_MC),
               1, tolerance = 1e-1)
  expect_equal(c_kern(h = h, d = d, kernel = 1),
               1 / ((2 * pi)^((d + 1) / 2) *
                      besselI(x = 1 / h^2, nu = (d - 1) / 2) *
                      exp(-1 / h^2) * h^(d - 1)))
})

test_that("Kernel constants for Epa", {
  expect_equal(c_kern(h = h, d = d, kernel = 2) *
                 mean(L(t = drop(1 - X_MC %*% x) / h^2, kernel = 2) / d_MC),
               1, tolerance = 1e-1)
})

test_that("Kernel constants for sfp", {
  for (k in c(1, 5, 10, 100)) {
    expect_equal(c_kern(h = h, d = d, kernel = 3, k = k) *
                   mean(L(t = drop(1 - X_MC %*% x) / h^2, kernel = 3, k = k) /
                          d_MC),
                 1, tolerance = 1e-1)
  }
})

test_that("Kernel vMF constant for the specific case d = 2", {
  expect_equal(c_kern(h = h, d = 2, kernel = 1),
               exp((log(1 / h^2) - log(2 * pi)) - log1p(-exp(-2 / h^2))),
               tolerance = 1e-10)
})

## Intrinsic constants

test_that("Intrinsic kernel constants", {
  for (kern in 1:3) {
    expect_equal(c_kern(h = h, d = d, kernel = kern, intrinsic = TRUE) *
                   mean(L(t = acos(pmax(pmin(X_MC %*% x, 1), -1))^2 / (2 * h^2),
                          kernel = kern) / d_MC),
                 1, tolerance = 1e-1)
  }
})

test_that("Intrinsic kernel constants become extrinsic constants for small h", {
  h_small <- 0.2
  for (kern in 1:3) {
    expect_equal(c_kern(h = h_small, d = d, kernel = kern, k = 1,
                        intrinsic = TRUE, inc_sfp = FALSE),
                 c_kern(h = h_small, d = d, kernel = kern, k = 1,
                        intrinsic = FALSE, inc_sfp = FALSE),
                 tolerance = 1e-1)
  }
})

## Spherically symmetric kernel constants

test_that("Asymptotic sphericaly symmetric kernel constants become product
          kernel constants for small h and r = 1", {
  h_small <- 0.1
  for (kern in 1:3) {
    expect_equal(c_kern(h = h_small, d = d, kernel = kern, k = 10,
                             kernel_type = 1, inc_sfp = FALSE),
                 c_kern(h = h_small, d = d, kernel = kern, k = 10,
                        kernel_type = 2, inc_sfp = FALSE),
                 tolerance = 1e-1)
  }
})

test_that("Kernel constants for sphericaly symmetric vMF", {
  h_small <- h_r
  ind_dj <- comp_ind_dj(d_r)
  tj_r <- sapply(1:r, function(j) {
    ind_j <- (ind_dj[j] + 1):ind_dj[j + 1]
    drop(1 - X_MC_r[, ind_j] %*% x_r[ind_j]) / h_small[j]^2
  })
  expect_equal(c_kern(h = h_small, d = d_r, kernel = 1, kernel_type = 2) *
                 mean(L(t = rowSums(tj_r), kernel = 1) / d_MC_r),
               1, tolerance = 1e-1)
})

test_that("Kernel constants for sphericaly symmetric Epa", {
  h_small <- h_r
  ind_dj <- comp_ind_dj(d_r)
  tj_r <- sapply(1:r, function(j) {
    ind_j <- (ind_dj[j] + 1):ind_dj[j + 1]
    drop(1 - X_MC_r[, ind_j] %*% x_r[ind_j]) / h_small[j]^2
  })
  expect_equal(c_kern(h = h_small, d = d_r, kernel = 2, kernel_type = 2) *
                 mean(L(t = rowSums(tj_r), kernel = 2) / d_MC_r),
               1, tolerance = 1e-1)
})

test_that("Kernel constants for sphericaly symmetric sfp", {
  h_small <- h_r
  ind_dj <- comp_ind_dj(d_r)
  tj_r <- sapply(1:r, function(j) {
    ind_j <- (ind_dj[j] + 1):ind_dj[j + 1]
    drop(1 - X_MC_r[, ind_j] %*% x_r[ind_j]) / h_small[j]^2
  })
  expect_equal(c_kern(h = h_small, d = d_r, kernel = 3, k = 10,
                        kernel_type = 2) *
                 mean(L(t = rowSums(tj_r), kernel = 3, k = 10) / d_MC_r),
               1, tolerance = 1e-1)
})

test_that("Coherency between asymptotic and exact kernel constants for
          sphericaly symmetric Epa for d = 2 and h common", {
  r <- 10
  d2 <- rep(2, r)
  h_small <- rep(0.2, r)
  h_small_eps <- h_small + seq(-0.01, 0.01, l = r)
  expect_equal(c_kern(h = h_small, d = d2, kernel = 2, kernel_type = 2),
               c_kern(h = h_small_eps, d = d2, kernel = 2, kernel_type = 2),
               tolerance = 1e-1)
})

test_that("Coherency between exact kernel constants for sphericaly symmetric
          Epa for d = 2 and h common, cases above sqrt(2) and below", {
  r <- 10
  d2 <- rep(2, r)
  h_sqrt_minus <- rep(sqrt(2) - 1e-10, r)
  h_sqrt_plus <- rep(sqrt(2) + 1e-10, r)
  expect_equal(c_kern(h = h_sqrt_minus, d = d2, kernel = 2, kernel_type = 2),
               c_kern(h = h_sqrt_plus, d = d2, kernel = 2, kernel_type = 2),
               tolerance = 1e-1)
})

test_that("Coherency between asymptotic and exact kernel constants for
          sphericaly symmetric sfp for d = 2 and h common", {
  r <- 10
  d2 <- rep(2, r)
  h_small <- rep(0.2, r)
  h_small_eps <- h_small + seq(-0.01, 0.01, l = r)
  expect_equal(c_kern(h = h_small, d = d2, kernel = 3, kernel_type = 2, k = 1),
               c_kern(h = h_small_eps, d = d2, kernel = 3, kernel_type = 2,
                      k = 1),
               tolerance = 1e-1)
})

## lambda_d

dd <- 1:5

test_that("Coherency between lambda_h, b_d, and v_d", {
  for (di in dd) {
    for (ki in 1:3) {
      expect_equal(di * b_d(kernel = ki, d = di),
                   lambda_h(d = di, kernel = ki, bias = TRUE, h = NULL) /
                     lambda_h(d = di, kernel = ki, h = NULL),
                   tolerance = 1e-6)
      expect_equal(v_d(kernel = ki, d = di),
                   lambda_h(d = di, kernel = ki, squared = TRUE, h = NULL) /
                     lambda_h(d = di, kernel = ki, h = NULL)^2,
                   tolerance = 1e-6)
    }
  }
})

test_that("Coherency between lambda_h(h = NULL) and lambda_h(h = small)", {
  for (di in 1) {
    for (ki in 1:3) {
      h <- 0.1
      expect_equal(lambda_h(d = di, kernel = ki, h = NULL),
                   lambda_h(d = di, kernel = ki, h = h,
                            abs.tol = 1e-15, rel.tol = 1e-15,
                            subdivisions = 1e3, stop.on.error = FALSE),
                   tolerance = 1e-2)
    }
  }
})

test_that("Coherency between lambda_h and lambda_vmf_h", {
  for (di in dd) {
    for (hi in seq(0.05, 2, l = 10)) {
      expect_equal(lambda_h(d = di, kernel = 1, h = hi),
                   lambda_vmf_h(d = di, h = hi),
                   tolerance = 1e-3)
    }
  }
})

## Kernels moments for product kernels

dd <- 1:5
L_vMF <- function(t) exp(-t)
L_Epa <- function(t) (1 - t) * (0 <= t & t <= 1)
L_sfp <- function(t, k) log1p(exp(k * (1 - t))) / log1p(exp(k))

test_that("Kernel moments for product vMF", {
  for (di in dd) {
    expect_equal(b_d(kernel = 1, d = di), b_d(kernel = L_vMF, d = di),
                 tolerance = 1e-6)
    expect_equal(v_d(kernel = 1, d = di), v_d(kernel = L_vMF, d = di),
                 tolerance = 1e-6)
  }
})

test_that("Kernel moments for product Epa", {
  for (di in dd) {
    expect_equal(b_d(kernel = 2, d = di), b_d(kernel = L_Epa, d = di),
                 tolerance = 1e-6)
    expect_equal(v_d(kernel = 2, d = di), v_d(kernel = L_Epa, d = di),
                 tolerance = 1e-6)
  }
})

test_that("Kernel moments for product softplus", {
  for (di in dd) {
    for (ki in c(1, 10, 100)) {
      expect_equal(b_d(kernel = 3, d = di, k = ki),
                   b_d(kernel = function(x) L_sfp(x, k = ki), d = di, k = ki),
                   tolerance = 1e-6)
      expect_equal(v_d(kernel = 3, d = di, k = ki),
                   v_d(kernel = function(x) L_sfp(x, k = ki), d = di, k = ki),
                   tolerance = 1e-6)
    }
  }
})

## Kernels moments for spherically symmetric kernels

dd <- 1:5
rr <- 1:5

test_that("Kernel moments for spherically symmetric vMF", {
  expect_equal(b_d(kernel = 1, d = dd, kernel_type = "sph"),
               b_d(kernel = L_vMF, d = dd, kernel_type = "sph"),
               tolerance = 1e-6)
  expect_equal(v_d(kernel = 1, d = dd, kernel_type = "sph"),
               v_d(kernel = L_vMF, d = dd, kernel_type = "sph"),
               tolerance = 1e-6)
})

test_that("Kernel moments for spherically symmetric Epa", {
  expect_equal(b_d(kernel = 2, d = dd, kernel_type = "sph"),
               b_d(kernel = L_Epa, d = dd, kernel_type = "sph"),
               tolerance = 1e-6)
  expect_equal(v_d(kernel = 2, d = dd, kernel_type = "sph"),
               v_d(kernel = L_Epa, d = dd, kernel_type = "sph"),
               tolerance = 1e-6)
})

test_that("Kernel moments for spherically symmetric softplus", {
  for (ki in c(1, 10, 100)) {
    expect_equal(b_d(kernel = 3, d = dd, k = ki, kernel_type = "sph"),
                 b_d(kernel = function(x) L_sfp(x, k = ki), d = dd, k = ki,
                     kernel_type = "sph"),
                 tolerance = 1e-5)
    expect_equal(v_d(kernel = 3, d = dd, k = ki, kernel_type = "sph"),
                 v_d(kernel = function(x) L_sfp(x, k = ki), d = dd, k = ki,
                     kernel_type = "sph"),
                 tolerance = 1e-5)
  }
})

test_that("Product and spherically symmetric moments equal for vMF with
          r = 1", {

  for (di in dd) {
    expect_equal(b_d(kernel = "1", d = di, kernel_type = "prod"),
                 b_d(kernel = "1", d = di, kernel_type = "sph"))
    expect_equal(v_d(kernel = "1", d = di, kernel_type = "prod"),
                 v_d(kernel = "1", d = di, kernel_type = "sph"))
  }

})

test_that("Product and spherically symmetric moments equal for Epa with
          r = 1", {

  for (di in dd) {
    expect_equal(b_d(kernel = "2", d = di, kernel_type = "prod"),
                 b_d(kernel = "2", d = di, kernel_type = "sph"))
    expect_equal(v_d(kernel = "2", d = di, kernel_type = "prod"),
                 v_d(kernel = "2", d = di, kernel_type = "sph"))
  }

})

test_that("Product and spherically symmetric moments equal for softplus with
          r = 1", {

  for (di in dd) {
    expect_equal(b_d(kernel = "3", d = di, kernel_type = "prod"),
                 b_d(kernel = "3", d = di, kernel_type = "sph"))
    expect_equal(v_d(kernel = "3", d = di, kernel_type = "prod"),
                 v_d(kernel = "3", d = di, kernel_type = "sph"))
  }

})

## Intrinsic kernel moments

lambda_d_tilde <- function(d, kernel = "1", bias = FALSE, squared = FALSE,
                           k = 10) {

  rotasym::w_p(p = d) *
    integrate(function(t) L(t = t^2 / 2, kernel = kernel, squared = squared,
                            k = k) * t^(d - 1 + 2 * bias),
              lower = 0, upper = Inf)$value

}
b_d_tilde <- function(d, kernel = "1", k = 10) {

  num <- 0.5 * lambda_d_tilde(d = d, kernel = kernel, bias = TRUE,
                              squared = FALSE, k = k)
  den <- d * lambda_d_tilde(d = d, kernel = kernel, bias = FALSE,
                            squared = FALSE, k = k)
  return(num / den)

}
v_d_tilde <- function(d, kernel = "1", k = 10) {

  num <- lambda_d_tilde(d = d, kernel = kernel, bias = FALSE, squared = TRUE,
                        k = k)
  den <- lambda_d_tilde(d = d, kernel = kernel, bias = FALSE, squared = FALSE,
                        k = k)^2
  return(num / den)

}
lambda_d <- function(d, kernel = "1", bias = FALSE, squared = FALSE,
                     k = 10) {
  2^(d / 2 - 1) * rotasym::w_p(p = d) *
    integrate(function(t) L(t = t, kernel = kernel, squared = squared,
                            k = k) * t^(d / 2 - !bias),
            lower = 0, upper = Inf)$value

}

test_that("Intrinsic and extrinsic moments match for J(s) := L(s^2/2) or
          L(s) := J(√(2s))", {
  for (di in 1:5) {
    for (kern in 1:3) {
      expect_equal(lambda_d_tilde(d = di, kernel = kern),
                   lambda_d(d = di, kernel = kern), tolerance = 1e-5)
      expect_equal(b_d_tilde(d = di, kernel = kern),
                   b_d(d = di, kernel = kern), tolerance = 1e-5)
      expect_equal(v_d_tilde(d = di, kernel = kern),
                   v_d(d = di, kernel = kern), tolerance = 1e-5)
    }
  }
})

## Kernels efficiencies

dd <- 1:5
rr <- 1:5

effic_prop7 <- function(d, r, k = 10, type = "vmf") {

  if (type == "vmf") {

    eff <- exp((d * r + 2) * log(2) + lgamma(d * r / 2 + 2) -
                 (d * r / 2 + 1) * log(d * r + 4))

  } else if (type == "sfp-s" || type == "sfp-p") {

    log_ct <- (d * r / 2 + 2) * log(2) + lgamma(d * r / 2 + 2) -
      (d * r / 2 + 1) * log(d * r + 4) - (d * r / 2) * log(k)
    if (type == "sfp-s") {

      eff <- exp(log_ct + lgamma(d * r / 2))
      Li_num <- -polylog_minus_exp_mu(s = d * r / 2 + 1, mu = k)
      Li_den <- -polylog_minus_exp_mu(s = d * r / 2 + 2, mu = k)
      J <- J_d_k(d = d * r, k = k)
      eff <- eff * (Li_num^(d * r / 2 + 2) / (J * Li_den^(d * r / 2)))

    } else if (type == "sfp-p") {

      eff <- exp(log_ct + r * lgamma(d / 2))
      Li_num <- -polylog_minus_exp_mu(s = d / 2 + 1, mu = k)
      Li_den <- -polylog_minus_exp_mu(s = d / 2 + 2, mu = k)
      J <- J_d_k(d = d, k = k)^r
      eff <- eff * (Li_num^(r * (d / 2 + 2)) / (J * Li_den^(d * r / 2)))

    }

  } else if (type == "epa-p") {

    eff <- exp(r * (d / 2 + 1) * log(d + 4) + lgamma(d * r / 2 + 2) -
                 (r - 1) * log(4) - (d * r / 2 + 1) * log(d * r + 4) -
                 r * lgamma(d / 2 + 2))

  }
  return(eff)

}

ratio_eff_sfp_prop7 <- function(d, r, k = 10) {

  # ratio = \frac{
  # \Gamma(d/2)^r J_{dr}(k)
  #   |Li|_{dr/2+2}^{dr/2}(-e^k) |Li|_{d/2+1}^{r(d/2+2)}(-e^k)
  # }{
  # \Gamma(dr/2) J_d^r(k)
  #   |Li|_{dr/2+1}^{dr/2+2}(-e^k) |Li|_{d/2+2}^{dr/2}(-e^k)
  # }
  log_num <- r * lgamma(d / 2) + log(J_d_k(d = d * r, k = k)) +
    (d * r / 2) * log(abs(polylog_minus_exp_mu(s = d * r / 2 + 2,
                                               mu = k))) +
    (r * (d / 2 + 2)) * log(abs(polylog_minus_exp_mu(s = d / 2 + 1, mu = k)))
  log_den <- lgamma(d * r / 2) + r * log(J_d_k(d = d, k = k)) +
    (d * r / 2 + 2) * log(abs(polylog_minus_exp_mu(s = d * r / 2 + 1,
                                                   mu = k))) +
    (d * r / 2) * log(abs(polylog_minus_exp_mu(s = d / 2 + 2, mu = k)))
  return(exp(log_num - log_den))

}

test_that("Efficiencies from effic_prop7() and eff_kern() coincide", {

  for (ri in rr) {
    for (di in dd) {
      expect_equal(effic_prop7(d = di, r = ri, type = "vmf"),
                   eff_kern(kernel = 1, kernel_type = "sph", d = di, r = ri))
      expect_equal(effic_prop7(d = di, r = ri, k = ki, type = "epa-p"),
                   eff_kern(kernel = 2, kernel_type = "prod", d = di, r = ri))
      for (ki in c(1, 10, 100)) {
        expect_equal(effic_prop7(d = di, r = ri, k = ki, type = "sfp-s"),
                     eff_kern(kernel = 3, d = di, r = ri, k = ki,
                              kernel_type = "sph"))
        expect_equal(effic_prop7(d = di, r = ri, k = ki, type = "sfp-p"),
                     eff_kern(kernel = 3, d = di, r = ri, k = ki,
                              kernel_type = "prod"))
        expect_equal(effic_prop7(d = di, r = ri, k = ki, type = "sfp-p"),
                     effic_prop7(d = di, r = ri, k = ki, type = "sfp-s") *
                       ratio_eff_sfp_prop7(d = di, r = ri, k = ki))

      }
    }
  }

})

## Kernels gradients

test_that("Kernel gradients for vMF and softplus", {
  for (kernel in c(1, 3)) {
    expect_equal(
      grad_L(x = x, y = y, h = h, kernel = kernel),
      numDeriv::grad(func = function(x)
        L(t = drop(1 - x %*% y) / h^2, kernel = kernel, deriv = 0), x = x),
      tolerance = 1e-3)
  }
})

test_that("Kernel gradients for Epa", {
  expect_equal(
    grad_L(x = x, y = y, h = h, kernel = 2),
    numDeriv::grad(func = function(x)
      L(t = drop(1 - x %*% y) / h^2, kernel = 2, deriv = 0), x = x),
    tolerance = 1e-3)
})

## Kernels Hessians

test_that("Kernel Hessians for vMF and sfp", {
  for (kernel in c(1, 3)) {
        expect_equal(
      hess_L(x = x, y = y, h = h, kernel = kernel),
      numDeriv::hessian(func = function(x)
        L(t = drop(1 - x %*% y) / h^2, kernel = kernel, deriv = 0), x = x,
        method.args = list(eps = 1e-10)),
      tolerance = 1e-3)
  }
})

## Kde gradients

test_that("Kde gradient for vMF and sfp", {
  for (kernel in c(1, 3)) {
    expect_equal(
      drop(-c_kern(h = h, d = d, kernel = kernel) / h^2 * colMeans(
        L(t = drop(1 - X %*% x) / h^2, deriv = 1, kernel = kernel) * X
      )),
      numDeriv::grad(func = function(x)
        drop(kde_polysph(x = rbind(x), X = X, norm_x = FALSE, h = h, d = d,
                         kernel = kernel)), x = x),
      tolerance = 1e-3)
  }
})

test_that("Kde projected gradient for vMF and sfp", {
  for (kernel in c(1, 3)) {
    expect_equal(
      drop(-c_kern(h = h, d = d, kernel = kernel) / h^2 * colMeans(
        L(t = drop(1 - X %*% x) / h^2, deriv = 1, kernel = kernel) * X
        ) %*% (diag(rep(1, d + 1)) - tcrossprod(x))),
      numDeriv::grad(func = function(x)
        drop(kde_polysph(x = rbind(x), X = X, norm_x = TRUE, h = h, d = d,
                         kernel = kernel)), x = x),
      tolerance = 1e-3)
  }
})

## Kde Hessians

test_that("Kde Hessian for vMF and sfp", {
  for (kernel in c(1, 3)) {
    expect_equal(
      c_kern(h = h, d = d, kernel = kernel) / h^4 *
             matrix(rowMeans(sapply(1:n, function(i)
               L(t = drop(1 - X[i, ] %*% x) / h^2, deriv = 2, kernel = kernel) *
                 tcrossprod(X[i, ]))), nrow = d + 1, ncol = d + 1),
      numDeriv::hessian(func = function(x)
        drop(kde_polysph(x = rbind(x), X = X, norm_x = FALSE, h = h, d = d,
                         kernel = kernel)), x = x,
        method.args = list(eps = 1e-10)),
      tolerance = 5e-2)
  }
})

## Kernel sampling

test_that("Kernel sampling coherency between Epa and sfp in d = 2", {
  set.seed(12414)
  expect_true(
    ks.test(r_g_kern(n = 100, d = 2, h = 1, kernel = 2, k = 10),
            r_g_kern(n = 100, d = 2, h = 1, kernel = 2))$p.value > 0.1)
  expect_true(
    ks.test(r_g_kern(n = 100, d = 2, h = 0.5, kernel = 3, k = 10),
            r_g_kern(n = 100, d = 2, h = 0.5, kernel = 2))$p.value > 0.1)
  expect_false(
    ks.test(r_g_kern(n = 100, d = 2, h = 0.5, kernel = 3, k = 1),
            r_g_kern(n = 100, d = 2, h = 0.5, kernel = 2))$p.value > 0.1)
})
