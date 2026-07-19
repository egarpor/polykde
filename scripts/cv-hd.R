
library(polykde)
library(testthat)
library(rotasym)
library(movMF)
stopifnot(packageVersion("polykde") >= "1.1.8")

## polykde:::fast_log_c_vMF() works fine for large dimensions
{

# movMF's normalizing constant
log_c_vMF_movMF <- function(p, kappa, spline = FALSE) {
  -(movMF:::lH(kappa = kappa, nu = p / 2 - 1) + rotasym::w_p(p = p, log = TRUE))
}

my_odd_besselI_scaled <- function(x, m) {

  # 8.467 in Table of integrals
  k <- 0:m
  log_fact_k <- lfactorial(m + k) - (lfactorial(k) + lfactorial(m - k))
  log_const <- -0.5 * log(2 * pi * x)
  exp_k <- ((-1)^k + (-1)^(m + 1) * exp(-2 * x)) / (2 * x)^k
  sum(exp(log_const + log_fact_k) * exp_k)


  exp_k_A <- (-1)^k / (2 * x)^k
  exp_k_B <- 1 / (2 * x)^k
  exp(log_const) * (
    sum(exp(log_fact_k) * exp_k_A) +
    (-1)^(m + 1) * exp(-2 * x) * sum(exp(log_fact_k) * exp_k_B)
  )



  # log_2x <- k * log(2 * x)
  # log_k_x <- log_fact_k - log_2x
  # log_A_pos <- polykde:::log_sum_exp(log_k_x[seq(0, m, by = 2)])
  # log_A_neg <- polykde:::log_sum_exp(log_k_x[seq(1, m, by = 2)])
  # log_A <- polykde:::log_diff_exp(log_p = log_A_pos, log_n = log_A_neg)
  # exp(log_A$log_abs) * log_A$sgn
  # log_B <- -2 * x + polykde:::log_sum_exp(-log_k_x)
  # stopifnot(log_A$sgn >= 0)
  # if (m %% 2 == 0) { # Even
  #   log_result <- polykde:::log_diff_exp(log_p = log_A$log_abs,
  #                                        log_n = log_B)
  #   log_result <- log_result$log_abs
  # } else { # Odd
  #   log_result <- polykde:::log_sum_exp(c(log_A$log_abs, log_B))
  # }
  # exp(log_const + log_result)


}

besselI_scaled_half <- function(x, m, log = TRUE) {
  gsl::bessel_il_scaled(x = x, l = m) * sqrt(2 * x / pi)
  if (!log) {

    # Recursive computation
    # I_nu_p1 = I_nu_m1 - (2 * nu / x) * I_nu
    I_minus_0.5 <- 0.5 * sqrt(2 / (pi  * x)) * (1 + exp(-2 * x))
    I_plus_0.5 <- 0.5 * sqrt(2 / (pi  * x)) * (1 - exp(-2 * x))
    I_nu_minus_1 <- I_minus_0.5
    I_nu <- I_plus_0.5

    for (n in 0:(m - 1)) {
      I_nu_plus_1 <- I_nu_minus_1 - ((2 * n + 1) / x) * I_nu
      I_nu_minus_1 <- I_nu
      I_nu <- I_nu_plus_1
    }
    return(I_nu)

  } else {

    # Log-Recursive computation
    # log(I_nu_p1) = log_diff(log(I_nu_m1), log(2 * nu / x) + log(I_nu))
    log_I_minus_0.5 <- log(0.5 * sqrt(2 / (pi  * x))) + log1p(exp(-2 * x))
    log_I_plus_0.5 <- log(0.5 * sqrt(2 / (pi  * x))) + log(-expm1(-2 * x))
    log_I_nu_minus_1 <- log_I_minus_0.5
    log_I_nu <- log_I_plus_0.5
    for (k in 0:(m - 1)) {

      log_I_nu_plus_1 <- polykde:::log_diff_exp(
        log_p = log_I_nu_minus_1,
        log_n = log((2 * k + 1) / x) + log_I_nu)
      # print(log_I_nu_plus_1)
      print(k)
      stopifnot(log_I_nu_plus_1$sgn >= 0) # Can't happen unless precision is lost
      log_I_nu_plus_1 <- log_I_nu_plus_1$log_abs
      log_I_nu_minus_1 <- log_I_nu
      log_I_nu <- log_I_nu_plus_1

    }
    return(log_I_nu)

  }

}

p <- 1 + ceiling(10^seq(0, 4, by = 0.25))
kappa <- seq(0, 50, l = 100)
flog <- flog2 <- matrix(NA, nrow = length(kappa), ncol = length(p))
for (i in seq_along(p)) {
  flog[, i] <- polykde:::fast_log_c_vMF(p = p[i], kappa = p[i] * kappa)
  flog2[, i] <- log_c_vMF_movMF(p = p[i], kappa = p[i] * kappa)
}

matplot(log10(kappa), flog, type = "l", lty = 1, lwd = 2, col = 2,
        xlab = expression(kappa), ylab = expression(log(c[vmf](p, kappa))),
        main = "Log normalizing constant of vMF distribution",
        ylim = c(-1e3, 1e3))
matlines(log10(kappa), flog2, type = "l", lty = 2, col = 1)

}

## Evaluation of the CV for large dimensions
{

# vMF CV loss
exact_loss_cv <- function(h, X, spline = FALSE, arcsinh = TRUE, alt = FALSE,
                          Xi_Xj_l = NULL) {

  # Dimensions
  d <- ncol(X) - 1
  n <- nrow(X)
  r <- 1

  # Precompute matrix with the lower triangular parts of the matrices
  # (X_{il}'X_{jl})_{ij}, l = 1, ..., r.
  # ||X_i - X_j|| = sqrt(2 * (1 - X_i'X_j))
  #  => ||X_i - X_j||^2 = 2 * (1 - X_i'X_j)
  #  => ||X_i - X_j||^2 / 2 = 1 - X_i'X_j
  #  => X_i'X_j = 1 - ||X_i - X_j||^2 / 2
  ind_dj <- comp_ind_dj(d = d)
  if (is.null(Xi_Xj_l)) {

    Xi_Xj_l <- sapply(seq_len(r), function(j) {

      d_ij <- dist(X[, (ind_dj[j] + 1):ind_dj[j + 1]], method = "euclidean",
                   diag = FALSE, upper = FALSE)
      1 - 0.5 * d_ij^2

    })
    Xi_Xj_l <- t(Xi_Xj_l) # For better column recycling later

  }
  norm_Xi_Xj_l <- sqrt(2 * (1 + Xi_Xj_l))

  # Precompute other fixed objects in the LSCV loss
  log_n <- log(n)
  log_2_n1 <- log(2 / (n - 1))
  n2 <- n * (n - 1) / 2

  # Loop over bandwidths
  cv_loss <- numeric(length(h))
  for (i in seq_along(h)) {

    # Log-constants
    h_pos <- h[i]
    h_pos2 <- 1 / h_pos^2
    log_c_h2 <- sum(polykde:::fast_log_c_vMF(p = d + 1, kappa = h_pos2,
                                             spline = spline))
    log_c_2h2 <- sum(polykde:::fast_log_c_vMF(p = d + 1, kappa = 2 * h_pos2,
                                              spline = spline))

    # Compute X_{il}'X_{jl} / h_l^2 and
    # \sum_l \log(c_vMF(||X_{il}'X_{jl}|| / h_l^2))
    Xi_Xj_l_h <- .colSums(Xi_Xj_l * h_pos2, m = r, n = n2)
    log_c_norm_Xi_Xj_l_h <- numeric(n2)
    for (l in seq_len(r)) {

      log_c_norm_Xi_Xj_l_h <- log_c_norm_Xi_Xj_l_h +
        polykde:::fast_log_c_vMF(p = d[l] + 1, kappa = norm_Xi_Xj_l[l, ] *
                                   h_pos2[l], spline = spline)

    }

    # Compute arcsinh(CV) or CV loss?
    if (arcsinh) {

      # CV terms with LogSumExp trick
      log_cv_1 <- 2 * log_c_h2 - log_c_2h2 - log_n
      log_cv_2a <- log(2) +
        polykde:::log_sum_exp(Xi_Xj_l_h + log_c_h2 + (log_2_n1 - log_n))
      log_cv_2b <- log(2) +
        polykde:::log_sum_exp(2 * (log_c_h2 - log_n) - log_c_norm_Xi_Xj_l_h)

      # Positive and negative terms
      log_A <- polykde:::log_sum_exp(c(log_cv_1, log_cv_2b))
      log_B <- log_cv_2a

      # asinh-difference
      log_diff <- polykde:::log_diff_exp(log_p = log_A, log_n = log_B)
      asinh_cv <- polykde:::asinh_log(log_diff$log_abs, sgn = log_diff$sgn)
      loss <- asinh_cv

    } else if (alt) {

      # Common factor
      log_C <- log(2) + 2 * (log_c_h2 - log_n)

      # Terms in the sum
      log_A_1 <- rep(-(log_c_2h2 + log(n - 1)) + log_C, n2)
      log_A_2 <- -log_c_norm_Xi_Xj_l_h + log_C
      log_B <- log(2) + Xi_Xj_l_h - (log(1 - 1 / n) + log_c_h2) + log_C

      # Sum positives
      log_A <- apply(cbind(log_A_1, log_A_2), 1, polykde:::log_sum_exp)

      # Difference with negatives
      log_diff <- polykde:::log_diff_exp(log_p = log_A, log_n = log_B)

      # Sum all positives and negatives
      if (any(log_diff$sgn > 0)) {

        log_plus <- polykde:::log_sum_exp(log_diff$log_abs[log_diff$sgn > 0])

      } else {

        log_plus <- -Inf

      }
      if (any(log_diff$sgn < 0)) {

        log_minus <- polykde:::log_sum_exp(log_diff$log_abs[log_diff$sgn < 0])

      } else {

        log_minus <- -Inf

      }

      # Differences
      log_diff <- polykde:::log_diff_exp(log_p = log_plus, log_n = log_minus)
      loss <- polykde:::asinh_log(log_abs = log_diff$log_abs,
                                  sgn = log_diff$sgn)

    } else {

      # CV loss
      cv_1 <- exp(2 * log_c_h2 - log_c_2h2 - log_n)
      cv_2 <- 2 * sum(exp(Xi_Xj_l_h + log_c_h2 + (log_2_n1 - log_n)) -
                        exp(2 * (log_c_h2 - log_n) -
                              log_c_norm_Xi_Xj_l_h))
      cv <- cv_1 - cv_2
      loss <- cv

    }
    cv_loss[i] <- ifelse(!is.finite(loss), 1e6, loss)

  }
  return(list("h" = h, "cv_loss" = cv_loss, "Xi_Xj_l" = Xi_Xj_l))

}

# Gaussian CV loss
exact_loss_cv_gauss <- function(h, X, arcsinh = FALSE) {

  n <- nrow(X)
  d <- ncol(X)
  # Sinv <- solve(var(x))
  # Sinv <- (Sinv+t(Sinv))/2
  # Sinv2 <- Sinv%*%Sinv
  # hnorm <- ((4*d)/(det(Sinv)^(1/2)*(2*sum(diag(Sinv2))+sum(Sinv)^2)*n))^(1/(d+4))
  difs <- dist(X)
  difs2 <- difs^2
  if (!arcsinh) {

    edifs <- exp(-difs2 / 2)
    factor2 <- 2 / ((n^2 - n) * (2 * pi)^(d / 2))
    RK <- (4 * pi)^(-d / 2)
    RKn <- RK / n

    lscv.nd.temp <- function(h) {

      # lscv1 <- (1 - 1 / n) * sum(edifs^(1 / (2 * h^2))) / (2 * h^2)^(d / 2)
      # lscv2 <- 2 * sum(edifs^(1 / h^2)) / (h^2)^(d / 2)
      # lscv <- RKn / h^d + factor2 * (lscv1 - lscv2)
      lscv1 <- RKn / h^d + factor2 * (1 - 1 / n) * sum(edifs^(1 / (2 * h^2))) / (2 * h^2)^(d / 2)
      lscv2 <- factor2 * 2 * sum(edifs^(1 / h^2)) / (h^2)^(d / 2)
      lscv <- lscv1 - lscv2
      return(c(lscv, lscv1, lscv2))

    }

  } else {

    log_edifs <- -difs2 / 2
    log_factor2 <- log(2) - log(n^2 - n) - (d / 2) * log(2 * pi)
    log_RKn <- -(d / 2) * log(4 * pi) - log(n)

    lscv.nd.temp <- function(h) {

      log_lscv1_1 <- log_RKn - d * log(h)
      log_lscv1_2 <- log_factor2 + log(1 - 1 / n) +
        polykde:::log_sum_exp((1 / (2 * h^2)) * log_edifs) -
        (d / 2) * log(2 * h^2)
      log_lscv1 <- polykde:::log_sum_exp(c(log_lscv1_1, log_lscv1_2))

      log_lscv2 <- log_factor2 + log(2) +
        polykde:::log_sum_exp((1 / h^2) * log_edifs) - (d / 2) * log(h^2)

      lscv <- polykde:::log_diff_exp(log_p = log_lscv1, log_n = log_lscv2)
      loss <- polykde:::asinh_log(log_abs = lscv$log_abs, sgn = lscv$sgn)
      return(c(loss, log_lscv1, log_lscv2))

    }

  }
  cv_loss <- sapply(h, lscv.nd.temp)

  # opt <- nlm(f=lscv.nd.temp, p=hnorm)
  return(list("h" = h,
              "cv_loss" = cv_loss[1, ],
              "cv_loss_1" = cv_loss[2, ],
              "cv_loss_2" = cv_loss[3, ]))
}

# vMF CV loss via Monte Carlo
mc_loss_cv <- function(h, X, spline = FALSE, arcsinh = TRUE,
                       Xi_Xj_l = NULL, kernel = 1, seed_mc = 42, M = 1e4) {

  # Dimensions
  d <- ncol(X) - 1
  n <- nrow(X)
  r <- 1

  # Constant Epa wrt unif
  m_h <- pmax(-1, 1 - h^2)
  F_d <- sapply(seq_along(d), function(i)
    drop(sphunif::p_proj_unif(x = m_h[i], p = d[i] + 1)))
  const_epa <- (1 - h^(-2)) * (1 - F_d) +
    exp(rotasym::w_p(p = d, log = TRUE) - rotasym::w_p(p = d + 1, log = TRUE)) *
    (1 - m_h^2)^(d / 2) / (d * h^2)

  # Loop over bandwidths
  cv_loss <- numeric(length(h))
  for (i in seq_along(h)) {

    # Sum part with LogSumExp trick
    h_pos <- h[i]
    log_cv_kde <- log_cv_kde_polysph(X = X, d = d, h = h_pos, kernel = kernel,
                                     normalized = FALSE, wrt_unif = TRUE) -
      log(n) + log(2) + log(const_epa[i]) #+ polykde:::fast_log_c_vMF(p = d + 1, kappa = 1 / h_pos^2,
                      #                           spline = spline) + 1 / h_pos^2
    log_sum_cv_kde <- polykde:::log_sum_exp(log_cv_kde)

    # h-dependent Monte Carlo
    set.seed(seed_mc)
    mc_kde_samp <- r_kde_polysph(n = M, X = X, d = d,
                                 h = h_pos, kernel = kernel)
    log_kde2_mc <- kde_polysph(x = mc_kde_samp, X = X, d = d,
                               h = h_pos, kernel = kernel, log = TRUE,
                               normalized = FALSE, wrt_unif = TRUE) - log(M) +
      rotasym::w_p(p = p, log = TRUE) + log(const_epa[i]) #+
      #polykde:::fast_log_c_vMF(p = d + 1, kappa = 1 / h_pos^2,
      #                         spline = spline) + 1 / h_pos^2
    log_int_kde2 <- polykde:::log_sum_exp(log_kde2_mc)

    # Compute arcsinh(CV) or CV loss?
    if (arcsinh) {

      # arcsinh-loss
      log_diff <- polykde:::log_diff_exp(log_p = log_int_kde2,
                                         log_n = log_sum_cv_kde)
      asinh_cv <- polykde:::asinh_log(log_diff$log_abs, sgn = log_diff$sgn)
      loss <- asinh_cv

    } else {

      # CV loss
      cv <- exp(log_int_kde2) - exp(log_sum_cv_kde)
      loss <- cv

    }
    cv_loss[i] <- ifelse(!is.finite(loss), 1e6, loss)

  }
  return(list("h" = h, "cv_loss" = cv_loss, "Xi_Xj_l" = Xi_Xj_l))

}

# CV curve under perfect orthogonality
cv_0 <- function(h, p, n) {

  l1 <- polykde:::fast_log_c_vMF(p = p, kappa = 1 / h^2)
  l2 <- polykde:::fast_log_c_vMF(p = p, kappa = 2 / h^2)
  l3 <- polykde:::fast_log_c_vMF(p = p, kappa = sqrt(2) / h^2)
  log_cv_1_1 <- 2 * l1 - log(n) - l2
  log_cv_1_2 <- log(1 - 1 / n) + 2 * l1 - l3
  log_cv_1 <- apply(cbind(log_cv_1_1, log_cv_1_2), 1, polykde:::log_sum_exp)
  log_cv_2 <- log(2) + l1
  ld <- polykde:::log_diff_exp(log_p = log_cv_1, log_n = log_cv_2)
  polykde:::asinh_log(log_abs = ld$log_abs, sgn = ld$sgn)

}

# vMF example
n <- 100
p <- 100
kappa <- 10
Y <- rotasym::r_vMF(n = n, mu = c(1, rep(0, p - 1)), kappa = kappa)
# Y <- rotasym::r_unif_sphere(n = n, p = p)

h_grid <- 10^seq(-2, 1, by = 0.01)
cv_vmf <- exact_loss_cv(h = h_grid, X = Y, arcsinh = TRUE)
cv_vmf2 <- mc_loss_cv(h = h_grid, X = Y, M = 1e3, arcsinh = TRUE, kernel = 1)
plot(cv_vmf$h, cv_vmf$cv_loss, type = "o",
     xlab = "Bandwidth", ylab = "CV loss (arcsinh)",
     main = "Exact CV loss for vMF sample in high dimensions")

plot(cv_vmf2$h, cv_vmf2$cv_loss, type = "o", col = 4)
rug(cv_vmf$h)
abline(v = cv_vmf$h[which.min(cv_vmf$cv_loss)], col = 2, lty = 2)
# abline(v = bw_cv$bw, col = 3, lty = 2)
# The plot barely changes between samples...

# Compare with CV_0
curve(cv_0(x, n = n, p = p), n = 1e3, add = TRUE, col = 3)

# Gaussian example
n <- 100
p <- 100
sigma <- 0.1
Z <- rotasym::r_vMF(n = n, mu = c(1, rep(0, p - 1)), kappa = kappa)
# Z <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = sigma * diag(p))
h_grid <- 10^seq(-3, 1, by = 0.01)
cv_gauss_arc <- exact_loss_cv_gauss(h = h_grid, X = Z, arcsinh = TRUE)
plot(cv_gauss_arc$h, cv_gauss_arc$cv_loss, type = "o", log = "x")
# cv_gauss <- exact_loss_cv_gauss(h = h_grid, X = Y, arcsinh = FALSE)
# max(abs(cv_gauss_arc$cv_loss - asinh(cv_gauss$cv_loss)), na.rm = TRUE)
# lines(cv_gauss$h, asinh(cv_gauss$cv_loss), col = 4)
abline(v = cv_gauss_arc$h[which.min(cv_gauss_arc$cv_loss)],col=2)

# However, if the data is on the sphere...
r_err <- 0
cv_gauss_sph <- exact_loss_cv_gauss(h = h_grid,
                                    X = Y * (1 + rnorm(length(Y), sd = 0)),
                                    # X = Z / sqrt(rowSums(Z^2)),
                                    arcsinh = TRUE)
plot(cv_gauss_sph$h, cv_gauss_sph$cv_loss, type = "o")
lines(cv_gauss_sph$h, cv_gauss_sph$cv_loss_1, col = 3, type = "l")
lines(cv_gauss_sph$h, cv_gauss_sph$cv_loss_2, col = 2, type = "l")

# Check CV_0
cv_0(h = 1, p = 10, n = 2)
exact_loss_cv(h = 1,
              X = rbind(c(1, rep(0, 9)), c(rep(0, 9), 1)),
              arcsinh = TRUE)$cv_loss

# Compare CV vMF vs Gauss -- very similar!

n <- 500
p <- 2000
kappa <- 100
Y <- rbind(rotasym::r_vMF(n = n, mu = c(1, rep(0, p - 1)), kappa = kappa),
           rotasym::r_vMF(n = n, mu = c(-1, rep(0, p - 1)), kappa = kappa))
h_grid <- seq(0.01, 0.25, by = 0.005)
cv_vmf <- exact_loss_cv(h = h_grid, X = Y, arcsinh = TRUE)
cv_gauss <- exact_loss_cv_gauss(h = h_grid, X = Y, arcsinh = TRUE)
plot(cv_vmf$h, cv_vmf$cv_loss, type = "o",
     xlab = "Bandwidth", ylab = "CV loss (arcsinh)",
     main = "Exact CV loss for vMF sample in high dimensions")
lines(cv_gauss$h, cv_gauss$cv_loss, type = "o", col = 4)
abline(v = cv_vmf$h[which.min(cv_vmf$cv_loss)], col = 1, lty = 2)
abline(v = cv_gauss$h[which.min(cv_gauss$cv_loss)], col = 4, lty = 2)
rug(cv_vmf$h)
legend("topright", legend = c("vMF CV", "Gauss CV"),
       col = c(1, 4), pch = 1)

}

## Evaluation of the kde for large dimensions
{

kde_hd <- function(x, X, d, h, kernel, ...) {

  dens <- kde_polysph(x = x, X = X, d = d, h = h, kernel = kernel,
                      normalized = FALSE, log = TRUE, wrt_unif = TRUE, ...)
  if (kernel == 1)  {

    dens <- dens +
      sum(polykde:::fast_log_c_vMF(p = d + 1, kappa = 1 / h^2) + 1 / h^2)

  } else if (kernel == 2) {

    # Constant Epa wrt unif
    m_h <- pmax(-1, 1 - h^2)
    F_d <- sapply(seq_along(d), function(i)
      drop(sphunif::p_proj_unif(x = m_h[i], p = d[i] + 1)))
    const_epa <- (1 - h^(-2)) * (1 - F_d) +
      exp(rotasym::w_p(p = d, log = TRUE) - rotasym::w_p(p = d + 1, log = TRUE)) *
      (1 - m_h^2)^(d / 2) / (d * h^2)
    dens <- dens + sum(log(const_epa))

  }
  dens

}

test_that("Equivalence between kde_hd() and kde_polysph() for low dimension", {

  p <- 100
  n <- 50
  kappa <- 50
  x <- rotasym::r_vMF(n = 10, mu = c(1, rep(0, p - 1)), kappa = kappa)
  X <- rotasym::r_vMF(n = n, mu = c(1, rep(0, p - 1)), kappa = kappa)
  h <- 1
  kde1 <- kde_hd(x = x, X = X, d = p - 1, h = h, kernel = 1)
  kde2 <- kde_polysph(x = x, X = X, d = p - 1, h = h,
                      kernel = 1, log = TRUE, wrt_unif = TRUE)
  expect_equal(kde1, kde2)

})

}

## classic3
{

library(ZILGM)
data("classic3", package = "ZILGM")

# Extract labels and words
rownames(classic3) <- paste0("doc_", seq_len(nrow(classic3)))
doc_labels <- as.factor(classic3$doc_name)
names(doc_labels) <- rownames(classic3)
doc <- classic3
doc$doc_name <- NULL
doc <- as.matrix(doc)

# Filter words that do not appear often, relatively
doc_prop <- doc / rowSums(doc)
col_sums <- colSums(doc_prop)
summary(col_sums)
filt <- col_sums >= quantile(col_sums, probs = 0.1) &
  col_sums <= quantile(col_sums, probs = 1)
mean(filt)
doc_filt <- doc[, filt]
doc_filt <- doc_filt[rowSums(doc_filt > 0) > 0, ]
dim(doc_filt)
doc_filt <- doc_filt[, 1:100] # because yes

# Normalize
X_filt <- doc_filt / sqrt(rowSums(doc_filt^2))
X_filt <- X_filt[!duplicated(X_filt), ] # Causing big isssues on bwds!!!
labels <- doc_labels[names(doc_labels) %in% rownames(X_filt)]
stopifnot(labels == doc_labels[rownames(X_filt)])
dim(X_filt)

# Split into train + test
set.seed(42)
train_ind <- sample(nrow(X_filt), size = 0.8 * nrow(X_filt))
X_train <- X_filt[train_ind, ]
X_test <- X_filt[-train_ind, ]
labels_train <- labels[train_ind]
labels_test <- labels[-train_ind]

# Split train by class
X_train_cisi <- X_train[rownames(X_train) %in%
                          names(labels_train)[labels_train == "cisi"], ]
X_train_cran <- X_train[rownames(X_train) %in%
                          names(labels_train)[labels_train == "cran"], ]
X_train_med  <- X_train[rownames(X_train) %in%
                          names(labels_train)[labels_train == "med"], ]

#
#
# h_grid <- 10^seq(-2, 1, length.out = 20)
# cv <- exact_loss_cv(h = h_grid, X = X,
#                     spline = FALSE, arcsinh = TRUE)
# plot(h_grid, cv$cv_loss, type = "l", log = "x")
# rug(h_grid)
#
# h_grid2 <- c(1e-5, 1e-4, 1e-3, 1e3)
# cv2 <- exact_loss_cv(h = h_grid2, X = X,
#                     spline = FALSE, arcsinh = TRUE,
#                     Xi_Xj_l = cv$Xi_Xj_l)
# plot(h_grid2, cv2$cv_loss, type = "l", log = "x")
# rug(h_grid2)


# Fits for each group

h_grid <- 10^seq(-1, 1, length.out = 200) #10^seq(-8, 1, length.out = 200)
# cv_cisi <- exact_loss_cv(h = h_grid, X = X_train_cisi, arcsinh = TRUE)
# cv_cran <- exact_loss_cv(h = h_grid, X = X_train_cran, arcsinh = TRUE)
# cv_med <- exact_loss_cv(h = h_grid, X = X_train_med, arcsinh = TRUE)
cv_cisi <- exact_loss_cv_gauss(h = h_grid, X = X_train_cisi, arcsinh = TRUE)
cv_cran <- exact_loss_cv_gauss(h = h_grid, X = X_train_cran, arcsinh = TRUE)
cv_med <- exact_loss_cv_gauss(h = h_grid, X = X_train_med, arcsinh = TRUE)

plot(cv_cisi$h, cv_cisi$cv_loss, type = "l", log = "x")
lines(cv_cran$h, cv_cran$cv_loss, col = 2)
lines(cv_med$h, cv_med$cv_loss, col = 3)
rug(h_grid)

# h_grid <- 10^seq(-1, 1, length.out = 100)
# cv_cisi <- mc_loss_cv(h = h_grid, X = X_train_cisi, arcsinh = TRUE, M = 1e3, kernel = 2)
# cv_cran <- mc_loss_cv(h = h_grid, X = X_train_cran, arcsinh = TRUE, M = 1e3, kernel = 2)
# cv_med <- mc_loss_cv(h = h_grid, X = X_train_med, arcsinh = TRUE, M = 1e3, kernel = 2)
#
# plot(cv_cisi$h, cv_cisi$cv_loss, type = "l", log = "x")
# lines(cv_cran$h, cv_cran$cv_loss, col = 2)
# lines(cv_med$h, cv_med$cv_loss, col = 3)
# rug(h_grid)

# Fit kde to each class
log_kde_cisi <- kde_hd(x = X_test, X = X_train_cisi, d = ncol(X_train_cisi) - 1,
                       kernel = 3, h = 5e-1)# h = cv_cisi$h[which.min(cv_cisi$cv_loss)])
log_kde_cran <- kde_hd(x = X_test, X = X_train_cran, d = ncol(X_train_cran) - 1,
                       kernel = 3, h = 5e-1)# h = cv_cran$h[which.min(cv_cran$cv_loss)])
log_kde_med <- kde_hd(x = X_test, X = X_train_med, d = ncol(X_train_med) - 1,
                      kernel = 3, h = 5e-1)# h = cv_med$h[which.min(cv_med$cv_loss)])
log_kde_matrix <- cbind(cisi = c(log_kde_cisi),
                        cran = c(log_kde_cran),
                        med = c(log_kde_med))
head(log_kde_matrix)
pred_labels_kde <- colnames(log_kde_matrix)[max.col(log_kde_matrix)]
pred_labels_kde <- as.factor(pred_labels_kde)

# Compute confusion matrices
conf_kde <- table(True = labels_test, Predicted = pred_labels_kde)
print(conf_kde)
sum(diag(conf_kde)) / sum(conf_kde)

# Fit single-component vMF to each class
fit_movMF_cisi <- movMF(X_train_cisi, k = 1)
fit_movMF_cran <- movMF(X_train_cran, k = 1)
fit_movMF_med <- movMF(X_train_med, k = 1)
log_movMF_cisi <- dmovMF(X_test, theta = fit_movMF_cisi$theta,
                         alpha = fit_movMF_cisi$alpha, log = TRUE)
log_movMF_cran <- dmovMF(X_test, theta = fit_movMF_cran$theta,
                         alpha = fit_movMF_cran$alpha, log = TRUE)
log_movMF_med <- dmovMF(X_test, theta = fit_movMF_med$theta,
                        alpha = fit_movMF_med$alpha, log = TRUE)
log_movMF_matrix <- cbind(cisi = log_movMF_cisi,
                          cran = log_movMF_cran,
                          med = log_movMF_med)
pred_labels_movMF <- colnames(log_movMF_matrix)[max.col(log_movMF_matrix)]
pred_labels_movMF <- as.factor(pred_labels_movMF)

# Compute confusion matrices
conf_movMF <- table(True = labels_test, Predicted = pred_labels_movMF)
print(conf_movMF)
sum(diag(conf_movMF)) / sum(conf_movMF)

}

## MNIST
{

# Load dataset
load("qmnist_nist.RData")

# Filter and split train
train_nist_sd <- apply(train_nist$px, 2, sd)
train_nist_fil <- train_nist
sel_cols <- train_nist_sd > 1
train_nist_fil$px <- train_nist$px[, sel_cols]
train_nist_fil <- train_nist_fil[rowSums(train_nist_fil$px) > 0, ]
train_d <- list()
for (i in 0:9) {
  i_ch <- as.character(i)
  train_d[[i_ch]] <- as.matrix(train_nist_fil$px[train_nist_fil$digit == i, ])
  train_d[[i_ch]] <- train_d[[i_ch]] / sqrt(rowSums(train_d[[i_ch]]^2))
}

# CV bandwidths
h_grid <- seq(0.005, 0.1, by = 0.005)
cv_d <- list()
for (i in 0:9) {
  cat("Computing CV for digit", i, "\n")
  i_ch <- as.character(i)
  cv_d[[i_ch]] <- exact_loss_cv(h = h_grid, X = train_d[[i_ch]],
                                spline = FALSE, arcsinh = TRUE)
  plot(cv_d[[i_ch]]$h, cv_d[[i_ch]]$cv_loss, type = "o", log = "x",
       col = i + 1)
  rug(h_grid)
}

# KDA predictions
kda_hd <- function(x_new, train_d, cv_d) {

  kda <- matrix(NA, nrow = nrow(x_new), ncol = length(train_d))
  for (i in seq_along(train_d)) {
    cat("Evaluating class", names(train_d)[i], "\n")
    kda[, i] <- c(kde_hd(x = x_new, X = train_d[[i]], kernel = 1,
                         d = ncol(train_d[[i]]) - 1,
                         h = cv_d[[i]]$h[which.min(cv_d[[i]]$cv_loss)]))
  }
  c(0:9)[apply(kda, 1, which.max)]

}
test_fil <- test_nist$px[, sel_cols]
test_fil <- test_fil / sqrt(rowSums(test_fil^2))
test_fil <- as.matrix(test_fil)
ind_test <- sample(nrow(test_fil), size = 5e3)
pred_labels_kda <- kda_hd(x_new = test_fil, train_d = train_d, cv_d = cv_d)

# Confusion matrix
conf_kda <- table(True = test_nist$digit, Predicted = pred_labels_kda)
conf_kda
sum(diag(conf_kda)) / sum(conf_kda)

# Same approach with movMF

mvmfda <- function(x_new, train_d) {

  da <- matrix(NA, nrow = nrow(x_new), ncol = length(train_d))
  for (i in seq_along(train_d)) {
    cat("Evaluating class", names(train_d)[i], "\n")
    fit <- list()
    for (k in 1:20) {
      fit[[k]] <- tryCatch(movMF::movMF(x = train_d[[i]], k = k),
                           error = function(e) NA)
    }
    bics <- sapply(fit, function(x) ifelse(any(is.na(x)), NA, BIC(x)))
    fit_best <- fit[[which.min(bics)]]
    cat("Best model with", which.min(bics), "components\n")
    da[, i] <- dmovMF(x_new, theta = fit_best$theta, alpha = fit_best$alpha,
                      log = TRUE)
  }
  c(0:9)[apply(da, 1, which.max)]

}

pred_labels_mv <- mvmfda(x_new = test_fil, train_d = train_d)
conf_mv <- table(True = test_nist$digit, Predicted = pred_labels_mv)
conf_mv
sum(diag(conf_mv)) / sum(conf_mv)

}

# Compare with movMF
# jiffe faces?

## useR_2008_abstracts
{

library("tm")
library("slam")
data("useR_2008_abstracts", package = "movMF")


abstracts_titles <-
  apply(useR_2008_abstracts[,c("Title", "Abstract")],
        1, paste, collapse = " ")
useR_2008_abstracts_corpus <- VCorpus(VectorSource(abstracts_titles))
useR_2008_abstracts_DTM <-
  DocumentTermMatrix(useR_2008_abstracts_corpus,
                     control = list(
                       tokenize = "MC",
                       stopwords = TRUE,
                       stemming = TRUE,
                       wordLengths = c(3, Inf)))

ColSums <- col_sums(useR_2008_abstracts_DTM > 0)
sort(ColSums, decreasing = TRUE)[1:10]



useR_2008_abstracts_DTM <-
  useR_2008_abstracts_DTM[, ColSums >= 5 & ColSums <= 90]
useR_2008_abstracts_DTM



useR_2008_abstracts_DTM <- weightTfIdf(useR_2008_abstracts_DTM)


m <- movMF(useR_2008_abstracts_DTM,
           k = 10, nruns = 20)

dim(useR_2008_abstracts_DTM)

X_DTM <- as.matrix(useR_2008_abstracts_DTM)
X_DTM <- X_DTM / sqrt(rowSums(X_DTM^2))

cv_DTM <- exact_loss_cv(h = seq(0.01, 0.5, by = 0.001), X = X_DTM,
                spline = FALSE, arcsinh = TRUE)
plot(cv_DTM$h, cv_DTM$cv_loss, type = "o")



###################################################
### code chunk number 17: movMF.Rnw:1512-1525
###################################################
set.seed(2008)
library("movMF")
Ks <- c(1:5, 10, 20)
splits <- sample(rep(1:10, length.out = nrow(useR_2008_abstracts_DTM)))
useR_2008_movMF <-
  lapply(Ks, function(k)
    sapply(1:10, function(s) {
      m <- movMF(useR_2008_abstracts_DTM[splits != s,],
                 k = k, nruns = 20)
      logLik(m, useR_2008_abstracts_DTM[splits == s,])
    }))
useR_2008_movMF_common <-
  lapply(Ks, function(k)
    sapply(1:10, function(s) {
      m <- movMF(useR_2008_abstracts_DTM[splits != s,],
                 k = k, nruns = 20,
                 kappa = list(common = TRUE))
      logLik(m, useR_2008_abstracts_DTM[splits == s,])
    }))
if(cache) {
  save(useR_2008_movMF, useR_2008_movMF_common,
       file = "movMF.rda")
} else {
  if(file.exists("movMF.rda")) file.remove("movMF.rda")
}


###################################################
### code chunk number 18: movMF.Rnw:1538-1550
###################################################
logLiks <- data.frame(logLik = c(unlist(useR_2008_movMF),
                                 unlist(useR_2008_movMF_common)),
                      K = c(rep(Ks, sapply(useR_2008_movMF, length)),
                            rep(Ks, sapply(useR_2008_movMF_common, length))),
                      Dataset = seq_len(length(useR_2008_movMF[[1]])),
                      Method = factor(rep(1:2, each = length(unlist(useR_2008_movMF))),
                                      1:2, c("free", "common")))
logLiks$logLik <- logLiks$logLik - rep(rep(with(logLiks, tapply(logLik, Dataset, mean)), length(Ks)), 2)
print(xyplot(logLik ~ K | Method, data = logLiks, groups = Dataset, type = "l", lty = 1,
             xlab = "Number of components", ylab = "Predictive log-likelihood",
             strip = strip.custom(factor.levels  =
                                    expression(paste("Free ", kappa), paste("Common ", kappa)))))


###################################################
### code chunk number 19: movMF.Rnw:1563-1566
###################################################
set.seed(2008)
best_model <- movMF(useR_2008_abstracts_DTM, k = 2, nruns = 20,
                    kappa = list(common = TRUE))


###################################################
### code chunk number 20: movMF.Rnw:1572-1574
###################################################
apply(coef(best_model)$theta, 1, function(x)
  colnames(coef(best_model)$theta)[order(x, decreasing = TRUE)[1:10]])


###################################################
### code chunk number 21: movMF.Rnw:1587-1594
###################################################
clustering <- predict(best_model)
keywords <- useR_2008_abstracts[, "Keywords"]
keywords <- sapply(keywords,
                   function(x) sapply(strsplit(x, ", ")[[1]], function(y)
                     strsplit(y, "-")[[1]][1]))
tab <- table(Keyword = unlist(keywords),
             Component = rep(clustering, sapply(keywords, length)))


###################################################
### code chunk number 22: movMF.Rnw:1600-1601
###################################################
(tab <- tab[rowSums(tab) > 8, ])


###################################################
### code chunk number 23: movMF.Rnw:1614-1618
###################################################
library("vcd")
mosaic(tab, shade = TRUE,
       labeling_args = list(rot_labels = 0, just_labels = c("center", "right"),
                            pos_varnames = c("left", "right"), rot_varnames = 0))



}
