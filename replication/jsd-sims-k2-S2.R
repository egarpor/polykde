
## Functions

# Set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Required libraries
library(foreach)
library(progressr)
library(future)
library(doFuture)
library(polykde)
library(rgl)
library(rotasym)
library(sphunif)
library(viridis)

# Nice sphere
sphere1.f <- function(x0 = 0, y0 = 0, z0 = 0, r = 1, n = 101, ...){
  f <- function(s,t){
    cbind(r * cos(t) * cos(s) + x0,
          r * sin(s) + y0,
          r * sin(t) * cos(s) + z0)
  }
  persp3d(f, slim = c(-pi / 2, pi / 2), tlim = c(0, 2 * pi), n = n, add = TRUE,
          ...)
}

## Samplers

dgp <- 1

### DGP 1: Simultaneous failure of mean and scatter tests

if (dgp == 1) {

# Sampler of the mixture f_a = (1 - a) / 2 * f_1 + (1 + a) / 2 * f_2
r_f_a <- function(n, d, kappa1 = 25, kappa2 = 50, nu1 = 0, F_inv_1 = NULL,
                  a = 0) {

  # Sample labels and sample sizes of each component
  stopifnot(length(d) == 1)
  p_a <- (1 - a) / 2
  stopifnot(p_a >= 0 && p_a <= 1)
  n1 <- sum(runif(n) < p_a)
  n2 <- n - n1

  # Subdivide samples of the second component
  n2 <- table(sample(x = 1:4, size = n2, replace = TRUE))

  # Samples from each component
  samp_1 <- switch((n1 >= 1) + 1, NULL,
                   rbind(r_alt(n = n1, p = d + 1, alt = "SC", kappa = kappa1,
                               nu = nu1, F_inv = F_inv_1)[, , 1]))
  samp_2_1 <- switch((n2[1] >= 1) + 1, NULL,
                     rbind(r_vMF(n = n2[1], mu = c(1, 0, 0), kappa = kappa2)))
  samp_2_2 <- switch((n2[2] >= 1) + 1, NULL,
                     rbind(r_vMF(n = n2[2], mu = c(-1, 0, 0), kappa = kappa2)))
  samp_2_3 <- switch((n2[3] >= 1) + 1, NULL,
                     rbind(r_vMF(n = n2[3], mu = c(0, 1, 0), kappa = kappa2)))
  samp_2_4 <- switch((n2[4] >= 1) + 1, NULL,
                     rbind(r_vMF(n = n2[4], mu = c(0, -1, 0), kappa = kappa2)))

  # Merging
  samp <- rbind(samp_1, samp_2_1, samp_2_2, samp_2_3, samp_2_4)
  return(samp[sample(n), , drop = FALSE])

}

# Sampler of H_a. Samples from the two groups defined by
# * n1 observations from f_1a = (1 - a) / 2 * f_1 + (1 + a) / 2 * f_2
# * n2 observations from f_2a = (1 + a) / 2 * f_1 + (1 - a) / 2 * f_2
r_hyp_a <- function(n1, n2, d, kappa1 = 25, kappa2 = 50, nu1 = 0,
                    F_inv_1 = NULL, a = 0) {

  X <- rbind(r_f_a(n = n1, d = d, kappa1 = kappa1, kappa2 = kappa2,
                   nu1 = nu1, F_inv_1 = F_inv_1, a = a),
             r_f_a(n = n1, d = d, kappa1 = kappa1, kappa2 = kappa2,
                   nu1 = nu1, F_inv_1 = F_inv_1, a = -a))
  labels <- rep(x = 1:2, times = c(n1, n2))
  ind <- sample(n1 + n2)
  X <- X[ind, , drop = FALSE]
  labels <- labels[ind]
  return(list("X" = X, "labels" = labels))

}

### DGP 2: Scatter fails, mean works

} else if (dgp == 2) {

# Sampler of H_a. Samples from the two groups defined by
# * n1 observations from f_1a = f_SC(nu = a)
# * n2 observations from f_2a = f_SC(nu = -a)
r_hyp_a <- function(n1, n2, d, kappa = 25, a = 0,
                    F_inv_a1 = NULL, F_inv_a2 = NULL) {

  stopifnot(a >= 0 && a <= 1)
  a <- ifelse(a == 1, 1 - 1e-15, a)
  X <- rbind(r_alt(n = n1, p = d + 1, alt = "SC", kappa = kappa, nu = a,
                   F_inv = F_inv_a1)[, , 1],
             r_alt(n = n2, p = d + 1, alt = "SC", kappa = kappa, nu = -a,
                   F_inv = F_inv_a2)[, , 1])
  labels <- rep(x = 1:2, times = c(n1, n2))
  ind <- sample(n1 + n2)
  X <- X[ind, , drop = FALSE]
  labels <- labels[ind]
  return(list("X" = X, "labels" = labels))

}

### DGP 3: Mean fails, scatter works

} else if (dgp == 3) {

# Sampler of the mixture f_a = (1 - a) / 2 * f_1 + (1 + a) / 2 * f_2
r_f_a <- function(n, d, kappa1 = 25, kappa2 = 150, nu2 = 0.75,
                  F_inv_kappa2 = NULL, a = 0) {

  # Sample labels and sample sizes of each component
  stopifnot(length(d) == 1)
  p_a <- (1 - a) / 2
  stopifnot(p_a >= 0 && p_a <= 1)
  n1 <- sum(runif(n) < p_a)
  n2 <- n - n1

  # Samples from each component
  samp_1 <- switch((n1 >= 1) + 1, NULL,
                   r_vMF(n = n1, mu = c(rep(0, d), 1), kappa = kappa1))
  samp_2 <- switch((n2 >= 1) + 1, NULL,
                   rbind(r_alt(n = n2, p = d + 1, alt = "SC", kappa = kappa2,
                               nu = nu2, F_inv = F_inv_kappa2)[, , 1]))

  # Merging
  samp <- rbind(samp_1, samp_2)
  return(samp[sample(n), , drop = FALSE])

}

# Sampler of H_a. Samples from the two groups defined by
# * n1 observations from f_1a = (1 - a) / 2 * f_1 + (1 + a) / 2 * f_2
# * n2 observations from f_2a = (1 + a) / 2 * f_1 + (1 - a) / 2 * f_2
r_hyp_a <- function(n1, n2, d, kappa1 = 25, kappa2 = 150, nu2 = 0.75,
                    F_inv_kappa2 = NULL, a = 0) {

  X <- rbind(r_f_a(n = n1, d = d, kappa1 = kappa1, kappa2 = kappa2, nu2 = nu2,
                   F_inv_kappa2 = F_inv_kappa2, a = a),
             r_f_a(n = n2, d = d, kappa1 = kappa1, kappa2 = kappa2, nu2 = nu2,
                   F_inv_kappa2 = F_inv_kappa2, a = -a))
  labels <- rep(x = 1:2, times = c(n1, n2))
  ind <- sample(n1 + n2)
  X <- X[ind, , drop = FALSE]
  labels <- labels[ind]
  return(list("X" = X, "labels" = labels))

}

}

# Visual check dispersion
open3d()
mfrow3d(nr = 1, nc = 3, sharedMouse = TRUE)
for (a in seq(0, 1, l = 3)) {

  samp <- switch(dgp,
                 r_hyp_a(n1 = 1e3, n2 = 1e3, d = 2, a = a),
                 r_hyp_a(n1 = 1e3, n2 = 1e3, d = 2, a = a),
                 r_hyp_a(n1 = 1e3, n2 = 1e3, d = 2, a = a))
  plot3d(samp$X, col = samp$labels + 2, xlim = c(-1, 1), ylim = c(-1, 1),
         zlim = c(-1, 1), xlab = "", ylab = "", zlab = "")

}

# Save figures
P <- matrix(c(0.888265907764435, -0.183141931891441, 0.421239256858826,
              0, 0.458672642707825, 0.402677863836288, -0.79213011264801, 0,
              -0.0245515555143356, 0.896833121776581, 0.441687434911728, 0,
              0, 0, 0, 1), nrow = 4, ncol = 4)
windowRect <- c(80, 80, 980, 980)
highlight <- c(88, 122, 146, 157)
for (a in seq(0, 1, l = 3)) {

  open3d()
  samp <- switch(dgp,
                 r_hyp_a(n1 = 1e3, n2 = 1e3, d = 2, a = a),
                 r_hyp_a(n1 = 1e3, n2 = 1e3, d = 2, a = a),
                 r_hyp_a(n1 = 1e3, n2 = 1e3, d = 2, a = a))
  plot3d(samp$X, col = c(2, 1)[samp$labels], xlim = c(-1, 1), ylim = c(-1, 1),
         zlim = c(-1, 1), size = 8, axes = FALSE, box = FALSE,
         lit = FALSE, xlab = "", ylab = "", zlab = "")
  sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 0.99, col = "lightblue",
            lit = TRUE, n = 201, alpha = 0.35)
  par3d(windowRect = windowRect, userMatrix = P, zoom = 0.4810174)
  snapshot3d(paste0("dgp_1_", a, ".png"), top = FALSE, webshot = FALSE,
             width = windowRect[3] - windowRect[1],
             height = windowRect[4] - windowRect[2])
  close3d()

}

# All tests are blind except JSD
hom_test_polysph(X = samp$X, d = 2, labels = samp$labels, type = "mean",
                 plot_boot = TRUE)
hom_test_polysph(X = samp$X, d = 2, labels = samp$labels, type = "scatter",
                 plot_boot = TRUE)
hom_test_polysph(X = samp$X, d = 2, labels = samp$labels, type = "jsd", B = 100,
                 plot_boot = TRUE)

## Simulations

# Define a progress bar
handlers(handler_progress(
  format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
  clear = FALSE))

# Monte Carlo setup
save <- TRUE
len_a <- 5
a <- seq(0, 1, l = len_a)
M <- 1e4
d <- 2
F_inv_1 <- F_inv_from_f(f = function(z) exp(-25 * z^2), p = d + 1)
F_inv_kappa2 <- F_inv_from_f(f = function(z) exp(-150 * (z - 0.75)^2),
                             p = d + 1)
F_inv_a1 <- F_inv_from_f(f = function(z) exp(-25 * z^2), p = d + 1)
F_inv_a2 <- F_inv_from_f(f = function(z) exp(-25 * z^2), p = d + 1)
n1 <- 100
n2 <- 100
type <- "jsd"
# type <- "mean"
# type <- "scatter"

# Average bandwidths
h_rot <- h_cv <- numeric(M)
progressr::with_progress({
  prog <- progressr::progressor(along = seq_len(M))
  for (i in seq_len(M)) {

    set.seed(1e6 + i)
    X <- switch(dgp,
                r_hyp_a(n1 = n1, n2 = n2, d = d, F_inv_1 = F_inv_1)$X,
                r_hyp_a(n1 = n1, n2 = n2, d = d, F_inv_a1 = F_inv_a1,
                        F_inv_a2 = F_inv_a2)$X,
                r_hyp_a(n1 = n1, n2 = n2, d = d,
                        F_inv_kappa2 = F_inv_kappa2)$X)
    h_rot[i] <- bw_rot_polysph(X = X, d = d)$bw
    h_cv[i] <- bw_cv_polysph(X = X, d = d, method = "L-BFGS-B",
                             bw0 = h_rot[i])$bw
    prog()

  }
})
(h_med_rot <- median(h_rot))
(h_med_cv <- median(h_cv))
h_med_rot / h_med_cv # 7.001773 (dgp = 1)

# Bandwidths for the tests (change)
(c <- 2^(-3:5)[5])
h_test <- c * h_med_cv

# File for saving & loading
file <- paste0("pow_S2_dgp", dgp, "_", type, "_c", sprintf("%.2f", c),
               "_M1e", log10(M))
# load(paste0(file, ".RData"))

# Warp-speed Monte Carlo
tests <- rep(list(), length = len_a)
options(future.rng.onMisuse = "ignore")
doFuture::registerDoFuture()
future::plan(future::multisession(), workers = 7)
for (j in seq_len(len_a)) {

  cat("\n", j, " / ", len_a, "\n")
  if (dgp == 2) {

    F_inv_a1 <- F_inv_from_f(f = function(z) exp(-25 * (z - a)^2), p = d + 1)
    F_inv_a2 <- F_inv_from_f(f = function(z) exp(-25 * (z + a)^2), p = d + 1)

  }
  progressr::with_progress({
    prog <- progressr::progressor(along = seq_len(M))
    tests[[j]] <- foreach(k = 1:M, .inorder = TRUE,
                          .packages = "polykde") %dopar% {

      if (requireNamespace("progressr", quietly = TRUE)) prog()
      set.seed(k, kind = "Mersenne-Twister")
      samp_hyp_a <- switch(dgp,
                           r_hyp_a(n1 = n1, n2 = n2, d = d, F_inv_1 = F_inv_1,
                                   a = a[j]),
                           r_hyp_a(n1 = n1, n2 = n2, d = d, a = a[j]),
                           r_hyp_a(n1 = n1, n2 = n2, d = d,
                                   F_inv_kappa2 = F_inv_kappa2, a = a[j]))
      hom_test_polysph(X = samp_hyp_a$X, d = d, labels = samp_hyp_a$labels,
                       type = type, h = h_test, B = 1, plot_boot = FALSE,
                       M = 1e4, cv_jsd = 1)

    }
  })

}

# Statistics (M x len_a)
stats <- sapply(tests, function(x) sapply(x, function(y)
  y$statistic))
stats_perm <- sapply(tests, function(x) sapply(x, function(y)
  y$statistic_perm[1]))
colnames(stats) <- colnames(stats_perm) <- a
boxplot(stats)

# Compute powers
alphas <- seq(0, 1, l = 101)
pow_ws <- t(sapply(alphas, function(alph, two = FALSE) {
  if (!two) {

    q_boot_alph_null_up <- quantile(x = stats_perm[, 1], probs = 1 - alph)
    return(rowMeans(t(stats) > q_boot_alph_null_up))

  } else {

    q_boot_alph2_null_up <- quantile(x = stats_perm[, 1], probs = 1 - alph / 2)
    q_boot_alph2_null_lo <- quantile(x = stats_perm[, 1], probs = alph / 2)
    return(rowMeans((t(stats) < q_boot_alph2_null_lo) |
                      (t(stats) > q_boot_alph2_null_up)))

  }
}))
rownames(pow_ws) <- alphas

# Plots
for (type in c("jsd", "mean", "scatter")) {
  for (c in 2^(-3:5)) {
    if (type != "jsd" && c != 1) {
      next
    } else {
      for (dgp in 1:3) {
        file <- paste0("pow_S2_dgp", dgp, "_", type, "_c", sprintf("%.2f", c),
                       "_M1e", log10(M))
        load(paste0(file, ".RData"))
if (save) pdf(file = paste0(file, ".pdf"), width = 6, height = 6)
matplot(alphas[-1], pow_ws[-1, ], type = "l", lty = 1, col = viridis(len_a),
        ylab = "Power", xlab = expression(alpha), lwd = 2, axes = FALSE,
        log = "xy", xlim = c(0.01, 1), ylim = c(0.01, 1), cex.lab = 1.25)
axis(1, at = c(0.01, 0.05, 0.10, 0.25, 0.5, 1), cex.axis = 1.25)
axis(2, at = c(0.01, 0.05, 0.10, 0.25, 0.5, 1), cex.axis = 1.25)
box()
abline(h = c(0.01, 1), lty = 3, col = "gray")
abline(v = c(0.01, 0.05, 0.10, 1), lty = 3, col = "gray")
abline(a = 0, b = 1, lty = 3)
lines(alphas, alphas - qnorm(0.975) * sqrt(alphas * (1 - alphas) / M),
      lty = 2)
lines(alphas, alphas + qnorm(0.975) * sqrt(alphas * (1 - alphas) / M), lty = 2)
points(rep(alphas[-1][5], len_a), pow_ws[-1, ][5, ], col = viridis(len_a),
       pch = 19)
segments(x0 = 1e-20, y0 = pow_ws[-1, ][5, ], x1 = 0.05, y1 = pow_ws[-1, ][5, ],
         col = viridis(len_a))
legend("bottomright", legend = paste("a =", round(rev(a), 2)), lwd = 2,
       col = viridis(len_a, direction = -1), bg = "white", cex = 1.5)
if (save) dev.off()
pow_ws
      }
    }
  }
}

# Save RData
if (save) save(list = ls(), file = paste0(file, ".RData"))
