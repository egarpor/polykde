
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

## Samplers

# Sampler of the mixture f_a = a * f_1 + b * f_2 + c * f_3
r_f_a <- function(n, d, r, kappa1 = 10, kappa2 = 10, kappa3 = 10,
                  a = rep(1 / 3, 3)) {

  # Sample labels and sample sizes of each component
  stopifnot(length(d) == 1)
  stopifnot(d >= 2)
  stopifnot(all(a >= 0 & a <= 1))
  stopifnot(abs(sum(a) - 1) < 1e-10 )
  labels <- sample(x = 1:3, size = n, replace = TRUE, prob = a)
  n1 <- sum(labels == 1)
  n2 <- sum(labels == 2)
  n3 <- n - n1 - n2

  # Generate means
  mu1 <- rep(c(1, 0, 0, rep(0, d - 2)), r)
  mu2 <- rep(c(0, 1, 0, rep(0, d - 2)), r)
  mu3 <- rep(c(0, 0, 1, rep(0, d - 2)), r)

  # Samples from each component
  samp_1 <- switch((n1 >= 1) + 1, NULL,
                   r_vmf_polysph(n = n1, d = rep(d, r), mu = mu1,
                                 kappa = rep(kappa1, r)))
  samp_2 <- switch((n2 >= 1) + 1, NULL,
                   r_vmf_polysph(n = n2, d = rep(d, r), mu = mu2,
                                 kappa = rep(kappa2, r)))
  samp_3 <- switch((n3 >= 1) + 1, NULL,
                   r_vmf_polysph(n = n3, d = rep(d, r), mu = mu3,
                                 kappa = rep(kappa3, r)))

  # Merging
  samp <- rbind(samp_1, samp_2, samp_3)
  return(samp[sample(n), , drop = FALSE])

}

# Sampler of H_a. Samples from the three groups defined by
# * n1 observations from f_1a = (1 + 2 * a) / 3 * f_1 + (1 - a) / 3 * (f_2 + f_3)
# * n2 observations from f_2a = (1 + 2 * a) / 3 * f_2 + (1 - a) / 3 * (f_1 + f_3)
# * n3 observations from f_3a = (1 + 2 * a) / 3 * f_3 + (1 - a) / 3 * (f_1 + f_2)
# Only the first two groups differ with moving a
r_hyp_a <- function(n1, n2, n3, d, r, kappa1 = 10, kappa2 = 10, kappa3 = 10,
                    a = 0) {

  X <- rbind(r_f_a(n = n1, d = d, r = r, kappa1 = kappa1, kappa2 = kappa2,
                   kappa3 = kappa3, a = (1 + c(2 * a, -a, -a)) / 3),
             r_f_a(n = n2, d = d, r = r, kappa1 = kappa1, kappa2 = kappa2,
                   kappa3 = kappa3, a = (1 + c(-a, 2 * a, -a)) / 3),
             r_f_a(n = n3, d = d, r = r, kappa1 = kappa1, kappa2 = kappa2,
                   kappa3 = kappa3, a = (1 + c(-a, -a, 2 * a)) / 3))
  labels <- rep(x = 1:3, times = c(n1, n2, n3))
  ind <- sample(n1 + n2 + n3)
  X <- X[ind, , drop = FALSE]
  labels <- labels[ind]
  return(list("X" = X, "labels" = labels))

}

# Visual check dispersion
open3d()
mfrow3d(nr = 2, nc = 2, sharedMouse = TRUE)
for (a in seq(0, 1, l = 4)) {

  samp <- r_hyp_a(n1 = 1e3, n2 = 1e3, n3 = 1e3, d = 2, r = 1,
                  kappa1 = 10, kappa2 = 10, kappa3 = 10, a = a)
  plot3d(samp$X, col = samp$labels + 1, xlim = c(-1, 1), ylim = c(-1, 1),
         zlim = c(-1, 1), xlab = "", ylab = "", zlab = "")

}

# All tests detect
hom_test_polysph(X = samp$X[samp$labels %in% c(1, 2), ], d = 2,
                 labels = samp$labels[samp$labels %in% c(1, 2)],
                 type = "mean", plot_boot = TRUE)
hom_test_polysph(X = samp$X[samp$labels %in% c(1, 2), ], d = 2,
                 labels = samp$labels[samp$labels %in% c(1, 2)],
                 type = "scatter", plot_boot = TRUE)
hom_test_polysph(X = samp$X[samp$labels %in% c(1, 2), ], d = 2,
                 labels = samp$labels[samp$labels %in% c(1, 2)],
                 type = "jsd", h = 0.25, B = 100, plot_boot = TRUE)

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
d <- 10
r <- 2
dd <- rep(d, r)
n1 <- 100
n2 <- 75
n3 <- 50
type <- "jsd"
# type <- "mean"
# type <- "scatter"

# Average bandwidths
h_rot <- h_cv <- matrix(nrow = M, ncol = r)
progressr::with_progress({
  prog <- progressr::progressor(along = seq_len(M))
  for (i in seq_len(M) {

    set.seed(1e6 + i)
    X <- r_hyp_a(n1 = n1, n2 = n2, n3 = n3, d = d, r = r)$X
    h_rot[i, ] <- bw_rot_polysph(X = X, d = dd)$bw
    h_cv[i, ] <- bw_cv_polysph(X = X, d = dd, method = "L-BFGS-B",
                               bw0 = h_rot[i, ])$par
    prog()

  }
})
(h_med_rot <- apply(h_rot, 2, median))
(h_med_cv <- apply(h_cv, 2, median))
mean(h_med_rot / h_med_cv) # 0.8925769

# Bandwidths for the tests
(c <- 2^(-3:5)[9])
h_test <- c * h_med_rot

# File for saving & loading
file <- paste0("pow_S10^2_", type, "_c", sprintf("%.2f", c), "_M1e", log10(M))
# load(paste0(file, ".RData"))

# Warp-speed Monte Carlo
tests <- rep(list(), length = len_a)
options(future.rng.onMisuse = "ignore")
doFuture::registerDoFuture()
future::plan(future::multisession(), workers = 7)
for (j in seq_len(len_a)) {

  cat("\n", j, " / ", len_a, "\n")
  progressr::with_progress({
    prog <- progressr::progressor(along = seq_len(M))
    tests[[j]] <- foreach(k = 1:M, .inorder = TRUE,
                          .packages = "polykde") %dopar% {

      if (requireNamespace("progressr", quietly = TRUE)) prog()
      set.seed(k, kind = "Mersenne-Twister")
      samp_hyp_a <- r_hyp_a(n1 = n1, n2 = n2, n3 = n3, d = d,
                            r = r, a = a[j])
      hom_test_polysph(X = samp_hyp_a$X, d = dd, labels = samp_hyp_a$labels,
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
for (type in c("jsd")) {
  for (c in 2^(-3:5)) {
    file <- paste0("pow_S10^2_", type, "_c", sprintf("%.2f", c), "_M1e",
                   log10(M))
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

# Save RData
if (save) save(list = ls(), file = paste0(file, ".RData"))
