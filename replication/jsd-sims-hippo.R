
# This script reproduces Figures 11a-1h in the SM. Comment/uncomment "type" to
# run different tests ("jsd"/"mean"/"scatter"). The "Simulation" block runs the
# simulations (parallelized and time-consuming), and the "Analysis" block
# generates the figures. Faster outputs can be obtained by reducing M = 1e4 to
# a smaller value, e.g., M = 5e2.

## Functions

# Required libraries
library(foreach)
library(progressr)
library(future)
library(doFuture)
library(polykde)
library(rgl)
library(viridis)

## Hippocampus data

# Load data
load("spokes.rda")

# Data transformation. dirs is an array c(n = 177, r = 168, d = 2 + 1).
# We need to transform it to the c(n, sum(d) + r) format that concatenates
# columns.
r <- ncol(dirs)
d <- rep(2, r)
n <- nrow(dirs)
X <- do.call(cbind, args = lapply(1:r, function(i) dirs[, i, ]))
stopifnot(ncol(X) == (sum(d) + r))
stopifnot(all(abs(rowSums(X^2) - r) < 1e-10))

## Bandwidths for initial estimation and separation, based on the vMF kernels

kernel_sim <- 1
kernel_type_sim <- 2
k_sfp_sim <- 10
h_mrot <- bw_mrot_polysph(X = X, d = d, kernel = kernel_sim, k = k_sfp_sim)
bw0 <- 3 * h_mrot
h_rot <- bw_rot_polysph(X = X, d = d, bw0 = bw0, kernel = kernel_sim,
                        kernel_type = kernel_type_sim, k = k_sfp, iterlim = 200,
                        steptol = 1e-10, gradtol = 1e-10)
h_rot$opt$minimum
mean(h_rot$bw - h_mrot)
mean(h_rot$bw - bw0)
h_rot <- h_rot$bw

## Density of hippocampus shapes

# Log-cross-validated kde evaluated at each point
h <- h_rot
log_kde <- log_cv_kde_polysph(X = X, d = d, h = h, wrt_unif = TRUE,
                              kernel = kernel_sim, kernel_type = kernel_type_sim,
                              k = k_sfp_sim)
sum(duplicated(log_kde))
sum(!is.finite(log_kde))

# Quantiles
qs_log_kde <- quantile(log_kde)

# Split sample into top 25% and bottom 25% density
X_inw <- X[log_kde > qs_log_kde[4], ]
X_out <- X[log_kde < qs_log_kde[2], ]

# Visualize groups
open3d()
mfrow3d(nr = 12, nc = 14, sharedMouse = TRUE)
col <- c("red", "white", "blue")[cut(log_kde, breaks = c(-Inf, qs_log_kde[2],
                                                         qs_log_kde[4], Inf))]
for (i in seq_len(ncol(dirs))) {
  col_i <- col # ids_labs + 1
  plot3d(dirs[, i, ], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
         col = col_i, axes = FALSE, box = FALSE, xlab = "",
         ylab = "", zlab = "")
  spheres3d(0, 0, 0, radius = 1, col = "lightblue", lit = FALSE,
            alpha = 0.5)
}

# Estimate bandwidths for the two groups
h_rot_joi <- bw_rot_polysph(X = rbind(X_inw, X_out), d = d, bw0 = 1.5 * h_rot,
                            kernel = kernel_sim, kernel_type = kernel_type_sim,
                            k = k_sfp_sim, iterlim = 1e3, steptol = 1e-10,
                            gradtol = 1e-10, print.level = 1)
h_rot_joi$opt$minimum
h_rot_inw <- bw_rot_polysph(X = X_inw, d = d, bw0 = 0.8 * h_rot, iterlim = 1e3,
                            kernel = kernel_sim, kernel_type = kernel_type_sim,
                            k = k_sfp_sim, steptol = 1e-10, gradtol = 1e-10,
                            print.level = 1)
h_rot_inw$opt$minimum
h_rot_out <- bw_rot_polysph(X = X_out, d = d, bw0 = 1.5 * h_rot, iterlim = 1e3,
                            kernel = kernel_sim, kernel_type = kernel_type_sim,
                            k = k_sfp_sim, steptol = 1e-10, gradtol = 1e-10,
                            print.level = 1)
h_rot_out$opt$minimum

## Samplers of H_a

# Sampler of the mixture f_a = (1 - a) / 2 * f_1 + (1 + a) / 2 * f_2
r_f_a <- function(n, X1, X2, d, h1, h2, a = 0) {

  stopifnot(a >= 0 && a <= 1)
  n1 <- sum(runif(n) < (1 - a) / 2)
  n2 <- n - n1
  samp_1 <- r_kde_polysph(n = n1, X = X1, d = d, h = h1, kernel = kernel_sim,
                          kernel_type = kernel_type_sim, k = k_sfp_sim)
  samp_2 <- r_kde_polysph(n = n2, X = X2, d = d, h = h2, kernel = kernel_sim,
                          kernel_type = kernel_type_sim, k = k_sfp_sim)
  samp <- rbind(samp_1, samp_2)
  return(samp[sample(n), , drop = FALSE])

}

# Sampler of H_a. Samples from the two groups defined by
# * n1 observations from f_1a = (1 - a) / 2 * f_1 + (1 + a) / 2 * f_2
# * n2 observations from f_2a = (1 + a) / 2 * f_1 + (1 - a) / 2 * f_2
r_hyp_a <- function(n1, n2, X1, X2, d, h1, h2, a = 0) {

  X <- rbind(r_f_a(n = n1, a = a, X1 = X1, X2 = X2, h1 = h1, h2 = h2, d = d),
             r_f_a(n = n2, a = a, X1 = X2, X2 = X1, h1 = h2, h2 = h1, d = d))
  labels <- rep(x = 1:2, times = c(n1, n2))
  ind <- sample(n1 + n2)
  X <- X[ind, , drop = FALSE]
  labels <- labels[ind]
  return(list("X" = X, "labels" = labels))

}

## Simulations

# Define a progress bar
handlers(handler_progress(
  format = "(:spin) [:bar] :percent Iter: :current/:total Rate: :tick_rate iter/sec ETA: :eta Elapsed: :elapsedfull",
  clear = FALSE))

# Monte Carlo setup
save <- TRUE
len_a <- 5
a <- seq(0, 1, l = len_a)
kernel_test <- 3
kernel_type_test <- 2
k_sfp_test <- 100
M <- 1e4
n1 <- 45
n2 <- 45
type <- "jsd"
# type <- "mean"
# type <- "scatter"

# Simulation bandwidths
h1_sim <- h_rot_joi$bw
h2_sim <- h_rot_joi$bw

# Average bandwidths
options(future.rng.onMisuse = "ignore")
doFuture::registerDoFuture()
future::plan(future::multisession(), workers = 8)
progressr::with_progress({
  prog <- progressr::progressor(along = seq_len(M))
  h_rot <- foreach(i = 1:M, .packages = "polykde", .combine = "rbind") %dopar% {

    if (requireNamespace("progressr", quietly = TRUE)) prog()
    set.seed(1e6 + i, kind = "Mersenne-Twister")
    X_med <- r_hyp_a(n1 = n1, n2 = n2, X1 = X_inw, X2 = X_out, d = d,
                     h1 = h1_sim, h2 = h2_sim)$X
    bw_rot_polysph(X = X_med, d = d, bw0 = 1.5 * h1_sim, kernel = kernel_test,
                   kernel_type = kernel_type_test, k = k_sfp_test)$bw

  }
})
progressr::with_progress({
  prog <- progressr::progressor(along = seq_len(M))
  h_cv <- foreach(i = 1:M, .packages = "polykde", .combine = "rbind",
                  .export = "h_rot") %dopar% {

    if (requireNamespace("progressr", quietly = TRUE)) prog()
    set.seed(1e6 + i, kind = "Mersenne-Twister")
    X_med <- r_hyp_a(n1 = n1, n2 = n2, X1 = X_inw, X2 = X_out, d = d,
                     h1 = h1_sim, h2 = h2_sim)$X
    bw_cv_polysph(X = X_med, d = d, method = "L-BFGS-B",
                  bw0 = h_rot[i, ], kernel = kernel_test,
                  kernel_type = kernel_type_test, k = k_sfp_test)$bw

  }
})
(h_med_rot <- apply(h_rot, 2, median))
(h_med_cv <- apply(h_cv, 2, median))
mean(h_med_rot / h_med_cv) # ratio: 0.9263784 for kernel = 3; 0.9218813 for kernel = 1

# Bandwidths for the tests (change)
(c <- 2^(-3:5)[5])
h_test <- c * h_med_rot

# File for saving & loading
file <- paste0("pow_hippo_", type, "_c", sprintf("%.2f", c), "_M1e", log10(M))
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
      samp_hyp_a <- r_hyp_a(n1 = n1, n2 = n2, a = a[j],
                            X1 = X_inw, X2 = X_out, d = d,
                            h1 = h1_sim, h2 = h2_sim)
      tryCatch(
        hom_test_polysph(X = samp_hyp_a$X, d = d, labels = samp_hyp_a$labels,
                         type = type, h = h_test, B = 1, plot_boot = FALSE,
                         M = 1e4, kernel = kernel_test,
                         kernel_type = kernel_type_test, k = k_sfp_test,
                         cv_jsd = 1),
        error = function(e) list(statistic = NA, statistic_perm = NA))

     }
  })

}

# Statistics (M x len_a)
stats <- sapply(tests, function(x) sapply(x, function(y)
  y$statistic))
stats_perm <- sapply(tests, function(x) sapply(x, function(y)
  y$statistic_perm[1]))
colnames(stats) <- colnames(stats_perm) <- a
tryCatch(boxplot(stats), error = function(e) NULL)

# Compute powers
alphas <- seq(0, 1, l = 101)
pow_ws <- t(sapply(alphas, function(alph, two = FALSE) {
  if (!two) {

    q_boot_alph_null_up <- quantile(x = stats_perm[, 1], probs = 1 - alph,
                                    na.rm = TRUE)
    return(rowMeans(t(stats) > q_boot_alph_null_up))

  } else {

    q_boot_alph2_null_up <- quantile(x = stats_perm[, 1], probs = 1 - alph / 2,
                                     na.rm = TRUE)
    q_boot_alph2_null_lo <- quantile(x = stats_perm[, 1], probs = alph / 2,
                                     na.rm = TRUE)
    return(rowMeans((t(stats) < q_boot_alph2_null_lo) |
                      (t(stats) > q_boot_alph2_null_up)))

  }
}))
rownames(pow_ws) <- alphas

# Save RData
if (save) save(list = ls(), file = paste0(file, ".RData"))

## Analysis

# Plots
M <- 1e4
save <- TRUE
for (type in c("jsd", "mean", "scatter")) {
  for (c in 2^(-3:5)) {
    if (type != "jsd" && c != 1) {
      next
    } else {
      file <- paste0("pow_hippo_", type, "_c", sprintf("%.2f", c), "_M1e",
                     log10(M))
      load(paste0(file, ".RData"))
      if (save) pdf(file = paste0(file, ".pdf"), width = 6, height = 6)
      matplot(alphas[-1], pow_ws[-1, ], type = "l", lty = 1,
              col = viridis(len_a), ylab = "Power", xlab = expression(alpha),
              lwd = 2, axes = FALSE, log = "xy", xlim = c(0.01, 1),
              ylim = c(0.01, 1), cex.lab = 1.25)
      axis(1, at = c(0.01, 0.05, 0.10, 0.25, 0.5, 1), cex.axis = 1.25)
      axis(2, at = c(0.01, 0.05, 0.10, 0.25, 0.5, 1), cex.axis = 1.25)
      box()
      abline(h = c(0.01, 1), lty = 3, col = "gray")
      abline(v = c(0.01, 0.05, 0.10, 1), lty = 3, col = "gray")
      abline(a = 0, b = 1, lty = 3)
      lines(alphas, alphas - qnorm(0.975) * sqrt(alphas * (1 - alphas) / M),
            lty = 2)
      lines(alphas, alphas + qnorm(0.975) * sqrt(alphas * (1 - alphas) / M),
            lty = 2)
      points(rep(alphas[-1][5], len_a), pow_ws[-1, ][5, ], col = viridis(len_a),
             pch = 19)
      segments(x0 = 1e-20, y0 = pow_ws[-1, ][5, ], x1 = 0.05,
               y1 = pow_ws[-1, ][5, ], col = viridis(len_a))
      legend("bottomright", legend = paste("a =", round(rev(a), 2)), lwd = 2,
             col = viridis(len_a, direction = -1), bg = "white", cex = 1.5)
      if (save) dev.off()
      pow_ws
    }
  }
}
