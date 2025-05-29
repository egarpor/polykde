
# This script produces Figure 4 in the paper.

## Functions

# Required libraries
library(foreach)
library(future)
library(progressr)
library(parallel)
library(doFuture)
library(doRNG)
library(polykde)
library(rotasym)
library(viridis)

# Define a progress bar
handlers(handler_progress(
  format = "(:spin) [:bar] :percent Iter: :current/:total Rate: :tick_rate iter/sec ETA: :eta Elapsed: :elapsedfull",
  clear = FALSE))

# Load data
data("hippocampus", package = "polykde")
attach(hippocampus)

# Data transformation. dirs is an array c(n = 177, r = 168, d = 2 + 1).
# We need to transform it to the c(n, sum(d) + r) format that concatenates
# columns.
r <- ncol(dirs)
d <- rep(2, r)
n <- nrow(dirs)
X <- do.call(cbind, args = lapply(1:r, function(i) dirs[, i, ]))
ind <- cumsum(c(1, d + 1))
stopifnot(ncol(X) == (sum(d) + r))
stopifnot(all(abs(rowSums(X^2) - r) < 1e-10))

## Homogeneity tests

## All the spokes together

# ROT bandwidths
k_sfp <- 100
h_mrot_jsd <- bw_mrot_polysph(X = X, d = d, kernel = 3, k = k_sfp)
h_rot_jsd <- bw_rot_polysph(X = X, d = d, bw0 = 2^(-1:2) %o% h_mrot_jsd,
                            kernel = 3, kernel_type = 2, k = k_sfp,
                            iterlim = 500, steptol = 1e-10, gradtol = 1e-10)

# Common parameters
B <- 5e3
ks <- 2^seq(-3, 5, by = 0.5)
cores <- 6
doFuture::registerDoFuture()
future::plan(future::multisession(), workers = cores)

# JSD
progressr::with_progress({
  prog <- progressr::progressor(along = seq_along(ks))
  test_jsd <- foreach(i = seq_along(ks), .combine = rbind,
                      .inorder = TRUE, .packages = c("polykde")) %dorng% {

    # Run test
    set.seed(42, kind = "Mersenne-Twister")
    res <- tryCatch(
      hom_test_polysph(X = X, d = d, labels = ids_labs, type = "jsd",
                       h = ks[i] * h_rot_jsd$bw, kernel = 3, kernel_type = 2,
                       k = k_sfp, B = B, cv_jsd = 1, plot_boot = FALSE),
      error = function(e) {print(e); NA})

    # Signal progress
    if (requireNamespace("progressr", quietly = TRUE)) {

      prog()

    }
    return(res)
  }
})

# Means
test_mean <- hom_test_polysph(X = X, d = d, labels = ids_labs, type = "mean",
                              B = B, plot_boot = FALSE)
test_mean

# Scatter
test_scat <- hom_test_polysph(X = X, d = d, labels = ids_labs, type = "scatter",
                              B = B, plot_boot = FALSE)
test_scat

# Save & load tests
save(list = ls(), file = "tests.RData")
# load("tests.RData")

# JSD -- Figure 4
pdf("pval_jsd.pdf", width = 10, height = 5)
plot(ks, unlist(test_jsd[, 2]), type = "o", pch = 16, ylim = c(0, 1),
     xlab = expression(c),
     ylab = expression(p * "-value"), axes = FALSE)
abline(v = 1, lty = 2)
axis(1, at = ks, labels = round(ks, 2)); axis(2, at = seq(0, 1, l = 11)); box()
abline(h = c(0.01, 0.05, 0.10), lty = 3)
dev.off()

## Each spoke separately

# Common parameters
B <- 5e3
ks <- 2^seq(-3, 5, by = 0.5)
cores <- 6
doFuture::registerDoFuture()
future::plan(future::multisession(), workers = cores)

# JSD
progressr::with_progress({
  prog <- progressr::progressor(along = 1:r)
  pval_jsd <- foreach(j = 1:r, .combine = rbind, .inorder = TRUE,
                      .packages = c("polykde")) %dorng% {

    # Run tests
    ind_j <- ind[j]:(ind[j + 1] - 1)
    res <- numeric(length(ks))
    for (i in seq_along(ks)) {
      set.seed(42)
      res[i] <- tryCatch(
        hom_test_polysph(X = X[, ind_j], d = d[j], labels = ids_labs,
                         type = "jsd", kernel = 3, kernel_type = 2, k = k_sfp,
                         h = ks[i] * h_mrot_jsd[j], B = B, cv_jsd = 1,
                         plot_boot = FALSE)$p.value,
        error = function(e) {print(e); NA})
    }

    # Signal progress
    if (requireNamespace("progressr", quietly = TRUE)) {

      prog()

    }
    return(res)

  }
})

# Save & load tests
save(list = ls(), file = "tests.RData")
# load("tests.RData")

# p-values traces
pdf("pval_each_jsd.pdf", width = 6, height = 6)
pval_jsd_fdr <- apply(pval_jsd, 2, p.adjust, method = "BY")
matplot(1:r, pval_jsd_fdr, type = "l", lty = 1, ylim = c(0, 1),
        xlab = "Spoke index", ylab = "FDR-corrected p-value", axes = FALSE,
        col = viridis(length(ks)))
axis(1, at = seq(0, r, by = 20)); axis(2, at = seq(0, 1, l = 11)); box()
abline(h = c(0.01, 0.05, 0.10), lty = 3)
dev.off()

# Mean and scatter tests
pval_mean <- pval_scat <- rep(NA, r)
for (j in seq_len(r)) {

  # Store p-values
  ind_j <- ind[j]:(ind[j + 1] - 1)
  pval_mean[j] <- hom_test_polysph(X = X[, ind_j], d = d[j], labels = ids_labs,
                                   type = "mean", B = B,
                                   plot_boot = FALSE)$p.value
  pval_scat[j] <- hom_test_polysph(X = X[, ind_j], d = d[j], labels = ids_labs,
                                   type = "scatter", B = B,
                                   plot_boot = FALSE)$p.value

  # Progress
  cat("\nj =", j, "\n")
  plot(pval_mean, type = "l", ylim = c(0, 1))
  lines(pval_scat, col = 2)
  abline(h = c(0.01, 0.05, 0.10), lty = 3)

}

# Save & load tests
save(list = ls(), file = "tests.RData")
# load("tests.RData")

# p-values traces
pdf("pval_mean_scat_fdr.pdf", width = 6, height = 6)
pval_mean_fdr <- p.adjust(pval_mean, method = "BY")
pval_scat_fdr <- p.adjust(pval_scat, method = "BY")
plot(pval_mean_fdr, type = "l", ylim = c(0, 1), xlab = "Spoke index",
     ylab = "FDR-corrected p-value", axes = FALSE)
lines(pval_scat_fdr, col = 2)
axis(1, at = seq(0, r, by = 20)); axis(2, at = seq(0, 1, l = 11)); box()
abline(h = c(0.01, 0.05, 0.10), lty = 3)
legend("topleft", lwd = 2, col = 1:2, bg = "white",
       legend = c("Equal location vectors", "Equal scatter matrices"))
dev.off()
