
# Required libraries
library(polykde)
library(foreach)
library(latex2exp)
library(doParallel)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Simulation

# Setting
r <- 2
d <- 2
dd <- rep(d, r)
stopifnot(all(dd == d))
mu <- rep(c(1, rep(0, d)), r)
kappa <- rep(5, r)
x <- rbind(mu)
kernel <- 1
# x <- r_vmf_polysph(n = 1, d = dd, mu = mu, kappa = kappa)

# Variance elements
f_x <- drop(kde_polysph(x = x, X = rbind(mu), d = dd, h = 1 / sqrt(kappa),
                        wrt_unif = FALSE))
vd <- v_d(d = d, kernel = kernel)

# Bias elements
bd <- b_d(d = d, kernel = kernel)
Hf <- function(x) grad_hess_kde_polysph(x = rbind(x), X = rbind(mu), d = dd,
                                        h = 1 / sqrt(kappa))$hess[1, , ]
Hf_x <- Hf(x)
nabla2_samp <- r_vmf_polysph(n = 1e5, d = dd, mu = mu, kappa = rep(0, r))
nabla2_x <- apply(nabla2_samp, 1, function(xi) sum(diag(Hf(x = xi)))^2)
(R <- rotasym::w_p(p = d + 1)^r * mean(nabla2_x))
(C <- (vd^r / bd^2 * (d * r) / (4 * R))^(1 / (d * r + 4)))

# # When r = 1, check R with:
# d^2 * DirStats::R_Psi_mixvmf(q = d, mu = rbind(mu), kappa = kappa, p = 1)

# Monte Carlo
M <- 1e4
(n <- 2^(7:17))
delta <- c(-2:2, 4)
kde <- array(NA, dim = c(M, length(n), length(delta)))
h <- matrix(NA, nrow = length(n), ncol = length(delta))
cl <- makeCluster(7)
registerDoParallel(cl)
for (j in seq_along(n)) {

  # Progress
  message("\nn = ", n[j])
  pb <- txtProgressBar(style = 3)

  for (k in seq_along(delta)) {

    # Bandwidths
    h[j, k] <- C * n[j]^(-1 / (4 + d * r + delta[k]))
    hh <- rep(h[j, k], r)

    # # Sequential Monte Carlo
    # for (i in seq_len(M)) {
    #
    #   set.seed(i, kind = "Mersenne-Twister")
    #   samp <- r_vmf_polysph(n = n[j], d = dd, mu = mu, kappa = kappa)
    #   kde[i, j, k] <- tryCatch(kde_polysph(x = x, X = samp, d = dd, h = hh,
    #                                        wrt_unif = FALSE, kernel = kernel),
    #                            error = function(e) NA)
    #   setTxtProgressBar(pb, value = ((k - 1) * M + i) / (M * length(delta)))
    #
    # }

    # Parallel Monte Carlo
    kde[, j, k] <- foreach::foreach(i = 1:M, .combine = "c",
                                    .inorder = TRUE, .multicombine = TRUE,
                                    .maxcombine = 1e3, .packages = "polykde"
                                    ) %dopar% {

      set.seed(i, kind = "Mersenne-Twister")
      samp <- r_vmf_polysph(n = n[j], d = dd, mu = mu, kappa = kappa)
      tryCatch(kde_polysph(x = x, X = samp, d = dd, h = hh,
                           wrt_unif = FALSE, kernel = kernel),
               error = function(e) NA)

    }
    setTxtProgressBar(pb, value = k / length(delta))

  }

}
stopCluster(cl)

# Asymptotic expectation, variance, and rate
asymp_exp <- f_x + bd * sum(diag(Hf_x)) * h^2
asymp_var <- (vd^r * f_x) / (n * (h^d)^r)
# bd_vmf_h <- (1 / d) * lambda_vmf_h(d = d, h = h, bias = TRUE) /
#   lambda_vmf_h(d = d, h = h, bias = FALSE)
# vd_vmf_h <- lambda_vmf_h(d = d, h = h, squared = TRUE) /
#   lambda_vmf_h(d = d, h = h, squared = FALSE)^2
# asymp_exp2 <- f_x + bd_vmf_h * sum(diag(Hf_x)) * h^2
# asymp_var2 <- (vd_vmf_h^r * f_x) / (n * (h^d)^r)
rate_n <- 1 / sqrt(asymp_var)

# Exact expectation and variance
exact_exp <- colMeans(kde)
exact_var <- apply(kde, 2:3, var)

# Zn statistics
zn_1 <- zn_2 <- array(NA, dim = c(M, length(n), length(delta)))
for (i in seq_len(M) {

  zn_1[i, , ] <- rate_n * (kde[i, , ] - exact_exp)    # Central
  zn_2[i, , ] <- rate_n * (kde[i, , ] - asymp_exp)   # Expansion

}

# Save all
save(list = ls(), file = paste0("kde-norm-kernel", kernel,
                                "-r", r, "-d", d, "-k", kappa[1],
                                "-M", log10(M), ".RData"))

## Analysis

# (r, d, kernel) = (1, 1, *): Perfect match for bias/var for delta = -2:2
# (r, d, kernel) = (1, 2, *): Perfect match for bias/var for delta = -1:2
# (r, d, kernel) = (2, 1, *): NO match for bias for delta = -2:2,
#                             but match for var for delta = -2:0 (larger n's).
#                             The kernel has an effect.
# (r, d, kernel) = (2, 2, *): NO match for bias for delta = -2:2,
#                             but match for var for delta = -2:0 (larger n's).
#                             The kernel has an effect.
# Problem: ^r must be removed on the bias moment!

# Load data
r <- 2
d <- 2
kernel <- 1
k <- 5
M <- 1e5
load(dir(pattern = paste0("kde-norm-kernel", kernel,
                          "-r", r, "-d", d, "-k", k, "-M", log10(M), "-*"))[1])

# # Subset to delta = -2:2, 4
# delta_small <- c(-2:2, 4)
# ind <- delta %in% delta_small
# zn_1 <- zn_1[ , , ind]
# zn_2 <- zn_2[ , , ind]
# exact_exp <- exact_exp[, ind]
# asymp_exp <- asymp_exp[, ind]
# exact_var <- exact_var[, ind]
# asymp_var <- asymp_var[, ind]
# h <- h[, ind]
# delta <- delta_small

# Plots of expectation, log(bias), and variance
col <- viridis::plasma(length(delta))
pdf("kde_exp.pdf", width = 6, height = 6)
matplot(log2(n), exact_exp, type = "l", lty = 1, xlab = "\u2113",
        ylab = "Expectation", col = col)
matlines(log2(n), asymp_exp, type = "l", lty = 2, col = col)
abline(h = f_x)
legend("bottomright", legend = TeX(paste("\\delta =", delta)), lwd = 2,
       col = col, bg = "white")
legend("topleft", legend = c("Exact", "Asymptotic"), lty = c(1, 2),
       bg = "white")
dev.off()
pdf("kde_logbias.pdf", width = 6, height = 6)
matplot(log2(n), log10(f_x - exact_exp), type = "l", lty = 1, xlab = "\u2113",
        ylab = TeX("$\\log_{10}(Bias)$"), col = col, ylim = c(-1.7, 0))
matlines(log2(n), log10(-bd * sum(diag(Hf_x)) * h^2),
         type = "l", lty = 2, col = col)
dev.off()
pdf("kde_sd.pdf", width = 6, height = 6)
matplot(log2(n), sqrt(exact_var), type = "l", lty = 1, xlab = "\u2113",
        ylab = "Standard deviation", col = col)
matlines(log2(n), sqrt(asymp_var), type = "l", lty = 2, col = col)
dev.off()

# Evolution of mean, variance, and KS statistic
mean_1 <- apply(zn_1, 2:3, mean)
sd_1 <- apply(zn_1, 2:3, sd)
mean_2 <- apply(zn_2, 2:3, mean)
sd_2 <- apply(zn_2, 2:3, sd)
stat_1 <- apply(zn_1, c(2, 3), function(x) ks.test(x, "pnorm")$statistic)
stat_2 <- apply(zn_2, c(2, 3), function(x) ks.test(x, "pnorm")$statistic)
colnames(mean_1) <- colnames(mean_2) <- colnames(sd_1) <- colnames(sd_2) <-
  colnames(stat_1) <- colnames(stat_2) <- delta

# Plots of the evolution of summary statistics
col <- viridis::plasma(length(delta))
pdf("zn2_logexp.pdf", width = 6, height = 6)
matplot(log2(n), log10(mean_2), type = "l", lty = 1, xlab = "\u2113",
        ylab = TeX("$\\log_{10}(Expectation)$"), col = col, axes = FALSE)
axis(1, at = 7:17); axis(2); box()
abline(h = 0, lty = 3)
legend("bottomleft", legend = TeX(paste("\\delta =", delta)), lwd = 2,
       col = col, bg = "white")
dev.off()
pdf("zn2_sd.pdf", width = 6, height = 6)
matplot(log2(n), sd_2, type = "l", lty = 1, xlab = "\u2113",
        ylab = "Standard deviation", col = col, ylim = c(0.5, 1), axes = FALSE)
abline(h = 1, lty = 3)
axis(1, at = 7:17); axis(2); box()
dev.off()
pdf("zn2_pval.pdf", width = 6, height = 6)
matplot(log2(n), stat_1, type = "l", lty = 2, xlab = "\u2113",
        ylab = "Kolmogorov-Smirnov statistic", col = col, ylim = c(0, 1),
        axes = FALSE)
matlines(log2(n), stat_2, type = "l", lty = 1, col = col)
legend("right", legend = TeX(c("$Z_{n,\\delta}^{(1)}$",
                               "$Z_{n,\\delta}^{(2)}$")),
       lty = c(1, 2), bg = "white")
axis(1, at = 7:17); axis(2); box()
dev.off()

# Density plots
col <- viridis::viridis(length(n))
for (k in seq_along(delta)) {

  # Zn1
  pdf(paste0("zn1_delta_", delta[k], ".pdf"), width = 6, height = 6)
  curve(dnorm(x), lwd = 2, n = 1e3, from = -5, to = 5, xlim = c(-4.5, 4.5),
        ylim = c(0, 0.55), xlab = expression(Z[list(n, delta)]^{(1)}),
        ylab = "Density", main = "")
  for (j in seq_along(n)) {
    lines(density(zn_1[, j, k], bw = "nrd", n = 1024, from = -5, to = 5),
          col = col[j])
  }
  curve(dnorm(x), lwd = 2, n = 1e3, from = -5, to = 5, add = TRUE)
  if (k == 1) legend("topleft", legend = paste("n =", n), lwd = 2, col = col,
                     bg = "white")
  dev.off()

  # Zn2
  pdf(paste0("zn2_delta_", delta[k], ".pdf"), width = 6, height = 6)
  curve(dnorm(x), lwd = 2, n = 200, from = -5, to = 5, xlim = c(-4.5, 4.5),
        ylim = c(0, 0.55), xlab = expression(Z[list(n, delta)]^{(2)}),
        ylab = "Density", main = "")
  for (j in seq_along(n)) {
    lines(density(zn_2[, j, k], bw = "nrd", n = 1024, from = -5, to = 5),
          col = col[j])
  }
  curve(dnorm(x), lwd = 2, n = 1e3, from = -5, to = 5, add = TRUE)
  if (k == 1) legend("topleft", legend = paste("n =", n), lwd = 2, col = col,
                     bg = "white")
  dev.off()

}
