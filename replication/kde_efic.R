
# Required libraries
library(polykde)
library(kableExtra)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Visualize kernels

# Plot kernels
pdf("kerns.pdf", width = 6, height = 6)
t <- seq(0, 2, l = 201)
k <- c(1, 2, 5, 10, 20)
col <- rev(viridis::viridis(2 * length(k) + 2))[seq_len(length(k) + 2)]
plot(t, L(t = t, kernel = 2), type = "l", ylab = expression(L(t)),
     xlab = expression(t), xlim = c(0, 2), ylim = c(0, 1), axes = FALSE,
     col = col[length(k) + 2])
for (ki in seq_along(k)) {
  lines(t, L(t = t, kernel = 3, k = k[ki]), type = "l", col = col[ki])
}
lines(t, L(t = t, kernel = 1), type = "l", col = 2)
abline(h = 0:1, lty = 2, col = "gray")
axis(1)
axis(2, at = seq(0, 1, l = 11))
box()
legend("topright", legend = expression("vMF",
                                       "sfp " * (upsilon == 1),
                                       "sfp " * (upsilon == 2),
                                       "sfp " * (upsilon == 5),
                                       "sfp " * (upsilon == 10),
                                       "sfp " * (upsilon == 20),
                                       "Epa"),
       col = c(2, col[1:length(k)], col[length(k) + 2]), lwd = 2, bg = "white")
dev.off()

# Varying r, for fixed d = 1,2
for (d in 1:2) {
  pdf(paste0("effic_r_d_", d, ".pdf"), width = 6, height = 6)
  rr <- 1:20
  k <- c(1, 5, 10, 20)
  col <- viridis::viridis(length(k), direction = -1)
  plot(rr, rep(1, length(rr)), type = "o", pch = 16, col = 4, ylab = "Efficiency",
       xlab = expression(r), ylim = c(0, 1), axes = FALSE)
  lines(rr, sapply(rr, function(r) {
    eff_kern(d = d, r = r, kernel = 1, kernel_type = "prod")
  }), col = 2, lty = 1, type = "o", pch = 16)
  for (ki in seq_along(k)) {
    lines(rr, sapply(rr, function(r) {
      eff_kern(d = d, r = r, kernel = 3, kernel_type = "sph", k = k[ki])
    }), col = col[ki], lty = 1, type = "o", pch = 16)
    lines(rr, sapply(rr, function(r) {
      eff_kern(d = d, r = r, kernel = 3, kernel_type = "prod", k = k[ki])
    }), col = col[ki], lty = 2, type = "o", pch = 16)
  }
  lines(rr, sapply(rr, function(r) {
    eff_kern(d = d, r = r, kernel = 2, kernel_type = "prod")
  }), col = 4, lty = 2, type = "o", pch = 16)
  axis(1, at = seq(1, max(rr), by = 2), labels = seq(1, max(rr), by = 2))
  axis(1, at = seq(2, max(rr), by = 2), labels = seq(2, max(rr), by = 2))
  axis(1, at = rr, labels = rr)
  axis(2, at = seq(0, 1, l = 11))
  box()
  abline(h = 0:1, lty = 2, col = "gray")
  legend(ifelse(d == 2, "topright", "bottomleft"),
         legend = expression("Epa",
                             "vMF",
                             "sfp " * (upsilon == 1),
                             "sfp " * (upsilon == 5),
                             "sfp " * (upsilon == 10),
                             "sfp " * (upsilon == 20),
                             "Sph. sym.",
                             "Product"),
         col = c(4, 2, col[seq_along(k)], 1, 1),
         lty = c(1, 1, rep(1, length(k)), 1, 2), lwd = 2, bg = "white")
  dev.off()
}

# Varying d, fixed r = 1
for (r in 1:2) {
  pdf(paste0("effic_d_r_", r, ".pdf"), width = 6, height = 6)
  dd <- 1:20
  k <- c(1, 5, 10, 20)
  col <- viridis::viridis(length(k), direction = -1)
  plot(dd, rep(1, length(dd)), col = 4, type = "o", pch = 16,
       xlab = expression(d), ylab = "Efficiency", ylim = c(0, 1), axes = FALSE)
  lines(dd, sapply(dd, function(d) eff_kern(d = d, r = r, kernel = 1,
                                            kernel_type = "prod")),
        col = 2, type = "o", pch = 16)
  for (ki in seq_along(k)) {
    lines(dd, sapply(dd, function(d) eff_kern(d = d, r = r, kernel = 3,
                                              kernel_type = "prod", k = k[ki])),
          col = col[ki], type = "o", pch = 16)
  }

  axis(1, at = seq(1, max(dd), by = 2), labels = seq(1, max(dd), by = 2))
  axis(1, at = seq(2, max(dd), by = 2), labels = seq(2, max(dd), by = 2))
  axis(1, at = dd, labels = dd); axis(2, at = seq(0, 1, l = 11)); box()
  abline(h = 0:1, lty = 2, col = "gray")
  legend("bottomleft", legend = expression("Epa",
                                           "vMF",
                                           "sfp " * (upsilon == 1),
                                           "sfp " * (upsilon == 5),
                                           "sfp " * (upsilon == 10),
                                           "sfp " * (upsilon == 20)),
         col = c(4, 2, col), lwd = 2, bg = "white")
  dev.off()
}

## Efficiency tables

dd <- c(1, 2, 3, 5, 10)
rr <- c(1, 2, 3, 5, 10)
eff_vmf <- c(sapply(rr, function(r) sapply(dd, function(d)
  eff_kern(d = d, r = r, kernel = "1", kernel_type = "sph"))))
eff_sfp_S_1 <- c(sapply(rr, function(r) sapply(dd, function(d)
  eff_kern(d = d, r = r, k = 1, kernel = "3", kernel_type = "sph"))))
eff_sfp_S_10 <- c(sapply(rr, function(r) sapply(dd, function(d)
  eff_kern(d = d, r = r, k = 10, kernel = "3", kernel_type = "sph"))))
eff_sfp_S_100 <- c(sapply(rr, function(r) sapply(dd, function(d)
  eff_kern(d = d, r = r, k = 100, kernel = "3", kernel_type = "sph"))))
eff_epa_P <- c(sapply(rr, function(r) sapply(dd, function(d)
  eff_kern(d = d, r = r, kernel = "2", kernel_type = "prod"))))
eff_sfp_P_1 <- c(sapply(rr, function(r) sapply(dd, function(d)
  eff_kern(d = d, r = r, k = 1, kernel = "3", kernel_type = "prod"))))
eff_sfp_P_10 <- c(sapply(rr, function(r) sapply(dd, function(d)
  eff_kern(d = d, r = r, k = 10, kernel = "3", kernel_type = "prod"))))
eff_sfp_P_100 <- c(sapply(rr, function(r) sapply(dd, function(d)
  eff_kern(d = d, r = r, k = 100, kernel = "3", kernel_type = "prod"))))

# Table
eff <- 100 * data.frame(eff_vmf,
                        eff_sfp_S_1, eff_sfp_S_10, eff_sfp_S_100,
                        eff_epa_P,
                        eff_sfp_P_1, eff_sfp_P_10, eff_sfp_P_100)
eff <- round(eff, 2)
eff <- cbind(expand.grid("d" = dd, "r" = rr)[, 2:1], eff)
knitr::kable(eff, format = "markdown")

# LaTeX table
kbl(eff, booktabs = TRUE, format = "latex")
