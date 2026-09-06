
# Required libraries
library(polykde)
library(movMF)
library(DirStats)
library(mlbench)
library(compositions)
library(parallel)
library(ggplot2)
stopifnot(packageVersion("polykde") >= "1.2.1")
stopifnot(packageVersion("DirStats") >= "1.0.0")

## Settings
{

# Number of stratified train/test splits and movMF random starts per fit
R <- 100
nruns <- 25

# Cores for parallelization
n_cores <- 12

# K grid: dense up to 10, then every third K up to the dataset's frontier
# K_top = 0.7 * (largest class) / 3, the EM-feasibility limit of a class fit
# (a component needs a few observations). The stricter one-parameter-per-
# observation cap n / (d + 2) truncates movMF below its accuracy optimum, so
# it is not used.
k_grid <- function(K_top) {

  if (K_top > 10) {

    sort(unique(c(1:10, seq(13, K_top, by = 3), K_top)))

  } else {

    seq_len(K_top)

  }

}

}

## Sphere embeddings
{

# Square-root map of a composition to the sphere positive orthant
sqrt_map <- function(X) {

  X <- as.matrix(X)
  sqrt(X / rowSums(X))

}

# L2 normalization: columns centered/scaled with the training statistics
# ctr and scl to avoid leakage
l2_map <- function(X, ctr, scl) {

  X <- sweep(sweep(as.matrix(X), 2, ctr, "-"), 2, scl, "/")
  X / sqrt(rowSums(X^2))

}

}

## Bandwidth selectors
{

# CV selector with vMF kernel and arcsinh trick
bw_cv <- function(X, d) {

  # suppressWarnings to silence the 1-d Nelder-Mead warning
  suppressWarnings(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV", exact_vmf = TRUE,
                  arcsinh = TRUE, spline = TRUE)$bw)

}

# ROT selector
bw_rot <- function(X, d) {

  bw_rot_polysph(X = X, d = d, kernel = 1)$bw

}

# EMI selector (depends on Monte Carlo, so seed is fixed)
bw_emi <- function(X, fit, seed) {

  set.seed(seed)
  bw_dir_emi(data = X, fit_mix = fit_to_mix(fit), optim = TRUE,
             plot_it = FALSE)$h_opt

}

# AMI selector
bw_ami <- function(X, fit) {

  bw_dir_ami(data = X, fit_mix = fit_to_mix(fit))

}

}

## Mixture fitting
{

# Convert a movMF fit (theta = kappa * mu, alpha) to the DirStats fit_mix format
fit_to_mix <- function(fit) {

  kap <- sqrt(rowSums(fit$theta^2))
  list(best_fit = list(mu_hat = fit$theta / kap, kappa_hat = kap,
                       p_hat = fit$alpha))

}


# movMF fit with k components; the single place where an EM failure is
# absorbed, into NULL (counted downstream)
fit_movmf <- function(X, k) {

  tryCatch(movMF(X, k = k, nruns = nruns), error = function(e) NULL)

}

}

## Classifiers
{

# kde log-density matrix (n_new x n_class): one common bandwidth h, or
# per-class selection via bwfun
kda_logdens <- function(Ztr, ytr, Znew, d, h = NULL, bwfun = NULL) {

  sapply(levels(ytr), function(l) {

    Xl <- Ztr[ytr == l, , drop = FALSE]
    hl <- if (is.null(h)) bwfun(Xl, d) else h
    kde_polysph(x = Znew, X = Xl, d = d, h = hl, kernel = 1, log = TRUE)

  })

}

# Assign each row to the class maximizing log-density plus log-prior
classify <- function(ld, logprior, cls) {

  cls[max.col(sweep(ld, 2, logprior, "+"), ties.method = "first")]

}

}

## Data splitting and embedding
{

# Remove coincident points before splitting: ties make the CV loss unbounded
# below (the leave-one-out density diverges as h -> 0). Only LetterRecognition
# is affected, via its integer features.
dedup_ds <- function(ds) {

  keep <- !duplicated(round(ds$X, 10))
  ds$X <- ds$X[keep, , drop = FALSE]
  ds$y <- droplevels(ds$y[keep])
  ds

}

# Stratified 70/30 split and leak-free embedding (l2 uses the training
# center/scale)
split_embed <- function(ds, seed, prop = 0.7) {

  X <- ds$X
  y <- ds$y
  set.seed(seed)
  tr <- unlist(lapply(levels(y), function(l) {

    i <- which(y == l)
    sample(i, round(prop * length(i)))

  }))
  if (ds$embed == "sqrt") {

    Ztr <- sqrt_map(X[tr, , drop = FALSE])
    Zte <- sqrt_map(X[-tr, , drop = FALSE])

  } else {

    ctr <- colMeans(X[tr, , drop = FALSE])
    scl <- apply(X[tr, , drop = FALSE], 2, sd)
    Ztr <- l2_map(X[tr, , drop = FALSE], ctr, scl)
    Zte <- l2_map(X[-tr, , drop = FALSE], ctr, scl)

  }
  list(Ztr = Ztr, ytr = droplevels(y[tr]), Zte = Zte, yte = y[-tr])

}

# One split: test accuracy of kde-CV (common and per-class bandwidth) and
# kde-ROT and, for each K in K_grid, of kde-EMI/kde-AMI with the pooled
# K-component mixture as reference and of the movMF classifier with per-class
# K-component fits (each class capped at its frontier). Every quantity is
# caught individually: a failure yields NA for it only, counted by run_ds().
sweep_split <- function(ds, seed, K_grid) {

  sp <- split_embed(ds, seed)
  d <- ds$d
  cls <- levels(sp$ytr)
  logprior <- log(as.numeric(table(sp$ytr)) / length(sp$ytr))
  acc <- function(ld) mean(classify(ld, logprior, cls) == sp$yte)
  acc_h <- function(h) acc(kda_logdens(sp$Ztr, sp$ytr, sp$Zte, d, h = h))
  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)

  # Mixture-free kdes: bandwidths (kept for later analysis) and accuracies
  h_cv <- safe(bw_cv(sp$Ztr, d))
  h_rot <- safe(bw_rot(sp$Ztr, d))
  flat <- c(CV_com = safe(acc_h(h_cv)),
            CV_cls = safe(acc(kda_logdens(sp$Ztr, sp$ytr, sp$Zte, d,
                                          bwfun = bw_cv))),
            ROT_com = safe(acc_h(h_rot)),
            ROT_cls = safe(acc(kda_logdens(sp$Ztr, sp$ytr, sp$Zte, d,
                                           bwfun = bw_rot))))

  # Pooled K-component mixtures, passed directly to EMI and AMI
  pooled <- lapply(K_grid, function(K) fit_movmf(sp$Ztr, k = K))
  h_emi <- sapply(seq_along(K_grid), function(i)
    safe(bw_emi(sp$Ztr, pooled[[i]], seed = seed * 1000L + K_grid[i])))
  h_ami <- sapply(pooled, function(f) safe(bw_ami(sp$Ztr, f)))
  emi <- sapply(h_emi, function(h) safe(acc_h(h)))
  ami <- sapply(h_ami, function(h) safe(acc_h(h)))

  # movMF classifier: per-class fits at the grid k's up to min(K, class
  # frontier), each distinct k fitted once per class; a class whose fit at
  # some k fails uses its largest successful k below it (dropped fits counted)
  Xc <- lapply(cls, function(l) sp$Ztr[sp$ytr == l, , drop = FALSE])
  cap <- sapply(Xc, function(X) max(1L, floor(nrow(X) / 3)))
  class_fits <- Map(function(X, kc) {

    ks <- unique(pmin(K_grid, kc))
    fits <- lapply(ks, function(k) fit_movmf(X, k = k))
    ok <- !sapply(fits, is.null)
    list(fits = fits[ok], k = ks[ok], n_fail = sum(!ok))

  }, Xc, cap)
  mv_acc <- function(k_of) safe(acc(sapply(class_fits, function(cf) {

    f <- cf$fits[[k_of(cf)]]
    dmovMF(sp$Zte, theta = f$theta, alpha = f$alpha, log = TRUE)

  })))
  mv <- sapply(K_grid, function(K)
    mv_acc(function(cf) which.max(cf$k * (cf$k <= K))))

  # movMF with K chosen by BIC (text-only, for the pooled vs per-class
  # comparison): once on the pooled fits, then per-class fits at that K, or
  # on each class's own fits
  ok <- !sapply(pooled, is.null)
  K_bic <- safe(K_grid[ok][which.min(sapply(pooled[ok], BIC))])
  mv_bic <- c(
    mv_bic_com = mv_acc(function(cf) which.max(cf$k * (cf$k <= K_bic))),
    mv_bic_cls = mv_acc(function(cf) which.min(sapply(cf$fits, BIC))))

  # Result: named accuracies, the selected bandwidths, and the number of
  # failed EM fits
  c(flat, setNames(emi, paste0("emi", K_grid)),
    setNames(ami, paste0("ami", K_grid)), setNames(mv, paste0("mv", K_grid)),
    mv_bic, h_CV = h_cv, h_ROT = h_rot,
    setNames(h_emi, paste0("hemi", K_grid)),
    setNames(h_ami, paste0("hami", K_grid)),
    n_fail = sum(!ok) + sum(sapply(class_fits, `[[`, "n_fail")))

}

}

## Datasets
{

# Llobregat-basin river hydrochemistry (compositions::Hydrochem): 14 chemical
# parts mapped to S^13 by the sqrt map, 4 river classes
load_hydrochem <- function() {

  data("Hydrochem", package = "compositions")
  parts <- c("H", "Na", "K", "Mg", "Ca", "Sr", "Ba", "NH4", "Cl", "NO3",
             "PO4", "SO4", "HCO3", "TOC")
  list(X = as.matrix(Hydrochem[, parts]), y = factor(Hydrochem$River),
       d = length(parts) - 1, embed = "sqrt", name = "Hydrochem")

}

# Deterding vowel recognition (mlbench::Vowel): the 9 LPC features V2:V10 (V1 is
# a speaker indicator) map to S^8, with 11 vowel classes
load_vowel <- function() {

  data("Vowel", package = "mlbench")
  feat <- paste0("V", 2:10)
  list(X = as.matrix(Vowel[, feat]), y = factor(Vowel$Class),
       d = length(feat) - 1, embed = "l2", name = "Vowel")

}

# Letter recognition (mlbench::LetterRecognition): 16 features map to S^15, with
# 26 classes; subsampled to n_sub rows for tractable repeated splits
load_letter <- function(n_sub = 5000, seed = 1) {

  data("LetterRecognition", package = "mlbench")
  set.seed(seed)
  i <- sample(nrow(LetterRecognition), min(n_sub, nrow(LetterRecognition)))
  feat <- setdiff(names(LetterRecognition), "lettr")
  list(X = as.matrix(LetterRecognition[i, feat]),
       y = droplevels(LetterRecognition$lettr[i]),
       d = length(feat) - 1, embed = "l2", name = "LetterRecognition")

}

}

## Experiments
{

# Run one dataset: de-duplicate, set the K grid from the largest class, sweep
# the splits in parallel, report the failure counts
run_ds <- function(ds) {

  ds <- dedup_ds(ds)
  K_grid <- k_grid(floor(0.7 * max(table(ds$y)) / 3))
  cat(ds$name, "n =", nrow(ds$X), " K_top =", max(K_grid), "\n",
      file = stderr())
  res <- mclapply(seq_len(R), function(s)
    tryCatch(sweep_split(ds, seed = s, K_grid = K_grid),
             error = function(e) NULL), mc.cores = n_cores)

  # A dead worker gives a non-numeric entry: keep it as an all-NA split
  bad <- !sapply(res, is.numeric)
  if (any(bad)) {

    stopifnot(!all(bad))
    res[bad] <- list(res[[which(!bad)[1]]] * NA)
    warning(sum(bad), " of ", R, " splits failed for ", ds$name, ".")

  }
  A <- do.call(rbind, res)
  nf <- colSums(is.na(A[, colnames(A) != "n_fail", drop = FALSE]))
  if (any(nf > 0)) {

    cat("  NA/", R, ": ", paste(names(nf)[nf > 0], nf[nf > 0], sep = "=",
                                collapse = ", "), "\n", sep = "",
        file = stderr())

  }
  cat("  dropped EM fits:", sum(A[, "n_fail"], na.rm = TRUE), "\n",
      file = stderr())
  list(name = ds$name, n = nrow(ds$X), d = ds$d, cl = nlevels(ds$y),
       K_grid = K_grid, acc = A)

}

# Each dataset is cached as soon as it finishes (cv-hd-<dataset>-R<R>.RData),
# so an interrupted run resumes and earlier R's stay reloadable
datasets <- list(load_hydrochem(), load_vowel(), load_letter())
res <- lapply(datasets, function(ds) {

  f <- sprintf("cv-hd-%s-R%d.RData", ds$name, R)
  if (!file.exists(f)) {

    z <- run_ds(ds)
    save(z, nruns, file = f)

  }
  load(f)
  z

})

}

## Summary
{

# Mean over the splits, NA unless at least half of them are available
mean_ok <- function(v) if (sum(!is.na(v)) >= R / 2) mean(v, na.rm = TRUE) else NA

# Numbers quoted in the text: mean (sd) accuracy, in %, of the mixture-free
# kdes and the maximum over K of each mixture curve
ms <- function(v) sprintf("%.1f (%.1f)", mean(v, na.rm = TRUE),
                          sd(v, na.rm = TRUE))
for (z in res) {

  A <- 100 * z$acc
  cat(z$name, ": kde-CV", ms(A[, "CV_com"]), " cls", ms(A[, "CV_cls"]),
      " kde-ROT", ms(A[, "ROT_com"]), " cls", ms(A[, "ROT_cls"]),
      " movMF-BIC", ms(A[, "mv_bic_com"]), " cls", ms(A[, "mv_bic_cls"]),
      "\n")
  for (pre in c("emi", "ami", "mv")) {

    m <- apply(A[, paste0(pre, z$K_grid), drop = FALSE], 2, mean_ok)
    cat("  ", pre, ": max ", sprintf("%.1f", max(m, na.rm = TRUE)), " at K = ",
        z$K_grid[which.max(m)], "\n", sep = "")

  }

}

}

## Final figure
{

# One row per (dataset, K, method) with the mean and sd over splits; the
# mixture-free methods are replicated along K
df <- do.call(rbind, lapply(res, function(z) {

  stat <- function(cols, method) {

    A <- z$acc[, cols, drop = FALSE]
    acc <- apply(A, 2, mean_ok)
    data.frame(dataset = z$name, K = z$K_grid, method = method, acc = acc,
               sd = ifelse(is.na(acc), NA, apply(A, 2, sd, na.rm = TRUE)),
               row.names = NULL)

  }
  nK <- length(z$K_grid)
  rbind(stat(rep("CV_com", nK), "kde-CV"), stat(rep("ROT_com", nK), "kde-ROT"),
        stat(paste0("mv", z$K_grid), "movMF"))

}))
df$dataset <- factor(df$dataset, levels = sapply(res, `[[`, "name"))
df$method <- factor(df$method, levels = c("kde-CV", "kde-ROT", "movMF"))
col <- c("kde-CV" = "#000000", "kde-ROT" = "#009E73", "movMF" = "#CC79A7")
lty <- c("kde-CV" = "dashed", "kde-ROT" = "dashed", "movMF" = "solid")
gg <- ggplot(df, aes(x = K, y = acc, color = method, linetype = method)) +
  geom_ribbon(aes(ymin = acc - sd, ymax = acc + sd, fill = method),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(data = df[df$method == "movMF", ], size = 0.9) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  scale_linetype_manual(values = lty) +
  scale_x_continuous(transform = "log10",
                     breaks = c(1, 2, 3, 5, 10, 20, 30, 50)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression("Number of mixture components" ~ K),
       y = "Test accuracy", color = "", linetype = "") +
  guides(color = guide_legend(nrow = 1), fill = "none") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor.x = element_blank())
ggsave("realdata_sweep.pdf", plot = gg, device = cairo_pdf, width = 12,
       height = 4, units = "in")

}
