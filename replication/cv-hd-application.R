
# This script reproduces the real-data application (Section 5) of the paper
# "Blessing of dimensionality in cross-validated bandwidth selection on the
# sphere" (Chacon, Garcia-Portugues, and Meilan-Vila). It compares kernel
# discriminant analysis with a cross-validated bandwidth against mixtures of
# von Mises-Fisher distributions (movMF) on three datasets embedded on the
# sphere: Hydrochem (geochemistry), Vowel (speech), and LetterRecognition
# (handwriting). It writes the results table (paper/tab_realdata.tex) and the
# accuracy-vs-dimension figure (paper/img/realdata_dim.pdf).
#
# Methods. Each class-conditional density has a tuning parameter set either on
# the whole training sample (com) or separately per class (cls), giving eight
# methods = {LSCV, ROT, movMF-BIC, movMF-AIC} x {com, cls}:
#   * kde-LSCV: least-squares cross-validation bandwidth, common or per-class.
#   * kde-ROT:  von Mises-Fisher rule-of-thumb plug-in bandwidth, common or
#               per-class.
#   * movMF:    number of mixture components chosen by BIC/AIC on the pooled data
#               (com, analogous to a common bandwidth) or per class (cls), then a
#               K-component vMF mixture refit to each class.
# The previous "median of the per-class selections" combine is not used.
#
# Tied observations. LSCV is unbounded below when the sample contains coincident
# points: for a point with a duplicate neighbour the leave-one-out density
# diverges as h -> 0, so CV(h) has no interior minimum and collapses to h_min.
# Each dataset is therefore de-duplicated (identical rows removed) before
# splitting. This is a no-op for the continuous datasets (Hydrochem, Vowel);
# LetterRecognition has integer features, hence exact ties after normalization.
#
# Tip: run with a single BLAS thread (e.g. OMP_NUM_THREADS=1) to avoid
# oversubscription across the parallel workers. Environment overrides:
#   CVAPP_R      number of stratified train/test splits (default 100)
#   CVAPP_CORES  number of parallel workers over the R splits (default 10)
#   CVAPP_IMG    output image directory
#   CVAPP_CACHE  cache directory

# Required libraries
library(polykde)
library(movMF)
library(mlbench)
library(compositions)
library(parallel)

## Settings

# Output paths (override with the environment variables if needed)
img_dir <- Sys.getenv("CVAPP_IMG",
  "/Users/Eduardo/GitHub/polykde/polykde/paper/img")
cache_dir <- Sys.getenv("CVAPP_CACHE",
  file.path(Sys.getenv("HOME"), ".cvapp-cache"))
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Number of stratified train/test splits and movMF fitting effort
R <- as.integer(Sys.getenv("CVAPP_R", "100"))
Kmax <- 12
nruns <- 5

# Minimum bandwidth safeguard. kappa = 1 / h^2, so h = 0.005 caps kappa at 4e4,
# below the ~1.1e5 where the vMF normalizing constant overflows; binds only on
# near-degenerate classes, never on the datasets reported here.
h_min <- 0.005

# Number of fork-based workers spread over the repeated splits
n_cores <- min(as.integer(Sys.getenv("CVAPP_CORES", "10")), detectCores())

## Sphere embeddings

# Square-root map of a composition to the positive orthant of S^{P-1}. Rows of X
# must be nonnegative (zeros allowed: the point lands on a sub-face).
sqrt_map <- function(X) {

  X <- as.matrix(X)
  sqrt(X / rowSums(X))

}

# L2 map of general features to S^{P-1}. Columns are centered/scaled with the
# training statistics ctr and scl (passed in) to avoid leakage.
l2_map <- function(X, ctr, scl) {

  X <- sweep(sweep(as.matrix(X), 2, ctr, "-"), 2, scl, "/")
  X / sqrt(rowSums(X^2))

}

## Density-based classification

# Least-squares cross-validation (LSCV, the selector analyzed in the paper) and
# rule-of-thumb (ROT) bandwidths, floored at h_min. Degenerate LSCV fits fall
# back to ROT. The argument X is whatever sample the bandwidth is selected on:
# the whole (de-duplicated) training set for the common bandwidth, or a single
# class for the per-class bandwidths.
bw_lscv <- function(X, d) {

  # The Bessel spline (fast, bit-exact vs besselI) is valid for nu = (p-2)/2 in
  # {0, 0.5, ..., 6}, i.e. d <= 13; d = 15 (nu = 7) must use the exact loss.
  h <- suppressWarnings(tryCatch(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV", exact_vmf = TRUE,
                  arcsinh = TRUE, spline = (d <= 13))$bw,
    error = function(e) bw_rot_polysph(X = X, d = d, kernel = 1)$bw))
  max(h_min, h)

}

bw_rot <- function(X, d) {

  max(h_min, bw_rot_polysph(X = X, d = d, kernel = 1)$bw)

}

# kde log-density matrix (n_new x n_class). When h is supplied, all classes
# share that single common bandwidth; otherwise each class gets its own
# bandwidth from bwfun applied to that class's data.
kda_logdens <- function(Ztr, ytr, Znew, d, h = NULL, bwfun = NULL) {

  cls <- levels(ytr)
  sapply(cls, function(l) {

    Xl <- Ztr[ytr == l, , drop = FALSE]
    hl <- if (is.null(h)) bwfun(Xl, d) else h
    kde_polysph(x = Znew, X = Xl, d = d, h = hl, kernel = 1, log = TRUE)

  })

}

# Fit the movMF path k = 1, ..., Kc on X (Kc capped by the sample size); returns
# the successful fits and their component counts.
mvmf_path <- function(X, Kmax = 12, nruns = 5) {

  Kc <- max(1L, min(Kmax, floor(nrow(X) / 5)))
  fits <- lapply(seq_len(Kc), function(k)
    tryCatch(movMF(X, k = k, nruns = nruns), error = function(e) NULL))
  ok <- which(!vapply(fits, is.null, logical(1)))
  list(fits = fits[ok], k = ok)

}

# movMF log-densities with the number of components chosen PER CLASS by BIC and
# by AIC on each class's own path. Returns the BIC and AIC n_new x n_class
# matrices.
mvmf_logdens <- function(Ztr, ytr, Znew, Kmax = 12, nruns = 5) {

  cls <- levels(ytr)
  dens <- function(f) dmovMF(Znew, theta = f$theta, alpha = f$alpha, log = TRUE)
  lb <- la <- matrix(NA_real_, nrow(Znew), length(cls),
                     dimnames = list(NULL, cls))
  for (j in seq_along(cls)) {

    p <- mvmf_path(Ztr[ytr == cls[j], , drop = FALSE], Kmax, nruns)
    lb[, j] <- dens(p$fits[[which.min(sapply(p$fits, BIC))]])
    la[, j] <- dens(p$fits[[which.min(sapply(p$fits, AIC))]])

  }
  list(BIC = lb, AIC = la)

}

# movMF log-densities with the number of components chosen ONCE on the pooled
# training sample by BIC and by AIC (analogous to the common LSCV bandwidth
# selected on the whole data), then a K-component movMF refit to each class (K
# capped by the class size). Returns the BIC and AIC n_new x n_class matrices.
mvmf_logdens_common <- function(Ztr, ytr, Znew, Kmax = 12, nruns = 5) {

  cls <- levels(ytr)
  dens <- function(f) dmovMF(Znew, theta = f$theta, alpha = f$alpha, log = TRUE)
  pp <- mvmf_path(Ztr, Kmax, nruns)                  # pooled path selects K
  K_bic <- pp$k[which.min(sapply(pp$fits, BIC))]
  K_aic <- pp$k[which.min(sapply(pp$fits, AIC))]
  per_class <- function(K) {

    m <- matrix(NA_real_, nrow(Znew), length(cls), dimnames = list(NULL, cls))
    for (j in seq_along(cls)) {

      Xl <- Ztr[ytr == cls[j], , drop = FALSE]
      k <- max(1L, min(K, floor(nrow(Xl) / 5)))
      f <- tryCatch(movMF(Xl, k = k, nruns = nruns), error = function(e)
        tryCatch(movMF(Xl, k = 1, nruns = nruns), error = function(e) NULL))
      m[, j] <- dens(f)

    }
    m

  }
  list(BIC = per_class(K_bic), AIC = per_class(K_aic))

}

# Assign each row to the class maximizing log-density plus log-prior.
classify <- function(ld, logprior, cls) {

  cls[max.col(sweep(ld, 2, logprior, "+"))]

}

# Remove coincident points (identical rows, rounded to 10 digits) from a loaded
# dataset before splitting. Duplicates make the LSCV loss unbounded below (the
# leave-one-out density diverges as h -> 0), so its minimizer collapses to h_min;
# this is a no-op for the continuous datasets (Hydrochem, Vowel), while
# LetterRecognition's integer features produce exact ties that are removed here.
dedup_ds <- function(ds) {

  keep <- !duplicated(round(ds$X, 10))
  ds$X <- ds$X[keep, , drop = FALSE]
  ds$y <- droplevels(ds$y[keep])
  ds

}

# One stratified split: misclassification error of each method on the test set.
# For each of LSCV, ROT, and movMF (BIC/AIC), the tuning parameter (bandwidth or
# number of components) is set either on the whole training sample (com) or per
# class (cls). De-duplication is done at the dataset level, before splitting.
one_split <- function(ds, seed, prop = 0.7, Kmax = 12, nruns = 5) {

  X <- ds$X
  y <- ds$y
  d <- ds$d
  cls <- levels(y)
  set.seed(seed)
  tr <- unlist(lapply(cls, function(l) {

    i <- which(y == l)
    sample(i, max(2L, round(prop * length(i))))

  }))
  ytr <- droplevels(y[tr])
  yte <- y[-tr]

  # Embed (leak-free): sqrt is row-wise; l2 uses the training center/scale
  if (ds$embed == "sqrt") {

    Ztr <- sqrt_map(X[tr, , drop = FALSE])
    Zte <- sqrt_map(X[-tr, , drop = FALSE])

  } else {

    ctr <- colMeans(X[tr, , drop = FALSE])
    scl <- apply(X[tr, , drop = FALSE], 2, sd)
    Ztr <- l2_map(X[tr, , drop = FALSE], ctr, scl)
    Zte <- l2_map(X[-tr, , drop = FALSE], ctr, scl)

  }

  logprior <- log(as.numeric(table(ytr)) / length(ytr))
  err <- function(ld) mean(classify(ld, logprior, levels(ytr)) != yte)

  # kde: common vs per-class bandwidth, for the LSCV and ROT selectors
  ld_lscv_c <- kda_logdens(Ztr, ytr, Zte, d, h = bw_lscv(Ztr, d))
  ld_lscv_p <- kda_logdens(Ztr, ytr, Zte, d, bwfun = bw_lscv)
  ld_rot_c <- kda_logdens(Ztr, ytr, Zte, d, h = bw_rot(Ztr, d))
  ld_rot_p <- kda_logdens(Ztr, ytr, Zte, d, bwfun = bw_rot)
  # movMF: number of components chosen on the pooled data (com) or per class (cls)
  mv_com <- mvmf_logdens_common(Ztr, ytr, Zte, Kmax = Kmax, nruns = nruns)
  mv_cls <- mvmf_logdens(Ztr, ytr, Zte, Kmax = Kmax, nruns = nruns)
  c(LSCV_com = err(ld_lscv_c),
    LSCV_cls = err(ld_lscv_p),
    ROT_com = err(ld_rot_c),
    ROT_cls = err(ld_rot_p),
    movMF_BIC_com = err(mv_com$BIC),
    movMF_BIC_cls = err(mv_cls$BIC),
    movMF_AIC_com = err(mv_com$AIC),
    movMF_AIC_cls = err(mv_cls$AIC))

}

# Repeated stratified splits in parallel over the seeds; returns an R x 8 matrix
# of test errors (rows are splits, columns are methods).
repeated_error <- function(ds, R = 100, prop = 0.7, Kmax = 12, nruns = 5) {

  na <- c(LSCV_com = NA, LSCV_cls = NA, ROT_com = NA, ROT_cls = NA,
          movMF_BIC_com = NA, movMF_BIC_cls = NA,
          movMF_AIC_com = NA, movMF_AIC_cls = NA)
  res <- mclapply(seq_len(R), function(s) {

    tryCatch(one_split(ds, seed = s, prop = prop, Kmax = Kmax, nruns = nruns),
             error = function(e) na)

  }, mc.cores = n_cores)
  do.call(rbind, res)

}

## Datasets

# Llobregat-basin river hydrochemistry (compositions::Hydrochem). The first P of
# the 14 chemical parts form a composition mapped to S^{P-1} by the sqrt map.
load_hydrochem <- function(P = 14) {

  data("Hydrochem", package = "compositions")
  parts <- c("H", "Na", "K", "Mg", "Ca", "Sr", "Ba", "NH4", "Cl", "NO3",
             "PO4", "SO4", "HCO3", "TOC")[seq_len(P)]
  list(X = as.matrix(Hydrochem[, parts]), y = factor(Hydrochem$River),
       d = P - 1, embed = "sqrt", name = "Hydrochem")

}

# Deterding vowel recognition (mlbench::Vowel): the 9 LPC features V2:V10 (V1 is
# a speaker indicator) map to S^8, with 11 vowel classes.
load_vowel <- function() {

  data("Vowel", package = "mlbench")
  feat <- paste0("V", 2:10)
  list(X = as.matrix(Vowel[, feat]), y = factor(Vowel$Class),
       d = length(feat) - 1, embed = "l2", name = "Vowel")

}

# Letter recognition (mlbench::LetterRecognition): 16 features map to S^15, with
# 26 classes; subsampled to n_sub rows for tractable repeated splits.
load_letter <- function(n_sub = 3000, seed = 1) {

  data("LetterRecognition", package = "mlbench")
  set.seed(seed)
  i <- sample(nrow(LetterRecognition), min(n_sub, nrow(LetterRecognition)))
  feat <- setdiff(names(LetterRecognition), "lettr")
  list(X = as.matrix(LetterRecognition[i, feat]),
       y = droplevels(LetterRecognition$lettr[i]),
       d = length(feat) - 1, embed = "l2", name = "LetterRecognition")

}

# Order the columns of X by decreasing between-class dispersion, so that nested
# prefixes give an informative sequence of dimensions for the L2 datasets.
feat_order <- function(X, y) {

  disp <- apply(X, 2, function(x) var(tapply(x, y, mean)))
  order(disp, decreasing = TRUE)

}

# Restrict a loaded L2 dataset to its first p features (ord order), yielding a
# lower-dimensional version with d = p - 1.
subset_ds <- function(ds, ord, p) {

  ds$X <- ds$X[, ord[seq_len(p)], drop = FALSE]
  ds$d <- p - 1
  ds

}

## Experiments

# Run (or reload) the repeated-split errors for each dataset over a grid of
# dimensions: Hydrochem via nested subcompositions; Vowel and LetterRecognition
# via nested feature subsets ordered by between-class dispersion.
results_file <- file.path(cache_dir, "cv-hd-application.RData")
if (file.exists(results_file)) {

  load(results_file)

} else {

  # Run one dataset config: de-duplicate before splitting, then repeated splits.
  # Records the analyzed (de-duplicated) sample size n and number of classes.
  run_ds <- function(ds) {

    ds <- dedup_ds(ds)
    cat(ds$name, "d =", ds$d, " n =", nrow(ds$X), "\n")
    list(d = ds$d, n = nrow(ds$X), cl = nlevels(ds$y),
         err = repeated_error(ds, R = R, Kmax = Kmax, nruns = nruns))

  }

  # Hydrochem: first P parts (d = P - 1)
  hydro <- lapply(c(8, 11, 14), function(P) run_ds(load_hydrochem(P = P)))

  # Vowel and LetterRecognition: nested feature subsets (d = p - 1), ordered by
  # between-class dispersion (computed once on the full data before de-dup)
  sweep_l2 <- function(ds, ps) {

    ord <- feat_order(ds$X, ds$y)
    lapply(ps, function(p) run_ds(subset_ds(ds, ord, p)))

  }
  vowel <- sweep_l2(load_vowel(), c(5, 7, 9))
  letter <- sweep_l2(load_letter(), c(4, 7, 10, 13, 16))

  save(hydro, vowel, letter, R, file = results_file)

}

## Results table

# Methods (common/per-class for each selector), colors, and pretty labels shared
# by the table and the figure
meth <- c("LSCV_com", "LSCV_cls", "ROT_com", "ROT_cls",
          "movMF_BIC_com", "movMF_BIC_cls", "movMF_AIC_com", "movMF_AIC_cls")
meth_lab <- c("kde-LSCV (com)", "kde-LSCV (cls)", "kde-ROT (com)",
              "kde-ROT (cls)", "movMF-BIC (com)", "movMF-BIC (cls)",
              "movMF-AIC (com)", "movMF-AIC (cls)")
acc_mean <- function(err) colMeans(1 - err, na.rm = TRUE)[meth]
acc_sd <- function(err) apply(1 - err, 2, sd, na.rm = TRUE)[meth]

# Pick a sweep entry (d, n, cl, err) at a given dimension
at_d <- function(sweep, d) sweep[[which(sapply(sweep, `[[`, "d") == d)]]

# One row per dataset at a representative dimension, with the analyzed
# (de-duplicated) sample size n, the sphere dimension d, and the class count
info <- lapply(list(list("Hydrochem", hydro, 10), list("Vowel", vowel, 8),
                    list("LetterRecognition", letter, 15)), function(u) {

  z <- at_d(u[[2]], u[[3]])
  list(name = u[[1]], n = z$n, d = z$d, cl = z$cl, err = z$err)

})

# Header: kde-LSCV/kde-ROT/movMF-BIC/movMF-AIC each in com and cls variants
head_cols <- paste(sprintf("%s$_{\\mathrm{com}}$ & %s$_{\\mathrm{cls}}$",
                           c("kde-LSCV", "kde-ROT", "movMF-BIC", "movMF-AIC"),
                           c("kde-LSCV", "kde-ROT", "movMF-BIC", "movMF-AIC")),
                   collapse = " & ")
tex <- c("\\begin{tabular}{lccccccccccc}", "\\toprule",
         paste0("Dataset & $n$ & $d$ & Classes & ", head_cols, " \\\\"),
         "\\midrule")
for (z in info) {

  m <- acc_mean(z$err)
  s <- acc_sd(z$err)
  best <- meth[which.max(m)]
  cell <- function(k) {

    v <- sprintf("%.1f\\,(%.1f)", 100 * m[k], 100 * s[k])
    if (meth[k] == best) paste0("\\textbf{", v, "}") else v

  }
  cells <- paste(sapply(seq_along(meth), cell), collapse = " & ")
  tex <- c(tex, sprintf("%s & %d & %d & %d & %s \\\\",
                        z$name, z$n, z$d, z$cl, cells))

}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(img_dir, "..", "tab_realdata.tex"))

## Figure: accuracy vs dimension

# Encode the four selector families by color and the common/per-class variant by
# line type (common = solid, per-class = dashed).
fam_col <- c(LSCV = "#1b7837", ROT = "#984ea3", movMF_BIC = "#e41a1c",
             movMF_AIC = "#ff7f00")
col <- setNames(unname(fam_col[sub("_(com|cls)$", "", meth)]), meth)
lty <- setNames(ifelse(grepl("_com$", meth), 1, 2), meth)
panels <- list(Hydrochem = hydro, Vowel = vowel,
               LetterRecognition = letter)
pdf(file.path(img_dir, "realdata_dim.pdf"), width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (nm in names(panels)) {

  sweep <- panels[[nm]]
  dd <- sapply(sweep, `[[`, "d")
  acc <- sapply(sweep, function(z) acc_mean(z$err))
  matplot(dd, t(acc), type = "n", ylim = range(acc), main = nm,
          xlab = "Sphere dimension d", ylab = "Test accuracy", xaxt = "n")
  axis(1, at = dd)
  for (k in meth) {

    lines(dd, acc[k, ], col = col[k], lwd = 2, lty = lty[k], type = "o",
          pch = 19, cex = 0.7)

  }
  if (nm == names(panels)[1]) {

    legend("bottomright", meth_lab, col = col, lwd = 2, lty = lty, pch = 19,
           bty = "n", cex = 0.8)

  }

}
dev.off()
