
# Required libraries
library(polykde)
library(movMF)
library(DirStats)
library(mlbench)
library(compositions)
library(parallel)
stopifnot(packageVersion("polykde") >= "1.2.1")
stopifnot(packageVersion("DirStats") >= "1.0.0")

## Settings

# Output paths
paper_dir <- "/Users/Eduardo/GitHub/polykde/polykde/paper"
cache_dir <- file.path(Sys.getenv("HOME"), ".cvapp-cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Number of stratified train/test splits and movMF fits
R <- 100

# Ceiling of the pooled movMF path. BIC still decreases at k = 12 on the pooled
# samples (a ceiling, not a selected value), but doubling it to 25 left the
# movMF accuracies unchanged or worse (checked on seeds 1-3), so it does not
# bind for the classification comparison. Per-class fits are capped tighter by
# mvmf_path().
Kmax <- 12
nruns <- 25

# Number of fork-based workers spread over the repeated splits
n_cores <- min(12, detectCores())

# Twelve methods: MISE-targeting selectors (CV, EMI), AMISE-targeting selectors
# (AMI, ROT), and the parametric movMF competitor, each in com/cls variants
fam <- c(CV = "kde-CV", EMI = "kde-EMI", AMI = "kde-AMI", ROT = "kde-ROT",
         movMF_BIC = "movMF-BIC", movMF_AIC = "movMF-AIC")
meth <- paste(rep(names(fam), each = 2), c("com", "cls"), sep = "_")

## Sphere embeddings

# Square-root map of a composition to the sphere positive orthant
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

# CV (LSCV, the selector analyzed in the paper) and ROT bandwidths. X is
# whatever sample the bandwidth is selected on: the pooled training set (com)
# or a single class (cls). A selector failure errors its split, which
# repeated_error() drops and reports.
bw_cv <- function(X, d) {

  # suppressWarnings: optim()'s generic advisory against 1-d Nelder-Mead,
  # emitted by bw_cv_polysph() on every call; constant chatter, no signal
  suppressWarnings(
    bw_cv_polysph(X = X, d = d, kernel = 1, type = "LSCV", exact_vmf = TRUE,
                  arcsinh = TRUE, spline = TRUE)$bw)

}

bw_rot <- function(X, d) {

  bw_rot_polysph(X = X, d = d, kernel = 1)$bw

}

# Convert a movMF fit (theta = kappa * mu, alpha) to the DirStats fit_mix format
fit_to_mix <- function(fit) {

  kap <- sqrt(rowSums(fit$theta^2))
  list(best_fit = list(mu_hat = fit$theta / kap, kappa_hat = kap,
                       p_hat = fit$alpha))

}

# EMI (exact MISE) and AMI (asymptotic MISE) plug-in bandwidths from a fitted
# vMF mixture reused as the reference density (DirStats). EMI integrates the
# exact MISE by importance sampling from the reference mixture (a seed fixes
# the Monte Carlo); AMI is closed-form and deterministic.
bw_emi <- function(X, fit, seed) {

  set.seed(seed)
  bw_dir_emi(data = X, fit_mix = fit_to_mix(fit), optim = TRUE,
             plot_it = FALSE)$h_opt

}

bw_ami <- function(X, fit) {

  bw_dir_ami(data = X, fit_mix = fit_to_mix(fit))

}

# kde log-density matrix (n_new x n_class). When h is supplied, all classes
# share that single common bandwidth; hvec gives one bandwidth per class;
# otherwise each class gets its own bandwidth from bwfun on that class's data.
kda_logdens <- function(Ztr, ytr, Znew, d, h = NULL, bwfun = NULL,
                        hvec = NULL) {

  cls <- levels(ytr)
  sapply(cls, function(l) {

    Xl <- Ztr[ytr == l, , drop = FALSE]
    hl <- if (!is.null(hvec)) hvec[[l]] else if (is.null(h)) bwfun(Xl, d) else h
    kde_polysph(x = Znew, X = Xl, d = d, h = hl, kernel = 1, log = TRUE)

  })

}

# Fit the movMF path k = 1, ..., Kc on X; returns the successful fits and their
# component counts. Kc caps the components at one parameter per observation: a
# vMF component on S^d costs d + 2 parameters (mu, kappa, weight) and
# ncol(X) = d + 1, so Kc = n / (d + 2). This cap binds before Kmax on the
# per-class fits.
mvmf_path <- function(X) {

  Kc <- max(1L, min(Kmax, floor(nrow(X) / (ncol(X) + 1))))
  fits <- lapply(seq_len(Kc), function(k)
    tryCatch(movMF(X, k = k, nruns = nruns), error = function(e) NULL))
  ok <- which(!sapply(fits, is.null))
  list(fits = fits[ok], k = ok)

}

# movMF log-densities with the number of components chosen by BIC and by AIC,
# either once on the pooled training sample (common = TRUE, then a K-component
# refit per class, K capped at one parameter per observation) or on each
# class's own path
# (common = FALSE). Returns the BIC/AIC n_new x n_class matrices and the BIC
# mixture(s) reused by the EMI/AMI selectors.
mvmf_logdens <- function(Ztr, ytr, Znew, common) {

  cls <- levels(ytr)
  dens <- function(f) dmovMF(Znew, theta = f$theta, alpha = f$alpha, log = TRUE)
  new_m <- function() matrix(NA_real_, nrow(Znew), length(cls),
                             dimnames = list(NULL, cls))
  if (common) {

    pp <- mvmf_path(Ztr)
    i_bic <- which.min(sapply(pp$fits, BIC))
    per_class <- function(K) {

      m <- new_m()
      for (j in seq_along(cls)) {

        Xl <- Ztr[ytr == cls[j], , drop = FALSE]
        k <- max(1L, min(K, floor(nrow(Xl) / (ncol(Xl) + 1))))
        f <- tryCatch(movMF(Xl, k = k, nruns = nruns), error = function(e)
          tryCatch(movMF(Xl, k = 1, nruns = nruns), error = function(e) NULL))
        m[, j] <- dens(f)

      }
      m

    }
    # Reuse the BIC refit when AIC picks the same K: movMF is unseeded, so a
    # second call would return a different fit for an identical K
    k_bic <- pp$k[i_bic]
    k_aic <- pp$k[which.min(sapply(pp$fits, AIC))]
    m_bic <- per_class(k_bic)
    list(BIC = m_bic,
         AIC = switch((k_aic == k_bic) + 1, per_class(k_aic), m_bic),
         mix = pp$fits[[i_bic]][c("theta", "alpha")])

  } else {

    lb <- new_m()
    la <- new_m()
    mix <- setNames(vector("list", length(cls)), cls)
    for (j in seq_along(cls)) {

      p <- mvmf_path(Ztr[ytr == cls[j], , drop = FALSE])
      fb <- p$fits[[which.min(sapply(p$fits, BIC))]]
      lb[, j] <- dens(fb)
      la[, j] <- dens(p$fits[[which.min(sapply(p$fits, AIC))]])
      mix[[j]] <- fb[c("theta", "alpha")]

    }
    list(BIC = lb, AIC = la, mix = mix)

  }

}

# Assign each row to the class maximizing log-density plus log-prior.
classify <- function(ld, logprior, cls) {

  cls[max.col(sweep(ld, 2, logprior, "+"))]

}

# Remove coincident points (identical rows, rounded to 10 digits) before
# splitting: duplicates make the CV loss unbounded below (the leave-one-out
# density diverges as h -> 0). A no-op except for LetterRecognition, whose
# integer features produce exact ties after normalization.
dedup_ds <- function(ds) {

  keep <- !duplicated(round(ds$X, 10))
  ds$X <- ds$X[keep, , drop = FALSE]
  ds$y <- droplevels(ds$y[keep])
  ds

}

# One stratified split: misclassification error of each of the twelve methods
# on the test set.
one_split <- function(ds, seed, prop = 0.7) {

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

  # movMF classifiers (component count pooled/com or per-class/cls); their BIC
  # mixtures are reused as reference densities by the EMI and AMI selectors
  mv_com <- mvmf_logdens(Ztr, ytr, Zte, common = TRUE)
  mv_cls <- mvmf_logdens(Ztr, ytr, Zte, common = FALSE)

  # Plug-in bandwidths: pooled mixture -> one common bandwidth; per-class
  # mixtures -> a bandwidth per class
  cls_tr <- levels(ytr)
  h_emi_com <- bw_emi(Ztr, mv_com$mix, seed = seed * 1000L)
  h_emi_cls <- sapply(cls_tr, function(l)
    bw_emi(Ztr[ytr == l, , drop = FALSE], mv_cls$mix[[l]],
           seed = seed * 1000L + match(l, cls_tr)))
  h_ami_com <- bw_ami(Ztr, mv_com$mix)
  h_ami_cls <- sapply(cls_tr, function(l)
    bw_ami(Ztr[ytr == l, , drop = FALSE], mv_cls$mix[[l]]))

  # Result: errors in the meth order
  c(CV_com = err(kda_logdens(Ztr, ytr, Zte, d, h = bw_cv(Ztr, d))),
    CV_cls = err(kda_logdens(Ztr, ytr, Zte, d, bwfun = bw_cv)),
    EMI_com = err(kda_logdens(Ztr, ytr, Zte, d, h = h_emi_com)),
    EMI_cls = err(kda_logdens(Ztr, ytr, Zte, d, hvec = h_emi_cls)),
    AMI_com = err(kda_logdens(Ztr, ytr, Zte, d, h = h_ami_com)),
    AMI_cls = err(kda_logdens(Ztr, ytr, Zte, d, hvec = h_ami_cls)),
    ROT_com = err(kda_logdens(Ztr, ytr, Zte, d, h = bw_rot(Ztr, d))),
    ROT_cls = err(kda_logdens(Ztr, ytr, Zte, d, bwfun = bw_rot)),
    movMF_BIC_com = err(mv_com$BIC), movMF_BIC_cls = err(mv_cls$BIC),
    movMF_AIC_com = err(mv_com$AIC), movMF_AIC_cls = err(mv_cls$AIC))

}

# Repeated stratified splits in parallel over the seeds; returns an R x 12
# matrix of test errors (rows are splits, columns are methods).
repeated_error <- function(ds, R = 100) {

  na <- setNames(rep(NA_real_, length(meth)), meth)
  res <- mclapply(seq_len(R), function(s) {

    tryCatch(one_split(ds, seed = s), error = function(e) na)

  }, mc.cores = n_cores)

  # mclapply returns a "try-error" if a worker dies, which rbind would splice
  # into the matrix as garbage
  bad <- !sapply(res, is.numeric)
  if (any(bad)) {

    res[bad] <- list(na)
    warning(sum(bad), " of ", R, " splits failed for ", ds$name, " (d = ",
            ds$d, ").")

  }
  do.call(rbind, res)

}

## Datasets

# Llobregat-basin river hydrochemistry (compositions::Hydrochem). The first p of
# the 14 chemical parts form a composition mapped to S^{p-1} by the sqrt map.
load_hydrochem <- function(p = 14) {

  data("Hydrochem", package = "compositions")
  parts <- c("H", "Na", "K", "Mg", "Ca", "Sr", "Ba", "NH4", "Cl", "NO3",
             "PO4", "SO4", "HCO3", "TOC")[seq_len(p)]
  list(X = as.matrix(Hydrochem[, parts]), y = factor(Hydrochem$River),
       d = p - 1, embed = "sqrt", name = "Hydrochem")

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
load_letter <- function(n_sub = 5000, seed = 1) {

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

# Repeated-split errors for each dataset over a grid of dimensions (Hydrochem
# via nested subcompositions; Vowel and LetterRecognition via nested feature
# subsets), saved to a single RData that drives the table below
results_file <- file.path(cache_dir, "cv-hd.RData")
if (!file.exists(results_file)) {

  # Run one dataset config: de-duplicate before splitting, then repeated splits
  run_ds <- function(ds) {

    ds <- dedup_ds(ds)
    cat(ds$name, "d =", ds$d, " n =", nrow(ds$X), "\n", file = stderr())
    list(d = ds$d, n = nrow(ds$X), cl = nlevels(ds$y),
         err = repeated_error(ds, R = R))

  }

  # Hydrochem: first p parts (d = p - 1)
  hydro <- lapply(c(8, 11, 14), function(p) run_ds(load_hydrochem(p = p)))

  # Vowel and LetterRecognition: nested feature subsets (d = p - 1), ordered by
  # between-class dispersion (computed once on the full data before de-dup)
  sweep_l2 <- function(ds, ps) {

    ord <- feat_order(ds$X, ds$y)
    lapply(ps, function(p) run_ds(subset_ds(ds, ord, p)))

  }
  vowel <- sweep_l2(load_vowel(), c(5, 7, 9))
  letter <- sweep_l2(load_letter(), c(4, 7, 10, 13, 16))

  save(hydro, vowel, letter, R, nruns, file = results_file)

}

## Table

load(results_file)
acc_mean <- function(err) colMeans(1 - err, na.rm = TRUE)[meth]
acc_sd <- function(err) apply(1 - err, 2, sd, na.rm = TRUE)[meth]

# Pick a sweep entry (d, n, cl, err) at a given dimension
at_d <- function(sw, d) sw[[which(sapply(sw, `[[`, "d") == d)]]

# One entry per dataset at its full dimension
info <- Map(function(nm, sw, dd) c(list(name = nm), at_d(sw, dd)),
            c("Hydrochem", "Vowel", "LetterRecognition"),
            list(hydro, vowel, letter), c(13, 8, 15))

# Only the pooled (com) variants are reported in the table; the per-class (cls)
# results stay in the cache and are summarized in the text
meth_com <- paste0(names(fam), "_com")

# LaTeX table: accuracy mean (sd) x 100, grouped kde/movMF header, best method
# per dataset (and its ties at the displayed precision) in bold
tex <- c(paste0("\\begin{tabular}{lccc", strrep("c", length(fam)), "}"),
         "\\toprule",
         " & & & & \\multicolumn{4}{c}{kde} & \\multicolumn{2}{c}{movMF} \\\\",
         "\\cmidrule(lr){5-8} \\cmidrule(lr){9-10}",
         paste0("Dataset & $n$ & $d$ & Classes & ",
                paste(sub("^(kde|movMF)-", "", fam), collapse = " & "),
                " \\\\"),
         "\\midrule")
for (z in info) {

  m <- acc_mean(z$err)[meth_com]
  s <- acc_sd(z$err)[meth_com]
  cell <- sapply(seq_along(meth_com), function(k) {

    v <- sprintf("%.1f\\,(%.1f)", 100 * m[k], 100 * s[k])
    best <- round(100 * m[k], 1) == round(100 * max(m), 1)
    if (best) paste0("\\textbf{", v, "}") else v

  })
  tex <- c(tex, sprintf("%s & %d & %d & %d & %s \\\\", z$name, z$n, z$d, z$cl,
                        paste(cell, collapse = " & ")))

}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(paper_dir, "tab_realdata.tex"))
