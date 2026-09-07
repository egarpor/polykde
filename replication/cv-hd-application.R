
# Required libraries
library(polykde)
library(movMF)
library(mlbench)
library(compositions)
library(parallel)
library(ggplot2)
library(patchwork)
stopifnot(packageVersion("polykde") >= "1.3.0")

## Sphere embeddings
{

# Square-root map of a composition to the sphere positive orthant
sqrt_map <- function(X) {

  X <- as.matrix(X)
  sqrt(X / rowSums(X))

}

# L2 normalization: columns centered/scaled with the training statistics
# center and scale to avoid leakage
l2_map <- function(X, center, scale) {

  X <- sweep(sweep(as.matrix(X), 2, center, "-"), 2, scale, "/")
  X / sqrt(rowSums(X^2))

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
  features <- paste0("V", 2:10)
  list(X = as.matrix(Vowel[, features]), y = factor(Vowel$Class),
       d = length(features) - 1, embed = "l2", name = "Vowel")

}

# Letter recognition (mlbench::LetterRecognition): 16 features map to S^15, with
# 26 classes
load_letter <- function() {

  data("LetterRecognition", package = "mlbench")
  features <- setdiff(names(LetterRecognition), "lettr")
  list(X = as.matrix(LetterRecognition[, features]),
       y = factor(LetterRecognition$lettr),
       d = length(features) - 1, embed = "l2", name = "LetterRecognition")

}

}

## Fit classifiers (bandwidth selectors and movMF)
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

# movMF fit with k components and n_runs_em random EM starts, capturing
# when the fits fail and measuring fitting time.
fit_movmf <- function(X, k, n_runs_em) {

  # Fit with timing and output capturing
  time <- system.time(
    msgs <- capture.output(fit <- tryCatch(
      movMF(X, k = k, nruns = n_runs_em, verbose = 1e9),
      error = function(e) NULL), type = "message"))[["elapsed"]]

  # Filter messages with error
  msgs <- trimws(msgs[startsWith(msgs, "Error")])

  # List with fit and information
  list(fit = fit, log = data.frame(
    fail = is.null(fit), n_runs_failed = length(msgs),
    iter = if (is.null(fit)) NA else fit$details$iter[["iter"]], time = time,
    msg = if (length(msgs)) names(which.max(table(msgs))) else ""))

}

}

## Classifiers
{

# kde log-density matrix (n_new x n_class): one kde per class on that class'
# training points, with h a common bandwidth (if scalar) or one bandwidth
# per class (if vector).
kda_log_dens <- function(X_train, y_train, X_new, d, h) {

  h <- rep_len(h, nlevels(y_train))
  sapply(seq_len(nlevels(y_train)), function(i) {

    X_lev <- X_train[y_train == levels(y_train)[i], , drop = FALSE]
    kde_polysph(x = X_new, X = X_lev, d = d, h = h[i], kernel = 1, log = TRUE)

  })

}

# Assign each row to the class maximizing log-density plus log-prior
classify <- function(log_dens, log_prior, classes) {

  # Ties are broken by the first level to avoid hidden randomization
  classes[max.col(sweep(log_dens, 2, log_prior, "+"), ties.method = "first")]

}

}

## Data splitting and embedding
{

# Remove coincident points before splitting
dedup_dataset <- function(dataset) {

  keep <- !duplicated(round(dataset$X, 10))
  dataset$X <- dataset$X[keep, , drop = FALSE]
  dataset$y <- droplevels(dataset$y[keep])
  dataset

}

# Stratified split with training proportion prop and leak-free data embedding
# (l2 uses the training center/scale)
split_embed <- function(dataset, seed, prop) {

  # Split into predictors and labels
  X <- dataset$X
  y <- dataset$y

  # Draw a proportion prop of the indices within each class
  set.seed(seed)
  ind_train <- unlist(lapply(levels(y), function(lev) {

    i <- which(y == lev)
    sample(i, round(prop * length(i)))

  }))

  # The l2 map needs a pre-standardization, the sqrt map does not
  if (dataset$embed == "sqrt") {

    X_train <- sqrt_map(X[ind_train, , drop = FALSE])
    X_test <- sqrt_map(X[-ind_train, , drop = FALSE])

  } else {

    center <- colMeans(X[ind_train, , drop = FALSE])
    scale <- apply(X[ind_train, , drop = FALSE], 2, sd)
    X_train <- l2_map(X[ind_train, , drop = FALSE], center, scale)
    X_test <- l2_map(X[-ind_train, , drop = FALSE], center, scale)

  }
  list(X_train = X_train, y_train = droplevels(y[ind_train]),
       X_test = X_test, y_test = y[-ind_train])

}

# Test accuracy and training time of kde-CV and kde-ROT (common and per-class
# bandwidths) and of movMF with r components per class, for each r in r_grid,
# and with r selected per class by BIC.
sweep_split <- function(dataset, seed, r_grid, n_runs_em, prop) {

  # Split and unpack
  split_data <- split_embed(dataset, seed, prop)
  X_train <- split_data$X_train
  y_train <- split_data$y_train
  X_test <- split_data$X_test
  y_test <- split_data$y_test
  d <- dataset$d
  classes <- levels(y_train)
  log_prior <- log(as.numeric(table(y_train)) / length(y_train))
  X_by_class <- lapply(classes, function(lev)
    X_train[y_train == lev, , drop = FALSE])

  ## KDE computations

  # Safe wrapper and timer (value plus elapsed seconds)
  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
  timed <- function(expr) {

    time <- system.time(value <- expr)[["elapsed"]]
    list(value = value, time = time)

  }

  # Common and per-class bandwidths
  cv_com <- timed(safe(bw_cv(X_train, d)))
  rot_com <- timed(safe(bw_rot(X_train, d)))
  cv_cls <- timed(safe(sapply(X_by_class, bw_cv, d = d)))
  rot_cls <- timed(safe(sapply(X_by_class, bw_rot, d = d)))

  # Test accuracy of the kda with a common or per-class bandwidths
  acc <- function(log_dens) {

    mean(classify(log_dens, log_prior, classes) == y_test)

  }
  acc_h <- function(h) safe(acc(kda_log_dens(X_train, y_train, X_test, d, h)))
  kde <- c(CV_com = acc_h(cv_com$value), CV_cls = acc_h(cv_cls$value),
           ROT_com = acc_h(rot_com$value), ROT_cls = acc_h(rot_cls$value))

  ## movMF computations

  # Exactly r components per class for each r in r_grid
  fits <- lapply(r_grid, function(r)
    lapply(X_by_class, fit_movmf, k = r, n_runs_em = n_runs_em))

  # Test accuracy at each r, NA if EM failed for any class
  acc_r <- sapply(fits, function(fits_r) {

    if (any(sapply(fits_r, function(f) is.null(f$fit)))) return(NA_real_)
    safe(acc(sapply(fits_r, function(f)
      dmovMF(X_test, theta = f$fit$theta, alpha = f$fit$alpha, log = TRUE))))

  })

  # movMF with r selected per class by BIC among its successful fits on the
  # grid (r = 1 never fails, so every class has a choice)
  bic <- sapply(fits, function(fits_r) sapply(fits_r, function(f)
    if (is.null(f$fit)) NA else BIC(f$fit)))
  j_bic <- apply(bic, 1, which.min)
  acc_bic <- safe(acc(sapply(seq_along(classes), function(c) {

    f <- fits[[j_bic[c]]][[c]]$fit
    dmovMF(X_test, theta = f$theta, alpha = f$alpha, log = TRUE)

  })))

  # EM log, one row per class and r, flagging the r selected by BIC
  em_log <- do.call(rbind, lapply(seq_along(r_grid), function(j)
    cbind(seed = seed, class = classes, r = r_grid[j],
          do.call(rbind, lapply(fits[[j]], `[[`, "log")))))
  em_log$bic_sel <- em_log$r == r_grid[j_bic][match(em_log$class, classes)]

  ## Final results

  # Accuracies, bandwidths, and training times
  t_r <- sapply(fits, function(fits_r) sum(sapply(fits_r,
                                                  function(f) f$log$time)))
  list(acc = c(kde, setNames(acc_r, paste0("movMF", r_grid)),
               movMF_bic_cls = acc_bic,
               h_CV = cv_com$value, h_ROT = rot_com$value,
               t_CV_com = cv_com$time, t_CV_cls = cv_cls$time,
               t_ROT_com = rot_com$time, t_ROT_cls = rot_cls$time,
               setNames(t_r, paste0("t_movMF", r_grid))),
       em_log = em_log)

}

}

## Experiments
{

# Number of stratified train/test splits and training proportion
M <- 100
prop <- 0.7

# Number of EM initializations per movMF fit
n_runs_em <- 25

# Cores for parallelization
n_cores <- 12

# r grid dense up to 10, then geometrically spaced
mixture_grid <- function(r_top) {

  k_geom <- round(10 * 1.25^(1:40))
  sort(unique(c(seq_len(min(10, r_top)), k_geom[k_geom < r_top], r_top)))

}

# Run one dataset: de-duplicate, set the r grid from the smallest class, sweep
# the splits in parallel, report the NA counts and the EM failures
run_dataset <- function(dataset, prop) {

  # Grid up to the cap floor(n_c / 3) of the smallest training class, so every
  # class can be asked for r components at every r on the grid
  dataset <- dedup_dataset(dataset)
  r_grid <- mixture_grid(floor(min(round(prop * table(dataset$y))) / 3))
  cat(dataset$name, "n =", nrow(dataset$X), "r_top =", max(r_grid), "\n",
      file = stderr())

  # One seed per split, run in parallel; an erroring split yields NULL
  res <- mclapply(seq_len(M), function(seed)
    tryCatch(sweep_split(dataset, seed = seed, r_grid = r_grid,
                         n_runs_em = n_runs_em, prop = prop),
             error = function(e) NULL), mc.cores = n_cores)

  # Accuracies (splits by rows) and EM logs; report the NA count per quantity
  acc_mat <- do.call(rbind, lapply(res, `[[`, "acc"))
  em_log <- do.call(rbind, lapply(res, `[[`, "em_log"))
  n_na <- colSums(is.na(acc_mat))
  if (any(n_na > 0)) {

    cat("  NA/", M, ": ", paste(names(n_na)[n_na > 0], n_na[n_na > 0],
                                sep = "=", collapse = ", "), "\n", sep = "",
        file = stderr())

  }

  # EM failures per r: share of failed fits (class level) and of failed runs,
  # share of best runs stopped by maxiter (movMF default 100), error messages
  fail_r <- tapply(em_log$fail, em_log$r, mean)
  runs_r <- tapply(em_log$n_runs_failed, em_log$r, sum) /
    (n_runs_em * table(em_log$r))
  cat("  EM failed fits per r (%):", round(100 * fail_r), "\n", file = stderr())
  cat("  EM failed runs per r (%):", round(100 * runs_r), "\n", file = stderr())
  cat("  best EM runs at maxiter (%):",
      round(100 * mean(em_log$iter >= 100, na.rm = TRUE)), "\n",
      file = stderr())
  print(table(em_log$msg[em_log$msg != ""]))
  list(name = dataset$name, n = nrow(dataset$X), d = dataset$d,
       n_classes = nlevels(dataset$y), r_grid = r_grid, acc = acc_mat,
       em_log = em_log)

}

# Each dataset is saved as soon as it finishes (cv-hd-<dataset>-M<M>.RData),
# so an interrupted run resumes and earlier M's stay reloadable
datasets <- list(load_hydrochem(), load_vowel(), load_letter())
res <- lapply(datasets, function(dataset) {

  results_file <- sprintf("cv-hd-%s-M%d.RData", dataset$name, M)
  if (!file.exists(results_file)) {

    res_dataset <- run_dataset(dataset, prop)
    save(res_dataset, n_runs_em, file = results_file)

  }
  load(results_file)
  res_dataset

})

}

## Paper numbers and figure
{

# Load the results saved by the experiments
M <- 100
res <- lapply(c("Hydrochem", "Vowel", "LetterRecognition"), function(name) {

  load(sprintf("cv-hd-%s-M%d.RData", name, M), envir = globalenv())
  res_dataset

})

# Mean (standard deviation) of a quantity across the splits
mean_sd <- function(v) sprintf("%.1f (%.1f)", mean(v, na.rm = TRUE),
                               sd(v, na.rm = TRUE))

# Numbers quoted in the text
for (res_dataset in res) {

  # Accuracies as percentages: kdes with common ("com") and per-class ("cls")
  # bandwidths
  acc_mat <- 100 * res_dataset$acc
  cat(res_dataset$name, ": kde-CV", mean_sd(acc_mat[, "CV_com"]),
      " cls", mean_sd(acc_mat[, "CV_cls"]),
      " kde-ROT", mean_sd(acc_mat[, "ROT_com"]),
      " cls", mean_sd(acc_mat[, "ROT_cls"]), "\n")

  # Best mean accuracy of movMF along the r grid, over the splits with all
  # classes fitted
  acc_r <- acc_mat[, paste0("movMF", res_dataset$r_grid), drop = FALSE]
  mean_acc_r <- colMeans(acc_r, na.rm = TRUE)
  i_best <- which.max(mean_acc_r)
  r_best <- res_dataset$r_grid[i_best]
  cat("  movMF: max ", sprintf("%.1f (%.1f)", mean_acc_r[i_best],
                                sd(acc_r[, i_best], na.rm = TRUE)),
      " at r = ", r_best, " over ", sum(!is.na(acc_r[, i_best])),
      " splits\n", sep = "")
  cat("  movMF: splits per r:", colSums(!is.na(acc_r)), "\n")

  # movMF with r selected per class by BIC, and the range of the selected r's
  em_log <- res_dataset$em_log
  cat("  movMF-BIC per class", mean_sd(acc_mat[, "movMF_bic_cls"]),
      " r_c in", range(em_log$r[em_log$bic_sel]), "\n")

  # Median training time per split (s) and EM failure share at the oracle r
  t_med <- apply(res_dataset$acc[, c("t_CV_com", "t_CV_cls", "t_ROT_com",
                                     "t_ROT_cls", paste0("t_movMF", r_best))],
                 2, median, na.rm = TRUE)
  cat("  training time (s):", paste(names(t_med), sprintf("%.2f", t_med),
                                    sep = "=", collapse = " "), "\n")
  cat("  EM failed fits at r = ", r_best, ": ",
      round(100 * mean(em_log$fail[em_log$r == r_best])), "%\n", sep = "")

}

# Figure: one row per (dataset, method, r) with the mean accuracy across the
# splits with all classes fitted and its standard deviation, plus the share of
# EM runs that error
acc_df <- do.call(rbind, lapply(res, function(res_dataset) {

  # Summary of the given accuracy columns, one per r
  stat_cols <- function(cols, method) {

    acc_cols <- res_dataset$acc[, cols, drop = FALSE]
    data.frame(dataset = res_dataset$name, r = res_dataset$r_grid,
               method = method, acc = colMeans(acc_cols, na.rm = TRUE),
               sd = apply(acc_cols, 2, sd, na.rm = TRUE), row.names = NULL)

  }

  # The kde methods do not depend on r: their column is recycled along the grid
  n_r <- length(res_dataset$r_grid)
  em_log <- res_dataset$em_log
  runs_r <- tapply(em_log$n_runs_failed, em_log$r, sum) /
    (n_runs_em * table(em_log$r))
  rbind(stat_cols(rep("CV_com", n_r), "kde-CV"),
        stat_cols(rep("ROT_com", n_r), "kde-ROT"),
        stat_cols(paste0("movMF", res_dataset$r_grid), "movMF"),
        data.frame(dataset = res_dataset$name, r = res_dataset$r_grid,
                   method = "EM failure rate", acc = as.numeric(runs_r),
                   sd = 0, row.names = NULL))

}))
acc_df <- acc_df[!is.na(acc_df$acc), ]

# Panel and legend orders
acc_df$dataset <- factor(acc_df$dataset, levels = sapply(res, `[[`, "name"))
acc_df$method <- factor(acc_df$method, levels = c("kde-CV", "kde-ROT", "movMF",
                                                  "EM failure rate"))

# Palette of Figure 1, with the r-free kdes dashed, movMF solid, and the EM
# failure rate dotted
col_method <- c("kde-CV" = "green3", "kde-ROT" = "blue", "movMF" = "red",
                "EM failure rate" = "grey40")
lty_method <- c("kde-CV" = "dashed", "kde-ROT" = "dashed", "movMF" = "solid",
                "EM failure rate" = "dotted")

# One panel per dataset, with mean +/- standard deviation ribbons and r in
# log10 scale. The dataset names are not drawn: they go in the LaTeX
# subcaptions of the figure.
panels <- lapply(levels(acc_df$dataset), function(lev) {

  acc_lev <- acc_df[acc_df$dataset == lev, ]
  ggplot(acc_lev, aes(x = r, y = acc, color = method, linetype = method)) +
    geom_ribbon(aes(ymin = acc - sd, ymax = acc + sd, fill = method),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.6) +
    geom_point(data = acc_lev[acc_lev$method == "movMF", ], size = 0.9) +
    scale_color_manual(values = col_method) +
    scale_fill_manual(values = col_method) +
    scale_linetype_manual(values = lty_method) +
    scale_x_continuous(transform = "log10",
                       breaks = c(1, 2, 3, 5, 10, 20, 50, 100, 200)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = expression("Number of mixture components" ~ r),
         y = "Test accuracy / EM failure rate", color = "", linetype = "") +
    guides(color = guide_legend(nrow = 1), fill = "none") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

})

# Three panels side by side, sharing a single legend at the bottom
plot_sweep <- wrap_plots(panels, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")
ggsave("realdata_sweep.pdf", plot = plot_sweep, device = cairo_pdf,
       width = 12, height = 4, units = "in")

}
