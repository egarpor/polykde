
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

## Bandwidth selectors and movMF fits
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

# movMF fit with k components and n_runs_em random EM starts, with error
# catching
fit_movmf <- function(X, k, n_runs_em) {

  tryCatch(movMF(X, k = k, nruns = n_runs_em), error = function(e) NULL)

}

}

## Classifiers
{

# kde log-density matrix (n_new x n_class). Two modes: common bandwidth h
# learned on the whole training set, or per-class bandwidths learned on each
# class' and computed with bw_fun(). In any case, the kda compares the kdes
# of each class' training points at the new data points.
kda_log_dens <- function(X_train, y_train, X_new, d, h = NULL, bw_fun = NULL) {

  # One column per class, each a kde fitted on that class' training points
  sapply(levels(y_train), function(lev) {

    X_lev <- X_train[y_train == lev, , drop = FALSE]
    h_lev <- if (is.null(h)) bw_fun(X_lev, d) else h
    kde_polysph(x = X_new, X = X_lev, d = d, h = h_lev, kernel = 1, log = TRUE)

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

# Stratified 70/30 split and leak-free data embedding (l2 uses the training
# center/scale)
split_embed <- function(dataset, seed, prop = 0.70) {

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

# Test accuracy of kde-CV and kde-ROT (common and per-class bandwidths) and of
# movMF at each r in r_grid (per-class r-component fits)
sweep_split <- function(dataset, seed, r_grid, n_runs_em) {

  # Split and unpack
  split_data <- split_embed(dataset, seed)
  X_train <- split_data$X_train
  y_train <- split_data$y_train
  X_test <- split_data$X_test
  y_test <- split_data$y_test
  d <- dataset$d
  classes <- levels(y_train)
  log_prior <- log(as.numeric(table(y_train)) / length(y_train))

  ## KDE computations

  # Safe wrapper
  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)

  # Compute the common kde bandwidths
  h_cv <- safe(bw_cv(X_train, d))
  h_rot <- safe(bw_rot(X_train, d))

  # Test accuracy functions for kde, with common and per-class variants
  acc <- function(log_dens) {

    mean(classify(log_dens, log_prior, classes) == y_test)

  }
  acc_h <- function(h) acc(kda_log_dens(X_train, y_train, X_test, d, h = h))
  acc_bw <- function(bw_fun) acc(kda_log_dens(X_train, y_train, X_test, d,
                                              bw_fun = bw_fun))
  kde <- c(CV_com = safe(acc_h(h_cv)), CV_cls = safe(acc_bw(bw_cv)),
           ROT_com = safe(acc_h(h_rot)), ROT_cls = safe(acc_bw(bw_rot)))

  ## movMF computations

  # Per-class fits at the distinct grid k's below the class cap n_c / 3. A fit
  # in which EM fails is dropped and counted (n_fail out of n_try).
  class_fits <- lapply(classes, function(lev) {

    X <- X_train[y_train == lev, , drop = FALSE]
    k_vals <- unique(pmin(r_grid, floor(nrow(X) / 3)))
    fits <- lapply(k_vals, function(k) fit_movmf(X, k = k,
                                                 n_runs_em = n_runs_em))
    ok <- !sapply(fits, is.null)
    list(fits = fits[ok], k = k_vals[ok], n_try = length(k_vals),
         n_fail = sum(!ok))

  })

  # Test accuracy at each shared r, each class using its largest fitted k <= r
  acc_r <- sapply(r_grid, function(r) safe(acc(sapply(class_fits,
    function(class_fit) {

      fit <- class_fit$fits[[max(which(class_fit$k <= r))]]
      dmovMF(X_test, theta = fit$theta, alpha = fit$alpha, log = TRUE)

    }))))

  ## Final results

  # Named accuracies, kde bandwidths, and attempted/failed movMF fit counts
  c(kde, setNames(acc_r, paste0("movMF", r_grid)),
    h_CV = h_cv, h_ROT = h_rot,
    n_try = sum(sapply(class_fits, `[[`, "n_try")),
    n_fail = sum(sapply(class_fits, `[[`, "n_fail")))

}

}

## Experiments
{

# Number of stratified train/test splits
M <- 100

# Number of EM initializations per movMF fit
n_runs_em <- 25

# Cores for parallelization
n_cores <- 12

# r grid dense up to 10, then geometrically spaced
mixture_grid <- function(r_top) {

  k_geom <- round(10 * 1.25^(1:40))
  sort(unique(c(seq_len(min(10, r_top)), k_geom[k_geom < r_top], r_top)))

}

# Run one dataset: de-duplicate, set the r grid from the largest class, sweep
# the splits in parallel, report the failure counts
run_dataset <- function(dataset) {

  dataset <- dedup_dataset(dataset)
  r_grid <- mixture_grid(floor(0.7 * max(table(dataset$y)) / 3))
  cat(dataset$name, "n =", nrow(dataset$X), " r_top =", max(r_grid), "\n",
      file = stderr())

  # One seed per split, run in parallel; an erroring split yields NULL
  res <- mclapply(seq_len(M), function(seed)
    tryCatch(sweep_split(dataset, seed = seed, r_grid = r_grid,
                         n_runs_em = n_runs_em),
             error = function(e) NULL), mc.cores = n_cores)

  # Splits by rows; report the NA count of each quantity, excluding the counts
  acc_mat <- do.call(rbind, res)
  n_na <- colSums(is.na(acc_mat[, !colnames(acc_mat) %in% c("n_try", "n_fail"),
                                drop = FALSE]))
  if (any(n_na > 0)) {

    cat("  NA/", M, ": ", paste(names(n_na)[n_na > 0], n_na[n_na > 0],
                                sep = "=", collapse = ", "), "\n", sep = "",
        file = stderr())

  }
  n_fail <- sum(acc_mat[, "n_fail"], na.rm = TRUE)
  n_try <- sum(acc_mat[, "n_try"], na.rm = TRUE)
  cat("  dropped EM fits: ", n_fail, " of ", n_try, " (",
      round(100 * n_fail / n_try), "%)\n", sep = "", file = stderr())
  list(name = dataset$name, n = nrow(dataset$X), d = dataset$d,
       n_classes = nlevels(dataset$y), r_grid = r_grid, acc = acc_mat)

}

# Each dataset is cached as soon as it finishes (cv-hd-<dataset>-M<M>.RData),
# so an interrupted run resumes and earlier M's stay reloadable
datasets <- list(load_hydrochem(), load_vowel(), load_letter())
res <- lapply(datasets, function(dataset) {

  cache_file <- sprintf("cv-hd-%s-M%d.RData", dataset$name, M)
  if (!file.exists(cache_file)) {

    res_dataset <- run_dataset(dataset)
    save(res_dataset, n_runs_em, file = cache_file)

  }
  load(cache_file)
  res_dataset

})

}

## Paper numbers and figure
{

# Mean over the splits, NA unless at least half of them are available
mean_ok <- function(v) {

  if (sum(!is.na(v)) >= M / 2) mean(v, na.rm = TRUE) else NA

}

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

  # Best mean accuracy of movMF along the r grid
  mean_acc_r <- apply(acc_mat[, paste0("movMF", res_dataset$r_grid),
                              drop = FALSE], 2, mean_ok)
  cat("  movMF: max ", sprintf("%.1f", max(mean_acc_r, na.rm = TRUE)),
      " at r = ", res_dataset$r_grid[which.max(mean_acc_r)], "\n", sep = "")

}

# Figure: one row per (dataset, method, r) with the mean accuracy across the
# splits and its standard deviation
acc_df <- do.call(rbind, lapply(res, function(res_dataset) {

  # Summary of the given accuracy columns, one per r
  stat_cols <- function(cols, method) {

    acc_cols <- res_dataset$acc[, cols, drop = FALSE]
    acc <- apply(acc_cols, 2, mean_ok)
    data.frame(dataset = res_dataset$name, r = res_dataset$r_grid,
               method = method, acc = acc,
               sd = ifelse(is.na(acc), NA,
                           apply(acc_cols, 2, sd, na.rm = TRUE)),
               row.names = NULL)

  }

  # The kde methods do not depend on r: their column is recycled along the grid
  n_r <- length(res_dataset$r_grid)
  rbind(stat_cols(rep("CV_com", n_r), "kde-CV"),
        stat_cols(rep("ROT_com", n_r), "kde-ROT"),
        stat_cols(paste0("movMF", res_dataset$r_grid), "movMF"))

}))

# Panel and legend orders
acc_df$dataset <- factor(acc_df$dataset, levels = sapply(res, `[[`, "name"))
acc_df$method <- factor(acc_df$method,
                        levels = c("kde-CV", "kde-ROT", "movMF"))

# Palette of Figure 1, with the r-free kdes dashed and movMF solid
col_method <- c("kde-CV" = "green3", "kde-ROT" = "blue", "movMF" = "red")
lty_method <- c("kde-CV" = "dashed", "kde-ROT" = "dashed", "movMF" = "solid")

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
         y = "Test accuracy", color = "", linetype = "") +
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
