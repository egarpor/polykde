
# Oracle-K sweep for the movMF competitor in the real-data application
# (Section 5): test accuracy of the movMF classifier at EVERY admissible
# number of components K, over the same R stratified splits as the main
# pipeline. "Admissible" means one parameter per observation within each
# class, K <= n_c / (d + 2), the entire range where the mixture is not
# overparameterized; each class uses min(K, its own cap) components. The sweep
# shows the accuracy-optimal K is interior to the explored range, so the
# ceilings of cv-hd-application.R do not bind for the comparison with the kde.
# Writes movmf-k-sweep.RData (cache) and paper/img/movmf_k_sweep.pdf.

# The main pipeline's cache must exist: sourcing the script below would
# otherwise trigger its full recompute inside this process
stopifnot(file.exists(file.path(Sys.getenv("HOME"), ".cvapp-cache",
                                "cv-hd.RData")))
source("/Users/Eduardo/GitHub/polykde/polykde/replication/cv-hd-application.R")

# Lift the pooled-path ceiling: the sweep explores the parameter cap alone
Kmax <- Inf

# Same split + embedding as one_split()
split_embed <- function(ds, seed, prop = 0.7) {

  X <- ds$X
  y <- ds$y
  cls <- levels(y)
  set.seed(seed)
  tr <- unlist(lapply(cls, function(l) {

    i <- which(y == l)
    sample(i, max(2L, round(prop * length(i))))

  }))
  ytr <- droplevels(y[tr])
  yte <- y[-tr]
  if (ds$embed == "sqrt") {

    Ztr <- sqrt_map(X[tr, , drop = FALSE])
    Zte <- sqrt_map(X[-tr, , drop = FALSE])

  } else {

    ctr <- colMeans(X[tr, , drop = FALSE])
    scl <- apply(X[tr, , drop = FALSE], 2, sd)
    Ztr <- l2_map(X[tr, , drop = FALSE], ctr, scl)
    Zte <- l2_map(X[-tr, , drop = FALSE], ctr, scl)

  }
  list(Ztr = Ztr, ytr = ytr, Zte = Zte, yte = yte)

}

# One split: accuracy of the movMF classifier at every shared K = 1, ..., K_top
# (each class truncated at its own parameter cap), plus the accuracy at the
# per-class BIC selection (the cls variant of the main pipeline)
sweep_split <- function(ds, seed, K_top) {

  sp <- split_embed(ds, seed)
  cls <- levels(sp$ytr)
  logprior <- log(as.numeric(table(sp$ytr)) / length(sp$ytr))
  err <- function(ld) mean(classify(ld, logprior, cls) != sp$yte)

  # Full path and test log-densities for each class, up to its own cap
  paths <- lapply(cls, function(l) {

    p <- mvmf_path(sp$Ztr[sp$ytr == l, , drop = FALSE])
    ld <- sapply(p$fits, function(f)
      dmovMF(sp$Zte, theta = f$theta, alpha = f$alpha, log = TRUE))
    list(k = p$k, ld = ld, bic = sapply(p$fits, BIC))

  })

  # Assemble the classifier at each shared K by indexing the stored densities:
  # each class uses its largest fitted k not exceeding K
  acc_K <- sapply(seq_len(K_top), function(K) {

    ld <- sapply(paths, function(p) p$ld[, which.max(p$k * (p$k <= K))])
    1 - err(ld)

  })

  # Per-class BIC selection (cls variant)
  ld_bic <- sapply(paths, function(p) p$ld[, which.min(p$bic)])
  c(acc_K, bic_cls = 1 - err(ld_bic))

}

## Sweep

sweep_file <- file.path(cache_dir, "movmf-k-sweep.RData")
if (!file.exists(sweep_file)) {

  datasets <- list(dedup_ds(load_hydrochem()), dedup_ds(load_vowel()),
                   dedup_ds(load_letter()))
  sweeps <- lapply(datasets, function(ds) {

    # Largest admissible K across classes: the x-range of this dataset's curve
    n_cls <- table(ds$y)
    K_top <- floor(0.7 * max(n_cls) / (ncol(ds$X) + 1))
    cat(ds$name, "K_top =", K_top, "\n", file = stderr())
    acc <- mclapply(seq_len(R), function(s)
      tryCatch(sweep_split(ds, seed = s, K_top = K_top),
               error = function(e) rep(NA_real_, K_top + 1)),
      mc.cores = n_cores)
    list(name = ds$name, d = ds$d, K_top = K_top, acc = do.call(rbind, acc))

  })
  save(sweeps, R, nruns, file = sweep_file)

}

## Figure: accuracy vs K, against the kde-CV reference

load(sweep_file)
load(results_file)

# kde-CV (com) accuracy at each dataset's full dimension, as reference line
acc_cv <- sapply(list(at_d(hydro, 13), at_d(vowel, 8), at_d(letter, 15)),
                 function(z) mean(1 - z$err[, "CV_com"], na.rm = TRUE))

library(ggplot2)
df <- do.call(rbind, Map(function(sw, cv) {

  a <- sw$acc[, seq_len(sw$K_top), drop = FALSE]
  data.frame(dataset = sw$name, K = seq_len(sw$K_top),
             acc = colMeans(a, na.rm = TRUE),
             sd = apply(a, 2, sd, na.rm = TRUE), cv = cv)

}, sweeps, acc_cv))
df$dataset <- factor(df$dataset, levels = sapply(sweeps, `[[`, "name"))
gg <- ggplot(df, aes(x = K, y = acc)) +
  geom_ribbon(aes(ymin = acc - sd, ymax = acc + sd), alpha = 0.2) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  geom_hline(aes(yintercept = cv), linetype = "dashed") +
  facet_wrap(~ dataset, nrow = 1, scales = "free") +
  scale_x_continuous(breaks = function(x) seq_len(floor(x[2]))) +
  labs(x = expression("Shared number of components" ~ K),
       y = "Test accuracy") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank())
ggsave(file.path(paper_dir, "img", "movmf_k_sweep.pdf"), plot = gg,
       device = cairo_pdf, width = 12, height = 4, units = "in")
