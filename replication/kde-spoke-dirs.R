
# This script produces Figure 3 in the paper and Figures 13a-13b in the SM.

## Preprocessing

# Required libraries
library(foreach)
library(progressr)
library(parallel)
library(optimParallel)
library(polykde)
library(rgl)

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
stopifnot(ncol(X) == (sum(d) + r))
stopifnot(all(abs(rowSums(X^2) - r) < 1e-10))

## ROT bandwidth selection

# Plug-in
kernel_type <- 2
kernel <- 3
k_sfp <- 100
bwd_mrot <- bw_mrot_polysph(X = X, d = d, kernel = kernel, k = k_sfp)
bwd_rot <- bw_rot_polysph(X = X, d = d, bw0 = 2^(-1:2) %o% bwd_mrot,
                          kernel = kernel, kernel_type = kernel_type,
                          k = k_sfp)

## Analysis of the cross-validated depths (avoid counting the contribution of
## a point to its own depth)

# Log-cross-validated kde evaluated at each point
h <- bwd_rot$bw
log_kde <- log_cv_kde_polysph(X = X, d = d, h = h, wrt_unif = TRUE,
                              kernel = kernel, kernel_type = kernel_type,
                              k = k_sfp)
sum(duplicated(log_kde))
sum(!is.finite(log_kde))

# All subjects without ordering (for illustration purposes)
open3d()
P <- matrix(c(-0.01408212, 0.6652657, -0.7464739, 0,
               0.80961037, 0.4457031,  0.3819424, 0,
               0.58679885,-0.5989745, -0.5448825, 0,
               0.00000000, 0.00000000, 0.00000000, 1),
            nrow = 4, ncol = 4, byrow = TRUE)
par3d(windowRect = c(0, 45, 1440, 898), zoom = 0.5050682, userMatrix = P)
mfrow3d(nr = 12, nc = 15, sharedMouse = TRUE)
for (i in seq_len(nrow(dirs))) {

  # Boundary points
  x <- bdry[i, , ]
  plot3d(x, col = "lightblue", axes = FALSE, box = FALSE,
         xlab = "", ylab = "", zlab = "",
         aspect = apply(x, 2, function(x) diff(range(x))))

  # Magically build the boundary
  ash <- alphashape3d::ashape3d(x, alpha = 6)
  m <- as.mesh3d(ash)
  shade3d(m, alpha = 0.5, col = "lightblue")
  addNormals(x = m)

}
snapshot3d("hippo_unsorted.png", top = FALSE, webshot = FALSE,
           width = 1440, height = 853)

# Rank all subjects from the highest depth to the lowest depth
col <- viridis::viridis(nrow(dirs), direction = -1)
open3d()
P <- matrix(c(-0.01408212,  0.6652657, -0.7464739, 0,
               0.80961037,  0.4457031,  0.3819424, 0,
               0.58679885, -0.5989745, -0.5448825, 0,
               0.00000000, 0.00000000, 0.00000000, 1),
            nrow = 4, ncol = 4, byrow = TRUE)
par3d(windowRect = c(0, 45, 1440, 898), zoom = 0.5050682, userMatrix = P)
mfrow3d(nr = 12, nc = 15, sharedMouse = TRUE)
inds <- order(log_kde, decreasing = TRUE)
for (i in seq_len(nrow(dirs))) {

  # Boundary points
  x <- bdry[inds[i], , ]
  plot3d(x, col = col[i], axes = FALSE, box = FALSE,
         xlab = "", ylab = "", zlab = "",
         aspect = apply(x, 2, function(x) diff(range(x))))

  # Magically build the boundary
  ash <- alphashape3d::ashape3d(x, alpha = 6)
  m <- as.mesh3d(ash)
  shade3d(m, alpha = 0.5, col = col[i])
  addNormals(x = m)

}
snapshot3d("hippo_sorted.png", top = FALSE, webshot = FALSE,
           width = 1440, height = 853)

# Rank all subjects from the highest depth to the lowest depth, colored
# for minority class -- Figure 3
col <- viridis::viridis(nrow(dirs), direction = -1)
open3d()
P <- matrix(c(-0.01408212,  0.6652657, -0.7464739, 0,
               0.80961037,  0.4457031,  0.3819424, 0,
               0.58679885, -0.5989745, -0.5448825, 0,
               0.00000000, 0.00000000, 0.00000000, 1),
            nrow = 4, ncol = 4, byrow = TRUE)
par3d(windowRect = c(0, 45, 1440, 898), zoom = 0.51, userMatrix = P)
mfrow3d(nr = 12, nc = 15, sharedMouse = TRUE)
inds <- order(log_kde, decreasing = TRUE)
for (i in seq_len(nrow(dirs))) {

  # Boundary points
  x <- bdry[inds[i], , ]
  col_i <- col[i] # ifelse(ids_labs[inds[i]], "red", col[i])
  plot3d(x, col = col_i,
         # axes = TRUE, box = TRUE, xlab = "x", ylab = "y", zlab = "z",
         axes = FALSE, box = FALSE, xlab = "", ylab = "", zlab = "",
         aspect = apply(x, 2, function(x) diff(range(x))))

  # Magically build the boundary
  ash <- alphashape3d::ashape3d(x, alpha = 6)
  m <- as.mesh3d(ash)
  shade3d(m, alpha = 0.5, col = col_i)
  addNormals(x = m)

  # Add a red point
  if (ids_labs[inds[i]]) {

    points3d(x = -125, y = -107.5, z = 49.5, radius = 1, col = "red", size = 7)

  }

}
snapshot3d("hippo_sorted_red.png", top = FALSE, webshot = FALSE,
           width = 1440, height = 853)

# 9 most outlying subjects -- Figure 13b
col <- viridis::viridis(nrow(dirs))
open3d()
P <- matrix(c(-0.02513775,  0.7123916, -0.7013319, 0,
               0.91663814,  0.2963825,  0.2682016, 0,
               0.39892703, -0.6361256, -0.6604558, 0,
               0.00000000,  0.0000000,  0.0000000, 1),
            nrow = 4, ncol = 4, byrow = TRUE)
par3d(windowRect = c(0, 45, 1440, 898), zoom = 0.5050682, userMatrix = P)
mfrow3d(nr = 3, nc = 3, sharedMouse = TRUE)
inds <- order(log_kde)
for (i in 9:1) {

  # Boundary points
  x <- bdry[inds[i], , ]
  plot3d(x, col = col[i], axes = FALSE, box = FALSE, xlab = "", ylab = "",
         zlab = "", aspect = apply(x, 2, function(x) diff(range(x))))

  # Magically build the boundary
  ash <- alphashape3d::ashape3d(x, alpha = 6)
  m <- as.mesh3d(ash)
  shade3d(m, alpha = 0.5, col = col[i])
  addNormals(x = m)

}
snapshot3d("hippo_most_out.png", top = FALSE, webshot = FALSE,
           width = 1440, height = 853)

# 9 most deep subjects -- Figure 13a
col <- viridis::viridis(nrow(dirs), direction = -1)
open3d()
P <- matrix(c(-0.02513775,  0.7123916, -0.7013319, 0,
               0.91663814,  0.2963825,  0.2682016, 0,
               0.39892703, -0.6361256, -0.6604558, 0,
               0.00000000,  0.0000000,  0.0000000, 1),
            nrow = 4, ncol = 4, byrow = TRUE)
par3d(windowRect = c(0, 45, 1440, 898), zoom = 0.5050682, userMatrix = P)
mfrow3d(nr = 3, nc = 3, sharedMouse = TRUE)
inds <- order(log_kde, decreasing = TRUE)
for (i in 9:1) {

  # Boundary points
  x <- bdry[inds[i], , ]
  plot3d(x, col = col[i], axes = FALSE, box = FALSE,
         xlab = "", ylab = "", zlab = "",
         aspect = apply(x, 2, function(x) diff(range(x))))

  # Magically build the boundary
  ash <- alphashape3d::ashape3d(x, alpha = 6)
  m <- as.mesh3d(ash)
  shade3d(m, alpha = 0.5, col = col[i])
  addNormals(x = m)

}
snapshot3d("hippo_most_in.png", top = FALSE, webshot = FALSE,
           width = 1440, height = 853)
