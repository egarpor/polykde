# Generates the polykde hexagonal sticker logo.
#
# Concept: a cluster of three S^2 globes, each painted with a multimodal
# von Mises-Fisher mixture density (a "KDE on the sphere"), coloured with
# viridis and framed in a dark-viridis hexagon. The interior is rendered with
# polykde's own density evaluator (d_vmf_polysph), so the logo literally shows
# what the package computes.
#
# Rendering note: the paper's sphere figures use rgl + snapshot3d, but rgl
# needs an interactive OpenGL context (RStudio / a real X11 GL window) that is
# not available under plain `Rscript`. To keep this script fully reproducible
# (no GL, no headless browser), the globes are drawn with a small base-graphics
# orthographic ray-caster: for every pixel of a sphere's front hemisphere we
# recover the surface normal, evaluate the mixture density there, map it to
# viridis and apply Lambert shading. An optional rgl path (render_rgl(), not
# run by default) is kept at the bottom for interactive use.
#
# Run from the package root:
#   Rscript logo/logo.R
# Output: logo/logo.png (master) and man/figures/logo.png (shipped mirror).
# The transparent interior render (the hex subplot) is written to a tempfile so
# only logo.png ships in the package.

library(polykde)
library(viridis)
library(hexSticker)

## ---- Shared logo standard (identical across egarpor packages) -------------
# Aller_Rg is bundled with (and auto-registered by) hexSticker, so no
# showtext/font_add setup is needed for the sticker() idiom.
FONT    <- "Aller_Rg"   # typeface for the package name and the GitHub URL
P_SIZE  <- 31.2         # package-name size (shared across packages)
U_SIZE  <- 9.0          # GitHub URL size (large enough to read)
U_X     <- 1.00         # GitHub URL position: along the lower-right hex edge
U_Y     <- 0.08
U_ANGLE <- 30
H_SIZE  <- 1.5          # hexagon border thickness
DPI     <- 600

set.seed(20251110)

dir.create("logo", showWarnings = FALSE)
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Sphere ray-caster ------------------------------------------------ ##

# Renders one S^2 globe as an npix x npix matrix of "#RRGGBBAA" colours, with
# the front hemisphere shaded and coloured by a vMF mixture density. Pixels
# outside the disk (and the anti-aliased rim) are transparent.
#
#   mu    : m x 3 matrix of unit-norm modal directions (view space; +z faces
#           the camera, so modes with z > 0 are visible).
#   kappa : length-m concentrations.
#   prop  : length-m mixture weights (need not sum to one; rescaled here).
#   light : length-3 light direction (view space), pointing toward the light.
#   gamma : tone curve on the density (<1 lifts mid densities).
# Evaluate the vMF mixture density at rows of `x` (unit vectors on S^2).
mix_dens <- function(x, mu, kappa, prop) {
  mu <- rbind(mu)
  prop <- prop / sum(prop)
  d <- numeric(nrow(x))
  for (k in seq_len(nrow(mu))) {
    d <- d + prop[k] * d_vmf_polysph(x = x, d = 2, mu = mu[k, ], kappa = kappa[k])
  }
  d
}

# Renders one translucent "crystal" globe carrying a KDE contour plot: the vMF
# mixture density on the front hemisphere is drawn as filled viridis level bands
# with thin isolines between them, over a glassy shading (mostly flat fill, a
# specular highlight and a Fresnel rim). The globe is semi-transparent (alpha <
# 1, more opaque at the rim) so overlapping globes blend like glass.
sphere_raster <- function(mu, kappa, prop, npix = 620, nlev = 9,
                          base_alpha = 0.60, light = c(-0.5, 0.62, 0.61),
                          gamma = 0.80, pal = viridis(256),
                          iso_w = 0.10, iso_dark = 0.45,
                          ambient = 0.72, diffuse = 0.34,
                          spec_str = 0.65, spec_pow = 34,
                          rim_str = 0.55, rim_pow = 3.2) {

  light <- light / sqrt(sum(light^2))

  # Pixel grid over [-1, 1]^2; rows go top (+y) to bottom (-y).
  u <- seq(-1, 1, length.out = npix)
  v <- seq(1, -1, length.out = npix)
  U <- matrix(u, npix, npix, byrow = TRUE)
  V <- matrix(v, npix, npix, byrow = FALSE)
  rr <- sqrt(U^2 + V^2)
  inside <- rr <= 1

  # Front-hemisphere surface points = unit normals on S^2 (view space).
  Z <- sqrt(pmax(0, 1 - U^2 - V^2))
  nx <- U[inside]; ny <- V[inside]; nz <- Z[inside]
  pts <- cbind(nx, ny, nz)

  # KDE contour plot: quantise the density into filled viridis level bands.
  dens <- mix_dens(pts, mu, kappa, prop)
  s <- ((dens - min(dens)) / (max(dens) - min(dens) + 1e-12))^gamma
  band <- pmin(nlev - 1, floor(s * nlev))
  band_col <- pal[pmax(1, pmin(256, 1 + round((band + 0.5) / nlev * 255)))]
  base_rgb <- col2rgb(band_col)                      # 3 x K

  # Thin isolines at the band boundaries.
  iso <- pmin(1, abs(s * nlev - round(s * nlev)) / iso_w)
  iso_fac <- iso_dark + (1 - iso_dark) * iso

  # Glassy shading: mostly flat (so bands read), a specular glint and a Fresnel
  # rim highlight toward cyan-white.
  lam <- pmax(0, nx * light[1] + ny * light[2] + nz * light[3])
  shade <- pmin(1, ambient + diffuse * lam)
  rz <- 2 * lam * nz - light[3]
  spec <- spec_str * pmax(0, rz)^spec_pow
  rim <- rim_str * (1 - nz)^rim_pow

  R <- pmin(255, round(base_rgb[1, ] * shade * iso_fac + 255 * spec + 205 * rim))
  G <- pmin(255, round(base_rgb[2, ] * shade * iso_fac + 255 * spec + 232 * rim))
  B <- pmin(255, round(base_rgb[3, ] * shade * iso_fac + 255 * spec + 250 * rim))

  # Translucent interior, more opaque at the rim (glass edge), anti-aliased.
  a <- pmin(1, base_alpha + (1 - base_alpha) * (1 - nz)^2)
  aa <- pmin(1, pmax(0, (1 - rr[inside]) / (2.5 / npix)))
  A <- round(255 * a * aa)

  cols <- character(npix * npix)
  cols[] <- "#00000000"
  cols[which(inside)] <- rgb(R, G, B, A, maxColorValue = 255)
  matrix(cols, npix, npix)
}

## ---- 2. Compose the three-globe cluster ---------------------------------- ##

# Unit-norm helper for hand-placed modes.
u3 <- function(...) { v <- c(...); v / sqrt(sum(v^2)) }

# Three globes forming a triangular cluster that reads as a product of spheres:
# a slightly larger front-centre one and two behind it. Modal directions are
# hand-placed on the visible (z > 0) face so each globe shows viridis hotspots.
globes <- list(
  # front, centre
  list(cx = 0.50, cy = 0.44, R = 0.255,
       mu = rbind(u3(0.15, 0.25, 0.95), u3(-0.55, -0.25, 0.8),
                  u3(0.55, -0.30, 0.78)),
       kappa = c(11, 15, 12), prop = c(0.4, 0.32, 0.28)),
  # upper-left
  list(cx = 0.33, cy = 0.62, R = 0.225,
       mu = rbind(u3(-0.15, 0.35, 0.92), u3(0.5, -0.2, 0.84)),
       kappa = c(12, 10), prop = c(0.55, 0.45)),
  # upper-right
  list(cx = 0.67, cy = 0.62, R = 0.225,
       mu = rbind(u3(0.1, 0.15, 0.98), u3(-0.5, 0.3, 0.81),
                  u3(0.35, -0.5, 0.79)),
       kappa = c(13, 11, 16), prop = c(0.4, 0.33, 0.27))
)

# Draw order: back (upper) globes first, front globe last so it overlaps.
draw_order <- c(2, 3, 1)

W <- 1600
subplot_png <- tempfile(fileext = ".png")   # interior render (not shipped)
png(subplot_png, width = W, height = W, bg = "transparent")
op <- par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
for (i in draw_order) {
  g <- globes[[i]]
  npix <- min(660, round(2 * g$R * W))
  ras <- sphere_raster(mu = g$mu, kappa = g$kappa, prop = g$prop, npix = npix)
  rasterImage(as.raster(ras),
              xleft = g$cx - g$R, ybottom = g$cy - g$R,
              xright = g$cx + g$R, ytop = g$cy + g$R,
              interpolate = TRUE)
}
par(op)
dev.off()

## ---- 3. Assemble the hex sticker ----------------------------------------- ##

# Dark-viridis palette.
hex_fill   <- "#241436"   # deep viridis-purple
hex_border <- "#48C16E"   # bright viridis green
txt_col    <- "#FFFFFF"
url_col    <- "#B9A9D6"

sticker(
  subplot  = subplot_png,
  s_x = 1.00, s_y = 0.98, s_width = 1.04, s_height = 1.04,
  package  = "polykde",
  p_x = 1.00, p_y = 0.44, p_color = txt_col, p_size = P_SIZE, p_family = FONT,
  h_fill   = hex_fill, h_color = hex_border, h_size = H_SIZE,
  spotlight = FALSE,
  url = "github.com/egarpor/polykde", u_color = url_col,
  u_x = U_X, u_y = U_Y, u_angle = U_ANGLE, u_size = U_SIZE, u_family = FONT,
  dpi = DPI,
  filename = "logo/logo.png"
)
file.copy("logo/logo.png", "man/figures/logo.png", overwrite = TRUE)

message("Wrote logo/logo.png and man/figures/logo.png")

## ---- Optional: native rgl render (interactive only) ---------------------- ##
# Run this in RStudio / an interactive session with a working OpenGL context to
# get the rgl-native version (lit translucent spheres). Not used by default
# because rgl cannot open a GL window under `Rscript` here.
render_rgl <- function(file = "logo-spheres-rgl.png") {
  stopifnot(requireNamespace("rgl", quietly = TRUE))
  rgl::open3d(useNULL = FALSE)
  rgl::bg3d("white")
  centers <- list(c(0, 0, 0), c(-1.6, 1.5, -0.6), c(1.6, 1.5, -0.6))
  radii   <- c(1, 0.72, 0.72)
  for (j in seq_along(centers)) {
    g <- globes[[j]]
    latt <- fib_latt(n = 8000)
    dens <- numeric(nrow(latt))
    for (k in seq_len(nrow(rbind(g$mu)))) {
      dens <- dens + g$prop[k] *
        d_vmf_polysph(x = latt, d = 2, mu = rbind(g$mu)[k, ], kappa = g$kappa[k])
    }
    cols <- viridis(256)[1 + round(255 * (dens - min(dens)) /
                                     diff(range(dens)))]
    P <- latt * radii[j]
    P <- sweep(P, 2, centers[[j]], `+`)
    rgl::spheres3d(centers[[j]], radius = radii[j] * 0.985,
                   col = "grey85", alpha = 0.25, lit = TRUE)
    rgl::plot3d(P, col = cols, size = 5, add = TRUE,
                axes = FALSE, box = FALSE, xlab = "", ylab = "", zlab = "")
  }
  rgl::par3d(zoom = 0.62)
  rgl::snapshot3d(file, width = 1600, height = 1600, webshot = FALSE)
  rgl::close3d()
}
