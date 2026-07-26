
# Hex sticker for the 'polykde' package: a cluster of three S^2 globes, each
# painted with a multimodal von Mises-Fisher mixture density (a KDE on the
# sphere) via polykde's own d_vmf_polysph. The globes are drawn with a base-
# graphics orthographic ray-caster (no rgl needed) and framed in a dark-viridis
# hexagon. An optional rgl path (render_rgl) is kept at the bottom.

library(polykde)
library(viridis)
library(hexSticker)

# Shared logo standards
font <- "Aller_Rg"
name_size <- 31.2
url_size <- 9.0
url_x <- 1.00
url_y <- 0.08
url_angle <- 30
hex_border <- 1.5
dpi <- 600
set.seed(20251110)

# Output directories
dir.create("logo", showWarnings = FALSE)
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

# Sphere ray-caster

# Evaluate the vMF mixture density at rows of x (unit vectors on S^2).
mix_dens <- function(x, mu, kappa, prop) {
  mu <- rbind(mu)
  prop <- prop / sum(prop)
  d <- numeric(nrow(x))
  for (k in seq_len(nrow(mu))) {
    d <- d + prop[k] *
      d_vmf_polysph(x = x, d = 2, mu = mu[k, ], kappa = kappa[k])
  }
  d
}

# Render one translucent "crystal" globe as an npix x npix "#RRGGBBAA" raster:
# the vMF mixture density on the front hemisphere as filled viridis level bands
# with thin isolines, over glassy shading (flat fill, specular glint, Fresnel
# rim). Semi-transparent (more opaque at the rim) so overlapping globes blend.
sphere_raster <- function(mu, kappa, prop, npix = 620, nlev = 9,
                          base_alpha = 0.60, light = c(-0.5, 0.62, 0.61),
                          gamma = 0.80, pal = viridis(256),
                          iso_w = 0.10, iso_dark = 0.45,
                          ambient = 0.72, diffuse = 0.34,
                          spec_str = 0.65, spec_pow = 34,
                          rim_str = 0.55, rim_pow = 3.2) {

  light <- light / sqrt(sum(light^2))

  # Pixel grid over [-1, 1]^2; rows go top (+y) to bottom (-y).
  grid_x <- seq(-1, 1, length.out = npix)
  grid_y <- seq(1, -1, length.out = npix)
  px_x <- matrix(grid_x, npix, npix, byrow = TRUE)
  px_y <- matrix(grid_y, npix, npix, byrow = FALSE)
  rr <- sqrt(px_x^2 + px_y^2)
  inside <- rr <= 1

  # Front-hemisphere surface points = unit normals on S^2 (view space).
  px_z <- sqrt(pmax(0, 1 - px_x^2 - px_y^2))
  nx <- px_x[inside]; ny <- px_y[inside]; nz <- px_z[inside]
  pts <- cbind(nx, ny, nz)

  # KDE contour plot: quantise the density into filled viridis level bands.
  dens <- mix_dens(pts, mu, kappa, prop)
  s <- ((dens - min(dens)) / (max(dens) - min(dens) + 1e-12))^gamma
  band <- pmin(nlev - 1, floor(s * nlev))
  band_col <- pal[pmax(1, pmin(256, 1 + round((band + 0.5) / nlev * 255)))]
  base_rgb <- col2rgb(band_col) # 3 x K

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

  red <- pmin(255, round(base_rgb[1, ] * shade * iso_fac +
                           255 * spec + 205 * rim))
  green <- pmin(255, round(base_rgb[2, ] * shade * iso_fac +
                             255 * spec + 232 * rim))
  blue <- pmin(255, round(base_rgb[3, ] * shade * iso_fac +
                            255 * spec + 250 * rim))

  # Translucent interior, more opaque at the rim (glass edge), anti-aliased.
  a <- pmin(1, base_alpha + (1 - base_alpha) * (1 - nz)^2)
  aa <- pmin(1, pmax(0, (1 - rr[inside]) / (2.5 / npix)))
  alpha <- round(255 * a * aa)

  cols <- character(npix * npix)
  cols[] <- "#00000000"
  cols[which(inside)] <- rgb(red, green, blue, alpha, maxColorValue = 255)
  matrix(cols, npix, npix)

}

# Compose the three-globe cluster

# Unit-norm helper for hand-placed modes.
u3 <- function(...) { v <- c(...); v / sqrt(sum(v^2)) }

# Three globes forming a triangular cluster (a "product of spheres"): a larger
# front-centre globe and two behind it. Modes are hand-placed on the visible
# (z > 0) face so each globe shows viridis hotspots.
globes <- list(
  # front, centre
  list(cx = 0.50, cy = 0.44, radius = 0.255,
       mu = rbind(u3(0.15, 0.25, 0.95), u3(-0.55, -0.25, 0.8),
                  u3(0.55, -0.30, 0.78)),
       kappa = c(11, 15, 12), prop = c(0.4, 0.32, 0.28)),
  # upper-left
  list(cx = 0.33, cy = 0.62, radius = 0.225,
       mu = rbind(u3(-0.15, 0.35, 0.92), u3(0.5, -0.2, 0.84)),
       kappa = c(12, 10), prop = c(0.55, 0.45)),
  # upper-right
  list(cx = 0.67, cy = 0.62, radius = 0.225,
       mu = rbind(u3(0.1, 0.15, 0.98), u3(-0.5, 0.3, 0.81),
                  u3(0.35, -0.5, 0.79)),
       kappa = c(13, 11, 16), prop = c(0.4, 0.33, 0.27))
)

# Draw order: back (upper) globes first, front globe last so it overlaps.
draw_order <- c(2, 3, 1)

# Render the three globes onto a transparent canvas
canvas_px <- 1600
subplot_png <- tempfile(fileext = ".png") # interior render (not shipped)
png(subplot_png, width = canvas_px, height = canvas_px, bg = "transparent")
op <- par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
for (i in draw_order) {
  g <- globes[[i]]
  npix <- min(660, round(2 * g$radius * canvas_px))
  ras <- sphere_raster(mu = g$mu, kappa = g$kappa, prop = g$prop, npix = npix)
  rasterImage(as.raster(ras),
              xleft = g$cx - g$radius, ybottom = g$cy - g$radius,
              xright = g$cx + g$radius, ytop = g$cy + g$radius,
              interpolate = TRUE)
}
par(op)
dev.off()

# Assemble the hex sticker

# Dark-viridis palette.
hex_fill <- "#5B4A87" # lighter viridis-violet
hex_border_col <- "#48C16E" # bright viridis green
txt_col <- "#FFFFFF"
url_col <- "#B9A9D6"

# Draw the hexagon, wordmark and URL, then mirror to man/figures
sticker(
  subplot = subplot_png,
  s_x = 1.00, s_y = 0.98, s_width = 1.04, s_height = 1.04,
  package = "polykde", p_x = 1.00, p_y = 0.44, p_size = name_size,
  p_color = txt_col, p_family = font,
  h_fill = hex_fill, h_color = hex_border_col, h_size = hex_border,
  spotlight = FALSE,
  url = "github.com/egarpor/polykde",
  u_x = url_x, u_y = url_y, u_angle = url_angle, u_size = url_size,
  u_color = url_col, u_family = font,
  dpi = dpi, filename = "logo/logo.png"
)
file.copy("logo/logo.png", "man/figures/logo.png", overwrite = TRUE)

# Optional: native rgl render (interactive only)
# Run in RStudio / an interactive session with a working OpenGL context to get
# the rgl-native version (lit translucent spheres). Not used by default because
# rgl cannot open a GL window under `Rscript` here.
render_rgl <- function(file = "logo-spheres-rgl.png") {
  stopifnot(requireNamespace("rgl", quietly = TRUE))
  rgl::open3d(useNULL = FALSE)
  rgl::bg3d("white")
  centers <- list(c(0, 0, 0), c(-1.6, 1.5, -0.6), c(1.6, 1.5, -0.6))
  radii <- c(1, 0.72, 0.72)
  for (j in seq_along(centers)) {
    g <- globes[[j]]
    latt <- fib_latt(n = 8000)
    dens <- numeric(nrow(latt))
    for (k in seq_len(nrow(rbind(g$mu)))) {
      dens <- dens + g$prop[k] *
        d_vmf_polysph(x = latt, d = 2, mu = rbind(g$mu)[k, ],
                      kappa = g$kappa[k])
    }
    cols <- viridis(256)[1 + round(255 * (dens - min(dens)) /
                                     diff(range(dens)))]
    pts <- latt * radii[j]
    pts <- sweep(pts, 2, centers[[j]], `+`)
    rgl::spheres3d(centers[[j]], radius = radii[j] * 0.985,
                   col = "grey85", alpha = 0.25, lit = TRUE)
    rgl::plot3d(pts, col = cols, size = 5, add = TRUE,
                axes = FALSE, box = FALSE, xlab = "", ylab = "", zlab = "")
  }
  rgl::par3d(zoom = 0.62)
  rgl::snapshot3d(file, width = 1600, height = 1600, webshot = FALSE)
  rgl::close3d()
}
