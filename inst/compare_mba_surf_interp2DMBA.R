# Compare MBA::mba.surf() with terramba::interp2DMBA()
#
# The comparison uses trend = NULL and m = n = 1 so that both functions
# receive the same observations and interpolation settings. The direct MBA
# result is converted to a terra::SpatRaster using the same orientation
# conversion as interp2DMBA().

library(MBA)
library(terra)
library(terramba)

set.seed(42)

# Irregular sample points from a known synthetic surface
n_points <- 120
xyz <- data.frame(
  x = runif(n_points, 0, 100),
  y = runif(n_points, 0, 70)
)
xyz$z <- with(
  xyz,
  0.04 * x - 0.03 * y +
    2.5 * sin(x / 14) * cos(y / 11) +
    rnorm(n_points, sd = 0.15)
)

# Identical settings for both approaches
dx <- 2
dy <- 2
buffer <- 0
h <- 6
m <- 1
n <- 1
xlim <- range(xyz$x)
ylim <- range(xyz$y)

# 1. terramba wrapper
r_terramba <- interp2DMBA(
  x = xyz,
  dx = dx,
  dy = dy,
  xlim = xlim,
  ylim = ylim,
  buffer = buffer,
  h = h,
  m = m,
  n = n,
  trend = NULL,
  output = "raster"
)

# 2. Direct MBA::mba.surf() call
# Reproduce the grid definition used by interp2DMBA().
para <- terramba:::.getbbox_nx_ny(
  xpos = seq(xlim[1], xlim[2], by = dx),
  ypos = seq(ylim[1], ylim[2], by = dy),
  dx = dx,
  dy = dy,
  buffer = buffer
)

mba_result <- MBA::mba.surf(
  xyz = xyz,
  no.X = para$nx,
  no.Y = para$ny,
  # Retain the deliberate mapping used inside interp2DMBA().
  n = m,
  m = n,
  extend = TRUE,
  h = h,
  b.box = para$bbox
)

# Convert the MBA image-format result exactly as interp2DMBA() does.
S  <- mba_result$xyz.est
SZ <- matrix(
  rev(as.vector(S$z)),
  nrow = para$nx,
  ncol = para$ny
)

r_mba <- terra::rast(
  ncols = para$nx,
  nrows = para$ny,
  xmin = para$bbox[1],
  xmax = para$bbox[2] + dx,
  ymin = para$bbox[3],
  # Retain the deliberate geometry used inside interp2DMBA().
  ymax = para$bbox[4] + dx
)
terra::values(r_mba) <- t(SZ[rev(seq_len(nrow(SZ))), , drop = FALSE])

# 3. Numerical comparison
r_difference <- r_terramba - r_mba
v1 <- terra::values(r_terramba, mat = FALSE)
v2 <- terra::values(r_mba, mat = FALSE)
valid <- is.finite(v1) & is.finite(v2)

comparison <- c(
  n_compared = sum(valid),
  mean_difference = mean(v1[valid] - v2[valid]),
  rmse = sqrt(mean((v1[valid] - v2[valid])^2)),
  maximum_absolute_difference = max(abs(v1[valid] - v2[valid]))
)

print(comparison)
print(terra::compareGeom(r_terramba, r_mba, stopOnError = FALSE))

# The differences should be zero, apart from possible floating-point noise.
stopifnot(
  isTRUE(terra::compareGeom(r_terramba, r_mba, stopOnError = FALSE)),
  comparison[["maximum_absolute_difference"]] < 1e-12
)

# 4. Visual comparison
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(1, 3), mar = c(3, 3, 3, 5))

plot(r_mba, main = "Direct MBA::mba.surf()")
points(xyz$x, xyz$y, pch = 16, cex = 0.35)

plot(r_terramba, main = "terramba::interp2DMBA()")
points(xyz$x, xyz$y, pch = 16, cex = 0.35)

plot(r_difference, main = "Difference")


