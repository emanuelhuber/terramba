test_that("interp2DMBA reproduces direct mba.surf output", {
  skip_if_not_installed("MBA")

  set.seed(42)
  xyz <- data.frame(x = runif(100, 0, 80), y = runif(100, 0, 50))
  xyz$z <- with(xyz, sin(x / 10) + cos(y / 8) + 0.01 * x)

  # Distinct resolutions and distinct m/n values protect the established,
  # deliberate grid geometry and m/n mapping from accidental changes.
  dx <- 2
  dy <- 4
  m <- 2
  n <- 1
  xlim <- range(xyz$x)
  ylim <- range(xyz$y)

  r_wrapper <- interp2DMBA(
    xyz, dx = dx, dy = dy, xlim = xlim, ylim = ylim,
    m = m, n = n, h = 6, buffer = 0, trend = NULL
  )

  para <- terramba:::.getbbox_nx_ny(
    xpos = seq(xlim[1], xlim[2], by = dx),
    ypos = seq(ylim[1], ylim[2], by = dy),
    dx = dx, dy = dy, buffer = 0
  )
  direct <- MBA::mba.surf(
    xyz, para$nx, para$ny,
    n = m, m = n, extend = TRUE, h = 6, b.box = para$bbox
  )$xyz.est
  z_matrix <- matrix(rev(as.vector(direct$z)), nrow = para$nx, ncol = para$ny)
  r_direct <- terra::rast(
    ncols = para$nx, nrows = para$ny,
    xmin = para$bbox[1], xmax = para$bbox[2] + dx,
    ymin = para$bbox[3], ymax = para$bbox[4] + dx
  )
  terra::values(r_direct) <- t(z_matrix[nrow(z_matrix):1, ])

  expect_true(terra::compareGeom(r_wrapper, r_direct, stopOnError = FALSE))
  expect_equal(
    terra::values(r_wrapper, mat = FALSE),
    terra::values(r_direct, mat = FALSE),
    tolerance = 1e-12
  )
})
