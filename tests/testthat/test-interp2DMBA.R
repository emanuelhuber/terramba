test_that("basic interpolation and matrix output work", {
  skip_if_not_installed("MBA")
  d <- expand.grid(x = seq(0, 10, by = 2), y = seq(0, 8, by = 2))
  d$z <- 2 + 0.5 * d$x - 0.2 * d$y
  expect_s4_class(interp2DMBA(d, dx = 1), "SpatRaster")
  expect_true(is.matrix(interp2DMBA(d, dx = 1, output = "matrix", trend = NULL)))
})

test_that("non-finite observations are removed", {
  skip_if_not_installed("MBA")
  d <- data.frame(x = c(0, 1, 2, Inf), y = c(0, 1, 0, 2), z = c(1, 2, 3, 4))
  expect_s4_class(interp2DMBA(d, dx = 1, trend = NULL), "SpatRaster")
  d$x[1:2] <- Inf
  expect_error(interp2DMBA(d, dx = 1, trend = NULL), "finite x, y, and z")
})

test_that("sf input must contain POINT geometries", {
  skip_if_not_installed("MBA")
  p <- sf::st_as_sf(data.frame(z = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1)), coords = c("x", "y"))
  expect_s4_class(interp2DMBA(p, dx = 0.5), "SpatRaster")
  polygon <- sf::st_sf(z = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0))))))
  expect_error(interp2DMBA(polygon, dx = 1), "POINT geometries")
})

test_that("failed linear trends give a clear error", {
  skip_if_not_installed("MBA")
  d <- data.frame(x = 1:5, y = 2 * (1:5), z = 1:5)
  expect_error(interp2DMBA(d, dx = 1, trend = "linear"), "linear trend could not be estimated")
})

test_that("a SpatVector mask is accepted directly", {
  skip_if_not_installed("MBA")
  d <- expand.grid(x = 0:4, y = 0:4)
  d$z <- d$x + d$y
  m <- terra::vect(terra::ext(0, 4, 0, 4), type = "polygons")
  expect_s4_class(interp2DMBA(d, dx = 1, msk = m, trend = NULL), "SpatRaster")
})
