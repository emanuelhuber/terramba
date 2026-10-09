#' 2D raster interpolation using multilevel B-splines
#'
#' Interpolates spatial point data to a regular grid using
#' [MBA::mba.surf()]. The result is returned as a [terra::SpatRaster] or a
#' matrix.
#'
#' @param x `(matrix[n,3]|sf)` A matrix or data frame whose first three
#' columns are x, y, and z, or an `sf` object with POINT geometries.
#' For an `sf` object, select the attribute to interpolate with `z`.
#' If `z = NULL`, the object must contain exactly one non-geometry column,
#' which is then used automatically.
#' @param dx `(numeric[1])` Grid resolution in the x-direction.
#' @param dy `(numeric[1])` Grid resolution in the y-direction. If missing,
#' `dx` is used.
#' @param xlim `(numeric[2]|NULL)` Optional x-limits. If `NULL`, they are
#'   determined from `extend`.
#' @param ylim `(numeric[2]|NULL)` Optional y-limits. If `NULL`, they are
#'   determined from `extend`.
#' @param buffer `(numeric[1])` Deprecated-compatible alias for `bufferDist`.
#'   Retained for backward compatibility.
#' @param extend Extent method: `"bbox"` (axis-aligned bounding box),
#'   `"obbox"` (oriented bounding box), `"chull"` (convex hull), or
#'   `"buffer"` (buffered convex hull). Polygon methods clip the output to
#'   their polygon.
#' @param bufferDist Non-negative buffer distance in coordinate units.
#'   For `"bbox"` it expands the bounds; for polygon methods it buffers the
#'   polygon. `"buffer"` requires a positive value. If `NULL`, `buffer` is used.
#' @param h `(numeric[1])` Exponent of the power term used by [MBA::mba.surf()].
#' @param m `(numeric[1])` Number of coefficients used in the B-spline
#' approximation. If both are `NULL`, they are selected from the aspect ratio.
#' @param n `(numeric[1])` initial size of the spline space in the hierarchical
#' construction along the x axis. If the rectangular domain is a square,
#' `n = m = 1` is recommended. If the x axis is `k` times the length of the
#' y axis, `n = 1, m = k` is recommended. The default is `n = 1`.
#' @param output Either `"raster"` or `"matrix"`.
#' @param msk Optional [terra::SpatVector] or `sf` object used to define the
#'   output extent and mask the result.
#' @param trend Either `"linear"` to remove and restore a linear trend, or
#'   `NULL`/`"none"` for no trend removal.
#' @param z Optional name or column index of the value attribute in an `sf`
#'   object.
#' @param crs Optional fallback coordinate reference system for raster output.
#'   Any value accepted by [terra::crs()] can be used, for example
#'   `"EPSG:2056"` or WKT. The argument is used only when no CRS can be
#'   inferred from `x` or `msk`. If a CRS can be inferred and `crs` is also
#'   supplied, `crs` is ignored with a warning. If both `x` and `msk` provide
#'   CRS metadata, their coordinate reference systems must match. The argument
#'   assigns metadata only and does not transform coordinates. It has no effect
#'   on `output = "matrix"`.
#'
#' @details
#' Rows with non-finite x, y, or z values are removed before interpolation.
#' At least three finite observations must remain.
#'
#' With `trend = "linear"`, a plane of the form `z ~ x + y` is fitted. MBA is
#' applied to the residuals and the fitted plane is subsequently added back.
#' An informative error is raised when the plane cannot be estimated.
#'
#' For compatibility with the original implementation, `m` and `n` are
#' deliberately mapped to `n = m` and `m = n` in [MBA::mba.surf()]. The raster
#' uses the package's established node-based grid convention.
#'
#' @return A [terra::SpatRaster] when `output = "raster"`, otherwise a list
#'   containing `x`, `y`, and `z` (the interpolated matrix).
#' @export
#'
#' @examples
#' if (requireNamespace("MBA", quietly = TRUE)) {
#'   d <- expand.grid(x = seq(0, 10, by = 2), y = seq(0, 8, by = 2))
#'   d$z <- 2 + 0.5 * d$x - 0.2 * d$y
#'   r <- interp2DMBA(d, dx = 1)
#' }
interp2DMBA <- function(x, dx, dy = NULL, xlim = NULL, ylim = NULL,
                        buffer = 0, extend = c("bbox", "obbox", "chull", "buffer"),
                        bufferDist = NULL, h = 6, m = NULL, n = NULL,
                        output = c("raster", "matrix"), msk = NULL,
                        trend = "linear", z = NULL,
                        crs = NULL) {
  if (missing(dx) && is.null(dy)) {
    stop("You must specify at least 'dx' or 'dy'.", call. = FALSE)
  }
  if (missing(dx)) dx <- dy
  if (is.null(dy)) dy <- dx

  check_resolution <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop(sprintf("'%s' must be a positive finite number.", name),
           call. = FALSE)
    }
  }
  check_resolution(dx, "dx")
  check_resolution(dy, "dy")

  if (xor(is.null(m), is.null(n))) {
    stop("Supply both 'm' and 'n', or leave both as NULL.", call. = FALSE)
  }
  output <- match.arg(output)
  extend <- match.arg(extend)
  if (identical(trend, "none")) trend <- NULL
  if (!is.null(trend) && !identical(trend, "linear")) {
    stop("'trend' must be 'linear', 'none', or NULL.", call. = FALSE)
  }

  # Keep the legacy `buffer` argument working. `bufferDist`, when supplied,
  # takes precedence. For bbox it expands the rectangle; for polygon methods
  # it buffers the selected polygon.
  if (is.null(bufferDist)) bufferDist <- buffer
  if (!is.numeric(bufferDist) || length(bufferDist) != 1L ||
      !is.finite(bufferDist) || bufferDist < 0) {
    stop("'bufferDist' (and legacy 'buffer') must be one non-negative finite number.",
         call. = FALSE)
  }
  if (extend == "buffer" && bufferDist <= 0) {
    stop("When 'extend = \"buffer\"', 'bufferDist' must be greater than zero.",
         call. = FALSE)
  }

  check_limits <- function(lim, name) {
    if (!is.null(lim) &&
        (!is.numeric(lim) || length(lim) != 2L || any(!is.finite(lim)) ||
         lim[1] == lim[2])) {
      stop(sprintf("'%s' must be NULL or two distinct finite numeric values.", name),
           call. = FALSE)
    }
    if (!is.null(lim)) sort(as.numeric(lim)) else NULL
  }
  xlim <- check_limits(xlim, "xlim")
  ylim <- check_limits(ylim, "ylim")

  has_crs <- function(value) {
    is.character(value) && length(value) == 1L &&
      !is.na(value) && nzchar(value)
  }

  is_sf <- inherits(x, "sf")
  is_spatvector <- inherits(msk, "SpatVector")

  if (is_sf) {
    geom_types <- unique(as.character(sf::st_geometry_type(x)))
    if (!length(geom_types) || !all(geom_types == "POINT")) {
      stop("'x' must contain only POINT geometries.", call. = FALSE)
    }
    geom_col <- attr(x, "sf_column")
    value_cols <- setdiff(names(x), geom_col)
    if (is.null(z)) {
      if (length(value_cols) != 1L) {
        stop("For an 'sf' object with multiple attributes, specify 'z'.",
             call. = FALSE)
      }
      z_col <- value_cols[[1L]]
    } else {
      z_col <- if (is.numeric(z)) names(x)[z] else z
      if (length(z_col) != 1L || is.na(z_col) || !z_col %in% value_cols) {
        stop("'z' must identify one non-geometry column in 'x'.", call. = FALSE)
      }
    }
    coords <- sf::st_coordinates(x)
    vals <- x[[z_col]]
    x_crs_obj <- sf::st_crs(x)
    x_crs <- if (!is.na(x_crs_obj)) x_crs_obj$wkt else NULL
  } else {
    x <- as.matrix(x)
    if (ncol(x) < 3L) {
      stop("'x' must contain at least three columns: x, y, and z.",
           call. = FALSE)
    }
    suppressWarnings(storage.mode(x) <- "double")
    coords <- x[, 1:2, drop = FALSE]
    vals <- x[, 3]
    x_crs <- NULL
  }

  keep <- is.finite(coords[, 1]) & is.finite(coords[, 2]) & is.finite(vals)
  coords <- coords[keep, 1:2, drop = FALSE]
  vals <- as.numeric(vals[keep])
  if (nrow(coords) < 3L) {
    stop("At least three points with finite x, y, and z values are required.",
         call. = FALSE)
  }
  xy <- cbind(coords, vals)

  msk_crs <- NULL
  if (!is.null(msk)) {
    if (is_spatvector) {
      mask_bbox <- as.vector(terra::ext(msk))
      msk_crs_value <- terra::crs(msk)
      if (has_crs(msk_crs_value)) msk_crs <- msk_crs_value
    } else if (inherits(msk, c("sf", "sfc"))) {
      bb <- sf::st_bbox(msk)
      mask_bbox <- unname(bb[c("xmin", "xmax", "ymin", "ymax")])
      msk_crs_obj <- sf::st_crs(msk)
      if (!is.na(msk_crs_obj)) msk_crs <- msk_crs_obj$wkt
    } else {
      stop("'msk' must be a 'SpatVector', 'sf', or 'sfc' object.",
           call. = FALSE)
    }
  } else {
    mask_bbox <- NULL
  }

  # Spatial inputs are authoritative. A supplied `crs` is only a fallback
  # when no CRS can be inferred from either `x` or `msk`.
  if (!is.null(x_crs) && !is.null(msk_crs) &&
      !terra::same.crs(x_crs, msk_crs)) {
    stop(
      "The coordinate reference systems of 'x' and 'msk' do not match. ",
      "Transform one of them before calling interp2DMBA().",
      call. = FALSE
    )
  }

  crs_source <- if (!is.null(x_crs)) {
    "'x'"
  } else if (!is.null(msk_crs)) {
    "'msk'"
  } else {
    NULL
  }
  r_crs <- if (!is.null(x_crs)) x_crs else msk_crs

  if (output == "raster") {
    if (!is.null(r_crs) && !is.null(crs)) {
      warning(
        sprintf(
          "'crs' was ignored because a coordinate reference system was inferred from %s.",
          crs_source
        ),
        call. = FALSE
      )
    } else if (is.null(r_crs) && !is.null(crs)) {
      crs_test <- terra::rast(nrows = 1, ncols = 1)
      r_crs <- tryCatch({
        terra::crs(crs_test) <- crs
        validated_crs <- terra::crs(crs_test)
        if (!has_crs(validated_crs)) stop("Invalid CRS.")
        validated_crs
      }, error = function(e) {
        stop(
          "'crs' is not a valid coordinate reference system for terra.",
          call. = FALSE
        )
      })
    }
  }

  # Construct a polygon for methods that need a non-rectangular/oriented
  # domain. These methods need sf; the bbox method does not.
  clip_polygon <- NULL
  if (extend != "bbox") {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required for extend = 'obbox', 'chull', or 'buffer'.",
           call. = FALSE)
    }
    point_sfc <- sf::st_as_sf(
      data.frame(.x = coords[, 1], .y = coords[, 2]),
      coords = c(".x", ".y"),
      crs = if (is_sf) sf::st_crs(x) else sf::st_crs(NA)
    )
    point_union <- sf::st_union(sf::st_geometry(point_sfc))
    clip_polygon <- if (extend == "obbox") {
      sf::st_minimum_rotated_rectangle(point_union)
    } else {
      sf::st_convex_hull(point_union)
    }
    if (is.null(clip_polygon) || all(sf::st_is_empty(clip_polygon))) {
      stop("Could not construct a valid interpolation polygon from the input points.",
           call. = FALSE)
    }
    if (extend == "buffer" || bufferDist > 0) {
      clip_polygon <- sf::st_buffer(clip_polygon, dist = bufferDist)
    }
    if (any(!sf::st_is_valid(clip_polygon))) {
      clip_polygon <- sf::st_make_valid(clip_polygon)
    }
    if (all(as.numeric(sf::st_area(clip_polygon)) <= 0)) {
      stop("The selected extent has zero area; check the observation coordinates.",
           call. = FALSE)
    }
    bb <- sf::st_bbox(clip_polygon)
    auto_xlim <- unname(bb[c("xmin", "xmax")])
    auto_ylim <- unname(bb[c("ymin", "ymax")])
  } else {
    if (!is.null(msk)) {
      auto_xlim <- mask_bbox[1:2]
      auto_ylim <- mask_bbox[3:4]
    } else {
      auto_xlim <- range(coords[, 1])
      auto_ylim <- range(coords[, 2])
    }
    if (bufferDist > 0) {
      auto_xlim <- auto_xlim + c(-bufferDist, bufferDist)
      auto_ylim <- auto_ylim + c(-bufferDist, bufferDist)
    }
  }

  # Explicit xlim/ylim independently override automatically derived limits.
  xlim_use <- if (is.null(xlim)) auto_xlim else xlim
  ylim_use <- if (is.null(ylim)) auto_ylim else ylim
  if (diff(xlim_use) <= 0 || diff(ylim_use) <= 0) {
    stop("The interpolation extent must have positive width and height.",
         call. = FALSE)
  }

  # Snap/expand the requested bounds to an integer number of grid intervals.
  para <- .getbbox_nx_ny(
    xpos = xlim_use, ypos = ylim_use, dx = dx, dy = dy, buffer = 0
  )

  if (is.null(m) && is.null(n)) {
    m <- n <- 1L
    ratio_x_y <- diff(para$bbox[3:4]) / diff(para$bbox[1:2])
    if (!is.finite(ratio_x_y) || ratio_x_y <= 0) {
      stop("Could not determine MBA refinement from the interpolation extent.",
           call. = FALSE)
    }
    if (ratio_x_y < 1) m <- max(1L, round(1 / ratio_x_y)) else
      n <- max(1L, round(ratio_x_y))
  }

  # Retain the original interpolation coordinate convention and MBA argument
  # mapping, which are relied on by existing regression tests.
  if (!is.null(trend)) {
    x_values <- seq(para$bbox[1], para$bbox[2] + dx, length.out = para$nx)
    y_values <- seq(para$bbox[3], para$bbox[4] + dx, length.out = para$ny)
    grid <- expand.grid(x = x_values, y = y_values)
    df <- data.frame(x = xy[, 1], y = xy[, 2], z = xy[, 3])
    fit <- tryCatch(stats::lm(z ~ x + y, data = df), error = function(e) NULL)
    if (is.null(fit) || fit$rank < 3L ||
        any(!is.finite(stats::coef(fit)))) {
      stop("The linear trend could not be estimated from the supplied points.",
           call. = FALSE)
    }
    predicted_z <- tryCatch(
      stats::predict(fit, newdata = grid),
      error = function(e) rep(NA_real_, nrow(grid))
    )
    if (length(predicted_z) != nrow(grid) || any(!is.finite(predicted_z))) {
      stop("The fitted linear trend could not be predicted on the output grid.",
           call. = FALSE)
    }
    r0 <- terra::rast(
      ncols = para$nx, nrows = para$ny,
      xmin = para$bbox[1], xmax = para$bbox[2] + dx,
      ymin = para$bbox[3], ymax = para$bbox[4] + dy
    )
    if (!is.null(r_crs)) terra::crs(r0) <- r_crs
    z_matrix <- matrix(predicted_z, nrow = para$ny, ncol = para$nx,
                       byrow = TRUE)
    terra::values(r0) <- z_matrix[nrow(z_matrix):1, ]
    xy[, 3] <- stats::residuals(fit)
  }

  S <- suppressWarnings(MBA::mba.surf(
    xy, para$nx, para$ny, n = m, m = n,
    extend = TRUE, h = h, b.box = para$bbox
  )$xyz.est)
  SZ <- matrix(rev(as.vector(S$z)), nrow = para$nx, ncol = para$ny)
  r <- terra::rast(
    ncols = para$nx, nrows = para$ny,
    xmin = para$bbox[1], xmax = para$bbox[2] + dx,
    ymin = para$bbox[3], ymax = para$bbox[4] + dy
  )
  if (!is.null(r_crs)) terra::crs(r) <- r_crs
  terra::values(r) <- t(SZ[nrow(SZ):1, ])
  if (!is.null(trend)) r <- r + r0

  # Apply geometric clipping for obbox/chull/buffer and then the optional
  # user-provided mask. Masking never changes the selected raster extent.
  if (!is.null(clip_polygon)) {
    clip_vector <- terra::vect(clip_polygon)
    if (!is.null(r_crs) && !has_crs(terra::crs(clip_vector))) {
      terra::crs(clip_vector) <- r_crs
    }
    r <- terra::mask(r, clip_vector)
  }
  if (!is.null(msk)) {
    mask_vector <- if (is_spatvector) msk else terra::vect(msk)
    r <- terra::mask(r, mask_vector)
  }

  if (output == "matrix") {
    return(list(
      x = terra::xFromCol(r, seq_len(terra::ncol(r))),
      y = rev(terra::yFromRow(r, seq_len(terra::nrow(r)))),
      z = terra::as.matrix(r, wide = TRUE)[
        terra::nrow(r):1, , drop = FALSE
      ]
    ))
  }
  r
}

.getbbox_nx_ny <- function(xpos, ypos, dx, dy, buffer = NULL) {
  xpos_rg <- range(xpos)
  ypos_rg <- range(ypos)
  if (is.null(buffer)) buffer <- min(diff(xpos_rg) * 0.05, diff(ypos_rg) * 0.05)
  if (!is.null(buffer) && min(buffer) > 0) {
    if (length(buffer) == 1L) buffer <- buffer[c(1, 1)]
    xpos_rg <- xpos_rg + c(-1, 1) * buffer
    ypos_rg <- ypos_rg + c(-1, 1) * buffer
  }
  bbox <- c(xpos_rg, ypos_rg)
  bbox_dx <- bbox[2] - bbox[1]
  bbox_dy <- bbox[4] - bbox[3]
  nx <- ceiling(bbox_dx / dx)
  ny <- ceiling(bbox_dy / dy)
  Dx <- (nx * dx - bbox_dx) / 2
  Dy <- (ny * dy - bbox_dy) / 2
  bbox[1:2] <- bbox[1:2] + c(-1, 1) * Dx
  bbox[3:4] <- bbox[3:4] + c(-1, 1) * Dy
  list(bbox = bbox, nx = nx + 1, ny = ny + 1)
}
