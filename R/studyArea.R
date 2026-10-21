#' The study area: ONE definition for the whole pipeline (DECISIONS.md 2026-10-21)
#'
#' Everything the pipeline does is restricted to a study area that is defined ONCE (`sharedStudyArea` in tools/sharedConfig.R) and handed to the modules as a
#' small list (the "spec"). Change the spec and the whole pipeline runs for another area: one German state, Germany, several countries, Europe, or a shapefile.
#' Predictor layers are cropped and masked to the area when they are made (before anything is extracted from them), response records are filtered by it,
#' models are fitted and predicted on what is left. Following the SpaDES convention there are two roles: `studyArea` (where the models are fitted and predicted:
#' the reporting area plus a buffer, per scale) and `studyAreaReporting` (what is reported: the area itself). Only the climate model is fitted on a larger area (Europe,
#' where the EBBA2 atlas is), see `applyStudyAreaToProcessed()`.
#'
#' The spec: `list(label = "Germany", countries = "DEU", regions = NULL, shapefile = NULL)`
#'   countries  ISO3 codes (GADM level 0); several countries are merged.
#'   regions    names of regions (GADM level 1, `NAME_1`, e.g. "Sachsen") of the country; NULL = the whole country.
#'   shapefile  path of a polygon file instead of GADM (lon/lat or EPSG:3035); the other fields are then ignored.
#'   label      name used in messages and stamps.
#'
#' The polygon is built from lon/lat vertices with the closed-form `lonLatToLAEA()` (never through a GDAL/PROJ transformation: the EPSG:3035 axis order of
#' terra/sf transformations has come out swapped on EVE, DECISIONS.md 2026-10-05) and labelled `crs`.
#'
#' @param spec List, see above.
#' @param gadmDir Character. Folder of the geodata cache (inputs/predictors/raw/gadm).
#' @param crs Character. Label of the result (EPSG:3035).
#' @return SpatVector (one dissolved polygon geometry), attribute `label`.
buildStudyArea <- function(spec, gadmDir, crs = "EPSG:3035") {
  if (is.null(spec) || is.null(spec$label)) stop("buildStudyArea(): the study-area spec must be a list with at least a label (tools/sharedConfig.R: sharedStudyArea).")
  toLAEA <- function(g) {                                                         # lon/lat polygons -> EPSG:3035, vertices by the closed form
    v <- terra::geom(g); xy <- lonLatToLAEA(v[, "x"], v[, "y"])
    terra::vect(cbind(v[, c("geom", "part")], x = xy[, "x"], y = xy[, "y"], hole = v[, "hole"]), type = "polygons", crs = crs)
  }
  if (!is.null(spec$shapefile)) {
    g <- terra::vect(spec$shapefile)
    area <- if (isTRUE(terra::is.lonlat(g))) toLAEA(g) else if (grepl("3035", terra::crs(g, describe = TRUE)$code %||% "")) g else
      stop("buildStudyArea(): the shapefile must be in lon/lat or EPSG:3035 (got ", terra::crs(g, describe = TRUE)$code, ").")
  } else {
    if (!requireNamespace("geodata", quietly = TRUE)) stop("buildStudyArea(): the 'geodata' package is needed for the national outline.")
    dir.create(gadmDir, recursive = TRUE, showWarnings = FALSE)
    parts <- list()
    for (cc in spec$countries) {
      g <- if (length(spec$regions)) {
        g1 <- geodata::gadm(country = cc, level = 1, path = gadmDir)
        keep <- tolower(g1$NAME_1) %in% tolower(spec$regions)
        if (!any(keep)) stop("buildStudyArea(): no region of ", cc, " named ", paste(spec$regions, collapse = ", "), " (GADM level 1 names: ", paste(sort(unique(g1$NAME_1)), collapse = ", "), ")")
        g1[keep, ]
      } else geodata::gadm(country = cc, level = 0, path = gadmDir)
      parts[[cc]] <- toLAEA(g)
    }
    if (!length(parts)) stop("buildStudyArea(): the spec has no countries and no shapefile.")
    area <- if (length(parts) == 1) parts[[1]] else do.call(rbind, unname(parts))
  }
  area <- terra::aggregate(area)                                                  # one geometry
  area$label <- spec$label
  e <- as.vector(terra::ext(area))
  if (!(e[1] > 0.9e6 && e[2] < 7.5e6 && e[3] > 0.9e6 && e[4] < 6.1e6))
    stop("buildStudyArea(): the study area is outside Europe's EPSG:3035 box (x ", round(e[1]), "-", round(e[2]), ", y ", round(e[3]), "-", round(e[4]), "): wrong spec or swapped axes.")
  area
}

#' The study area widened by `bufferM` metres (0 = the area itself)
studyAreaWorking <- function(area, bufferM = 0) if (!isTRUE(bufferM > 0)) area else terra::aggregate(terra::buffer(area, width = bufferM))

#' A short text that changes whenever the study area (or its buffer) changes: label, extent, number of vertices, buffer.
#' Used in the stamps of cropped layers and in cache keys, so a layer made for one area is never reused for another.
studyAreaFingerprint <- function(area, bufferM = 0) {
  v <- terra::geom(area)
  paste0(area$label[1], "|", paste(round(as.vector(terra::ext(area))), collapse = ","), "|", nrow(v), "|", round(sum(v[, "x"] + v[, "y"])), "|buffer", bufferM)
}

#' Crop a raster to the study area (extent, whole cells) and set the cells outside it to NA
#' @param area SpatVector from `buildStudyArea()` (use `studyAreaWorking()` for a buffer).
maskToStudyArea <- function(r, area) terra::mask(terra::crop(r, terra::ext(area), snap = "out"), area, touches = TRUE)

#' Which point records (EPSG:3035 x/y) lie inside the study area (plus `bufferM`)?
#'
#' A point-in-polygon test of millions of records against a detailed outline takes tens of minutes (measured: 3 million points, more than 10 min), so the area is
#' rasterised once on a grid of `res` m (cells touched by it count as inside, as in the masks of the layers: `touches = TRUE`) and the records look up their cell
#' (3 million points: 1.5 s). The border is therefore accurate to about one cell. A buffer keeps records on the coast whose centroid lies in the water
#' (2 of 2,637 MhB route centroids lie outside the German outline).
#'
#' @param xy Matrix or data.frame with the columns x and y (EPSG:3035, m).
#' @param area SpatVector from `buildStudyArea()`.
#' @param bufferM Numeric. Widen the area by this many metres.
#' @param res Numeric. Cell size (m) of the lookup grid.
#' @return Logical vector, TRUE = inside.
insideStudyArea <- function(xy, area, bufferM = 0, res = 200) {
  xy <- as.matrix(xy[, c("x", "y")])
  if (!nrow(xy)) return(logical(0))
  wk <- studyAreaWorking(area, bufferM)
  mk <- terra::rasterize(wk, terra::rast(terra::ext(wk), resolution = res, crs = terra::crs(wk)), field = 1, background = NA, touches = TRUE)
  !is.na(terra::extract(mk, xy)[, 1])
}

#' Crop and mask every processed predictor layer to the study area (before anything is extracted from it)
#'
#' Walks `processedRoot/scale_*/*.tif`. Each layer is cropped to the extent of the study area plus the buffer of its scale and the cells outside the (buffered) area are set to
#' NA; the file is replaced (written to `.part` first). A sidecar `<file>.studyarea` holds the fingerprint, so a layer is only redone when the area or the buffer changes.
#' Layers are derived products (the raw downloads are never touched). Buffers: climate scale `climatePredictBufferM` (one coarse cell, so the resampling to the fine grids has
#' neighbours at the border), all other scales `recordsBufferM`. The one climate layer the climate model is TRAINED on (the window ending in `ebba2TrainingYear`) is left alone:
#' the climate model is fitted on all of Europe (where the EBBA2 atlas is), only the prediction windows follow the study area.
#'
#' @param processedRoot Character. `predictors/processed`.
#' @param area SpatVector from `buildStudyArea()`.
#' @param recordsBufferM,climatePredictBufferM Numeric. Buffers (m).
#' @param climateScaleLabel Character. Folder of the climate scale (`scaleLabel(climateResolutionM)`).
#' @param ebba2TrainingYear Integer. Last year of the climate training window.
#' @param dryRun Logical. Only report.
#' @return Invisibly, a data.frame of the layers handled (file, action).
applyStudyAreaToProcessed <- function(processedRoot, area, recordsBufferM = 1000, climatePredictBufferM = 50000, climateScaleLabel = "scale_50",
                                      ebba2TrainingYear = 2017, dryRun = FALSE) {
  dirs <- list.dirs(processedRoot, recursive = FALSE, full.names = TRUE); dirs <- dirs[grepl("^scale_", basename(dirs))]
  rows <- list(); wkCache <- list()
  for (d in dirs) {
    isClimate <- basename(d) == climateScaleLabel
    bufferM <- if (isClimate) climatePredictBufferM else recordsBufferM
    key <- as.character(bufferM); if (is.null(wkCache[[key]])) wkCache[[key]] <- studyAreaWorking(area, bufferM)
    wk <- wkCache[[key]]; stamp <- studyAreaFingerprint(area, bufferM)
    for (f in list.files(d, pattern = "[.]tif$", full.names = TRUE)) {
      sf <- paste0(f, ".studyarea")
      if (isClimate && grepl(paste0("bioclim_[0-9]{4}-", ebba2TrainingYear, "_"), basename(f))) { rows[[length(rows) + 1]] <- data.frame(file = f, action = "training window, left as it is"); next }
      if (file.exists(sf) && identical(readLines(sf, warn = FALSE)[1], stamp)) { rows[[length(rows) + 1]] <- data.frame(file = f, action = "current"); next }
      if (dryRun) { rows[[length(rows) + 1]] <- data.frame(file = f, action = "would crop and mask"); next }
      r <- terra::rast(f); dt <- terra::datatype(r)[1]; nm <- names(r)
      out <- maskToStudyArea(r, wk); names(out) <- nm
      part <- paste0(f, ".part")
      terra::writeRaster(out, part, filetype = "GTiff", overwrite = TRUE, datatype = dt, gdal = c("COMPRESS=DEFLATE"))
      rm(r, out); invisible(gc())
      file.remove(f); file.rename(part, f); writeLines(stamp, sf)
      rows[[length(rows) + 1]] <- data.frame(file = f, action = "cropped and masked")
      message("study area (", area$label[1], ", +", bufferM, " m): ", basename(d), "/", basename(f))
    }
  }
  invisible(do.call(rbind, rows))
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
