# Regional index (10 / 20 / 50 km grids) with uncertainty: per-replicate CELL MEANS of the 200 m predictions.
#
# This module only produces the sufficient statistics: for every species, grid size, year and replicate, the mean probability of
# occurrence of the 200 m pixels inside each coarse cell. The pixels are masked with the German national outline (GADM level 0, the
# same source and the same rule as the baseline regional index, `aggregateSpeciesToGrid()`) BEFORE they are averaged, so a border cell
# holds German pixels only. The index (100 x cell mean / baseline cell mean, geometric mean over species) and the percentile intervals
# are computed by runIndex_Monitor::computeRegionalIndexUncertainty(), as for the national indices. No smoothing here or there:
# smoothing is applied only to the final raw result.
#
# Steps (tools/runUncertaintyTask.R): `regionband` (species x band: partial sums per coarse cell) -> `regionassemble` (species:
# add the bands up -> `regional_means_<km>km.rds`) -> `regionindex` (once: runIndex_Monitor).

#' lon/lat (degrees, GRS80/WGS84) -> ETRS89-LAEA Europe (EPSG:3035) easting/northing (m), closed form (Snyder 1987; EPSG Guidance Note 7-2)
#'
#' A copy of dataPrep_Monitor's `lonLatToLAEA()` so this module does not depend on another one. Pure R, NO GDAL/PROJ: on EVE the
#' EPSG:3035 axis order of terra/sf transformations depends on which geo library initialised first and came out SWAPPED in some
#' sessions (DECISIONS.md 2026-10-05; again 2026-10-08 in a plain terra job, where the German outline got x and y exchanged and
#' `crop()` stopped with "extents do not overlap"). This formula cannot flip.
uncLonLatToLAEA <- function(lon, lat) {
  a <- 6378137; f <- 1 / 298.257222101; e2 <- 2 * f - f^2; e <- sqrt(e2)
  qf <- function(phi) { s <- sin(phi); (1 - e2) * (s / (1 - e2 * s^2) - (1 / (2 * e)) * log((1 - e * s) / (1 + e * s))) }
  phi <- lat * pi / 180; lam <- lon * pi / 180; phi0 <- 52 * pi / 180; lam0 <- 10 * pi / 180
  q <- qf(phi); qp <- qf(pi / 2); q0 <- qf(phi0)
  Rq <- a * sqrt(qp / 2); beta <- asin(q / qp); beta0 <- asin(q0 / qp)
  m1 <- cos(phi0) / sqrt(1 - e2 * sin(phi0)^2); D <- a * m1 / (Rq * cos(beta0))
  B <- Rq * sqrt(2 / (1 + sin(beta0) * sin(beta) + cos(beta0) * cos(beta) * cos(lam - lam0)))
  cbind(x = 4321000 + B * D * cos(beta) * sin(lam - lam0),
        y = 3210000 + (B / D) * (cos(beta0) * sin(beta) - sin(beta0) * cos(beta) * cos(lam - lam0)))
}

#' German national outline (GADM level 0, as in the baseline) in EPSG:3035, built from the lon/lat vertices with `uncLonLatToLAEA()`
#' (never through a GDAL/PROJ transformation, see above) and labelled with the CRS of the predictions. Stops if the result is not in Germany's box.
#' @param gadmDir folder of the geodata cache (inputs/predictors/raw/gadm); @param crs CRS of the rasters it will be used with.
uncOutlineLAEA <- function(gadmDir, crs) {
  if (!requireNamespace("geodata", quietly = TRUE))
    stop("The regional step needs the 'geodata' package (national outline); the baseline regional index needs it too.")
  dir.create(gadmDir, recursive = TRUE, showWarnings = FALSE)
  g <- geodata::gadm(country = "DEU", level = 0, path = gadmDir)                 # lon/lat
  v <- terra::geom(g); xy <- uncLonLatToLAEA(v[, "x"], v[, "y"])
  out <- terra::vect(cbind(v[, c("geom", "part")], x = xy[, "x"], y = xy[, "y"], hole = v[, "hole"]), type = "polygons", crs = crs)
  e <- as.vector(terra::ext(out))
  if (!(e[1] > 3.9e6 && e[2] < 4.8e6 && e[3] > 2.6e6 && e[4] < 3.7e6))
    stop("German outline outside Germany's EPSG:3035 box: x ", round(e[1]), "-", round(e[2]), ", y ", round(e[3]), "-", round(e[4]))
  out
}

#' German outline for the regional step, in the CRS `crs` of the predictions (inputs/predictors/raw/gadm cache, as in the baseline)
uncGermanyBoundary <- function(cfg, crs) uncOutlineLAEA(file.path(cfg$inputRoot, "predictors", "raw", "gadm"), crs)

#' Coarse grid of one cell size: the 200 m window cropped to the outline, then aggregated by the same factor as the baseline
#' (`aggregateSpeciesToGrid()`: crop to the boundary, `aggregate(fact = round(cellSizeM / res))`), so the cells are the same.
uncRegionalGrid <- function(win, bProj, cellSizeM) {
  fact <- max(1L, as.integer(round(cellSizeM / terra::res(win)[1])))
  list(template = terra::aggregate(terra::crop(win, bProj), fact = fact), fact = fact)
}

#' Partial sums for one species and one band: for every year, per coarse cell and replicate, the sum of the predictions of the German
#' pixels of this band and their number. Written per run label under `regional/<label>/band<kk>.rds`; finished files are skipped.
#'
#' @param win window template (`uncBands(cfg, sp)$window`); @param bProj outline in the CRS of `win`.
uncRegionalBand <- function(cfg, sp, band, win, bProj, cellSizesM) {
  labs <- uncRunLabels(cfg, sp)
  if (!length(labs)) { message(sp, ": no replicate predictions found -- nothing to do"); return(invisible(NULL)) }
  outOf <- function(lab) file.path(uncSpDir(cfg, sp, "regional", lab), sprintf("band%02d.rds", band$k))
  todo <- labs[!vapply(labs, function(l) file.exists(outOf(l)), logical(1))]
  if (!length(todo)) return(invisible("cached"))

  tmpl <- band$template
  ger <- terra::rasterize(bProj, tmpl, field = 1, touches = TRUE)         # the same rule as terra::mask(x, outline) (touches = TRUE)
  gerIdx <- which(!is.na(terra::values(ger)[, 1])); rm(ger)
  lookup <- list(); nGer <- list()
  for (cs in cellSizesM) {
    g <- uncRegionalGrid(win, bProj, cs)
    cg <- rep(NA_integer_, terra::ncell(tmpl))
    if (length(gerIdx)) cg[gerIdx] <- as.integer(terra::cellFromXY(g$template, terra::xyFromCell(tmpl, gerIdx)))
    lookup[[as.character(cs)]] <- cg
    nGer[[as.character(cs)]] <- tabulate(cg[!is.na(cg)], nbins = terra::ncell(g$template))
  }
  for (lab in todo) {
    perSize <- lapply(cellSizesM, function(cs) list(nGerman = nGer[[as.character(cs)]], perYear = list()))
    names(perSize) <- as.character(cellSizesM)
    ids <- NULL; yearsDone <- integer(0)
    for (yr in cfg$outYears) {
      f <- file.path(uncSpDir(cfg, sp, "pred", lab), sprintf("%d_band%02d.rds", yr, band$k))
      if (!file.exists(f)) next
      p <- uncReadInt16(f); ids <- p$ids; yearsDone <- c(yearsDone, yr)
      for (cs in cellSizesM) {
        cg <- lookup[[as.character(cs)]][p$idx]; keep <- !is.na(cg)
        if (!any(keep)) next
        P <- p$P[keep, , drop = FALSE]; ok <- !is.na(P); P[!ok] <- 0
        S <- rowsum(P, cg[keep], reorder = TRUE); N <- rowsum(ok * 1, cg[keep], reorder = TRUE)
        perSize[[as.character(cs)]]$perYear[[as.character(yr)]] <- list(cells = as.integer(rownames(S)), S = S, N = N)
      }
      rm(p)
    }
    uncSaveRDS(list(species = sp, band = band$k, label = lab, ids = ids, years = yearsDone, bySize = perSize), outOf(lab))
    message(sp, " band ", band$k, " (", lab, "): regional partial sums for ", length(yearsDone), " years, ",
            length(cellSizesM), " grids, ", length(gerIdx), " German pixels")
  }
  invisible("done")
}

#' Add the band pieces up per species: `regional/regional_means_<km>km.rds` with the mean probability per coarse cell, year and replicate
#' (array cells x years x replicates, NA where a cell has no data), the German pixel count per cell and the grid geometry.
uncRegionalAssemble <- function(cfg, sp, win, bProj, cellSizesM) {
  labs <- uncRunLabels(cfg, sp); years <- cfg$outYears
  if (!length(labs)) { message(sp, ": no replicate predictions found"); return(invisible(NULL)) }
  natRows <- list()                                   # German national area means (all German pixels), from the first grid's sums
  for (cs in cellSizesM) {
    km <- cs / 1000; g <- uncRegionalGrid(win, bProj, cs); nc <- terra::ncell(g$template)
    pieces <- list(); idsAll <- integer(0); nGerAll <- NULL
    for (lab in labs) {
      fs <- file.path(uncSpDir(cfg, sp, "regional", lab), sprintf("band%02d.rds", seq_len(cfg$nBands)))
      if (!all(file.exists(fs))) { warning(sp, " ", lab, ": ", sum(!file.exists(fs)), " regional band pieces missing -- run skipped", call. = FALSE); next }
      parts <- lapply(fs, readRDS)
      ids <- parts[[1]]$ids
      S <- array(0, c(nc, length(years), length(ids))); N <- S; nGer <- numeric(nc)
      for (pt in parts) {
        b <- pt$bySize[[as.character(cs)]]; nGer <- nGer + b$nGerman
        for (yi in seq_along(years)) {
          e <- b$perYear[[as.character(years[yi])]]
          if (is.null(e)) next
          S[e$cells, yi, ] <- S[e$cells, yi, ] + e$S; N[e$cells, yi, ] <- N[e$cells, yi, ] + e$N
        }
      }
      if (cs == cellSizesM[1]) {                                  # national mean over the German pixels = sum of the cell sums / number of pixels
        Sy <- apply(S, c(2, 3), sum); Ny <- apply(N, c(2, 3), sum)
        nr <- data.frame(species = sp, replicate = rep(ids, each = length(years)), year = rep(years, times = length(ids)),
                         areaMean = as.vector(Sy / Ny), nCells = as.vector(Ny))
        natRows[[lab]] <- nr[nr$nCells > 0, ]
      }
      m <- S / N; m[N == 0] <- NA_real_
      pieces[[lab]] <- m; idsAll <- c(idsAll, ids); if (is.null(nGerAll)) nGerAll <- nGer
    }
    if (!length(pieces)) next
    mean3 <- array(unlist(pieces, use.names = FALSE), c(nc, length(years), length(idsAll)))
    keep <- which(nGerAll > 0)
    out <- list(species = sp, cellSizeM = cs, fact = g$fact, years = years, ids = idsAll, cells = keep,
                nGerman = nGerAll[keep], mean = mean3[keep, , , drop = FALSE],
                grid = list(nrow = terra::nrow(g$template), ncol = terra::ncol(g$template), ext = as.vector(terra::ext(g$template)),
                            crs = terra::crs(g$template)))
    f <- file.path(uncSpDir(cfg, sp, "regional"), sprintf("regional_means_%dkm.rds", km))
    uncSaveRDS(out, f)
    message(sp, ": ", km, " km grid: ", length(keep), " cells x ", length(years), " years x ", length(idsAll), " replicates -> ", basename(f))
  }
  if (length(natRows)) {
    nat <- do.call(rbind, natRows)
    utils::write.csv(nat, file.path(uncSpDir(cfg, sp), "area_mean_replicates_germany.csv"), row.names = FALSE)
    message(sp, ": German-only area means of ", length(unique(nat$replicate)), " replicates x ", length(unique(nat$year)), " years -> area_mean_replicates_germany.csv")
  }
  invisible(TRUE)
}

#' Germany-only COPIES of the finished uncertainty maps of one species (`maps/` -> `maps_germany/`) or of the community maps
#' (`community/` -> `community_germany/`, `sp = NULL`): same files, pixels outside the German outline set to NA. The originals stay.
#' TEMPORARY (DECISIONS.md 2026-10-08): the proper fix is to cut all inputs to the study area at the source.
uncMaskMaps <- function(cfg, sp = NULL, outline = NULL) {
  src <- if (is.null(sp)) file.path(uncRoot(cfg), "community") else uncSpDir(cfg, sp, "maps")
  dst <- if (is.null(sp)) file.path(uncRoot(cfg), "community_germany") else file.path(uncSpDir(cfg, sp), "maps_germany")
  fs <- list.files(src, pattern = "[.]tif$", full.names = TRUE)
  if (!length(fs)) { message("No maps to mask in ", src); return(invisible(NULL)) }
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  r1 <- terra::rast(fs[1])
  if (is.null(outline)) outline <- uncGermanyBoundary(cfg, terra::crs(r1))
  mk <- terra::rasterize(outline, r1[[1]], field = 1, touches = TRUE)
  for (f in fs) {
    out <- file.path(dst, basename(f))
    if (file.exists(out)) next
    r <- terra::rast(f); m <- if (isTRUE(terra::compareGeom(r, mk, stopOnError = FALSE))) mk else terra::rasterize(outline, r[[1]], field = 1, touches = TRUE)
    part <- sub("[.]tif$", ".part.tif", out)
    terra::writeRaster(terra::mask(r, m), part, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
    file.rename(part, out)
  }
  message("Germany-only maps: ", length(fs), " files -> ", dst)
  invisible(dst)
}
