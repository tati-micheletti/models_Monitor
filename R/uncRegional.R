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

#' German national outline (GADM level 0), in the coordinate system `crs`. Cached under inputs/predictors/raw/gadm (as in the baseline).
uncGermanyBoundary <- function(cfg, crs) {
  if (!requireNamespace("geodata", quietly = TRUE))
    stop("The regional step needs the 'geodata' package (national outline); the baseline regional index needs it too.")
  d <- file.path(cfg$inputRoot, "predictors", "raw", "gadm"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  terra::project(geodata::gadm(country = "DEU", level = 0, path = d), crs)
}

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
  invisible(TRUE)
}
