#' Combined multi-species index, per grid cell, per year -- a regional
#' rather than national version of `computeCombinedIndex()`
#'
#' Same geometric-mean logic as the national combined index
#' (`SBI_t = exp(mean(log(index_i,t)))` across species `i`), but computed
#' independently at every cell of a `cellSizeM` grid instead of once for
#' the whole country -- i.e. baseline-normalization AND the cross-species
#' geometric mean both happen cell-by-cell, using `terra` raster algebra
#' (each operation below acts on every cell simultaneously). The result is
#' a map of "how has the multi-species index changed here" instead of one
#' number for all of Germany.
#'
#' A species is excluded from a given cell's geometric mean wherever its
#' own baseline-year value there is too close to zero to divide by
#' (`minBaseline`) -- the per-cell analogue of a species just not being
#' present at that location; it doesn't break the calculation, that
#' species' signal there is simply absent, same principle as
#' `computeCombinedIndex()` dropping a species missing a whole year.
#'
#' @param species Character vector of Latin species names.
#' @param years Integer vector of years.
#' @param baselineYear Integer. Index baseline (must be in `years`).
#' @param metaDir Character. `metaModel()`'s output directory.
#' @param cellSizeM Numeric. Grid cell size in meters (e.g. 20000, 50000).
#' @param minBaseline Numeric. Per-cell baseline values below this are
#'   treated as "species effectively absent here" and excluded from that
#'   cell's geometric mean, rather than producing a huge/unstable ratio.
#' @return SpatRaster, one layer per year (named by year), at `cellSizeM`
#'   resolution.
computeGriddedCombinedIndex <- function(species, years, baselineYear, metaDir,
                                         cellSizeM, minBaseline = 1e-6) {
  baseKey <- as.character(baselineYear)

  speciesGrids <- lapply(species, aggregateSpeciesToGrid, years = years,
                         metaDir = metaDir, cellSizeM = cellSizeM)
  names(speciesGrids) <- species
  speciesGrids <- speciesGrids[!sapply(speciesGrids, is.null)]

  indexGrids <- lapply(names(speciesGrids), function(sp) {
    g <- speciesGrids[[sp]]
    if (!baseKey %in% names(g)) {
      warning(sp, ": missing baseline year ", baselineYear, " -- excluded from gridded index")
      return(NULL)
    }
    base <- terra::ifel(g[[baseKey]] < minBaseline, NA, g[[baseKey]])
    100 * g / base
  })
  names(indexGrids) <- names(speciesGrids)
  indexGrids <- indexGrids[!sapply(indexGrids, is.null)]

  if (length(indexGrids) == 0) stop("No species had a valid baseline year -- cannot compute gridded index.")

  allYearKeys <- sort(unique(unlist(lapply(indexGrids, names))))

  combinedLayers <- lapply(allYearKeys, function(yk) {
    layersForYear <- Filter(Negate(is.null), lapply(indexGrids, function(g) {
      if (yk %in% names(g)) g[[yk]] else NULL
    }))
    stackForYear <- terra::rast(layersForYear)
    exp(mean(log(stackForYear), na.rm = TRUE))
  })

  combined <- terra::rast(combinedLayers)
  names(combined) <- allYearKeys
  combined
}
