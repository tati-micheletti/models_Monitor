# Test: a grid size added later is ADDED to a finished regional band piece, the sizes already in it stay identical (2026-11-26).
# From the repo root: Rscript modules/models_Monitor/tests/test_uncRegionalBandSizes.R
suppressMessages(library(terra))
src <- readLines("modules/models_Monitor/R/uncRegional.R")
grab <- function(name) { i <- grep(paste0("^", name, " <- function"), src); j <- i + which(src[(i + 1):length(src)] == "}")[1]; eval(parse(text = src[i:j]), envir = globalenv()) }
for (f in c("uncRegionalGrid", "uncRegionalPieceComplete", "uncRegionalBand")) grab(f)
chk <- function(label, ok) { cat(if (isTRUE(ok)) "ok   " else "FAIL ", label, "\n"); if (!isTRUE(ok)) quit(status = 1) }

tmp <- tempfile(); dir.create(tmp)
cfg <- list(outYears = 2017:2018, nBands = 1)
uncRunLabels <- function(cfg, sp) "L1"
uncSpDir <- function(cfg, sp, kind, lab = NULL) { p <- file.path(tmp, sp, kind, lab); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
uncSaveRDS <- function(x, f) { dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); saveRDS(x, f) }
win <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 4000, ymin = 0, ymax = 4000, crs = "EPSG:3035")   # 200 m pixels
bProj <- vect(cbind(c(0, 4000, 4000, 0), c(0, 0, 4000, 4000)), type = "polygons", crs = "EPSG:3035")
band <- list(k = 1L, template = win)
set.seed(1); nc <- ncell(win)
for (yr in cfg$outYears) { dir.create(uncSpDir(cfg, "sp", "pred", "L1"), recursive = TRUE, showWarnings = FALSE); saveRDS(list(ids = 1:3, idx = seq_len(nc), P = matrix(runif(nc * 3), nc, 3)), file.path(uncSpDir(cfg, "sp", "pred", "L1"), sprintf("%d_band01.rds", yr))) }
uncReadInt16 <- function(f) readRDS(f)
piece <- file.path(uncSpDir(cfg, "sp", "regional", "L1"), "band01.rds")

r1 <- uncRegionalBand(cfg, "sp", band, win, bProj, c(1000)); a <- readRDS(piece)
chk("first run computes the 1000 m grid", identical(r1, "done") && identical(names(a$bySize), "1000") && identical(a$years, 2017:2018))
chk("same sizes again: cached", identical(uncRegionalBand(cfg, "sp", band, win, bProj, c(1000)), "cached"))
chk("a size the piece lacks makes it NOT complete", !uncRegionalPieceComplete(piece, cfg$outYears, c(1000, 2000)) && uncRegionalPieceComplete(piece, cfg$outYears, 1000))
r2 <- uncRegionalBand(cfg, "sp", band, win, bProj, c(1000, 2000)); b <- readRDS(piece)
chk("asking for 1000 and 2000 computes and ADDS the 2000 m grid", identical(r2, "done") && setequal(names(b$bySize), c("1000", "2000")))
chk("the 1000 m grid already in the piece is identical to before", identical(b$bySize[["1000"]], a$bySize[["1000"]]))
chk("the 2000 m sums add up to the 1000 m sums (same pixels)", isTRUE(all.equal(sum(unlist(lapply(b$bySize[["2000"]]$perYear[["2017"]]$S, sum))), sum(unlist(lapply(b$bySize[["1000"]]$perYear[["2017"]]$S, sum))))))
chk("a run asking for only 2000 (a subset) is cached and leaves the piece alone", identical(uncRegionalBand(cfg, "sp", band, win, bProj, c(2000)), "cached"))
cfg$outYears <- 2017:2019; saveRDS(readRDS(file.path(uncSpDir(cfg, "sp", "pred", "L1"), "2018_band01.rds")), file.path(uncSpDir(cfg, "sp", "pred", "L1"), "2019_band01.rds"))
r3 <- uncRegionalBand(cfg, "sp", band, win, bProj, c(1000)); d <- readRDS(piece)
chk("a new year recomputes the piece AND keeps the 2000 m grid it already had", identical(r3, "done") && setequal(names(d$bySize), c("1000", "2000")) && identical(d$years, 2017:2019))
chk("the sums of the old years are unchanged by the recompute", identical(d$bySize[["1000"]]$perYear[["2017"]], a$bySize[["1000"]]$perYear[["2017"]]))
cat("all passed\n")
