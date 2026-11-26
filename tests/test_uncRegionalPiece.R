# Test of uncRegionalPieceComplete() (2026-11-26; from the repo root: Rscript modules/models_Monitor/tests/test_uncRegionalPiece.R)
src <- readLines("modules/models_Monitor/R/uncRegional.R"); i <- grep("^uncRegionalPieceComplete <- function", src); j <- i + which(src[(i + 1):length(src)] == "}")[1]
eval(parse(text = src[i:j]))
chk <- function(label, ok) { cat(if (isTRUE(ok)) "ok   " else "FAIL ", label, "\n"); if (!isTRUE(ok)) quit(status = 1) }
d <- tempfile(); dir.create(d); f <- file.path(d, "band01.rds"); saveRDS(list(years = 2017:2025), f)
chk("a piece with years 2017-2025 is complete for 2017-2025 (cached as before)", uncRegionalPieceComplete(f, 2017:2025))
chk("the same piece is NOT complete for 2017-2026 (the new year makes it stale)", !uncRegionalPieceComplete(f, 2017:2026))
chk("a missing file is not complete", !uncRegionalPieceComplete(file.path(d, "none.rds"), 2017:2025))
writeLines("not an rds", g <- file.path(d, "bad.rds")); chk("an unreadable file is not complete (no error)", !uncRegionalPieceComplete(g, 2017:2025))
saveRDS(list(years = integer(0)), h <- file.path(d, "empty.rds")); chk("a piece without years is not complete", !uncRegionalPieceComplete(h, 2017:2025))
saveRDS(list(years = 2017:2025, bySize = list(`10000` = 1, `20000` = 1)), k <- file.path(d, "sized.rds"))
chk("sizes present: complete", uncRegionalPieceComplete(k, 2017:2025, c(10000, 20000)))
chk("a requested size the piece lacks (1000 m): NOT complete", !uncRegionalPieceComplete(k, 2017:2025, c(1000, 10000)))
cat("all passed\n")
