# Local test of the ensemble workflow on the mini baseline (Alauda arvensis, outputs/utest): model jobs -> ensembles of chosen member sets.
#   Rscript modules/models_Monitor/tests/ensemble/test_pipeline.R <scales, e.g. climate,landscape> <years, e.g. 2020:2021> [--algos glm,gam,rf,nn]
# Writes next to the BRT outputs of outputs/utest (new file names only; nothing existing is touched).
args <- commandArgs(trailingOnly = TRUE)
scales <- strsplit(if (length(args) >= 1) args[1] else "climate", ",")[[1]]
rng <- function(x) { p <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]]); if (length(p) == 1) p else p[1]:p[2] }
years <- rng(if (length(args) >= 2) args[2] else "2020:2021")
algos <- if ("--algos" %in% args) strsplit(args[which(args == "--algos") + 1], ",")[[1]] else c("glm", "gam", "rf")
Sys.setenv(BIRDMONITOR_RUNNAME = "utest", BIRDMONITOR_SPECIES = "Alauda arvensis", BIRDMONITOR_UNC_REPS = "0:3", BIRDMONITOR_UNC_YEARS = paste0(min(years), ":", max(years)))
suppressMessages({library(terra); library(gbm); library(glmnet); library(ranger); library(mgcv); library(nnet)})
source("modules/models_Monitor/tests/uncertainty/helper.R")
cfg <- testCfg(); sp <- cfg$species
ok <- function(label, cond) { cat(if (isTRUE(cond)) "ok   " else "FAIL ", label, "\n"); if (!isTRUE(cond)) quit(status = 1) }

# naming rule and registry
ok("ensTag: default set = 'ens'", ensTag(c("brt", "glm", "gam", "rf")) == "ens")
ok("ensTag: any other set is named by its members in canonical order", ensTag(c("gam", "brt")) == "ens_brt-gam" && ensTag("rf") == "ens_rf")
ok("registry has the families", all(c("glm", "gam", "rf", "nn") %in% names(algoRegistry())))

for (scale in scales) {
  cat("\n==== ", scale, " ====\n")
  for (a in c("brt", algos)) algoScaleRun(cfg, sp, scale, a, threads = 2L, years = years)
  for (a in algos) {
    ok(paste(a, "model, perf, cvpred, oofhab saved"), all(file.exists(vapply(c("model", "perf", "cvpred", "oofhab"), function(k) ensFile(cfg, sp, scale, a, k), ""))))
    ok(paste(a, "maps written for", length(years), "years"), all(file.exists(vapply(years, function(y) ensFile(cfg, sp, scale, a, "pred", y), ""))))
  }
  ok("brt member: out-of-fold files rebuilt from the fold models", all(file.exists(vapply(c("cvpred", "oofhab"), function(k) ensFile(cfg, sp, scale, "brt", k), ""))))
  tbl <- uncTrainingTable(cfg, sp, scale); ctx <- uncOofContext(cfg, sp)

  # several ensembles from the SAME fitted members: the full set, and subsets ("one out")
  sets <- list(c("brt", algos), "brt", algos, c("brt", algos[1]))
  for (members in sets) {
    ensembleScaleRun(cfg, sp, scale, members, years = years)
    tag <- ensTag(members)
    ok(paste(tag, ": perf/cvpred/oofhab saved"), all(file.exists(vapply(c("perf", "cvpred", "oofhab"), function(k) ensFile(cfg, sp, scale, tag, k), ""))))
    ok(paste(tag, ": cvpred has one value per training row; oofhab one per habitat record"),
       length(readRDS(ensFile(cfg, sp, scale, tag, "cvpred"))) == nrow(tbl) && length(readRDS(ensFile(cfg, sp, scale, tag, "oofhab"))) == length(ctx$rows))
    yr <- years[1]
    v <- sapply(vapply(members, function(a) ensFile(cfg, sp, scale, a, "pred", yr), ""), function(f) values(rast(f)[["mean_prob"]])[, 1])
    e <- values(rast(ensFile(cfg, sp, scale, tag, "pred", yr))[["mean_prob"]])[, 1]
    both <- stats::complete.cases(v) & !is.na(e)
    ok(paste(tag, ": map = mean of its members (", sum(both), "cells )"), max(abs(rowMeans(as.matrix(v))[both] - e[both])) < 1e-5)
    ok(paste(tag, ": NA where any member is NA"), all(is.na(e[!stats::complete.cases(v)])))
    if (length(members) > 1L) {
      sdv <- values(rast(ensFile(cfg, sp, scale, tag, "sd", yr))[[1]])[, 1]
      ok(paste(tag, ": SD layer = standard deviation across the members"), max(abs(sdv[both] - apply(as.matrix(v)[both, , drop = FALSE], 1, sd))) < 1e-5)
    }
  }
  ok("a single-member ensemble equals that model", isTRUE(all.equal(readRDS(ensFile(cfg, sp, scale, "ens_brt", "cvpred")), readRDS(ensFile(cfg, sp, scale, "brt", "cvpred")))))
  print(readRDS(ensFile(cfg, sp, scale, ensTag(c("brt", algos)), "members")), digits = 3)
}
cat("\nALL ENSEMBLE PIPELINE TESTS PASSED for", paste(scales, collapse = ", "), "\n")
