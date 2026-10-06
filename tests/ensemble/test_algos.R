# Local test of the alternative algorithms (GLM, GAM, RF) on the mini baseline's training tables (Alauda arvensis).
#   Rscript modules/models_Monitor/tests/ensemble/test_algos.R [scale: climate|landscape|habitat, default climate]
args <- commandArgs(trailingOnly = TRUE); scale <- if (length(args)) args[1] else "climate"
Sys.setenv(BIRDMONITOR_RUNNAME = "utest", BIRDMONITOR_SPECIES = "Alauda arvensis", BIRDMONITOR_UNC_REPS = "0:3")
suppressMessages({library(terra); library(gbm); library(glmnet); library(randomForest); library(mgcv)})
source("modules/models_Monitor/tests/uncertainty/helper.R")
cfg <- testCfg(); sp <- cfg$species
d <- uncTrainingTable(cfg, sp, scale); brt <- uncMainModel(cfg, sp, scale); predSel <- brt$var.names
cat(scale, ":", nrow(d), "records,", length(predSel), "predictors, folds:", paste(sort(unique(d$foldID)), collapse = ","), "\n")
ok <- function(label, cond) { cat(if (isTRUE(cond)) "ok   " else "FAIL ", label, "\n"); if (!isTRUE(cond)) quit(status = 1) }

res <- list(); cvs <- list()
for (a in c("glm", "gam", "rf")) {
  t0 <- Sys.time(); f <- algoFit(a, d, predSel, seed = 1L); tf <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  p <- algoPredict(f, d)
  ok(paste(a, "predicts probabilities for every row"), length(p) == nrow(d) && all(p >= 0 & p <= 1) && !anyNA(p))
  t0 <- Sys.time(); cv <- blockCVPredictAlgo(f, d, d$foldID); tc <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  ok(paste(a, "block-CV covers every row"), !anyNA(cv))
  pf <- evalSDM(d$occurrence, cv)
  cat(sprintf("     %s: %d predictors kept, fit %.1fs, CV %.1fs | block-CV AUC %.3f TSS %.3f D2 %.3f %s\n", a, length(f$predSel), tf, tc, pf$AUC, pf$TSS, pf$D2, f$note))
  res[[a]] <- f; cvs[[a]] <- cv
}
# refit on the same data reproduces the model's predictions (glm/gam deterministic)
for (a in c("glm", "gam")) ok(paste(a, "algoRefit on the same data = same predictions"), max(abs(algoPredict(algoRefit(res[[a]], d), d) - algoPredict(res[[a]], d))) < 1e-8)
# BRT block-CV for comparison, and the ensemble
cvB <- blockCVPredictBRT(d, predSel, "occurrence", d$foldID, brt); attr(cvB, "foldModels") <- NULL
ens <- algoEnsembleMean(list(cvs$glm, cvs$gam, cvs$rf, cvB))
pB <- evalSDM(d$occurrence, cvB); pE <- evalSDM(d$occurrence, ens)
cat(sprintf("     BRT block-CV AUC %.3f TSS %.3f | ensemble (glm+gam+rf+brt) AUC %.3f TSS %.3f D2 %.3f\n", pB$AUC, pB$TSS, pE$AUC, pE$TSS, pE$D2))
ok("ensemble mean is the arithmetic mean", isTRUE(all.equal(ens, rowMeans(cbind(cvs$glm, cvs$gam, cvs$rf, cvB)))))
ok("ensemble is NA where an algorithm is NA", is.na(algoEnsembleMean(list(c(0.2, NA), c(0.4, 0.5)))[2]))
cat("ALL ALGORITHM TESTS PASSED\n")
