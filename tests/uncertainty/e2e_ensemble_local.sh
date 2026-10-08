#!/bin/bash
# Local end-to-end rehearsal of the uncertainty workflow with the ENSEMBLE inside every replicate (members brt, glm, nn: the cheap ones).
# Needs: step0_mini_baseline.R and tests/ensemble/test_meta_ens.R first (they make the baseline ensemble, which replicate 0 must reproduce).
#   bash modules/models_Monitor/tests/uncertainty/e2e_ensemble_local.sh
set -euo pipefail
cd "$(dirname "$0")/../../../.."
R="${RSCRIPT:-/c/Users/michelet/AppData/Local/Programs/R/R-4.6.0/bin/Rscript.exe}"
export BIRDMONITOR_RUNNAME=utest BIRDMONITOR_SPECIES="Alauda arvensis" BIRDMONITOR_UNC_REPS=0:2 BIRDMONITOR_UNC_MEMBERS=brt,glm,nn \
       BIRDMONITOR_UNC_YEARS=2024:2025 BIRDMONITOR_UNC_BANDS=40 BIRDMONITOR_UNC_BASELINE=2024 BIRDMONITOR_UNC_REPBATCH=3
unset BIRDMONITOR_UNC_TAG || true
t() { echo "=== $(date +%H:%M:%S) $*"; }
t preflight;   $R tools/runUncertaintyTask.R --step preflight
for s in fit coarse oof ridge; do t $s; $R tools/runUncertaintyTask.R --step $s --index 1; done
t bandpredict; seq 1 40 | xargs -P 6 -I{} $R tools/runUncertaintyTask.R --step bandpredict --index {} > /dev/null
t summarize;   seq 1 40 | xargs -P 6 -I{} $R tools/runUncertaintyTask.R --step summarize --index {} > /dev/null
t assemble;    $R tools/runUncertaintyTask.R --step assemble --index 1
t community;   seq 1 40 | xargs -P 6 -I{} $R tools/runUncertaintyTask.R --step community --index {} > /dev/null
t assembleAll; $R tools/runUncertaintyTask.R --step assembleAll
t "DONE e2e ensemble"
