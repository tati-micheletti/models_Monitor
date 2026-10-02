#!/bin/bash -l
# SLURM job array: fit + cross-validate + predict the German landscape BRT,
# one task per species (11 total -- keep this in sync with length(sharedSpecies) in sharedConfig.R, i.e. the number of include=X rows in data/speciesCanonical.csv). Locally this ranged ~16-71 min per
# species (worst case ~71 min). Walltime below also covers the case where
# this task turns out to be the one that finishes last for its species and
# self-triggers that species' meta-model (~144-183 min -- see
# checkAndRunMeta logic in models_Monitor.R).
#
# TODO before the first submission, fill in (both can be looked up on EVE itself):
#   - #SBATCH --partition=<PARTITION>   (run `sinfo -s` on a login node and pick a partition that suits a
#     single-core, <=16GB, <=5h job)
#   - the `module load` line(s) below, with the exact R module name from `module spider R`
#
# Submit from the birdMonitor repo ROOT (not from inside this cluster/ folder):
#   sbatch modules/models_Monitor/cluster/eve_array_landscape.sbatch

#SBATCH --job-name=birdmon-landscape
#SBATCH --array=1-11
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=13G
#SBATCH --time=05:00:00
#SBATCH --output=logs/landscape_%A_%a.out
#SBATCH --error=logs/landscape_%A_%a.err

set -euo pipefail
mkdir -p logs

# EVE_R_MODULE is set at submit time (see cluster/submit_eve_pipeline.sh) -- the exact
# name from `module spider R`. Nothing to edit in this file.
if [ -n "${EVE_R_MODULE:-}" ]; then module load ${EVE_R_MODULE}; fi   # may hold several names: "GCC/15.3.0 R/4.6.1" (unquoted on purpose)

Rscript tools/runClusterTask.R \
  --scale landscape \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --run-name "${BIRDMONITOR_RUNNAME:-test4}" \
  --repo-root "${SLURM_SUBMIT_DIR}"
