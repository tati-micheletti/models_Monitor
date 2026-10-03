#!/bin/bash -l
# SLURM job array: fit + cross-validate + predict the European climate BRT,
# one task per species (11 total -- keep this in sync with length(sharedSpecies) in sharedConfig.R, i.e. the number of include=X rows in data/speciesCanonical.csv). Cheapest of the 3 scales by far --
# locally this took ~10-15 seconds PER SPECIES. Walltime below is sized for
# the case where this task also turns out to be the one that finishes last
# for its species and self-triggers that species' meta-model (~144-183 min
# -- see checkAndRunMeta logic in models_Monitor.R), not for europe's own
# (near-instant) cost alone.
#
# TODO before the first submission, fill in (both can be looked up on EVE itself):
#   - #SBATCH --partition=<PARTITION>   (run `sinfo -s` on a login node and pick a partition that suits a
#     single-core, <=16GB, <=4h job)
#   - the `module load` line(s) below, with the exact R module name from `module spider R`
#
# Submit from the birdMonitor repo ROOT (not from inside this cluster/ folder):
#   sbatch modules/models_Monitor/cluster/eve_array_europe.sbatch

#SBATCH --job-name=birdmon-europe
#SBATCH --array=1-11
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=13G
#SBATCH --time=04:00:00
#SBATCH --output=logs/europe_%A_%a.out
#SBATCH --error=logs/europe_%A_%a.err

set -euo pipefail
mkdir -p logs

# EVE_R_MODULE is set at submit time (see cluster/submit_eve_pipeline.sh) -- the exact
# name from `module spider R`. Nothing to edit in this file.
# Module list is set HERE (not only by the submit script) so the job never depends on
# EVE forwarding environment variables. Override at submit time with EVE_R_MODULE=... if needed.
EVE_R_MODULE="${EVE_R_MODULE:-GCC/13.3.0 OpenMPI/5.0.5 R/4.5.1 GDAL/3.10.3 CMake ImageMagick/7.1.1-38 UDUNITS/2.2.28}"
module load ${EVE_R_MODULE}
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR on $(hostname): Rscript not found after: module load ${EVE_R_MODULE}" >&2
  module list >&2
  exit 1
fi

# $TMPDIR on EVE is a RAM disk (half of the node's RAM, shared with the job itself),
# so terra's large temporary rasters would exhaust memory. Use /work instead and
# clean up when the job ends (EVE wiki, Storage page).
export TMPDIR="/work/${USER}/birdMonitor/tmp/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
trap 'rm -rf "${TMPDIR}"' EXIT

Rscript tools/runClusterTask.R \
  --scale europe \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --run-name "${BIRDMONITOR_RUNNAME:-test4}" \
  --repo-root "${SLURM_SUBMIT_DIR}"
