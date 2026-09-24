#!/bin/bash
# SLURM job array: fit + cross-validate + predict the European climate BRT,
# one task per species (12 total -- keep this in sync with length(sharedSpecies) in sharedConfig.R). Cheapest of the 3 scales by far --
# locally this took ~10-15 seconds PER SPECIES. Walltime below is sized for
# the case where this task also turns out to be the one that finishes last
# for its species and self-triggers that species' meta-model (~144-183 min
# -- see checkAndRunMeta logic in models_Monitor.R), not for europe's own
# (near-instant) cost alone.
#
# TODO once your WKDV ticket answer arrives, fill in:
#   - #SBATCH --partition=<PARTITION>   (ask which partition suits a
#     single-core, <=16GB, <=4h job)
#   - the `module load` line(s) below, with the exact R module name EVE uses
#
# Submit from the birdMonitor repo ROOT (not from inside this cluster/ folder):
#   sbatch modules/models_Monitor/cluster/eve_array_europe.sbatch

#SBATCH --job-name=birdmon-europe
#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/europe_%A_%a.out
#SBATCH --error=logs/europe_%A_%a.err

set -euo pipefail
mkdir -p logs

# module load R/<version>   # uncomment once the exact module name is known

Rscript tools/runClusterTask.R \
  --scale europe \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --run-name test1 \
  --repo-root "${SLURM_SUBMIT_DIR}"
