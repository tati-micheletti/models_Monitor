#!/bin/bash
# SLURM job array: fit + cross-validate + predict the German landscape BRT,
# one task per species (12 total -- keep this in sync with length(sharedSpecies) in sharedConfig.R). Locally this ranged ~16-71 min per
# species (worst case ~71 min). Walltime below also covers the case where
# this task turns out to be the one that finishes last for its species and
# self-triggers that species' meta-model (~144-183 min -- see
# checkAndRunMeta logic in models_Monitor.R).
#
# TODO once your WKDV ticket answer arrives, fill in:
#   - #SBATCH --partition=<PARTITION>   (ask which partition suits a
#     single-core, <=16GB, <=5h job)
#   - the `module load` line(s) below, with the exact R module name EVE uses
#
# Submit from the birdMonitor repo ROOT (not from inside this cluster/ folder):
#   sbatch modules/models_Monitor/cluster/eve_array_landscape.sbatch

#SBATCH --job-name=birdmon-landscape
#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=05:00:00
#SBATCH --output=logs/landscape_%A_%a.out
#SBATCH --error=logs/landscape_%A_%a.err

set -euo pipefail
mkdir -p logs

# module load R/<version>   # uncomment once the exact module name is known

Rscript tools/runClusterTask.R \
  --scale landscape \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --run-name test1 \
  --repo-root "${SLURM_SUBMIT_DIR}"
