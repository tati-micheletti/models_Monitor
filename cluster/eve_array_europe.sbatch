#!/bin/bash
# SLURM job array: fit + cross-validate + predict the European climate BRT,
# one task per species (11 total). Cheapest of the 3 scales by far --
# locally this took ~10-15 seconds PER SPECIES, so the resource requests
# below are deliberately generous, not tight.
#
# TODO once your WKDV ticket answer arrives, fill in:
#   - #SBATCH --partition=<PARTITION>   (ask which partition suits a
#     single-core, <=4GB, <15min job)
#   - the `module load` line(s) below, with the exact R module name EVE uses
#
# Submit from the birdMonitor repo ROOT (not from inside this cluster/ folder):
#   sbatch modules/models_Monitor/cluster/eve_array_europe.sbatch

#SBATCH --job-name=birdmon-europe
#SBATCH --array=1-11
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/europe_%A_%a.out
#SBATCH --error=logs/europe_%A_%a.err

set -euo pipefail
mkdir -p logs

# module load R/<version>   # uncomment once the exact module name is known

Rscript modules/models_Monitor/R/runClusterTask.R \
  --scale europe \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --run-name test1 \
  --repo-root "${SLURM_SUBMIT_DIR}"
