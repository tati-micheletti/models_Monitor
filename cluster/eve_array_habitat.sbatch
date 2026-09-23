#!/bin/bash
# SLURM job array: fit + cross-validate + predict the German habitat BRT,
# one task per species (11 total). This is the dominant cost of the whole
# pipeline -- locally it ranged ~57-154 min per species (worst case ~2h34m,
# ~15.7h total sequentially across all 11). The walltime/mem below give real
# margin on top of the worst observed case.
#
# TODO once your WKDV ticket answer arrives, fill in:
#   - #SBATCH --partition=<PARTITION>   (ask which partition suits a
#     single-core, <=16GB, <=4h job)
#   - the `module load` line(s) below, with the exact R module name EVE uses
#
# Submit from the birdMonitor repo ROOT (not from inside this cluster/ folder):
#   sbatch modules/models_Monitor/cluster/eve_array_habitat.sbatch

#SBATCH --job-name=birdmon-habitat
#SBATCH --array=1-11
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/habitat_%A_%a.out
#SBATCH --error=logs/habitat_%A_%a.err

set -euo pipefail
mkdir -p logs

# module load R/<version>   # uncomment once the exact module name is known

Rscript modules/models_Monitor/R/runClusterTask.R \
  --scale habitat \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --run-name test1 \
  --repo-root "${SLURM_SUBMIT_DIR}"
