#!/bin/bash
# NOT NEEDED IN NORMAL OPERATION. The meta-model now self-triggers: whichever
# of the europe/habitat/landscape arrays actually finishes last for a given
# species runs that species' meta-model itself, in the same job (see
# checkAllScalesReady() in models_Monitor.R) -- no manual --dependency
# chaining needed. Locally the meta step ranged ~144-183 min/species (worst
# case ~3h03m), comparable to or larger than the habitat-BRT bottleneck, so
# whichever scale task self-triggers it should keep this script's own
# walltime/memory budget, not the scale array's.
#
# This script exists only to catch STRAGGLERS: if a species' scale task
# failed (so meta never self-triggered for it), submit this for just that
# species' index after fixing the underlying failure -- e.g.:
#   sbatch --array=7 modules/models_Monitor/cluster/eve_array_meta.sbatch
# It safely no-ops (not an error) for any species whose 3 scales aren't
# all actually done yet.
#
# TODO once your WKDV ticket answer arrives, fill in:
#   - #SBATCH --partition=<PARTITION>   (ask which partition suits a
#     single-core, <=16GB, <=4h job)
#   - the `module load` line(s) below, with the exact R module name EVE uses
#
# Submit from the birdMonitor repo ROOT (not from inside this cluster/
# folder), restricted to the specific species that need it (--array=<index>,
# e.g. "3" or "3,7,9" -- see straggler example above).

#SBATCH --job-name=birdmon-meta
#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/meta_%A_%a.out
#SBATCH --error=logs/meta_%A_%a.err

set -euo pipefail
mkdir -p logs

# module load R/<version>   # uncomment once the exact module name is known

Rscript tools/runClusterTask.R \
  --scale meta \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --run-name test1 \
  --repo-root "${SLURM_SUBMIT_DIR}"
