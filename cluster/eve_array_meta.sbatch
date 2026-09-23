#!/bin/bash
# SLURM job array: combine the europe/habitat/landscape scales via the ridge
# regression meta-model and predict onto the full period, one task per
# species (11 total). Locally this ranged ~144-183 min per species (worst
# case ~3h03m) -- comparable to or larger than the habitat-BRT bottleneck,
# so don't assume the fitting stage is the only expensive one.
#
# STAGE DEPENDENCY: each species' meta task needs that species' own
# europe/habitat/landscape tasks already finished (it reads their output
# directories). This array must be submitted with a dependency on the three
# scale arrays, e.g.:
#   jid_eu=$(sbatch --parsable modules/models_Monitor/cluster/eve_array_europe.sbatch)
#   jid_ha=$(sbatch --parsable modules/models_Monitor/cluster/eve_array_habitat.sbatch)
#   jid_la=$(sbatch --parsable modules/models_Monitor/cluster/eve_array_landscape.sbatch)
#   sbatch --dependency=afterok:${jid_eu}:${jid_ha}:${jid_la} \
#     modules/models_Monitor/cluster/eve_array_meta.sbatch
#
# TODO once your WKDV ticket answer arrives, fill in:
#   - #SBATCH --partition=<PARTITION>   (ask which partition suits a
#     single-core, <=16GB, <=4h job)
#   - the `module load` line(s) below, with the exact R module name EVE uses
#
# Submit from the birdMonitor repo ROOT (not from inside this cluster/ folder,
# and not on its own -- see the dependency example above):
#   sbatch modules/models_Monitor/cluster/eve_array_meta.sbatch

#SBATCH --job-name=birdmon-meta
#SBATCH --array=1-11
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
