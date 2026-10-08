# Uncertainty (option B): confidence layers from a spatial-block bootstrap of the BRTs

What this adds to the baseline run: for every species, **maps of how uncertain the predicted probability of
occurrence is**, **intervals on the index time series**, and **per-pixel confidence layers for change and trend** — all
from the same models, refitted on resampled data. Nothing in the baseline results is touched; everything is written under
`outputs/<runName>/uncertainty/` (the index intervals also next to the annual report, in `outputs/<runName>/annual_report/`).

Plan and decisions: `improvements.md` item 15 and `DECISIONS.md` (2026-10-06). Plain-language explanation for readers
outside the project: the "uncertainty" document (A) of the document set.

## The method in six lines

0. **Replicate 0** is the main models themselves, pushed through exactly the same machinery. It must reproduce the baseline
   meta-model map (`parity_check.txt` says whether it does); it is never part of an interval.
1. Cut the training data of each scale (climate, landscape, habitat) into **spatial blocks** (size = the cross-validation
   block size of that scale, `cv_spatial_autocor()` rule).
2. One **replicate** = draw blocks *with replacement* until as many blocks as the original; all records of a drawn block
   enter (twice if drawn twice). Redraw, deterministically, if a draw has fewer than 10 presences or absences.
3. Refit each scale's BRT on the draw with the **fixed hyperparameters** of the main model (trees, learning rate, bag
   fraction, tree complexity). Hyperparameters are not re-tuned (as in the Boreal Avian Modelling Project).
4. Fit the **ridge meta-model** of that replicate (weights never negative) on the replicate's own scale predictions at
   the habitat records (training rows = the same habitat draw). The predictions are **out-of-fold**: each record is
   predicted by a refit of the replicate's BRT that never saw the record's block (the replicate's draw minus the record's
   block-CV fold, same fixed hyperparameters), so the weights are honest, as in the baseline (`uncOofSpecies()`, step `oof`).
   Replicate 0 uses the saved fold models of the main BRTs, i.e. it equals the baseline.
5. Predict **all years** with the replicate's three BRTs + ridge, resampling the coarse scales to the 200 m grid exactly
   like the baseline.
6. Compute everything derived (change, trend, area mean, index) **inside each replicate**, then take percentiles
   across replicates. This is what carries spatial *and* temporal uncertainty into the derived layers.

Seeds are fixed per species / scale / replicate id (`stableSeed()` of a key), logged in `replicate_log.csv`, and the
caller's RNG state is restored. Adding replicates later (ids 51–100) never changes replicates 1–50.

## Outputs (per species, under `outputs/<runName>/uncertainty/<Species_name>/`)

| File | What it is |
|---|---|
| `maps/<sp>_unc_<year>.tif` | 5 layers per year: `mean`, `sd`, `lwr` (5th pct), `upr` (95th pct), `width` (= upr − lwr) of the probability of occurrence across replicates. **`width` is the confidence layer: wide = uncertain.** |
| `maps/<sp>_unc_change_vsBaseline.tif`, `…_vs5YearsAgo.tif`, `…_vsLastYear.tif` | change in probability, current year minus reference year, computed per replicate. 7 layers: `deltaMean`, `deltaSd`, `deltaLwr`, `deltaUpr`, `deltaWidth`, **`shareDecrease`** (share of replicates with a decrease), **`shareIncrease`**. Reference years: 2005, current−5, current−1; current = 2025. |
| `maps/<sp>_unc_trend_per_decade.tif` | per-pixel linear trend over all mapped years (change in probability per decade), per replicate. 7 layers: `slopeMean`, `slopeSd`, **`slopeLwr`, `slopeUpr`** (the 90% interval of the slope), `slopeWidth`, `shareDecrease`, `shareIncrease`. |
| `species_index_uncertainty.csv` (annual_report/, all species) | per species and year: mean, median, 90% interval of the species' area-mean probability and of the **index** (100 × area mean ÷ baseline-year area mean). |
| `area_mean_replicates.csv` | the raw numbers behind it: area mean of every replicate and year (long format). |
| `replicate_log.csv` | one row per fitted model: seeds, attempts, blocks, rows drawn, presences/absences, hyperparameters, git commits, time. |
| `oof/oof_<run>.rds` | out-of-fold scale predictions at the habitat records, per replicate (the inputs the weights are trained on). |
| `ridge/ridge_<run>.rds` | ridge coefficients per replicate and how many records each used. |

Across species (`outputs/<runName>/uncertainty/`):

| File | What it is |
|---|---|
| `combined_index_uncertainty.csv` (annual_report/) | the combined multi-species indices per year, each with a 90% interval: **SBI** (geometric mean of the species indices) and **Chain** (LPI-style) computed inside each replicate, **from the same replicate in every species**; and **Analytical** and **MSI** (the two DDA/PECBMS methods), which use the SD of each species index across replicates as their standard error. |
| `combined_index_replicates.csv` | the per-replicate series behind it. |
| `community/richness_expected_unc_<year>.tif` | expected species richness (sum of probabilities) with the 5-layer summary. |
| `community/community_meanDeltaP_unc_<comparison>.tif` | mean change in probability over species, 7 layers as above. |
| `community/community_turnoverBC_unc_<comparison>.tif` | **Bray-Curtis turnover** between the two years, computed per replicate over all species and then summarised: 5 layers (`mean`, `sd`, `lwr`, `upr`, `width`; no share layers, a dissimilarity has no sign). 0 = the expected community is unchanged, 1 = no species in common. See the note below. |

Why a turnover layer: expected richness (sum of probabilities) and `community_meanDeltaP` (mean of the signed changes) both cancel under a swap of
species. Example: 5 species at 0.8 in 2005; in 2025 four are lost and four others gained -> richness identical to "nothing changed", mean change
0. Bray-Curtis, BC = sum_s |p_2025 - p_ref| / sum_s (p_ref + p_2025) = 1 - 2 sum_s min(p_ref, p_2025) / (sum_s p_ref + sum_s p_2025), is 0.8 for the swap
and 0 for no change. The probabilities of occurrence stand in for abundances; BC also rises when richness changes without any swap
(gain/loss of species are both turnover and nestedness in Baselga's terms; they are not separated here). Baseline layer: `turnoverBC`, the 5th layer
of `change_<comparison>_community.tif`.

How to read the confidence layers: the interval is a **percentile interval of the fitted probability surface** (how much
the model would change if the training data had been different), **not** a prediction interval for observed occurrences.
`shareDecrease = 0.95` means 95% of the replicates predict a decrease at that pixel.

## Running it on EVE

After the baseline (prep → model arrays → index) has finished:

```bash
cd ~/projects/birdMonitor
git pull
bash --login cluster/setup_uncertainty.sh          # once: installs matrixStats, lists the key packages
BIRDMONITOR_RUNNAME=test4 bash --login cluster/submit_eve_uncertainty.sh
```

The script first runs a pre-flight check (every input must exist; it lists what is missing and stops), then submits the
whole chain with dependencies. Useful settings (environment variables in front of the command):

| Variable | Default | Use |
|---|---|---|
| `BIRDMONITOR_UNC_REPS` | `0:50` | replicate ids of this run; **0 is the main models (no resampling), run through the same machinery as a built-in check** (`parity_check.txt` per species) and left out of every interval. **To add replicates: `51:100`** (new run label; old files stay). |
| `BIRDMONITOR_UNC_YEARS` | all prediction years | e.g. `2005,2020:2025` for a quick first pass |
| `BIRDMONITOR_UNC_BANDS` | `16` | the country is cut into this many horizontal bands (more = smaller jobs) |
| `UNC_THROTTLE` | none | maximum band tasks running at the same time |
| `UNC_AFTER` | none | job id(s) to wait for first |
| `BIRDMONITOR_UNC_BLOCKMULT` | `1` | multiplies the resampling block size (sensitivity test: 0.5 and 2) |

## First run on EVE: a cheap timing/test run (do this before the real one)

Replicate 0 plus three replicates, one year, written to a separate folder (`uncertainty_timing`) that can never mix with real replicates.
It runs every step on real data, so it also shows any problem the local tests could not:

```bash
BIRDMONITOR_RUNNAME=test4 BIRDMONITOR_UNC_TAG=timing BIRDMONITOR_UNC_REPS=0:3 BIRDMONITOR_UNC_YEARS=2025   UNC_THROTTLE=40 bash --login cluster/submit_eve_uncertainty.sh
```

Afterwards, `sacct -j <jobid> --format=JobID,JobName,Elapsed,MaxRSS,State` per step gives the real minutes per replicate and
year; the cost of the real run is then (minutes per replicate-year) × replicates × years. The covariate cache this run
builds is re-used by the real run.

## Cost and storage (B = 50, 11 species, 21 years) — estimates, to be replaced by the timing run

Prediction of the 200 m habitat model dominates. About 770 core-hours in total; wall time = core-hours ÷ cores used.
Per-replicate predictions are stored as 16-bit integers (resolution 3.3 × 10⁻⁵): roughly 8–12 GB per species, ~100 GB for
all species, under `…/<species>/pred/<run label>/`. They are what makes adding replicates later possible; delete them
once the final replicate set is fixed if space is needed.

## Regional index (10 / 20 / 50 km) with uncertainty

The regional index of the baseline (`computeRegionalIndex()`) is a per-cell version of the combined index. The uncertainty version does the same
inside every replicate and then summarises over replicates, on the RAW grid (smoothing is applied only to the final raw result, never inside
replicates; the smoothed baseline map therefore has no interval of its own).

1. **Masking before averaging.** The 200 m pixels are masked with the German outline (GADM level 0, as the baseline) BEFORE they are averaged into
   the coarse cells, so a border cell holds German pixels only. Nothing in the models or predictions is cropped.
2. **Cell means per replicate** (`R/uncRegional.R`, steps `regionband` per species x band and `regionassemble` per species): for each species, year and
   replicate, the mean probability of the German pixels of every coarse cell, stored as `<species>/regional/regional_means_<km>km.rds`.
3. **Index per replicate** (`computeRegionalIndexUncertainty()` in runIndex_Monitor, step `regionindex`): species index = 100 x cell mean in year t /
   cell mean in the baseline year (a species is left out of a cell where its baseline mean is below `minBaseline`, 10^-6 as in the baseline), combined
   index = geometric mean over species. **The floor (and a possible DDA-style cap) must be the same in the baseline and here; the decision is open
   (DECISIONS.md 2026-10-08).**
4. **Summaries over the replicates** (replicate 0 is never part of an interval): per year mean, sd, lwr, upr, width; change between years (3 comparisons
   as for the species maps) and trend per decade with `shareDecrease` / `shareIncrease` (cold and hot spots); a coverage layer (share of the cell's pixels
   inside the outline); a parity check of replicate 0 against the baseline `regional_index_<km>km_raw.tif` (`regional_parity_<km>km.txt`).

| File (`uncertainty_<tag>/regional/`) | What it is |
|---|---|
| `regional_index_<km>km_unc_<year>.tif` | 5 layers: `mean`, `sd`, `lwr`, `upr`, `width` of the combined regional index |
| `regional_index_<km>km_unc_change_<vsBaseline\|vs5YearsAgo\|vsLastYear>.tif` | 7 layers: `deltaMean`, `deltaSd`, `deltaLwr`, `deltaUpr`, `deltaWidth`, `shareDecrease`, `shareIncrease` |
| `regional_index_<km>km_unc_trend_per_decade.tif` | 7 layers: `slopeMean` ... `slopeWidth`, `shareDecrease`, `shareIncrease` (index points per decade) |
| `regional_index_<km>km_coverage.tif` | share of the cell inside the German outline (border cells < 1) |
| `regional_index_<km>km_replicates.rds`, `regional_parity_<km>km.txt` | the per-replicate index array behind the maps; the parity check |

Run on EVE after the band stage: `BIRDMONITOR_UNC_TAG=<tag> UNC_AFTER=<band job id> bash --login cluster/submit_eve_regional.sh` (nothing is refitted;
it reads the stored replicate predictions). Local test: replicate 0, built from its own pixels, reproduces the baseline regional machinery to 0.002
index points (10 km), including the border cells.

## Germany-only national means and maps (temporary, DECISIONS.md 2026-10-08)

The prediction window is the bounding box of Germany; for seven species it is predicted entirely, so `area_mean_replicates.csv` (all pixels) averages neighbouring
countries too. `regionassemble` therefore also writes `<species>/area_mean_replicates_germany.csv` (German pixels only: sum of the cell sums / number of pixels,
per replicate and year); `assembleAll` with `BIRDMONITOR_UNC_AREAMEAN=area_mean_replicates_germany.csv BIRDMONITOR_INDEX_TAG=germany` writes the species and
combined index intervals from it to `annual_report_germany/`. The `maskmaps` step writes Germany-only COPIES of the finished maps (`<species>/maps_germany/`,
`community_germany/`); the originals stay. All of this is chained in `cluster/submit_eve_regional.sh`. The proper fix (cut every input to the study-area outline so
nothing is predicted outside Germany) is planned after the Steering Consortium meeting.

## How to read the numbers — limits that belong in the methods

- The `mean` layer is the mean over replicates, not the baseline map (bootstrap averaging shifts it slightly). Use the
  baseline map as the point estimate and these layers as its uncertainty; replicate 0 shows they describe the same model.
- Percentiles from 50 replicates are noisy in the tails (the 5th and 95th percentile each rest on about 2–3 replicates); the
  widths will tighten and stabilise when replicates are added (B = 100).
- Years before the habitat training years (2005–2021 here) are hindcasts that assume the 2022–2025 relationships hold;
  the intervals do **not** include that extrapolation uncertainty, and per-pixel "trends" over 2005–2025 are driven by
  how the covariates changed, not by independent evidence about the birds.
- As in the baseline, each scale's BRT is evaluated at its own training records when the meta-model is fitted
  (in-sample suitability); the bootstrap reproduces that choice, it does not correct for it.
- `shareDecrease` / `shareIncrease` are shares of replicates, not probabilities in a Bayesian sense, and have no threshold:
  a change of −0.0001 counts as a decrease. Read them together with `deltaLwr` / `deltaUpr`.
- The climate scale is resampled in blocks of the cross-validation block size, capped at 1500 km: about 13 blocks for Europe.

## What is NOT covered (to state in the methods)

1. Hyperparameter-tuning uncertainty (settings fixed from the main fit).
2. Error in the predictor layers (land use, climate, DEM treated as exact).
3. Structural / model-form uncertainty (one algorithm per scale; the algorithm ensemble is planned, improvements.md item 7/15).
4. Survey-design bias, and spatial/temporal autocorrelation beyond what block resampling captures.
5. The thinning randomness of the occurrence data (thinning is part of the data preparation, not resampled here).
6. The SMOOTHED regional index maps (smoothing is applied to the final raw result only; the raw regional index has intervals, see above).

## Files

Code (SpaDES events of **models_Monitor**, selected with `runScale`): `R/uncCommon.R` (config, paths, bootstrap draw/fit),
`R/uncFit.R` (replicate fits, covariate cache, coarse predictions), `R/uncBand.R` (band prediction, ridge), `R/uncSummarize.R`
(summaries, assembly), `R/uncSim.R` (the events), `R/uncPreflight.R`. The index intervals are an event of **runIndex_Monitor**
(`computeIndexUncertainty`, parameter `uncertaintyDir`). Entry point `tools/runUncertaintyTask.R`; SLURM scripts
`cluster/eve_unc_*.sbatch`, chain `cluster/submit_eve_uncertainty.sh`; local tests `tests/uncertainty/`.
