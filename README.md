# Somatic compartment-network simulations of module-integrity erosion

This repository contains the individual-based simulations supporting the
somatic (within-body) channel of module-integrity erosion described in
the manuscript, together with the generator of the manuscript's
conceptual landscape figure (Fig. 2, a pure schematic with no simulation
output). The model follows the irreversible loss-of-function
mutation of a bottleneck core gene in a mycelial compartment network:
a 12 × 12 lattice of compartments (periodic boundaries) holding K = 8
nuclei each (N_tot = 1,152), nuclear migration through hyphal fusion at
four network architectures (C0 fragmented; C1 local fusion; C2 local
plus rare long-range anastomosis; C3 strong long-range fusion),
compartment turnover (γ = 0.25 per nuclear generation) implementing
between-compartment selection through donor weighting 1 − S·x_d
(floored at 0.02), and an optional selfish within-compartment advantage
w of mutant nuclei. Each nuclear generation applies within-compartment
sampling, migration, turnover and mutation, in that order; migration
pools a four-neighbour local contribution (weight m_loc) with a
whole-body contribution (weight m_glob). Mutation, migration, turnover
and selection act on a somatic timescale and are deliberately decoupled
from any population-level process.

The simulations are conceptual demonstrations of mechanism, not
parameterised fits to any particular taxon. All parameter values are
given in the figure captions in `figures/`.

## Contents

| Path | Contents |
| --- | --- |
| `code/somatic_ibm.py` | Simulation engine (CLI: `calib`, `selresp`, `fixprob`, `erosion`, `disconnect`, `snapshots`, `batch`, `smoke`) |
| `code/make_figures.py` | Builds the manuscript figures and `somatic_ibm_summary.json` from `data/ibm_results/` |
| `code/verify_numbers.py` | Checks every number quoted in the figure captions against `somatic_ibm_summary.json` and the sensitivity summary |
| `code/sensitivity_scan.py` | One-factor-at-a-time sensitivity scan producing `data/sensitivity/` |
| `code/make_fig_landscape.py` | Draws the conceptual landscape-continuum figure (Fig. 2); no data dependency |
| `data/ibm_results/` | Archived simulation output (45 `.npz` files; key-level dictionary in `data/README.md`) |
| `data/sensitivity/` | Sensitivity scan output (120 jobs; protocol and result tables in `data/sensitivity/README.md`) |
| `figures/` | Manuscript figures in editable vector formats (SVG and PDF) with LaTeX caption sources and the numerical summary |

## Scope and limitations

- Parameter values (lattice size, K, γ, m_loc, m_glob, μ, S, w) were
  chosen to keep the four architectures distinguishable within tractable
  run times; none is estimated from data.
- Quantitative emphasis is on endpoint contrasts between the fragmented
  (C0) and strongly mixed (C3) architectures. The intermediate
  architectures C1 and C2 exchange rank between selective regimes, so no
  monotone connectivity effect should be inferred from their ordering.
- Panel a of the main figure is an illustrative schematic; panels b–e
  are simulation output. Panel b displays one representative replicate
  per architecture (the replicate closest to the ensemble-mean burden at
  the common checkpoint t = 1,800); the sectoring statistics quoted in
  its caption are medians across all polymorphic replicates
  (0.02 < x̄ < 0.98) at that checkpoint, with completed or nearly clean
  replicates excluded and their number reported.
- Fixation probabilities are resolved down to two events in R = 50,000
  replicates (4 × 10⁻⁵); conditions with no observed fixation are
  plotted at that detection limit as downward triangles.
- Erosion-time medians are computed over completed replicates only and
  are right-censored at the run horizon for up to 22.5% of replicates
  (costly regime, C0); at R = 200 they also carry appreciable
  Monte-Carlo error (seed-to-seed standard deviations of a few per cent
  to about 15%; see the provenance notes below).
- The donor-weighting rule 1 − S·x_d and the selfish advantage w are
  phenomenological devices chosen for transparency, not derived
  mechanistically.
- A one-factor-at-a-time sensitivity scan (each ecological parameter
  halved and doubled) is archived in `data/sensitivity/`: the direction
  of every endpoint contrast is preserved in all scanned variants,
  while the multiplicative factors are baseline-specific (see
  "Sensitivity analysis" below).
- The landscape-continuum figure (Fig. 2) is an illustrative schematic
  with no data dependency: basin positions, depths and widths are
  chosen for legibility, and the topology of its panel b is drawn, not
  measured. In panel b the dashed ring around each node is drawn at the
  basin-width radius used for the corresponding Gaussian basin of
  panel a, so the two panels stay geometrically consistent.

## Requirements

Python 3 (tested with 3.12), NumPy and Matplotlib; pinned versions in
`requirements.txt`.

```
pip install -r requirements.txt
```

## Reproducing the archived results

From the repository root:

```
python code/somatic_ibm.py batch --outdir data/ibm_results
python code/make_figures.py
python code/make_fig_landscape.py
python code/verify_numbers.py
```

`batch` executes the full job list (45 jobs: 4 neutral calibrations,
4 selection-response calibrations, 20 fixation-probability runs at
R = 50,000 replicates each, 12 erosion runs, 1 disconnection experiment
and 4 snapshot runs) and writes one `.npz` per job. Existing files are
skipped, so with the archived data in place the batch simply confirms
that all 45 outputs are present; to regenerate the archive from scratch,
first move the `.npz` files out of `data/ibm_results/`, or point
`--outdir` at an empty directory and copy the results back afterwards.
`make_figures.py` rebuilds the figures and the numerical summary in
`figures/` from `data/ibm_results/`; `make_fig_landscape.py` redraws the
conceptual Fig. 2 (no data input, runs in seconds); `verify_numbers.py`
cross-checks every number quoted in the captions and should print
`PASS`.

Runtime (single core, Linux, NumPy 2.1.3): a complete batch of all 45
jobs takes about 9 minutes. The 20 fixation-probability jobs
(R = 50,000 each) take 5–8 seconds apiece; the erosion and
disconnection jobs, which run to a 30,000-generation horizon for
censored replicates, account for the remainder.

## Sensitivity analysis

`code/sensitivity_scan.py` reruns the three headline endpoint contrasts
under one-factor-at-a-time perturbations: each ecological parameter
(`m_loc`, `m_glob`, `gamma`, `K`, `mu`, `w`) halved and doubled around
its baseline, evaluated on the fragmented (C0) and strongly mixed (C3)
architectures (fixation additionally at the intermediate S = 0.01; the
selection-response calibration at all four levels). The scan uses the
released deterministic seeding, so it is bit-reproducible: re-running
it regenerates `data/sensitivity/` exactly (about 30 minutes on a
single core).

```
python code/sensitivity_scan.py --outdir data/sensitivity
```

Findings (full tables in `data/sensitivity/README.md`): the direction
of every endpoint contrast is preserved in all scanned variants —
mixing raises the fixation probability of a deleterious loss at both
tested values of S, fragmentation delays complete core-gene loss in
both selective regimes, and selection is more attenuated under mixing
(F_eff(C0) > F_eff(C3) everywhere). At S = 0.02 the fixation contrast
is resolved at 95% in every variant, with C3/C0 ratios of 4.2–23.0 and
conservative ratio lower bounds above 1.75. The costly-regime
erosion-time contrast stays within 1.8–2.7, the engine-regime contrast
within 3.0–7.5, and the selection-attenuation factor
F_eff(C0)/F_eff(C3) within 1.5–4.0. The multiplicative factors quoted
in the captions are therefore baseline-specific magnitudes rather than
robust constants; only the directions and orders of magnitude are
invariant across the scanned parameter space.

## Reproducibility and provenance

- **Deterministic seeds.** Every job seeds NumPy's PCG64 generator from
  its own parameters through a SHA-256 based hash, so re-running the
  same job with the released code reproduces its output exactly on any
  machine.
- **Status of the archived data.** The files in `data/ibm_results/` are
  the canonical results used for the figures and the numbers quoted in
  the manuscript. They predate the seed change above (they were produced
  with Python's per-process randomised string hashing), so a full
  regeneration yields statistically equivalent but not bit-identical
  results. Spot checks with the released code: the calibration
  N_e,var agrees to within 1% (554 versus 553 for C0); fixation counts
  agree within Poisson error (for example 41 versus 46 fixed replicates
  at C0, S = 0, and 26 versus 30 at C3, S = 0.02); the time-averaged
  sectoring statistic Fhat agrees to within 0.005 in every regime; and
  in the disconnection experiment the connected-arm median reproduces
  to within 1% (4,675 versus 4,636), while the split time t* and the
  disconnected-arm median — statistics dominated by the slow, partly
  censored arm — moved by 10–20% in a regenerated run (t* = 622,
  disconnected median 8,059, against archived 751 and 9,228). The
  medians of the right-skewed, partly censored erosion-time
  distributions are the least stable statistic at R = 200:
  seed-to-seed standard deviations run from a few per cent to about 15%,
  and regenerated costly-regime medians for C0 span roughly
  8,800–11,400 generations against an archived 8,146. The numbers in
  the captions refer to the archived data, and `verify_numbers.py` pins
  them to `somatic_ibm_summary.json`.
- **Regenerated snapshot family.** The four `snap_*.npz` files — the
  C1 run added for the four-architecture panel b and regenerated C0,
  C2, C3 runs — were produced with the released deterministic seeding
  (unlike the families above), so re-running the snapshot jobs
  reproduces them byte-for-byte, and every panel-b statistic is exactly
  reproducible from `make_figures.py`. The regenerated grids are
  statistically equivalent to the earlier three-level archive; the
  panel-b medians (F = 1.00, 0.75, 0.65, 0.29 for C0–C3 at
  t = 1,800) agree with the independent calibration F_struct values
  (1.00, 0.76, 0.63, 0.30) within 0.03.
- **Restored `selresp`.** The selection-response calibration
  (`run_selresp`) was absent from the working script at archival time
  and has been restored for this release. The archived `sel_*.npz` were
  produced by an earlier implementation of the same measurement (an OLS
  slope of logit(x̄) on time over polymorphic checkpoints, converted to
  F_eff = −slope/(γS)) whose exact window and inclusion details were not
  preserved. The restored estimator regenerates F_eff values of 0.61,
  0.50, 0.46, 0.26 for C0–C3, against archived values of 0.55, 0.41,
  0.43, 0.22 (differences of up to 0.09; the cross-seed standard
  deviation of the restored estimator is about 0.03). The ordering and
  the C0/C3 contrast are preserved, but the near-tie between C1 and C2
  in the archive (0.41 versus 0.43) reverses in regenerated sets (0.50
  versus 0.46). The archived values remain the calibration used in the
  figures.
- **Batch configurations** in `build_jobs` reproduce the archived run
  configurations, including the disconnection trigger (median x̄ = 0.05)
  and the snapshot checkpoint times.

## Data

See `data/README.md` for a key-level dictionary of every `.npz` family
(producer command, keys, shapes and semantics).

## Figures

Both simulation figures (Fig. 1, Fig. S1) and the conceptual landscape
figure (Fig. 2) are delivered as editable vector files: SVG, with all
text preserved as editable text objects, and PDF, with embedded Type 42
(CID TrueType) fonts and text preserved as text. Fig. 2 is drawn at
exactly 7.2 × 4.05 in (183 mm, double column) without tight cropping, so
the physical size survives export. Raster PNG previews are regenerable
by the drawing scripts and are excluded from version control.
On-chart text follows an identity-only rule (which line, which event);
exact numbers and statistical detail live in the captions, whose LaTeX
sources accompany the figures. Panel letters are lowercase (a–e),
following the caption convention of the manuscript. The font stack is
Arial / Helvetica with a DejaVu Sans fallback on systems without Arial.

## Licence and citation

Code and data are released under the MIT License, copyright 2026
Ruizhe Xu (see `LICENSE`). If you use this code or the archived data,
please cite the repository:

```
Xu, R. (2026). Somatic compartment-network simulations of
module-integrity erosion (Version 1.0.0). Computer software.
MIT License.
```

Citation metadata is also provided in `CITATION.cff`.
