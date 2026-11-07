# Sensitivity analysis: `data/sensitivity/`

A one-factor-at-a-time (OFAT) robustness scan of the endpoint contrasts
quoted in the figure captions. Each ecological parameter was halved and
doubled around its baseline value; every other parameter was held fixed.
The scan asks whether the *direction* and order of magnitude of the
C0-versus-C3 contrasts survive, not whether the multiplicative factors
themselves are stable.

Producer: `python code/sensitivity_scan.py --outdir data/sensitivity`
(run from the repository root; about 30 minutes on a single core).
Unlike the archived results in `data/ibm_results/` (which predate the
seed change), the scan uses the released deterministic seeding and is
bit-reproducible: re-running it regenerates every `.npz` exactly.

## Protocol

| Factor | Low | Baseline | High | Applies to |
| --- | --- | --- | --- | --- |
| `m_loc` | 0.04 | 0.08 | 0.16 | architectures with local fusion (C1–C3; C0 is defined by m_loc = 0 and stays fragmented) |
| `m_glob` | ×0.5 | C2 0.03, C3 0.25 | ×2 | architectures with long-range fusion (C2–C3) |
| `gamma` | 0.125 | 0.25 | 0.5 | all levels |
| `K` | 4 | 8 | 16 | all levels (N_tot = 576 / 1152 / 2304) |
| `mu` | 1.5e-4 | 3e-4 | 6e-4 | erosion experiment |
| `w` | 0.01 | 0.02 | 0.04 | engine erosion regime |

Three experiments, endpoint architectures C0 and C3 throughout
(C1/C2 rank non-monotonically between regimes, so no intermediate
ordering is tested):

1. **fix** — single-mutant fixation probability (S = 0.01 and 0.02),
   R = 100,000 (two blocks of 50,000), Tmax = 40,000.
2. **ero** — median time to complete core-gene loss under the costly
   (S = 0.02, w = 0) and engine (S = 0, w = 0.02) regimes, R = 200,
   Tmax = 30,000 (60,000 for K = 16, where erosion is slower).
3. **sel** — selection-response calibration F_eff (S = 0.05, R = 400,
   T = 400, Tw = 60) at all four levels, as used for the diffusion
   overlay of panel c.

## Files

- `sensitivity_scan.csv` — one row per job (120 rows): counts, medians,
  censoring, P_fix with 95% Wilson intervals, slopes and F_eff.
- `sensitivity_summary.json` — per-variant ratios and headline ranges,
  as quoted in the repository README.
- `npz/<experiment>_<variant>_<level>[__<S/regime>].npz` — per-job raw
  statistics; `ero_*` files additionally store the per-replicate
  `t_ero` array.

## Results

**Fixation contrast C3/C0 (panel c).** Fixation of the same deleterious
loss is more likely under strong mixing than under fragmentation in
every variant and at both values of S:

| variant | fixed C0 (S=0.01) | fixed C3 (S=0.01) | ratio | fixed C0 (S=0.02) | fixed C3 (S=0.02) | ratio |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | 40 | 72 | 1.80 | 11 | 46 | 4.18 |
| m_loc 0.04 | 31 | 70 | 2.26 | 11 | 46 | 4.18 |
| m_loc 0.16 | 42 | 85 | 2.02 | 5 | 66 | 13.2 |
| m_glob ×0.5 | 40 | 77 | 1.92 | 9 | 43 | 4.78 |
| m_glob ×2 | 39 | 70 | 1.79 | 8 | 70 | 8.75 |
| gamma 0.125 | 35 | 86 | 2.46 | 13 | 60 | 4.62 |
| gamma 0.5 | 32 | 49 | 1.53 | 5 | 54 | 10.8 |
| K 4 | 76 | 158 | 2.08 | 22 | 105 | 4.77 |
| K 16 | 22 | 25 | 1.14 | 1 | 23 | 23.0 |

Counts are out of R = 100,000. At S = 0.02 the ratio lies between 4.2
and 23.0 across variants, and the conservative lower bound of every
ratio (C3 Wilson low divided by C0 Wilson high) exceeds 1.75, so the
contrast is individually resolved at 95% in every variant. At S = 0.01
the direction is preserved everywhere but the individual ratios
(1.1–2.5) are not all resolved at 95% at this replicate count. The
K = 16 / S = 0.02 cell rests on a single C0 fixation event and is
reported with that caveat. The scan baseline ratio (4.18) agrees with
the archived 5.0 (6 versus 30 events at R = 50,000) within Monte-Carlo
error.

**Erosion-time contrast C0/C3 (panel d).** Fragmentation always delays
complete core-gene loss:

| variant | costly C0 | costly C3 | ratio | engine C0 | engine C3 | ratio |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | 10885 (22.4% cens.) | 4368 (0.5%) | 2.49 | 2535 | 637.5 | 3.98 |
| m_loc 0.04 | 9801 (19.9%) | 4347 | 2.25 | 3302.5 | 629 | 5.25 |
| m_loc 0.16 | 9483 (23.0%) | 4379.5 | 2.17 | 3057 | 574.5 | 5.32 |
| m_glob ×0.5 | 10309 (20.4%) | 4560 (0.5%) | 2.26 | 2411.5 | 645.5 | 3.74 |
| m_glob ×2 | 10537 (17.0%) | 3840 | 2.74 | 2898.5 | 545.5 | 5.31 |
| gamma 0.125 | 7787 (9.0%) | 4408 | 1.77 | 3557.5 | 587.5 | 6.06 |
| gamma 0.5 | 11524.5 (31.9%) | 5567 (2.0%) | 2.07 | 2550 | 617.5 | 4.13 |
| K 4 | 11695 (3.5%) | 4275 | 2.74 | 2880 | 813 | 3.54 |
| K 16 | 13267 (3.0%) | 4862.5 | 2.73 | 3035.5 | 535.5 | 5.67 |
| mu 1.5e-4 | 14666 (45.4%) | 7660 (7.5%) | 1.91 | 5906 | 844 | 7.00 |
| mu 6e-4 | 5722 (1.5%) | 2416 | 2.37 | 1754.5 | 498 | 3.52 |
| w 0.01 | — | — | — | 2881.5 | 952.5 | 3.03 |
| w 0.04 | — | — | — | 2681.5 | 359.5 | 7.46 |

Medians are over completed replicates only (R = 200); percentages in
parentheses are replicates right-censored at the horizon, so heavily
censored cells (mu = 1.5e-4, costly regime) are lower bounds. The
costly-regime ratio stays within 1.77–2.74 and the engine-regime ratio
within 3.03–7.46 across all variants, against archived baseline values
of 2.08 and 4.91; the scan baseline medians themselves sit inside the
seed-to-seed span documented in the repository README (regenerated
costly C0 medians roughly 8,800–11,400 against an archived 8,146).

**Selection attenuation (panel c overlay).** The effective sectoring
coefficient is always larger under fragmentation, F_eff(C0) >
F_eff(C3), with C0/C3 between 1.47 and 3.95 (baseline 2.74):

| variant | F_eff C0 | F_eff C1 | F_eff C2 | F_eff C3 | C0/C3 |
| --- | --- | --- | --- | --- | --- |
| baseline | 0.594 | 0.486 | 0.434 | 0.217 | 2.74 |
| m_loc 0.04 | 0.583 | 0.493 | 0.523 | 0.257 | 2.26 |
| m_loc 0.16 | 0.608 | 0.416 | 0.334 | 0.191 | 3.18 |
| m_glob ×0.5 | 0.600 | 0.476 | 0.487 | 0.370 | 1.62 |
| m_glob ×2 | 0.593 | 0.441 | 0.382 | 0.156 | 3.79 |
| gamma 0.125 | 0.649 | 0.480 | 0.499 | 0.289 | 2.25 |
| gamma 0.5 | 0.566 | 0.403 | 0.380 | 0.278 | 2.04 |
| K 4 | 0.628 | 0.516 | 0.621 | 0.427 | 1.47 |
| K 16 | 0.582 | 0.353 | 0.328 | 0.147 | 3.95 |

C0 is insensitive to the migration factors by construction; the
attenuation is driven by the fused architectures. The C1/C2 ordering
flips in some variants (for example K = 4), consistent with their
near-tie in the archived calibration.

## Conclusions

Across every variant scanned (nine architecture factors for fixation
and attenuation; eleven for erosion; plus the engine-regime w factors):

- the direction of each endpoint contrast is preserved without
  exception: mixing raises the fixation probability of a deleterious
  loss, fragmentation delays complete core-gene loss in both regimes,
  and selection is more attenuated under mixing;
- at S = 0.02 the fixation contrast is resolved at 95% in every
  variant (conservative ratio lower bounds above 1.75);
- the multiplicative factors quoted in the captions (about five-fold
  for fixation at S = 0.02 and for the engine-regime erosion time,
  about two-fold for the costly regime) are baseline-specific: across
  variants the fixation ratio spans 4.2–23.0 (S = 0.02) and the
  engine-regime erosion ratio 3.0–7.5. The captions' wording
  ("roughly five-fold") should be read with that spread in mind.

The scan perturbs one factor at a time and leaves combinations
(for example simultaneous changes of migration and turnover) untested;
joint sweeps remain future work.
