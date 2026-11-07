# Data dictionary: `data/ibm_results/`

All files are NumPy `.npz` archives written by `code/somatic_ibm.py`.
Throughout, the state variable is n[r, i, j] = number of mutant
(loss-of-function core-gene) nuclei in compartment (i, j) of replicate r;
x = n/K is the mutant fraction of a compartment; x̄ is the whole-body
mutant frequency; and F = Var(x_c) / [x̄(1 − x̄)] is the sectoring
coefficient across the 144 compartments.

Common scalar keys in every family: `level` (connectivity architecture,
where applicable), `N_tot` = 1152 (census number of nuclei), `M` = 144
(compartments), `K` = 8 (nuclei per compartment), `gamma` = 0.25
(compartment turnover per nuclear generation).

Connectivity architectures (migration parameters in `code/somatic_ibm.py`,
`LEVELS`):

| Level | m_loc | m_glob | Interpretation |
| --- | --- | --- | --- |
| C0 | 0.00 | 0.000 | fragmented, no nuclear migration through fusion |
| C1 | 0.08 | 0.000 | local hyphal fusion only (stepping-stone) |
| C2 | 0.08 | 0.030 | local fusion plus rare long-range anastomosis |
| C3 | 0.08 | 0.250 | strong long-range fusion (small-world) |

## `calib_<level>.npz` — neutral calibration (4 files)

Producer: `python code/somatic_ibm.py calib --level <lv> --R 2000 --T 300 --Tw 60`

Neutral runs from x0 = 0.5. Keys: `m_loc`, `m_glob`; `Ne_var`,
`Ne_var_sd` — variance-effective population size from the one-generation
increments of x̄ (theory: x̄(1 − x̄)/Var(Δx̄)); `F_struct` — mean sectoring
coefficient at quasi-equilibrium, averaged over polymorphic replicates
(0.05 < x̄ < 0.95) from generation Tw onwards.

## `sel_<level>.npz` — selection-response calibration (4 files)

Producer: `python code/somatic_ibm.py selresp --level <lv> --S 0.05 --R 400 --T 400 --Tw 60`

Selection response at fixed S = 0.05 from x0 = 0.5. Keys: `slope` —
median across replicates of the OLS slope of logit(x̄) on t over
polymorphic checkpoints (t ≥ Tw, 0.05 < x̄ < 0.95); `slope_sd` — robust
scale (1.4826 × MAD) of the per-replicate slopes; `n_fit` — number of
contributing replicates; `F_eff` = −slope/(γS) and `F_eff_sd` — the
effective sectoring coefficient, i.e. the F that reproduces the measured
response if between-compartment selection acts as s_eff = −γ·F_eff·S.
Provenance: the estimator was restored for the release; see the
repository README ("Restored `selresp`").

## `fix_<level>_S<x>.npz` — fixation probability (20 files)

Producer: `python code/somatic_ibm.py fixprob --level <lv> --S <x> --w 0 --R 50000 --Tmax 40000`

A single mutant nucleus introduced at a uniformly random position; runs
to absorption. Keys: `fixed`, `lost`, `censored` — fate counts over R
replicates; `P_fix` — fixed/absorbed (conditional on absorption within
Tmax); `P_fix_all` — fixed/R; `Fhat` — unweighted mean F over polymorphic
checkpoints; `Fsel` — x̄(1 − x̄)-weighted mean F; `t_fix_mean` — mean
generation of fixation (fixed replicates only).

## `ero_<level>_S<x>_w<y>.npz` — erosion time (12 files)

Producer: `python code/somatic_ibm.py erosion --level <lv> --S <x> --w <y> --mu 3e-4 --R 200 --Tmax 30000`

Recurrent mutation (μ = 3 × 10⁻⁴; note that μ is stored in the file but
not encoded in the file name) from an all-wild-type start. Keys:
`n_eroded` — replicates reaching complete core-gene loss within Tmax;
`t_ero` — int32 array, generation of complete loss per replicate
(−1 = censored at Tmax); `Fhat` — mean F over polymorphic checkpoints.
Regimes archived: (S, w) = (0.02, 0) costly core; (0, 0) short-circuited
selection; (0, 0.02) selfish saving.

## `disc_S0.02.npz` — paired-arm disconnection experiment (1 file)

Producer: `python code/somatic_ibm.py disconnect --S 0.02 --mu 3e-4 --R 200 --thr 0.05 --T1max 30000 --T2 28000 --nrec 280`

Phase 1 erodes the tissue under strong mixing (C3) until the median
replicate crosses x̄ = `thr`; the split occurs at generation `t_star`.
Arm A then remains connected (C3), arm B is disconnected (m = 0); both
run to complete erosion or the horizon T2. Keys: `thr`, `t_star`, `T2`;
`t_rec`, `xA`, `xB` — 281 checkpoints of the ensemble-mean mutant
frequency of arms A and B (times relative to t_star); `F_B` — sectoring
coefficient of arm B at the same checkpoints; `tA`, `tB` — int32 arrays,
absolute generation of complete erosion per replicate (−1 = censored).

## `snap_<level>.npz` — snapshot grids (4 files)

Producer: `python code/somatic_ibm.py snapshots --level <lv> --S 0.004 --mu 3e-4 --R 40 --times 0 1200 1800 2400 3000 3600 4200 4800 6000`

Keys: `times` — checkpoint generations; `g<t>` — int16 array of shape
(R, 12, 12): mutant-nucleus counts per compartment at generation t
(0 ≤ count ≤ K). Levels archived: C0, C1, C2, C3. This family was
generated with the released deterministic seeding (see the repository
README, "Regenerated snapshot family"): the C1 run was added when the
main figure's panel b was redesigned to cover all four architectures,
and the C0–C3 runs regenerate byte-for-byte with the released code.

## `somatic_ibm_summary.json` (in `figures/`)

Derived summary of every number quoted in the figure captions, written
by `code/make_figures.py`; validated against the captions by
`code/verify_numbers.py`.
