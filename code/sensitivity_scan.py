#!/usr/bin/env python3
"""One-factor-at-a-time (OFAT) sensitivity scan for the somatic IBM.

Repository entry point: produces data/sensitivity/ (120 jobs, about 30
minutes on a single core). The archived results in data/sensitivity/
were produced with this script and the released deterministic seeding,
so a re-run regenerates them exactly. Protocol and result tables are
documented in data/sensitivity/README.md.

Question: are the endpoint contrasts quoted in the figure captions robust
to halving and doubling the ecological parameters?

Experiments (all run with the released deterministic seeding, so they are
bit-reproducible; contrast with the pre-release archive in data/ibm_results):
  sel  selection-response calibration F_eff per connectivity level
       (S=0.05, R=400, T=400, Tw=60), all four levels C0-C3
  fix  single-mutant fixation probability, endpoints C0 vs C3,
       S in {0.01, 0.02}, R=100,000 (two blocks of 50,000), Tmax=40,000
  ero  time to complete core-gene loss, endpoints C0 vs C3, regimes
       costly (S=0.02, w=0) and engine (S=0, w=0.02), R=200,
       Tmax=30,000 (60,000 for K=16)

Variants (baseline values in parentheses):
  m_loc   0.04 / 0.16     (0.08; local fusion, affects C1-C3; C0 has no
                          fusion by definition and stays at m_loc = 0)
  m_glob  x0.5 / x2       (C2 0.03, C3 0.25; long-range anastomosis)
  gamma   0.125 / 0.5     (0.25; compartment turnover)
  K       4 / 16          (8; nuclei per compartment, N_tot 576 / 2304)
  mu      1.5e-4 / 6e-4   (3e-4; erosion experiment only)
  w       0.01 / 0.04     (0.02; engine regime only)

Outputs, all under --outdir:
  npz/<job>.npz         per-job raw statistics (counts, medians, F_eff)
  sensitivity_scan.csv  one row per job
  sensitivity_summary.json  headline robustness numbers

Usage (from the repository root):
  python code/sensitivity_scan.py --outdir data/sensitivity
  python code/sensitivity_scan.py --outdir data/sensitivity --only mloc
  python code/sensitivity_scan.py --outdir data/sensitivity --aggregate-only
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from somatic_ibm import generation, seed_of  # noqa: E402

D = 12                     # lattice side (compartments), as in the manuscript runs
LVLS = ["C0", "C1", "C2", "C3"]
BASE_M = {"C0": (0.00, 0.000), "C1": (0.08, 0.000),
          "C2": (0.08, 0.030), "C3": (0.08, 0.250)}
BASE = dict(K=8, gamma=0.25, mu=3e-4)

VARIANTS = {
    "baseline": {},
    "mloc_lo": {"m_loc": 0.04}, "mloc_hi": {"m_loc": 0.16},
    "mglob_lo": {"m_glob_scale": 0.5}, "mglob_hi": {"m_glob_scale": 2.0},
    "gamma_lo": {"gamma": 0.125}, "gamma_hi": {"gamma": 0.5},
    "K_lo": {"K": 4}, "K_hi": {"K": 16},
    "mu_lo": {"mu": 1.5e-4}, "mu_hi": {"mu": 6e-4},
    "w_lo": {"w": 0.01}, "w_hi": {"w": 0.04},
}
ARCH = ["baseline", "mloc_lo", "mloc_hi", "mglob_lo", "mglob_hi",
        "gamma_lo", "gamma_hi", "K_lo", "K_hi"]


def params(variant):
    """Resolve (m_loc, m_glob) per level + scalars for a variant."""
    v = VARIANTS[variant]
    K = v.get("K", BASE["K"])
    gam = v.get("gamma", BASE["gamma"])
    mu = v.get("mu", BASE["mu"])
    w = v.get("w", 0.02)
    m = {}
    for lv in LVLS:
        mloc, mglob = BASE_M[lv]
        # C0 is *defined* by the absence of fusion (m_loc = 0): the local-
        # fusion factor perturbs only architectures that have fusion.
        if mloc > 0:
            mloc = v.get("m_loc", mloc)
        if "m_glob_scale" in v and mglob > 0:
            mglob = mglob * v["m_glob_scale"]
        m[lv] = (mloc, mglob)
    return dict(K=K, gamma=gam, mu=mu, w=w, m=m)


# ---------------------------------------------------------------- runners
def run_fix(variant, level, S, R=100000, blocks=2, Tmax=40000):
    p = params(variant)
    K, gam = p["K"], p["gamma"]
    m_loc, m_glob = p["m"][level]
    N_tot = K * D * D
    fixed = lost = cens = 0
    Rb = R // blocks
    for b in range(blocks):
        rng = np.random.default_rng(
            seed_of("sensfix", variant, level, S, b))
        n = np.zeros((Rb, D, D), dtype=np.int16)
        ii = rng.integers(0, D, Rb)
        jj = rng.integers(0, D, Rb)
        n[np.arange(Rb), ii, jj] = 1
        alive = np.arange(Rb)
        fate = np.full(Rb, -1, np.int8)
        t = 0
        while len(alive) and t < Tmax:
            n = generation(n, K, m_loc, m_glob, gam, S, 0.0, 0.0, rng)
            tot = n.sum(axis=(1, 2))
            done = (tot == 0) | (tot == N_tot)
            if done.any():
                fate[alive[done]] = (tot[done] == N_tot).astype(np.int8)
            alive = alive[~done]
            n = n[~done]
            t += 1
        fixed += int((fate == 1).sum())
        lost += int((fate == 0).sum())
        cens += int((fate == -1).sum())
    return dict(experiment="fix", variant=variant, level=level, S=S, R=R,
                Tmax=Tmax, fixed=fixed, lost=lost, censored=cens)


def run_ero(variant, level, S, w, R=200):
    p = params(variant)
    K, gam, mu = p["K"], p["gamma"], p["mu"]
    m_loc, m_glob = p["m"][level]
    Tmax = 30000 if K == 8 else 60000
    N_tot = K * D * D
    rng = np.random.default_rng(seed_of("sensero", variant, level, S, w, R))
    n = np.zeros((R, D, D), dtype=np.int16)
    alive = np.arange(R)
    t_ero = np.full(R, -1, np.int32)
    t = 0
    while len(alive) and t < Tmax:
        n = generation(n, K, m_loc, m_glob, gam, S, w, mu, rng)
        done = n.sum(axis=(1, 2)) == N_tot
        if done.any():
            t_ero[alive[done]] = t + 1
        alive = alive[~done]
        n = n[~done]
        t += 1
    er = t_ero >= 0
    med = float(np.median(t_ero[er])) if er.any() else float("nan")
    return dict(experiment="ero", variant=variant, level=level, S=S, w=w,
                mu=mu, R=R, Tmax=Tmax, n_eroded=int(er.sum()),
                censored_pct=100.0 * (1.0 - er.mean()), median=med,
                t_ero=t_ero)


def run_sel(variant, level, S=0.05, R=400, T=400, Tw=60):
    p = params(variant)
    K, gam = p["K"], p["gamma"]
    m_loc, m_glob = p["m"][level]
    rng = np.random.default_rng(seed_of("senssel", variant, level, R))
    n = rng.binomial(K, 0.5, size=(R, D, D)).astype(np.int16)
    x = n.astype(np.float64) / K
    traj = np.empty((T + 1, R))
    traj[0] = x.mean(axis=(1, 2))
    for t in range(T):
        n = generation(n, K, m_loc, m_glob, gam, S, 0.0, 0.0, rng)
        traj[t + 1] = (n.astype(np.float64) / K).mean(axis=(1, 2))
    ts = np.arange(Tw, T + 1)
    slopes = []
    for r in range(R):
        xx = traj[ts, r]
        ok = (xx > 0.05) & (xx < 0.95)
        if ok.sum() >= 2:
            q = np.log(xx[ok] / (1.0 - xx[ok]))
            tt = ts[ok].astype(np.float64)
            tt = tt - tt.mean()
            slopes.append(float((tt * (q - q.mean())).sum() / (tt * tt).sum()))
    slopes = np.asarray(slopes)
    med = float(np.median(slopes))
    return dict(experiment="sel", variant=variant, level=level, S=S, R=R,
                T=T, Tw=Tw, K=K, gamma=gam, m_loc=m_loc, m_glob=m_glob,
                n_fit=int(slopes.size), slope=med,
                F_eff=float(-med / (gam * S)))


# ---------------------------------------------------------------- jobs
def build_jobs():
    jobs = []
    for v in ARCH:                                   # sel: 4 levels
        for lv in LVLS:
            jobs.append(("sel", v, dict(level=lv), "sel_%s_%s" % (v, lv)))
    for v in ARCH:                                   # fix: endpoints, 2 S
        for lv in ["C0", "C3"]:
            for S in [0.01, 0.02]:
                jobs.append(("fix", v, dict(level=lv, S=S),
                             "fix_%s_%s_S%g" % (v, lv, S)))
    ero_regimes = [("costly", 0.02, 0.0), ("engine", 0.0, 0.02)]
    for v in ARCH + ["mu_lo", "mu_hi"]:
        for reg, S, w in ero_regimes:
            for lv in ["C0", "C3"]:
                jobs.append(("ero", v, dict(level=lv, S=S, w=w),
                             "ero_%s_%s_%s" % (v, lv, reg)))
    for v in ["w_lo", "w_hi"]:
        for lv in ["C0", "C3"]:
            jobs.append(("ero", v, dict(level=lv, S=0.0, w=VARIANTS[v]["w"]),
                         "ero_%s_%s_engine" % (v, lv)))
    return jobs


def wilson(k, n, z=1.959963985):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    den = 1.0 + z * z / n
    ctr = (p + z * z / (2.0 * n)) / den
    half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return ctr - half, ctr + half


def npz_path(outdir, jid):
    return os.path.join(outdir, "npz", jid + ".npz")


def save(outdir, jid, res):
    path = npz_path(outdir, jid)
    arr = {k: v for k, v in res.items() if isinstance(v, np.ndarray)}
    sca = {k: v for k, v in res.items() if not isinstance(v, np.ndarray)}
    np.savez(path, **sca, **arr)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--only", default=None,
                    help="run only jobs whose id contains this substring")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--aggregate-only", action="store_true")
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, "npz"), exist_ok=True)

    if not args.aggregate_only:
        for kind, variant, kw, jid in build_jobs():
            if args.only and args.only not in jid:
                continue
            path = npz_path(args.outdir, jid)
            if os.path.exists(path) and not args.force:
                continue
            if kind == "sel":
                res = run_sel(variant, **kw)
                print("sel %-18s %-3s F_eff=%.3f (n=%d)"
                      % (variant, kw["level"], res["F_eff"], res["n_fit"]),
                      flush=True)
            elif kind == "fix":
                res = run_fix(variant, kw["level"], kw["S"])
                print("fix %-18s %-3s S=%.2g fixed=%d lost=%d cens=%d"
                      % (variant, kw["level"], kw["S"], res["fixed"],
                         res["lost"], res["censored"]), flush=True)
            else:
                res = run_ero(variant, kw["level"], kw["S"], kw["w"])
                print("ero %-18s %-3s S=%.2g w=%.2g median=%s eroded=%d/%d"
                      % (variant, kw["level"], kw["S"], kw["w"],
                         ("%.0f" % res["median"]) if res["median"] == res["median"] else "NA",
                         res["n_eroded"], res["R"]), flush=True)
            save(args.outdir, jid, res)

    aggregate(args.outdir)


def aggregate(outdir):
    rows, fixd, erod, seld = [], {}, {}, {}
    for path in sorted(glob.glob(os.path.join(outdir, "npz", "*.npz"))):
        z = dict(np.load(path, allow_pickle=False))
        jid = os.path.basename(path)[:-4]
        exp = jid.split("_")[0]
        if exp == "fix":
            k, n = int(z["fixed"]), int(z["R"])
            lo, hi = wilson(k, n)
            pfix = k / max(n - int(z["censored"]), 1)
            row = dict(experiment="fix", variant=z["variant"],
                       level=str(z["level"]), S=float(z["S"]), R=n,
                       fixed=k, lost=int(z["lost"]),
                       censored=int(z["censored"]), P_fix=pfix,
                       ci_lo=max(lo, 0.0), ci_hi=hi)
            fixd[(str(z["variant"]), str(z["level"]), float(z["S"]))] = row
        elif exp == "ero":
            row = dict(experiment="ero", variant=str(z["variant"]),
                       level=str(z["level"]), S=float(z["S"]), w=float(z["w"]),
                       mu=float(z["mu"]), R=int(z["R"]), Tmax=int(z["Tmax"]),
                       n_eroded=int(z["n_eroded"]),
                       censored_pct=float(z["censored_pct"]),
                       median=float(z["median"]))
            erod[(str(z["variant"]), str(z["level"]),
                  float(z["S"]), float(z["w"]))] = row
        else:
            row = dict(experiment="sel", variant=str(z["variant"]),
                       level=str(z["level"]), S=float(z["S"]), R=int(z["R"]),
                       n_fit=int(z["n_fit"]), slope=float(z["slope"]),
                       F_eff=float(z["F_eff"]))
            seld[(str(z["variant"]), str(z["level"]))] = row
        rows.append(row)

    cols = ["experiment", "variant", "level", "S", "w", "mu", "R", "Tmax",
            "fixed", "lost", "censored", "P_fix", "ci_lo", "ci_hi",
            "n_eroded", "censored_pct", "median", "n_fit", "slope", "F_eff"]
    with open(os.path.join(outdir, "sensitivity_scan.csv"), "w",
              newline="") as f:
        wtr = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        wtr.writeheader()
        for r in rows:
            wtr.writerow({k: r.get(k, "") for k in cols})

    # ---------------- headline robustness numbers
    summary = dict(protocol=dict(
        design="one-factor-at-a-time; each factor halved and doubled "
               "around the baseline",
        factors={k: v for k, v in VARIANTS.items()},
        fix=dict(levels=["C0", "C3"], S=[0.01, 0.02], R=100000,
                 blocks=2, Tmax=40000),
        ero=dict(levels=["C0", "C3"],
                 regimes={"costly": "S=0.02, w=0",
                          "engine": "S=0, w=0.02"},
                 R=200, Tmax="30000 (60000 for K=16)"),
        sel=dict(levels=LVLS, S=0.05, R=400, T=400, Tw=60)))

    fix_out = {}
    for S in [0.01, 0.02]:
        ratios, sign, undef, lo_ratio = {}, True, [], []
        for v in ARCH:
            c0 = fixd.get((v, "C0", S))
            c3 = fixd.get((v, "C3", S))
            if not (c0 and c3):
                continue
            if c3["fixed"] <= c0["fixed"]:
                sign = False
            if c0["fixed"] >= 1 and c3["fixed"] >= 1:
                ratios[v] = c3["P_fix"] / c0["P_fix"]
                # conservative ratio lower bound: C3 CI low / C0 CI high
                lo_ratio.append(c3["ci_lo"] / c0["ci_hi"])
            else:
                undef.append(v)
        fix_out["S=%g" % S] = dict(
            ratio={k: round(v, 2) for k, v in sorted(ratios.items())},
            ratio_min=round(min(ratios.values()), 2) if ratios else None,
            ratio_max=round(max(ratios.values()), 2) if ratios else None,
            ratio_ci_lower_min=round(min(lo_ratio), 2) if lo_ratio else None,
            sign_preserved_all_variants=sign,
            undefined_because_zero_events=undef)
    summary["fix"] = fix_out

    ero_out = {}
    w_eng = {"w_lo": VARIANTS["w_lo"]["w"], "w_hi": VARIANTS["w_hi"]["w"]}
    for reg, S, w in [("costly", 0.02, 0.0), ("engine", 0.0, 0.02)]:
        ratios = {}
        vs = ARCH + ["mu_lo", "mu_hi"]
        if reg == "engine":
            vs = vs + ["w_lo", "w_hi"]
        for v in vs:
            ww = w_eng.get(v, w)          # w variants use their own w
            c0 = erod.get((v, "C0", S, ww))
            c3 = erod.get((v, "C3", S, ww))
            if not (c0 and c3):
                continue
            if c0["median"] == c0["median"] and c3["median"] == c3["median"]:
                ratios[v] = c0["median"] / c3["median"]
        ero_out[reg] = dict(
            ratio={k: round(v, 2) for k, v in sorted(ratios.items())},
            ratio_min=round(min(ratios.values()), 2) if ratios else None,
            ratio_max=round(max(ratios.values()), 2) if ratios else None,
            sign_preserved_all_variants=bool(
                all(v > 1.0 for v in ratios.values())))
    summary["ero"] = ero_out

    att, order = {}, True
    for v in ARCH:
        f0 = seld.get((v, "C0"))
        f3 = seld.get((v, "C3"))
        if not (f0 and f3):
            continue
        if f0["F_eff"] <= f3["F_eff"]:
            order = False
        att[v] = f0["F_eff"] / f3["F_eff"]
    summary["sel"] = dict(
        F_eff_C0_over_C3={k: round(v, 2) for k, v in sorted(att.items())},
        attenuation_min=round(min(att.values()), 2) if att else None,
        attenuation_max=round(max(att.values()), 2) if att else None,
        C0_gt_C3_all_variants=order,
        F_eff_by_variant={v: {lv: round(seld[(v, lv)]["F_eff"], 3)
                              for lv in LVLS if (v, lv) in seld}
                          for v in ARCH})

    with open(os.path.join(outdir, "sensitivity_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print("aggregated %d rows -> sensitivity_scan.csv / "
          "sensitivity_summary.json" % len(rows))
    print(json.dumps({k: summary[k] for k in ["fix", "ero"]}, indent=1)[:1500])


if __name__ == "__main__":
    main()
