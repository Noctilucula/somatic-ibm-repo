#!/usr/bin/env python3
"""
Two-level somatic evolution IBM on a mycelial compartment network.

Levels:  nuclei (within compartments, complementation + selfish selection w)
         compartments (lattice d x d, turnover gamma, between-selection S)
Nuclei migrate through hyphal fusion: local neighbour exchange (m_loc)
plus long-range anastomosis modelled as mixing with the global pool (m_glob).

State: n[r, i, j] = number of loss-of-function mutant nuclei in compartment
(i, j) of replicate r.  K nuclei per compartment, N_tot = K * d * d.

One generation (vectorised over replicates):
  1. within-compartment Wright-Fisher with selection w on mutants
  2. migration (pooled binomial resampling from the mixed pool)
  3. compartment turnover: a fraction gamma of compartments is replaced by a
     copy of a neighbour, donor drawn with weight W_d = 1 - S * x_d
  4. irreversible mutation of wild-type nuclei at rate mu

Subcommands (all resumable, npz outputs):
  calib      neutral variance-effective size N_e,var per connectivity level
  selresp    selection-response calibration: effective sectoring F_eff
  fixprob    single-mutant fixation probability vs S, per connectivity level
  erosion    recurrent mutation, time to complete core-gene loss
  disconnect paired-arm disconnection experiment
  snapshots  illustrative compartment grids at checkpoint times
"""
import argparse
import hashlib
import os
import time

import numpy as np

# ---------------------------------------------------------------- connectivity
# name: (m_loc, m_glob)
LEVELS = {
    "C0": (0.00, 0.000),   # fragmented: no nuclear migration through fusion
    "C1": (0.08, 0.000),   # local hyphal fusion only (stepping-stone)
    "C2": (0.08, 0.030),   # local + rare long-range anastomosis
    "C3": (0.08, 0.250),   # strong long-range fusion (small-world)
}

DEF = dict(d=12, K=8, gamma=0.25)   # compartment turnover rate per nuclear gen


def seed_of(*parts):
    """Deterministic, platform-independent seed from the run parameters.

    SHA-256 based (stable across processes and machines).  Note: the
    archived results in ibm_results/ were produced with an earlier
    per-process randomised string-hash seeding, so re-running a job
    regenerates statistically equivalent but not bit-identical output.
    """
    txt = "|".join(str(p) for p in parts)
    return int(hashlib.sha256(txt.encode("utf-8")).hexdigest()[:8], 16) \
        % (2**31 - 1)


# ---------------------------------------------------------------- one generation
def generation(n, K, m_loc, m_glob, gamma, S, w, mu, rng):
    """Advance one generation. Returns new state array (int16)."""
    R, d, _ = n.shape

    # 1. within-compartment WF with selfish advantage w for mutants
    nf = n.astype(np.float64)
    p = nf * (1.0 + w)
    p = p / (p + (K - nf))
    np.clip(p, 0.0, 1.0, out=p)
    n = rng.binomial(K, p).astype(np.int16)

    # 2. migration through the fused network (pooled resampling)
    m_tot = m_loc + m_glob
    if m_tot > 0:
        x = n.astype(np.float64) / K
        nb = 0.25 * (np.roll(x, 1, 1) + np.roll(x, -1, 1)
                     + np.roll(x, 1, 2) + np.roll(x, -1, 2))
        glo = x.mean(axis=(1, 2), keepdims=True)
        x_mix = (1.0 - m_tot) * x + m_loc * nb + m_glob * glo
        n = rng.binomial(K, x_mix).astype(np.int16)

    # 3. compartment turnover (local growth replaces dead compartments)
    if gamma > 0:
        mask = rng.random(n.shape) < gamma
        if mask.any():
            r_idx, i_idx, j_idx = np.nonzero(mask)
            m = len(r_idx)
            x = n.astype(np.float64) / K
            im, ip = (i_idx - 1) % d, (i_idx + 1) % d
            jm, jp = (j_idx - 1) % d, (j_idx + 1) % d
            xN = np.stack((x[r_idx, i_idx, jm], x[r_idx, i_idx, jp],
                           x[r_idx, im, j_idx], x[r_idx, ip, j_idx]))
            W = np.clip(1.0 - S * xN, 0.02, None)
            Wn = W / W.sum(axis=0)
            cum = np.cumsum(Wn, axis=0)
            don = (rng.random(m) < cum).argmax(axis=0)
            x_donor = xN[don, np.arange(m)]
            n[r_idx, i_idx, j_idx] = rng.binomial(K, x_donor)

    # 4. irreversible loss-of-function mutation
    if mu > 0:
        n += rng.binomial(K - n, mu).astype(np.int16)

    return n


def freq_stats(n, K):
    """Return xbar, F (per replicate). F = Var(x_c)/[xbar(1-xbar)]."""
    x = n.astype(np.float64) / K
    xb = x.mean(axis=(1, 2))
    var = x.var(axis=(1, 2))
    den = xb * (1.0 - xb)
    F = np.where(den > 1e-12, var / np.maximum(den, 1e-12), np.nan)
    return xb, F


# ---------------------------------------------------------------- subcommands
def run_calib(args):
    """Neutral runs from x0=0.5: variance-effective size N_e,var and the
    structural sectoring coefficient F (equilibrium, selection-weighted)."""
    d, K, gamma = DEF["d"], DEF["K"], DEF["gamma"]
    m_loc, m_glob = LEVELS[args.level]
    R, T = args.R, args.T
    rng = np.random.default_rng(seed_of("calib", args.level, R))
    n = rng.binomial(K, 0.5, size=(R, d, d)).astype(np.int16)
    traj = np.empty((T + 1, R))
    Fs = []
    traj[0] = freq_stats(n, K)[0]
    for t in range(T):
        n = generation(n, K, m_loc, m_glob, gamma, 0.0, 0.0, 0.0, rng)
        xb, F = freq_stats(n, K)
        traj[t + 1] = xb
        if t >= args.Tw:
            poly = (xb > 0.05) & (xb < 0.95)
            if poly.any():
                Fs.append(float(np.nanmean(F[poly])))
    inc = traj[1:args.Tw + 1] - traj[:args.Tw]
    var_inc = inc.var(axis=1)
    xbar = traj[:args.Tw].mean(axis=1)
    ne = xbar * (1 - xbar) / np.maximum(var_inc, 1e-15)
    out = dict(level=args.level, m_loc=m_loc, m_glob=m_glob,
               Ne_var=float(ne.mean()), Ne_var_sd=float(ne.std()),
               F_struct=float(np.mean(Fs)) if Fs else float("nan"),
               N_tot=K * d * d, M=d * d, K=K, gamma=gamma)
    np.savez(args.out, **out)
    print("calib %s Ne_var=%.1f (sd %.1f) F_struct=%.3f N_tot=%d"
          % (args.level, out["Ne_var"], out["Ne_var_sd"], out["F_struct"],
             out["N_tot"]))


def run_selresp(args):
    """Selection-response calibration at fixed S from x0 = 0.5.

    Per replicate, the OLS slope of logit(xbar) on t is estimated over
    polymorphic checkpoints (t >= Tw, 0.05 < xbar < 0.95); replicates
    with at least two such checkpoints contribute.  The reported slope is
    the median across replicates (robust against fast-loss trajectories),
    and slope_sd the robust scale (1.4826 x MAD).  The slope is converted
    into an effective sectoring coefficient
        F_eff = -slope / (gamma * S),
    i.e. the F that reproduces the measured response if between-
    compartment selection acts as s_eff = -gamma * F_eff * S.

    Provenance note: this function was restored for the code release;
    the archived sel_*.npz were produced by an earlier implementation of
    the same measurement whose window/inclusion details were not
    preserved.  Regenerated values agree with the archive to within
    Monte-Carlo error (see the repository README).
    """
    if args.S <= 0:
        raise SystemExit("selresp requires S > 0")
    d, K, gamma = DEF["d"], DEF["K"], DEF["gamma"]
    m_loc, m_glob = LEVELS[args.level]
    R, T = args.R, args.T
    rng = np.random.default_rng(seed_of("sel", args.level, args.S, R))
    n = rng.binomial(K, 0.5, size=(R, d, d)).astype(np.int16)
    traj = np.empty((T + 1, R))
    traj[0] = freq_stats(n, K)[0]
    for t in range(T):
        n = generation(n, K, m_loc, m_glob, gamma, args.S, 0.0, 0.0, rng)
        traj[t + 1] = freq_stats(n, K)[0]
    ts = np.arange(args.Tw, T + 1)
    slopes = []
    for r in range(R):
        x = traj[ts, r]
        ok = (x > 0.05) & (x < 0.95)
        if ok.sum() >= 2:
            q = np.log(x[ok] / (1.0 - x[ok]))
            tt = ts[ok].astype(np.float64)
            tt = tt - tt.mean()
            slopes.append(float((tt * (q - q.mean())).sum() / (tt * tt).sum()))
    slopes = np.asarray(slopes)
    med = float(np.median(slopes))
    mad = float(np.median(np.abs(slopes - med))) * 1.4826
    out = dict(level=args.level, S=args.S, R=R, T=T,
               slope=med, slope_sd=mad,
               n_fit=int(slopes.size),
               F_eff=float(-med / (gamma * args.S)),
               F_eff_sd=float(mad / (gamma * args.S)),
               N_tot=K * d * d, M=d * d, K=K, gamma=gamma)
    np.savez(args.out, **out)
    print("selresp %s S=%.4g -> slope=%.5f (mad %.5f, n=%d), F_eff=%.3f"
          % (args.level, args.S, out["slope"], out["slope_sd"],
             out["n_fit"], out["F_eff"]))


def run_fixprob(args):
    """Single mutant introduced; run to absorption; record fate + mean F."""
    d, K, gamma = DEF["d"], DEF["K"], DEF["gamma"]
    m_loc, m_glob = LEVELS[args.level]
    R, Tmax = args.R, args.Tmax
    rng = np.random.default_rng(seed_of("fix", args.level, args.S, R))
    n = np.zeros((R, d, d), dtype=np.int16)
    ii = rng.integers(0, d, R)
    jj = rng.integers(0, d, R)
    n[np.arange(R), ii, jj] = 1          # exactly one mutant per replicate
    N_tot = K * d * d
    alive = np.arange(R)
    fate = np.full(R, -1, np.int8)        # -1 censored, 0 lost, 1 fixed
    t_abs = np.full(R, -1, np.int32)
    Fsum = np.zeros(R); Fcnt = np.zeros(R)
    FWsum = np.zeros(R); Wsum = np.zeros(R)
    t = 0
    while len(alive) and t < Tmax:
        n = generation(n, K, m_loc, m_glob, gamma, args.S, args.w,
                       0.0, rng)
        tot = n.sum(axis=(1, 2))
        xb, F = freq_stats(n, K)
        good = np.isfinite(F)
        Fsum[alive[good]] += F[good]
        Fcnt[alive[good]] += 1
        wgt = xb * (1.0 - xb)
        goodw = good & (wgt > 1e-9)
        FWsum[alive[goodw]] += (F[goodw] * wgt[goodw])
        Wsum[alive[goodw]] += wgt[goodw]
        done = (tot == 0) | (tot == N_tot)
        if done.any():
            idx = alive[done]
            fate[idx] = (tot[done] == N_tot).astype(np.int8)
            t_abs[idx] = t + 1
        keep = ~done
        alive = alive[keep]
        n = n[keep]
        if t % 2000 == 0 and t > 0:
            print("  t=%d alive=%d" % (t, len(alive)), flush=True)
        t += 1
    absorbed = fate >= 0
    fixed = int((fate == 1).sum())
    lost = int((fate == 0).sum())
    cens = int((fate == -1).sum())
    Fhat = float(np.mean(Fsum[Fcnt > 0] / np.maximum(Fcnt[Fcnt > 0], 1))) \
        if (Fcnt > 0).any() else float("nan")
    Fsel = float(np.mean(FWsum[Wsum > 0] / np.maximum(Wsum[Wsum > 0], 1e-30))) \
        if (Wsum > 0).any() else float("nan")
    out = dict(level=args.level, S=args.S, w=args.w, R=R, Tmax=Tmax,
               fixed=fixed, lost=lost, censored=cens,
               P_fix=fixed / max(int(absorbed.sum()), 1),
               P_fix_all=fixed / R, Fhat=Fhat, Fsel=Fsel,
               t_fix_mean=float(t_abs[fate == 1].mean()) if fixed else np.nan,
               N_tot=N_tot, M=d * d, K=K, gamma=gamma)
    np.savez(args.out, **out)
    print("fixprob %s S=%.4g -> P_fix=%.6f (fixed %d / lost %d / cens %d), "
          "Fsel=%.3f" % (args.level, args.S, out["P_fix"], fixed, lost, cens,
                         Fsel))


def run_erosion(args):
    """Recurrent mutation; time to complete loss of the wild-type core gene."""
    d, K, gamma = DEF["d"], DEF["K"], DEF["gamma"]
    m_loc, m_glob = LEVELS[args.level]
    R, Tmax = args.R, args.Tmax
    rng = np.random.default_rng(seed_of("ero", args.level, args.S, args.w,
                                        args.mu, R))
    n = np.zeros((R, d, d), dtype=np.int16)
    N_tot = K * d * d
    alive = np.arange(R)
    t_ero = np.full(R, -1, np.int32)
    Fsum = np.zeros(R); Fcnt = np.zeros(R)
    FWsum = np.zeros(R); Wsum = np.zeros(R)
    t = 0
    while len(alive) and t < Tmax:
        n = generation(n, K, m_loc, m_glob, gamma, args.S, args.w,
                       args.mu, rng)
        tot = n.sum(axis=(1, 2))
        xb, F = freq_stats(n, K)
        good = np.isfinite(F) & (tot > 0) & (tot < N_tot)
        Fsum[alive[good]] += F[good]
        Fcnt[alive[good]] += 1
        wgt = xb * (1.0 - xb)
        goodw = good & (wgt > 1e-9)
        FWsum[alive[goodw]] += (F[goodw] * wgt[goodw])
        Wsum[alive[goodw]] += wgt[goodw]
        done = tot == N_tot
        if done.any():
            t_ero[alive[done]] = t + 1
        keep = ~done
        alive = alive[keep]
        n = n[keep]
        t += 1
    eroded = t_ero >= 0
    Fhat = float(np.mean(Fsum[Fcnt > 0] / np.maximum(Fcnt[Fcnt > 0], 1))) \
        if (Fcnt > 0).any() else float("nan")
    out = dict(level=args.level, S=args.S, w=args.w, mu=args.mu, R=R,
               Tmax=Tmax, n_eroded=int(eroded.sum()),
               t_ero=t_ero, Fhat=Fhat, N_tot=N_tot, M=d * d, K=K, gamma=gamma)
    np.savez(args.out, **out)
    med = float(np.median(t_ero[eroded])) if eroded.any() else float("nan")
    print("erosion %s S=%.3g w=%.3g mu=%.1g -> eroded %d/%d, median t=%s, "
          "Fhat=%.3f" % (args.level, args.S, args.w, args.mu,
                         eroded.sum(), R,
                         ("%.0f" % med) if med == med else "NA", Fhat))


def run_disconnect(args):
    """Paired erosion-rate experiment (E3-style discriminator).  Phase 1:
    erosion under C3 (recurrent mutation) until the median replicate crosses
    thr.  Arms: A stays connected (C3), B is disconnected (m = 0).  Both run
    to complete erosion; completion times and trajectories are recorded."""
    d, K, gamma = DEF["d"], DEF["K"], DEF["gamma"]
    R = args.R
    rng = np.random.default_rng(seed_of("disc", args.S, args.mu, R))
    n = np.zeros((R, d, d), dtype=np.int16)
    N_tot = K * d * d
    m3 = LEVELS["C3"]
    t_star = -1
    for t in range(args.T1max):
        n = generation(n, K, m3[0], m3[1], gamma, args.S, 0.0, args.mu, rng)
        xb = n.sum(axis=(1, 2)) / N_tot
        if np.median(xb) >= args.thr:
            t_star = t + 1
            break
    nA = n.copy()
    nB = n.copy()
    tA = np.full(R, -1, np.int32)
    tB = np.full(R, -1, np.int32)
    actA = np.arange(R)
    actB = np.arange(R)
    T2 = args.T2
    sub = max(1, T2 // args.nrec)
    rec_t, rec_xA, rec_xB, rec_FB = [], [], [], []
    for t2 in range(T2 + 1):
        if t2 % sub == 0:
            rec_t.append(t2)
            sA = nA.sum() + N_tot * (R - len(actA))
            sB = nB.sum() + N_tot * (R - len(actB))
            rec_xA.append(sA / (N_tot * R))
            rec_xB.append(sB / (N_tot * R))
            _, FB = freq_stats(nB, K)
            rec_FB.append(float(np.nanmean(FB)))
        if t2 < T2:
            if len(actA):
                nA = generation(nA, K, m3[0], m3[1], gamma, args.S, 0.0,
                                args.mu, rng)
                done = nA.sum(axis=(1, 2)) == N_tot
                if done.any():
                    tA[actA[done]] = t_star + t2 + 1
                    actA = actA[~done]
                    nA = nA[~done]
            if len(actB):
                nB = generation(nB, K, 0.0, 0.0, gamma, args.S, 0.0,
                                args.mu, rng)
                done = nB.sum(axis=(1, 2)) == N_tot
                if done.any():
                    tB[actB[done]] = t_star + t2 + 1
                    actB = actB[~done]
                    nB = nB[~done]
    out = dict(S=args.S, mu=args.mu, R=R, thr=args.thr, T2=T2, t_star=t_star,
               t_rec=np.array(rec_t), xA=np.array(rec_xA), xB=np.array(rec_xB),
               F_B=np.array(rec_FB), tA=tA, tB=tB,
               N_tot=N_tot, M=d * d, K=K, gamma=gamma)
    np.savez(args.out, **out)
    mA = np.median(tA[tA >= 0]) if (tA >= 0).any() else np.nan
    mB = np.median(tB[tB >= 0]) if (tB >= 0).any() else np.nan
    print("disconnect: t*=%d, erosion median A=%s (n=%d) B=%s (n=%d)"
          % (t_star, mA, (tA >= 0).sum(), mB, (tB >= 0).sum()))


def run_snapshots(args):
    """R replicate grids per level at checkpoint times (figure picks the
    median-state replicate per time)."""
    d, K, gamma = DEF["d"], DEF["K"], DEF["gamma"]
    m_loc, m_glob = LEVELS[args.level]
    R = args.R
    rng = np.random.default_rng(seed_of("snap", args.level, args.S, args.mu,
                                        R))
    n = np.zeros((R, d, d), dtype=np.int16)
    times = sorted(set([0] + list(args.times)))
    grids = {}
    N_tot = K * d * d
    for t in range(max(times) + 1):
        if t in times:
            grids[t] = n.astype(np.int16).copy()
        n = generation(n, K, m_loc, m_glob, gamma, args.S, 0.0, args.mu, rng)
    out = dict(level=args.level, S=args.S, mu=args.mu, R=R,
               times=np.array(times), N_tot=N_tot)
    for t in times:
        out["g%d" % t] = grids[t]
    np.savez_compressed(args.out, **out)
    means = [float(grids[t].sum() / (N_tot * R)) for t in times]
    print("snapshots %s: mean x per time %s" % (args.level,
          ["%.3f" % m for m in means]))


# ---------------------------------------------------------------- batch driver
def build_jobs(outdir):
    jobs = []
    levels = ["C0", "C1", "C2", "C3"]
    for lv in levels:
        jobs.append(("calib", dict(level=lv, R=2000, T=300, Tw=60),
                     os.path.join(outdir, "calib_%s.npz" % lv)))
        jobs.append(("selresp", dict(level=lv, S=0.05, R=400, T=400, Tw=60),
                     os.path.join(outdir, "sel_%s.npz" % lv)))
    Sgrid = [0.0, 0.01, 0.02, 0.05, 0.1]
    for lv in levels:
        for S in Sgrid:
            jobs.append(("fixprob", dict(level=lv, S=S, w=0.0, R=50000,
                                         Tmax=40000),
                         os.path.join(outdir, "fix_%s_S%.4g.npz" % (lv, S))))
    regimes = [dict(S=0.02, w=0.0), dict(S=0.0, w=0.0), dict(S=0.0, w=0.02)]
    for reg in regimes:
        for lv in levels:
            jobs.append(("erosion", dict(level=lv, S=reg["S"], w=reg["w"],
                                          mu=3e-4, R=200, Tmax=30000),
                         os.path.join(outdir, "ero_%s_S%.3g_w%.3g.npz"
                                      % (lv, reg["S"], reg["w"]))))
    jobs.append(("disconnect", dict(S=0.02, mu=3e-4, R=200, thr=0.05,
                                    T1max=30000, T2=28000, nrec=280),
                 os.path.join(outdir, "disc_S0.02.npz")))
    for lv in levels:
        jobs.append(("snapshots", dict(level=lv, S=0.004, mu=3e-4, R=40,
                                       times=[0, 1200, 1800, 2400, 3000,
                                              3600, 4200, 4800, 6000]),
                     os.path.join(outdir, "snap_%s.npz" % lv)))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["calib", "fixprob", "erosion",
                                    "selresp", "disconnect", "snapshots",
                                    "batch", "smoke"])
    ap.add_argument("--level", default="C0")
    ap.add_argument("--S", type=float, default=0.0)
    ap.add_argument("--w", type=float, default=0.0)
    ap.add_argument("--mu", type=float, default=0.0)
    ap.add_argument("--R", type=int, default=1000)
    ap.add_argument("--T", type=int, default=80)
    ap.add_argument("--Tw", type=int, default=60)
    ap.add_argument("--Tmax", type=int, default=40000)
    ap.add_argument("--thr", type=float, default=0.4)
    ap.add_argument("--T1max", type=int, default=30000)
    ap.add_argument("--T2", type=int, default=12000)
    ap.add_argument("--nrec", type=int, default=240)
    ap.add_argument("--times", type=int, nargs="+",
                    default=[100, 800, 2500, 5000])
    ap.add_argument("--out", default="/tmp/ibm_out.npz")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--max-wall", type=float, default=420.0)
    args = ap.parse_args()

    if args.cmd == "smoke":
        args.level, args.S, args.mu, args.R, args.Tmax = \
            "C2", 0.01, 0.0, 300, 5000
        args.out = "/tmp/smoke_fix.npz"
        run_fixprob(args)
        args.level, args.S, args.mu, args.R, args.Tmax = \
            "C3", 0.0, 3e-4, 100, 4000
        args.out = "/tmp/smoke_ero.npz"
        run_erosion(args)
        return

    if args.cmd == "calib":
        run_calib(args)
    elif args.cmd == "selresp":
        run_selresp(args)
    elif args.cmd == "fixprob":
        run_fixprob(args)
    elif args.cmd == "erosion":
        run_erosion(args)
    elif args.cmd == "disconnect":
        run_disconnect(args)
    elif args.cmd == "snapshots":
        run_snapshots(args)
    elif args.cmd == "batch":
        outdir = args.outdir or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "ibm_results")
        os.makedirs(outdir, exist_ok=True)
        jobs = build_jobs(outdir)
        t0 = time.time()
        done = 0
        for kind, kw, path in jobs:
            if os.path.exists(path):
                done += 1
                continue
            if time.time() - t0 > args.max_wall:
                print("wall limit reached; %d/%d jobs done, resumable"
                      % (done, len(jobs)))
                return
            a = argparse.Namespace(cmd=kind, **dict(kw, out=path,
                                                     max_wall=args.max_wall))
            print("[job %s %s]" % (kind, path), flush=True)
            if kind == "calib":
                run_calib(a)
            elif kind == "selresp":
                run_selresp(a)
            elif kind == "fixprob":
                run_fixprob(a)
            elif kind == "erosion":
                run_erosion(a)
            elif kind == "disconnect":
                run_disconnect(a)
            elif kind == "snapshots":
                run_snapshots(a)
            done += 1
        print("batch complete: %d/%d jobs" % (done, len(jobs)))


if __name__ == "__main__":
    main()
