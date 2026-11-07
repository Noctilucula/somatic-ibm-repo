#!/usr/bin/env python3
"""Cross-check every number quoted in the two caption .tex files against
somatic_ibm_summary.json (the numbers produced by make_figures.py from
the archived data in data/ibm_results), and pin the headline sensitivity
numbers quoted in README.md to data/sensitivity/sensitivity_summary.json.
Exits non-zero on any mismatch; prints "PASS n checks" on success."""
import json
import os
import sys

D = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 os.pardir, "figures")
S = json.load(open(os.path.join(D, "somatic_ibm_summary.json")))
cap = open(os.path.join(D, "Fig_SomaticIBM_main_caption.tex"),
           encoding="utf-8").read()
capS = open(os.path.join(D, "Fig_SomaticIBM_snapshots_caption.tex"),
            encoding="utf-8").read()
capL = open(os.path.join(D, "Fig_LandscapeContinuum_caption.tex"),
            encoding="utf-8").read()

lv = S["levels"]
ok, bad = [], []


def chk(name, cond):
    (ok if cond else bad).append(name)


# panel a / c: F_eff and Ne_var rounded
chk("F_eff 0.55/0.41/0.43/0.22",
    [round(lv[k]["F_eff"], 2) for k in ["C0", "C1", "C2", "C3"]]
    == [0.55, 0.41, 0.43, 0.22])
chk("Ne_var 553/550/508/480",
    [round(lv[k]["Ne_var"]) for k in ["C0", "C1", "C2", "C3"]]
    == [553, 550, 508, 480])
chk("caption quotes F_eff list", "0.55" in cap and "553" in cap)

# panel b: common-checkpoint structure statistics
LVL4 = ["C0", "C1", "C2", "C3"]
pb = S["panelB"]
chk("panel-b F medians 1.00/0.75/0.65/0.29",
    [round(pb[k]["F_med"], 2) for k in LVL4] == [1.00, 0.75, 0.65, 0.29])
chk("panel-b common checkpoint t=1800, n_poly 19/22/21/20",
    all(pb[k]["t"] == 1800 for k in LVL4)
    and [pb[k]["n_poly"] for k in LVL4] == [19, 22, 21, 20])
chk("panel-b C0 replicate all-or-nothing, C3 intermediate",
    abs(pb["C0"]["rep_frac0"] + pb["C0"]["rep_frac8"] - 0.993) < 0.004
    and abs(pb["C3"]["rep_frac_int"] - 0.889) < 0.004)
chk("panel-b F medians monotone and near F_struct",
    all(pb[a]["F_med"] > pb[b_]["F_med"]
        for a, b_ in zip(LVL4, LVL4[1:]))
    and all(abs(pb[k]["F_med"] - lv[k]["F_struct"]) < 0.03 for k in LVL4))
chk("caption quotes panel-b medians, ns and Supplementary Fig. 1",
    "0.29" in cap and "n=19" in cap and "Supplementary Fig.~1" in cap)

# panel c: five-fold at S=0.02, detection limit, neutral expectation
r = S["P_fix"]["C3|S=0.02"]["p"] / S["P_fix"]["C0|S=0.02"]["p"]
chk("C3/C0 fixation ratio = 5.0", abs(r - 5.0) < 0.05)
chk("detection limit 2/R = 4e-5",
    abs(2.0 / S["P_fix"]["C0|S=0"]["R"] - 4e-5) < 1e-12)
chk("neutral 1/N_tot ~ 8.7e-4",
    abs(1.0 / S["N_tot"] - 8.7e-4) < 3e-5)

# panel d: erosion medians and censoring
e = S["erosion"]
costly = [int(e["costly|%s" % k]["median"] + 0.5)
          for k in ["C0", "C1", "C2", "C3"]]
chk("costly 8146/6756/7470/3919", costly == [8146, 6756, 7470, 3919])
neut = [int(e["neutral|%s" % k]["median"] + 0.5)
        for k in ["C0", "C1", "C2", "C3"]]
chk("neutral 2855/3436/3415/3347", neut == [2855, 3436, 3415, 3347])
eng = [int(e["engine|%s" % k]["median"] + 0.5)
      for k in ["C0", "C1", "C2", "C3"]]
chk("engine 2905/1293/971/591", eng == [2905, 1293, 971, 591])
cen = [round(100.0 * (e["costly|%s" % k]["R"] - e["costly|%s" % k]["eroded"])
              / e["costly|%s" % k]["R"], 6)
       for k in ["C0", "C1", "C2", "C3"]]
chk("censored 22.5/6.5/2.5/0 %", cen == [22.5, 6.5, 2.5, 0.0])
chk("engine ratio ~5x C0->C3", abs(eng[0] / eng[3] - 4.91) < 0.05)

# panel e: disconnection
d = S["disconnect"]
chk("t*=751, thr=0.05, mean at t* ~0.16",
    d["t_star"] == 751 and abs(d["thr"] - 0.05) < 1e-12
    and abs(d["xbar_at_tstar"] - 0.16) < 0.01)
chk("absolute medians A=4636 (198/200), B=9228 (163/200)",
    (round(d["median_A"]), round(d["median_B"]), d["nA"], d["nB"])
    == (4636, 9228, 198, 163))
chk("post-t* medians 3885/8477",
    (round(d["median_A_post"]), round(d["median_B_post"])) == (3885, 8477))
chk("post-t* ratio > 2", d["median_B_post"] / d["median_A_post"] > 2.0)
chk("caption quotes 0.05 trigger and both time bases",
    "0.05" in cap and "3{,}885" in cap and "8{,}477" in cap)
chk("neutral clock 1/mu = 3333", round(1 / S["mu"]) == 3333)

# S1 caption
chk("S1: R=40, times list, S=0.004",
    "R=40" in capS and "1{,}200" in capS and "6{,}000" in capS
    and "0.004" in capS)

# final numbering; methodological statements live in the repository
chk("main caption numbered Fig.~1", "Fig.~1\\,|" in cap)
chk("no stale Fig.~2 reference in main caption", "Fig.~2" not in cap)
chk("S1 caption cross-references Fig.~1a", "Fig.~1a" in capS)
chk("no stale Fig.~2 reference in S1 caption", "Fig.~2" not in capS)
chk("no [TO-METHODS] marker left", "[TO-METHODS]" not in cap)
chk("update-rule mechanics moved out of the legend", "replenished" not in cap)
chk("scope statement moved out of the legend",
    "conceptual demonstrations" not in cap)

# sensitivity scan: headline numbers quoted in the README
SEN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   os.pardir, "data", "sensitivity",
                   "sensitivity_summary.json")
sen = json.load(open(SEN))
chk("sens: fixation contrast direction preserved (both S)",
    sen["fix"]["S=0.01"]["sign_preserved_all_variants"]
    and sen["fix"]["S=0.02"]["sign_preserved_all_variants"])
chk("sens: fix S=0.02 ratio range 4.18-23.0 (README: 4.2-23.0)",
    sen["fix"]["S=0.02"]["ratio_min"] == 4.18
    and sen["fix"]["S=0.02"]["ratio_max"] == 23.0)
chk("sens: fix S=0.02 conservative CI lower bounds >= 1.75",
    sen["fix"]["S=0.02"]["ratio_ci_lower_min"] == 1.75)
chk("sens: erosion contrast direction preserved (both regimes)",
    sen["ero"]["costly"]["sign_preserved_all_variants"]
    and sen["ero"]["engine"]["sign_preserved_all_variants"])
chk("sens: costly ratio range 1.77-2.74 (README: 1.8-2.7)",
    sen["ero"]["costly"]["ratio_min"] == 1.77
    and sen["ero"]["costly"]["ratio_max"] == 2.74)
chk("sens: engine ratio range 3.03-7.46 (README: 3.0-7.5)",
    sen["ero"]["engine"]["ratio_min"] == 3.03
    and sen["ero"]["engine"]["ratio_max"] == 7.46)
chk("sens: F_eff(C0) > F_eff(C3) in every variant",
    sen["sel"]["C0_gt_C3_all_variants"])
chk("sens: attenuation range 1.47-3.95 (README: 1.5-4.0)",
    sen["sel"]["attenuation_min"] == 1.47
    and sen["sel"]["attenuation_max"] == 3.95)
chk("sens: 120-job scan CSV archived",
    sum(1 for _ in open(SEN.replace("sensitivity_summary.json",
                                    "sensitivity_scan.csv"))) == 121)

# Fig. 2 landscape caption (pure schematic: no numbers to cross-check)
flatL = " ".join(capL.split())
chk("Fig2 caption numbered Fig.~2", "Fig.~2\\,|" in capL)
chk("Fig2 caption has panels a and b",
    "\\textbf{a}," in capL and "\\textbf{b}," in capL)
chk("Fig2 caption ring wording (single ring, basin width)",
    "a dashed ring whose radius follows the basin width" in flatL)
chk("Fig2 caption describes all four stages",
    all(k in flatL for k in ["curved edges", "straighten and vanish",
                            "red crossed dotted line",
                            "only the single deep basin remains"]))
chk("Fig2 caption British spelling",
    "coloured" in capL and "colored" not in capL)
for f in ["Fig_LandscapeContinuum.svg", "Fig_LandscapeContinuum.pdf"]:
    p = os.path.join(D, f)
    chk("Fig2 figure present & non-trivial: %s" % f,
        os.path.exists(p) and os.path.getsize(p) > 10000)

# placeholders / leftovers
for name, txt in [("main caption", cap), ("S1 caption", capS),
                  ("Fig2 caption", capL)]:
    for p in ["TODO", "XXX", "【", "____", "placeholder"]:
        chk("%s free of %r" % (name, p), p.lower() not in txt.lower())

print("PASS %d checks" % len(ok))
if bad:
    print("FAILED:", bad)
    sys.exit(1)
