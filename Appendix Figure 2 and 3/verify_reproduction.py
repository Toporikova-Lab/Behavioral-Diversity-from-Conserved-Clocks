# =============================================================================
# verify_reproduction.py  --  did this machine reproduce the published numbers?
#
# Collects every number that appears in the response to Reviewer 1, in Supplement
# S1 and in Discussion 4.4, straight out of the files the pipeline just wrote, and
# compares them with reference_values.csv (the values from the run of 10 Sep 2026).
#
# The simulation is deterministic: every random draw uses a fixed seed, so a rerun
# on 500 individuals should match the reference to the last digit. Differences
# larger than the tolerances below mean something really is different, not noise.
#
#   python verify_reproduction.py                  compare against the reference
#   python verify_reproduction.py --make-reference  (re)write reference_values.csv
# =============================================================================

# %% imports and configuration
import numpy as np
import pandas as pd
import os
import sys

try:
    BASE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    BASE = os.getcwd()

SP = ["L. cornutus", "A. pennsylvanica", "S. grossa"]
PUB_SIM = {"L. cornutus": 0.09, "A. pennsylvanica": 0.06, "S. grossa": 0.75}   # beta in the code
I_SIM = 0.7                                    # lux 700 / 1000; manuscript uses I = 1
PUB_MS = {s: PUB_SIM[s] * I_SIM for s in SP}   # the same beta on the manuscript scale
CIRC = {"L. cornutus": "Larinioides_circular_stats_v2.csv",
        "A. pennsylvanica": "Agelenopsis_circular_stats_v2.csv",
        "S. grossa": "Steatoda_circular_stats_v2.csv"}
COMP = {"L. cornutus": "LC_spider_analysis_comprehensive_with_LD_split.csv",
        "A. pennsylvanica": "Ag_spider_analysis_comprehensive_with_LD_split.csv",
        "S. grossa": "Sg_spider_analysis_comprehensive_with_LD_split.csv"}
# absolute tolerance per kind of quantity; a deterministic rerun should give 0
TOL = {"percent": 0.1, "hours": 0.02, "param": 0.002, "ratio": 0.005}


# %% read a results file whether or not it is gzipped
def read_results(track):
    for name in ("results_individual.csv.gz", "results_individual.csv"):
        p = os.path.join(BASE, track, name)
        if os.path.exists(p):
            return pd.read_csv(p)
    raise FileNotFoundError(f"no results_individual.csv[.gz] in {track}")


# %% every checked quantity, computed from the files in this folder
def collect():
    rows = []

    def add(check, value, kind):
        rows.append(dict(check=check, value=value, kind=kind))

    # --- 1. observed values: data only, must match exactly on any machine
    for sp in SP:
        c = pd.read_csv(os.path.join(BASE, "spider_data", CIRC[sp]))
        sig = c.is_significant
        add(f"observed | {sp} | % significant rhythm in LD", 100 * sig.mean(), "percent")
        add(f"observed | {sp} | % in 22.8-25.2 h given a rhythm", 100 * c.is_entrained[sig].mean(), "percent")
        add(f"observed | {sp} | % entrained in LD", 100 * (sig & c.is_entrained).mean(), "percent")
        d = pd.read_csv(os.path.join(BASE, "spider_data", COMP[sp]))
        dd = d[(d.Condition == "DD") & (d.Mean_Activity > 0)]
        rh = dd[dd.Period_p_value < 0.05]
        add(f"observed | {sp} | % rhythmic in DD", 100 * len(rh) / len(dd), "percent")
        add(f"observed | {sp} | mean free-running period (h)", rh.Period_hours.mean(), "hours")
        add(f"observed | {sp} | SD of free-running period (h)", rh.Period_hours.std(ddof=1), "hours")

    # --- 2. calibration
    fit = pd.read_csv(os.path.join(BASE, "calibration_fit.csv"))
    for _, r in fit[fit.variant == "fitted"].iterrows():
        for col, kind in [("k_bar", "param"), ("sigma_scale", "param"), ("sigma_asym", "param")]:
            add(f"calibration | {r.species} | {col}", r[col], kind)
        add(f"calibration | {r.species} | fitted % rhythmic in DD", 100 * r.f_sim, "percent")

    # --- 3. main run, at the published beta
    res = read_results("track1_main")
    fitted = res[res.variant == "fitted"]
    for sp in SP:
        m = fitted[(fitted.species == sp) & np.isclose(fitted.beta, PUB_SIM[sp])]
        add(f"track1 fitted | {sp} | % rhythmic in DD", 100 * m.rhythmic_DD.mean(), "percent")
        add(f"track1 fitted | {sp} | % entrained (spider criterion)", 100 * m.entrained_LS.mean(), "percent")
        add(f"track1 fitted | {sp} | % phase-locked (exact)", 100 * m.entrained_true.mean(), "percent")
        add(f"track1 fitted | {sp} | mean light/dark activity ratio", m.light_dark_ratio.mean(), "ratio")

    # --- 4. the beta sweep, S. grossa, the curve the response quotes
    sg = fitted[fitted.species == "S. grossa"].groupby("beta").entrained_LS.mean() * 100
    for b in [0.0, 0.75, 2.0]:
        add(f"track1 fitted | S. grossa | % entrained at simulation beta {b}", sg.loc[b], "percent")
    add("track1 fitted | S. grossa | beta sweep is monotonic", float(np.all(np.diff(sg.values) >= -0.2)), "param")

    # --- 5. the supplement numbers, on the manuscript beta scale
    for sp in SP:
        g = fitted[fitted.species == sp].groupby("beta")
        x = g.entrained_LS.mean().index.values * I_SIM
        for col, lab in [("entrained_LS", "entrained (spider criterion)"), ("entrained_true", "phase-locked")]:
            y = getattr(g, col).mean().values * 100
            add(f"supplement | {sp} | % {lab} at manuscript beta {PUB_MS[sp]:.3f}",
                float(np.interp(PUB_MS[sp], x, y)), "percent")
        y = g.entrained_LS.mean().values * 100
        obs = 100 * (lambda c: (c.is_significant & c.is_entrained).mean())(
            pd.read_csv(os.path.join(BASE, "spider_data", CIRC[sp])))
        need = float(np.interp(obs, y, x)) if y.max() >= obs else np.nan
        add(f"supplement | {sp} | manuscript beta needed for the observed %", need, "param")

    return pd.DataFrame(rows)


# %% make the reference, or check against it
REF = os.path.join(BASE, "reference_values.csv")

if "--make-reference" in sys.argv:
    collect().to_csv(REF, index=False)
    print(f"wrote {REF} with {len(collect())} reference values")
else:
    got = collect()
    ref = pd.read_csv(REF)
    cmp = ref.merge(got, on="check", suffixes=("_reference", "_yours"), how="outer")
    cmp["tolerance"] = cmp.kind_reference.map(TOL)
    cmp["difference"] = cmp.value_yours - cmp.value_reference
    both_nan = cmp.value_yours.isna() & cmp.value_reference.isna()   # e.g. a beta the grid never reaches
    one_nan = (cmp.value_yours.isna() ^ cmp.value_reference.isna())
    cmp["result"] = np.where(both_nan, "match (both undefined)",
                     np.where(one_nan, "MISSING",
                      np.where(cmp.difference.abs() <= cmp.tolerance, "match", "DIFFERS")))
    cmp = cmp[["check", "value_reference", "value_yours", "difference", "tolerance", "result"]]
    cmp.to_csv(os.path.join(BASE, "verification_report.csv"), index=False)

    n_bad = (~cmp.result.str.startswith("match")).sum()
    exact = int(((cmp.difference.abs() < 1e-9) | both_nan).sum())
    pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200)
    print(cmp.round(4).to_string(index=False))
    print(f"\n{len(cmp)} checks: {len(cmp) - n_bad} match, {n_bad} do not. "
          f"{exact} are identical to the last digit.")
    if n_bad == 0:
        print("REPRODUCED. Every number in the response letter, Supplement S1 and Discussion 4.4\n"
              "was regenerated on this machine from the raw spider CSVs.")
    else:
        print("Some checks differ. Rows marked MISSING usually mean a stage has not been run yet;\n"
              "rows marked DIFFERS with a small difference usually mean N was reduced for a quick test.")
    print("\nwrote verification_report.csv")
