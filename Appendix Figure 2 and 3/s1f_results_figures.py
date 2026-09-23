# =============================================================================
# s1f_results_figures.py   --   Simulation S1, step F: comparisons and figures
#
# Reads results_individual.csv from the folder it sits in and produces the
# comparisons listed in the plan: within species (DD -> LD, paired), between
# species at published beta (vs observed), and across the beta grid (paired,
# same individuals at every beta). Observed values are recomputed from the
# spider CSVs in ../spider_data every time.
# =============================================================================

# %% imports and inputs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os

try:
    OUT = os.path.dirname(os.path.abspath(__file__))
except NameError:
    OUT = os.getcwd()
DATA = os.path.join(OUT, "spider_data") if os.path.isdir(os.path.join(OUT, "spider_data")) \
    else os.path.join(os.path.dirname(OUT), "spider_data")

res = pd.read_csv(os.path.join(OUT, "results_individual.csv.gz" if os.path.exists(os.path.join(OUT, "results_individual.csv.gz")) else "results_individual.csv"))
VARIANT = os.environ.get("S1_PLOT_VARIANT", res.variant.iloc[0])
res = res[res.variant == VARIANT].copy()
PUB = {"L. cornutus": 0.09, "A. pennsylvanica": 0.06, "S. grossa": 0.75}
COL = {"L. cornutus": "#1f77b4", "A. pennsylvanica": "#ff7f0e", "S. grossa": "#2ca02c"}
CIRC = {"L. cornutus": "Larinioides_circular_stats_v2.csv", "A. pennsylvanica": "Agelenopsis_circular_stats_v2.csv",
        "S. grossa": "Steatoda_circular_stats_v2.csv"}
COMP = {"L. cornutus": "LC_spider_analysis_comprehensive_with_LD_split.csv",
        "A. pennsylvanica": "Ag_spider_analysis_comprehensive_with_LD_split.csv",
        "S. grossa": "Sg_spider_analysis_comprehensive_with_LD_split.csv"}
res["published"] = res.apply(lambda r: np.isclose(r.beta, PUB[r.species]), axis=1)
pub = res[res.published]


# %% observed values, recomputed
def observed(sp):
    c = pd.read_csv(os.path.join(DATA, CIRC[sp]))
    d = pd.read_csv(os.path.join(DATA, COMP[sp]))
    ent = c.is_significant & c.is_entrained
    dd = d[(d.Condition == "DD") & (d.Mean_Activity > 0)]
    ldl = d[d.Condition == "LD_L"].Mean_Activity.mean(); ldd = d[d.Condition == "LD_D"].Mean_Activity.mean()
    ang = 2 * np.pi * c.loc[ent, "mean_phase_ZT"] / 24
    return dict(pct_entrained=100 * ent.mean(), n_ld=len(c),
                tau_LD_sd=c.loc[ent, "period_hr"].std(ddof=1),
                phase_ZT=(np.degrees(np.arctan2(np.sin(ang).sum(), np.cos(ang).sum())) % 360) / 15,
                vector_r=c.vector_strength.mean(),
                pct_rhythmic_DD=100 * (dd.Period_p_value < 0.05).mean(), n_dd=len(dd),
                light_dark=ldl / ldd)

obs = pd.DataFrame({sp: observed(sp) for sp in PUB}).T
print("observed (recomputed from spider_data):")
print(obs.round(3).to_string())


def wilson(k, n):
    lo, hi = stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0, stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1
    return 100 * lo, 100 * hi


# %% table 1: species at published beta, model vs observed
rows = []
for sp in PUB:
    m = pub[pub.species == sp]
    e = m[m.entrained_LS]
    rows.append(dict(species=sp, n=len(m),
                     model_pct_rhythmic_DD=100 * m.rhythmic_DD.mean(), obs_pct_rhythmic_DD=obs.loc[sp, "pct_rhythmic_DD"],
                     model_pct_entrained_LS=100 * m.entrained_LS.mean(), obs_pct_entrained=obs.loc[sp, "pct_entrained"],
                     model_pct_entrained_true=100 * m.entrained_true.mean(),
                     model_tau_LD_sd=e.tau_LS_LD.std(), obs_tau_LD_sd=obs.loc[sp, "tau_LD_sd"],
                     model_phase_ZT=e.phase_ZT_LD.mean(), obs_phase_ZT=obs.loc[sp, "phase_ZT"],
                     model_vector_r=m.vector_r_LD.mean(), obs_vector_r=obs.loc[sp, "vector_r"],
                     model_light_dark=m.light_dark_ratio.mean(), obs_light_dark=obs.loc[sp, "light_dark"],
                     model_change_DD_to_LD=100 * (m.entrained_LS.mean() - m.rhythmic_DD.mean()),
                     obs_change_DD_to_LD=obs.loc[sp, "pct_entrained"] - obs.loc[sp, "pct_rhythmic_DD"]))
t1 = pd.DataFrame(rows)
t1.to_csv(os.path.join(OUT, "table1_species_vs_observed.csv"), index=False)
print("\nTable 1: model at published beta vs observed")
print(t1.round(2).T.to_string())

# %% within species: paired DD -> LD (McNemar) and entrainment vs distance from 24 h
rows = []
for sp in PUB:
    m = pub[pub.species == sp]
    b = ((m.rhythmic_DD) & (~m.entrained_LS)).sum(); c = ((~m.rhythmic_DD) & (m.entrained_LS)).sum()
    p_mcnemar = stats.binomtest(int(min(b, c)), int(b + c), 0.5).pvalue if b + c > 0 else np.nan
    x = np.abs(m.tau_true_DD - 24).values; yv = m.entrained_LS.values.astype(int)
    ok = np.isfinite(x)
    # logistic regression by simple Newton steps (no statsmodels dependency needed)
    X = np.c_[np.ones(ok.sum()), x[ok]]; w = np.zeros(2)
    for _ in range(50):
        pr = 1 / (1 + np.exp(-X @ w)); g = X.T @ (yv[ok] - pr); H = (X * (pr * (1 - pr))[:, None]).T @ X
        w += np.linalg.solve(H + 1e-9 * np.eye(2), g)
    rows.append(dict(species=sp, rhythmic_DD_only=int(b), entrained_LD_only=int(c), mcnemar_p=p_mcnemar,
                     logit_slope_per_h=w[1], odds_ratio_per_h=np.exp(w[1]),
                     median_abs_dtau_entrained=np.nanmedian(x[yv == 1]) if (yv == 1).any() else np.nan,
                     median_abs_dtau_not=np.nanmedian(x[yv == 0]) if (yv == 0).any() else np.nan))
t2 = pd.DataFrame(rows); t2.to_csv(os.path.join(OUT, "table2_within_species.csv"), index=False)
print("\nTable 2: within species, paired DD->LD and entrainment vs |tau_DD - 24|")
print(t2.round(3).to_string(index=False))

# %% across beta: paired comparisons
rows = []
for sp in PUB:
    d = res[res.species == sp]
    wide = d.pivot(index="id", columns="beta", values="entrained_LS")
    for beta in wide.columns:
        k = int(wide[beta].sum()); n = len(wide)
        lo, hi = wilson(k, n)
        # paired difference vs beta = 0 (same individuals)
        diff = (wide[beta].astype(int) - wide[0.0].astype(int))
        rows.append(dict(species=sp, beta=beta, pct_entrained=100 * k / n, ci_lo=lo, ci_hi=hi,
                         paired_diff_vs_beta0=100 * diff.mean(), paired_se=100 * diff.std(ddof=1) / np.sqrt(n),
                         gained=int((diff > 0).sum()), lost=int((diff < 0).sum())))
t3 = pd.DataFrame(rows); t3.to_csv(os.path.join(OUT, "table3_beta_paired.csv"), index=False)
print("\nTable 3: % entrained across beta (paired vs beta=0)")
print(t3.round(1).to_string(index=False))

# key questions
sg = t3[t3.species == "S. grossa"].set_index("beta"); lc = t3[t3.species == "L. cornutus"].set_index("beta")
print(f"\nS. grossa: beta=0 -> {sg.loc[0.0,'pct_entrained']:.1f}%,  beta=0.75 -> {sg.loc[0.75,'pct_entrained']:.1f}%  "
      f"(paired diff {sg.loc[0.75,'paired_diff_vs_beta0']:+.1f} +- {sg.loc[0.75,'paired_se']:.1f})")
print(f"S. grossa best over the grid: {sg.pct_entrained.max():.1f}% at beta={sg.pct_entrained.idxmax()};  "
      f"L. cornutus at its published beta: {lc.loc[0.09,'pct_entrained']:.1f}%")

# %% cross-species beta swap: each population at each species' published beta
swap = np.full((3, 3), np.nan)
for i, sp in enumerate(PUB):
    for j, sp2 in enumerate(PUB):
        d = res[(res.species == sp) & np.isclose(res.beta, PUB[sp2])]
        swap[i, j] = 100 * d.entrained_LS.mean()
t4 = pd.DataFrame(swap, index=[f"{s} population" for s in PUB], columns=[f"beta of {s}" for s in PUB])
t4.to_csv(os.path.join(OUT, "table4_beta_swap.csv"))
print("\nTable 4: % entrained, each population run with each species' beta")
print(t4.round(1).to_string())

# %% figure 1: entrainment vs beta (single panel: % entrained, spider criterion)
fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))
for sp in PUB:
    d = t3[t3.species == sp]
    ax.plot(d.beta, d.pct_entrained, "o-", color=COL[sp], label=sp)
    ax.axvline(PUB[sp], color=COL[sp], ls=":", lw=1)
ax.set_title("Entrainment percentage vs beta for three spider species", fontsize=16.5)
ax.set_xlabel("Light sensitivity beta", fontsize=15)
ax.set_ylabel("Entrainment percentage", fontsize=15)
ax.set_xscale("symlog", linthresh=0.1)
ax.tick_params(axis="both", labelsize=15)
ax.grid(alpha=0.3)
ax.legend(fontsize=12, loc="lower right")
ax.set_ylim(0, 105)
ax.set_xlim(left=-0.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_entrainment_vs_beta.png"), dpi=1000)

# %% figure 2: tau_DD vs tau_LD at published beta
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
for a, sp in zip(ax, PUB):
    m = pub[pub.species == sp]
    a.scatter(m.tau_true_DD, m.tau_LS_LD, c=np.where(m.entrained_LS, COL[sp], "0.6"), s=12, alpha=0.7)
    a.axhline(24, color="k", lw=1); a.axhspan(22.8, 25.2, color="k", alpha=0.07)
    a.plot([16, 34], [16, 34], "k:", lw=1)
    a.set_title(sp, fontsize=13.5)
    a.set_xlabel("FRP in DD (hr)", fontsize=15); a.set_xlim(16, 34); a.set_ylim(14, 35)
    a.tick_params(axis="both", labelsize=15)
ax[0].set_ylabel("Period in LD (hr)", fontsize=15)
fig.suptitle("Period in LD vs FRP in DD for three spider species", fontsize=18)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_tauDD_vs_tauLD.png"), dpi=1000)

# %% figure 3: phase, vector strength, light/dark ratio, amplitude vs beta, with observed
fig, ax = plt.subplots(1, 4, figsize=(15, 3.8))
for sp in PUB:
    g = res[res.species == sp]
    ge = g[g.entrained_LS].groupby("beta")
    ga = g.groupby("beta")
    ax[0].plot(ge.phase_ZT_LD.mean().index, ge.phase_ZT_LD.mean().values, "o-", color=COL[sp], label=sp)
    ax[0].plot(PUB[sp], obs.loc[sp, "phase_ZT"], "*", color=COL[sp], ms=14, mec="k")
    ax[1].plot(ga.vector_r_LD.mean().index, ga.vector_r_LD.mean().values, "o-", color=COL[sp])
    ax[1].plot(PUB[sp], obs.loc[sp, "vector_r"], "*", color=COL[sp], ms=14, mec="k")
    ax[2].plot(ga.light_dark_ratio.mean().index, ga.light_dark_ratio.mean().values, "o-", color=COL[sp])
    ax[2].plot(PUB[sp], obs.loc[sp, "light_dark"], "*", color=COL[sp], ms=14, mec="k")
    ax[3].plot(ga.A_LD.mean().index, ga.A_LD.mean().values / ga.A_DD.mean().values, "o-", color=COL[sp])
for a, lab in zip(ax, ["mean phase ZT of entrained (h; dark from 12)", "vector strength", "light / dark activity ratio", "A_LD / A_DD"]):
    a.set_ylabel(lab); a.set_xscale("symlog", linthresh=0.1); a.set_xlabel("beta"); a.grid(alpha=0.3)
ax[0].axhline(12, color="k", lw=0.8, ls="--"); ax[0].legend(fontsize=7)
fig.suptitle("Phase, concentration and masking readouts across beta (stars = observed at published beta)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_phase_masking_vs_beta.png"), dpi=150)

# %% figure 4: DD -> LD change, model vs observed
fig, ax = plt.subplots(figsize=(7, 4))
x = np.arange(3); wdt = 0.2
for k, (lab, col_m, col_o) in enumerate([("rhythmic DD", "model_pct_rhythmic_DD", "obs_pct_rhythmic_DD"),
                                         ("entrained LD", "model_pct_entrained_LS", "obs_pct_entrained")]):
    ax.bar(x + (k - 1) * wdt * 2 + 0.0, t1[col_m], wdt, color=["0.3", "0.6"][k], label=f"model {lab}")
    ax.bar(x + (k - 1) * wdt * 2 + wdt, t1[col_o], wdt, color=["#d62728", "#ff9896"][k], label=f"observed {lab}")
ax.set_xticks(x); ax.set_xticklabels(t1.species); ax.set_ylabel("% of animals"); ax.legend(fontsize=8, ncol=2)
ax.set_title("DD rhythmicity (a fit) and LD entrainment (a prediction), model vs observed", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig4_DD_to_LD.png"), dpi=150)

# %% figure 5: beta swap heatmap
fig, ax = plt.subplots(figsize=(6, 4))
im = ax.imshow(swap, cmap="viridis", vmin=0, vmax=100)
for i in range(3):
    for j in range(3):
        ax.text(j, i, f"{swap[i, j]:.0f}%", ha="center", va="center", color="w" if swap[i, j] < 60 else "k", fontsize=11)
ax.set_xticks(range(3)); ax.set_xticklabels([f"beta = {PUB[s]}\n({s})" for s in PUB], fontsize=8)
ax.set_yticks(range(3)); ax.set_yticklabels([f"{s}\npopulation" for s in PUB], fontsize=8)
plt.colorbar(im, label="% entrained (spider criterion)")
ax.set_title("Is it beta or the gate? Each population with each species' beta", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig5_beta_swap.png"), dpi=150)
print("\nwrote tables 1-4 and figures 1-5")
