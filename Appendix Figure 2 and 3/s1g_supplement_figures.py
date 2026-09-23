# =============================================================================
# s1g_supplement_figures.py  --  figures for the reviewer response / Supplement S1
#
# beta is shown on the MANUSCRIPT scale. The manuscript defines I(t) = 1 in light,
# so its beta is the fractional rise of k2 in light. The simulation code used
# I = 0.7 (700 lux / 1000), so beta_manuscript = 0.7 * beta_simulation.
# All other quantities are read from track1_main/ (variant "fitted").
# =============================================================================

# %% inputs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

OUT = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
TRACK = os.path.join(OUT, "track1_main")
FIG = os.path.join(OUT, "supplement_figures"); os.makedirs(FIG, exist_ok=True)
I_SIM = 0.7

def track_file(name):
    """The sweep writes plain CSVs; the shipped copies are gzipped. Accept either."""
    gz = os.path.join(TRACK, name + ".gz")
    return gz if os.path.exists(gz) else os.path.join(TRACK, name)

res = pd.read_csv(track_file("results_individual.csv"))
res = res[res.variant == "fitted"].copy()
res["beta_ms"] = res.beta * I_SIM
prof = pd.read_csv(track_file("profiles_hourly.csv")); prof = prof[prof.variant == "fitted"]
ref_mean = pd.read_csv(os.path.join(OUT, "spider_mean_profiles.csv"))

SP = ["L. cornutus", "A. pennsylvanica", "S. grossa"]
PUB = {"L. cornutus": 0.09, "A. pennsylvanica": 0.06, "S. grossa": 0.75}          # manuscript beta
OBS = {"L. cornutus": 90.6, "A. pennsylvanica": 63.4, "S. grossa": 48.4}
COL = {"L. cornutus": "#0072B2", "A. pennsylvanica": "#E69F00", "S. grossa": "#009E73"}
MK = {"L. cornutus": "o", "A. pennsylvanica": "s", "S. grossa": "^"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def pct_curve(sp, col):
    g = res[res.species == sp].groupby("beta_ms")[col].mean() * 100
    return g.index.values, g.values


def at_published(sp, col):
    x, y = pct_curve(sp, col)
    return np.interp(PUB[sp], x, y)


# %% Figure S1-1: entrainment vs beta (single panel, spider-scoring criterion; matches fig1_entrainment_vs_beta.png)
fig, ax = plt.subplots(1, 1, figsize=(6.5, 4.6))
for sp in SP:
    x, y = pct_curve(sp, "entrained_LS")
    ax.plot(x, y, marker=MK[sp], ms=5, lw=1.8, color=COL[sp], label=sp)
    ax.axvline(PUB[sp], color=COL[sp], ls=":", lw=1)
ax.set_xscale("symlog", linthresh=0.05); ax.set_xlim(-0.005, 1.6); ax.set_ylim(0, 102)
ax.set_xticks([0, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5]); ax.set_xticklabels(["0", "0.05", "0.1", "0.25", "0.5", "1", "1.5"])
ax.grid(alpha=0.3, lw=0.6)
ax.set_title("Entrainment percentage vs beta for three spider species", fontsize=15)
ax.set_xlabel("Light sensitivity beta", fontsize=13.5)
ax.set_ylabel("Entrainment percentage", fontsize=13.5)
ax.tick_params(axis="both", labelsize=13.5)
ax.legend(fontsize=11.25, frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "FigS1_1_entrainment_vs_beta.png"), dpi=1000)

# %% Figure S1-2: which clocks entrain? free-running period vs period in LD, at the fitted beta
# (title/axis wording matches fig2_tauDD_vs_tauLD.png)
fig, ax = plt.subplots(1, 3, figsize=(10, 3.5), sharey=True)
for a, sp in zip(ax, SP):
    b_sim = min(res.beta.unique(), key=lambda b: abs(b * I_SIM - PUB[sp]))
    m = res[(res.species == sp) & np.isclose(res.beta, b_sim)]
    a.axhspan(22.8, 25.2, color="0.88", zorder=0)
    a.plot([16, 34], [16, 34], ":", color="0.5", lw=1)
    a.axhline(24, color="0.3", lw=0.8)
    e = m.entrained_LS
    a.scatter(m.tau_true_DD[~e], m.tau_LS_LD[~e], s=9, color="0.65", alpha=0.7, label="not entrained")
    a.scatter(m.tau_true_DD[e], m.tau_LS_LD[e], s=9, color=COL[sp], alpha=0.85, label="entrained")
    a.set_xlim(16, 34); a.set_ylim(14, 35)
    a.set_xlabel("FRP in DD (hr)", fontsize=13.5)
    a.set_title(sp, fontsize=13.5)
    a.tick_params(axis="both", labelsize=13.5)
ax[0].set_ylabel("Period in LD (hr)", fontsize=13.5)
ax[0].legend(fontsize=11.25, frameon=False, loc="upper left")
fig.suptitle("Period in LD vs FRP in DD for three spider species", fontsize=18)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "FigS1_2_which_clocks_entrain.png"), dpi=1000)

# %% Figure S1-3: model and spider daily waveforms
H = [f"h{k}" for k in range(24)]
fig, ax = plt.subplots(2, 3, figsize=(10, 5.2), sharex=True, sharey=True)
for c, sp in enumerate(SP):
    b_sim = min(res.beta.unique(), key=lambda b: abs(b * I_SIM - PUB[sp]))
    for r, cond in enumerate(["DD", "LD"]):
        m = prof[(prof.species == sp) & (prof.condition == cond)]
        if cond == "LD":
            m = m[np.isclose(m.beta, b_sim)]
        P = m[H].values; mm = P.mean(0)
        for row in P[::12]:
            ax[r, c].plot(range(24), row / row.max(), color="0.82", lw=0.5)
        ax[r, c].plot(range(24), mm / mm.max(), color="k", lw=2, label="model, population mean")
        s = ref_mean[(ref_mean.species == sp) & (ref_mean.condition == cond)].mean_activity.values
        ax[r, c].plot(range(24), s / s.max(), color=COL[sp], lw=2, label="spiders, population mean")
        if cond == "LD":
            ax[r, c].axvspan(12, 24, color="0.92", zorder=0)
        ax[r, c].set_title(f"{sp}   {cond}", fontsize=9, loc="left")
ax[0, 0].set_ylabel("activity (normalised)"); ax[1, 0].set_ylabel("activity (normalised)")
ax[1, 1].set_xlabel("hour of cycle   (DD: peak aligned to hour 12;   LD: ZT, dark phase shaded)")
ax[0, 0].legend(fontsize=7.5, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "FigS1_3_waveforms.png"), dpi=300)

# %% Figure S1-4: how the observed entrainment percentage decomposes (data only)
rows = []
for sp, f in [("L. cornutus", "Larinioides"), ("A. pennsylvanica", "Agelenopsis"), ("S. grossa", "Steatoda")]:
    c = pd.read_csv(os.path.join(OUT, "spider_data", f"{f}_circular_stats_v2.csv"))
    sig = c.is_significant
    rows.append(dict(species=sp, n=len(c), significant=100 * sig.mean(),
                     in_window_given_sig=100 * c.is_entrained[sig].mean(), entrained=100 * (sig & c.is_entrained).mean()))
dec = pd.DataFrame(rows); dec.to_csv(os.path.join(FIG, "observed_decomposition.csv"), index=False)
fig, ax = plt.subplots(figsize=(6.2, 3.4))
x = np.arange(3); w = 0.26
for k, (col, lab, shade) in enumerate([("significant", "significant rhythm (any period)", 0.75),
                                       ("in_window_given_sig", "period within 22.8 to 25.2 h, given a rhythm", 0.5),
                                       ("entrained", "both: scored as entrained", 0.0)]):
    ax.bar(x + (k - 1) * w, dec[col], w, color=[COL[s] for s in SP], alpha=1 - shade * 0.7, edgecolor="k", lw=0.5, label=lab)
    for xi, v in zip(x + (k - 1) * w, dec[col]):
        ax.text(xi, v + 1.5, f"{v:.0f}", ha="center", fontsize=7.5)
ax.set_xticks(x); ax.set_xticklabels([f"{s}\n(n = {n})" for s, n in zip(dec.species, dec.n)]); ax.set_ylim(0, 110)
ax.set_ylabel("spiders in LD (%)"); ax.legend(fontsize=7.5, frameon=False, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "FigS1_4_observed_decomposition.png"), dpi=300)

# %% numbers used in the text, on the manuscript beta scale
print("model % entrained at the fitted beta (manuscript scale), interpolated on the grid:")
for sp in SP:
    print(f"  {sp:<18} beta {PUB[sp]:.2f}: spider criterion {at_published(sp,'entrained_LS'):.1f}%   "
          f"phase-locked {at_published(sp,'entrained_true'):.1f}%   observed {OBS[sp]}")
print("\nbeta (manuscript scale) needed to reach the observed percentage:")
for sp in SP:
    x, y = pct_curve(sp, "entrained_LS")
    need = np.interp(OBS[sp], y, x) if y.max() >= OBS[sp] else np.nan
    print(f"  {sp:<18} {need:.2f}   (grid max {y.max():.1f}% at beta {x[np.argmax(y)]:.2f})")
print("\nS. grossa across the grid (manuscript beta -> %):")
x, y = pct_curve("S. grossa", "entrained_LS"); print("  " + "  ".join(f"{a:.2f}:{b:.0f}" for a, b in zip(x, y)))
print("\nobserved decomposition:"); print(dec.round(1).to_string(index=False))
