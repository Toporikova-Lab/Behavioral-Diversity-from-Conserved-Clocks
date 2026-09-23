# =============================================================================
# s1e_shape_qc.py   --   Simulation S1, step E: waveform quality control
#
# Checks every simulated spider's activity profile and compares the population's
# waveform shape to the real spiders' aligned profiles, so that the beta effect is
# judged on oscillators that look like spider locomotion, not on artefacts.
#
# Per individual: flags for non-oscillatory clock, damped output, extremely sharp
# edges, multi-peaked days, and shape outside the spider range.
# Per population: shape statistics against the spider reference, at the published
# beta and across the beta grid.
# =============================================================================

# %% imports and inputs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

try:
    OUT = os.path.dirname(os.path.abspath(__file__))
except NameError:
    OUT = os.getcwd()

def find(name):
    """The reference CSVs live beside this script or one level up (the track folders)."""
    here = os.path.join(OUT, name)
    return here if os.path.exists(here) else os.path.join(os.path.dirname(OUT), name)

res = pd.read_csv(os.path.join(OUT, "results_individual.csv.gz" if os.path.exists(os.path.join(OUT, "results_individual.csv.gz")) else "results_individual.csv"))
prof = pd.read_csv(os.path.join(OUT, "profiles_hourly.csv.gz" if os.path.exists(os.path.join(OUT, "profiles_hourly.csv.gz")) else "profiles_hourly.csv"))
ref = pd.read_csv(find("spider_shape_reference.csv"))
ref_mean = pd.read_csv(find("spider_mean_profiles.csv"))
PUB = {"L. cornutus": 0.09, "A. pennsylvanica": 0.06, "S. grossa": 0.75}
COL = {"L. cornutus": "#1f77b4", "A. pennsylvanica": "#ff7f0e", "S. grossa": "#2ca02c"}
H = [f"h{k}" for k in range(24)]
SHAPE = ["fwhm_h", "peak6_frac", "max_jump", "peak_to_mean", "n_peaks"]

VARIANT = os.environ.get("S1_PLOT_VARIANT", res.variant.iloc[0])
res = res[res.variant == VARIANT].copy()
prof = prof[prof.variant == VARIANT].copy()
res["published"] = res.apply(lambda r: np.isclose(r.beta, PUB[r.species]), axis=1)
pub = res[res.published].copy()

# %% per-individual flags at the published beta
def spider_range(sp, cond, metric, lo=2.5, hi=97.5):
    x = ref[(ref.species == sp) & (ref.condition == cond)][metric].dropna()
    return np.percentile(x, lo), np.percentile(x, hi)

flags = pub[["species", "id", "beta"]].copy()
flags["non_oscillatory_clock"] = (pub.n_peaks_prot_DD < 2) | (pub.A_prot_DD < 1e-3)
flags["arrhythmic_activity_DD"] = ~pub.rhythmic_DD
flags["damped_DD"] = pub.damping_DD < 0.8
flags["damped_LD"] = pub.damping_LD < 0.8
flags["beating_LD"] = pub.amp_cv_LD > 0.25          # day-to-day amplitude wobble: relative coordination
flags["sharp_edge_LD"] = pub.edge_sharp_LD > 0.5   # > half the amplitude within 0.1 h
flags["sharp_edge_DD"] = pub.edge_sharp_DD > 0.5
flags["multi_peak_LD"] = pub.n_peaks_LD > 2          # clock peak + masking edge is 2; more is odd
for m in ["fwhm_h", "peak6_frac", "peak_to_mean"]:
    for cond in ["DD", "LD"]:
        lo, hi = zip(*[spider_range(sp, cond, m) for sp in pub.species])
        val = pub[f"{m}_{cond}"].values
        flags[f"{m}_{cond}_outside_spider_range"] = (val < np.array(lo)) | (val > np.array(hi))
flags.to_csv(os.path.join(OUT, "qc_flags.csv"), index=False)

print("per-individual QC flags at published beta (% of individuals):")
print((flags.drop(columns=["id", "beta"]).groupby("species").mean() * 100).round(1).T.to_string())

# %% population shape vs spider reference, published beta
rows = []
for sp in PUB:
    for cond in ["DD", "LD"]:
        m = pub[pub.species == sp]
        s = ref[(ref.species == sp) & (ref.condition == cond)]
        for metric in SHAPE:
            rows.append(dict(species=sp, condition=cond, metric=metric,
                             model_mean=m[f"{metric}_{cond}"].mean(), model_sd=m[f"{metric}_{cond}"].std(),
                             spider_mean=s[metric].mean(), spider_sd=s[metric].std(), n_spiders=len(s)))
        rows.append(dict(species=sp, condition=cond, metric="vector_r",
                         model_mean=m.vector_r_LD.mean() if cond == "LD" else np.nan,
                         model_sd=m.vector_r_LD.std() if cond == "LD" else np.nan,
                         spider_mean=s.vector_r.mean(), spider_sd=s.vector_r.std(), n_spiders=len(s)))
cmp = pd.DataFrame(rows)
cmp["z_of_model_mean"] = (cmp.model_mean - cmp.spider_mean) / cmp.spider_sd
cmp.to_csv(os.path.join(OUT, "qc_shape_comparison.csv"), index=False)
print("\nmodel vs spider waveform shape (z = model mean in spider SD units):")
print(cmp.round(2).to_string(index=False))

# %% figure 1: waveforms, model vs spider
fig, ax = plt.subplots(2, 3, figsize=(13, 6.5), sharex=True)
for c, sp in enumerate(PUB):
    for r, cond in enumerate(["DD", "LD"]):
        m = prof[(prof.species == sp) & (prof.condition == cond)]
        if cond == "LD":
            m = m[np.isclose(m.beta, PUB[sp])]
        P = m[H].values
        P = P / np.maximum(P.max(axis=1, keepdims=True), 1e-9)
        for row in P[::max(1, len(P) // 40)]:
            ax[r, c].plot(range(24), row, color="0.8", lw=0.6)
        mm = m[H].values.mean(axis=0)
        ax[r, c].plot(range(24), mm / mm.max(), "k", lw=2.5, label=f"model mean (n={len(m)})")
        s = ref_mean[(ref_mean.species == sp) & (ref_mean.condition == cond)].mean_activity.values
        ax[r, c].plot(range(24), s / s.max(), color="r", lw=2.5, label="spider mean")
        if cond == "LD":
            ax[r, c].axvspan(12, 24, color="0.9", zorder=0)
        ax[r, c].set_title(f"{sp}   {cond}", fontsize=10)
ax[0, 0].legend(fontsize=8)
ax[0, 0].set_ylabel("normalised activity"); ax[1, 0].set_ylabel("normalised activity")
ax[1, 1].set_xlabel("hour   (DD: folded at own period, peak at 12;   LD: ZT, dark shaded)")
fig.suptitle(f"[{VARIANT}] Waveform QC at published beta: simulated spiders (grey, black mean) vs real spiders (red)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "qc_waveforms.png"), dpi=150)

# %% figure 2: shape distributions, model vs spider
fig, ax = plt.subplots(2, 3, figsize=(13, 6))
for c, (metric, xlabel) in enumerate([("fwhm_h", "FWHM of hourly profile (h)"),
                                      ("peak6_frac", "share of activity in best 6 h"),
                                      ("peak_to_mean", "peak / mean")]):
    for r, cond in enumerate(["DD", "LD"]):
        for sp in PUB:
            mv = pub[pub.species == sp][f"{metric}_{cond}"].dropna()
            sv = ref[(ref.species == sp) & (ref.condition == cond)][metric].dropna()
            bins = np.linspace(min(mv.min(), sv.min()), max(mv.max(), sv.max()), 20)
            ax[r, c].hist(mv, bins=bins, density=True, histtype="step", lw=2, color=COL[sp], label=f"{sp} model")
            ax[r, c].hist(sv, bins=bins, density=True, histtype="stepfilled", alpha=0.25, color=COL[sp], label=f"{sp} spiders")
        ax[r, c].set_title(cond, fontsize=10)
    ax[1, c].set_xlabel(xlabel)
ax[0, 0].set_ylabel("density"); ax[1, 0].set_ylabel("density")
ax[0, 0].legend(fontsize=6, ncol=2)
fig.suptitle("Shape distributions at published beta: model (lines) vs spiders (filled)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "qc_shape_distributions.png"), dpi=150)

# %% figure 3: does beta change the population waveform?
g = res.groupby(["species", "beta"]).agg(
    fwhm=("fwhm_h_LD", "mean"), peak6=("peak6_frac_LD", "mean"), p2m=("peak_to_mean_LD", "mean"),
    vec=("vector_r_LD", "mean"), amp=("A_LD", "mean"), damp=("damping_LD", "mean"),
    cv=("amp_cv_LD", "mean"), edge=("edge_sharp_LD", "mean")).reset_index()
panels = [("fwhm", "FWHM in LD (h)"), ("peak6", "best-6 h share"), ("p2m", "peak / mean"),
          ("vec", "vector strength"), ("amp", "activity amplitude A_LD"), ("cv", "day-to-day amplitude CV")]
fig, ax = plt.subplots(2, 3, figsize=(13, 6), sharex=True)
for a, (col, lab) in zip(ax.ravel(), panels):
    for sp in PUB:
        d = g[g.species == sp]
        a.plot(d.beta, d[col], "o-", color=COL[sp], label=sp)
        a.axvline(PUB[sp], color=COL[sp], ls=":", lw=1)
    a.set_ylabel(lab); a.set_xscale("symlog", linthresh=0.1); a.grid(alpha=0.3)
for a in ax[1]:
    a.set_xlabel("beta (symlog)")
ax[0, 0].legend(fontsize=8)
fig.suptitle("Population waveform shape in LD across the beta grid (dotted = published beta)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "qc_shape_vs_beta.png"), dpi=150)

# %% verdict
print("\nVERDICT")
bad = flags[["non_oscillatory_clock", "damped_DD", "damped_LD", "sharp_edge_LD", "sharp_edge_DD"]].mean()
for k, val in bad.items():
    print(f"  {k:<24} {100 * val:5.1f}%  {'OK' if val < 0.05 else 'INVESTIGATE'}")
worst = cmp.loc[cmp.z_of_model_mean.abs().idxmax()]
print(f"  largest shape mismatch: {worst.species} {worst.condition} {worst.metric}: "
      f"model {worst.model_mean:.2f} vs spiders {worst.spider_mean:.2f} +- {worst.spider_sd:.2f} (z = {worst.z_of_model_mean:.1f})")
print("wrote qc_flags.csv, qc_shape_comparison.csv, qc_waveforms.png, qc_shape_distributions.png, qc_shape_vs_beta.png")
