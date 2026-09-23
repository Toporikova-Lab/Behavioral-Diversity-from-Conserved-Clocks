# =============================================================================
# s1_spider_reference.py   --   Simulation S1, reference yardstick
#
# Per-animal waveform shape statistics from the 24 h aligned activity profiles
# (Figure 3 folder). These define what "looks like spider locomotion" means when
# the model waveforms are judged in s1e_shape_qc.py.
#
# DD profiles are peak-aligned to hour 12 by the Figure 3 script; LD profiles are
# ZT-aligned (ZT0 = lights on, dark from ZT12). The same shape metrics are applied
# to hourly-binned model output, aligned the same way.
# =============================================================================

# %% imports and paths
import numpy as np
import pandas as pd
import os

try:
    OUT = os.path.dirname(os.path.abspath(__file__))
except NameError:
    OUT = os.getcwd()

# folder holding the *_aligned_activity.csv files; edit if you move them
PROFILE_DIR = os.path.join(OUT, "spider_profiles")

FILES = {("L. cornutus", "DD"): "Larinioides_dd_aligned_activity.csv",
         ("L. cornutus", "LD"): "Larinioides_ld_aligned_activity.csv",
         ("A. pennsylvanica", "DD"): "Agelenopsis_dd_aligned_activity.csv",
         ("A. pennsylvanica", "LD"): "Agelenopsis_ld_aligned_activity.csv",
         ("S. grossa", "DD"): "Steatoda_dd_aligned_activity.csv",
         ("S. grossa", "LD"): "Steatoda_ld_aligned_activity.csv"}


# %% shape metrics on a 24-bin hourly profile
def shape_metrics(p):
    """Shape statistics of one 24-bin circular activity profile. Returns a dict."""
    p = np.asarray(p, float)
    if not np.isfinite(p).all() or p.sum() <= 0:
        return None
    lo, hi = p.min(), p.max()
    amp = hi - lo
    half = lo + 0.5 * amp
    fwhm = int((p > half).sum())                        # hours above half-maximum
    rolled = np.array([np.roll(p, -s)[:6].sum() for s in range(24)])
    peak6 = rolled.max() / p.sum()                      # share of activity in best 6 h
    jump = np.abs(np.diff(np.r_[p, p[0]])).max() / amp if amp > 0 else np.nan
    ang = 2 * np.pi * np.arange(24) / 24
    r = np.hypot((p * np.cos(ang)).sum(), (p * np.sin(ang)).sum()) / p.sum()
    # local maxima with prominence above 20% of amplitude, circular
    left, right = np.roll(p, 1), np.roll(p, -1)
    npk = int(((p > left) & (p >= right) & (p - lo > 0.2 * amp)).sum())
    return dict(fwhm_h=fwhm, peak6_frac=peak6, max_jump=jump,
                peak_to_mean=hi / p.mean(), vector_r=r, n_peaks=npk)


# %% run over every animal
rows, means = [], []
for (sp, cond), fname in FILES.items():
    d = pd.read_csv(os.path.join(PROFILE_DIR, fname))
    hours = d.iloc[:, 0].values
    for col in d.columns[1:]:
        m = shape_metrics(d[col].values)
        if m is None:
            continue
        rows.append(dict(species=sp, condition=cond, animal=col, **m))
    prof = d.iloc[:, 1:].mean(axis=1).values
    means.append(pd.DataFrame(dict(species=sp, condition=cond, hour=hours, mean_activity=prof)))

ref = pd.DataFrame(rows)
ref.to_csv(os.path.join(OUT, "spider_shape_reference.csv"), index=False)
pd.concat(means).to_csv(os.path.join(OUT, "spider_mean_profiles.csv"), index=False)

# %% summary table: what the model has to look like
summary = (ref.groupby(["species", "condition"])
              [["fwhm_h", "peak6_frac", "max_jump", "peak_to_mean", "vector_r", "n_peaks"]]
              .agg(["mean", "std"]).round(3))
print(summary.to_string())
print("\nn animals per group:")
print(ref.groupby(["species", "condition"]).size().to_string())
