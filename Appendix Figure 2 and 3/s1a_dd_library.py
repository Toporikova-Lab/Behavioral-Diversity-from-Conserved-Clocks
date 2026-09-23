# =============================================================================
# s1a_dd_library.py   --   Simulation S1, step A
#
# One DD pass over 20,000 random (k1,k2,k3) draws, stored as a reusable library.
# Every later calibration reweights this file instead of running ODEs again.
#
# DD only, so I = 0 and beta / alpha / L_baseline are inert by construction.
# The core (x,y,z) does not depend on the output gate, so both gates are carried
# alongside a single core integration.
#
# Metrics are recorded at each species' OWN observed DD record length, because
# record length differs nearly twofold between species and drives what the
# periodogram can resolve.
# =============================================================================

# %% imports, constants, and run size
import numpy as np
import pandas as pd
import os

# __file__ is not always defined when running cell-by-cell in Spyder
try:
    OUT = os.path.dirname(os.path.abspath(__file__))
except NameError:
    OUT = os.getcwd()
print("outputs will be written to:", OUT)

rng = np.random.default_rng(20260909)

# ---- fixed model constants --------------------------------------------------
v, K, n_hill, gamma, L_amp = 0.84, 1.0, 12, 12.0, 2.5

TAU_H_LC,  KAPPA_LC,  YTH_LC  = 2.0, 30.0, 0.5      # L. cornutus gate
TAU_H_APSG, KAPPA_APSG, YTH_APSG = 3.0, 10.0, 0.3   # A. pennsylvanica == S. grossa gate

# observed DD record lengths, days (Valid_Minutes/1440 in the comprehensive CSVs)
REC_LC, REC_AP, REC_SG = 8.45, 4.76, 7.46

# N = 2000 runs in about 15 s and is enough to check the script works.
# N = 20000 is the real library. Same seed, so the first 2000 rows are identical.
N = 20000
DT = 0.01
EQ_DAYS = 10.0

# %% draw the proposal population
# ---- proposal distribution --------------------------------------------------
# Broad and bounded, so importance weights stay bounded for any plausible target,
# but not so broad that most draws land in the dead region. k_geo uniform in log
# over 0.10-0.30 (tau roughly 14-42 h); asymmetry normal with sigma 0.7, wider
# than any sigma_asym we expect to fit but inside the region where the oscillator
# still runs.
SIG_PROP = 0.7
log_kgeo = rng.uniform(np.log(0.10), np.log(0.30), N)
k_geo = np.exp(log_kgeo)
a = rng.normal(0.0, SIG_PROP, (N, 3))
a = a - a.mean(axis=1, keepdims=True)
k1, k2, k3 = (k_geo * np.exp(a[:, 0]), k_geo * np.exp(a[:, 1]), k_geo * np.exp(a[:, 2]))
asym = np.sqrt((a ** 2).sum(axis=1))


# %% model definitions
def deriv(x, y, z, h_lc, h_ap):
    """Vectorised RHS for the whole population. DD, so no light and no masking."""
    dx = v * (K ** n_hill / (K ** n_hill + z ** n_hill)) - k1 * x
    dy = v * x - k2 * y
    dz = v * y - k3 * z
    dh_lc = (1.0 / (1.0 + np.exp(np.clip(KAPPA_LC * (y - YTH_LC), -60, 60))) - h_lc) / TAU_H_LC
    dh_ap = (1.0 / (1.0 + np.exp(np.clip(KAPPA_APSG * (y - YTH_APSG), -60, 60))) - h_ap) / TAU_H_APSG
    return dx, dy, dz, dh_lc, dh_ap


def rk4_step(s, dt):
    """One RK4 step on the packed state tuple s = (x, y, z, h_lc, h_ap)."""
    d1 = deriv(*s)
    s2 = tuple(si + 0.5 * dt * di for si, di in zip(s, d1))
    d2 = deriv(*s2)
    s3 = tuple(si + 0.5 * dt * di for si, di in zip(s, d2))
    d3 = deriv(*s3)
    s4 = tuple(si + dt * di for si, di in zip(s, d3))
    d4 = deriv(*s4)
    return tuple(si + dt / 6.0 * (a1 + 2 * a2 + 2 * a3 + a4)
                 for si, a1, a2, a3, a4 in zip(s, d1, d2, d3, d4))


def activity(y, h, y_th):
    """Locomotor output before masking. L_baseline is additive and cancels in max-mean."""
    return L_amp * (1.0 / (1.0 + np.exp(np.clip(-gamma * (y - y_th), -60, 60)))) * h


# %% 10-day equilibration   (slow cell, about 1 min at N = 20000)
# ---- 10-day equilibration ---------------------------------------------------
state = (np.full(N, 0.5), np.full(N, 0.5), np.full(N, 0.5), np.full(N, 0.5), np.full(N, 0.5))

y_hi = np.full(N, -np.inf)
y_lo = np.full(N, np.inf)
for i in range(int(EQ_DAYS * 24 / DT)):
    state = rk4_step(state, DT)
    if i * DT > (EQ_DAYS - 3) * 24:              # judge oscillation on the last 3 days
        y_hi = np.maximum(y_hi, state[1])
        y_lo = np.minimum(y_lo, state[1])

eq_flat = (y_hi - y_lo) < 1e-4
print(f"equilibration done. {eq_flat.sum()} of {N} draws show no oscillation in the "
      f"last 3 equilibration days; their final state is still used as the initial "
      f"condition, as required.")

# %% DD assay   (slow cell, about 1 min at N = 20000)
# ---- DD assay with running metrics and on-the-fly peak detection -------------
n_steps = int(REC_LC * 24 / DT)
checkpoints = {"AP": int(REC_AP * 24 / DT), "SG": int(REC_SG * 24 / DT), "LC": n_steps - 1}

acc = dict(y_max=np.full(N, -np.inf), y_sum=np.zeros(N),
           lc_max=np.full(N, -np.inf), lc_sum=np.zeros(N),
           ap_max=np.full(N, -np.inf), ap_sum=np.zeros(N),
           y_min=np.full(N, np.inf),
           npk=np.zeros(N, int), pk_first=np.zeros(N), pk_last=np.zeros(N))
MIN_PEAK_GAP = 2.0     # hours; below this a "peak" is numerical jitter, not a cycle
snap = {}

y2 = state[1].copy()
state = rk4_step(state, DT)
y1 = state[1].copy()

for i in range(2, n_steps):
    state = rk4_step(state, DT)
    x, y, z, h_lc, h_ap = state
    L_lc = activity(y, h_lc, YTH_LC)
    L_ap = activity(y, h_ap, YTH_APSG)

    acc["y_max"] = np.maximum(acc["y_max"], y);   acc["y_sum"] += y
    acc["y_min"] = np.minimum(acc["y_min"], y)
    acc["lc_max"] = np.maximum(acc["lc_max"], L_lc); acc["lc_sum"] += L_lc
    acc["ap_max"] = np.maximum(acc["ap_max"], L_ap); acc["ap_sum"] += L_ap

    tp = (i - 1) * DT
    hit = (y1 > y2) & (y1 >= y)                  # interior maximum of the protein
    hit &= (acc["npk"] == 0) | (tp - acc["pk_last"] > MIN_PEAK_GAP)
    if hit.any():
        acc["pk_first"] = np.where(hit & (acc["npk"] == 0), tp, acc["pk_first"])
        acc["pk_last"] = np.where(hit, tp, acc["pk_last"])
        acc["npk"] += hit

    if i in checkpoints.values():
        for key, idx in checkpoints.items():
            if i == idx:
                snap[key] = {k: (val.copy() if hasattr(val, "copy") else val)
                             for k, val in acc.items()}
                snap[key]["nstep"] = i
    y2, y1 = y1, y

# %% assemble, save, and summarise
# ---- assemble the library ---------------------------------------------------
df = pd.DataFrame(dict(k1=k1, k2=k2, k3=k3, k_geo=k_geo, asym=asym,
                       a1=a[:, 0], a2=a[:, 1], a3=a[:, 2], sigma_proposal=SIG_PROP))

for key, gate_max, gate_sum in [("LC", "lc_max", "lc_sum"),
                                ("AP", "ap_max", "ap_sum"),
                                ("SG", "ap_max", "ap_sum")]:
    s = snap[key]
    m = s["nstep"] - 2
    df[f"amp_prot_{key}"] = s["y_max"] - s["y_sum"] / m
    df[f"amp_act_{key}"] = s[gate_max] - s[gate_sum] / m
    swings = (s["y_max"] - s["y_min"]) > 1e-3
    df[f"tau_{key}"] = np.where((s["npk"] >= 2) & swings,
                                (s["pk_last"] - s["pk_first"]) / np.maximum(s["npk"] - 1, 1),
                                np.nan)
    df[f"npeaks_{key}"] = s["npk"]
    # published DD classification: protein must peak, activity amplitude must clear 0.4
    df[f"rhythmic_{key}"] = (s["npk"] >= 2) & swings & (df[f"amp_act_{key}"] > 0.4)

df.to_csv(os.path.join(OUT, "dd_library.csv"), index=False)

print(f"\nlibrary written: {len(df)} rows")
print(df[["tau_LC", "tau_AP", "tau_SG"]].describe().round(2).to_string())
for key, rec in [("LC", REC_LC), ("AP", REC_AP), ("SG", REC_SG)]:
    print(f"  gate/record {key} ({rec} d): rhythmic in {100 * df[f'rhythmic_{key}'].mean():.1f}% "
          f"of the proposal, activity amplitude median {df[f'amp_act_{key}'].median():.3f}")

# %% sanity checks -- these must all pass before the library is used
# 1. A. pennsylvanica and S. grossa share a gate and a core, so their periods differ
#    only through record length. Median difference should be a small fraction of an hour.
d_tau = (df.tau_AP - df.tau_SG).abs()
print(f"\ncheck 1  median |tau_AP - tau_SG| = {d_tau.median():.4f} h   (same core and gate, "
      f"differing only in record length; expect well under 0.1 h)")

# 2. Period must follow the scale of k with about 6% scatter.
ok = df.tau_LC.notna()
fit = np.polyfit(np.log(df.k_geo[ok]), np.log(df.tau_LC[ok]), 1)
resid = np.log(df.tau_LC[ok]) - np.polyval(fit, np.log(df.k_geo[ok]))
print(f"check 2  tau = {np.exp(fit[1]):.3f} * k_geo^{fit[0]:.3f}, "
      f"residual SD {100 * resid.std():.1f}%   (expect ~3.5, ~-1.03, ~7%)")

# 3. No period may be shorter than the peak refractory. Jitter peaks would show here.
print(f"check 3  shortest period {df.tau_LC.min():.2f} h   (must be well above "
      f"{MIN_PEAK_GAP} h; a value near it means jitter is being counted as cycles)")
