# =============================================================================
# s1b_calibrate.py   --   Simulation S1, step B
#
# Fit the three population parameters per species by REWEIGHTING dd_library.csv,
# so no ODE is solved here. Targets, all from the observed DD records:
#     fraction rhythmic,  mean tau (rhythmic animals only),  SD tau (rhythmic only)
# Three targets, three parameters:  k_bar, sigma_scale, sigma_asym.
#
# =============================================================================

# %% imports and paths
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import os

try:
    OUT = os.path.dirname(os.path.abspath(__file__))
except NameError:
    OUT = os.getcwd()
DATA = os.path.join(OUT, "spider_data")


SPECIES = {"L. cornutus": ("LC", "LC_spider_analysis_comprehensive_with_LD_split.csv"),
           "A. pennsylvanica": ("AP", "Ag_spider_analysis_comprehensive_with_LD_split.csv"),
           "S. grossa": ("SG", "Sg_spider_analysis_comprehensive_with_LD_split.csv")}

# %% targets, recomputed from the CSVs every time (never typed in)
def dd_targets(path):
    """Observed DD rhythmic fraction and period statistics, rhythmic animals only."""
    d = pd.read_csv(path)
    dd = d[d.Condition == "DD"]
    # animals with no beam crossings at all (Mean_Activity == 0) cannot inform a clock
    # model; they are dropped from the denominator. Only L. cornutus DD has any (11 of 45).
    n_empty = int((dd.Mean_Activity <= 0).sum())
    dd = dd[dd.Mean_Activity > 0]
    rhythmic = dd[dd.Period_p_value < 0.05]
    n, nr = len(dd), len(rhythmic)
    f = nr / n
    tau, sd = rhythmic.Period_hours.mean(), rhythmic.Period_hours.std(ddof=1)
    return dict(n=n, n_rhythmic=nr, n_excluded_zero_activity=n_empty, f=f, se_f=np.sqrt(f * (1 - f) / n),
                tau=tau, se_tau=sd / np.sqrt(nr),
                sd=sd, se_sd=sd / np.sqrt(2 * (nr - 1)))

targets = {sp: dd_targets(os.path.join(DATA, fn)) for sp, (_, fn) in SPECIES.items()}
print("targets from the DD records (Condition == DD, Period_p_value < 0.05):")
for sp, t in targets.items():
    print(f"  {sp:<18} rhythmic {t['n_rhythmic']}/{t['n']} = {100*t['f']:.1f}%   "
          f"tau {t['tau']:.2f} +- {t['sd']:.2f} h   (excluded {t['n_excluded_zero_activity']} zero-activity animals)")

# %% the library and the proposal it was drawn from
lib = pd.read_csv(os.path.join(OUT, "dd_library.csv"))
A2 = (lib[["a1", "a2", "a3"]].values ** 2).sum(axis=1)
LOG_KGEO = np.log(lib.k_geo.values)
SIG_PROP = float(lib.sigma_proposal.iloc[0])
LOG_RANGE = np.log(0.30) - np.log(0.10)          # uniform-in-log proposal for k_geo


def importance_weights(k_bar, sigma_scale, sigma_asym):
    """Normalised weights turning the proposal draws into draws from the candidate."""
    # asymmetry: demeaned 3-vector lives on a 2-D plane, hence the 2 in the normaliser
    lw = -A2 / (2 * sigma_asym ** 2) + A2 / (2 * SIG_PROP ** 2) + 2 * np.log(SIG_PROP / sigma_asym)
    # scale: lognormal candidate over uniform-in-log proposal
    lw += -(LOG_KGEO - np.log(k_bar)) ** 2 / (2 * sigma_scale ** 2) \
          - np.log(sigma_scale * np.sqrt(2 * np.pi)) + np.log(LOG_RANGE)
    w = np.exp(lw - lw.max())
    return w / w.sum()


def weighted_stats(key, w):
    """Rhythmic fraction and tau mean/SD over rhythmic draws, under weights w."""
    r = lib[f"rhythmic_{key}"].values.astype(float)
    tau = np.nan_to_num(lib[f"tau_{key}"].values)
    f = (w * r).sum()
    wr = w * r
    if wr.sum() <= 0:
        return f, np.nan, np.nan, 1 / (w ** 2).sum()
    wr = wr / wr.sum()
    mu = (wr * tau).sum()
    sd = np.sqrt((wr * (tau - mu) ** 2).sum())
    return f, mu, sd, 1 / (w ** 2).sum()


def loss(params, key, t):
    """Sum of squared z-scores of the three targets. Large penalty off-domain."""
    k_bar, sigma_scale, sigma_asym = params
    if not (0.08 < k_bar < 0.35 and 0.005 < sigma_scale < 0.5 and 0.02 < sigma_asym < 1.2):
        return 1e6
    f, mu, sd, _ = weighted_stats(key, importance_weights(k_bar, sigma_scale, sigma_asym))
    if not np.isfinite(mu):
        return 1e6
    return ((f - t["f"]) / t["se_f"]) ** 2 + ((mu - t["tau"]) / t["se_tau"]) ** 2 \
           + ((sd - t["sd"]) / t["se_sd"]) ** 2


# %% coarse grid, then local refinement, per species
K_GRID = np.linspace(0.12, 0.24, 25)
SS_GRID = np.linspace(0.02, 0.30, 15)
SA_GRID = np.linspace(0.05, 0.80, 16)

fits, surfaces = [], {}
for sp, (key, _) in SPECIES.items():
    t = targets[sp]
    grid = np.full((len(K_GRID), len(SS_GRID), len(SA_GRID)), np.nan)
    for i, kb in enumerate(K_GRID):
        for j, ss in enumerate(SS_GRID):
            for l, sa in enumerate(SA_GRID):
                grid[i, j, l] = loss((kb, ss, sa), key, t)
    i, j, l = np.unravel_index(np.nanargmin(grid), grid.shape)
    x0 = (K_GRID[i], SS_GRID[j], SA_GRID[l])
    res = minimize(loss, x0, args=(key, t), method="Nelder-Mead",
                   options=dict(xatol=1e-4, fatol=1e-4, maxiter=2000))
    kb, ss, sa = res.x
    f, mu, sd, ess = weighted_stats(key, importance_weights(kb, ss, sa))
    surfaces[sp] = (grid, (i, j, l))
    fits.append(dict(species=sp, key=key, variant="fitted",
                     k_bar=kb, sigma_scale=ss, sigma_asym=sa, loss=res.fun, ess=ess,
                     f_sim=f, tau_sim=mu, sd_sim=sd,
                     f_obs=t["f"], tau_obs=t["tau"], sd_obs=t["sd"], n_obs=t["n"],
                     n_excluded_zero_activity=t["n_excluded_zero_activity"]))


fit = pd.DataFrame(fits)
fit.to_csv(os.path.join(OUT, "calibration_fit.csv"), index=False)

print("\nfitted population parameters (loss is a sum of three squared z-scores):")
cols = ["species", "variant", "k_bar", "sigma_scale", "sigma_asym", "loss", "ess",
        "f_sim", "f_obs", "tau_sim", "tau_obs", "sd_sim", "sd_obs"]
print(fit[cols].round(3).to_string(index=False))

# %% loss-surface figure: two slices through each optimum
fig, ax = plt.subplots(3, 2, figsize=(9, 10))
for row, (sp, (grid, (i, j, l))) in enumerate(surfaces.items()):
    z1 = np.log10(grid[:, j, :].T)          # k_bar x sigma_asym at best sigma_scale
    z2 = np.log10(grid[i, :, :].T)          # sigma_scale x sigma_asym at best k_bar
    im = ax[row, 0].contourf(K_GRID, SA_GRID, z1, levels=20, cmap="viridis_r")
    ax[row, 0].plot(K_GRID[i], SA_GRID[l], "r+", ms=12, mew=2)
    ax[row, 1].contourf(SS_GRID, SA_GRID, z2, levels=20, cmap="viridis_r")
    ax[row, 1].plot(SS_GRID[j], SA_GRID[l], "r+", ms=12, mew=2)
    ax[row, 0].set_ylabel(f"{sp}\n$\\sigma_{{asym}}$")
    plt.colorbar(im, ax=ax[row, 1], label="log10 loss")
ax[2, 0].set_xlabel("$\\bar{k}$")
ax[2, 1].set_xlabel("$\\sigma_{scale}$")
ax[0, 0].set_title("k_bar vs sigma_asym (sigma_scale at optimum)", fontsize=9)
ax[0, 1].set_title("sigma_scale vs sigma_asym (k_bar at optimum)", fontsize=9)
fig.suptitle("Calibration loss surfaces (red cross = grid optimum before refinement)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "calibration_loss_surfaces.png"), dpi=150)
print("\nwrote calibration_fit.csv and calibration_loss_surfaces.png")
