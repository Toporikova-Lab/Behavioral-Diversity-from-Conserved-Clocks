# =============================================================================
# gate_shape_scan.py  --  does any output-gate setting reproduce the spider DD
# waveform width? Representative clock per species = symmetric k at fitted k_bar.
# DD only, so beta / alpha / L_base are inert except L_base's effect on peak/mean.
# =============================================================================
import numpy as np, pandas as pd, itertools, os
from scipy.integrate import odeint
from scipy.signal import find_peaks

v, K, n, gamma, L_amp = 0.84, 1.0, 12, 12.0, 2.5
fit = pd.read_csv("calibration_fit.csv"); fit = fit[fit.variant == "fitted"].set_index("species")
ref = pd.read_csv("spider_shape_reference.csv")
ref = ref[ref.condition == "DD"].groupby("species")[["fwhm_h", "peak6_frac", "peak_to_mean"]].mean()
L_BASE = {"L. cornutus": 0.0, "A. pennsylvanica": 0.2, "S. grossa": 0.15}
PUB = {"L. cornutus": (2.0, 30, 0.5), "A. pennsylvanica": (3.0, 10, 0.3), "S. grossa": (3.0, 10, 0.3)}


def run(k, tau_h, kappa, y_th, L_base, days=12):
    def rhs(s, t):
        x, y, z, h = s
        return [v * (K**n / (K**n + z**n)) - k * x, v * x - k * y, v * y - k * z,
                (1 / (1 + np.exp(np.clip(kappa * (y - y_th), -60, 60))) - h) / tau_h]
    s = odeint(rhs, [0.5, 0.5, 0.5, 0.5], np.linspace(0, 240, 2401))[-1]
    t = np.linspace(0, days * 24, days * 240 + 1)
    s = odeint(rhs, s, t)
    y, h = s[:, 1], s[:, 3]
    L = L_base + L_amp / (1 + np.exp(-gamma * (y - y_th))) * h
    return t, y, L


def hourly_shape(t, L, tau):
    phase = ((t + 1e-6) % tau) / tau
    b = np.minimum((phase * 24).astype(int), 23)
    p = np.bincount(b, weights=L, minlength=24) / np.maximum(np.bincount(b, minlength=24), 1)
    lo, hi = p.min(), p.max(); amp = hi - lo
    fwhm = int((p > lo + 0.5 * amp).sum())
    peak6 = max(np.roll(p, -s)[:6].sum() for s in range(24)) / p.sum()
    return fwhm, peak6, hi / p.mean()


rows = []
for sp in fit.index:
    k = fit.loc[sp, "k_bar"]
    t, y, L = run(k, *PUB[sp], L_BASE[sp])
    pk = t[find_peaks(y, distance=20)[0]]; tau = np.diff(pk).mean()
    print(f"{sp}: k_bar={k:.3f} tau={tau:.2f} h   y range {y.min():.2f}-{y.max():.2f}   "
          f"spider DD target fwhm {ref.loc[sp,'fwhm_h']:.1f} h, peak6 {ref.loc[sp,'peak6_frac']:.2f}, p2m {ref.loc[sp,'peak_to_mean']:.1f}")
    for tau_h, kappa, y_th in itertools.product([0.5, 1.0, 1.5, 2.0, 3.0], [10, 30, 60], [0.25, 0.3, 0.4, 0.5, 0.6]):
        t, y, L = run(k, tau_h, kappa, y_th, L_BASE[sp])
        fw, p6, p2m = hourly_shape(t, L, tau)
        rows.append(dict(species=sp, tau_h=tau_h, kappa=kappa, y_th=y_th, fwhm=fw, peak6=p6, p2m=p2m,
                         A_L=L.max() - L.mean(), published=(tau_h, kappa, y_th) == PUB[sp]))
d = pd.DataFrame(rows); d.to_csv("gate_shape_scan.csv", index=False)
for sp in fit.index:
    g = d[d.species == sp].copy()
    tg = ref.loc[sp]
    g["dist"] = ((g.fwhm - tg.fwhm_h) / 2.5) ** 2 + ((g.peak6 - tg.peak6_frac) / 0.1) ** 2 + ((g.p2m - tg.peak_to_mean) / 1.5) ** 2
    print(f"\n{sp}  published row, then 6 closest to spider DD shape (A_L must stay > 0.4):")
    print(pd.concat([g[g.published], g[~g.published].nsmallest(6, "dist")]).round(3).to_string(index=False))
