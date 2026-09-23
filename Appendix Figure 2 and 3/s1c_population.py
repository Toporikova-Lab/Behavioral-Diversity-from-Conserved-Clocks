# =============================================================================
# s1c_population.py   --   Simulation S1, step C
#
# Draw the simulated spiders. N per species, from the fitted (k_bar, sigma_scale,
# sigma_asym) in calibration_fit.csv. The SAME underlying standard-normal draws
# are used for every species and variant, so individual i sits at the same
# quantile of each population and between-species comparisons are paired.
# Nobody is excluded: draws whose clock does not oscillate are the arrhythmic
# spiders and stay in every denominator.
# =============================================================================

# %% imports
import numpy as np
import pandas as pd
import os

try:
    OUT = os.path.dirname(os.path.abspath(__file__))
except NameError:
    OUT = os.getcwd()

N = 500
SEED = 20260910
rng = np.random.default_rng(SEED)

# %% common random numbers, drawn once
z_scale = rng.normal(0, 1, N)              # scale direction
z_asym = rng.normal(0, 1, (N, 3))          # asymmetry direction, demeaned below
z_asym -= z_asym.mean(axis=1, keepdims=True)
phase_offset = rng.uniform(0, 24, N)       # extra DD hours before lights-on

# %% build every population from the same draws
fit = pd.read_csv(os.path.join(OUT, "calibration_fit.csv"))
rows = []
for _, p in fit.iterrows():
    k_geo = p.k_bar * np.exp(p.sigma_scale * z_scale)
    a = p.sigma_asym * z_asym
    k = k_geo[:, None] * np.exp(a)
    rows.append(pd.DataFrame(dict(
        id=np.arange(N), species=p.species, key=p.key, variant=p.variant,
        k1=k[:, 0], k2=k[:, 1], k3=k[:, 2], k_geo=k_geo,
        asym=np.sqrt((a ** 2).sum(axis=1)), phase_offset_h=phase_offset,
        k_bar=p.k_bar, sigma_scale=p.sigma_scale, sigma_asym=p.sigma_asym, seed=SEED)))

pop = pd.concat(rows, ignore_index=True)
pop.to_csv(os.path.join(OUT, "population.csv"), index=False)

print(f"population.csv: {len(pop)} rows = {N} individuals x {pop.groupby(['species','variant']).ngroups} populations")
print(pop.groupby(["species", "variant"])[["k_geo", "asym"]].agg(["mean", "std"]).round(3).to_string())
