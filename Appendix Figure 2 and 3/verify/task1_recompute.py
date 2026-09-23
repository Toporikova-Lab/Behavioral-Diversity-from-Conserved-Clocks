"""Task 1: recompute the four claims directly from results_individual.csv."""
import numpy as np, pandas as pd
CSV = "/home/claude/S1/track1_main/results_individual.csv"
df = pd.read_csv(CSV)
df = df[df.variant == "fitted"].copy()
PUB = {"L. cornutus": 0.09, "A. pennsylvanica": 0.06, "S. grossa": 0.75}
print("rows (fitted):", len(df), "| n per species x beta:")
print(df.groupby(["species", "beta"]).size().unstack().to_string())

# Claim 4 first: definition check
recomputed = (df.tau_LS_LD.between(22.8, 25.2) & (df.A_LD > 0.4))
mism = (recomputed != df.entrained_LS.astype(bool)).sum()
print(f"\nClaim 4  entrained_LS == (22.8<=tau_LS_LD<=25.2 & A_LD>0.4): mismatches = {mism} of {len(df)}")

# Claim 1
print("\nClaim 1  % entrained_LS at published beta:")
for sp, b in PUB.items():
    sub = df[(df.species == sp) & np.isclose(df.beta, b)]
    print(f"  {sp:18s} beta={b:5.2f} n={len(sub)}  pct={100*sub.entrained_LS.mean():.1f}  (flag beta_is_published all True: {sub.beta_is_published.all()})")

# Claim 2
sg = df[df.species == "S. grossa"].groupby("beta").entrained_LS.mean().mul(100)
print("\nClaim 2  S. grossa % entrained_LS by beta:")
print(sg.round(1).to_string())
print("  monotonic non-decreasing:", bool(np.all(np.diff(sg.values) >= 0)),
      "| strictly increasing:", bool(np.all(np.diff(sg.values) > 0)),
      f"| beta=0: {sg.iloc[0]:.1f}, beta=2: {sg.iloc[-1]:.1f}")

# Claim 3
print("\nClaim 3  % rhythmic_DD by species and beta:")
r = df.groupby(["species", "beta"]).rhythmic_DD.mean().mul(100).unstack()
print(r.round(2).to_string())
print("  identical across betas per species:", {sp: bool(r.loc[sp].nunique() == 1) for sp in r.index})
# stronger check: DD columns identical per id across betas
ddcols = [c for c in df.columns if c.endswith("_DD")]
nun = df.groupby(["species", "id"])[ddcols].nunique().max()
print("  max distinct values of any *_DD column within an individual across betas:", int(nun.max()))
