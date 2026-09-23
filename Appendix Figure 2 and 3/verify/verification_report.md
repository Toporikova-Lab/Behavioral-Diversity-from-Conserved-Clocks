# Independent verification of S1 beta sweep (track1_main), 2026-09-10

Inputs read (not modified): `/home/claude/S1/track1_main/results_individual.csv`,
`/home/claude/S1/track1_main/population.csv`, `/home/claude/S1/track1_main/s1d_beta_sweep.py` (read for protocol only).
Scripts and outputs in `/home/claude/S1/verify/`: `task1_recompute.py` (+ `task1_output.txt`),
`task2_resim_odeint.py` (+ `task2_output.txt`, `task2_resim_results.csv`).

## Task 1 - claims recomputed from results_individual.csv (variant = fitted, 500 ids x 12 betas x 3 species = 18000 rows)

| Claim | Recomputed | Verdict |
|---|---|---|
| 1. % entrained_LS at published beta | L. cornutus (0.09) **64.0**; A. pennsylvanica (0.06) **21.4**; S. grossa (0.75) **60.4** (n = 500 each; `beta_is_published` flag consistent) | AGREE, exact |
| 2. S. grossa % entrained_LS rises monotonically with beta | 16.2, 16.6, 17.0, 21.4, 24.4, 37.4, 44.4, 55.0, 60.4, 63.4, 68.8, 71.8 for beta = 0, .03, .06, .09, .12, .2, .3, .5, .75, 1, 1.5, 2. Strictly increasing at every step; beta=0 -> 16.2, beta=2 -> 71.8 | AGREE, exact |
| 3. % rhythmic_DD identical across betas | A. pennsylvanica 95.8, L. cornutus 88.2, S. grossa 69.2 at all 12 betas; every `*_DD` column has exactly 1 distinct value per individual across betas | AGREE |
| 4. entrained_LS == (22.8 <= tau_LS_LD <= 25.2) & (A_LD > 0.4) | 0 mismatches in 18000 rows | AGREE |

Caveat on claim 3: it holds by construction, not as an independent measurement. The script integrates
segments 0-2 (all dark) once per individual and tiles the end state across the beta grid, so DD columns
are copied, not recomputed per beta. That is legitimate because beta enters the model only through
`k2_eff = k2*(1+beta*I)` and I = 0 in DD, so beta is inert in DD analytically (the script's test T1 also
checks this bitwise). The claim is true but should be described as a property of the model, not a finding.

## Task 2 - independent re-simulation with scipy.integrate.odeint (S. grossa, beta = 0.75)

Method: model written from the equations given (v=0.84, K=1, n=12, gamma=12, L_amp=2.5, tau_m=2, I=0.7;
S. grossa alpha=0.25, L_base=0.15, tau_h=3, kappa=10, y_th=0.3), individual k1,k2,k3 from the CSV and
phase_offset_h from population.csv. Protocol: 10 d DD from [0.5,0.5,0.5,0.5,0] -> 7.46 d DD assay ->
extra DD of phase_offset_h (rounded to 0.01 h) -> LD 12:12 (lights on when t mod 24 < 12, t=0 at LD
onset) 15 d transient -> 6.36 d LD assay. Because the light input is discontinuous, I integrated
piecewise per constant-light segment (12 h blocks), carrying the state across boundaries; inside each
segment odeint ran on a 0.01 h output grid with hmax=0.01, rtol=1e-10, atol=1e-12. Lomb-Scargle:
astropy, frequencies linspace(1/35, 1/14, 2000), period at argmax, on the 0.1 h subsample the pipeline
uses (also on the full 0.01 h grid - identical result). A_LD = max(L) - mean(L).

| id | case | tau_true_DD csv / mine | tau_LS_LD csv / mine | A_LD csv / mine | entrained_LS csv / mine |
|---|---|---|---|---|---|
| 0 | entrained True | 21.700 / 21.699 | **23.936 / 23.936** | 0.8656 / 0.8658 | True / True |
| 2 | False, tau_DD < 22 | 18.475 / 18.468 | **16.708 / 16.708** | 0.6808 / 0.6812 | False / False |
| 113 | False, tau_DD > 30 | 31.875 / 31.883 | **15.077 / 15.077** | 0.7921 / 0.7921 | False / False |
| 201 | bonus: False, tau_DD > 30 | 34.125 / 34.148 | **30.896 / 30.896** | 0.7738 / 0.7736 | False / False |

tau_LS_DD (21.718, 18.497, 32.043, 34.204) and A_DD also agree to 4 decimals; tau_true_LD
(protein-peak interval) agrees to <0.01 h. All four entrained_LS calls reproduce. Residual differences
(<0.01 h in peak-interval periods, <0.0005 in amplitudes) are consistent with RK4 dt=0.01 vs adaptive
odeint and with my A_LD being taken on the 0.01 h grid. The exact agreement in tau_LS_* reflects the
discrete 2000-point frequency grid: both runs land on the same grid point.

Side check: population.csv carries no `L_amp` column, so the script's per-individual L_amp branch is not
taken and the constant 2.5 applies - confirmed by the amplitude agreement above.

## Concerns / notes (none invalidate the claims)

1. Claim 3 is tautological given the tiling design (see caveat above); wording in the manuscript should
   reflect that beta cannot affect DD by construction of the model.
2. For S. grossa the LD assay is only 6.36 d (152.6 h). The Lomb-Scargle period resolution near 24 h is
   roughly 24^2/152.6 ~ 3.8 h, wider than the 2.4 h entrainment window [22.8, 25.2]. Peak location is
   finer than resolution, so the criterion is usable, but tau_LS_LD values near the window edges are not
   sharply determined. This is inherited from the spiders' pipeline, so it is consistent with the data
   analysis, just worth keeping in mind (e.g. the entrained id 0 reads 23.94 h rather than 24.00 h).
3. Non-entrained long-tau individuals (ids 113, 201) show LS peaks at ~15 h and ~31 h respectively,
   i.e. the periodogram picks a sub-harmonic or the free-running period from a relative-coordination
   pattern; either way the "not entrained" call is correct and reproduced.
4. No discrepancy found between results_individual.csv, results_summary.csv/run.log, and the
   independent simulation.

**Overall: all four claims verified; independent odeint re-simulation reproduces the CSV for 4/4
individuals, including the entrainment calls.**
