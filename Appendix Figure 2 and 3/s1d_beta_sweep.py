# =============================================================================
# s1d_beta_sweep.py   --   Simulation S1, step D
#
# For every simulated spider in population.csv, in this order:
#   0  DD equilibration, 10 d (oscillation check, note printed if absent)
#   1  DD assay at that species' observed DD record length
#   2  extra DD of U(0,24) h so lights-on meets each clock at a random phase
#   3  LD 12:12 transient, 15 d, discarded
#   4  LD assay at that species' observed LD record length   <- primary outcome
#   5  DD release, 10 d, reported separately
# Segments 0-2 are dark, where beta is inert, so they run once per individual and
# the end state is tiled across the whole beta grid. The beta comparison is thus
# paired down to the individual clock and its phase at lights-on.
#
# Measurements per individual per beta: amplitudes, tau by peak interval (exact)
# and by Lomb-Scargle (the spiders' pipeline), the two entrainment calls, phase
# and vector strength (activity-weighted, as in the Figure 1 script), light/dark
# ratio, and waveform shape (width, rise/fall, edge sharpness, damping, peaks/day).
# =============================================================================

# %% imports and paths
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import os, time

# astropy is used ONLY by test T6, which checks the fast periodogram below against
# astropy's. The simulation itself needs numpy/scipy/pandas/matplotlib only.
try:
    from astropy.timeseries import LombScargle
    HAVE_ASTROPY = True
except ImportError:
    HAVE_ASTROPY = False

try:
    OUT = os.path.dirname(os.path.abspath(__file__))
except NameError:
    OUT = os.getcwd()

# %% fixed constants, species table, protocol lengths
v, K, n_hill, gamma, L_amp, tau_m = 0.84, 1.0, 12, 12.0, 2.5, 2.0
LUX = 700.0
I_ON = LUX / 1000.0
DT = 0.01                 # integration step, h
STORE_EVERY = 10          # store every 0.1 h
STORE_DT = DT * STORE_EVERY

SPECIES = {  # published values, unchanged; DD/LD record lengths from Valid_Minutes/1440
    "L. cornutus":      dict(beta=0.09, alpha=0.90, L_base=0.00, tau_h=2.0, kappa=30.0, y_th=0.5, rec_dd=8.45, rec_ld=15.49),
    "A. pennsylvanica": dict(beta=0.06, alpha=0.10, L_base=0.20, tau_h=3.0, kappa=10.0, y_th=0.3, rec_dd=4.76, rec_ld=5.57),
    "S. grossa":        dict(beta=0.75, alpha=0.25, L_base=0.15, tau_h=3.0, kappa=10.0, y_th=0.3, rec_dd=7.46, rec_ld=6.36),
}
BETA_GRID = np.array([0.0, 0.03, 0.06, 0.09, 0.12, 0.20, 0.30, 0.50, 0.75, 1.0, 1.5, 2.0])
EQ_DAYS, TRANSIENT_DAYS, RELEASE_DAYS = 10.0, 15.0, 10.0
AMP_THRESHOLD = 0.4       # published activity-amplitude criterion
VARIANTS = os.environ.get("S1_VARIANTS", "fitted").split(",")
N_LIMIT = int(os.environ.get("S1_N_LIMIT", "0")) or None    # set for a quick test run


# %% model
def light(t, ld):
    """Light input I. ld is a boolean array (per element); lights on for ZT 0-12."""
    return np.where(ld & ((t % 24.0) < 12.0), I_ON, 0.0)


def deriv(s, t, p, ld):
    """Vectorised RHS. p holds per-element parameter arrays; ld marks LD elements."""
    x, y, z, h, M = s
    I = light(t, ld)
    k2_eff = p["k2"] * (1.0 + p["beta"] * I)
    dx = v * (K ** n_hill / (K ** n_hill + z ** n_hill)) - p["k1"] * x
    dy = v * x - k2_eff * y
    dz = v * y - p["k3"] * z
    h_inf = 1.0 / (1.0 + np.exp(np.clip(p["kappa"] * (y - p["y_th"]), -60, 60)))
    dh = (h_inf - h) / p["tau_h"]
    dM = ((I > 0).astype(float) - M) / tau_m
    return (dx, dy, dz, dh, dM)


def rk4_step(s, t, p, ld):
    d1 = deriv(s, t, p, ld)
    s2 = tuple(a + 0.5 * DT * b for a, b in zip(s, d1)); d2 = deriv(s2, t + 0.5 * DT, p, ld)
    s3 = tuple(a + 0.5 * DT * b for a, b in zip(s, d2)); d3 = deriv(s3, t + 0.5 * DT, p, ld)
    s4 = tuple(a + DT * b for a, b in zip(s, d3));       d4 = deriv(s4, t + DT, p, ld)
    return tuple(a + DT / 6.0 * (b1 + 2 * b2 + 2 * b3 + b4)
                 for a, b1, b2, b3, b4 in zip(s, d1, d2, d3, d4))


def activity(s, p):
    """Locomotor output L including baseline and masking."""
    x, y, z, h, M = s
    m_gate = 1.0 / (1.0 + np.exp(np.clip(-gamma * (y - p["y_th"]), -60, 60)))
    return (p["L_base"] + p["L_amp"] * m_gate * h) * (1.0 - p["alpha"] * M)


def integrate(s, t0, hours, p, ld, store=False):
    """Advance by `hours`. Returns final state, and (t, L, y) arrays if store."""
    n = int(round(hours / DT))
    n_store = n // STORE_EVERY
    L_out = np.empty((n_store, len(s[0])), np.float32) if store else None
    y_out = np.empty((n_store, len(s[0])), np.float32) if store else None
    t = t0
    for i in range(n):
        s = rk4_step(s, t, p, ld)
        t = t0 + (i + 1) * DT
        if store and (i + 1) % STORE_EVERY == 0:
            j = (i + 1) // STORE_EVERY - 1
            L_out[j] = activity(s, p)
            y_out[j] = s[1]
    if store:
        t_out = t0 + STORE_DT * np.arange(1, n_store + 1)
        return s, t_out, L_out, y_out
    return s


def params_for(pop_rows, sp, beta=None):
    """Per-element parameter dict for a population slice of one species."""
    c = SPECIES[sp]
    n = len(pop_rows)
    p = dict(k1=pop_rows.k1.values, k2=pop_rows.k2.values, k3=pop_rows.k3.values)
    for name in ["alpha", "L_base", "tau_h", "kappa", "y_th"]:
        p[name] = np.full(n, c[name])
    p["beta"] = np.full(n, c["beta"] if beta is None else beta)
    # output strength: published constant unless population.csv carries a per-individual L_amp
    p["L_amp"] = pop_rows.L_amp.values if "L_amp" in pop_rows else np.full(n, L_amp)
    return p


# %% analysis helpers (one individual at a time)
def peak_times(t, x, min_gap_h=2.0, prom_frac=0.2):
    """Times of prominent peaks, at least min_gap_h apart."""
    amp = x.max() - x.min()
    if amp < 1e-6:
        return np.array([])
    idx, _ = find_peaks(x, distance=int(min_gap_h / STORE_DT), prominence=prom_frac * amp)
    return t[idx]


LS_FREQS = np.linspace(1 / 35.0, 1 / 14.0, 2000)      # the spiders' grid: 14-35 h, 2000 points


def lomb_scargle_period(t, x):
    """The spiders' pipeline (astropy, floating mean): period at the periodogram peak."""
    if x.max() - x.min() < 1e-6:
        return np.nan
    power = LombScargle(t, x).power(LS_FREQS)
    return 1.0 / LS_FREQS[np.argmax(power)]


def batch_ls_period(t, X, chunk=1500):
    """Same floating-mean Lomb-Scargle as astropy's default, for every column of X at
    once. Exact least squares of [cos, sin, 1] at each frequency; the argmax is the
    period. Checked against lomb_scargle_period() in the test suite."""
    w = 2 * np.pi * LS_FREQS
    C = np.cos(np.outer(w, t)); S = np.sin(np.outer(w, t))          # (n_f, n_t)
    n_t = len(t)
    G = np.empty((len(w), 3, 3))
    G[:, 0, 0] = (C * C).sum(1); G[:, 1, 1] = (S * S).sum(1); G[:, 2, 2] = n_t
    G[:, 0, 1] = G[:, 1, 0] = (C * S).sum(1)
    G[:, 0, 2] = G[:, 2, 0] = C.sum(1); G[:, 1, 2] = G[:, 2, 1] = S.sum(1)
    Ginv = np.linalg.inv(G)
    out = np.full(X.shape[1], np.nan)
    for a in range(0, X.shape[1], chunk):
        Y = X[:, a:a + chunk].astype(np.float64)
        Y = Y - Y.mean(axis=0)
        b = np.stack([C @ Y, S @ Y, np.repeat(Y.sum(axis=0)[None, :], len(w), 0)], axis=1)  # (n_f, 3, n)
        red = np.einsum("fin,fij,fjn->fn", b, Ginv, b)                 # chi2 reduction
        flat = (X[:, a:a + chunk].max(0) - X[:, a:a + chunk].min(0)) < 1e-6
        idx = np.argmax(red, axis=0)
        out[a:a + chunk] = np.where(flat, np.nan, 1.0 / LS_FREQS[idx])
    return out


def circular_stats(t, x):
    """Activity-weighted mean phase (h, ZT) and vector strength, as in Figure 1."""
    if x.sum() <= 0:
        return np.nan, np.nan
    ang = 2 * np.pi * (t % 24.0) / 24.0
    c, s = (x * np.cos(ang)).sum(), (x * np.sin(ang)).sum()
    r = np.hypot(c, s) / x.sum()
    return (np.arctan2(s, c) % (2 * np.pi)) / (2 * np.pi) * 24.0, r


def daily_peak_locking(t, x):
    """Per-day argmax times -> (slope of peak time on day, Rayleigh r of peak phases)."""
    n_days = int(np.floor((t[-1] - t[0] + STORE_DT) / 24.0))
    if n_days < 3:
        return np.nan, np.nan
    tp = []
    for d in range(n_days):
        m = (t >= t[0] + 24 * d) & (t < t[0] + 24 * (d + 1))
        tp.append(t[m][np.argmax(x[m])])
    tp = np.array(tp)
    slope = np.polyfit(np.arange(n_days), tp, 1)[0]
    ang = 2 * np.pi * (tp % 24.0) / 24.0
    r = np.hypot(np.cos(ang).sum(), np.sin(ang).sum()) / n_days
    return slope, r


def fold_profile(t, x, period, n_bins=24, align_peak=False):
    """Mean waveform folded at `period` into n_bins; optionally peak at bin n_bins//2."""
    # tiny offset keeps samples that sit exactly on a bin edge from being split by float error
    phase = ((t - t[0] + 1e-6) % period) / period
    b = np.minimum((phase * n_bins).astype(int), n_bins - 1)
    prof = np.bincount(b, weights=x, minlength=n_bins) / np.maximum(np.bincount(b, minlength=n_bins), 1)
    if align_peak:
        prof = np.roll(prof, n_bins // 2 - int(np.argmax(prof)))
    return prof


def shape_metrics(p):
    """Shape statistics of one 24-bin circular profile (same code as spider reference)."""
    p = np.asarray(p, float)
    lo, hi = p.min(), p.max(); amp = hi - lo
    if amp < 1e-6 or p.sum() <= 0:
        return dict(fwhm_h=np.nan, peak6_frac=np.nan, max_jump=np.nan, peak_to_mean=np.nan, n_peaks=0)
    fwhm = int((p > lo + 0.5 * amp).sum())
    peak6 = max(np.roll(p, -s)[:6].sum() for s in range(24)) / p.sum()
    jump = np.abs(np.diff(np.r_[p, p[0]])).max() / amp
    left, right = np.roll(p, 1), np.roll(p, -1)
    npk = int(((p > left) & (p >= right) & (p - lo > 0.2 * amp)).sum())
    return dict(fwhm_h=fwhm, peak6_frac=peak6, max_jump=jump, peak_to_mean=hi / p.mean(), n_peaks=npk)


def fine_shape(t, x, period):
    """Rise/fall times, edge sharpness and damping from the 0.1 h trace."""
    w = fold_profile(t, x, period, n_bins=int(round(period / STORE_DT)), align_peak=True)
    lo, hi = w.min(), w.max(); amp = hi - lo
    if amp < 1e-6:
        return dict(rise_h=np.nan, fall_h=np.nan, edge_sharp=np.nan, fwhm_fine_h=np.nan, damping=np.nan)
    ipk = int(np.argmax(w))
    above10 = np.where(w > lo + 0.1 * amp)[0]; above90 = np.where(w > lo + 0.9 * amp)[0]
    rise = (above90[above90 <= ipk].min() - above10[above10 <= ipk].min()) * STORE_DT if (above90 <= ipk).any() and (above10 <= ipk).any() else np.nan
    fall = (above10[above10 >= ipk].max() - above90[above90 >= ipk].max()) * STORE_DT if (above90 >= ipk).any() and (above10 >= ipk).any() else np.nan
    edge = np.abs(np.diff(np.r_[w, w[0]])).max() / amp            # largest change per 0.1 h, in amplitudes
    fwhm = (w > lo + 0.5 * amp).sum() * STORE_DT
    n_days = int(np.floor((t[-1] - t[0]) / 24.0))
    if n_days >= 2:
        daily = np.array([np.ptp(x[(t >= t[0] + 24 * d) & (t < t[0] + 24 * (d + 1))]) for d in range(n_days)])
        damping = daily[-1] / max(daily[0], 1e-9)
        amp_cv = daily.std() / max(daily.mean(), 1e-9)
    else:
        damping, amp_cv = np.nan, np.nan
    return dict(rise_h=rise, fall_h=fall, edge_sharp=edge, fwhm_fine_h=fwhm, damping=damping, amp_cv=amp_cv)


# %% built-in tests (must all pass before anything else is believed)
def test_suite(pop):
    print("running built-in tests")
    sub = pop[(pop.species == "S. grossa") & (pop.variant == pop.variant.iloc[0])].head(40)
    # T3: k = 0.20 symmetric gives tau_DD = 20.90 h
    one = sub.head(1).copy(); one[["k1", "k2", "k3"]] = 0.20
    p = params_for(one, "S. grossa")
    s = tuple(np.full(1, val) for val in (0.5, 0.5, 0.5, 0.5, 0.0))
    s = integrate(s, 0.0, 240.0, p, np.zeros(1, bool))
    s, t, L, y = integrate(s, 0.0, 240.0, p, np.zeros(1, bool), store=True)
    tau = np.diff(peak_times(t, y[:, 0])).mean()
    print(f"  T3  tau at k=0.20: {tau:.2f} h (expect 20.90)")
    assert abs(tau - 20.90) < 0.05
    # T1: beta is inert in DD (bitwise)
    p0, p2 = params_for(sub, "S. grossa", 0.0), params_for(sub, "S. grossa", 2.0)
    s = tuple(np.full(len(sub), val) for val in (0.5, 0.5, 0.5, 0.5, 0.0))
    a = integrate(s, 0.0, 100.0, p0, np.zeros(len(sub), bool))
    b = integrate(s, 0.0, 100.0, p2, np.zeros(len(sub), bool))
    print(f"  T1  max |state(beta=0) - state(beta=2)| in DD: {max(np.abs(u - w).max() for u, w in zip(a, b)):.1e} (expect 0)")
    assert all(np.array_equal(u, w) for u, w in zip(a, b))
    # T2: A. pennsylvanica and S. grossa identical in DD for the same k draws
    pa, ps = params_for(sub, "A. pennsylvanica"), params_for(sub, "S. grossa")
    a = integrate(s, 0.0, 100.0, pa, np.zeros(len(sub), bool)); b = integrate(s, 0.0, 100.0, ps, np.zeros(len(sub), bool))
    core_same = all(np.array_equal(a[i], b[i]) for i in range(4))
    print(f"  T2  A.p and S.g share core and gate in DD: {core_same} (expect True; L differs only by L_base)")
    assert core_same
    # T5: at beta = 0 the clock in LD equals the clock in DD (light acts through alpha only)
    p0 = params_for(sub, "L. cornutus", 0.0)
    a = integrate(s, 0.0, 100.0, p0, np.ones(len(sub), bool)); b = integrate(s, 0.0, 100.0, p0, np.zeros(len(sub), bool))
    print(f"  T5  max |y_LD - y_DD| at beta=0: {np.abs(a[1] - b[1]).max():.1e} (expect 0)")
    assert np.array_equal(a[1], b[1])
    # T6: batch Lomb-Scargle reproduces astropy's peak period
    if not HAVE_ASTROPY:
        print("  T6  SKIPPED - astropy not installed (pip install astropy). The simulation\n"
              "      does not need it; only this cross-check does.")
        print("  all other tests passed\n")
        return
    tt = 360 + STORE_DT * np.arange(1, 1500)
    X = np.stack([(1 + np.sin(2 * np.pi * tt / per)) ** 3 + 0.2 * np.cos(2 * np.pi * tt / 24) for per in (20.5, 23.1, 24.0, 27.7)], 1)
    ref = np.array([lomb_scargle_period(tt, X[:, i]) for i in range(4)])
    got = batch_ls_period(tt, X)
    print(f"  T6  batch vs astropy Lomb-Scargle, max |dperiod| = {np.abs(ref - got).max():.4f} h (expect < 0.02)")
    assert np.abs(ref - got).max() < 0.02
    print("  all tests passed\n")


# %% main protocol for one species x variant
def run_population(pop_sp, sp, variant):
    c = SPECIES[sp]; n = len(pop_sp); p = params_for(pop_sp, sp)
    dark = np.zeros(n, bool)
    t_start = time.time()

    # 0. equilibration
    s = tuple(np.full(n, val) for val in (0.5, 0.5, 0.5, 0.5, 0.0))
    s = integrate(s, 0.0, (EQ_DAYS - 3) * 24, p, dark)
    s, t_eq, L_eq, y_eq = integrate(s, 0.0, 3 * 24, p, dark, store=True)
    n_flat = sum(len(peak_times(t_eq, y_eq[:, i])) < 2 for i in range(n))
    print(f"  equilibration: {n_flat} of {n} clocks do not oscillate; their end state is used anyway")

    # 1. DD assay
    s, t_dd, L_dd, y_dd = integrate(s, 0.0, c["rec_dd"] * 24, p, dark, store=True)

    # 2. phase offset: capture each individual's state at its own offset time
    off_idx = np.round(pop_sp.phase_offset_h.values / DT).astype(int)
    start = [a.copy() for a in s]
    t = 0.0
    for i in range(1, int(24 / DT) + 1):
        s = rk4_step(s, t, p, dark); t = i * DT
        m = off_idx == i
        for a, b in zip(start, s):
            a[m] = b[m]
    start = tuple(start)

    # 3-5. tile across beta, LD transient, LD assay, release
    nb = len(BETA_GRID)
    big = tuple(np.tile(a, nb) for a in start)
    pb = {k: np.tile(a, nb) for k, a in p.items()}
    pb["beta"] = np.repeat(BETA_GRID, n)
    lit = np.ones(n * nb, bool)
    big = integrate(big, 0.0, TRANSIENT_DAYS * 24, pb, lit)
    t0 = TRANSIENT_DAYS * 24
    big, t_ld, L_ld, y_ld = integrate(big, t0, c["rec_ld"] * 24, pb, lit, store=True)
    t1 = t_ld[-1]
    big, t_rel, L_rel, y_rel = integrate(big, t1, RELEASE_DAYS * 24, pb, np.zeros(n * nb, bool), store=True)
    print(f"  integration done in {time.time() - t_start:.0f} s; Lomb-Scargle in batch")
    tau_ls_dd = batch_ls_period(t_dd, L_dd)
    tau_ls_ld = batch_ls_period(t_ld, L_ld)
    print(f"  periodograms done at {time.time() - t_start:.0f} s; per-individual analysis")

    rows, profiles = [], []
    dd_cache = {}
    milestones = {max(1, round(n * f)) for f in (0.25, 0.5, 0.75)}
    for j in range(n * nb):
        i, beta = j % n, BETA_GRID[j // n]
        if i == 0:
            print(f"  beta = {beta:g}: testing {n} oscillators")
        if i not in dd_cache:
            dd_cache[i] = analyse_dd(t_dd, L_dd[:, i], y_dd[:, i], tau_ls_dd[i])
        d = dict(species=sp, variant=variant, id=int(pop_sp.id.iloc[i]), beta=beta,
                 beta_is_published=bool(np.isclose(beta, c["beta"])),
                 k1=pop_sp.k1.iloc[i], k2=pop_sp.k2.iloc[i], k3=pop_sp.k3.iloc[i],
                 k_geo=pop_sp.k_geo.iloc[i], asym=pop_sp.asym.iloc[i])
        d.update(dd_cache[i])
        d.update(analyse_ld(t_ld, L_ld[:, j], y_ld[:, j], tau_ls_ld[j]))
        d.update(analyse_release(t_rel, L_rel[:, j]))
        rows.append(d)
        profiles.append(dict(species=sp, variant=variant, id=d["id"], beta=beta, condition="LD",
                             **{f"h{k}": val for k, val in enumerate(fold_profile(t_ld, L_ld[:, j], 24.0))}))
        if j < n:
            per = dd_cache[i]["tau_true_DD"] if np.isfinite(dd_cache[i]["tau_true_DD"]) else 24.0
            profiles.append(dict(species=sp, variant=variant, id=d["id"], beta=np.nan, condition="DD",
                                 **{f"h{k}": val for k, val in enumerate(fold_profile(t_dd, L_dd[:, i], per, align_peak=True))}))
        if (i + 1) in milestones:
            print(f"    beta = {beta:g}: {round(100 * (i + 1) / n)}% done ({i + 1} of {n} oscillators)")
    print(f"  {sp} / {variant} finished in {time.time() - t_start:.0f} s")
    return rows, profiles


def analyse_dd(t, L, y, tau_ls):
    pk_y, pk_L = peak_times(t, y), peak_times(t, L)
    A = L.max() - L.mean()
    tau = np.diff(pk_L).mean() if len(pk_L) >= 2 else np.nan
    tau_prot = np.diff(pk_y).mean() if len(pk_y) >= 2 else np.nan
    # exact period = protein peak interval (one clean peak per cycle); the L-peak interval
    # is kept as a diagnostic because some gates give a shoulder peak within a cycle
    d = dict(A_DD=A, A_prot_DD=y.max() - y.mean(), n_peaks_prot_DD=len(pk_y),
             tau_true_DD=tau_prot, tau_prot_DD=tau_prot, tau_Lpeaks_DD=tau,
             tau_LS_DD=tau_ls,
             rhythmic_DD=bool(len(pk_y) >= 2 and A > AMP_THRESHOLD))
    per = tau_prot if np.isfinite(tau_prot) else 24.0
    d.update({f"{k}_DD": val for k, val in shape_metrics(fold_profile(t, L, per, align_peak=True)).items()})
    d.update({f"{k}_DD": val for k, val in fine_shape(t, L, per).items()})
    return d


def analyse_ld(t, L, y, tau_ls):
    pk_L, pk_y = peak_times(t, L), peak_times(t, y)
    A = L.max() - L.mean()
    tau = np.diff(pk_L).mean() if len(pk_L) >= 2 else np.nan
    tau_prot = np.diff(pk_y).mean() if len(pk_y) >= 2 else np.nan
    slope, r_lock = daily_peak_locking(t, L)
    phase, vec = circular_stats(t, L)
    zt = t % 24.0
    d = dict(A_LD=A, tau_true_LD=tau_prot, tau_Lpeaks_LD=tau, tau_LS_LD=tau_ls, peak_slope_LD=slope, lock_r_LD=r_lock,
             entrained_LS=bool(np.isfinite(tau_ls) and 22.8 <= tau_ls <= 25.2 and A > AMP_THRESHOLD),
             entrained_true=bool(np.isfinite(slope) and abs(slope - 24) < 0.05 and r_lock > 0.9 and A > AMP_THRESHOLD),
             phase_ZT_LD=phase, vector_r_LD=vec,
             light_dark_ratio=L[zt < 12].mean() / max(L[zt >= 12].mean(), 1e-9),
             peaks_per_day_LD=len(pk_L) / ((t[-1] - t[0]) / 24.0))
    d.update({f"{k}_LD": val for k, val in shape_metrics(fold_profile(t, L, 24.0)).items()})
    d.update({f"{k}_LD": val for k, val in fine_shape(t, L, 24.0).items()})
    return d


def analyse_release(t, L):
    pk = peak_times(t, L)
    return dict(A_release=L.max() - L.mean(),
                tau_release=np.diff(pk).mean() if len(pk) >= 2 else np.nan,
                first_peak_release_h=(pk[0] - t[0]) if len(pk) else np.nan)


# %% run everything
if __name__ == "__main__":
    pop = pd.read_csv(os.path.join(OUT, "population.csv"))
    if N_LIMIT:
        pop = pop[pop.id < N_LIMIT]
        print(f"TEST MODE: {N_LIMIT} individuals per population")
    test_suite(pop)

    all_rows, all_profiles = [], []
    for variant in VARIANTS:
        for sp in SPECIES:
            sel = pop[(pop.species == sp) & (pop.variant == variant)].reset_index(drop=True)
            print(f"{sp} / {variant}: {len(sel)} individuals x {len(BETA_GRID)} betas")
            rows, profs = run_population(sel, sp, variant)
            all_rows += rows; all_profiles += profs
            pd.DataFrame(all_rows).to_csv(os.path.join(OUT, "results_individual.csv"), index=False)
            pd.DataFrame(all_profiles).to_csv(os.path.join(OUT, "profiles_hourly.csv"), index=False)

    res = pd.DataFrame(all_rows)
    ent = res[res.entrained_LS]
    summary = res.groupby(["species", "variant", "beta"]).agg(
        n=("id", "size"), pct_rhythmic_DD=("rhythmic_DD", lambda s: 100 * s.mean()),
        pct_entrained_LS=("entrained_LS", lambda s: 100 * s.mean()),
        pct_entrained_true=("entrained_true", lambda s: 100 * s.mean()),
        tau_DD_mean=("tau_true_DD", "mean"), tau_DD_sd=("tau_true_DD", "std"),
        tau_LS_LD_mean=("tau_LS_LD", "mean"), tau_LS_LD_sd=("tau_LS_LD", "std"),
        A_DD=("A_DD", "mean"), A_LD=("A_LD", "mean"), phase_ZT=("phase_ZT_LD", "mean"),
        vector_r=("vector_r_LD", "mean"), light_dark=("light_dark_ratio", "mean"),
        fwhm_LD=("fwhm_h_LD", "mean"), peak6_LD=("peak6_frac_LD", "mean"),
        max_jump_LD=("max_jump_LD", "mean"), edge_sharp_LD=("edge_sharp_LD", "mean"),
        damping_LD=("damping_LD", "mean"), peaks_per_day=("peaks_per_day_LD", "mean")).reset_index()
    ent_sd = ent.groupby(["species", "variant", "beta"]).tau_LS_LD.std().rename("tau_LD_sd_entrained").reset_index()
    summary = summary.merge(ent_sd, how="left")
    summary.to_csv(os.path.join(OUT, "results_summary.csv"), index=False)
    print("\nsummary at published beta:")
    print(summary[summary.apply(lambda r: np.isclose(r.beta, SPECIES[r.species]["beta"]), axis=1)]
          .round(2).to_string(index=False))
