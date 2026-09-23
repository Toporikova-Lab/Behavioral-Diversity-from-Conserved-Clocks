"""Task 2: independent re-simulation of selected S. grossa individuals with scipy odeint.
Written from the model equations only (no code reused from s1d_beta_sweep.py).
Light is piecewise constant, so the integration is done segment by segment (each
segment has constant I); inside a segment odeint runs on a 0.01 h output grid with
hmax=0.01 and tight tolerances. State is carried across segment boundaries."""
import numpy as np, pandas as pd
from scipy.integrate import odeint
from scipy.signal import find_peaks
from astropy.timeseries import LombScargle

RES = "/home/claude/S1/track1_main/results_individual.csv"
POP = "/home/claude/S1/track1_main/population.csv"
v, K, n, gamma, L_amp, tau_m = 0.84, 1.0, 12, 12.0, 2.5, 2.0
I_ON = 0.7
SG = dict(alpha=0.25, L_base=0.15, tau_h=3.0, kappa=10.0, y_th=0.3)
DT = 0.01
REC_DD, REC_LD, EQ_D, TRANS_D = 7.46, 6.36, 10.0, 15.0
FREQS = np.linspace(1 / 35.0, 1 / 14.0, 2000)


def rhs(s, t, k1, k2, k3, beta, I):
    x, y, z, h, M = s
    k2_eff = k2 * (1 + beta * I)
    dx = v * K**n / (K**n + z**n) - k1 * x
    dy = v * x - k2_eff * y
    dz = v * y - k3 * z
    h_inf = 1.0 / (1.0 + np.exp(SG["kappa"] * (y - SG["y_th"])))
    dh = (h_inf - h) / SG["tau_h"]
    dM = ((1.0 if I > 0 else 0.0) - M) / tau_m
    return [dx, dy, dz, dh, dM]


def activity(S):
    x, y, z, h, M = S.T
    m_gate = 1.0 / (1.0 + np.exp(-gamma * (y - SG["y_th"])))
    return (SG["L_base"] + L_amp * m_gate * h) * (1 - SG["alpha"] * M)


def run_segment(s0, t0, hours, I, k):
    """Integrate `hours` with constant light I. Returns (t, states) excluding t0 sample."""
    npts = int(round(hours / DT))
    t = t0 + DT * np.arange(npts + 1)
    S = odeint(rhs, s0, t, args=(*k, I), rtol=1e-10, atol=1e-12, hmax=DT)
    return t[1:], S[1:]


def run_dd(s0, t0, hours, k):
    return run_segment(s0, t0, hours, 0.0, k)


def run_ld(s0, t0, hours, k):
    """LD 12:12 with lights on when (t mod 24) < 12, t measured from LD onset (t0=0)."""
    ts, Ss, s, t = [], [], np.array(s0), t0
    end = t0 + hours
    while t < end - 1e-9:
        zt = t % 24.0
        seg_end = t + (12.0 - zt if zt < 12.0 else 24.0 - zt)
        seg_end = min(seg_end, end)
        I = I_ON if zt < 12.0 else 0.0
        tt, S = run_segment(s, t, seg_end - t, I, k)
        ts.append(tt); Ss.append(S); s = S[-1]; t = seg_end
    return np.concatenate(ts), np.concatenate(Ss)


def ls_period(t, x):
    p = LombScargle(t, x).power(FREQS)
    return 1.0 / FREQS[np.argmax(p)]


def peak_tau(t, x):
    amp = x.max() - x.min()
    idx, _ = find_peaks(x, distance=int(2.0 / DT), prominence=0.2 * amp)
    return np.diff(t[idx]).mean() if len(idx) >= 2 else np.nan


def simulate(k1, k2, k3, beta, offset_h):
    k = (k1, k2, k3, beta)
    s = [0.5, 0.5, 0.5, 0.5, 0.0]
    _, S = run_dd(s, 0.0, EQ_D * 24, k)                       # 10 d equilibration
    t_dd, S_dd = run_dd(S[-1], 0.0, REC_DD * 24, k)           # DD assay 7.46 d
    _, S = run_dd(S_dd[-1], 0.0, round(offset_h / DT) * DT, k)  # extra DD, phase offset
    _, S = run_ld(S[-1], 0.0, TRANS_D * 24, k)               # LD transient 15 d
    t_ld, S_ld = run_ld(S[-1], TRANS_D * 24, REC_LD * 24, k)  # LD assay 6.36 d
    L_dd, L_ld = activity(S_dd), activity(S_ld)
    # analyse on the 0.1 h subsample used by the pipeline, and on the full 0.01 h grid
    sub = slice(9, None, 10)
    return dict(tau_true_DD=peak_tau(t_dd, S_dd[:, 1]),
                A_DD=L_dd.max() - L_dd.mean(),
                tau_LS_DD=ls_period(t_dd[sub], L_dd[sub]),
                A_LD=L_ld.max() - L_ld.mean(),
                tau_LS_LD=ls_period(t_ld[sub], L_ld[sub]),
                tau_LS_LD_fullgrid=ls_period(t_ld, L_ld),
                tau_true_LD=peak_tau(t_ld, S_ld[:, 1]))


if __name__ == "__main__":
    res = pd.read_csv(RES)
    pop = pd.read_csv(POP)
    sg = res[(res.species == "S. grossa") & (res.variant == "fitted") & np.isclose(res.beta, 0.75)].set_index("id")
    off = pop[(pop.species == "S. grossa") & (pop.variant == "fitted")].set_index("id").phase_offset_h
    IDS = {0: "entrained_LS True", 2: "False, tau_DD<22", 113: "False, tau_DD>30", 201: "bonus: False, tau_DD>30, tau_LS_LD~31"}
    rows = []
    for i, lab in IDS.items():
        r = sg.loc[i]
        mine = simulate(r.k1, r.k2, r.k3, float(r.beta), off.loc[i])
        ent = (22.8 <= mine["tau_LS_LD"] <= 25.2) and mine["A_LD"] > 0.4
        rows.append(dict(id=i, case=lab, k1=r.k1, k2=r.k2, k3=r.k3, phase_offset_h=off.loc[i],
                         csv_tau_true_DD=r.tau_true_DD, my_tau_true_DD=mine["tau_true_DD"],
                         csv_A_DD=r.A_DD, my_A_DD=mine["A_DD"],
                         csv_tau_LS_DD=r.tau_LS_DD, my_tau_LS_DD=mine["tau_LS_DD"],
                         csv_tau_LS_LD=r.tau_LS_LD, my_tau_LS_LD=mine["tau_LS_LD"], my_tau_LS_LD_fullgrid=mine["tau_LS_LD_fullgrid"],
                         csv_A_LD=r.A_LD, my_A_LD=mine["A_LD"],
                         csv_tau_true_LD=r.tau_true_LD, my_tau_true_LD=mine["tau_true_LD"],
                         csv_entrained_LS=bool(r.entrained_LS), my_entrained_LS=ent, match=bool(r.entrained_LS) == ent))
        print(f"id {i:4d} [{lab}] done")
    out = pd.DataFrame(rows)
    out.to_csv("/home/claude/S1/verify/task2_resim_results.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.round(4).T.to_string())
