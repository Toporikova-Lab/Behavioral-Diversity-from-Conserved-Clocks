# =============================================================================
# s1h_methods_figures.py  --  two extra figures for the collaborator document
#   M1  protocol schematic: what each oscillator is run through, in order
#   M2  calibration check: simulated vs observed free-running behaviour in DD
# beta is plotted on the MANUSCRIPT scale (beta_ms = 0.7 * beta_simulation).
# =============================================================================

# %% inputs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
import os

OUT = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
FIG = os.path.join(OUT, "supplement_figures"); os.makedirs(FIG, exist_ok=True)
I_SIM = 0.7


def track_file(track, name):
    p = os.path.join(OUT, track, name + ".gz")
    return p if os.path.exists(p) else os.path.join(OUT, track, name)


SP = ["L. cornutus", "A. pennsylvanica", "S. grossa"]
COL = {"L. cornutus": "#0072B2", "A. pennsylvanica": "#E69F00", "S. grossa": "#009E73"}
MK = {"L. cornutus": "o", "A. pennsylvanica": "s", "S. grossa": "^"}
PUB_MS = {"L. cornutus": 0.09, "A. pennsylvanica": 0.06, "S. grossa": 0.75}
REC = {"L. cornutus": (8.45, 15.49), "A. pennsylvanica": (4.76, 5.57), "S. grossa": (7.46, 6.36)}
COMP = {"L. cornutus": "LC_spider_analysis_comprehensive_with_LD_split.csv",
        "A. pennsylvanica": "Ag_spider_analysis_comprehensive_with_LD_split.csv",
        "S. grossa": "Sg_spider_analysis_comprehensive_with_LD_split.csv"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

res1 = pd.read_csv(track_file("track1_main", "results_individual.csv"))

# %% Figure M0 -- pipeline diagram: what each stage consumes and produces
INK, MUTED, ACC, ACCFC = "#1a1a1a", "#5f5f5f", "#0072B2", "#EAF3FA"
TITLE_H, LINE_H, PAD = 0.055, 0.034, 0.022


def box(ax, xc, ytop, w, title, lines, fc="#FCFCFA", ec="0.55", tc=INK, mono=False):
    """Box centred on xc with its top at ytop, height set by the number of body lines."""
    h = TITLE_H + LINE_H * len(lines) + PAD
    x, y = xc - w / 2, ytop - h
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.018",
                                facecolor=fc, edgecolor=ec, linewidth=1.0, zorder=2))
    ax.text(xc, ytop - 0.030, title, ha="center", va="center", fontsize=8.6, weight="bold",
            color=tc, zorder=3, family="monospace" if mono else None)
    for i, ln in enumerate(lines):
        ax.text(xc, ytop - TITLE_H - LINE_H * (i + 0.5) + 0.004, ln, ha="center", va="center",
                fontsize=7.5, color=MUTED, zorder=3)
    return dict(xc=xc, top=ytop, bot=y, left=x, right=x + w)


def arrow(ax, p0, p1, color="0.45", ls="-", rad=0.0, lw=1.1):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=11, color=color,
                                 linewidth=lw, linestyle=ls, connectionstyle=f"arc3,rad={rad}",
                                 zorder=1, shrinkA=2, shrinkB=2))


fig, ax = plt.subplots(figsize=(12.6, 8.3))
ax.set_xlim(0, 1); ax.set_ylim(-0.11, 1.01); ax.axis("off")

ax.text(0.012, 0.985, "Pipeline: from the recorded spiders to the figures",
        fontsize=12, weight="bold", va="top", color=INK)
ax.text(0.012, 0.945, "Blue: measured data, never regenerated.  Grey: computed by the script named in the box.  "
                      "Every seed is fixed, so a rerun reproduces each number exactly.",
        fontsize=8.2, va="top", color=MUTED)

# --- row 1: what goes in
IN_Y = 0.895
b_data = box(ax, 0.135, IN_Y, 0.245, "Measured behaviour",
             ["3 period files, one row per spider", "per condition: period, p, activity,", "record length"],
             fc=ACCFC, ec=ACC, tc=ACC)
b_circ = box(ax, 0.415, IN_Y, 0.245, "Measured entrainment",
             ["3 LD circular-statistics files:", "period, phase, vector strength,", "entrained yes/no per animal"],
             fc=ACCFC, ec=ACC, tc=ACC)
b_wave = box(ax, 0.695, IN_Y, 0.245, "Measured waveforms",
             ["6 aligned 24 h profiles,", "one column per animal,", "DD and LD"],
             fc=ACCFC, ec=ACC, tc=ACC)
b_par = box(ax, 0.918, IN_Y, 0.155, "Table 2",
            ["\u03b2, \u03b1, L_base, \u03c4_h,", "\u03ba, y_th per species,", "held fixed"],
            fc="#F4F2EC", ec="0.6")

# --- row 2: the four computational stages
ST_Y, SW = 0.545, 0.215
st = []
for xc, title, lines in [
        (0.135, "A   s1a_dd_library.py", ["20,000 random (k\u2081,k\u2082,k\u2083)", "run 10 d + 8.45 d in darkness"]),
        (0.385, "B   s1b_calibrate.py", ["reweight the library to fit", "k\u0304, \u03c3_scale, \u03c3_asym per species"]),
        (0.635, "C   s1c_population.py", ["draw 500 individuals per", "species"]),
        (0.885, "D   s1d_beta_sweep.py", ["6-segment protocol for every", "individual at 12 values of \u03b2"])]:
    st.append(box(ax, xc, ST_Y, SW, title, lines, mono=False))
for i in range(3):
    arrow(ax, (st[i]["right"], ST_Y - 0.062), (st[i + 1]["left"], ST_Y - 0.062))

# --- feeds from the measured data into the stages that consume them
arrow(ax, (0.135, b_data["bot"]), (0.135, st[0]["top"]), color=ACC, ls="--")
ax.text(0.148, (b_data["bot"] + st[0]["top"]) / 2, "record lengths", fontsize=7.4, color=ACC,
        style="italic", ha="left", va="center")
arrow(ax, (0.245, b_data["bot"]), (0.360, st[1]["top"]), color=ACC, ls="--", rad=-0.12)
ax.text(0.312, (b_data["bot"] + st[1]["top"]) / 2 + 0.012, "calibration targets:\n% rhythmic, \u03c4 mean, \u03c4 SD",
        fontsize=7.4, color=ACC, style="italic", ha="center", va="center")
arrow(ax, (0.885, b_par["bot"]), (0.885, st[3]["top"]), color="0.5", ls="--")

# --- row 3: what each stage writes
for s_, txt in zip(st, ["dd_library.csv\n20,000 rows",
                        "calibration_fit.csv\n3 parameters \u00d7 3 species",
                        "population.csv\n500 \u00d7 3 species",
                        "results_individual.csv  18,000 rows\nprofiles_hourly.csv"]):
    arrow(ax, (s_["xc"], s_["bot"]), (s_["xc"], s_["bot"] - 0.045))
    ax.text(s_["xc"], s_["bot"] - 0.052, txt, ha="center", va="top", fontsize=7.2,
            color=INK, family="monospace", linespacing=1.6)

# --- row 4: analysis, then verification, stacked full width so no arrow crosses a box
AN_Y = 0.215
b_an = box(ax, 0.515, AN_Y, 0.93, "E   s1e_shape_qc.py,  s1f_results_figures.py,  s1g / s1h figure scripts",
           ["quality-control every simulated waveform, compare it with the recorded spiders, compute the \u03b2 response,",
            "and draw every figure and table that appears in the response letter and in Supplement S1"])
b_ver = box(ax, 0.515, b_an["bot"] - 0.035, 0.93, "Verification",
            ["six model tests run before every sweep (\u03b2 inert in darkness, shared gate, known period, periodogram cross-check);",
             "55 published numbers rechecked by verify_reproduction.py; a clean-room rerun from the measured data alone",
             "reproduced all 18,000 result rows and every true/false flag bit for bit"],
            fc="#F4F2EC", ec="0.6")
arrow(ax, (0.885, st[3]["bot"] - 0.105), (0.885, b_an["top"]))
arrow(ax, (0.515, b_an["bot"]), (0.515, b_ver["top"]))

# measured entrainment and waveforms feed the comparison only: route below the input
# row and down the left margin, so no arrow passes through a computational box
XL, YH = 0.022, b_circ["bot"] - 0.030
for bx in (b_circ, b_wave):
    ax.plot([bx["left"] + 0.02, bx["left"] + 0.02], [bx["bot"], YH], color=ACC, ls="--", lw=1.1, zorder=1)
# horizontal run, broken where it would cross the "record lengths" arrow
for x0, x1 in [(b_wave["left"] + 0.02, 0.148), (0.122, XL)]:
    ax.plot([x0, x1], [YH, YH], color=ACC, ls="--", lw=1.1, zorder=1)
ax.plot([XL, XL], [YH, AN_Y - 0.055], color=ACC, ls="--", lw=1.1, zorder=1)
arrow(ax, (XL, AN_Y - 0.055), (b_an["left"], AN_Y - 0.055), color=ACC, ls="--")
ax.text(XL + 0.009, (YH + AN_Y) / 2 - 0.03,
        "observed entrainment percentages, phases and waveforms,\nused for comparison only, never for fitting",
        fontsize=7.4, color=ACC, style="italic", ha="center", va="center", rotation=90)

fig.savefig(os.path.join(FIG, "FigM0_pipeline.png"), dpi=300, bbox_inches="tight")

# %% Figure M1 -- protocol schematic
SEGMENTS = [("0  equilibration", 10.0, "DD", "discarded"),
            ("1  DD assay", None, "DD", "free-running period,\namplitude, waveform"),
            ("2  phase offset", 0.5, "DD", "discarded"),
            ("3  LD transient", 15.0, "LD", "discarded"),
            ("4  LD assay", None, "LD", "entrainment, period, phase,\namplitude, masking, waveform"),
            ("5  DD release", 10.0, "DD", "reported separately")]

fig, ax = plt.subplots(figsize=(11, 4.4))
BADGE = ["0", "1", "2", "3", "4", "5"]
y0 = 0
for row, sp in enumerate(SP):
    x = 0.0
    for j, (name, dur, cond, _) in enumerate(SEGMENTS):
        d = dur if dur is not None else (REC[sp][0] if j == 1 else REC[sp][1])
        if cond == "DD":
            ax.add_patch(Rectangle((x, y0), d, 0.66, facecolor="0.25", edgecolor="w", lw=0.8))
        else:                                   # draw the 12:12 cycle explicitly
            for k in range(int(np.ceil(d)) + 1):
                w_l = min(0.5, max(0.0, d - k))
                if w_l > 0:
                    ax.add_patch(Rectangle((x + k, y0), w_l, 0.66, facecolor="#FFE9A8", edgecolor="none"))
                w_d = min(0.5, max(0.0, d - k - 0.5))
                if w_d > 0:
                    ax.add_patch(Rectangle((x + k + 0.5, y0), w_d, 0.66, facecolor="0.25", edgecolor="none"))
            ax.add_patch(Rectangle((x, y0), d, 0.66, facecolor="none", edgecolor="w", lw=0.8))
        if row == 0:                            # number badges on the top row only
            ax.text(x + d / 2, y0 + 0.95, BADGE[j], ha="center", va="bottom", fontsize=9,
                    weight="bold", color="0.2",
                    bbox=dict(boxstyle="circle,pad=0.22", fc="w", ec="0.5", lw=0.8))
        x += d
    ax.text(-0.8, y0 + 0.33, sp, ha="right", va="center", fontsize=9, color=COL[sp], style="italic")
    ax.text(x + 0.8, y0 + 0.33, f"{x:.1f} d total", ha="left", va="center", fontsize=8, color="0.35")
    y0 -= 1.15

# per-species durations of the two segments that differ, and the key below the bars
ax.text(0, -3.95, "segment", fontsize=8, weight="bold", color="0.2")
ax.text(9.5, -3.95, "duration", fontsize=8, weight="bold", color="0.2")
ax.text(19.5, -3.95, "light", fontsize=8, weight="bold", color="0.2")
ax.text(25, -3.95, "what is recorded", fontsize=8, weight="bold", color="0.2")
rows_txt = [
    ("0  equilibration", "10 d", "DD", "nothing; end state seeds segment 1"),
    ("1  DD assay", "8.45 / 4.76 / 7.46 d", "DD", "free-running period, amplitude, waveform, DD rhythmicity"),
    ("2  phase offset", "U(0, 24) h per individual", "DD", "nothing; randomises the phase at lights-on"),
    ("3  LD transient", "15 d", "LD 12:12", "nothing; discarded so only the steady state is scored"),
    ("4  LD assay", "15.49 / 5.57 / 6.36 d", "LD 12:12", "entrainment, both period estimators, phase, vector strength,"),
    ("", "", "", "amplitude, light/dark ratio, waveform shape"),
    ("5  DD release", "10 d", "DD", "period and amplitude after release; reported separately"),
]
yy = -4.5
for a_, b_, c_, d_ in rows_txt:
    ax.text(0, yy, a_, fontsize=8, color="0.15")
    ax.text(9.5, yy, b_, fontsize=8, color="0.35")
    ax.text(19.5, yy, c_, fontsize=8, color="0.35")
    ax.text(25, yy, d_, fontsize=8, color="0.35")
    yy -= 0.52
ax.text(9.5, yy + 0.1, "durations listed L. cornutus / A. pennsylvanica / S. grossa; segments 1 and 4",
        fontsize=7.5, color="0.5", style="italic")
ax.text(9.5, yy - 0.42, "match that species' own mean recording length", fontsize=7.5, color="0.5", style="italic")

ax.add_patch(Rectangle((46, 1.55), 1.4, 0.34, facecolor="0.25")); ax.text(47.9, 1.72, "darkness", va="center", fontsize=8)
ax.add_patch(Rectangle((53, 1.55), 1.4, 0.34, facecolor="#FFE9A8", ec="0.6", lw=0.5)); ax.text(54.9, 1.72, "light", va="center", fontsize=8)
ax.set_xlim(-9, 64); ax.set_ylim(yy - 1.0, 2.3); ax.axis("off")
ax.set_title("Simulation protocol per individual oscillator (bars to scale, in days)", fontsize=10, loc="left", x=0.0)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "FigM1_protocol.png"), dpi=300)

# %% Figure M2 -- calibration check
fig, ax = plt.subplots(1, 4, figsize=(13, 3.2))
bins = np.arange(14, 40, 1.5)
for c, sp in enumerate(SP):
    m = res1[(res1.variant == "fitted") & (res1.species == sp)].drop_duplicates("id")
    d = pd.read_csv(os.path.join(OUT, "spider_data", COMP[sp]))
    dd = d[(d.Condition == "DD") & (d.Mean_Activity > 0)]
    obs = dd[dd.Period_p_value < 0.05].Period_hours
    sim = m.tau_true_DD[m.rhythmic_DD]
    ax[c].hist(sim, bins=bins, density=True, color=COL[sp], alpha=0.35, label=f"model (n = {len(sim)})")
    ax[c].hist(obs, bins=bins, density=True, histtype="step", lw=2, color="0.15", label=f"spiders (n = {len(obs)})")
    ax[c].axvline(24, color="0.5", lw=0.8, ls=":")
    ax[c].set_title(f"{sp}\nmodel {sim.mean():.2f} ± {sim.std(ddof=1):.2f} h,  "
                    f"spiders {obs.mean():.2f} ± {obs.std(ddof=1):.2f} h", fontsize=8, loc="left")
    ax[c].set_xlabel("free-running period in DD (h)")
    ax[c].legend(fontsize=7, frameon=False)
ax[0].set_ylabel("density")

w = 0.35; xx = np.arange(3)
sim_f = [100 * res1[(res1.variant == "fitted") & (res1.species == s)].drop_duplicates("id").rhythmic_DD.mean() for s in SP]
obs_f = []
for sp in SP:
    d = pd.read_csv(os.path.join(OUT, "spider_data", COMP[sp]))
    dd = d[(d.Condition == "DD") & (d.Mean_Activity > 0)]
    obs_f.append(100 * (dd.Period_p_value < 0.05).mean())
ax[3].bar(xx - w / 2, sim_f, w, color=[COL[s] for s in SP], alpha=0.45, label="model")
ax[3].bar(xx + w / 2, obs_f, w, color=[COL[s] for s in SP], edgecolor="0.15", lw=1.2, label="spiders")
for i, (a, b) in enumerate(zip(sim_f, obs_f)):
    ax[3].text(i - w / 2, a + 1.5, f"{a:.1f}", ha="center", fontsize=7)
    ax[3].text(i + w / 2, b + 1.5, f"{b:.1f}", ha="center", fontsize=7)
ax[3].set_xticks(xx); ax[3].set_xticklabels(["L.c", "A.p", "S.g"]); ax[3].set_ylim(0, 112)
ax[3].set_ylabel("rhythmic in DD (%)"); ax[3].legend(fontsize=7, frameon=False, loc="lower right")
ax[3].set_title("calibration targets", fontsize=8, loc="left")
fig.suptitle("Calibration check: the three fitted parameters per species reproduce the free-running data they were fitted to",
             fontsize=10, x=0.01, ha="left")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "FigM2_calibration_check.png"), dpi=300)

print("wrote FigM1_protocol.png, FigM2_calibration_check.png")
for sp in SP:
    m = res1[(res1.variant == "fitted") & (res1.species == sp)].drop_duplicates("id")
    print(f"  {sp:<18} simulated DD: {100*m.rhythmic_DD.mean():5.1f}% rhythmic, "
          f"tau {m.tau_true_DD[m.rhythmic_DD].mean():.2f} +- {m.tau_true_DD[m.rhythmic_DD].std(ddof=1):.2f} h")
