# -*- coding: utf-8 -*-
"""
Supplementary Figure S2 and Supplementary Table S2
Individual actogram and Lomb-Scargle periodogram atlas.

Behavioral Diversity from Conserved Clocks (Robb et al., J. Theor. Biol.)

For every spider in every lighting condition (DD, LD 12:12, LL) the script
draws a double-plotted actogram beside its Lomb-Scargle periodogram, all on
one common activity scale, and writes:

  Supplementary_Figure_S2_atlas.pdf            the atlas itself
  Supplementary_Table_S2_excluded_recordings.csv   recordings removed by the
                                               empty-day quality filter
  periodogram_summary.csv                      one row per retained recording
  individual_actograms/                        one PNG per recording (optional)

HOW TO RUN
  1. Set DATA_DIR below to the repository's Data folder on your computer.
  2. Run the whole file in Spyder (F5), or cell by cell (Ctrl+Enter).
  Outputs are written to an "output" folder next to this script.

Requires: numpy, pandas, matplotlib, astropy  (see README.md)
"""

# %% ---------------------------------------------------------- user settings
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")      # draws off screen; comment out to see figures in Spyder
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.timeseries import LombScargle

# >>> EDIT THIS: folder that contains the Agelenopsis, Larenioides and Steatoda folders
DATA_DIR = r"C:\path\to\Behavioral-Diversity-from-Conserved-Clocks\Data"

# Outputs go to an "output" folder next to this script
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else os.getcwd()
OUT_DIR = os.path.join(HERE, "output")

# If DATA_DIR was not edited, fall back to ../Data (the layout of the repository)
if not os.path.isdir(DATA_DIR):
    DATA_DIR = os.path.abspath(os.path.join(HERE, "..", "Data"))
if not os.path.isdir(DATA_DIR):
    raise FileNotFoundError("Data folder not found. Set DATA_DIR at the top of the script.")
os.makedirs(OUT_DIR, exist_ok=True)
print("data from :", DATA_DIR)
print("output to :", OUT_DIR)


# %% ------------------------------------------------------ analysis settings
BIN_MIN       = 6            # minutes per actogram bin
PERIOD_RANGE  = (14., 35.)   # h, periodogram window used in the manuscript
FAP_LEVEL     = 0.05         # false-alarm level for the dashed threshold line
ACT_SCALE     = 6.0          # crossings per bin drawn at full row height (all panels)
ROWS_PER_PAGE = 4
MIN_DAYS      = 3            # skip recordings shorter than this
SAVE_PNGS     = True         # one PNG per individual as well as the PDF

# Quality filter: a recording is excluded if the spider made zero beam crossings
# on this many whole days (a whole day = at least 20 h of recording).
# Set to None to keep every recording.
MAX_EMPTY_DAYS = 3

# Species folder name (first 3 letters, case insensitive) -> name used in the paper
SPECIES_BY_PREFIX = {"lar": "Larinioides cornutus",
                     "age": "Agelenopsis pennsylvanica",
                     "ste": "Steatoda grossa"}
SPECIES_ORDER = ["Larinioides cornutus", "Agelenopsis pennsylvanica", "Steatoda grossa"]
COND_ORDER    = ["DD", "LD", "LL"]

# Columns that are not spiders, and animals that are not part of this study
TIME_COLS   = ("datetime", "Date_Time")
DROP_COLS   = {"datetime", "Date_Time", "Light"}
DROP_PREFIX = ("At",)        # At32F..At43F belong to a different species

# CT0 for DD and LL recordings = clock hour of lights-on in the last LD cycle
# the animals experienced. Housing lights came on at 08:00.
CT0_DEFAULT = 8.0
CT0_BY_FILE = {}             # e.g. {"file name.csv": 7.5} if a batch differed


# %% ------------------------------------------------------------ file helpers
def species_of(folder):
    """Paper species name from a data subfolder name (tolerates spelling)."""
    return SPECIES_BY_PREFIX.get(folder[:3].lower(), folder)


def is_monitor_file(path):
    """True for TriKinetics activity exports (first column is a timestamp)."""
    first = pd.read_csv(path, nrows=0).columns[0]
    return first in TIME_COLS


def condition_of(path):
    """DD / LD / LL from the file name."""
    name = os.path.basename(path).upper()
    for c in COND_ORDER:
        if c in name:
            return c
    return "??"


def load_monitor(path):
    """Read one TriKinetics export -> (activity DataFrame, light Series)."""
    df = pd.read_csv(path)
    t = pd.to_datetime(df[df.columns[0]], format="mixed")
    df = df.set_index(pd.DatetimeIndex(t))
    light = df["Light"].astype(float) if "Light" in df.columns else pd.Series(0., df.index)
    keep = [c for c in df.columns if c not in DROP_COLS and not c.startswith(DROP_PREFIX)]
    act = df[keep].apply(pd.to_numeric, errors="coerce").fillna(0.)
    return act, light


def load_lab_stats(data_dir):
    """(spider, condition) -> (period, p) from *_with_LD_split.csv in DATA_DIR."""
    out = {}
    for f in glob.glob(os.path.join(data_dir, "*", "*with_LD_split*.csv")):
        d = pd.read_csv(f)
        for _, r in d.iterrows():
            key = (str(r.Spider_ID).lower().lstrip("s"), str(r.Condition))
            out[key] = (r.get("Period_hours", np.nan), r.get("Period_p_value", np.nan))
    return out


# %% ---------------------------------------------------- time and binning helpers
def bin_series(s, light, bin_min):
    """Sum counts into bin_min bins; return binned counts and mean light."""
    rule = "%dmin" % bin_min
    a = s.resample(rule).sum()
    l = light.resample(rule).mean().reindex(a.index).fillna(0.)
    return a, l


def zt_and_day(index, ref_hour):
    """Time stamps -> (phase 0-24 h from ref_hour, day number)."""
    t0 = index[0].normalize() + pd.Timedelta(hours=ref_hour)
    if t0 > index[0]:
        t0 -= pd.Timedelta(days=1)
    hrs = (index - t0).total_seconds() / 3600.
    return hrs % 24., np.floor(hrs / 24.).astype(int)


def lights_on_hour(light, index):
    """Clock hour of the first dark-to-light transition in an LD recording."""
    v = light.values
    on = np.where((v[1:] > .5) & (v[:-1] <= .5))[0] + 1
    if len(on) == 0:
        return CT0_DEFAULT
    t = index[on[0]]
    return t.hour + t.minute / 60.


def full_day_list(index_dt, day_of_bin):
    """Day numbers that cover at least 20 h of recording."""
    n = pd.Series(day_of_bin, index=index_dt).groupby(day_of_bin).size()
    return list(n[n >= 20 * 60].index)


def empty_days(series, day_of_bin, full_days):
    """Whole days with zero counts: (total, longest consecutive run)."""
    tot = series.groupby(day_of_bin).sum()
    z = sorted(d for d in full_days if tot.get(d, 0) == 0)
    run = best = 0
    prev = None
    for d in z:
        run = run + 1 if prev is not None and d == prev + 1 else 1
        best = max(best, run)
        prev = d
    return len(z), best


# %% ------------------------------------------------------------ Lomb-Scargle
def ls_period(counts_1min):
    """Lomb-Scargle on 1-min counts -> period grid, power, tau, p, significant."""
    y = counts_1min.values.astype(float)
    t = (counts_1min.index - counts_1min.index[0]).total_seconds() / 3600.
    if len(y) < 100 or y.std() == 0:
        return None, None, np.nan, np.nan, False
    ls = LombScargle(t, y)
    periods = np.linspace(PERIOD_RANGE[0], PERIOD_RANGE[1], 800)
    power = ls.power(1. / periods)
    tau = periods[np.argmax(power)]
    p = ls.false_alarm_probability(power.max(), method="baluev",
                                   minimum_frequency=1. / PERIOD_RANGE[1],
                                   maximum_frequency=1. / PERIOD_RANGE[0])
    return periods, power, tau, float(p), bool(p < FAP_LEVEL)


def fap_line(counts_1min):
    """Power level corresponding to FAP_LEVEL, for the dashed threshold."""
    y = counts_1min.values.astype(float)
    t = (counts_1min.index - counts_1min.index[0]).total_seconds() / 3600.
    ls = LombScargle(t, y)
    return ls.false_alarm_level(FAP_LEVEL, method="baluev",
                                minimum_frequency=1. / PERIOD_RANGE[1],
                                maximum_frequency=1. / PERIOD_RANGE[0])


# %% ------------------------------------------------------------ drawing helpers
def shade_dark(ax, cond, light_b, phase):
    """Grey bands where the lights were off."""
    if cond == "LL":
        return
    if cond == "DD":
        ax.axvspan(0, 48, color="0.90", zorder=0)
        return
    dark = phase[light_b.values < 0.5]
    if len(dark) == 0:
        return
    lo, hi = dark.min(), dark.max()
    for off in (0., 24.):
        ax.axvspan(lo + off, min(hi + off, 48.), color="0.85", zorder=0)


def draw_actogram(ax, act_b, light_b, cond, ref_hour, title):
    """Double-plotted actogram, common activity scale, ZT or CT x-axis."""
    phase, day = zt_and_day(act_b.index, ref_hour)
    n_days = int(day.max()) + 1
    h = np.clip(act_b.values / ACT_SCALE, 0, 1) * 0.85
    x = np.concatenate([phase, phase + 24.])
    row = np.concatenate([day, day - 1])
    hh = np.concatenate([h, h])
    keep = (row >= 0) & (row < n_days) & (hh > 0)
    ax.vlines(x[keep], -row[keep], -row[keep] + hh[keep], lw=0.45, color="k")
    shade_dark(ax, cond, light_b, phase)
    ax.set_xlim(0, 48); ax.set_ylim(-n_days + 0.05, 1.0)
    ax.set_xticks(np.arange(0, 49, 6))
    ax.set_yticks(-np.arange(0, n_days) + 0.45)
    ax.set_yticklabels(np.arange(1, n_days + 1), fontsize=5)
    ax.tick_params(axis="x", labelsize=6, length=2)
    ax.set_title(title, fontsize=7, loc="left", pad=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def draw_periodogram(ax, periods, power, tau, p, sig, thr, lab=None):
    """LS periodogram with the FAP threshold and the peak period printed."""
    if periods is None:
        ax.text(.5, .5, "no data", ha="center", va="center", fontsize=6)
        ax.set_axis_off(); return
    ax.plot(periods, power, lw=0.8, color="k")
    ax.axhline(thr, ls="--", lw=0.7, color="crimson")
    ax.axvline(24, ls=":", lw=0.6, color="0.5")
    if sig:
        ax.plot([tau], [power.max()], "o", ms=3, color="crimson")
    txt = r"$\tau$ = %.2f h" % tau + ("\np = %.1e" % p if p < 1e-3 else "\np = %.3f" % p)
    right = tau < 0.5 * (PERIOD_RANGE[0] + PERIOD_RANGE[1])   # keep text off the peak
    xt, ha = (.97, "right") if right else (.03, "left")
    box = dict(fc="white", ec="none", alpha=.75, pad=.6)
    ax.text(xt, .93, txt, transform=ax.transAxes, ha=ha, va="top", fontsize=5.5,
            color="crimson" if sig else "0.4", bbox=box)
    if lab is not None and np.isfinite(lab[0]):
        ax.text(xt, .74, "lab pipeline\n%.2f h, p = %.3g" % lab, transform=ax.transAxes,
                ha=ha, va="top", fontsize=5, color="steelblue", bbox=box)
    ax.set_xlim(*PERIOD_RANGE)
    ax.set_xticks([16, 20, 24, 28, 32])
    ax.tick_params(labelsize=5.5, length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_pair(fig, gs, k, r, act, light, lab):
    """One individual: actogram (left 3/4) and periodogram (right 1/4)."""
    s = act[r.spider]
    a_b, l_b = bin_series(s, light, BIN_MIN)
    periods, power, tau, p, sig = ls_period(s)
    thr = fap_line(s) if periods is not None else np.nan
    ax1 = fig.add_subplot(gs[k, 0:3])
    ax2 = fig.add_subplot(gs[k, 3])
    mean_cpm = a_b.sum() / ((a_b.index[-1] - a_b.index[0]).total_seconds() / 60.)
    title = "%s   |   %s   |   %s   |   mean %.2f crossings/min" % (
        r.spider, r.cond, r.expt, mean_cpm)
    draw_actogram(ax1, a_b, l_b, r.cond, r.ref_hour, title)
    draw_periodogram(ax2, periods, power, tau, p, sig, thr, lab=lab)
    ax1.set_xlabel("ZT (h)" if r.cond == "LD" else "CT (h)", fontsize=6, labelpad=1)
    ax1.set_ylabel("day", fontsize=6)
    ax2.set_xlabel("period (h)", fontsize=6, labelpad=1)
    ax2.set_ylabel("LS power", fontsize=6)


# %% ------------------------------------------------------ find the data files
LAB = load_lab_stats(DATA_DIR)
print("loaded %d lab-pipeline period entries" % len(LAB))

all_csv = sorted(glob.glob(os.path.join(DATA_DIR, "*", "*.csv")))
files = [f for f in all_csv if is_monitor_file(f)]
print("%d activity files found (%d other csv files ignored)" % (len(files), len(all_csv) - len(files)))
for f in files:
    print("   %-26s %s  %s" % (os.path.basename(os.path.dirname(f)), condition_of(f),
                               os.path.basename(f)))


# %% ------------------------------------ index every recording, apply quality filter
records = []
for path in files:
    cond = condition_of(path)
    act, light = load_monitor(path)
    span_days = (act.index[-1] - act.index[0]).total_seconds() / 86400.
    if span_days < MIN_DAYS:
        print("skipped (too short): %s" % os.path.basename(path)); continue
    if cond == "LD":
        ref = lights_on_hour(light, act.index)
    else:
        ref = CT0_BY_FILE.get(os.path.basename(path), CT0_DEFAULT)
    phase, day = zt_and_day(act.index, ref)
    day_s = pd.Series(day, index=act.index)
    full = full_day_list(act.index, day_s)
    for col in act.columns:
        n_empty, n_consec = empty_days(act[col], day_s, full)
        records.append(dict(species=species_of(os.path.basename(os.path.dirname(path))),
                            cond=cond, spider=col, file=os.path.basename(path),
                            path=path, ref_hour=ref,
                            expt=os.path.basename(path).replace(".csv", ""),
                            days=round(span_days, 2), n_full_days=len(full),
                            empty_days=n_empty, max_consec_empty=n_consec))
index_all = pd.DataFrame(records)
print("\n%d individual recordings in %d files" % (len(index_all), len(files)))

if MAX_EMPTY_DAYS is None:
    bad = pd.Series(False, index=index_all.index)
else:
    bad = index_all.empty_days >= MAX_EMPTY_DAYS
index, dropped = index_all[~bad].copy(), index_all[bad].copy()
dropped.drop(columns="path").to_csv(
    os.path.join(OUT_DIR, "Supplementary_Table_S2_excluded_recordings.csv"), index=False)

print("\nexcluded (%s or more empty days): %d" % (MAX_EMPTY_DAYS, len(dropped)))
print(dropped.groupby(["species", "cond"]).size().to_string())
print("\nretained: %d" % len(index))
print(index.groupby(["species", "cond"]).size().to_string())


# %% ------------------------------------------- period statistics for every spider
cache = {p: load_monitor(p) for p in index.path.unique()}
rows = []
for _, r in index.iterrows():
    s = cache[r.path][0][r.spider]
    periods, power, tau, p, sig = ls_period(s)
    lab = LAB.get((r.spider.lower().lstrip("s"), r.cond), (np.nan, np.nan))
    rows.append(dict(species=r.species, cond=r.cond, spider=r.spider, expt=r.expt,
                     days=r.days, n_full_days=r.n_full_days, empty_days=r.empty_days,
                     mean_counts_per_min=round(float(s.mean()), 4),
                     tau_h=round(tau, 3) if np.isfinite(tau) else np.nan,
                     LS_p=p, rhythmic=sig,
                     lab_period_h=lab[0], lab_p=lab[1]))
stats = pd.DataFrame(rows)
stats.to_csv(os.path.join(OUT_DIR, "periodogram_summary.csv"), index=False)
print(stats.groupby(["species", "cond"]).rhythmic.agg(["size", "sum", "mean"]).to_string())


# %% ------------------------------------------------------ PDF: title and summary
def title_page(pdf):
    """Cover page with the figure legend."""
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(.5, .70, "Supplementary Figure S2", ha="center", fontsize=22)
    fig.text(.5, .655, "Individual actograms and Lomb-Scargle periodograms",
             ha="center", fontsize=13)
    cap = ("Every spider recorded in this study is shown once per lighting condition.\n\n"
           "Left: double-plotted actogram of raw locomotor activity, %d-min bins. Bar height is\n"
           "proportional to beam crossings and uses the SAME scale in every panel (full row height\n"
           "= %.0f crossings per %d-min bin), so amplitudes may be compared directly between\n"
           "individuals, species and conditions. Grey shading marks darkness.\n\n"
           "LD panels are plotted on Zeitgeber Time (ZT0 = lights-on). DD and LL panels are plotted\n"
           "on Circadian Time (CT0 = projected lights-on of the last light-dark cycle experienced).\n\n"
           "Right: Lomb-Scargle periodogram over %g-%g h. Red dashed line is the p = %.2f\n"
           "false-alarm level; the peak period and its p value are printed in each panel.\n\n"
           "Recordings with %s or more whole days of zero beam crossings were excluded, in every\n"
           "lighting condition. %d of %d recordings were removed on this criterion; they are listed\n"
           "in Supplementary Table S2.\n"
           % (BIN_MIN, ACT_SCALE, BIN_MIN, PERIOD_RANGE[0], PERIOD_RANGE[1], FAP_LEVEL,
              MAX_EMPTY_DAYS, len(dropped), len(index_all)))
    fig.text(.12, .28, cap, ha="left", va="bottom", fontsize=9.5, linespacing=1.6)
    pdf.savefig(fig); plt.close(fig)


def strip_panel(ax, col, ylab, logy=False):
    """One dot per individual, grouped by species and condition."""
    xs, labels = [], []
    for i, sp in enumerate(SPECIES_ORDER):
        for j, cond in enumerate(COND_ORDER):
            d = stats[(stats.species == sp) & (stats.cond == cond)][col].dropna()
            x = i * 4 + j
            xs.append(x); labels.append(cond)
            if len(d) == 0:
                continue
            jit = (np.random.RandomState(0).rand(len(d)) - .5) * .5
            ax.plot(x + jit, d, "o", ms=3, mfc="none", mec="0.35", mew=.7)
            ax.plot([x - .35, x + .35], [d.median()] * 2, "-", lw=2, color="crimson")
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(ylab, fontsize=9)
    if logy:
        ax.set_yscale("log")
    for i, sp in enumerate(SPECIES_ORDER):
        ax.text(i * 4 + 1, 1.03, sp, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=8.5, style="italic")
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)


def summary_page(pdf):
    """Activity level and period of every individual, as strip plots."""
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 11))
    fig.subplots_adjust(left=.13, right=.95, top=.90, bottom=.10, hspace=.42)
    strip_panel(axes[0], "mean_counts_per_min", "mean activity (crossings/min)", logy=True)
    axes[0].set_title("Activity level of every individual", fontsize=11, loc="left", pad=22)
    strip_panel(axes[1], "tau_h", "Lomb-Scargle peak period (h)")
    axes[1].axhline(24, ls=":", lw=.8, color="0.5")
    axes[1].set_title("Peak period of every individual", fontsize=11, loc="left", pad=22)
    fig.text(.13, .035,
             "Each circle is one individual recording; red bars are medians. Top panel is on a log\n"
             "scale. These distributions show how representative any single actogram is of its species.",
             fontsize=8.5, va="bottom")
    pdf.savefig(fig); plt.close(fig)


# %% ------------------------------------------------------ PDF: individual panels
def section_page(pdf, text, sub):
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(.5, .56, text, ha="center", fontsize=20, style="italic")
    fig.text(.5, .50, sub, ha="center", fontsize=13)
    pdf.savefig(fig); plt.close(fig)


def species_condition_pages(pdf, sp, cond):
    """Section page plus ROWS_PER_PAGE individuals per page."""
    sub = index[(index.species == sp) & (index.cond == cond)]
    if len(sub) == 0:
        return
    section_page(pdf, sp, "%s   -   n = %d individuals" % (cond, len(sub)))
    sub = sub.sort_values(["expt", "spider"]).reset_index(drop=True)
    for start in range(0, len(sub), ROWS_PER_PAGE):
        fig = plt.figure(figsize=(8.5, 11))
        gs = fig.add_gridspec(ROWS_PER_PAGE, 4, hspace=.55, wspace=.45,
                              left=.09, right=.96, top=.94, bottom=.06)
        for k, (_, r) in enumerate(sub.iloc[start:start + ROWS_PER_PAGE].iterrows()):
            act, light = cache[r.path]
            lab = LAB.get((r.spider.lower().lstrip("s"), r.cond))
            panel_pair(fig, gs, k, r, act, light, lab)
        fig.suptitle("%s  -  %s" % (sp, cond), fontsize=10, y=.975)
        pdf.savefig(fig); plt.close(fig)
    print("rendered %s %s (%d panels)" % (sp, cond, len(sub)))


# %% ------------------------------------------------------------ write the PDF
pdf_path = os.path.join(OUT_DIR, "Supplementary_Figure_S2_atlas.pdf")
pdf = PdfPages(pdf_path)
title_page(pdf)
summary_page(pdf)
for sp in SPECIES_ORDER:
    for cond in COND_ORDER:
        species_condition_pages(pdf, sp, cond)
pdf.close()
print("\nwrote", pdf_path)


# %% ------------------------------------------ optional: one PNG per individual
if SAVE_PNGS:
    png_root = os.path.join(OUT_DIR, "individual_actograms")
    for _, r in index.iterrows():
        d = os.path.join(png_root, r.species.replace(" ", "_"), r.cond)
        os.makedirs(d, exist_ok=True)
        act, light = cache[r.path]
        lab = LAB.get((r.spider.lower().lstrip("s"), r.cond))
        fig = plt.figure(figsize=(7.5, 3.2))
        gs = fig.add_gridspec(1, 4, wspace=.45, left=.09, right=.96, top=.86, bottom=.16)
        panel_pair(fig, gs, 0, r, act, light, lab)
        fig.savefig(os.path.join(d, "%s_%s.png" % (r.spider, r.expt)), dpi=200)
        plt.close(fig)
    print("wrote per-individual PNGs to", png_root)
