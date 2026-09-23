"""
Fig6_generate_for_article_v2.py
================================
Generates Fig6.png combining five panels:

  A  – Binary raster plot (17 days, 15-min bins)       [from Fig6_binary_raster_for_article.py]
  B  – T2 masking: Light vs Dark activity comparison   [from Fig6_T2_analysis_for_article.py]
  C  – T2 masking: Masking index distribution          [from Fig6_T2_analysis_for_article.py]
  D  – T2 period distribution (Lomb-Scargle)           [from Fig6_T2_vs_DD_for_article.py]
  E  – T2 vs DD paired period comparison               [from Fig6_T2_vs_DD_for_article.py]

Output: Fig6.png at 1000 DPI, saved in the same folder as this script.
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from pathlib import Path
from scipy import stats
from astropy.timeseries import LombScargle


# =============================================================================
# PATHS & GLOBAL SETTINGS
# =============================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT  = r"C:\Users\robb\Dropbox\Polo\Dan\Research\Projects\Spider Circadian Oscillator\Masking-article-figures-and-data"

RASTER_CSV = os.path.join(DATA_ROOT, r"Figure 6\Lc_Ag_1204-1223 2025 Monitor1.csv")
T2_CSV     = os.path.join(DATA_ROOT, r"Figure 6\Lc_Ag_1204-1223 2025 Monitor1_T2.csv")
DD_CSV     = os.path.join(DATA_ROOT, r"Spider-data-for-figures\Larinioides\LC_spider_analysis_comprehensive_with_LD_split.csv")

OUTPUT_FILE = os.path.join(SCRIPT_DIR, "Fig6.png")
DPI_OUTPUT  = 1000

# Raster settings
N_DAYS        = 17
BIN_MINUTES   = 15
RASTER_SPIDER = "LcF95"

# Period analysis settings
MIN_PERIOD_HOURS      = 14.0
MAX_PERIOD_HOURS      = 35.0
MIN_ACTIVITY_BINS     = 10
SIGNIFICANCE_THRESHOLD = 0.05

# Colors for T2 vs DD panels
GREEN_DARK  = '#228B22'
GREEN_LIGHT = '#90EE90'


# =============================================================================
# PANEL A – BINARY RASTER
# =============================================================================

def aggregate_data(df, bin_minutes):
    """Resample to bin_minutes and binarize (1 if any activity, else 0)."""
    df_agg = df.copy() if bin_minutes == 1 else df.resample(f'{bin_minutes}min').sum()
    spider_cols = [col for col in df_agg.columns if col.upper() != 'LIGHT']
    for col in spider_cols:
        df_agg[col] = (df_agg[col] > 0).astype(int)
    return df_agg


def get_day_boundaries(df):
    """Return list of N_DAYS Timestamps to plot (excludes first/last day, pads if needed)."""
    unique_days = sorted(df.index.normalize().unique())
    if len(unique_days) > 2:
        unique_days = unique_days[1:-1]
    if len(unique_days) >= N_DAYS:
        return unique_days[-N_DAYS:]
    first_day = unique_days[0]
    pad = [first_day - pd.Timedelta(days=N_DAYS - len(unique_days) - i)
           for i in range(N_DAYS - len(unique_days))]
    pad.extend(unique_days)
    return pad


def prepare_day_data(df, day, spider_id):
    """Extract hours, activity, and light arrays for one day."""
    day_start = pd.Timestamp(day)
    day_end   = day_start + pd.Timedelta(days=1)
    day_data  = df[(df.index >= day_start) & (df.index < day_end)]
    if day_data.empty:
        return None, None, None
    hours    = (day_data.index - day_start).total_seconds() / 3600
    activity = day_data[spider_id].values
    light    = day_data['Light'].values if 'Light' in day_data.columns else None
    return hours, activity, light


def add_dark_shading(ax, hours, light):
    """Shade dark periods light blue."""
    dark_start = None
    for h, is_dark in zip(hours, light == 0):
        if is_dark and dark_start is None:
            dark_start = h
        elif not is_dark and dark_start is not None:
            ax.axvspan(dark_start, h, alpha=0.3, color='lightblue', zorder=0)
            dark_start = None
    if dark_start is not None:
        ax.axvspan(dark_start, 24, alpha=0.3, color='lightblue', zorder=0)


def format_raster_axis(ax, day_index):
    """Tick/spine formatting for one raster row."""
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.set_yticks([])
    ax.grid(True, alpha=0.3, axis='x')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if day_index == N_DAYS - 1:
        ax.set_xticks([0, 6, 12, 18, 24])
        ax.set_xticklabels(['00:00', '06:00', '12:00', '18:00', '24:00'], fontsize=8)
        ax.set_xlabel('Time of Day', fontsize=10)
    else:
        ax.set_xticks([])


def plot_raster_panel(raster_axes, csv_path):
    """Fill Panel A: draw one binary raster row per axis in raster_axes."""
    df     = pd.read_csv(csv_path, index_col=0, parse_dates=True)
    df_agg = aggregate_data(df, BIN_MINUTES)
    days   = get_day_boundaries(df_agg)

    for i, (ax, day) in enumerate(zip(raster_axes, days)):
        hours, activity, light = prepare_day_data(df_agg, day, RASTER_SPIDER)
        if hours is None:
            ax.set_xlim(0, 24)
            ax.set_ylim(0, 1.05)
        else:
            bar_width = BIN_MINUTES / 60.0
            ax.bar(hours, activity, width=bar_width, color='black', alpha=0.8, align='edge')
            if light is not None:
                add_dark_shading(ax, hours, light)
            ax.set_xlim(0, 24)
            ax.set_ylim(0, 1.05)
        format_raster_axis(ax, i)
        ax.set_ylabel(f'{i}', fontsize=8, rotation=0, ha='right', va='center')


# =============================================================================
# PANELS B & C – T2 MASKING ANALYSIS
# =============================================================================

def load_t2_masking_data(filepath):
    df      = pd.read_csv(filepath, index_col=0, parse_dates=True)
    spiders = [col for col in df.columns if col.lower() != 'light']
    return df, spiders


def calculate_masking(df, spiders):
    """Compute mean light/dark activity and masking index per spider."""
    light   = df['Light']
    results = []
    for spider in spiders:
        activity  = pd.to_numeric(df[spider], errors='coerce').fillna(0)
        val_light = activity[light == 1].mean()
        val_dark  = activity[light == 0].mean()
        total     = val_dark + val_light
        mi = (val_dark - val_light) / total if total > 0 else np.nan
        results.append({'spider_id': spider,
                        'mean_light': val_light,
                        'mean_dark':  val_dark,
                        'masking_index': mi})
    return pd.DataFrame(results)


def run_masking_statistics(results):
    """Wilcoxon signed-rank (light vs dark) and one-sample t-test (MI vs 0)."""
    light = results['mean_light'].values
    dark  = results['mean_dark'].values
    mi    = results['masking_index'].dropna().values
    stat1, p1 = stats.wilcoxon(light, dark)
    stat2, p2 = stats.ttest_1samp(mi, 0)
    return pd.DataFrame([
        {'test': 'Wilcoxon signed-rank', 'comparison': 'Light vs Dark',
         'statistic': stat1, 'p_value': p1, 'n': len(light)},
        {'test': 'One-sample t-test', 'comparison': 'Masking Index vs 0',
         'statistic': stat2, 'p_value': p2, 'n': len(mi),
         'mean': np.mean(mi), 'std': np.std(mi), 'sem': np.std(mi)/np.sqrt(len(mi))},
    ])


def plot_masking_panels(ax_B, ax_C, results):
    """Fill Panels B and C."""
    light_vals = results['mean_light'].values
    dark_vals  = results['mean_dark'].values
    mi_vals    = results['masking_index'].dropna().values
    n          = len(light_vals)
    np.random.seed(42)

    # Panel B: paired Light vs Dark
    bp = ax_B.boxplot([light_vals, dark_vals], positions=[1, 2], widths=0.5,
                      patch_artist=True, showfliers=False)
    bp['boxes'][0].set(facecolor='#FFD700', alpha=0.5)
    bp['boxes'][1].set(facecolor='#2C3E50', alpha=0.5)

    j1 = np.random.normal(0, 0.06, n)
    j2 = np.random.normal(0, 0.06, n)
    ax_B.scatter(1 + j1, light_vals, c='#FFD700', edgecolor='k', s=40, alpha=0.8, zorder=3)
    ax_B.scatter(2 + j2, dark_vals,  c='#2C3E50', edgecolor='w', s=40, alpha=0.8, zorder=3)
    for i in range(n):
        ax_B.plot([1+j1[i], 2+j2[i]], [light_vals[i], dark_vals[i]],
                  color='gray', alpha=0.3, lw=0.8)

    ymax = max(light_vals.max(), dark_vals.max()) * 1.1
    ax_B.set(xticks=[1, 2], xticklabels=['Light', 'Dark'],
             xlim=[0.4, 2.6], ylim=[0, ymax])
    ax_B.set_ylabel('Mean Crossings/Min', fontsize=10)

    # Panel C: Masking Index
    bp2 = ax_C.boxplot([mi_vals], positions=[1], widths=0.4,
                       patch_artist=True, showfliers=False)
    bp2['boxes'][0].set(facecolor='#3498DB', alpha=0.5)

    jm     = np.random.normal(0, 0.05, len(mi_vals))
    colors = ['#27AE60' if m > 0 else '#E74C3C' for m in mi_vals]
    ax_C.scatter(1 + jm, mi_vals, c=colors, edgecolor='k', s=40, alpha=0.8, zorder=3)
    ax_C.axhline(0, color='red', ls='--', lw=2, alpha=0.7)

    mi_max = max(abs(mi_vals.min()), abs(mi_vals.max())) * 1.15
    ax_C.set(xticks=[1], xticklabels=[''],
             xlim=[0.4, 1.6], ylim=[-mi_max, mi_max])
    ax_C.set_ylabel('Masking Index', fontsize=10)
    ax_C.text(1, -0.75, r'$MI = \frac{Dark - Light}{Dark + Light}$',
              ha='center', fontsize=12)


# =============================================================================
# PANELS D & E – T2 vs DD PERIOD COMPARISON
# =============================================================================

def load_t2_activity(filepath):
    """Load T2 activity CSV and add an 'hours' column."""
    data         = pd.read_csv(filepath)
    datetime_col = data.columns[0]
    data['datetime'] = pd.to_datetime(data[datetime_col])
    data['hours']    = (data['datetime'] - data['datetime'].iloc[0]).dt.total_seconds() / 3600
    exclude     = ['datetime', 'Light', 'light', datetime_col]
    spider_cols = [col for col in data.columns if col not in exclude]
    return data, spider_cols


def load_dd_data(filepath):
    df = pd.read_csv(filepath)
    return df[df['Condition'] == 'DD'].set_index('Spider_ID')


def extract_base_id(spider_id):
    """Strip sex indicator (M/F) from spider ID."""
    s = str(spider_id).strip()
    return s[:-1] if s and s[-1].upper() in ['M', 'F'] else s


def find_matching_id(query_id, candidate_ids):
    """Fuzzy match spider ID across naming conventions."""
    query_norm = re.sub(r'[_\-\s]', '', str(query_id).lower())
    for cand in candidate_ids:
        if re.sub(r'[_\-\s]', '', str(cand).lower()) == query_norm:
            return cand
    query_num = re.search(r'\d+', query_id)
    if query_num:
        for cand in candidate_ids:
            cand_num = re.search(r'\d+', str(cand))
            if cand_num and cand_num.group() == query_num.group():
                return cand
    return None


def analyze_periodicity(time_hours, activity):
    """Lomb-Scargle periodogram; returns dict with period, power, fap, status."""
    activity = np.asarray(activity, dtype=float)
    if np.count_nonzero(activity) < MIN_ACTIVITY_BINS:
        return {'period_hours': np.nan, 'power': np.nan, 'fap': np.nan,
                'status': 'insufficient_activity'}
    frequencies = np.linspace(1/MAX_PERIOD_HOURS, 1/MIN_PERIOD_HOURS, 2000)
    ls          = LombScargle(time_hours, activity)
    power       = ls.power(frequencies)
    peak_idx    = np.argmax(power)
    peak_period = 1 / frequencies[peak_idx]
    fap         = ls.false_alarm_probability(power[peak_idx])
    return {'period_hours': peak_period, 'power': power[peak_idx], 'fap': fap,
            'status': 'significant' if fap < SIGNIFICANCE_THRESHOLD else 'non_significant'}


def analyze_all_spiders(t2_data, spider_cols, dd_df):
    """Run period analysis for all T2 spiders and match to DD results."""
    time_hours = t2_data['hours'].values
    results    = []
    for t2_id in spider_cols:
        base_id   = extract_base_id(t2_id)
        dd_match  = find_matching_id(base_id, dd_df.index)
        t2_result = analyze_periodicity(time_hours, t2_data[t2_id].values)
        row = {
            'Spider_ID':       t2_id,
            'T2_Period_hours': t2_result['period_hours'],
            'T2_Power':        t2_result['power'],
            'T2_FAP':          t2_result['fap'],
            'T2_Status':       t2_result['status'],
            'DD_Period_hours': dd_df.loc[dd_match, 'Period_hours']     if dd_match else np.nan,
            'DD_Amplitude':    dd_df.loc[dd_match, 'Period_Amplitude'] if dd_match else np.nan,
            'DD_p_value':      dd_df.loc[dd_match, 'Period_p_value']   if dd_match else np.nan,
        }
        results.append(row)
    return pd.DataFrame(results)


def run_period_statistics(results_df):
    """Paired t-test and Wilcoxon on spiders with significant T2 periods and DD data."""
    sig = results_df.dropna(subset=['T2_Period_hours', 'DD_Period_hours'])
    sig = sig[sig['T2_Status'] == 'significant']
    if len(sig) < 3:
        return pd.DataFrame()
    dd, t2 = sig['DD_Period_hours'].values, sig['T2_Period_hours'].values
    stat1, p1   = stats.ttest_rel(dd, t2)
    stat2, p2   = stats.wilcoxon(dd, t2)
    n_decreased = np.sum(t2 < dd)
    return pd.DataFrame([
        {'test': 'Paired t-test', 'statistic': stat1, 'p_value': p1, 'n': len(dd),
         'mean_DD': np.mean(dd), 'std_DD': np.std(dd),
         'mean_T2': np.mean(t2), 'std_T2': np.std(t2)},
        {'test': 'Wilcoxon signed-rank', 'statistic': stat2, 'p_value': p2, 'n': len(dd)},
        {'test': 'Direction count', 'n_decreased': n_decreased, 'n_total': len(dd),
         'proportion_decreased': n_decreased / len(dd)},
    ])


def plot_t2_dd_panels(ax_D, ax_E, results_df):
    """Fill Panels D and E."""
    t2_valid   = results_df.dropna(subset=['T2_Period_hours'])
    paired_sig = results_df.dropna(subset=['T2_Period_hours', 'DD_Period_hours'])
    paired_sig = paired_sig[paired_sig['T2_Status'] == 'significant']
    np.random.seed(42)

    # Panel D: T2 period distribution
    sig_mask = t2_valid['T2_Status'] == 'significant'
    all_t2   = t2_valid['T2_Period_hours'].values

    bp = ax_D.boxplot([all_t2], positions=[1], widths=0.5, patch_artist=True, showfliers=False)
    bp['boxes'][0].set(facecolor='lightblue', alpha=0.5)
    bp['medians'][0].set(color='navy', linewidth=2)

    for mask, fill in [(sig_mask, GREEN_DARK), (~sig_mask, 'none')]:
        periods = t2_valid.loc[mask, 'T2_Period_hours'].values
        if len(periods) > 0:
            jitter = np.random.uniform(-0.15, 0.15, len(periods))
            label  = 'Significant' if fill != 'none' else 'Non-significant'
            ax_D.scatter(1 + jitter, periods, facecolors=fill, s=60,
                         edgecolors=GREEN_DARK, linewidths=1.5, label=label, zorder=3)

    ax_D.axhline(y=24, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax_D.set(ylabel='Peak Period (hours)',
             xlim=(0.4, 1.6), ylim=(MIN_PERIOD_HOURS - 1, MAX_PERIOD_HOURS + 1),
             xticks=[1], xticklabels=['T2'])
    ax_D.legend(loc='upper right', fontsize=8)

    # Panel E: paired T2 vs DD
    if len(paired_sig) > 0:
        dd  = paired_sig['DD_Period_hours'].values
        t2p = paired_sig['T2_Period_hours'].values
        idx = np.argsort(t2p)
        dd, t2p = dd[idx], t2p[idx]

        for i in range(len(dd)):
            ax_E.plot([0, 1], [t2p[i], dd[i]], '-', color='gray', alpha=0.5, lw=1.2)

        jit = 0.02
        ax_E.scatter(np.random.uniform(-jit, jit, len(t2p)), t2p,
                     c=GREEN_DARK, s=60, zorder=3)
        ax_E.scatter(1 + np.random.uniform(-jit, jit, len(dd)), dd,
                     c=GREEN_LIGHT, edgecolors=GREEN_DARK, s=60, zorder=3)
        ax_E.hlines(np.mean(t2p), -0.15,  0.15, colors='black', linewidths=2, zorder=4)
        ax_E.hlines(np.mean(dd),   0.85,  1.15, colors='black', linewidths=2, zorder=4)

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=GREEN_DARK,
                   markersize=9, label='T2'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=GREEN_LIGHT,
                   markeredgecolor=GREEN_DARK, markersize=9, label='DD'),
        ]
        ax_E.legend(handles=legend_elements, loc='upper right', fontsize=8)

    ax_E.axhline(y=24, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax_E.set(ylabel='Period (hours)',
             xlim=(-0.3, 1.3), ylim=(MIN_PERIOD_HOURS - 1, MAX_PERIOD_HOURS + 1),
             xticks=[0, 1], xticklabels=['T2', 'DD'])


# =============================================================================
# MAIN – ASSEMBLE FIGURE
# =============================================================================

def main():
    # --- Load and analyze data ---
    print("Loading raster data...")
    # (raster CSV read inside plot_raster_panel)

    print("Loading T2 masking data...")
    t2_mask_df, mask_spiders = load_t2_masking_data(T2_CSV)
    masking_results = calculate_masking(t2_mask_df, mask_spiders)
    masking_stats   = run_masking_statistics(masking_results)

    print(f"  Light: {masking_results['mean_light'].mean():.3f} ± {masking_results['mean_light'].std():.3f}")
    print(f"  Dark:  {masking_results['mean_dark'].mean():.3f} ± {masking_results['mean_dark'].std():.3f}")
    print(f"  MI:    {masking_results['masking_index'].mean():.3f} ± {masking_results['masking_index'].std():.3f}")

    print("Loading T2 activity and DD data...")
    t2_activity_df, t2_spiders = load_t2_activity(T2_CSV)
    dd_df          = load_dd_data(DD_CSV)
    period_results = analyze_all_spiders(t2_activity_df, t2_spiders, dd_df)
    period_stats   = run_period_statistics(period_results)

    sig_df = period_results[period_results['T2_Status'] == 'significant'].dropna(subset=['DD_Period_hours'])
    print(f"  Significant with DD match: {len(sig_df)}")
    if len(sig_df) > 0:
        print(f"  DD: {sig_df['DD_Period_hours'].mean():.1f} ± {sig_df['DD_Period_hours'].std():.1f} h")
        print(f"  T2: {sig_df['T2_Period_hours'].mean():.1f} ± {sig_df['T2_Period_hours'].std():.1f} h")

    # Save CSV outputs alongside the script
    masking_results.to_csv(os.path.join(SCRIPT_DIR, "Fig6_masking_results.csv"), index=False)
    masking_stats.to_csv(os.path.join(SCRIPT_DIR, "Fig6_statistics.csv"), index=False)
    period_results.to_csv(os.path.join(SCRIPT_DIR, "Fig6_T2_vs_DD_period_results.csv"), index=False)
    if len(period_stats) > 0:
        period_stats.to_csv(os.path.join(SCRIPT_DIR, "Fig6_T2_vs_DD_statistics.csv"), index=False)

    # --- Build figure layout ---
    # Left column: Panel A (raster, N_DAYS rows)
    # Right column: 2×2 grid of Panels B, C (top) and D, E (bottom)
    fig = plt.figure(figsize=(14, 10), facecolor='white')

    outer = gridspec.GridSpec(
        1, 2, figure=fig,
        width_ratios=[3.5, 10],
        wspace=0.25,
        left=0.06, right=0.97, top=0.97, bottom=0.06,
    )

    raster_gs   = gridspec.GridSpecFromSubplotSpec(N_DAYS, 1, subplot_spec=outer[0], hspace=0.01)
    raster_axes = [fig.add_subplot(raster_gs[i]) for i in range(N_DAYS)]

    right_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[1],
                                                hspace=0.20, wspace=0.25)
    ax_B = fig.add_subplot(right_gs[0, 0])
    ax_C = fig.add_subplot(right_gs[0, 1])
    ax_D = fig.add_subplot(right_gs[1, 0])
    ax_E = fig.add_subplot(right_gs[1, 1])

    # --- Fill panels ---
    plot_raster_panel(raster_axes, RASTER_CSV)
    plot_masking_panels(ax_B, ax_C, masking_results)
    plot_t2_dd_panels(ax_D, ax_E, period_results)

    # --- Panel labels ---
    raster_axes[0].text(-0.20, 1.08, 'A', fontsize=21, fontweight='bold',
                        ha='left', va='top', transform=raster_axes[0].transAxes)
    for ax, label in [(ax_B, 'B'), (ax_C, 'C'), (ax_D, 'D'), (ax_E, 'E')]:
        ax.text(-0.18, 1.12, label, fontsize=21, fontweight='bold',
                ha='left', va='top', transform=ax.transAxes)

    # --- Save ---
    print(f"\nSaving {OUTPUT_FILE} at {DPI_OUTPUT} DPI...")
    plt.savefig(OUTPUT_FILE, dpi=DPI_OUTPUT, format='png',
                bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()
    print("Done!")


if __name__ == "__main__":
    main()