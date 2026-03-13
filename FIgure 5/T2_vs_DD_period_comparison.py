"""
T2 vs DD Period Comparison Analysis
====================================
Compares circadian periods in constant darkness (DD) vs T2 light cycles.

OUTPUT FILES (saved in same folder as this script):
    {prefix}_period_results.csv  - Per-spider period data
    {prefix}_statistics.csv      - Statistical test results  
    {prefix}_figure.png          - Two-panel figure

USAGE: Set file paths in CONFIGURATION section, then run (F5 in Spyder)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from scipy import stats
from astropy.timeseries import LombScargle
import re
import os


# =============================================================================
# CONFIGURATION
# =============================================================================

T2_ACTIVITY_FILE = r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\figures\Figure 6 masking experiment\Lc_Ag_1204-1223 2025 Monitor1_T2.csv"
DD_RESULTS_FILE = r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\Data analysis\Larenioides\LC_spider_analysis_comprehensive_with_LD_split.csv"
OUTPUT_PREFIX = "T2_vs_DD_comparison"

# Analysis parameters
MIN_PERIOD_HOURS = 14.0
MAX_PERIOD_HOURS = 35.0
MIN_ACTIVITY_BINS = 10
SIGNIFICANCE_THRESHOLD = 0.05

# Figure colors
GREEN_DARK = '#228B22'
GREEN_LIGHT = '#90EE90'


# =============================================================================
# DATA LOADING
# =============================================================================

def load_t2_activity(filepath):
    """Load T2 activity CSV and convert datetime to hours."""
    data = pd.read_csv(filepath)
    datetime_col = data.columns[0]
    data['datetime'] = pd.to_datetime(data[datetime_col])
    data['hours'] = (data['datetime'] - data['datetime'].iloc[0]).dt.total_seconds() / 3600
    
    exclude = ['datetime', 'Light', 'light', datetime_col]
    spider_cols = [col for col in data.columns if col not in exclude]
    return data, spider_cols


def load_dd_data(filepath):
    """Load DD results and filter to DD condition."""
    df = pd.read_csv(filepath)
    return df[df['Condition'] == 'DD'].set_index('Spider_ID')


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def extract_base_id(spider_id):
    """Remove sex indicator (M/F) from spider ID."""
    s = str(spider_id).strip()
    return s[:-1] if s and s[-1].upper() in ['M', 'F'] else s


def find_matching_id(query_id, candidate_ids):
    """Match spider ID by normalized string or number."""
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
    """Run Lomb-Scargle periodogram on activity data."""
    activity = np.asarray(activity, dtype=float)
    
    if np.count_nonzero(activity) < MIN_ACTIVITY_BINS:
        return {'period_hours': np.nan, 'power': np.nan, 'fap': np.nan, 
                'status': 'insufficient_activity'}
    
    frequencies = np.linspace(1/MAX_PERIOD_HOURS, 1/MIN_PERIOD_HOURS, 2000)
    ls = LombScargle(time_hours, activity)
    power = ls.power(frequencies)
    
    peak_idx = np.argmax(power)
    peak_period = 1 / frequencies[peak_idx]
    fap = ls.false_alarm_probability(power[peak_idx])
    
    return {
        'period_hours': peak_period,
        'power': power[peak_idx],
        'fap': fap,
        'status': 'significant' if fap < SIGNIFICANCE_THRESHOLD else 'non_significant'
    }


# =============================================================================
# ANALYSIS
# =============================================================================

def analyze_all_spiders(t2_data, spider_cols, dd_df):
    """Analyze T2 periodicity for all spiders and match to DD data."""
    time_hours = t2_data['hours'].values
    results = []
    
    for t2_id in spider_cols:
        # Find matching DD spider
        base_id = extract_base_id(t2_id)
        dd_match = find_matching_id(base_id, dd_df.index)
        
        # Analyze T2 period
        t2_result = analyze_periodicity(time_hours, t2_data[t2_id].values)
        
        # Build result row
        row = {
            'Spider_ID': t2_id,
            'T2_Period_hours': t2_result['period_hours'],
            'T2_Power': t2_result['power'],
            'T2_FAP': t2_result['fap'],
            'T2_Status': t2_result['status'],
            'DD_Period_hours': dd_df.loc[dd_match, 'Period_hours'] if dd_match else np.nan,
            'DD_Amplitude': dd_df.loc[dd_match, 'Period_Amplitude'] if dd_match else np.nan,
            'DD_p_value': dd_df.loc[dd_match, 'Period_p_value'] if dd_match else np.nan
        }
        results.append(row)
    
    return pd.DataFrame(results)


def run_statistics(results_df):
    """Run paired t-test and Wilcoxon on significant spiders."""
    sig = results_df.dropna(subset=['T2_Period_hours', 'DD_Period_hours'])
    sig = sig[sig['T2_Status'] == 'significant']
    
    if len(sig) < 3:
        return pd.DataFrame()
    
    dd, t2 = sig['DD_Period_hours'].values, sig['T2_Period_hours'].values
    stat1, p1 = stats.ttest_rel(dd, t2)
    stat2, p2 = stats.wilcoxon(dd, t2)
    n_decreased = np.sum(t2 < dd)
    
    return pd.DataFrame([
        {'test': 'Paired t-test', 'statistic': stat1, 'p_value': p1, 'n': len(dd),
         'mean_DD': np.mean(dd), 'std_DD': np.std(dd),
         'mean_T2': np.mean(t2), 'std_T2': np.std(t2)},
        {'test': 'Wilcoxon signed-rank', 'statistic': stat2, 'p_value': p2, 'n': len(dd)},
        {'test': 'Direction count', 'n_decreased': n_decreased, 'n_total': len(dd),
         'proportion_decreased': n_decreased / len(dd)}
    ])


# =============================================================================
# VISUALIZATION
# =============================================================================

def create_figure(results_df, output_path):
    """Create two-panel figure: T2 distribution (left) + paired comparison (right)."""
    t2_valid = results_df.dropna(subset=['T2_Period_hours'])
    paired_sig = results_df.dropna(subset=['T2_Period_hours', 'DD_Period_hours'])
    paired_sig = paired_sig[paired_sig['T2_Status'] == 'significant']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    np.random.seed(42)
    
    # --- Panel A: T2 Distribution ---
    sig_mask = t2_valid['T2_Status'] == 'significant'
    all_t2 = t2_valid['T2_Period_hours'].values
    
    bp = ax1.boxplot([all_t2], positions=[1], widths=0.5, patch_artist=True, showfliers=False)
    bp['boxes'][0].set(facecolor='lightblue', alpha=0.5)
    bp['medians'][0].set(color='navy', linewidth=2)
    
    # Significant (filled) and non-significant (open) points
    for mask, fill in [(sig_mask, GREEN_DARK), (~sig_mask, 'none')]:
        periods = t2_valid.loc[mask, 'T2_Period_hours'].values
        if len(periods) > 0:
            jitter = np.random.uniform(-0.15, 0.15, len(periods))
            label = 'Significant' if fill != 'none' else 'Non-significant'
            ax1.scatter(1 + jitter, periods, facecolors=fill, s=80,
                        edgecolors=GREEN_DARK, linewidths=1.5, label=label, zorder=3)
    
    ax1.axhline(y=24, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax1.set(ylabel='Peak Period (hours)', xlim=(0.4, 1.6), 
            ylim=(MIN_PERIOD_HOURS - 1, MAX_PERIOD_HOURS + 1), xticks=[1], xticklabels=['T2'])
    ax1.legend(loc='upper right', fontsize=9)
    
    # --- Panel B: Paired T2 vs DD ---
    if len(paired_sig) > 0:
        dd = paired_sig['DD_Period_hours'].values
        t2 = paired_sig['T2_Period_hours'].values
        idx = np.argsort(t2)
        dd, t2 = dd[idx], t2[idx]
        
        # Lines and points
        for i in range(len(dd)):
            ax2.plot([0, 1], [t2[i], dd[i]], '-', color='gray', alpha=0.5, lw=1.2)
        
        jit = 0.02
        ax2.scatter(np.random.uniform(-jit, jit, len(t2)), t2, c=GREEN_DARK, s=80, zorder=3)
        ax2.scatter(1 + np.random.uniform(-jit, jit, len(dd)), dd, c=GREEN_LIGHT, 
                    edgecolors=GREEN_DARK, s=80, zorder=3)
        
        # Mean lines
        ax2.hlines(np.mean(t2), -0.15, 0.15, colors='black', linewidths=2, zorder=4)
        ax2.hlines(np.mean(dd), 0.85, 1.15, colors='black', linewidths=2, zorder=4)
        
        # Legend
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=GREEN_DARK, markersize=10, label='T2'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=GREEN_LIGHT,
                   markeredgecolor=GREEN_DARK, markersize=10, label='DD')
        ]
        ax2.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    ax2.axhline(y=24, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax2.set(ylabel='Period (hours)', xlim=(-0.3, 1.3), 
            ylim=(MIN_PERIOD_HOURS - 1, MAX_PERIOD_HOURS + 1), xticks=[0, 1], xticklabels=['T2', 'DD'])
    
    plt.tight_layout()
    plt.savefig(f'{output_path}_figure.png', dpi=300, bbox_inches='tight')
    plt.show()


# =============================================================================
# MAIN
# =============================================================================

def run_analysis():
    """Run complete T2 vs DD comparison analysis."""
    # Get output directory (same folder as this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Load data
    print("Loading data...")
    t2_data, spider_cols = load_t2_activity(T2_ACTIVITY_FILE)
    dd_df = load_dd_data(DD_RESULTS_FILE)
    print(f"  T2: {len(spider_cols)} spiders, DD: {len(dd_df)} spiders")
    
    # Analyze
    print("Analyzing periods...")
    results_df = analyze_all_spiders(t2_data, spider_cols, dd_df)
    
    # Summary
    sig_df = results_df[results_df['T2_Status'] == 'significant'].dropna(subset=['DD_Period_hours'])
    print(f"  Significant with DD match: {len(sig_df)}")
    if len(sig_df) > 0:
        print(f"  DD: {sig_df['DD_Period_hours'].mean():.1f} ± {sig_df['DD_Period_hours'].std():.1f} h")
        print(f"  T2: {sig_df['T2_Period_hours'].mean():.1f} ± {sig_df['T2_Period_hours'].std():.1f} h")
    
    # Save outputs
    print("Saving results...")
    out_path = os.path.join(script_dir, OUTPUT_PREFIX)
    
    results_df.to_csv(f"{out_path}_period_results.csv", index=False)
    
    stats_df = run_statistics(results_df)
    if len(stats_df) > 0:
        stats_df.to_csv(f"{out_path}_statistics.csv", index=False)
    
    create_figure(results_df, out_path)
    print("Done!")
    
    return results_df, stats_df


# Run analysis
results_df, stats_df = run_analysis()
