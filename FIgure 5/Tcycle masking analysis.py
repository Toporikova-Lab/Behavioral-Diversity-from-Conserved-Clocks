"""
T-cycle Masking Analysis
========================
Analyzes light/dark activity differences in spider locomotor activity data
collected under ultradian T-cycles (e.g., T2 = 1h light, 1h dark).

MASKING INDEX FORMULA:
    MI = (Dark - Light) / (Dark + Light)
    
    Interpretation:
      MI > 0  : More active in dark (negative masking / light suppression)
      MI < 0  : More active in light (positive masking / light activation)  
      MI = 0  : No difference between light and dark activity

INPUT FILE FORMAT:
    CSV with datetime index, 'Light' column (0/1), and spider activity columns.
    Example:
        datetime,Light,LcF8,LcF23,AgM5
        2025-01-16 11:00,1,0,5,2
        2025-01-16 11:01,0,3,8,1

OUTPUT FILES (saved in same folder as input):
    {base}_masking_results.csv : Per-animal data
        Columns: spider_id, mean_light, mean_dark, masking_index
        
    {base}_statistics.csv : Statistical tests for manuscript reporting
        Columns: test, comparison, statistic, p_value, n, mean, std, sem
        
    {base}_masking_figure.png : Two-panel figure (300 dpi)

STATISTICAL TESTS:
    1. Wilcoxon signed-rank test: Paired comparison of light vs dark activity
    2. One-sample t-test: Whether masking index differs from zero

USAGE:
    1. Set DATA_FILE path below (line 72)
    2. Run script (F5 in Spyder)

================================================================================
COMMON MODIFICATIONS (for reviewer requests)
================================================================================

CHANGE FIGURE SIZE:
    Line 127: figsize=(10, 5)  # (width, height) in inches

CHANGE COLORS:
    Lines 134-135: Box colors for light/dark
        '#FFD700' = gold (light phase)
        '#2C3E50' = dark blue (dark phase)
    Line 152: Box color for masking index
        '#3498DB' = blue
    Line 155: Point colors for masking index
        '#27AE60' = green (positive MI)
        '#E74C3C' = red (negative MI)

CHANGE POINT SIZE:
    Lines 139, 140, 156: s=50  # scatter point size

CHANGE JITTER (horizontal spread of points):
    Lines 137-138: np.random.normal(0, 0.06, n)  # 0.06 = spread width
    Line 154: np.random.normal(0, 0.05, ...)

CHANGE Y-AXIS LABELS:
    Line 148: ax1.set_ylabel('Mean Crossings/Min')
    Line 161: ax2.set_ylabel('Masking Index')

MOVE/RESIZE FORMULA:
    Line 164: ax2.text(1, -0.75, ...)  # x=1 (centered), y=-0.75
    Line 165: fontsize=11

ADD GRID LINES:
    After line 148: ax1.yaxis.grid(True, ls='--', alpha=0.4)
    After line 161: ax2.yaxis.grid(True, ls='--', alpha=0.4)

CHANGE FIGURE DPI:
    Line 168: dpi=300

ADD SIGNIFICANCE BARS:
    See previous version of this script or ask Claude for code.

================================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

# =============================================================================
# CONFIGURATION - Set your file path here
# =============================================================================
DATA_FILE = r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\figures\Figure 6 masking experiment\T2 analysis\Lc_Ag_1204-1223 2025 Monitor1_T2.csv"


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def load_data(filepath):
    """
    Load activity CSV and auto-detect T-cycle duration.
    
    Parameters
    ----------
    filepath : str
        Path to CSV file with datetime index, 'Light' column, and spider columns
    
    Returns
    -------
    df : DataFrame
        Activity data with datetime index
    spiders : list
        Column names for individual spiders
    cycle_name : str
        Detected cycle name (e.g., 'T2', 'T24')
    """
    df = pd.read_csv(filepath, index_col=0, parse_dates=True)
    spiders = [col for col in df.columns if col.lower() != 'light']
    
    # Detect cycle duration from light-dark transitions
    changes = df['Light'].diff().fillna(0) != 0
    change_times = df.index[changes].tolist()
    if len(change_times) >= 2:
        durations = [(change_times[i] - change_times[i-1]).total_seconds()/3600 
                     for i in range(1, len(change_times))]
        half_period = np.mean(durations)
    else:
        half_period = 12  # Default to T24 if can't detect
    
    cycle_name = f'T{round(half_period * 2)}'
    print(f"Loaded: {len(df)} timepoints, {len(spiders)} spiders, {cycle_name} cycle")
    return df, spiders, cycle_name


def calculate_masking(df, spiders):
    """
    Calculate masking index for each spider.
    
    For each spider, computes mean activity during light and dark phases,
    then calculates MI = (Dark - Light) / (Dark + Light).
    
    Parameters
    ----------
    df : DataFrame
        Activity data with 'Light' column (0=dark, 1=light)
    spiders : list
        Column names for individual spiders
    
    Returns
    -------
    DataFrame with columns: spider_id, mean_light, mean_dark, masking_index
    """
    light = df['Light']
    results = []
    
    for spider in spiders:
        activity = pd.to_numeric(df[spider], errors='coerce').fillna(0)
        val_light = activity[light == 1].mean()
        val_dark = activity[light == 0].mean()
        total = val_dark + val_light
        mi = (val_dark - val_light) / total if total > 0 else np.nan
        results.append({'spider_id': spider, 'mean_light': val_light, 
                        'mean_dark': val_dark, 'masking_index': mi})
    
    return pd.DataFrame(results)


def run_statistics(results):
    """
    Run statistical tests for manuscript reporting.
    
    Tests performed:
        1. Wilcoxon signed-rank: paired light vs dark activity
        2. One-sample t-test: masking index vs zero
    
    Parameters
    ----------
    results : DataFrame
        Output from calculate_masking()
    
    Returns
    -------
    DataFrame with test results including p-values, means, std, sem
    """
    light = results['mean_light'].values
    dark = results['mean_dark'].values
    mi = results['masking_index'].dropna().values
    
    # Wilcoxon signed-rank: Light vs Dark (paired, non-parametric)
    stat1, p1 = stats.wilcoxon(light, dark)
    
    # One-sample t-test: Masking Index vs 0
    stat2, p2 = stats.ttest_1samp(mi, 0)
    
    return pd.DataFrame([
        {'test': 'Wilcoxon signed-rank', 'comparison': 'Light vs Dark',
         'statistic': stat1, 'p_value': p1, 'n': len(light)},
        {'test': 'One-sample t-test', 'comparison': 'Masking Index vs 0',
         'statistic': stat2, 'p_value': p2, 'n': len(mi),
         'mean': np.mean(mi), 'std': np.std(mi), 'sem': np.std(mi)/np.sqrt(len(mi))}
    ])


def plot_masking(results, save_path=None):
    """
    Create two-panel figure showing masking analysis results.
    
    Panel A (left): Paired comparison of activity in light vs dark
    Panel B (right): Distribution of masking index values
    
    Parameters
    ----------
    results : DataFrame
        Output from calculate_masking()
    save_path : str or Path, optional
        If provided, saves figure to this path
    """
    # --- Figure setup ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))  # FIGURE SIZE
    
    light_vals = results['mean_light'].values
    dark_vals = results['mean_dark'].values
    mi_vals = results['masking_index'].dropna().values
    n = len(light_vals)
    
    np.random.seed(42)  # For reproducible jitter
    
    # =========================================================================
    # PANEL A: Light vs Dark activity comparison
    # =========================================================================
    
    # Box plots
    bp = ax1.boxplot([light_vals, dark_vals], positions=[1, 2], widths=0.5,
                     patch_artist=True, showfliers=False)
    bp['boxes'][0].set(facecolor='#FFD700', alpha=0.5)  # LIGHT BOX COLOR (gold)
    bp['boxes'][1].set(facecolor='#2C3E50', alpha=0.5)  # DARK BOX COLOR (dark blue)
    
    # Scatter points with jitter
    j1 = np.random.normal(0, 0.06, n)  # JITTER SPREAD for light
    j2 = np.random.normal(0, 0.06, n)  # JITTER SPREAD for dark
    ax1.scatter(1 + j1, light_vals, c='#FFD700', edgecolor='k', s=50, alpha=0.8, zorder=3)
    ax1.scatter(2 + j2, dark_vals, c='#2C3E50', edgecolor='w', s=50, alpha=0.8, zorder=3)
    
    # Paired lines connecting same individual
    for i in range(n):
        ax1.plot([1+j1[i], 2+j2[i]], [light_vals[i], dark_vals[i]], 
                 color='gray', alpha=0.3, lw=0.8)
    
    # Axis formatting
    ymax = max(light_vals.max(), dark_vals.max()) * 1.1  # Y-AXIS PADDING
    ax1.set(xticks=[1, 2], xticklabels=['Light', 'Dark'], xlim=[0.4, 2.6], ylim=[0, ymax])
    ax1.set_ylabel('Mean Crossings/Min')  # Y-AXIS LABEL
    
    # =========================================================================
    # PANEL B: Masking Index distribution
    # =========================================================================
    
    # Box plot
    bp2 = ax2.boxplot([mi_vals], positions=[1], widths=0.4, patch_artist=True, showfliers=False)
    bp2['boxes'][0].set(facecolor='#3498DB', alpha=0.5)  # MI BOX COLOR (blue)
    
    # Scatter points colored by sign
    jm = np.random.normal(0, 0.05, len(mi_vals))  # JITTER SPREAD
    colors = ['#27AE60' if m > 0 else '#E74C3C' for m in mi_vals]  # GREEN=positive, RED=negative
    ax2.scatter(1 + jm, mi_vals, c=colors, edgecolor='k', s=50, alpha=0.8, zorder=3)
    
    # Zero reference line
    ax2.axhline(0, color='red', ls='--', lw=2, alpha=0.7)
    
    # Axis formatting
    mi_max = max(abs(mi_vals.min()), abs(mi_vals.max())) * 1.15  # Y-AXIS PADDING
    ax2.set(xticks=[1], xticklabels=[''], xlim=[0.4, 1.6], ylim=[-mi_max, mi_max])
    ax2.set_ylabel('Masking Index')  # Y-AXIS LABEL
    
    # Formula annotation
    ax2.text(1, -0.75, r'$MI = \frac{Dark - Light}{Dark + Light}$',  # FORMULA POSITION
             ha='center', fontsize=14)  # FORMULA FONT SIZE
    
    # Save and display
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')  # FIGURE DPI
        print(f"Saved: {save_path}")
    plt.show()


# =============================================================================
# RUN ANALYSIS
# =============================================================================

# Load data
data_path = Path(DATA_FILE)
df, spiders, cycle_name = load_data(DATA_FILE)

# Calculate masking index for each spider
results = calculate_masking(df, spiders)
stats_df = run_statistics(results)

# Print summary to console
print(f"\nLight: {results['mean_light'].mean():.3f} ± {results['mean_light'].std():.3f}")
print(f"Dark:  {results['mean_dark'].mean():.3f} ± {results['mean_dark'].std():.3f}")
print(f"Masking Index: {results['masking_index'].mean():.3f} ± {results['masking_index'].std():.3f}")

# Save results to CSV files
out_dir = data_path.parent
base = data_path.stem
results.to_csv(out_dir / f"{base}_masking_results.csv", index=False)
stats_df.to_csv(out_dir / f"{base}_statistics.csv", index=False)
print(f"\nSaved: {base}_masking_results.csv")
print(f"Saved: {base}_statistics.csv")

# Generate and save figure
plot_masking(results, out_dir / f"{base}_masking_figure.png")
print("\nDone!")