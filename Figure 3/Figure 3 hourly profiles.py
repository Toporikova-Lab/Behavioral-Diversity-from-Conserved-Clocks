#!/usr/bin/env python3
"""
Complete Figure 2 Pipeline: Data Processing + Visualization

This script:
1. Processes raw locomotor activity data for all three spider species
2. Generates aligned hourly activity CSV files for DD, LD, and LL conditions
3. Creates the 9-panel Figure 2 showing all species and conditions

Usage: Configure DATA_FOLDERS and CONFIGS sections, then run.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION - All three species
# =============================================================================

# Data folder paths for each species
DATA_FOLDERS = {
    'Larinioides': r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\Data analysis\Larenioides",
    'Agelenopsis': r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\Data analysis\Agelenopsis",
    'Steatoda': r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\Data analysis\Steatoda"
}

# Configuration for each species
CONFIGS = {
    'Larinioides': {
        'species_name': 'Larinioides cornutus',
        'dd_activity_files': [
            "LC 1006-1104 2025 Monitor2_DD.csv",
            "LC 01162025 Monitor1_DD.csv"
        ],
        'ld_activity_files': [
            "LC 1006-1104 2025 Monitor2_LD.csv",
            "LC 10302024 Monitor2 _LD.csv"
        ],
        'll_activity_files': [
            "LC 12092024 Monitor1_LL.csv",
            "LC 12092024 Monitor2_LL.csv"
        ],
        'period_data_file': "LC_spider_analysis_comprehensive_with_LD_split.csv",
        'output_prefix': 'Larinioides',
        'ld_period_tolerance': 0.05,
    },
    'Agelenopsis': {
        'species_name': 'Agelenopsis pennsylvanica',
        'dd_activity_files': [
            "Ag 0825-0906 2025 Monitor1_DD.csv"
        ],
        'ld_activity_files': [
            "Ag 0825-0906 2025 Monitor1_LD.csv",
            "Ag 1015-1104 2025 Monitor1_LD.csv"
        ],
        'll_activity_files': [
            "Ag 1015-1104 2025 Monitor1_LL.csv"
        ],
        'period_data_file': "Ag_spider_analysis_comprehensive.csv",
        'output_prefix': 'Agelenopsis',
        'ld_period_tolerance': 0.05,
    },
    'Steatoda': {
        'species_name': 'Steatoda grossa',
        'dd_activity_files': [
            "StA DD 01182024.csv",
            "StB DD 01082024.csv"
        ],
        'ld_activity_files': [
            "StA 04032024_LD.csv",
            "StB 1-12  09232024_LD.csv",
            "StB 13-23 09232024_LD.csv"
        ],
        'll_activity_files': [
            "StA LL 12272023.csv",
            "StB LL 081620242.csv"
        ],
        'period_data_file': "Sg_spider_analysis_comprehensive_with_LD_split.csv",
        'output_prefix': 'Steatoda',
        'ld_period_tolerance': 0.05,
    }
}

# Display settings
COLORS = {'DD': '#2E86AB', 'LD': '#A23B72', 'LL': '#C73E1D'}
LD_DARK_START = 12  # ZT hour when dark phase begins
LD_DARK_END = 24    # ZT hour when dark phase ends

# =============================================================================
# DATA PROCESSING FUNCTIONS
# =============================================================================

def get_full_path(species, filename):
    """Join species data folder with filename."""
    return str(Path(DATA_FOLDERS[species]) / filename)


def find_datetime_column(df):
    """Auto-detect datetime column in dataframe."""
    common_names = ['datetime', 'date', 'time', 'timestamp', 'date_time']
    
    for col in df.columns:
        if col.lower() in common_names:
            return col
    
    try:
        pd.to_datetime(df[df.columns[0]].head(10))
        return df.columns[0]
    except:
        raise ValueError(f"Could not find datetime column. Columns: {list(df.columns)}")


def detect_light_schedule(df):
    """
    Detect light-on time from the Light column to properly align ZT0.
    
    Args:
        df: DataFrame with 'Light' column (1=on, 0=off) and 'time_minutes' column
    
    Returns:
        float: minute offset to first lights-on (for ZT0 alignment)
    """
    if 'Light' not in df.columns:
        print("      WARNING: No 'Light' column found, assuming ZT0 at start")
        return 0.0
    
    light = df['Light'].values
    time_min = df['time_minutes'].values
    
    # Find lights-on transitions (0 -> 1)
    light_diff = np.diff(light)
    lights_on_indices = np.where(light_diff == 1)[0] + 1
    
    # Determine first lights-on time for ZT0 alignment
    if len(lights_on_indices) > 0:
        first_lights_on_min = time_min[lights_on_indices[0]]
        print(f"      Detected lights-on at minute {first_lights_on_min:.1f}")
        return first_lights_on_min
    elif light[0] == 1:
        # Data starts during light phase
        print(f"      Data starts in light phase, using start as ZT0")
        return 0.0
    else:
        print("      WARNING: Could not detect lights-on, assuming ZT0 at start")
        return 0.0


def load_activity_files(species, file_list, detect_light=False):
    """
    Load activity data from multiple CSV files.
    
    Args:
        species: Species name
        file_list: List of activity data files
        detect_light: If True, detect and return ZT0 offset from Light column
    
    Returns:
        tuple: (all_data dict, zt0_offset in minutes)
        all_data: {spider_id: {'time_minutes': array, 'activity': array, 'total_minutes': float}}
        zt0_offset: minutes to first lights-on (only for LD condition)
    """
    all_data = {}
    zt0_offset = 0.0
    
    for filename in file_list:
        filepath = get_full_path(species, filename)
        print(f"    Loading: {filename}")
        
        try:
            df = pd.read_csv(filepath)
        except FileNotFoundError:
            print(f"      WARNING: File not found, skipping")
            continue
        
        # Parse datetime
        datetime_col = find_datetime_column(df)
        df['datetime'] = pd.to_datetime(df[datetime_col])
        df = df.sort_values('datetime').reset_index(drop=True)
        df['time_minutes'] = (df['datetime'] - df['datetime'].min()).dt.total_seconds() / 60
        
        # Detect light schedule if requested (for LD condition)
        if detect_light and zt0_offset == 0.0:  # Only detect once
            zt0_offset = detect_light_schedule(df)
        
        # Extract spider columns
        exclude = {'datetime', 'time_minutes', 'Light', datetime_col}
        spider_cols = [c for c in df.columns if c not in exclude]
        
        print(f"      Duration: {df['time_minutes'].max()/60:.1f}h, Spiders: {len(spider_cols)}")
        
        for spider_id in spider_cols:
            all_data[spider_id] = {
                'time_minutes': df['time_minutes'].values,
                'activity': df[spider_id].values,
                'total_minutes': df['time_minutes'].max()
            }
    
    return all_data, zt0_offset


def compute_subjective_hourly_activity(time_minutes, activity, period_hours, zt0_offset=0.0):
    """
    Compute mean activity for each of 24 subjective hours.
    Normalizes by number of timepoints to make results duration-independent.
    
    Args:
        time_minutes: Array of time values in minutes
        activity: Array of activity counts
        period_hours: Circadian period in hours
        zt0_offset: Offset in minutes to align ZT0 (for LD condition)
    
    Returns:
        Array of 24 mean activity values (one per subjective hour)
    """
    period_minutes = period_hours * 60
    
    hourly_activity = np.zeros(24)
    hourly_counts = np.zeros(24)
    
    for t, act in zip(time_minutes, activity):
        # Adjust time by ZT0 offset (important for LD to align to lights-on)
        adjusted_t = t - zt0_offset
        if adjusted_t < 0:
            continue  # Skip data before ZT0
        
        phase = (adjusted_t % period_minutes) / period_minutes
        hour_bin = int(phase * 24) % 24
        hourly_activity[hour_bin] += act
        hourly_counts[hour_bin] += 1
    
    # Compute mean activity per hour (not sum)
    with np.errstate(divide='ignore', invalid='ignore'):
        result = hourly_activity / hourly_counts
        result[hourly_counts == 0] = 0
    
    return result


def align_peak_to_hour_12(hourly_activity):
    """Circularly shift activity so peak is at hour 12."""
    peak_hour = np.argmax(hourly_activity)
    shift = 12 - peak_hour
    return np.roll(hourly_activity, shift)


def analyze_condition(species, file_list, period_df, condition_name,
                      use_fixed_period=False, period_tolerance=0.05):
    """
    Analyze locomotor activity for one condition.
    
    Args:
        species: Species name
        file_list: List of activity data files
        period_df: DataFrame with period information
        condition_name: 'DD', 'LD', or 'LL'
        use_fixed_period: If True, use 24h period for all spiders (LD condition)
        period_tolerance: For LD, only include spiders within this % of 24h
    
    Returns:
        dict: {spider_id: array of 24 hourly activities}
    """
    print(f"\n  Analyzing {condition_name}:")
    
    # Load activity data - detect light schedule for LD condition
    activity_data, zt0_offset = load_activity_files(species, file_list, 
                                                     detect_light=(condition_name == 'LD'))
    if not activity_data:
        print(f"    No activity files loaded for {condition_name}")
        return {}
    
    results = {}
    
    for spider_id in activity_data.keys():
        # Find period info for this spider
        spider_period = period_df[period_df['Spider_ID'] == spider_id]
        
        if len(spider_period) == 0:
            continue
        
        spider_period = spider_period.iloc[0]
        
        # Check periodicity significance: p-value < 0.05
        period_p = spider_period.get('Period_p_value', 1.0)
        if period_p >= 0.05:
            continue
        
        # Also check quality designation
        quality = spider_period.get('Period_Quality', '')
        period_hours = spider_period.get('Period_hours', np.nan)
        
        if pd.isna(period_hours):
            continue
        
        # For LD, also check entrainment
        if use_fixed_period:
            if not is_entrained(period_hours, period_tolerance):
                continue
            period_to_use = 24.0  # Fixed 24h for LD
        else:
            period_to_use = period_hours
        
        # Compute hourly activity with proper ZT0 offset for LD
        data = activity_data[spider_id]
        hourly_act = compute_subjective_hourly_activity(
            data['time_minutes'],
            data['activity'],
            period_to_use,
            zt0_offset=zt0_offset
        )
        
        # For DD and LL, align peak to hour 12
        # For LD, keep ZT alignment (no peak alignment)
        if not use_fixed_period:
            hourly_act = align_peak_to_hour_12(hourly_act)
        
        results[spider_id] = hourly_act
    
    print(f"    Included: {len(results)} spiders with p < 0.05")
    return results


def is_entrained(period_hours, tolerance=0.05):
    """Check if period is within tolerance of 24 hours."""
    min_period = 24 * (1 - tolerance)
    max_period = 24 * (1 + tolerance)
    return min_period <= period_hours <= max_period


def save_results_csv(results, output_path):
    """Save aligned activity results to CSV."""
    if not results:
        print(f"  No data to save for {output_path}")
        return
    
    df = pd.DataFrame({'subjective_hour': range(24)})
    for spider_id, activity in results.items():
        df[spider_id] = activity
    
    df.to_csv(output_path, index=False)
    print(f"  Saved: {output_path}")


def detect_significant_peaks(results, species_name, condition_name, alpha=0.05):
    """
    Identify statistically significant peaks in activity data.
    
    Uses one-sample t-test to determine if activity at each hour is 
    significantly elevated above the overall mean activity level.
    Applies Bonferroni correction for 24 multiple comparisons.
    
    Args:
        results: dict {spider_id: array of 24 hourly activities}
        species_name: Name of species
        condition_name: 'DD', 'LD', or 'LL'
        alpha: Significance level (default 0.05)
    
    Returns:
        list of dicts with peak information
    """
    if not results or len(results) == 0:
        return []
    
    # Convert to array: rows = spiders, columns = hours
    activities = np.array(list(results.values()))
    n_spiders = activities.shape[0]
    
    # Calculate overall mean activity for each spider (across all 24 hours)
    overall_means = np.mean(activities, axis=1)
    
    # Bonferroni correction for 24 comparisons
    corrected_alpha = alpha / 24
    
    peak_results = []
    
    for hour in range(24):
        hour_activities = activities[:, hour]
        
        # Test if this hour's activity is significantly above overall mean
        # One-sample t-test: H0 = hour activity equals overall mean
        # H1 = hour activity > overall mean (one-tailed)
        t_stat, p_value_two_tailed = stats.ttest_rel(hour_activities, overall_means)
        
        # Convert to one-tailed (testing if hour > overall mean)
        if t_stat > 0:  # Activity at this hour is elevated
            p_value = p_value_two_tailed / 2
        else:  # Activity at this hour is below mean
            p_value = 1.0  # Not a peak
        
        # Record results
        mean_activity = np.mean(hour_activities)
        sem_activity = stats.sem(hour_activities)
        
        result = {
            'species': species_name,
            'condition': condition_name,
            'hour': hour,  # Subjective hour for DD/LL, ZT for LD
            'mean_activity': mean_activity,
            'sem': sem_activity,
            'p_value': p_value,
            'significant_bonferroni': p_value < corrected_alpha,
            'n_spiders': n_spiders,
            'bonferroni_threshold': corrected_alpha
        }
        
        peak_results.append(result)
    
    return peak_results


# =============================================================================
# FIGURE GENERATION FUNCTIONS
# =============================================================================

def load_activity_data(filepath):
    """
    Load aligned activity data from CSV file.
    
    Returns:
        tuple: (hours, mean_activity, sem_activity, n_spiders) or None
    """
    if filepath is None or not Path(filepath).exists():
        return None
    
    try:
        df = pd.read_csv(filepath)
        hours = df.iloc[:, 0].values
        activity_data = df.iloc[:, 1:].values
        mean_activity = np.nanmean(activity_data, axis=1)
        sem_activity = stats.sem(activity_data, axis=1, nan_policy='omit')
        n_spiders = activity_data.shape[1]
        return hours, mean_activity, sem_activity, n_spiders
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def plot_activity_panel(ax, data, condition, species_name, show_ylabel=True, show_xlabel=True):
    """Plot activity data on a single panel."""
    color = COLORS[condition]
    
    # Set up panel
    ax.set_xlim([-0.5, 23.5])
    ax.set_ylim([0, 1.6])
    ax.set_xticks(range(0, 24, 2))
    ax.grid(False)
    
    # Add background shading for darkness periods (light blue)
    if condition == 'DD':
        # DD: entire cycle is dark
        ax.axvspan(-0.5, 23.5, alpha=0.3, color='lightblue', zorder=0)
    elif condition == 'LD':
        # LD: only ZT12-24 is dark
        ax.axvspan(LD_DARK_START, LD_DARK_END, alpha=0.3, color='lightblue', zorder=0)
    # LL: no shading (constant light)
    
    # Plot data if available
    if data is not None:
        hours, mean_act, sem_act, n_spiders = data
        
        ax.fill_between(hours, mean_act - sem_act, mean_act + sem_act,
                        alpha=0.3, color=color, zorder=2)
        ax.plot(hours, mean_act, '-o', markersize=5, linewidth=2,
               color=color, zorder=3)
        
        # Reference line
        if condition == 'LD':
            ax.axvline(x=LD_DARK_START, color='black', linestyle='-',
                      linewidth=1.5, alpha=0.7, zorder=4)
        else:
            ax.axvline(x=12, color='gray', linestyle='--',
                      linewidth=1, alpha=0.5, zorder=4)
        
        ax.set_title(condition, fontsize=14, fontweight='bold', color=color)
    else:
        ax.text(12, 0.8, 'No data available', ha='center', va='center',
               fontsize=11, color='gray', style='italic')
        ax.set_title(condition, fontsize=14, fontweight='bold', color=color)
    
    # Y-axis label (only leftmost column)
    if show_ylabel:
        ax.set_ylabel('Mean Activity\n(crossings/min)', fontsize=10)
    else:
        ax.set_yticklabels([])
    
    # X-axis label (only on bottom row - LL condition)
    if show_xlabel:
        if condition == 'LD':
            ax.set_xlabel('ZT', fontsize=10)
        else:
            ax.set_xlabel('Subjective Hour', fontsize=10)
    else:
        # Don't remove tick labels, just don't add xlabel
        pass


def create_figure_2(output_files, output_path='figure2_9panel.png'):
    """
    Create 9-panel Figure 2.
    
    Args:
        output_files: dict mapping (species, condition) -> filepath
        output_path: where to save the figure
    """
    print("\n" + "=" * 70)
    print("CREATING FIGURE 2")
    print("=" * 70)
    
    # Species order: left to right
    species_order = ['Larinioides', 'Agelenopsis', 'Steatoda']
    # Condition order: top to bottom
    condition_order = ['DD', 'LD', 'LL']
    
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    
    # Plot each panel: rows are conditions, columns are species
    for row_idx, condition in enumerate(condition_order):
        for col_idx, species in enumerate(species_order):
            ax = axes[row_idx, col_idx]
            
            # Load data
            filepath = output_files.get((species, condition))
            data = load_activity_data(filepath)
            
            # Determine if we show labels
            show_ylabel = (col_idx == 0)  # Only leftmost column
            show_xlabel = (row_idx == 2)  # Only bottom row (LL)
            
            # Plot
            plot_activity_panel(ax, data, condition, species, 
                              show_ylabel=show_ylabel,
                              show_xlabel=show_xlabel)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.08,
                       wspace=0.15, hspace=0.25)
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved: {output_path}")
    plt.close()


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    print("=" * 70)
    print("FIGURE 2 COMPLETE PIPELINE")
    print("Data Processing + Figure Generation + Peak Detection")
    print("=" * 70)
    
    # Track output files for figure generation
    output_files = {}
    
    # Collect all peak detection results
    all_peak_results = []
    
    # Process each species
    for species in ['Larinioides', 'Agelenopsis', 'Steatoda']:
        print(f"\n{'=' * 70}")
        print(f"PROCESSING: {species}")
        print("=" * 70)
        
        config = CONFIGS[species]
        
        # Load period data
        print("\n  Loading period data...")
        try:
            period_df = pd.read_csv(get_full_path(species, config['period_data_file']))
        except FileNotFoundError:
            print(f"  ERROR: Period data file not found for {species}")
            continue
        
        # Summary
        for cond in ['DD', 'LD', 'LL']:
            total = len(period_df[period_df['Condition'] == cond])
            sig_p = len(period_df[(period_df['Condition'] == cond) &
                                 (period_df['Period_p_value'] < 0.05)])
            highly_sig = len(period_df[(period_df['Condition'] == cond) &
                                      (period_df['Period_Quality'] == 'Highly Significant')])
            print(f"    {cond}: {total} total, {sig_p} with p<0.05, {highly_sig} highly significant")
        
        # Analyze each condition
        dd_periods = period_df[period_df['Condition'] == 'DD']
        ld_periods = period_df[period_df['Condition'] == 'LD']
        ll_periods = period_df[period_df['Condition'] == 'LL']
        
        # DD: individual periods, peak aligned
        dd_results = analyze_condition(
            species, config['dd_activity_files'], dd_periods, "DD",
            use_fixed_period=False
        )
        
        # LD: fixed 24h, ZT aligned
        ld_results = analyze_condition(
            species, config['ld_activity_files'], ld_periods, "LD",
            use_fixed_period=True,
            period_tolerance=config['ld_period_tolerance']
        )
        
        # LL: individual periods, peak aligned
        ll_results = analyze_condition(
            species, config['ll_activity_files'], ll_periods, "LL",
            use_fixed_period=False
        )
        
        # Save CSV files
        print("\n  Saving CSV files:")
        prefix = config['output_prefix']
        
        dd_path = f'{prefix}_dd_aligned_activity.csv'
        ld_path = f'{prefix}_ld_aligned_activity.csv'
        ll_path = f'{prefix}_ll_aligned_activity.csv'
        
        save_results_csv(dd_results, dd_path)
        save_results_csv(ld_results, ld_path)
        save_results_csv(ll_results, ll_path)
        
        # Track files for figure
        output_files[(species, 'DD')] = dd_path if dd_results else None
        output_files[(species, 'LD')] = ld_path if ld_results else None
        output_files[(species, 'LL')] = ll_path if ll_results else None
        
        # Perform peak detection
        print("\n  Detecting significant peaks:")
        for condition, results in [('DD', dd_results), ('LD', ld_results), ('LL', ll_results)]:
            if results:
                peaks = detect_significant_peaks(results, species, condition)
                all_peak_results.extend(peaks)
                
                # Count significant peaks
                n_sig = sum(1 for p in peaks if p['significant_bonferroni'])
                print(f"    {condition}: {n_sig} significant peaks (Bonferroni corrected, α=0.05/24)")
    
    # Save peak detection results
    print("\n" + "=" * 70)
    print("SAVING PEAK DETECTION RESULTS")
    print("=" * 70)
    
    if all_peak_results:
        # Save all peak data
        peak_df = pd.DataFrame(all_peak_results)
        peak_df.to_csv('peak_detection_all_timepoints.csv', index=False)
        print("\nSaved: peak_detection_all_timepoints.csv")
        print(f"  Contains: {len(peak_df)} timepoints across all species/conditions")
        
        # Save only significant peaks
        sig_peaks = peak_df[peak_df['significant_bonferroni']]
        sig_peaks_sorted = sig_peaks.sort_values(['species', 'condition', 'hour'])
        sig_peaks_sorted.to_csv('peak_detection_significant_only.csv', index=False)
        print("\nSaved: peak_detection_significant_only.csv")
        print(f"  Contains: {len(sig_peaks_sorted)} significant peaks")
        
        # Print summary
        print("\nSignificant Peak Summary:")
        for species in ['Larinioides', 'Agelenopsis', 'Steatoda']:
            print(f"\n  {species}:")
            for condition in ['DD', 'LD', 'LL']:
                cond_peaks = sig_peaks[(sig_peaks['species'] == species) & 
                                       (sig_peaks['condition'] == condition)]
                if len(cond_peaks) > 0:
                    peak_hours = cond_peaks['hour'].values
                    peak_pvals = cond_peaks['p_value'].values
                    print(f"    {condition}: Peaks at hours {list(peak_hours)}")
                    for h, p in zip(peak_hours, peak_pvals):
                        print(f"      Hour {h}: p = {p:.2e}")
                else:
                    print(f"    {condition}: No significant peaks")
    
    # Create Figure 2
    create_figure_2(output_files, 'figure2_9panel.png')
    
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()