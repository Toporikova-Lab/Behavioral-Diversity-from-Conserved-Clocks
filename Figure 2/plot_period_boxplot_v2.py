


"""
Circadian Period Analysis - Three Species Comparison
====================================================

Generates period distribution plots for three spider species across lighting
conditions (DD, LD, LL) with percentage labels and statistical analysis.

Significance Criteria:
    - Significant: p < 0.05
    - Non-significant: p >= 0.05

Visual Coding:
    - Significant individuals: filled with condition color, black outline
    - Non-significant individuals: small open circles
    - Box plots: semi-transparent condition colors

Output:
    - PNG figure with three horizontal subplots
    - CSV file with pairwise statistical comparisons

Author: Toporikova Lab
Version: 3.2
Date: January 2025
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from itertools import combinations
import os

# =============================================================================
# CONFIGURATION
# =============================================================================


# Data files for each species (order: left to right in figure)
DATA_FILES = {
    'Larinioides': r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\Data analysis\Larenioides\LC_spider_analysis_comprehensive_with_LD_split.csv",
    'Agelenopsis': r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\Data analysis\Agelenopsis\Ag_spider_analysis_comprehensive_with_LD_split.csv",
    'Steatoda'   : r"C:\Users\toporikovan\Box\students research\Stochastic Spiders\Steatoda Paper\Comparative actibity analysis\Data analysis\Steatoda\Sg_spider_analysis_comprehensive_with_LD_split.csv"
}


SPECIES_ORDER = ['Larinioides', 'Agelenopsis', 'Steatoda']
OUTPUT_FILENAME = 'three_species_period_comparison'
STATS_FILENAME = 'three_species_statistics'

# Plot settings
FIGURE_SIZE = (15, 4)
Y_AXIS_RANGE = (15, 33)
DPI = 300

# Conditions
CONDITIONS = ['DD', 'LD', 'LL']
CONDITION_COLORS = {
    'DD': '#1f77b4',
    'LD': '#9467bd',
    'LL': '#d62728'
}

# Significance criteria
SIGNIFICANCE_ALPHA = 0.05  # p < 0.05 = significant, p >= 0.05 = non-significant

# Visual parameters
MARKER_SIZE_SIGNIFICANT = 80
MARKER_SIZE_NONSIG = 25  # Very small for non-significant
JITTER_AMOUNT = 0.15
RANDOM_SEED = 42
BOXPLOT_ALPHA = 0.4  # Semi-transparent boxplots

# =============================================================================
# DATA LOADING
# =============================================================================

def load_and_prepare_data(filepath):
    """
    Load circadian period data for one species.
    
    Args:
        filepath: Path to CSV file
        
    Returns:
        DataFrame with added 'is_significant' column based on p < 0.05
    """
    df = pd.read_csv(filepath)
    df = df[df['Condition'].isin(CONDITIONS)].copy()
    df = df.dropna(subset=['Period_hours'])
    
    # Determine significance based on p-value threshold
    df['is_significant'] = df['Period_p_value'] < SIGNIFICANCE_ALPHA
    
    return df


def calculate_condition_stats(df):
    """
    Calculate summary statistics for each condition.
    
    Args:
        df: DataFrame with period data
        
    Returns:
        Dictionary with n_total, n_significant, percent_significant for each condition
    """
    stats_dict = {}
    for condition in CONDITIONS:
        subset = df[df['Condition'] == condition]
        n_total = len(subset)
        n_sig = subset['is_significant'].sum()
        pct_sig = (n_sig / n_total * 100) if n_total > 0 else 0
        
        stats_dict[condition] = {
            'n_total': n_total,
            'n_significant': n_sig,
            'percent_significant': pct_sig
        }
    return stats_dict

# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

def perform_pairwise_tests(df):
    """
    Perform pairwise Welch's t-tests with Bonferroni correction.
    
    Compares only spiders with significant periodicity.
    
    Args:
        df: DataFrame with period data
        
    Returns:
        Dictionary mapping (condition1, condition2) to p-value
    """
    results = {}
    condition_data = {}
    
    # Extract significant periods for each condition
    for condition in CONDITIONS:
        sig_periods = df[(df['Condition'] == condition) & 
                        df['is_significant']]['Period_hours'].values
        if len(sig_periods) > 0:
            condition_data[condition] = sig_periods
    
    # Perform pairwise comparisons
    pairs = list(combinations(condition_data.keys(), 2))
    n_comparisons = len(pairs)
    
    for cond1, cond2 in pairs:
        data1 = condition_data[cond1]
        data2 = condition_data[cond2]
        
        if len(data1) >= 2 and len(data2) >= 2:
            _, p_val = stats.ttest_ind(data1, data2, equal_var=False)
            p_val_corrected = min(p_val * n_comparisons, 1.0)
            results[(cond1, cond2)] = p_val_corrected
    
    return results


def format_significance(p_value):
    """Convert p-value to significance marker."""
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_boxplots(ax, df, positions):
    """
    Add semi-transparent box plots for significant periods.
    
    Args:
        ax: Matplotlib axis
        df: DataFrame with period data
        positions: List of x-axis positions
    """
    box_data = []
    box_positions = []
    box_colors = []
    
    for i, condition in enumerate(CONDITIONS):
        sig_periods = df[(df['Condition'] == condition) & 
                        df['is_significant']]['Period_hours'].values
        if len(sig_periods) > 0:
            box_data.append(sig_periods)
            box_positions.append(positions[i])
            box_colors.append(CONDITION_COLORS[condition])
    
    if len(box_data) > 0:
        bp = ax.boxplot(box_data, positions=box_positions, widths=0.6,
                       patch_artist=True, showfliers=False,
                       boxprops=dict(linewidth=1.5, edgecolor='black'),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5),
                       medianprops=dict(linewidth=2, color='black'))
        
        # Apply semi-transparent color fills
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(BOXPLOT_ALPHA)


def plot_scatter_points(ax, df, positions):
    """
    Add individual spider data points with jitter.
    
    Significant spiders: filled with condition color, black outline
    Non-significant spiders: very small open circles
    
    Args:
        ax: Matplotlib axis
        df: DataFrame with period data
        positions: List of x-axis positions
    """
    np.random.seed(RANDOM_SEED)
    
    for i, condition in enumerate(CONDITIONS):
        subset = df[df['Condition'] == condition]
        condition_color = CONDITION_COLORS[condition]
        
        # Plot significant individuals - filled with condition color, black outline
        sig_periods = subset[subset['is_significant']]['Period_hours'].values
        if len(sig_periods) > 0:
            jitter = np.random.uniform(-JITTER_AMOUNT, JITTER_AMOUNT, len(sig_periods))
            x_pos = positions[i] + jitter
            
            ax.scatter(x_pos, sig_periods, c=condition_color, s=MARKER_SIZE_SIGNIFICANT,
                      edgecolors='black', linewidths=1.5, alpha=1.0, zorder=3)
        
        # Plot non-significant individuals - small open circles
        nonsig_periods = subset[~subset['is_significant']]['Period_hours'].values
        if len(nonsig_periods) > 0:
            jitter = np.random.uniform(-JITTER_AMOUNT, JITTER_AMOUNT, len(nonsig_periods))
            x_pos = positions[i] + jitter
            
            ax.scatter(x_pos, nonsig_periods, c='white', s=MARKER_SIZE_NONSIG,
                      edgecolors='gray', linewidths=0.8, alpha=0.6, zorder=3)


def add_percent_labels(ax, condition_stats, positions):
    """
    Add percentage of significant spiders above subplot.
    
    Args:
        ax: Matplotlib axis
        condition_stats: Dictionary with condition statistics
        positions: List of x-axis positions
    """
    label_y = 33.3
    
    for i, condition in enumerate(CONDITIONS):
        pct = condition_stats[condition]['percent_significant']
        ax.text(positions[i], label_y, f"{pct:.0f}%",
               ha='center', va='bottom', fontsize=12, fontweight='bold')


def setup_subplot(ax, is_leftmost):
    """
    Configure subplot appearance.
    
    Args:
        ax: Matplotlib axis
        is_leftmost: Boolean, whether this is the leftmost subplot
    """
    positions = list(range(len(CONDITIONS)))
    
    # Add reference line at 24 hours
    ax.axhline(y=24, color='gray', linestyle='--', 
              linewidth=1.3, alpha=0.9, zorder=0)
    
    # Set x-axis
    ax.set_xticks(positions)
    ax.set_xticklabels(CONDITIONS, fontsize=12)
    
    # Set y-axis (label only on leftmost subplot)
    if is_leftmost:
        ax.set_ylabel('Period (hours)', fontsize=12)
        
        # Add legend for significance markers
        legend_elements = [
            ax.scatter([], [], c='dimgray', edgecolors='black',
                      s=MARKER_SIZE_SIGNIFICANT, linewidths=1.5,
                      label='Significant'),
            ax.scatter([], [], c='white', edgecolors='gray',
                      s=MARKER_SIZE_NONSIG, linewidths=0.8, alpha=0.6,
                      label='Non-significant')
        ]
        ax.legend(handles=legend_elements, loc='lower left', 
                 framealpha=0.9, fontsize=9)
    
    # Extend y-axis to accommodate percentage labels at 33.5
    ax.set_ylim(Y_AXIS_RANGE[0], 34.5)

# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

def save_statistics_to_csv(all_stats, output_dir):
    """
    Save statistical analysis results to CSV file.
    
    Args:
        all_stats: Dictionary containing stats for all species
        output_dir: Directory to save CSV file
    """
    rows = []
    
    for species in SPECIES_ORDER:
        stats_data = all_stats[species]
        condition_stats = stats_data['condition_stats']
        test_results = stats_data['test_results']
        df = stats_data['df']
        
        # Add condition statistics
        for condition in CONDITIONS:
            cond_stat = condition_stats[condition]
            all_periods = df[df['Condition'] == condition]['Period_hours']
            
            row = {
                'Species': species,
                'Comparison': condition,
                'Total_N': cond_stat['n_total'],
                'Significant_N': cond_stat['n_significant'],
                'Percent_Significant': cond_stat['percent_significant'],
                'Mean_Period_All': all_periods.mean() if len(all_periods) > 0 else np.nan,
                'SD_Period_All': all_periods.std() if len(all_periods) > 0 else np.nan
            }
            
            # Add significant-only statistics
            if cond_stat['n_significant'] > 0:
                sig_periods = df[(df['Condition'] == condition) & 
                               df['is_significant']]['Period_hours']
                row['Mean_Period_Significant'] = sig_periods.mean()
                row['SD_Period_Significant'] = sig_periods.std()
            else:
                row['Mean_Period_Significant'] = np.nan
                row['SD_Period_Significant'] = np.nan
            
            rows.append(row)
        
        # Add pairwise comparisons
        for (cond1, cond2), p_val in sorted(test_results.items()):
            comparison_name = f"{cond1} vs {cond2}"
            sig_marker = format_significance(p_val)
            
            row = {
                'Species': species,
                'Comparison': comparison_name,
                'Total_N': '',
                'Significant_N': '',
                'Percent_Significant': '',
                'Mean_Period_All': '',
                'SD_Period_All': '',
                'Mean_Period_Significant': '',
                'SD_Period_Significant': '',
                'P_Value': p_val,
                'Significance': sig_marker
            }
            rows.append(row)
    
    # Create DataFrame and save
    stats_df = pd.DataFrame(rows)
    output_path = os.path.join(output_dir, f'{STATS_FILENAME}.csv')
    stats_df.to_csv(output_path, index=False)
    print(f"Statistics saved: {output_path}")


def print_species_summary(species_name, df, condition_stats, test_results):
    """
    Print summary statistics to console.
    
    Args:
        species_name: Name of species
        df: DataFrame with period data
        condition_stats: Dictionary with condition statistics
        test_results: Dictionary with test results
    """
    print(f"\n{'='*70}")
    print(f"{species_name.upper()}")
    print('='*70)
    
    for condition in CONDITIONS:
        stats_cond = condition_stats[condition]
        all_periods = df[df['Condition'] == condition]['Period_hours']
        
        print(f"\n{condition}: {stats_cond['n_significant']}/{stats_cond['n_total']} "
              f"({stats_cond['percent_significant']:.1f}%) significant")
        
        if len(all_periods) > 0:
            print(f"  All: {all_periods.mean():.2f} ± {all_periods.std():.2f} h")
        
        if stats_cond['n_significant'] > 0:
            sig_periods = df[(df['Condition'] == condition) & 
                           df['is_significant']]['Period_hours']
            print(f"  Significant only: {sig_periods.mean():.2f} ± {sig_periods.std():.2f} h")
    
    print("\nPairwise comparisons (Welch's t-test, Bonferroni corrected):")
    for (cond1, cond2), p_val in sorted(test_results.items()):
        sig_marker = format_significance(p_val)
        print(f"  {cond1} vs {cond2}: p = {p_val:.4f} {sig_marker}")

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def create_three_species_plot():
    """
    Create three-species comparison figure and save statistics.
    
    Returns:
        fig, axes: Matplotlib figure and axes objects
    """
    fig, axes = plt.subplots(1, 3, figsize=FIGURE_SIZE, sharey=True)
    all_stats = {}
    
    # Process each species
    for idx, species in enumerate(SPECIES_ORDER):
        ax = axes[idx]
        data_file = DATA_FILES[species]
        
        # Load and analyze data
        df = load_and_prepare_data(data_file)
        condition_stats = calculate_condition_stats(df)
        test_results = perform_pairwise_tests(df)
        
        # Store for CSV output
        all_stats[species] = {
            'df': df,
            'condition_stats': condition_stats,
            'test_results': test_results
        }
        
        # Create plots
        positions = list(range(len(CONDITIONS)))
        plot_boxplots(ax, df, positions)
        plot_scatter_points(ax, df, positions)
        add_percent_labels(ax, condition_stats, positions)
        setup_subplot(ax, is_leftmost=(idx == 0))
        
        # Print summary
        print_species_summary(species, df, condition_stats, test_results)
    
    # Adjust subplot spacing - increase wspace for more white space between panels
    plt.subplots_adjust(wspace=0.3, left=0.06, right=0.98, top=0.92, bottom=0.1)
    
    # Save outputs
    output_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Save figure
    fig_path = os.path.join(output_dir, f'{OUTPUT_FILENAME}.png')
    fig.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    print(f"\n{'='*70}")
    print(f"Figure saved: {fig_path}")
    
    # Save statistics
    save_statistics_to_csv(all_stats, output_dir)
    print('='*70 + "\n")
    
    return fig, axes

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == '__main__':
    # Verify all data files exist
    missing_files = [f"{sp}: {fp}" for sp, fp in DATA_FILES.items() 
                     if not os.path.exists(fp)]
    
    if missing_files:
        print("ERROR: Missing data files:")
        for file in missing_files:
            print(f"  - {file}")
    else:
        fig, axes = create_three_species_plot()
        plt.close(fig)