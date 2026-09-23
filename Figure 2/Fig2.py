"""
Fig2.py
==========
Generates Fig2_v3.png at 1000 DPI combining six panels:

  A  – Larinioides cornutus raster plots (DD | LD | LL)
  B  – Agelenopsis pennsylvanica raster plots (DD | LD | LL)
  C  – Steatoda grossa raster plots (DD | LD | LL)
  D  – Larinioides cornutus period distribution (DD | LD | LL)
  E  – Agelenopsis pennsylvanica period distribution (DD | LD | LL)
  F  – Steatoda grossa period distribution (DD | LD | LL)

Layout: 2-row × 3-column grid.  Top row = rasters, bottom row = box plots.
Each column corresponds to one species.

Output: Fig2.png saved in the same folder as this script.
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from itertools import combinations


# =============================================================================
# PATHS & GLOBAL SETTINGS
# =============================================================================

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_FOLDER = r"C:\Users\robb\Dropbox\Polo\Dan\Research\Projects\Spider Circadian Oscillator\Masking-article-figures-and-data\Spider-data-for-figures"

LC_PATH = BASE_FOLDER + r"\Larinioides"
AG_PATH = BASE_FOLDER + r"\Agelenopsis"
SG_PATH = BASE_FOLDER + r"\Steatoda"

OUTPUT_FILE = os.path.join(SCRIPT_DIR, "Fig2.png")
DPI_OUTPUT  = 1000

# -------------------------------------------------------------------
# Raster settings
# -------------------------------------------------------------------
N_DAYS = 6   # days per raster column

# Each entry: (species_label, panel_label_raster, panel_label_box,
#              [(csv_path, spider_id, condition_label), ...],
#              period_csv_path)
PANELS = [
    ('Larinioides cornutus', 'A', 'D',
     [(LC_PATH + r"\LC 1006-1104 2025 Monitor2_DD.csv", 'Lc87',  'DD'),
      (LC_PATH + r"\LC 1006-1104 2025 Monitor2_LD.csv", 'Lc14',  'LD'),
      (LC_PATH + r"\LC 12092024 Monitor1_LL.csv",       'LcF8',  'LL')],
     LC_PATH + r"\LC_spider_analysis_comprehensive_with_LD_split.csv"),

    ('Agelenopsis pennsylvanica', 'B', 'E',
     [(AG_PATH + r"\Ag 0918-0926 2025 Monitor2_DD.csv", 'Ag11F', 'DD'),
      (AG_PATH + r"\Ag 0825-0906 2025 Monitor1_LD.csv", 'Ag4F',  'LD'),
      (AG_PATH + r"\Ag 1015-1104 2025 Monitor1_LL.csv", 'Ag45F', 'LL')],
     AG_PATH + r"\Ag_spider_analysis_comprehensive_with_LD_split.csv"),

    ('Steatoda grossa', 'C', 'F',
     [(SG_PATH + r"\StA DD 01182024.csv",       'SgA15', 'DD'),
      (SG_PATH + r"\StB 13-23 09232024_LD.csv", 'SgB13', 'LD'),
      (SG_PATH + r"\StB LL 081620242.csv",      'SgB8',  'LL')],
     SG_PATH + r"\Sg_spider_analysis_comprehensive_with_LD_split.csv"),
]

# -------------------------------------------------------------------
# Box-plot / period settings
# -------------------------------------------------------------------
CONDITIONS       = ['DD', 'LD', 'LL']
CONDITION_COLORS = {'DD': '#1f77b4', 'LD': '#9467bd', 'LL': '#d62728'}
SIGNIFICANCE_ALPHA    = 0.05
MARKER_SIZE_SIGNIFICANT = 80
MARKER_SIZE_NONSIG      = 25
JITTER_AMOUNT  = 0.15
RANDOM_SEED    = 42
BOXPLOT_ALPHA  = 0.4
Y_AXIS_RANGE   = (15, 33)


# =============================================================================
# PANEL A–C  –  RASTER HELPERS
# =============================================================================

def get_days_to_plot(df):
    """Return the last N_DAYS full days (skip first and last partial days)."""
    unique_days = sorted(set(df.index.date))
    return unique_days[-(N_DAYS + 1):-1]


def draw_raster_column(axes, df, spider_id):
    """
    Draw one raster column (N_DAYS rows) into a list of axes.

    Parameters
    ----------
    axes      : list of N_DAYS Axes, top to bottom
    df        : DataFrame with datetime index, 'Light' column, spider column
    spider_id : column name of the spider to plot
    """
    df = df.copy()
    df['date'] = df.index.date
    df['hour'] = df.index.hour + df.index.minute / 60.0
    days_to_plot = get_days_to_plot(df)

    for i, (ax, day) in enumerate(zip(axes, days_to_plot)):
        day_data = df[df['date'] == day]
        hours    = day_data['hour'].values
        activity = (day_data[spider_id] > 0).astype(int).values

        # Dark-period shading
        if 'Light' in day_data.columns:
            dark_mask = (day_data['Light'] == 0).values
            if dark_mask.any():
                ax.fill_between(hours, -0.1, 1.1, where=dark_mask,
                                step='pre', color='lightblue', alpha=0.3, zorder=0)

        # Activity
        ax.fill_between(hours, 0, activity, step='pre',
                        color='black', alpha=0.9, zorder=1)

        ax.set_xlim(0, 24)
        ax.set_ylim(0, 1.0)
        ax.set_ylabel(f'{i}', fontsize=11, rotation=0, ha='right', va='center')
        ax.tick_params(axis='y', left=False, labelleft=False)
        ax.set_yticks([])
        ax.grid(True, alpha=0.3, axis='x')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if i == N_DAYS - 1:
            ax.set_xticks([0, 6, 12, 18, 24])
            ax.set_xticklabels(['0', '6', '12', '18', '24'], fontsize=9)
        else:
            ax.set_xticks([])


# =============================================================================
# PANELS D–F  –  BOX-PLOT HELPERS
# =============================================================================

def load_period_data(filepath):
    """Load and prepare circadian period data for one species."""
    df = pd.read_csv(filepath)
    df = df[df['Condition'].isin(CONDITIONS)].copy()
    df = df.dropna(subset=['Period_hours'])
    df['is_significant'] = df['Period_p_value'] < SIGNIFICANCE_ALPHA
    return df


def calculate_condition_stats(df):
    """Return dict with n_total, n_significant, percent_significant per condition."""
    out = {}
    for cond in CONDITIONS:
        sub    = df[df['Condition'] == cond]
        n_tot  = len(sub)
        n_sig  = int(sub['is_significant'].sum())
        out[cond] = {
            'n_total':             n_tot,
            'n_significant':       n_sig,
            'percent_significant': (n_sig / n_tot * 100) if n_tot > 0 else 0,
        }
    return out


def perform_pairwise_tests(df):
    """Pairwise Welch's t-tests with Bonferroni correction (significant spiders only)."""
    cond_data = {}
    for cond in CONDITIONS:
        periods = df[(df['Condition'] == cond) & df['is_significant']]['Period_hours'].values
        if len(periods) > 0:
            cond_data[cond] = periods

    pairs = list(combinations(cond_data.keys(), 2))
    n_comp = len(pairs)
    results = {}
    for c1, c2 in pairs:
        if len(cond_data[c1]) >= 2 and len(cond_data[c2]) >= 2:
            _, p = stats.ttest_ind(cond_data[c1], cond_data[c2], equal_var=False)
            results[(c1, c2)] = min(p * n_comp, 1.0)
    return results


def plot_boxplots(ax, df, positions):
    """Semi-transparent box plots for significant periods."""
    box_data, box_pos, box_colors = [], [], []
    for i, cond in enumerate(CONDITIONS):
        periods = df[(df['Condition'] == cond) & df['is_significant']]['Period_hours'].values
        if len(periods) > 0:
            box_data.append(periods)
            box_pos.append(positions[i])
            box_colors.append(CONDITION_COLORS[cond])

    if box_data:
        bp = ax.boxplot(box_data, positions=box_pos, widths=0.6,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(linewidth=1.5, edgecolor='black'),
                        whiskerprops=dict(linewidth=1.5),
                        capprops=dict(linewidth=1.5),
                        medianprops=dict(linewidth=2, color='black'))
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(BOXPLOT_ALPHA)


def plot_scatter_points(ax, df, positions):
    """Individual spider points with jitter; filled = significant, open = non-significant."""
    np.random.seed(RANDOM_SEED)
    for i, cond in enumerate(CONDITIONS):
        sub   = df[df['Condition'] == cond]
        color = CONDITION_COLORS[cond]

        sig_periods = sub[sub['is_significant']]['Period_hours'].values
        if len(sig_periods) > 0:
            jit = np.random.uniform(-JITTER_AMOUNT, JITTER_AMOUNT, len(sig_periods))
            ax.scatter(positions[i] + jit, sig_periods, c=color,
                       s=MARKER_SIZE_SIGNIFICANT, edgecolors='black',
                       linewidths=1.5, alpha=1.0, zorder=3)

        nonsig = sub[~sub['is_significant']]['Period_hours'].values
        if len(nonsig) > 0:
            jit = np.random.uniform(-JITTER_AMOUNT, JITTER_AMOUNT, len(nonsig))
            ax.scatter(positions[i] + jit, nonsig, c='white',
                       s=MARKER_SIZE_NONSIG, edgecolors='gray',
                       linewidths=0.8, alpha=0.6, zorder=3)


def add_percent_labels(ax, condition_stats, positions):
    """Percentage of significant spiders printed above each box."""
    for i, cond in enumerate(CONDITIONS):
        pct = condition_stats[cond]['percent_significant']
        ax.text(positions[i], 33.3, f'{pct:.0f}%',
                ha='center', va='bottom', fontsize=12, fontweight='bold')


def setup_box_subplot(ax, is_leftmost):
    """Configure box-plot subplot appearance."""
    positions = list(range(len(CONDITIONS)))
    ax.axhline(y=24, color='gray', linestyle='--', linewidth=1.3, alpha=0.9, zorder=0)
    ax.set_xticks(positions)
    ax.set_xticklabels(CONDITIONS, fontsize=12)
    if is_leftmost:
        ax.set_ylabel('Period (hours)', fontsize=12)
        legend_elements = [
            ax.scatter([], [], c='dimgray', edgecolors='black',
                       s=MARKER_SIZE_SIGNIFICANT, linewidths=1.5, label='Significant'),
            ax.scatter([], [], c='white', edgecolors='gray',
                       s=MARKER_SIZE_NONSIG, linewidths=0.8, alpha=0.6,
                       label='Non-significant'),
        ]
        ax.legend(handles=legend_elements, loc='lower left', framealpha=0.9, fontsize=9)
    ax.set_ylim(Y_AXIS_RANGE[0], 34.5)


# =============================================================================
# MAIN – ASSEMBLE FIGURE
# =============================================================================

def main():
    n_species = len(PANELS)
    n_cond    = 3

    # ------------------------------------------------------------------
    # Figure and outer layout: 2 rows (rasters / box plots) × 3 species
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 12), facecolor='white')

    outer = gridspec.GridSpec(
        2, n_species, figure=fig,
        height_ratios=[1.5, 1],
        wspace=0.15, hspace=0.20,
        left=0.05, right=0.985, top=0.92, bottom=0.06,
    )

    # ------------------------------------------------------------------
    # Build axes
    # ------------------------------------------------------------------
    # axes_raster[s][c] = list of N_DAYS axes for species s, condition c
    axes_raster = []
    # axes_box[s] = single ax for species s period box plot
    axes_box = []

    for s in range(n_species):
        # --- raster inner grid: N_DAYS rows × 3 condition columns ---
        raster_gs = gridspec.GridSpecFromSubplotSpec(
            N_DAYS, n_cond, subplot_spec=outer[0, s],
            hspace=0.02, wspace=0.20,
        )
        species_axes = []
        for c in range(n_cond):
            col_axes = [fig.add_subplot(raster_gs[d, c]) for d in range(N_DAYS)]
            species_axes.append(col_axes)
        axes_raster.append(species_axes)

        # --- box-plot axes ---
        axes_box.append(fig.add_subplot(outer[1, s]))

    # Share y-axis across box plots
    for s in range(1, n_species):
        axes_box[s].sharey(axes_box[0])

    # ------------------------------------------------------------------
    # Fill panels
    # ------------------------------------------------------------------
    for s, (species_name, label_raster, label_box, conditions, period_csv) in enumerate(PANELS):

        # ---- Raster panels (A / B / C) ----
        for c, (csv_path, spider_id, cond_label) in enumerate(conditions):
            print(f"  Loading raster {species_name} {cond_label}: {spider_id}")
            df_raster = pd.read_csv(csv_path, index_col=0, parse_dates=True)
            draw_raster_column(axes_raster[s][c], df_raster, spider_id)

            # Condition label (italic) above the top row of this column
            axes_raster[s][c][0].text(
                0.5, 1.12, f'$\\it{{{cond_label}}}$',
                ha='center', va='bottom', fontsize=15,
                transform=axes_raster[s][c][0].transAxes, clip_on=False,
            )

        # Species name (italic) centered above the middle raster column
        axes_raster[s][1][0].text(
            0.5, 1.38, species_name,
            ha='center', va='bottom', fontsize=17, style='italic',
            transform=axes_raster[s][1][0].transAxes, clip_on=False,
        )

        # ---- Box-plot panels (D / E / F) ----
        print(f"  Loading period data: {species_name}")
        df_period     = load_period_data(period_csv)
        cond_stats    = calculate_condition_stats(df_period)

        positions = list(range(n_cond))
        plot_boxplots(axes_box[s], df_period, positions)
        plot_scatter_points(axes_box[s], df_period, positions)
        add_percent_labels(axes_box[s], cond_stats, positions)
        setup_box_subplot(axes_box[s], is_leftmost=(s == 0))

    # ------------------------------------------------------------------
    # Panel labels – placed in figure coordinates so A/D, B/E, C/F
    # share the same x position and are therefore vertically aligned.
    # ------------------------------------------------------------------
    fig.canvas.draw()
    LABEL_X_OFFSET = 0.013   # leftward gap from the axis left edge (figure fraction)
    LABEL_Y_OFFSET = 0.012   # upward gap from the axis top edge (figure fraction)

    for s, (_, label_raster, label_box, _, _) in enumerate(PANELS):
        pos_raster = axes_raster[s][0][0].get_position()
        pos_box    = axes_box[s].get_position()

        # Both labels share the same x (left edge of the leftmost raster column)
        x_fig = pos_raster.x0 - LABEL_X_OFFSET

        fig.text(x_fig, pos_raster.y1 + LABEL_Y_OFFSET, label_raster,
                 ha='left', va='bottom', fontsize=21, fontweight='bold')
        fig.text(x_fig, pos_box.y1 + LABEL_Y_OFFSET, label_box,
                 ha='left', va='bottom', fontsize=21, fontweight='bold')

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    print(f"\nSaving {OUTPUT_FILE} at {DPI_OUTPUT} DPI...")
    plt.savefig(OUTPUT_FILE, dpi=DPI_OUTPUT, format='png',
                bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()
    print("Done!")

if __name__ == "__main__":
    main()