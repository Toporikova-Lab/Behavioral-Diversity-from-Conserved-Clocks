"""
Fig5-v3.py
==========
Combined Figure 5: Experimental raster plots (Panel A) and model simulations (Panel B).

Panel A: Experimental actogram raster plots for LcF45 and LcF44
Panel B: Model simulation raster plots (Standard, No entrainment, No masking)

Outputs:
  Fig5.png  — 1000 dpi PNG saved in the same directory as this script

Run in Spyder IDE.  Edit USER CONFIGURATION below, then press Run.
"""

import pandas as pd
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import os

# =============================================================================
# USER CONFIGURATION
# =============================================================================
DATAFILES_PATH = r"C:\Users\robb\Dropbox\Polo\Dan\Research\Projects\Spider Circadian Oscillator\Masking-article-figures-and-data\Spider-data-for-figures\Larinioides\pre-processing"
TXT_FILE     = DATAFILES_PATH + r"\LC LD-DD 01162025\LC 01162025 Monitor1.txt"
EXCEL_FILE   = DATAFILES_PATH + r"\LC LD-DD 01162025\Monitors log 01162025.xlsx"
CSV_CHANNELS = DATAFILES_PATH + r"\monitor channels info.csv"
SHEET_NAME   = "Monitor 1"
N_DAYS       = 12              # rows per experimental raster

TARGET_SPIDERS = ['LcF45', 'LcF44']   # left-to-right order in panel A

# =============================================================================
# MODEL PROTOCOL SETTINGS
# =============================================================================
LD_DAYS_SIM  = 30   # days simulated in LD (transients decay here)
LD_DAYS_SHOW = 4    # last N LD days shown
DD_DAYS      = 8    # days in constant darkness (all shown)
DT           = 0.05 # integration time step (hours)
EQUIL_DAYS   = 10   # pre-protocol LD equilibration (not plotted)

# =============================================================================
# FONT SIZES  — adjust all text here
# =============================================================================
FS_TITLE      = 9    # panel titles (spider ID / model name + subtitle)
FS_XLABEL     = 8    # x-axis label "Zeitgeber Time (h)"
FS_TICK       = 7    # x-axis and y-axis tick labels
FS_PANEL      = 13   # bold A / B panel labels
FS_LD_DD      = 8    # LD / DD margin labels
FS_LEGEND     = 8    # legend text

# =============================================================================
# PHASE-ARROW COORDINATES  (data units: x = ZT hours, y = day number)
# Each entry: (tail_zt, tail_day, head_zt, head_day)
# tail = plain end (late DD, bottom of plot); head = arrowhead (early DD, top)
# Measured from Fig5-article.png; adjust here if needed.
# =============================================================================
# Experimental panels (Panel A)
ARROW_LCF45 = (21.5, 11, 10.5,  4)   # phase delay  (FRP > 24 h)
ARROW_LCF44 = ( 5, 11, 9.0,  4)   # phase advance (FRP < 24 h)

# Model panels (Panel B); LD/DD boundary is at day = LD_DAYS_SHOW = 4
ARROW_STANDARD      = (11.8, 11.2, 15.0,  3.4)
ARROW_NO_ENTRAIN    = ( 1.0, 11.2,  5.0,  3.4)
ARROW_NO_MASK       = (9.0, 11.6, 14.0,  3.2)

# Double-headed black arrow in the No-entrainment panel
# (horizontal, shows phase gap between free-running and entrained clocks)
DBL_ARROW_NO_ENTRAIN = (4.8, 13, 3.5)   # (zt_left, zt_right, day_level)

# =============================================================================
# MODEL PARAMETERS — Larinioides cornutus
# =============================================================================
LARINIOIDES_BASE = dict(
    v1=0.84, v2=0.84, v3=0.84,
    k1=0.18, k2=0.18, k3=0.18,
    n=12,
    K=1.0,
    sigmoid_steepness=12,
    L_amplitude=2.5,
    L_baseline=1.0,
    tau_m=2.0,
    tau_h=2.0,
    h_steepness=30,
    y_threshold=0.5,
    x0=0.1, y0=0.2, z0=1.5, h0=1.0, M0=0.0,
)

VARIANTS = [
    dict(
        title=r'$\it{L. cornutus}$' + '  Standard',
        subtitle='beta = 0.09,  alpha = 0.9',
        light_sensitivity=0.09,
        masking_strength=0.9,
    ),
    dict(
        title='No entrainment',
        subtitle='beta = 0,  alpha = 0.9',
        light_sensitivity=0.0,
        masking_strength=0.9,
    ),
    dict(
        title='No masking',
        subtitle='beta = 0.09,  alpha = 0',
        light_sensitivity=0.09,
        masking_strength=0.0,
    ),
]

# =============================================================================
# EXPERIMENTAL DATA LOADING
# =============================================================================

def load_raw_data(txt_file, channel_csv, excel_file, sheet_name):
    """Parse Trikinetics .txt file, rename sensor columns to animal IDs."""
    info    = pd.read_csv(channel_csv, index_col=0)
    animals = pd.read_excel(excel_file, sheet_name=sheet_name)

    df = pd.read_csv(txt_file, delimiter='\t', header=0,
                     names=info['Channel Info'])
    df['datetime'] = pd.to_datetime(df['Date'] + ' ' + df['Time'])
    df = df.set_index('datetime').drop(['Date', 'Time'], axis=1)

    keep = [c for c in df.columns if c.startswith('s') or c == 'Light']
    df   = df[keep]

    for _, row in animals.iterrows():
        ch = f"s{row['Channel']}"
        if ch in df.columns:
            animal_id = f"{row['Specie abbreviation']}{str(row['Subject ID']).replace(' ', '')}"
            df = df.rename(columns={ch: animal_id})

    return df


def find_first_light_onset(df):
    """Return first datetime where Light transitions 0 to 1."""
    if 'Light' not in df.columns:
        return df.index[0]
    light = df['Light']
    if light.iloc[0] == 1:
        return df.index[0]
    transitions = light.diff().fillna(0)
    on_times    = df.index[transitions > 0]
    return on_times[0] if len(on_times) > 0 else df.index[0]


# =============================================================================
# MODEL EQUATIONS
# =============================================================================

def light_LD(t):
    """LD 12:12: light on ZT0-12, off ZT12-24."""
    return 1.0 if (t % 24) < 12 else 0.0


def odes(state, t, p, light_fn):
    x, y, z, h, M = state
    light  = light_fn(t)
    k2_eff = p['k2'] * (1.0 + p['light_sensitivity'] * light)

    dx = p['v1'] * p['K']**p['n'] / (p['K']**p['n'] + z**p['n']) - p['k1'] * x
    dy = p['v2'] * x - k2_eff * y
    dz = p['v3'] * y - p['k3'] * z

    h_inf = 1.0 / (1.0 + np.exp(p['h_steepness'] * (y - p['y_threshold'])))
    dh    = (h_inf - h) / p['tau_h']

    M_target = 1.0 if light > 0 else 0.0
    dM = (M_target - M) / p['tau_m']

    return [dx, dy, dz, dh, dM]


def compute_activity(state, p):
    x, y, z, h, M = state.T
    m = 1.0 / (1.0 + np.exp(-p['sigmoid_steepness'] * (y - p['y_threshold'])))
    L = (p['L_baseline'] + p['L_amplitude'] * m * h) * (1.0 - p['masking_strength'] * M)
    return y, L


def equilibrate(p):
    t_eq = np.arange(0, EQUIL_DAYS * 24, DT)
    ic   = [p['x0'], p['y0'], p['z0'], p['h0'], p['M0']]
    sol  = odeint(odes, ic, t_eq, args=(p, light_LD))
    return sol[-1, :]


def run_experiment(variant):
    p = {**LARINIOIDES_BASE,
         'light_sensitivity': variant['light_sensitivity'],
         'masking_strength':  variant['masking_strength']}

    ic = equilibrate(p)

    t_ld   = np.arange(0, LD_DAYS_SIM * 24, DT)
    sol_ld = odeint(odes, ic, t_ld, args=(p, light_LD))

    t_dd   = np.arange(0, DD_DAYS * 24, DT)
    sol_dd = odeint(odes, sol_ld[-1, :], t_dd, args=(p, lambda t: 0.0))

    spd      = int(24 / DT)
    trim     = (LD_DAYS_SIM - LD_DAYS_SHOW) * spd
    sol_show = sol_ld[trim:]
    t_show   = t_ld[trim:]

    state_all = np.vstack([sol_show, sol_dd])
    t_all     = np.concatenate([t_show, t_show[-1] + DT + t_dd])

    y_all, L_all = compute_activity(state_all, p)
    return t_all, y_all, L_all


# =============================================================================
# DRAWING — experimental raster (single-axes style)
# =============================================================================

def draw_experimental_raster(ax, df, spider_id, zt0, n_days, show_yaxis):
    """
    Draw one experimental actogram on a single Axes.
    Each day occupies a y-band of height 1; day 0 at top (y-axis inverted).
    """
    for i in range(n_days):
        day_start = zt0 + pd.Timedelta(days=i)
        day_end   = day_start + pd.Timedelta(hours=24)
        mask      = (df.index >= day_start) & (df.index < day_end)
        day_data  = df[mask]

        if day_data.empty or spider_id not in day_data.columns:
            continue

        hours    = (day_data.index - day_start).total_seconds() / 3600.0
        activity = (day_data[spider_id] > 0).astype(float)

        # Dark-period shading
        if 'Light' in day_data.columns:
            dark = (day_data['Light'] == 0).values
            ax.fill_between(hours, i, i + 1, where=dark,
                            step='pre', color='lightblue', alpha=0.45, zorder=0)

        # Activity bars
        ax.fill_between(hours, i, i + activity, step='pre',
                        color='black', alpha=0.9, zorder=1)

    ax.set_xlim(0, 24)
    ax.set_ylim(0, n_days)
    ax.invert_yaxis()
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_xticklabels(['ZT0', 'ZT6', 'ZT12', 'ZT18', 'ZT24'])
    ax.tick_params(axis='x', labelsize=FS_TICK)
    ax.set_xlabel('Zeitgeber Time (h)', fontsize=FS_XLABEL)

    if show_yaxis:
        ax.set_yticks(np.arange(n_days) + 0.5)
        ax.set_yticklabels([str(d) for d in range(n_days)], fontsize=FS_TICK)
    else:
        ax.set_yticks([])


# =============================================================================
# DRAWING — model raster
# =============================================================================

def draw_model_raster(ax, t_all, y_all, L_all, total_days, show_yaxis):
    """Draw one model raster panel."""
    spd    = int(24 / DT)
    y_norm = np.max(y_all) if np.max(y_all) > 0 else 1.0
    L_norm = np.max(L_all) if np.max(L_all) > 0 else 1.0

    for day in range(total_days):
        i0 = day * spd
        i1 = min(i0 + spd, len(t_all))
        if i0 >= len(t_all):
            break

        t_seg = t_all[i0:i1] % 24
        y_seg = y_all[i0:i1]
        L_seg = L_all[i0:i1]

        # Blue shading: scotophase in LD; full row in DD
        if day < LD_DAYS_SHOW:
            ax.axhspan(day, day + 1, xmin=12/24, xmax=1.0,
                       color='lightblue', alpha=0.45, zorder=0)
        else:
            ax.axhspan(day, day + 1,
                       color='lightblue', alpha=0.45, zorder=0)

        ax.plot(t_seg, day + 1 - y_seg / y_norm * 0.85,
                color='#999999', linewidth=0.9, zorder=1,
                solid_capstyle='round')
        ax.plot(t_seg, day + 1 - L_seg / L_norm * 0.85,
                color='#C0392B', linewidth=2.0, zorder=2,
                solid_capstyle='round')

    # LD/DD boundary
    ax.axhline(LD_DAYS_SHOW, color='black', linewidth=1.2,
               linestyle='--', alpha=0.7, zorder=3)

    ax.set_xlim(0, 24)
    ax.set_ylim(0, total_days)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.tick_params(axis='x', labelsize=FS_TICK)
    ax.invert_yaxis()
    ax.set_xlabel('Zeitgeber Time (h)', fontsize=FS_XLABEL)

    if show_yaxis:
        ax.set_yticks(np.arange(total_days) + 0.5)
        ax.set_yticklabels([str(d) for d in range(total_days)], fontsize=FS_TICK)
    else:
        ax.set_yticks([])


def add_panel_label(ax, label, x=-0.14, y=1.12):
    """Add bold panel label (A, B) above the top-left of axes."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight='bold', va='bottom', ha='left',
            clip_on=False)


def draw_phase_arrow(ax, tail_zt, tail_day, head_zt, head_day):
    """
    Draw a blue dashed arrow from (tail_zt, tail_day) to (head_zt, head_day).
    Uses data coordinates (ZT on x, day on y).
    The arrowhead lands at the head point.
    """
    ax.annotate('',
                xy=(head_zt, head_day),
                xytext=(tail_zt, tail_day),
                xycoords='data', textcoords='data',
                arrowprops=dict(
                    arrowstyle='->',
                    color='blue',
                    lw=1.5,
                    linestyle='dashed',
                    mutation_scale=10,
                ),
                zorder=5)


def draw_double_arrow(ax, zt_left, zt_right, day_level):
    """
    Draw a horizontal double-headed black arrow at the given day level,
    spanning from zt_left to zt_right.
    """
    ax.annotate('',
                xy=(zt_right, day_level),
                xytext=(zt_left, day_level),
                xycoords='data', textcoords='data',
                arrowprops=dict(
                    arrowstyle='<->',
                    color='black',
                    lw=2.5,
                    mutation_scale=10,
                ),
                zorder=5)


def add_ld_dd_labels(ax, total_days):
    """Add LD / DD condition labels to the right margin of an axes."""
    ld_frac = (LD_DAYS_SHOW / 2) / total_days
    dd_frac = (LD_DAYS_SHOW + DD_DAYS / 2) / total_days
    ax.text(1.02, 1 - ld_frac, 'LD', transform=ax.transAxes,
            va='center', ha='left', fontsize=FS_LD_DD, fontweight='bold')
    ax.text(1.02, 1 - dd_frac, 'DD', transform=ax.transAxes,
            va='center', ha='left', fontsize=FS_LD_DD, fontweight='bold')


# =============================================================================
# FIGURE ASSEMBLY
# =============================================================================

def build_figure():
    total_model_days = LD_DAYS_SHOW + DD_DAYS   # 12, matches N_DAYS

    # Load experimental data
    print("Loading experimental data...")
    df  = load_raw_data(TXT_FILE, CSV_CHANNELS, EXCEL_FILE, SHEET_NAME)
    zt0 = find_first_light_onset(df)
    print(f"  ZT0 detected at: {zt0}")

    # Run model simulations
    sim_results = []
    for v in VARIANTS:
        print(f"  Simulating: {v['title']} ...")
        sim_results.append(run_experiment(v))

    # Figure layout: 5 equal-width columns (2 experimental + 3 model)
    n_exp   = len(TARGET_SPIDERS)   # 2
    n_model = len(VARIANTS)         # 3
    n_cols  = n_exp + n_model       # 5

    fig_h = total_model_days * 0.42 + 2.5 + 0.8
    fig_w = n_cols * 2.0   # equal column width of 2 inches each

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = gridspec.GridSpec(1, n_cols, figure=fig,
                            left=0.08, right=0.95,
                            top=0.85, bottom=0.20,
                            wspace=0.15)

    # Panel A: experimental rasters
    exp_arrows = [ARROW_LCF45, ARROW_LCF44]
    for ci, (spider_id, arrow) in enumerate(zip(TARGET_SPIDERS, exp_arrows)):
        ax = fig.add_subplot(gs[ci])
        draw_experimental_raster(ax, df, spider_id,
                                 zt0, N_DAYS,
                                 show_yaxis=(ci == 0))
        ax.set_title(spider_id, fontsize=FS_TITLE, fontweight='bold', pad=6)
        draw_phase_arrow(ax, *arrow)
        if ci == 0:
            add_panel_label(ax, 'A')

    # Panel B: model simulations
    model_arrows = [ARROW_STANDARD, ARROW_NO_ENTRAIN, ARROW_NO_MASK]
    for mi, (variant, (t_all, y_all, L_all), arrow) in enumerate(
            zip(VARIANTS, sim_results, model_arrows)):
        ci = n_exp + mi
        ax = fig.add_subplot(gs[ci])
        draw_model_raster(ax, t_all, y_all, L_all, total_model_days,
                          show_yaxis=False)
        ax.set_title(f"{variant['title']}\n{variant['subtitle']}",
                     fontsize=FS_TITLE, fontweight='bold', pad=6, linespacing=1.5)
        draw_phase_arrow(ax, *arrow)
        if mi == 1:   # No-entrainment panel: add double-headed arrow
            draw_double_arrow(ax, *DBL_ARROW_NO_ENTRAIN)
        if mi == 0:
            add_panel_label(ax, 'B')
        if mi == n_model - 1:
            add_ld_dd_labels(ax, total_model_days)

    # Shared legend (model panels only)
    line_y = Line2D([0], [0], color='#999999', linewidth=0.9,
                    label='Circadian protein $y$')
    line_L = Line2D([0], [0], color='#C0392B', linewidth=2.0,
                    label='Locomotor activity $L$')
    fig.legend(handles=[line_y, line_L], loc='lower center', ncol=2,
               fontsize=FS_LEGEND, frameon=True, bbox_to_anchor=(0.73, 0.07),
               framealpha=0.9)

    return fig


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == '__main__':
    print('Building Figure 5 v3 ...')
    fig = build_figure()

    save_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(save_dir, 'Fig5.png')
    fig.savefig(out_path, dpi=1000, format='png', bbox_inches='tight')
    print(f'\nSaved -> {out_path}')
    plt.show()
