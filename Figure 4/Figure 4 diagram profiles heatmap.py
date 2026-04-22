"""
Combined figure for the three-species spider circadian model.

Layout
------
  Panel A (top-left):    Dual-pathway model schematic.
  Panel B (bottom-left): LD/DD amplitude ratio heatmap across
                         (light_sensitivity, masking_strength) parameter space,
                         with species positions overlaid.
  Panel C (right, 3x3):  Simulated activity traces (species columns,
                         DD / LD / LL rows).

Designed to run in Spyder on Windows. The figure is saved as PNG and SVG
in the same directory as this script.

To speed things up, reduce HEATMAP_RESOLUTION (at the cost of a blockier
heatmap). The rest of the simulation is fast.
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle
from matplotlib.colors import LinearSegmentedColormap


# ============================================================================
# PARAMETERS
# ============================================================================

SHARED = {
    'v1': 0.84, 'v2': 0.84, 'v3': 0.84,
    'k1': 0.18, 'k2': 0.18, 'k3': 0.18,
    'n': 12, 'K': 1.0,
    'sigmoid_steepness': 12, 'L_amplitude': 2.5, 'tau_m': 2.0,
    'x0': 0.1, 'y0': 0.2, 'z0': 1.5, 'h0': 1.0, 'M0': 0.0,
}

SPECIES = {
    'L. cornutus': {
        **SHARED, 'light_sensitivity': 0.09, 'masking_strength': 0.9,
        'L_baseline': 0.0, 'tau_h': 2.0, 'h_steepness': 30,
        'y_threshold': 0.5,
        'color': '#D62828', 'marker': 'o',
    },
    'A. pennsylvanica': {
        **SHARED, 'light_sensitivity': 0.06, 'masking_strength': 0.1,
        'L_baseline': 0.2, 'tau_h': 3.0, 'h_steepness': 10,
        'y_threshold': 0.3,
        'color': '#003049', 'marker': 's',
    },
    'S. grossa': {
        **SHARED, 'light_sensitivity': 0.75, 'masking_strength': 0.25,
        'L_baseline': 0.15, 'tau_h': 3.0, 'h_steepness': 10,
        'y_threshold': 0.3,
        'color': '#F77F00', 'marker': '^',
    },
}

COND_COLORS = {'DD': '#2E86AB', 'LD': '#A23B72', 'LL': '#C73E1D'}
CONDITIONS = ['DD', 'LD', 'LL']

HEATMAP_RESOLUTION = 15  # 20-25 gives a smoother heatmap but is slower


# ============================================================================
# MODEL
# ============================================================================

def light_signal(t, condition):
    """Environmental light: 0 (dark) or 1 (light)."""
    if condition == 'DD':
        return 0.0
    if condition == 'LL':
        return 1.0
    return 1.0 if (t % 24) < 12 else 0.0  # LD 12:12, light first


def model_odes(state, t, p, condition):
    """Goodwin oscillator + locomotor gate + masking signal."""
    x, y, z, h, M = state
    light = light_signal(t, condition)
    k2e = p['k2'] * (1 + p['light_sensitivity'] * light)

    dx = p['v1'] * p['K']**p['n'] / (p['K']**p['n'] + z**p['n']) - p['k1'] * x
    dy = p['v2'] * x - k2e * y
    dz = p['v3'] * y - p['k3'] * z

    h_inf = 1.0 / (1.0 + np.exp(p['h_steepness'] * (y - p['y_threshold'])))
    dh = (h_inf - h) / p['tau_h']

    dM = ((1.0 if light > 0 else 0.0) - M) / p['tau_m']
    return [dx, dy, dz, dh, dM]


def equilibrate(p, condition, eq_days=10, verbose=True):
    """Run a 10-day equilibration; return final state, or ICs if flat."""
    t = np.arange(0, eq_days * 24, 0.1)
    init = [p['x0'], p['y0'], p['z0'], p['h0'], p['M0']]
    state = odeint(model_odes, init, t, args=(p, condition))
    y_last = state[-int(2 * 24 / 0.1):, 1]
    if np.max(y_last) - np.min(y_last) < 0.01:
        if verbose:
            print(f"  Note: no oscillations during {condition} "
                  f"equilibration; using initial conditions.")
        return init
    return state[-1].tolist()


def simulate(p, condition, sim_days=32, dt=0.01, verbose=True):
    """Equilibrate + simulate. Returns (t, L) for the activity output."""
    init = equilibrate(p, condition, verbose=verbose)
    t = np.arange(0, sim_days * 24, dt)
    state = odeint(model_odes, init, t, args=(p, condition))
    x, y, z, h, M = state.T
    m = 1.0 / (1.0 + np.exp(-p['sigmoid_steepness']
                            * (y - p['y_threshold'])))
    L = ((p['L_baseline'] + p['L_amplitude'] * m * h)
         * (1.0 - p['masking_strength'] * M))
    return t, L


def extract_cycle(t, L, condition, day=30):
    """One 24h cycle. LD is ZT-aligned; DD/LL are peak-centered at hour 12."""
    dt = t[1] - t[0]
    ppd = int(round(24 / dt))
    if condition == 'LD':
        i0 = day * ppd
        return np.arange(0, 24, dt), L[i0:i0 + ppd]
    i0 = (day - 1) * ppd
    L_ext = L[i0:i0 + 3 * ppd]
    peak_in_mid = int(np.argmax(L_ext[ppd:2 * ppd]))
    s = ppd + peak_in_mid - int(round(12 / dt))
    return np.arange(0, 24, dt), L_ext[s:s + ppd]


# ============================================================================
# PARAMETER SPACE (HEATMAP)
# ============================================================================

def ld_dd_ratio(p, dt=0.1, sim_days=20):
    """LD/DD amplitude ratio for a single parameter set."""
    try:
        _, L_dd = simulate(p, 'DD', sim_days=sim_days, dt=dt, verbose=False)
        _, L_ld = simulate(p, 'LD', sim_days=sim_days, dt=dt, verbose=False)
        i0 = int(15 * 24 / dt)
        amp_dd = np.max(L_dd[i0:]) - np.min(L_dd[i0:])
        amp_ld = np.max(L_ld[i0:]) - np.min(L_ld[i0:])
        return amp_ld / amp_dd if amp_dd > 0.01 else np.nan
    except Exception:
        return np.nan


def compute_heatmap(resolution=HEATMAP_RESOLUTION):
    """Sweep (light_sensitivity, masking_strength); return grid of LD/DD."""
    print(f"Computing {resolution}x{resolution} heatmap "
          f"({resolution**2} parameter combinations)...")
    base = {**SHARED, 'L_baseline': 0.1, 'tau_h': 2.5,
            'h_steepness': 10, 'y_threshold': 0.3}
    light_vals = np.linspace(0.0, 1.0, resolution)
    mask_vals = np.linspace(0.0, 1.0, resolution)
    grid = np.full((resolution, resolution), np.nan)
    for i, mv in enumerate(mask_vals):
        if i % 3 == 0:
            print(f"  row {i + 1}/{resolution}")
        for j, lv in enumerate(light_vals):
            grid[i, j] = ld_dd_ratio({**base,
                                      'light_sensitivity': lv,
                                      'masking_strength': mv})
    return light_vals, mask_vals, grid


# ============================================================================
# SCHEMATIC HELPERS
# ============================================================================

def _arrow(ax, p0, p1, color='black', lw=1.3):
    """Stimulating arrow (p0 -> p1)."""
    ax.annotate('', xy=p1, xytext=p0,
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=lw, mutation_scale=14))


def _tbar(ax, p0, p1, color='black', lw=1.3, bar=0.20):
    """Inhibitory connector: line from p0 to p1 with T-bar at p1."""
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    length = max(np.hypot(dx, dy), 1e-6)
    ux, uy = dx / length, dy / length
    p1s = (p1[0] - 0.15 * ux, p1[1] - 0.15 * uy)
    ax.annotate('', xy=p1s, xytext=p0,
                arrowprops=dict(arrowstyle='-', color=color, lw=lw))
    px, py = -uy, ux
    ax.plot([p1s[0] + bar * px, p1s[0] - bar * px],
            [p1s[1] + bar * py, p1s[1] - bar * py],
            color=color, lw=lw + 0.3, solid_capstyle='round')


def _draw_sun(ax, center, size=0.32):
    """Stylized yellow sun icon at `center`."""
    cx, cy = center
    ax.add_patch(Circle(center, size, facecolor='#FFD93D',
                        edgecolor='#E8A317', linewidth=1.2, zorder=3))
    for ang in np.arange(0, 2 * np.pi, np.pi / 4):
        x0 = cx + (size + 0.06) * np.cos(ang)
        y0 = cy + (size + 0.06) * np.sin(ang)
        x1 = cx + (size + 0.24) * np.cos(ang)
        y1 = cy + (size + 0.24) * np.sin(ang)
        ax.plot([x0, x1], [y0, y1], color='#E8A317', linewidth=1.4)


# ============================================================================
# PANEL DRAWING
# ============================================================================

def draw_schematic(ax):
    """Panel A: the dual-pathway model schematic."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_aspect('equal')
    ax.axis('off')

    # Clock: dotted red ellipse containing X, Y, Z with feedback loop
    ax.add_patch(Ellipse((2.4, 3.6), 3.4, 3.1, fill=False,
                         edgecolor='#D62828', linewidth=1.8, linestyle=':'))
    ax.text(2.4, 5.35, 'Clock', fontsize=11, fontweight='bold',
            color='#D62828', ha='center', style='italic')
    for name, (x, y) in [('X', (1.3, 4.1)),
                         ('Y', (3.3, 3.8)),
                         ('Z', (2.4, 2.3))]:
        ax.text(x, y, name, fontsize=14, fontweight='bold',
                ha='center', va='center')
    _arrow(ax, (1.45, 4.05), (3.1, 3.85))    # X -> Y
    _arrow(ax, (3.2, 3.55), (2.6, 2.5))      # Y -> Z
    _tbar(ax, (2.2, 2.5), (1.3, 3.85))       # Z -| X

    # Locomotion: dotted blue rounded rectangle containing L
    ax.add_patch(FancyBboxPatch((6.0, 3.05), 1.5, 1.3,
                                boxstyle="round,pad=0.1",
                                fill=False, edgecolor='#1F77B4',
                                linewidth=1.8, linestyle=':'))
    ax.text(6.75, 4.55, 'Locomotion', fontsize=11, fontweight='bold',
            color='#1F77B4', ha='center', style='italic')
    ax.text(6.75, 3.7, 'L', fontsize=16, fontweight='bold',
            ha='center', va='center')

    # Y -> L (both stimulation and gated inhibition)
    _arrow(ax, (3.55, 4.0), (6.05, 4.0))
    _tbar(ax, (3.55, 3.55), (6.05, 3.55))

    # Light (sun) and its two pathways
    _draw_sun(ax, (6.75, 1.3))
    ax.text(6.75, 0.45, 'Light', fontsize=11, fontweight='bold',
            color='#E8A317', ha='center', style='italic')
    _tbar(ax, (6.35, 1.6), (3.5, 3.4))       # Light -| Y (entrainment)
    _tbar(ax, (6.75, 1.7), (6.75, 3.0))      # Light -| L (masking)

    # Legend
    _arrow(ax, (0.3, 0.75), (1.1, 0.75))
    ax.text(1.25, 0.75, 'stimulation', fontsize=8, va='center')
    _tbar(ax, (0.3, 0.3), (1.1, 0.3), bar=0.10)
    ax.text(1.25, 0.3, 'inhibition', fontsize=8, va='center')

    ax.text(-0.03, 5.85, 'A', fontsize=15, fontweight='bold')


def draw_heatmap(ax, light_vals, mask_vals, grid):
    """Panel B: LD/DD amplitude ratio heatmap with species overlay."""
    cmap = LinearSegmentedColormap.from_list(
        'ratio', ['#2E86AB', '#FFFFFF', '#D62828'], N=100)
    im = ax.contourf(light_vals, mask_vals, grid,
                     levels=np.linspace(0.5, 1.5, 21),
                     cmap=cmap, extend='both')
    cs = ax.contour(light_vals, mask_vals, grid,
                    levels=[0.8, 1.0, 1.2],
                    colors='black', linewidths=[1, 1.5, 1],
                    linestyles=['--', '-', '--'])
    ax.clabel(cs, inline=True, fontsize=8, fmt='%.2f')

    for name, p in SPECIES.items():
        ax.plot(p['light_sensitivity'], p['masking_strength'],
                marker=p['marker'], markersize=11,
                markerfacecolor=p['color'],
                markeredgecolor='white', markeredgewidth=1.6,
                label=name, linestyle='None', zorder=10)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Light sensitivity\n(entrainment)', fontsize=10)
    ax.set_ylabel('Masking strength\n(suppression)', fontsize=10)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.92,
              handlelength=1.0, borderpad=0.4)

    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03,
                        ticks=[0.5, 0.75, 1.0, 1.25, 1.5])
    cbar.set_label('LD/DD ratio', fontsize=9, labelpad=2)
    cbar.ax.tick_params(labelsize=8)

    ax.text(-0.24, 1.02, 'B', transform=ax.transAxes,
            fontsize=15, fontweight='bold')


def _shade_activity_bg(ax, condition):
    """Add dark-phase shading to an activity subplot."""
    if condition == 'DD':
        ax.axvspan(0, 24, alpha=0.12, color='gray', zorder=0)
    elif condition == 'LD':
        ax.axvspan(12, 24, alpha=0.15, color='gray', zorder=0)
        ax.axvline(12, color='black', linestyle='--',
                   linewidth=1, alpha=0.5)


def draw_activity_grid(axes, results):
    """Panel C: 3x3 grid (species columns, DD/LD/LL rows)."""
    species_names = list(SPECIES.keys())
    for col, name in enumerate(species_names):
        for row, cond in enumerate(CONDITIONS):
            ax = axes[row][col]
            t, L = results[name][cond]
            ax.plot(t, L, color=COND_COLORS[cond], linewidth=2.0)
            _shade_activity_bg(ax, cond)

            ax.set_xlim(0, 24)
            ax.set_ylim(0, 1.5)
            ax.set_xticks([0, 6, 12, 18, 24])
            ax.grid(alpha=0.25, linewidth=0.5)

            if row == 0:
                ax.set_title(name, fontsize=11, fontstyle='italic',
                             fontweight='bold', pad=8)
            if col == 0:
                ax.text(0.03, 0.93, cond, transform=ax.transAxes,
                        fontsize=13, fontweight='bold',
                        color=COND_COLORS[cond],
                        ha='left', va='top',
                        bbox=dict(facecolor='white', edgecolor='none',
                                  alpha=0.75, pad=1.5))
                ax.set_ylabel('Mean activity\n(crossings/min)', fontsize=9)
            else:
                ax.set_yticklabels([])

            xl = ('Zeitgeber time (ZT)' if cond == 'LD'
                  else 'Hour (peak-centered)')
            ax.set_xlabel(xl, fontsize=9)

    axes[0][0].text(-0.50, 1.18, 'C', transform=axes[0][0].transAxes,
                    fontsize=15, fontweight='bold')


# ============================================================================
# MAIN
# ============================================================================

def run_all_simulations():
    """Simulate every species under every condition. Returns a nested dict."""
    print("Simulating species...")
    results = {name: {} for name in SPECIES}
    for name, p in SPECIES.items():
        print(f"  {name}")
        for cond in CONDITIONS:
            t_full, L_full = simulate(p, cond, sim_days=32, dt=0.01)
            results[name][cond] = extract_cycle(t_full, L_full, cond)
    return results


def build_figure(results, light_vals, mask_vals, grid):
    """Assemble the combined multi-panel figure."""
    fig = plt.figure(figsize=(17, 9))
    gs = GridSpec(3, 4, figure=fig,
                  width_ratios=[1.5, 1, 1, 1],
                  height_ratios=[1, 1, 1],
                  hspace=0.55, wspace=0.42,
                  left=0.05, right=0.97, top=0.94, bottom=0.08)

    ax_sch = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[1:, 0])
    act_axes = [[fig.add_subplot(gs[r, c]) for c in (1, 2, 3)]
                for r in range(3)]

    draw_schematic(ax_sch)
    draw_heatmap(ax_heat, light_vals, mask_vals, grid)
    draw_activity_grid(act_axes, results)
    return fig


def main():
    results = run_all_simulations()
    light_vals, mask_vals, grid = compute_heatmap()

    fig = build_figure(results, light_vals, mask_vals, grid)

    for ext in ('png', 'svg'):
        fname = f'Figure4 combined.{ext}'
        fig.savefig(fname, dpi=300, bbox_inches='tight')
        print(f"Saved: {fname}")

    plt.show()
    return fig


if __name__ == '__main__':
    main()
