"""
Fig1-Appendix: Combined Figure with Panels A–E
==================================================

Layout (2×2 outer grid):
  Top-left  [0,0]: A (raster LD 1:1) + B (raster DD) side by side
  Top-right [0,1]: C (internal model variables, 3×2)
  Bot-left  [1,0]: D (Period vs. Light Intensity)
  Bot-right [1,1]: E (Traces of activity L vs time)

Vertical alignment guaranteed by shared column boundaries:
  A, D share left edge of column 0
  C, E share left and right edges of column 1

Panels A, B, C share the same top edge (same outer row).
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.timeseries import LombScargle

# =============================================================================
# EXPERIMENTAL PROTOCOL CONFIGURATION  (Panels A–C)
# =============================================================================

EQUILIBRATION_DAYS   = 20
LD12_BASELINE_DAYS   = 7
TEST_CONDITION_DAYS  = 20
RASTER_LD12_DAYS     = 5
MECHANISM_DAY        = 1
FIRST_PULSE_ZT       = 0

# =============================================================================
# LIGHT-INTENSITY SWEEP CONFIGURATION  (Panels D–E)
# =============================================================================

LIGHT_MIN   = 0.0
LIGHT_MAX   = 10.0
LIGHT_STEP  = 0.5
EXAMPLE_LIGHT_INTENSITIES = [0.0, 0.5, 1.0, 2.0]
EQUIL_DAYS_DE  = 10
ANALYSIS_DAYS  = 10
PERIOD_MIN     = 20.0
PERIOD_MAX     = 24.0
SIGNIFICANCE_THRESHOLD = 0.05

# =============================================================================
# SHARED MODEL PARAMETERS
# =============================================================================

PARAMS = {
    'v1': 0.84, 'v2': 0.84, 'v3': 0.84,
    'k1': 0.18, 'k2': 0.18, 'k3': 0.18,
    'n': 12, 'K': 1.0,
    'sigmoid_steepness': 12,
    'light_sensitivity': 0.09,
    'masking_strength': 0.9,
    'L_baseline': 0.0,
    'L_amplitude': 2.5,
    'tau_h': 2.0,
    'h_steepness': 30,
    'y_threshold': 0.5,
    'tau_m': 2.0,
    'x0': 0.1, 'y0': 0.2, 'z0': 1.5, 'h0': 1.0, 'M0': 0.0
}

ENVIRONMENTAL_LIGHT = 1.0

# =============================================================================
# MODEL – PANELS A–C
# =============================================================================

def light_function_AB(t, test_condition, transition_time):
    if test_condition == 'LD12' or t < transition_time:
        return ENVIRONMENTAL_LIGHT if (t % 24) < 12 else 0.0
    if test_condition == 'DD':
        return 0.0
    if test_condition == 'LD1':
        time_since = t - transition_time
        if time_since < FIRST_PULSE_ZT:
            return 0.0
        return ENVIRONMENTAL_LIGHT if ((time_since - FIRST_PULSE_ZT) % 2) < 1 else 0.0
    return 0.0


def model_odes_AB(state, t, params, test_condition, transition_time):
    x, y, z, h, M = state
    light = light_function_AB(t, test_condition, transition_time)
    k2_eff = params['k2'] * (1 + params['light_sensitivity'] * light)
    dx = (params['v1'] * params['K']**params['n'] /
          (params['K']**params['n'] + z**params['n']) - params['k1'] * x)
    dy = params['v2'] * x - k2_eff * y
    dz = params['v3'] * y - params['k3'] * z
    h_inf = 1 / (1 + np.exp(params['h_steepness'] * (y - params['y_threshold'])))
    dh = (h_inf - h) / params['tau_h']
    M_target = 1.0 if light > 0 else 0.0
    dM = (M_target - M) / params['tau_m']
    return [dx, dy, dz, dh, dM]


def calculate_activity_AB(state, params):
    x, y, z, h, M = state.T
    m = 1 / (1 + np.exp(-params['sigmoid_steepness'] * (y - params['y_threshold'])))
    L_circ = params['L_baseline'] + params['L_amplitude'] * m * h
    return L_circ * (1 - params['masking_strength'] * M)


def run_equilibration_AB(params):
    print(f"Running {EQUILIBRATION_DAYS}-day LD 12:12 equilibration...")
    t_eq = np.arange(0, EQUILIBRATION_DAYS * 24, 0.1)
    state_eq = odeint(model_odes_AB,
                      [params['x0'], params['y0'], params['z0'],
                       params['h0'], params['M0']],
                      t_eq, args=(params, 'LD12', 1e6))
    return state_eq[-1, :]


def run_experiment_AB(initial_state, params, test_condition):
    print(f"  Running {test_condition} experiment...")
    total_hours = (LD12_BASELINE_DAYS + TEST_CONDITION_DAYS) * 24
    transition_time = LD12_BASELINE_DAYS * 24
    t = np.arange(0, total_hours, 0.1)
    state = odeint(model_odes_AB, initial_state, t,
                   args=(params, test_condition, transition_time))
    return t, state


def prepare_raster_data(t, L_output):
    dt = 0.1
    spd = int(24 / dt)
    ld12_start_day = LD12_BASELINE_DAYS - RASTER_LD12_DAYS
    combined_L = np.concatenate([
        L_output[ld12_start_day * spd : LD12_BASELINE_DAYS * spd],
        L_output[LD12_BASELINE_DAYS * spd : (LD12_BASELINE_DAYS + TEST_CONDITION_DAYS) * spd]
    ])
    total_days = RASTER_LD12_DAYS + TEST_CONDITION_DAYS
    return combined_L.reshape(total_days, spd), total_days


def extract_mechanism_day(t, state):
    dt = 0.1
    spd = int(24 / dt)
    start = LD12_BASELINE_DAYS * spd + (MECHANISM_DAY - 1) * spd
    t_day = t[start:start + spd] - t[start]
    return t_day, state[start:start + spd, :]


def calc_mechanism_vars(t, state, params):
    x, y, z, h, M = state.T
    m = 1 / (1 + np.exp(-params['sigmoid_steepness'] * (y - params['y_threshold'])))
    L_circ = params['L_baseline'] + params['L_amplitude'] * m * h
    L_out  = L_circ * (1 - params['masking_strength'] * M)
    return {'t': t, 'y': y, 'm': m, 'h': h, 'M': M, 'L': L_out}

# =============================================================================
# MODEL – PANELS D–E
# =============================================================================

def light_function_DE(t, light_intensity):
    if t < FIRST_PULSE_ZT:
        return 0.0
    return light_intensity if ((t - FIRST_PULSE_ZT) % 2) < 1 else 0.0


def model_odes_DE(state, t, params, light_intensity):
    x, y, z, h, M = state
    light = light_function_DE(t, light_intensity)
    k2_eff = params['k2'] * (1 + params['light_sensitivity'] * light)
    dx = (params['v1'] * params['K']**params['n'] /
          (params['K']**params['n'] + z**params['n']) - params['k1'] * x)
    dy = params['v2'] * x - k2_eff * y
    dz = params['v3'] * y - params['k3'] * z
    h_inf = 1 / (1 + np.exp(params['h_steepness'] * (y - params['y_threshold'])))
    dh = (h_inf - h) / params['tau_h']
    M_target = 1.0 if light > 0 else 0.0
    dM = (M_target - M) / params['tau_m']
    return [dx, dy, dz, dh, dM]


def calculate_activity_DE(state, params):
    x, y, z, h, M = state.T
    m = 1 / (1 + np.exp(-params['sigmoid_steepness'] * (y - params['y_threshold'])))
    L_circ = params['L_baseline'] + params['L_amplitude'] * m * h
    return L_circ * (1 - params['masking_strength'] * M)


def calculate_period_lombscargle(t, signal, min_period, max_period):
    signal_centered = signal - np.mean(signal)
    freq = np.linspace(1.0 / max_period, 1.0 / min_period, 10000)
    ls = LombScargle(t, signal_centered)
    power = ls.power(freq)
    peak_idx = np.argmax(power)
    peak_period = 1.0 / freq[peak_idx]
    fap = ls.false_alarm_probability(power[peak_idx])
    return peak_period, power[peak_idx], fap


def run_simulation_DE(light_intensity, params, return_trace=False):
    total_hours = (EQUIL_DAYS_DE + ANALYSIS_DAYS) * 24
    t = np.arange(0, total_hours, 0.1)
    state = odeint(model_odes_DE,
                   [params['x0'], params['y0'], params['z0'],
                    params['h0'], params['M0']],
                   t, args=(params, light_intensity))
    L = calculate_activity_DE(state, params)
    spd = int(24 / 0.1)
    start = EQUIL_DAYS_DE * spd
    t_an = t[start:] - t[start]
    L_an = L[start:]
    period, power, fap = calculate_period_lombscargle(t_an, L_an, PERIOD_MIN, PERIOD_MAX)
    rhythmic = fap < SIGNIFICANCE_THRESHOLD
    if return_trace:
        return period, power, rhythmic, t_an[:2*spd], L_an[:2*spd]
    return period, power, rhythmic

# =============================================================================
# PLOTTING HELPERS
# =============================================================================

def create_raster_plot(raster_data, total_days, test_condition, ax):
    """Plot raster onto provided axes."""
    dt = 0.1
    time_hours = np.arange(0, 24, dt)
    y_max = np.mean(raster_data) + 3 * np.std(raster_data)

    for day in range(total_days):
        y_offset = (total_days - day - 1) * (y_max + 0.1)
        if day < RASTER_LD12_DAYS:
            ax.fill_between([0, 12], y_offset, y_offset + y_max,
                            facecolor='lightblue', alpha=0.7, zorder=0, linewidth=0)
        else:
            if test_condition == 'LD1':
                ld1_day = day - RASTER_LD12_DAYS
                hrs_since_start = ld1_day * 24
                for hour in range(0, 24):
                    total_h = hrs_since_start + hour
                    if total_h >= FIRST_PULSE_ZT and ((total_h - FIRST_PULSE_ZT) % 2) < 1:
                        ax.fill_between([hour, hour + 1],
                                        y_offset, y_offset + y_max,
                                        facecolor='lightblue', alpha=0.7,
                                        zorder=0, linewidth=0)
            elif test_condition == 'DD':
                ax.fill_between([0, 24], y_offset, y_offset + y_max,
                                facecolor='whitesmoke', alpha=0.5, zorder=0, linewidth=0)

        ax.fill_between(time_hours, y_offset, y_offset + raster_data[day, :],
                        color='black', linewidth=0, zorder=2)
        ax.axhline(y_offset, color='black', linewidth=0.5, alpha=0.5, zorder=3)

    ax.set_xlim(0, 24)
    ax.set_ylim(0, total_days * (y_max + 0.1))
    ax.set_xlabel('Time of Day (hours)', fontsize=12, fontweight='bold')
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_xticklabels(['0', '6', '12', '18', '24'], fontsize=9)

    y_tick_pos = [(total_days - d - 0.5) * (y_max + 0.1) for d in range(total_days)]
    ax.set_yticks(y_tick_pos)
    ax.set_yticklabels(list(range(total_days)), fontsize=9)

    cond_label = 'LD 1:1' if test_condition == 'LD1' else 'DD'
    ax.set_title(f'LD 12:12 → {cond_label}', fontsize=13, fontweight='bold', pad=6)


def plot_mechanism_panels(data_ld1, data_dd, axes):
    """
    Plot internal-variable comparison into a 3×2 array of axes.
    Rows: Clock Protein / Activity Components / Activity
    Cols: LD 1:1 / DD
    """
    hours_since_ld1_start = (MECHANISM_DAY - 1) * 24

    def add_shading(ax, condition):
        if condition == 'LD1':
            for hour in range(0, 24):
                total_h = hours_since_ld1_start + hour
                if total_h >= FIRST_PULSE_ZT and ((total_h - FIRST_PULSE_ZT) % 2) >= 1:
                    ax.axvspan(hour, hour + 1, alpha=0.3, color='gray',
                               zorder=0, linewidth=0)
        else:
            ax.axvspan(0, 24, alpha=0.3, color='gray', zorder=0, linewidth=0)

    # Row 0: Clock Protein (y)
    for col, (data, cond, title) in enumerate([
            (data_ld1, 'LD1', f'LD 1:1 (Day {MECHANISM_DAY})'),
            (data_dd,  'DD',  f'DD (Day {MECHANISM_DAY})')]):
        ax = axes[0][col]
        add_shading(ax, cond)
        ax.plot(data['t'], data['y'], color='#2E86AB', linewidth=1.8, zorder=3,
                label='y')
        ax.axhline(PARAMS['y_threshold'], color='gray', linestyle='--',
                   linewidth=1.2, alpha=0.7, zorder=2, label='threshold')
        ax.fill_between(data['t'], PARAMS['y_threshold'], data['y'],
                        where=(data['y'] > PARAMS['y_threshold']),
                        color='lightblue', alpha=0.5, zorder=1,
                        label='y > threshold')
        ax.set_xlim(0, 24)
        ax.set_ylim(0, 0.80)
        ax.set_title(title, fontsize=9, fontweight='bold')
        if col == 0:
            ax.set_ylabel('Clock Protein', fontsize=9, fontweight='bold')
            ax.legend(loc='upper left', fontsize=8)
        else:
            ax.set_yticklabels([])
            ax.legend(loc='upper left', fontsize=8)
        ax.tick_params(labelbottom=False)

    # Row 1: Activity Components (m, h)
    for col, (data, cond) in enumerate([(data_ld1, 'LD1'), (data_dd, 'DD')]):
        ax = axes[1][col]
        add_shading(ax, cond)
        ax.fill_between(data['t'], 0, data['m'], color='#F4A460', alpha=0.6,
                        zorder=2, label='m')
        ax.fill_between(data['t'], 0, data['h'], color='#CD5C5C', alpha=0.6,
                        zorder=1, label='h')
        ax.set_xlim(0, 24)
        ax.set_ylim(0, 1.0)
        if col == 0:
            ax.set_ylabel('Activity\nComponents', fontsize=9, fontweight='bold')
            ax.legend(loc='upper right', fontsize=8)
        else:
            ax.set_yticklabels([])
            ax.legend(loc='upper right', fontsize=8)
        ax.tick_params(labelbottom=False)

    # Row 2: Activity (L, M)
    for col, (data, cond) in enumerate([(data_ld1, 'LD1'), (data_dd, 'DD')]):
        ax = axes[2][col]
        add_shading(ax, cond)
        ax.plot(data['t'], data['M'], color='#20B2AA', linewidth=1.8,
                linestyle='--', alpha=0.8, zorder=2, label='M')
        ax.plot(data['t'], data['L'], color='#FF4500', linewidth=2.0, zorder=3,
                label='L')
        if col == 0:
            ax.text(0.02, 0.95, r'L = (m × h) × (1 − α × M)',
                    transform=ax.transAxes, fontsize=8, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            ax.set_ylabel('Activity', fontsize=9, fontweight='bold')
            ax.legend(loc='upper right', fontsize=8)
        else:
            ax.set_yticklabels([])
            ax.legend(loc='upper right', fontsize=8)
        ax.set_xlim(0, 24)
        ax.set_ylim(0, 1.5)
        ax.set_xlabel('Zeitgeber Time (ZT)', fontsize=9, fontweight='bold')


def plot_period_vs_light(light_intensities, periods, rhythmic_array, ax):
    """Panel D."""
    rhythmic_mask   = rhythmic_array
    arrhythmic_mask = ~rhythmic_array

    ax.scatter(light_intensities[rhythmic_mask], periods[rhythmic_mask],
               s=80, c='#2E86AB', edgecolor='black', linewidth=1.5,
               zorder=3, alpha=0.9)
    if np.any(arrhythmic_mask):
        ax.scatter(light_intensities[arrhythmic_mask],
                   np.full(np.sum(arrhythmic_mask), PERIOD_MIN),
                   s=80, marker='x', c='red', linewidth=2.5, zorder=3)
    if np.sum(rhythmic_mask) > 1:
        ax.plot(light_intensities[rhythmic_mask], periods[rhythmic_mask],
                color='#2E86AB', linewidth=2.5, alpha=0.5, zorder=2)
    if rhythmic_mask[0]:
        ax.axhline(periods[0], color='gray', linestyle='--',
                   linewidth=1.5, alpha=0.7, zorder=1)

    ax.set_xlabel('Light Intensity', fontsize=11, fontweight='bold')
    ax.set_ylabel('Period (hours)', fontsize=11, fontweight='bold')
    ax.set_xlim(LIGHT_MIN - 0.1, LIGHT_MAX + 0.1)
    ax.set_ylim(PERIOD_MIN - 1, PERIOD_MAX + 1)
    ax.set_title('Period vs. Light Intensity', fontsize=11, fontweight='bold', pad=6)
    ax.grid(True, alpha=0.3, linewidth=0.5)


def plot_example_traces(example_traces, ax):
    """Panel E."""
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(EXAMPLE_LIGHT_INTENSITIES)))
    for i, light in enumerate(sorted(EXAMPLE_LIGHT_INTENSITIES)):
        if light in example_traces:
            t_tr, L_tr = example_traces[light]
            ax.plot(t_tr, L_tr, color=colors[i], linewidth=2.0,
                    label=f'Light = {light:.1f}', alpha=0.8)
    for hour in range(0, 48, 2):
        ax.axvspan(hour, hour + 1, alpha=0.15, color='lightblue',
                   zorder=0, linewidth=0)
    ax.set_xlabel('Time (hours)', fontsize=11, fontweight='bold')
    ax.set_ylabel('L (activity)', fontsize=11, fontweight='bold')
    ax.set_xlim(0, 48)
    ax.set_title('Traces of activity L vs time', fontsize=11, fontweight='bold', pad=6)
    ax.legend(loc='upper right', bbox_to_anchor=(0.80, 1.0),
              fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linewidth=0.5)

# =============================================================================
# MAIN EXECUTION
# =============================================================================

print("=" * 70)
print("RUNNING SIMULATIONS FOR PANELS A–C (LD 1:1 vs DD rasters/mechanism)")
print("=" * 70)

initial_state = run_equilibration_AB(PARAMS)
results_AB = {}
full_data_AB = {}

for cond in ['LD1', 'DD']:
    t, state = run_experiment_AB(initial_state, PARAMS, cond)
    L_out = calculate_activity_AB(state, PARAMS)
    raster_data, total_days = prepare_raster_data(t, L_out)
    t_day, state_day = extract_mechanism_day(t, state)
    mech = calc_mechanism_vars(t_day, state_day, PARAMS)
    results_AB[cond]   = {'raster_data': raster_data, 'total_days': total_days}
    full_data_AB[cond] = mech

print("\n" + "=" * 70)
print("RUNNING SIMULATIONS FOR PANELS D–E (Period vs light intensity)")
print("=" * 70)

light_intensities = np.arange(LIGHT_MIN, LIGHT_MAX + LIGHT_STEP, LIGHT_STEP)
periods = []
rhythmic_list = []
example_traces = {}

for light in light_intensities:
    if light in EXAMPLE_LIGHT_INTENSITIES:
        period, power, rhythmic, t_tr, L_tr = run_simulation_DE(
            light, PARAMS, return_trace=True)
        example_traces[light] = (t_tr, L_tr)
    else:
        period, power, rhythmic = run_simulation_DE(light, PARAMS)
    periods.append(period if rhythmic else np.nan)
    rhythmic_list.append(rhythmic)
    status = "Rhythmic" if rhythmic else "Arrhythmic"
    print(f"  Light {light:.1f}: Period = {period:.2f}h ({status})")

periods        = np.array(periods)
rhythmic_array = np.array(rhythmic_list)

# =============================================================================
# BUILD COMBINED FIGURE
# =============================================================================

print("\nBuilding combined figure...")

fig = plt.figure(figsize=(18, 14))

# -----------------------------------------------------------------------
# Outer 2×2 grid
#   Column 0 (left):  A+B rasters (top), D period plot (bottom)
#   Column 1 (right): C mechanism  (top), E trace plot  (bottom)
#
# This guarantees:
#   • A, B, C share the same top edge  (same outer row)
#   • A left edge == D left edge        (same outer column 0)
#   • C left/right edges == E left/right edges  (same outer column 1)
# -----------------------------------------------------------------------

gs_outer = gridspec.GridSpec(
    2, 2, figure=fig,
    width_ratios=[1, 1],
    height_ratios=[3, 1.5],
    hspace=0.18,       # vertical gap between top and bottom rows
    wspace=0.20        # horizontal gap between left and right columns
)

# ---- Top-left: A and B (two rasters, side by side) ----
gs_AB = gridspec.GridSpecFromSubplotSpec(
    1, 2, subplot_spec=gs_outer[0, 0],
    wspace=0.12        # tight gap between A and B
)
ax_A = fig.add_subplot(gs_AB[0])
ax_B = fig.add_subplot(gs_AB[1])

# ---- Top-right: C (3×2 mechanism sub-grid) ----
gs_C = gridspec.GridSpecFromSubplotSpec(
    3, 2, subplot_spec=gs_outer[0, 1],
    hspace=0.06,
    wspace=0.05
)
axes_C = [[fig.add_subplot(gs_C[r, c]) for c in range(2)] for r in range(3)]

# ---- Bottom-left: D ----
ax_D = fig.add_subplot(gs_outer[1, 0])

# ---- Bottom-right: E ----
ax_E = fig.add_subplot(gs_outer[1, 1])

# ---- Populate panels ----

create_raster_plot(results_AB['LD1']['raster_data'],
                   results_AB['LD1']['total_days'], 'LD1', ax_A)
ax_A.set_ylabel('Day', fontsize=12, fontweight='bold')

create_raster_plot(results_AB['DD']['raster_data'],
                   results_AB['DD']['total_days'], 'DD', ax_B)
ax_B.set_yticklabels([])

plot_mechanism_panels(full_data_AB['LD1'], full_data_AB['DD'], axes_C)

plot_period_vs_light(light_intensities, periods, rhythmic_array, ax_D)

plot_example_traces(example_traces, ax_E)

# -----------------------------------------------------------------------
# Panel labels A–E
# Placed well above each panel using figure-level coordinates so that
# A, B, and C are at the same absolute height regardless of inner layout.
# -----------------------------------------------------------------------

# Force layout so we can read axes positions
fig.tight_layout(pad=1.0)
fig.canvas.draw()

def add_panel_label(ax, label, fig,
                    x_offset=-0.08, y_offset=0.04,
                    fontsize=16):
    """
    Place a bold panel label in figure coordinates, anchored to the
    top-left corner of `ax`, shifted up by y_offset (figure fraction).
    """
    pos = ax.get_position()          # Bbox in figure fraction [0,1]
    fig.text(pos.x0 + x_offset * pos.width,
             pos.y1 + y_offset,
             label,
             fontsize=fontsize, fontweight='bold',
             va='bottom', ha='left',
             transform=fig.transFigure)

# A and B: anchor to their own axes
add_panel_label(ax_A, 'A', fig, x_offset=-0.10, y_offset=0.018)
add_panel_label(ax_B, 'B', fig, x_offset=-0.08, y_offset=0.018)

# C: anchor to top-left sub-axis of the mechanism block
add_panel_label(axes_C[0][0], 'C', fig, x_offset=-0.10, y_offset=0.018)

# D and E
add_panel_label(ax_D, 'D', fig, x_offset=-0.08, y_offset=0.018)
add_panel_label(ax_E, 'E', fig, x_offset=-0.06, y_offset=0.018)

# ---- Save ----
output_file = 'Fig1-Appendix.png'
plt.savefig(output_file, dpi=1000, format='png', bbox_inches='tight')
plt.show()
