"""
Larinioides Masking Experiment: LD 1:1 vs DD Comparison
========================================================

Protocol (both conditions):
- EQUILIBRATION_DAYS of LD 12:12 (equilibration)
- LD12_BASELINE_DAYS of LD 12:12 (experimental baseline)
- TEST_CONDITION_DAYS of test condition (LD 1:1 or DD)

Left subplot: LD 1:1 (ultradian masking test)
Right subplot: DD (free-running control)

Raster plots show the last RASTER_LD12_DAYS of LD 12:12 and all days of test.
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# =============================================================================
# EXPERIMENTAL PROTOCOL CONFIGURATION
# =============================================================================

# Control experimental durations here (in days)
EQUILIBRATION_DAYS = 20  # Initial LD 12:12 equilibration (fixed)
LD12_BASELINE_DAYS = 7   # LD 12:12 baseline after equilibration
TEST_CONDITION_DAYS = 20 # Duration of LD 1:1 or DD test condition

# Raster plot display settings
RASTER_LD12_DAYS = 5     # Number of LD12 baseline days to show in raster

# Mechanism figure settings
MECHANISM_DAY = 1        # Which day of test condition to plot (1 to TEST_CONDITION_DAYS)
                         # 1 = first day, 20 = last day, etc.

# LD 1:1 light pulse timing control
FIRST_PULSE_ZT = 0       # ZT time (0-24) for first light pulse in LD 1:1
                         # 0 = light pulse starts at transition
                         # 12 = light pulse starts 12h after transition (in dark phase)

# =============================================================================
# PARAMETERS
# =============================================================================

# Larinioides parameters
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
# MODEL EQUATIONS
# =============================================================================

def light_function(t, test_condition, transition_time):
    """
    Light intensity as function of time.
    
    Parameters:
    - test_condition: 'LD12', 'LD1' (ultradian), or 'DD' (constant darkness)
    - transition_time: time when protocol switches from LD12 to test
    """
    if test_condition == 'LD12' or t < transition_time:
        # LD 12:12 (baseline or equilibration)
        time_in_cycle = t % 24
        return ENVIRONMENTAL_LIGHT if time_in_cycle < 12 else 0.0
    else:
        # Test condition
        if test_condition == 'DD':
            return 0.0
        elif test_condition == 'LD1':
            # LD 1:1 (ultradian) - wait FIRST_PULSE_ZT hours before first pulse
            time_since_transition = t - transition_time
            
            # No light until FIRST_PULSE_ZT hours have passed
            if time_since_transition < FIRST_PULSE_ZT:
                return 0.0
            
            # After delay, start 1:1 cycles
            time_in_cycle = (time_since_transition - FIRST_PULSE_ZT) % 2
            return ENVIRONMENTAL_LIGHT if time_in_cycle < 1 else 0.0
    return 0.0

def model_odes(state, t, params, test_condition, transition_time):
    """Goodwin oscillator with locomotor gating and masking."""
    x, y, z, h, M = state
    
    # Get current light level
    light = light_function(t, test_condition, transition_time)
    
    # Light-dependent degradation (entrainment)
    k2_eff = params['k2'] * (1 + params['light_sensitivity'] * light)
    
    # Circadian oscillator
    dx = (params['v1'] * params['K']**params['n'] / 
          (params['K']**params['n'] + z**params['n']) - 
          params['k1'] * x)
    dy = params['v2'] * x - k2_eff * y
    dz = params['v3'] * y - params['k3'] * z
    
    # Locomotor gate (opens when y is low)
    h_inf = 1 / (1 + np.exp(params['h_steepness'] * 
                           (y - params['y_threshold'])))
    dh = (h_inf - h) / params['tau_h']
    
    # Masking signal
    M_target = 1.0 if light > 0 else 0.0
    dM = (M_target - M) / params['tau_m']
    
    return [dx, dy, dz, dh, dM]

# =============================================================================
# SIMULATION FUNCTIONS
# =============================================================================

def run_equilibration(params):
    """Run LD 12:12 equilibration for configured duration."""
    print(f"Running {EQUILIBRATION_DAYS}-day LD 12:12 equilibration...")
    
    t_eq = np.arange(0, EQUILIBRATION_DAYS * 24, 0.1)
    state_eq = odeint(
        model_odes,
        [params['x0'], params['y0'], params['z0'], 
         params['h0'], params['M0']],
        t_eq,
        args=(params, 'LD12', 1e6)  # Large transition time = always LD12
    )
    
    return state_eq[-1, :]

def run_experiment(initial_state, params, test_condition):
    """
    Run experimental protocol:
    - LD12_BASELINE_DAYS of LD 12:12
    - TEST_CONDITION_DAYS of test condition (LD1 or DD)
    """
    print(f"Running experimental protocol with {test_condition}...")
    print(f"  Days 0-{LD12_BASELINE_DAYS}: LD 12:12")
    print(f"  Days {LD12_BASELINE_DAYS}-{LD12_BASELINE_DAYS + TEST_CONDITION_DAYS}: {test_condition}")
    
    # Total time
    total_hours = (LD12_BASELINE_DAYS + TEST_CONDITION_DAYS) * 24
    transition_time = LD12_BASELINE_DAYS * 24  # Switch to test condition
    
    t = np.arange(0, total_hours, 0.1)
    
    state = odeint(
        model_odes,
        initial_state,
        t,
        args=(params, test_condition, transition_time)
    )
    
    return t, state

def calculate_activity(state, params):
    """Calculate locomotor activity from model state."""
    x, y, z, h, M = state.T
    
    # Gate modulation
    m = 1 / (1 + np.exp(-params['sigmoid_steepness'] * 
                        (y - params['y_threshold'])))
    
    # Circadian component
    L_circadian = params['L_baseline'] + params['L_amplitude'] * m * h
    
    # Apply masking
    L_output = L_circadian * (1 - params['masking_strength'] * M)
    
    return L_output

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def prepare_raster_data(t, L_output, transition_time):
    """
    Prepare data for raster plot.
    Shows last RASTER_LD12_DAYS of LD 12:12 + all TEST_CONDITION_DAYS.
    """
    dt = 0.1
    samples_per_hour = int(1 / dt)
    samples_per_day = 24 * samples_per_hour
    
    # Define plotting window
    ld12_start_day = LD12_BASELINE_DAYS - RASTER_LD12_DAYS  # Last N days of LD12
    test_start_day = LD12_BASELINE_DAYS  # Start of test condition
    test_end_day = LD12_BASELINE_DAYS + TEST_CONDITION_DAYS
    
    # Extract LD12 days (last RASTER_LD12_DAYS)
    ld12_start_idx = ld12_start_day * samples_per_day
    ld12_end_idx = test_start_day * samples_per_day
    
    # Extract test condition days
    test_start_idx = test_start_day * samples_per_day
    test_end_idx = test_end_day * samples_per_day
    
    # Combine data
    combined_L = np.concatenate([
        L_output[ld12_start_idx:ld12_end_idx],
        L_output[test_start_idx:test_end_idx]
    ])
    
    # Reshape into days x hours
    total_days = RASTER_LD12_DAYS + TEST_CONDITION_DAYS
    raster_data = combined_L.reshape(total_days, samples_per_day)
    
    return raster_data, total_days

def create_raster_plot(raster_data, total_days, test_condition, ax):
    """Create raster plot for one condition."""
    samples_per_day = raster_data.shape[1]
    dt = 0.1
    time_hours = np.arange(0, 24, dt)
    
    # Calculate y-axis range (mean + 3 SD)
    mean_activity = np.mean(raster_data)
    std_activity = np.std(raster_data)
    y_max = mean_activity + 3 * std_activity
    
    for day in range(total_days):
        # Vertical offset for each day
        y_offset = (total_days - day - 1) * (y_max + 0.1)
        
        # Draw shading based on light schedule
        if day < RASTER_LD12_DAYS:
            # LD 12:12 days - light from ZT 0-12 (blue), dark from ZT 12-24 (white)
            # Light period (0-12 hours)
            ax.fill_between([0, 12], y_offset, y_offset + y_max,
                           facecolor='lightblue', alpha=0.7, zorder=0, linewidth=0)
            # Dark period (12-24 hours) - leave white
        else:
            # Test condition days
            if test_condition == 'LD1':
                # LD 1:1 - determine where pulses occur based on FIRST_PULSE_ZT
                # Calculate which day of LD 1:1 this is (0-indexed)
                ld1_day = day - RASTER_LD12_DAYS
                hours_since_ld1_start = ld1_day * 24
                
                # Draw light pulses for each hour of this day
                for hour in range(0, 24):
                    total_hours = hours_since_ld1_start + hour
                    
                    # Check if we're past the initial delay
                    if total_hours >= FIRST_PULSE_ZT:
                        # Calculate position in 2-hour cycle
                        time_in_cycle = (total_hours - FIRST_PULSE_ZT) % 2
                        if time_in_cycle < 1:
                            ax.fill_between([hour, hour + 1], y_offset, y_offset + y_max,
                                           facecolor='lightblue', alpha=0.7, zorder=0, linewidth=0)
            elif test_condition == 'DD':
                # DD - all dark, use light gray background for visibility
                ax.fill_between([0, 24], y_offset, y_offset + y_max,
                               facecolor='whitesmoke', alpha=0.5, zorder=0, linewidth=0)
        
        # Plot activity on top
        activity = raster_data[day, :]
        ax.fill_between(time_hours, y_offset, y_offset + activity,
                       color='black', linewidth=0, zorder=2)
        
        # Add horizontal line at bottom of each day
        ax.axhline(y_offset, color='black', linewidth=0.5, alpha=0.5, zorder=3)
    
    # Formatting
    ax.set_xlim(0, 24)
    ax.set_ylim(0, total_days * (y_max + 0.1))
    ax.set_xlabel('Time of Day (hours)', fontsize=11, fontweight='bold')
    
    # Set x-ticks
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_xticklabels(['00:00', '06:00', '12:00', '18:00', '24:00'])
    
    # Set y-ticks (reversed day numbers)
    y_tick_positions = [(total_days - day - 0.5) * (y_max + 0.1) 
                        for day in range(total_days)]
    y_tick_labels = list(range(total_days))
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels(y_tick_labels)
    
    # Add condition label
    condition_label = 'LD 1:1' if test_condition == 'LD1' else 'DD'
    title = f'LD 12:12 → {condition_label}'
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    
    return y_max

def extract_mechanism_day_data(t, state, test_condition):
    """Extract data from specified day (MECHANISM_DAY) of test condition for mechanism analysis."""
    dt = 0.1
    samples_per_day = int(24 / dt)
    
    # Skip the LD12 baseline days - MECHANISM_DAY refers to day of TEST condition
    # Day 1 of test = day (LD12_BASELINE_DAYS + 1) overall
    baseline_samples = LD12_BASELINE_DAYS * samples_per_day
    
    # Calculate which samples to extract based on MECHANISM_DAY
    # Add baseline_samples to skip the LD12 days
    start_sample = baseline_samples + (MECHANISM_DAY - 1) * samples_per_day
    end_sample = start_sample + samples_per_day
    
    # Validate day number (must be within test condition days)
    if MECHANISM_DAY < 1 or MECHANISM_DAY > TEST_CONDITION_DAYS:
        raise ValueError(f"MECHANISM_DAY must be between 1 and {TEST_CONDITION_DAYS}")
    
    # Get the specified day
    t_day = t[start_sample:end_sample]
    t_day = t_day - t_day[0]  # Reset to 0-24 hours
    state_day = state[start_sample:end_sample, :]
    
    return t_day, state_day

def calculate_mechanism_variables(t, state, params, test_condition):
    """Calculate all internal model variables."""
    x, y, z, h, M = state.T
    
    # Gate modulation
    m = 1 / (1 + np.exp(-params['sigmoid_steepness'] * 
                        (y - params['y_threshold'])))
    
    # Activity output
    L_circadian = params['L_baseline'] + params['L_amplitude'] * m * h
    L_output = L_circadian * (1 - params['masking_strength'] * M)
    
    return {'t': t, 'y': y, 'm': m, 'h': h, 'M': M, 'L': L_output}

def create_mechanism_figure(data_ld1, data_dd):
    """Create side-by-side comparison of internal variables (Figure 5 style)."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
    
    # Calculate total hours since LD 1:1 started for light shading
    # MECHANISM_DAY = 1 means first day (0 complete days elapsed)
    hours_since_ld1_start = (MECHANISM_DAY - 1) * 24
    
    # Helper function to add light/dark shading
    def add_light_shading(ax, condition='LD1'):
        if condition == 'LD1':
            # Add gray shading for dark periods, white for light
            for hour in range(0, 24):
                total_hours = hours_since_ld1_start + hour
                if total_hours >= FIRST_PULSE_ZT:
                    time_in_cycle = (total_hours - FIRST_PULSE_ZT) % 2
                    # Dark phase (hour 1-2 of each 2h cycle)
                    if time_in_cycle >= 1:
                        ax.axvspan(hour, hour + 1, alpha=0.3, color='gray', 
                                  zorder=0, linewidth=0)
        else:  # DD
            # All gray for constant darkness
            ax.axvspan(0, 24, alpha=0.3, color='gray', zorder=0, linewidth=0)
    
    # =========================================================================
    # PANEL 1: Clock Protein (y) with threshold
    # =========================================================================
    
    # Left: LD 1:1
    ax = axes[0, 0]
    add_light_shading(ax, 'LD1')
    
    # Plot y
    ax.plot(data_ld1['t'], data_ld1['y'], color='#2E86AB', linewidth=2, 
            label='y', zorder=3)
    
    # Add threshold line
    ax.axhline(PARAMS['y_threshold'], color='gray', linestyle='--', 
              linewidth=1.5, alpha=0.7, label='threshold', zorder=2)
    
    # Shade area above threshold
    ax.fill_between(data_ld1['t'], PARAMS['y_threshold'], data_ld1['y'], 
                    where=(data_ld1['y'] > PARAMS['y_threshold']),
                    color='lightblue', alpha=0.5, zorder=1, label='y > threshold')
    
    # Formatting
    ax.set_xlim(0, 24)
    ax.set_ylabel('Clock Protein', fontsize=11, fontweight='bold')
    ax.set_ylim(0, 0.80)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_title(f'LD 1:1 (Day {MECHANISM_DAY})', fontsize=12, fontweight='bold')
    
    # Right: DD
    ax = axes[0, 1]
    add_light_shading(ax, 'DD')
    
    # Plot y
    ax.plot(data_dd['t'], data_dd['y'], color='#2E86AB', linewidth=2, 
            label='y', zorder=3)
    
    # Add threshold line
    ax.axhline(PARAMS['y_threshold'], color='gray', linestyle='--', 
              linewidth=1.5, alpha=0.7, zorder=2)
    
    # Shade area above threshold
    ax.fill_between(data_dd['t'], PARAMS['y_threshold'], data_dd['y'], 
                    where=(data_dd['y'] > PARAMS['y_threshold']),
                    color='lightblue', alpha=0.5, zorder=1)
    
    # Formatting
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 0.80)
    ax.set_yticklabels([])
    ax.set_title(f'DD (Day {MECHANISM_DAY})', fontsize=12, fontweight='bold')
    
    # =========================================================================
    # PANEL 2: Activity Components (m and h)
    # =========================================================================
    
    # Left: LD 1:1
    ax = axes[1, 0]
    add_light_shading(ax, 'LD1')
    
    # Fill areas under curves
    ax.fill_between(data_ld1['t'], 0, data_ld1['m'], 
                    color='#F4A460', alpha=0.6, label='m', zorder=2)
    ax.fill_between(data_ld1['t'], 0, data_ld1['h'], 
                    color='#CD5C5C', alpha=0.6, label='h', zorder=1)
    
    # Formatting
    ax.set_xlim(0, 24)
    ax.set_ylabel('Activity Components', fontsize=11, fontweight='bold')
    ax.set_ylim(0, 1.0)
    ax.legend(loc='upper right', fontsize=9)
    
    # Right: DD
    ax = axes[1, 1]
    add_light_shading(ax, 'DD')
    
    # Fill areas under curves
    ax.fill_between(data_dd['t'], 0, data_dd['m'], 
                    color='#F4A460', alpha=0.6, zorder=2)
    ax.fill_between(data_dd['t'], 0, data_dd['h'], 
                    color='#CD5C5C', alpha=0.6, zorder=1)
    
    # Formatting
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 1.0)
    ax.set_yticklabels([])
    
    # =========================================================================
    # PANEL 3: Activity (L and M)
    # =========================================================================
    
    # Left: LD 1:1
    ax = axes[2, 0]
    add_light_shading(ax, 'LD1')
    
    # Plot M (dashed)
    ax.plot(data_ld1['t'], data_ld1['M'], color='#20B2AA', linewidth=2, 
            linestyle='--', label='M', alpha=0.8, zorder=2)
    
    # Plot L (solid)
    ax.plot(data_ld1['t'], data_ld1['L'], color='#FF4500', linewidth=2.5, 
            label='L', zorder=3)
    
    # Add formula as text
    ax.text(0.02, 0.95, r'L = (m × h) × (1 - α × M)', 
           transform=ax.transAxes, fontsize=10, 
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Formatting
    ax.set_xlim(0, 24)
    ax.set_ylabel('Activity', fontsize=11, fontweight='bold')
    ax.set_ylim(0, 1.5)
    ax.set_xlabel('Zeitgeber Time (ZT)', fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    
    # Right: DD
    ax = axes[2, 1]
    add_light_shading(ax, 'DD')
    
    # Plot M (dashed)
    ax.plot(data_dd['t'], data_dd['M'], color='#20B2AA', linewidth=2, 
            linestyle='--', alpha=0.8, zorder=2)
    
    # Plot L (solid)
    ax.plot(data_dd['t'], data_dd['L'], color='#FF4500', linewidth=2.5, 
            zorder=3)
    
    # Formatting
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 1.5)
    ax.set_yticklabels([])
    ax.set_xlabel('Zeitgeber Time (ZT)', fontsize=11, fontweight='bold')
    
    # Overall title
    fig.suptitle('Larinioides: Internal Model Variables', 
                 fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    return fig

# =============================================================================
# MAIN EXECUTION
# =============================================================================

print("=" * 70)
print("LARINIOIDES MASKING EXPERIMENT: LD 1:1 vs DD")
print("=" * 70)
print(f"Configuration:")
print(f"  Equilibration: {EQUILIBRATION_DAYS} days LD 12:12")
print(f"  Baseline: {LD12_BASELINE_DAYS} days LD 12:12")
print(f"  Test: {TEST_CONDITION_DAYS} days")
print(f"  First LD 1:1 pulse at ZT: {FIRST_PULSE_ZT}")
print("=" * 70)

# Step 1: Equilibration (shared for both conditions)
initial_state = run_equilibration(PARAMS)
print(f"Equilibration complete. Final state: {initial_state}\n")

# Step 2: Run both experiments and save full data
results = {}
full_data = {}

for test_condition in ['LD1', 'DD']:
    print()
    t, state = run_experiment(initial_state, PARAMS, test_condition)
    L_output = calculate_activity(state, PARAMS)
    transition_time = LD12_BASELINE_DAYS * 24
    raster_data, total_days = prepare_raster_data(t, L_output, transition_time)
    
    # Extract last day for mechanism analysis
    t_last, state_last = extract_mechanism_day_data(t, state, test_condition)
    mech_data = calculate_mechanism_variables(t_last, state_last, PARAMS, test_condition)
    
    results[test_condition] = {
        'raster_data': raster_data,
        'total_days': total_days
    }
    full_data[test_condition] = mech_data

# Step 3: Create side-by-side raster plots
print("\nCreating raster comparison figure...")
fig1, axes = plt.subplots(1, 2, figsize=(14, 12))

# Left subplot: LD 1:1
create_raster_plot(results['LD1']['raster_data'], 
                   results['LD1']['total_days'], 
                   'LD1', axes[0])
axes[0].set_ylabel('Day', fontsize=11, fontweight='bold')

# Right subplot: DD
create_raster_plot(results['DD']['raster_data'], 
                   results['DD']['total_days'], 
                   'DD', axes[1])
axes[1].set_yticklabels([])  # Remove y-tick labels from right subplot

# Overall title
fig1.suptitle('Larinioides: Masking (LD 1:1) vs Free-Running (DD)', 
             fontsize=14, fontweight='bold', y=0.995)

plt.tight_layout(rect=[0, 0, 1, 0.99])

# Save raster figure
filename1 = 'larinioides_masking_raster.png'
plt.savefig(filename1, dpi=300, bbox_inches='tight')
print(f"Saved: {filename1}")

# Step 4: Create mechanism analysis figure
print("\nCreating mechanism analysis figure...")
fig2 = create_mechanism_figure(full_data['LD1'], full_data['DD'])

# Save mechanism figure
filename2 = 'larinioides_mechanism_variables.png'
plt.savefig(filename2, dpi=300, bbox_inches='tight')
print(f"Saved: {filename2}")

print("=" * 70)
print("Both figures created successfully!")
print("=" * 70)

plt.show()