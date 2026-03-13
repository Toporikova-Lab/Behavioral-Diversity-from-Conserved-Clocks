"""
Three-Species Spider Circadian Activity Model
==============================================

Model Summary:
--------------
This model simulates locomotor activity patterns in three spider species using
a Goodwin oscillator coupled with a locomotor gating mechanism and light masking.

Core Mechanism:
1. Circadian Clock: Goodwin oscillator (x, y, z) with ~24h period
2. Locomotor Gate: Activity-permissive gate (h) controlled by clock protein y
3. Light Entrainment: Light increases degradation of clock protein y
4. Masking: Light directly suppresses locomotor output (M variable)
5. Activity Output: L = (L_baseline + L_amplitude × m × h) × (1 - masking × M)

Species Differences (4 key parameters):
---------------------------------------

1. LIGHT_SENSITIVITY (combined entrainment and perception):
   - How strongly light affects the circadian clock
   - Larinioides: 0.09 (low) - clock minimally affected by light
   - Agelenopsis: 0.06 (very low) - clock nearly insensitive to light
   - Steatoda: 0.75 (high) - clock highly responsive to light
   - Effect: Higher values → stronger entrainment → higher LD amplitude

2. MASKING_STRENGTH:
   - How strongly light suppresses activity (acute behavioral response)
   - Larinioides: 0.9 (strong) - strongly avoids light
   - Agelenopsis: 0.1 (weak) - minimal light avoidance
   - Steatoda: 0.25 (weak) - minimal light avoidance
   - Effect: Higher values → stronger suppression → lower LD amplitude

3. BASELINE_ACTIVITY (L_baseline):
   - Constant activity level independent of circadian clock
   - Larinioides: 0.0 (none) - strictly clock-dependent
   - Agelenopsis: 0.2 (high) - always ready to respond
   - Steatoda: 0.15 (moderate) - some opportunism
   - Effect: Higher values → higher mean activity in all conditions

4. PEAK_WIDTH (tau_h, h_steepness):
   - How broad vs. narrow the daily activity peak is
   - Larinioides: narrow, sharp peaks (tau_h=2, h_steepness=30)
   - Agelenopsis: wide, smooth peaks (tau_h=3, h_steepness=10)
   - Steatoda: wide, smooth peaks (tau_h=3, h_steepness=10)
   - Effect: Wider peaks → extended activity window

Resulting Behaviors:
--------------------
LARINIOIDES (Daily Web-Builder):
- LD/DD ratio: 0.81 (LD < DD due to strong masking)
- Mechanism: Strong masking dominates weak entrainment
- Narrow, sharp activity peaks; no baseline activity
- Nearly arrhythmic in constant light

AGELENOPSIS (Opportunistic Hunter):
- LD/DD ratio: 0.95 (LD ≈ DD, balanced)
- Mechanism: Very weak entrainment balances weak masking
- Wide activity peaks with high baseline
- Strong rhythm maintained in all conditions

STEATODA (Light-Responsive Opportunist):
- LD/DD ratio: 1.33 (LD > DD due to strong entrainment)
- Mechanism: Strong entrainment dominates weak masking
- Wide activity peaks with moderate baseline
- Moderate rhythm in constant light
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# =============================================================================
# GLOBAL SETTINGS
# =============================================================================

# Environmental light intensity (same for all species)
ENVIRONMENTAL_LIGHT = 1.0

# Plotting colors for conditions
COLORS = {'DD': '#2E86AB', 'LD': '#A23B72', 'LL': '#C73E1D'}

# =============================================================================
# SPECIES PARAMETERS
# =============================================================================

# Shared circadian clock parameters (identical across species)
SHARED_CLOCK_PARAMS = {
    'v1': 0.84, 'v2': 0.84, 'v3': 0.84,
    'k1': 0.18, 'k2': 0.18, 'k3': 0.18,
    'n': 12, 'K': 1.0,
    'sigmoid_steepness': 12,
    'L_amplitude': 2.5,
    'tau_m': 2.0,
    'x0': 0.1, 'y0': 0.2, 'z0': 1.5, 'h0': 1.0, 'M0': 0.0
}

# Species-specific parameters
LARINIOIDES = {
    **SHARED_CLOCK_PARAMS,
    'species_name': 'Larinioides',
    'light_sensitivity': 0.09,     # Low - weak clock entrainment
    'masking_strength': 0.9,       # Strong - avoids light
    'L_baseline': 0.0,             # No baseline activity
    'tau_h': 2.0,                  # Fast gate recovery
    'h_steepness': 30,             # Sharp gate threshold
    'y_threshold': 0.5,            # High activity threshold
}

AGELENOPSIS = {
    **SHARED_CLOCK_PARAMS,
    'species_name': 'Agelenopsis',
    'light_sensitivity': 0.06,     # Very low - minimal clock entrainment
    'masking_strength': 0.1,       # Weak - little light avoidance
    'L_baseline': 0.2,             # High baseline activity
    'tau_h': 3.0,                  # Slow gate recovery
    'h_steepness': 10,             # Gradual gate threshold
    'y_threshold': 0.3,            # Low activity threshold
}

STEATODA = {
    **SHARED_CLOCK_PARAMS,
    'species_name': 'Steatoda',
    'light_sensitivity': 0.75,     # High - strong clock entrainment
    'masking_strength': 0.25,      # Weak - little light avoidance
    'L_baseline': 0.15,            # Moderate baseline activity
    'tau_h': 3.0,                  # Slow gate recovery
    'h_steepness': 10,             # Gradual gate threshold
    'y_threshold': 0.3,            # Low activity threshold
}

# =============================================================================
# MODEL FUNCTIONS
# =============================================================================

def light_function(t, light_type):
    """
    Environmental light intensity as function of time.
    Returns 0 (dark) or ENVIRONMENTAL_LIGHT (light).
    """
    if light_type == 'DD':
        return 0.0
    elif light_type == 'LL':
        return ENVIRONMENTAL_LIGHT
    elif light_type == 'LD':
        time_in_cycle = t % 24
        return ENVIRONMENTAL_LIGHT if time_in_cycle < 12 else 0.0
    return 0.0

def model_odes(state, t, params, light_type):
    """
    Differential equations for circadian clock with locomotor gating.
    
    Variables:
    x, y, z: Clock proteins (Goodwin oscillator)
    h: Locomotor gate (activity-permissive state)
    M: Masking signal (follows light on/off)
    """
    x, y, z, h, M = state
    
    # Environmental light
    light = light_function(t, light_type)
    
    # Light-dependent degradation of y (entrainment mechanism)
    k2_eff = params['k2'] * (1 + params['light_sensitivity'] * light)
    
    # Goodwin oscillator (circadian clock)
    dx = (params['v1'] * params['K']**params['n'] / 
          (params['K']**params['n'] + z**params['n']) - 
          params['k1'] * x)
    dy = params['v2'] * x - k2_eff * y
    dz = params['v3'] * y - params['k3'] * z
    
    # Locomotor gate (opens when y is low)
    h_inf = 1 / (1 + np.exp(params['h_steepness'] * 
                           (y - params['y_threshold'])))
    dh = (h_inf - h) / params['tau_h']
    
    # Masking signal (binary: 1 when light on, 0 when light off)
    M_target = 1.0 if light > 0 else 0.0
    dM = (M_target - M) / params['tau_m']
    
    return [dx, dy, dz, dh, dM]

def simulate_species(params, light_type):
    """
    Simulate one species under one light condition.
    Returns time array and activity output.
    """
    # Equilibration (10 days)
    eq_t = np.arange(0, 10 * 24, 0.01)
    eq_state = odeint(
        model_odes,
        [params['x0'], params['y0'], params['z0'], 
         params['h0'], params['M0']],
        eq_t,
        args=(params, light_type)
    )
    
    # Main simulation (32 days for stability)
    t = np.arange(0, 32 * 24, 0.01)
    state = odeint(
        model_odes, eq_state[-1, :], t,
        args=(params, light_type)
    )
    
    # Calculate activity output
    x, y, z, h, M = state.T
    m = 1 / (1 + np.exp(-params['sigmoid_steepness'] * 
                        (y - params['y_threshold'])))
    L_circadian = params['L_baseline'] + params['L_amplitude'] * m * h
    L_output = L_circadian * (1 - params['masking_strength'] * M)
    
    return t, L_output

def extract_cycle_for_plotting(t, L_out, condition, start_day=30):
    """
    Extract one 24-hour cycle for plotting.
    - LD: Natural ZT 0-24 alignment
    - DD/LL: Peak-centered at hour 12
    """
    dt = t[1] - t[0]
    samples_per_day = int(24 / dt)
    
    if condition == 'LD':
        # Natural ZT alignment for LD
        start = start_day * samples_per_day
        end = start + samples_per_day
        
        L_cycle = L_out[start:end]
        t_plot = np.arange(0, 24, dt)
        L_plot = L_cycle
        
        peak_time = t_plot[np.argmax(L_cycle)]
    else:
        # Peak-centered for DD and LL (avoids wrap-around artifacts)
        start = (start_day - 1) * samples_per_day
        end = (start_day + 2) * samples_per_day
        
        L_extended = L_out[start:end]
        
        # Find peak in middle day
        middle_start = samples_per_day
        middle_end = 2 * samples_per_day
        L_middle = L_extended[middle_start:middle_end]
        
        peak_idx_middle = np.argmax(L_middle)
        peak_time = peak_idx_middle * dt
        
        # Extract 24 hours centered on peak
        peak_idx_extended = middle_start + peak_idx_middle
        extract_start = peak_idx_extended - int(12 / dt)
        extract_end = extract_start + samples_per_day
        
        L_plot = L_extended[extract_start:extract_end]
        t_plot = np.arange(0, 24, dt)
    
    return t_plot, L_plot, peak_time

# =============================================================================
# MAIN SIMULATION
# =============================================================================

def run_simulation():
    """Run simulation for all species and conditions."""
    
    print("=" * 80)
    print("THREE-SPECIES SPIDER CIRCADIAN MODEL")
    print("=" * 80)
    print(f"Environmental light intensity: {ENVIRONMENTAL_LIGHT}")
    print("=" * 80)
    
    species_list = [LARINIOIDES, AGELENOPSIS, STEATODA]
    conditions = ['DD', 'LD', 'LL']
    all_results = {}
    
    for params in species_list:
        species_name = params['species_name']
        print(f"\n{species_name}:")
        print(f"  Light sensitivity: {params['light_sensitivity']}")
        print(f"  Masking strength:  {params['masking_strength']}")
        print(f"  Baseline activity: {params['L_baseline']}")
        
        species_results = {}
        
        for condition in conditions:
            t, L_out = simulate_species(params, condition)
            t_plot, L_plot, peak_time = extract_cycle_for_plotting(
                t, L_out, condition)
            
            amp = np.max(L_plot) - np.min(L_plot)
            
            if condition == 'LD':
                print(f"  {condition}: Amp={amp:.3f}, Peak at ZT {peak_time:.1f}")
            else:
                print(f"  {condition}: Amp={amp:.3f}")
            
            species_results[condition] = {
                't': t_plot,
                'L': L_plot,
                'color': COLORS[condition],
                'amplitude': amp
            }
        
        # Calculate LD/DD ratio
        ratio = (species_results['LD']['amplitude'] / 
                species_results['DD']['amplitude'])
        print(f"  LD/DD ratio: {ratio:.3f}")
        
        all_results[species_name] = species_results
    
    return all_results

# =============================================================================
# PLOTTING
# =============================================================================

def create_figure(all_results):
    """Create publication-ready 3x3 comparison figure."""
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 10))
    
    species_names = ['Larinioides', 'Agelenopsis', 'Steatoda']
    conditions = ['DD', 'LD', 'LL']
    
    for col_idx, species_name in enumerate(species_names):
        species_data = all_results[species_name]
        
        for row_idx, condition in enumerate(conditions):
            ax = axes[row_idx, col_idx]
            result = species_data[condition]
            
            # Plot trace
            ax.plot(result['t'], result['L'], 
                   color=result['color'], linewidth=2.5, alpha=0.9)
            
            # Add shading (only dark periods)
            if condition == 'DD':
                ax.axvspan(0, 24, alpha=0.15, color='gray', zorder=0)
            elif condition == 'LD':
                ax.axvspan(12, 24, alpha=0.15, color='gray', zorder=0)
                ax.axvline(12, color='black', linestyle='--', 
                          linewidth=1.5, alpha=0.6)
            
            # Formatting
            ax.set_xlim(0, 24)
            ax.set_ylim(0, 1.5)
            ax.set_xticks([0, 6, 12, 18, 24])
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
            
            # Y-axis (left column only)
            if col_idx == 0:
                ax.set_ylabel('Mean Activity\n(crossings/min)', 
                             fontsize=10, fontweight='bold')
            else:
                ax.set_yticklabels([])
            
            # X-axis labels
            if condition == 'LD':
                ax.set_xlabel('Zeitgeber Time (ZT)', 
                             fontsize=9, fontweight='bold')
            else:
                ax.set_xlabel('Hour (aligned to peak)', 
                             fontsize=9, fontweight='bold')
            
            # Condition label (left side)
            if col_idx == 0:
                ax.text(-0.25, 0.5, condition, 
                       transform=ax.transAxes,
                       fontsize=12, fontweight='bold',
                       color=result['color'],
                       va='center', ha='right')
            
            # Species title (top row)
            # if row_idx == 0:
            #     ax.set_title(species_name, fontsize=13, 
            #                 fontweight='bold', pad=10)
    
    plt.tight_layout(rect=[0.03, 0, 1, 0.97])
    # fig.suptitle('Three Spider Species: Model Activity Patterns',
    #              fontsize=14, fontweight='bold', y=0.995)
    
    return fig

# =============================================================================
# RUN AND SAVE
# =============================================================================

if __name__ == "__main__":
    # Run simulation
    results = run_simulation()
    
    # Create and save figure
    print("\n" + "=" * 80)
    print("Creating figure...")
    fig = create_figure(results)
    
    filename = 'spider_model_comparison.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved: {filename}")
    print("=" * 80)
    
    plt.show()
