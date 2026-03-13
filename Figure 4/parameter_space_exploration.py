"""
Parameter Space Exploration: LD/DD Ratio Heatmaps
==================================================

This script explores how different parameter combinations affect the LD/DD
amplitude ratio. We create pairwise heatmaps showing all combinations of
the 4 key parameters, with the three spider species overlaid as points.

Goal: Identify which parameter pairs best explain species differences
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# =============================================================================
# SPECIES PARAMETERS (from final model)
# =============================================================================

SPECIES_DATA = {
    'Larinioides': {
        'light_sensitivity': 0.09,
        'masking_strength': 0.9,
        'baseline_activity': 0.0,
        'tau_h': 2.0,
        'color': '#D62828',  # Red
        'marker': 'o'
    },
    'Agelenopsis': {
        'light_sensitivity': 0.06,
        'masking_strength': 0.1,
        'baseline_activity': 0.2,
        'tau_h': 3.0,
        'color': '#003049',  # Dark blue
        'marker': 's'
    },
    'Steatoda': {
        'light_sensitivity': 0.75,
        'masking_strength': 0.25,
        'baseline_activity': 0.15,
        'tau_h': 3.0,
        'color': '#F77F00',  # Orange
        'marker': '^'
    }
}

# Base parameters (shared across all species)
BASE_PARAMS = {
    'v1': 0.84, 'v2': 0.84, 'v3': 0.84,
    'k1': 0.18, 'k2': 0.18, 'k3': 0.18,
    'n': 12, 'K': 1.0,
    'sigmoid_steepness': 12,
    'L_amplitude': 2.5,
    'tau_m': 2.0,
    'h_steepness': 10,
    'y_threshold': 0.3,
    'x0': 0.1, 'y0': 0.2, 'z0': 1.5, 'h0': 1.0, 'M0': 0.0
}

# =============================================================================
# MODEL SIMULATION FUNCTIONS
# =============================================================================

def light_function(t, light_type):
    """Environmental light intensity."""
    if light_type == 'DD':
        return 0.0
    elif light_type == 'LL':
        return 1.0
    elif light_type == 'LD':
        return 1.0 if (t % 24) < 12 else 0.0
    return 0.0

def model_odes(state, t, params, light_type):
    """Model ODEs."""
    x, y, z, h, M = state
    light = light_function(t, light_type)
    
    k2_eff = params['k2'] * (1 + params['light_sensitivity'] * light)
    
    dx = (params['v1'] * params['K']**params['n'] / 
          (params['K']**params['n'] + z**params['n']) - 
          params['k1'] * x)
    dy = params['v2'] * x - k2_eff * y
    dz = params['v3'] * y - params['k3'] * z
    
    h_inf = 1 / (1 + np.exp(params['h_steepness'] * 
                           (y - params['y_threshold'])))
    dh = (h_inf - h) / params['tau_h']
    
    M_target = 1.0 if light > 0 else 0.0
    dM = (M_target - M) / params['tau_m']
    
    return [dx, dy, dz, dh, dM]

def calculate_ld_dd_ratio(params):
    """
    Calculate LD/DD amplitude ratio for given parameters.
    Returns ratio, or np.nan if simulation fails.
    """
    try:
        # Equilibration
        eq_t = np.arange(0, 10 * 24, 0.1)
        eq_state = odeint(
            model_odes,
            [params['x0'], params['y0'], params['z0'], 
             params['h0'], params['M0']],
            eq_t,
            args=(params, 'DD')
        )
        
        # DD simulation
        t_dd = np.arange(0, 20 * 24, 0.1)
        state_dd = odeint(model_odes, eq_state[-1, :], t_dd, 
                         args=(params, 'DD'))
        
        # LD simulation
        t_ld = np.arange(0, 20 * 24, 0.1)
        state_ld = odeint(model_odes, eq_state[-1, :], t_ld, 
                         args=(params, 'LD'))
        
        # Calculate activity
        for state, condition in [(state_dd, 'DD'), (state_ld, 'LD')]:
            x, y, z, h, M = state.T
            m = 1 / (1 + np.exp(-params['sigmoid_steepness'] * 
                                (y - params['y_threshold'])))
            L_circ = params['baseline_activity'] + params['L_amplitude'] * m * h
            L_out = L_circ * (1 - params['masking_strength'] * M)
            
            # Use last 5 days for amplitude calculation
            start_idx = int(15 * 24 / 0.1)
            L_last = L_out[start_idx:]
            amp = np.max(L_last) - np.min(L_last)
            
            if condition == 'DD':
                amp_dd = amp
            else:
                amp_ld = amp
        
        # Return ratio
        if amp_dd > 0.01:  # Avoid division by very small numbers
            return amp_ld / amp_dd
        else:
            return np.nan
            
    except:
        return np.nan

# =============================================================================
# PARAMETER SPACE EXPLORATION
# =============================================================================

def create_heatmap(param1_name, param1_range, param2_name, param2_range, 
                   fixed_params, resolution=20):
    """
    Create a heatmap for two parameters.
    
    Parameters:
    -----------
    param1_name : str - Name of parameter on x-axis
    param1_range : tuple - (min, max) for parameter 1
    param2_name : str - Name of parameter on y-axis
    param2_range : tuple - (min, max) for parameter 2
    fixed_params : dict - Fixed values for other parameters
    resolution : int - Grid resolution
    
    Returns:
    --------
    param1_vals, param2_vals, ratio_grid
    """
    param1_vals = np.linspace(param1_range[0], param1_range[1], resolution)
    param2_vals = np.linspace(param2_range[0], param2_range[1], resolution)
    
    ratio_grid = np.zeros((resolution, resolution))
    
    print(f"\nCalculating heatmap: {param1_name} vs {param2_name}")
    print(f"Resolution: {resolution}x{resolution} = {resolution**2} points")
    
    for i, p2_val in enumerate(param2_vals):
        if i % 5 == 0:
            print(f"  Progress: {i}/{resolution} rows")
        for j, p1_val in enumerate(param1_vals):
            # Build parameter dict
            params = fixed_params.copy()
            params[param1_name] = p1_val
            params[param2_name] = p2_val
            
            # Calculate ratio
            ratio_grid[i, j] = calculate_ld_dd_ratio(params)
    
    print(f"  Complete!")
    return param1_vals, param2_vals, ratio_grid

# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_parameter_space():
    """
    Create a figure with multiple pairwise parameter space heatmaps.
    """
    
    # Define parameter pairs to explore
    param_pairs = [
        # Most important: opposing mechanisms
        {
            'param1': 'light_sensitivity',
            'param1_range': (0.0, 1.0),
            'param1_label': 'Light Sensitivity\n(Entrainment)',
            'param2': 'masking_strength',
            'param2_range': (0.0, 1.0),
            'param2_label': 'Masking Strength\n(Suppression)',
            'fixed': {'baseline_activity': 0.1, 'tau_h': 2.5}
        },
        # Light sensitivity vs baseline
        {
            'param1': 'light_sensitivity',
            'param1_range': (0.0, 1.0),
            'param1_label': 'Light Sensitivity',
            'param2': 'baseline_activity',
            'param2_range': (0.0, 0.3),
            'param2_label': 'Baseline Activity\n(Opportunism)',
            'fixed': {'masking_strength': 0.4, 'tau_h': 2.5}
        },
        # Masking vs baseline
        {
            'param1': 'masking_strength',
            'param1_range': (0.0, 1.0),
            'param1_label': 'Masking Strength',
            'param2': 'baseline_activity',
            'param2_range': (0.0, 0.3),
            'param2_label': 'Baseline Activity',
            'fixed': {'light_sensitivity': 0.3, 'tau_h': 2.5}
        },
        # Peak width vs light sensitivity
        {
            'param1': 'light_sensitivity',
            'param1_range': (0.0, 1.0),
            'param1_label': 'Light Sensitivity',
            'param2': 'tau_h',
            'param2_range': (1.0, 5.0),
            'param2_label': 'Gate Recovery (tau_h)\n(Peak Width)',
            'fixed': {'masking_strength': 0.4, 'baseline_activity': 0.1}
        },
    ]
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()
    
    # Custom colormap: blue (LD<DD) -> white (LD=DD) -> red (LD>DD)
    colors = ['#2E86AB', '#FFFFFF', '#D62828']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('ratio', colors, N=n_bins)
    
    for idx, pair in enumerate(param_pairs):
        ax = axes[idx]
        
        # Build fixed parameters
        fixed_params = BASE_PARAMS.copy()
        fixed_params.update(pair['fixed'])
        
        # Calculate heatmap
        p1_vals, p2_vals, ratio_grid = create_heatmap(
            pair['param1'], pair['param1_range'],
            pair['param2'], pair['param2_range'],
            fixed_params,
            resolution=20
        )
        
        # Plot heatmap
        im = ax.contourf(p1_vals, p2_vals, ratio_grid,
                        levels=np.linspace(0.5, 1.5, 21),
                        cmap=cmap, extend='both')
        
        # Add contour lines
        contours = ax.contour(p1_vals, p2_vals, ratio_grid,
                             levels=[0.8, 1.0, 1.2],
                             colors='black', linewidths=[1, 2, 1],
                             linestyles=['--', '-', '--'])
        ax.clabel(contours, inline=True, fontsize=9, fmt='%.2f')
        
        # Plot species points
        for species_name, species_data in SPECIES_DATA.items():
            x_val = species_data[pair['param1']]
            y_val = species_data[pair['param2']]
            
            ax.plot(x_val, y_val,
                   marker=species_data['marker'],
                   markersize=15,
                   markerfacecolor=species_data['color'],
                   markeredgecolor='white',
                   markeredgewidth=2,
                   label=species_name,
                   zorder=10)
        
        # Formatting
        ax.set_xlabel(pair['param1_label'], fontsize=11, fontweight='bold')
        ax.set_ylabel(pair['param2_label'], fontsize=11, fontweight='bold')
        ax.set_title(f"Panel {chr(65+idx)}: {pair['param1_label'].split()[0]} vs {pair['param2_label'].split()[0]}",
                    fontsize=12, fontweight='bold', pad=10)
        
        if idx == 0:
            ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
        
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('LD/DD Ratio', fontsize=10, fontweight='bold')
        cbar.ax.axhline(1.0, color='black', linewidth=2, linestyle='-')
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.98])
    fig.suptitle('Parameter Space Exploration: How Parameters Affect LD/DD Ratio',
                 fontsize=14, fontweight='bold', y=0.995)
    
    return fig

# =============================================================================
# RUN
# =============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("PARAMETER SPACE EXPLORATION")
    print("=" * 80)
    print("\nGenerating pairwise parameter space heatmaps...")
    print("This may take a few minutes...")
    
    fig = plot_parameter_space()
    
    filename = 'parameter_space_exploration.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"\n{'=' * 80}")
    print(f"Saved: {filename}")
    print("=" * 80)
    
    plt.show()
