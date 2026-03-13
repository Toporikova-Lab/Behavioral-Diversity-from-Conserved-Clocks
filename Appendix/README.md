# Larinioides Masking Model - Appendix Figures

This directory contains two scripts for generating Appendix figures demonstrating the masking mechanism in Larinioides circadian behavior.

## Overview

**Purpose**: Demonstrate that light pulses in LD 1:1 fragment activity through masking while the circadian clock maintains smooth oscillations.

**Model**: Goodwin oscillator with separate masking pathway
- Clock: x, y, z variables (circadian oscillator)
- Gate: h variable (locomotor gate controlled by clock)
- Masking: M variable (direct light suppression)
- Output: L = (m × h) × (1 - α × M)

**Key Finding**: Same circadian clock produces different activity patterns through masking modulation.

---

## Script 1: Larinioides_masking_experiment.py

### What It Does

Creates **two figures** comparing LD 1:1 (ultradian light pulses) vs DD (constant darkness):

1. **Raster plots** (actograms): Shows full time course of activity
2. **Mechanism figure**: Shows internal model variables for one day

### Key Findings

**LD 1:1 condition:**
- Clock (y) oscillates smoothly despite 12 light pulses per day
- Activity (L) is fragmented by masking (M)
- Period ~23.6h (reduced from DD by ~0.3h)

**DD condition:**
- Clock (y) oscillates smoothly
- Activity (L) shows single smooth peak
- Period ~23.9h (natural free-running)

**Conclusion**: Clock maintains timekeeping; masking controls time-telling.

### Main Control Parameters

Located at the **top of script** in EXPERIMENTAL PROTOCOL CONFIGURATION section:

```python
# Duration controls
EQUILIBRATION_DAYS = 20    # Initial LD 12:12 equilibration
LD12_BASELINE_DAYS = 7     # LD 12:12 baseline before test
TEST_CONDITION_DAYS = 20   # Duration of LD 1:1 or DD test

# Mechanism figure
MECHANISM_DAY = 1          # Which day to plot (1-20)
                          # 1 = first day, 20 = last day

# Light timing
FIRST_PULSE_ZT = 0        # When first LD 1:1 pulse starts (0-24)
                          # 0 = at transition, 12 = 12h later
```

### How to Use

1. **Default run** (Day 1 mechanism):
   ```bash
   python Larinioides_masking_experiment.py
   ```

2. **View different day** (e.g., steady state):
   ```python
   MECHANISM_DAY = 20  # Change in script
   ```

3. **Change light phase**:
   ```python
   FIRST_PULSE_ZT = 12  # Start pulses in dark phase
   ```

### Output Files

- `larinioides_masking_raster.png` - Side-by-side actograms (LD 1:1 vs DD)
- `larinioides_mechanism_variables.png` - 3-panel mechanism figure
  - Top: Clock protein (y) with threshold
  - Middle: Activity components (m and h)
  - Bottom: Masking (M) and activity output (L)

### Figure Interpretation

**Mechanism Figure - 3 Panels:**

**Panel 1 (Clock Protein y):**
- Blue line: Clock protein level
- Dashed gray: Threshold (0.5)
- Light blue fill: Region where gate can open (y > threshold)
- Shows: Clock oscillates smoothly in both conditions

**Panel 2 (Activity Components):**
- Sandy brown (m): Gate modulation by clock
- Pink (h): Locomotor gate
- Overlap: Potential for activity
- Shows: Components driven by clock

**Panel 3 (Activity Output):**
- Teal dashed (M): Masking signal
- Orange solid (L): Final activity output
- Formula: L = (m × h) × (1 - α × M)
- Shows: LD 1:1 fragmented by M, DD smooth

**Key Observation**: 
- Left (LD 1:1): M oscillates 12×/day → L fragmented
- Right (DD): M near zero → L smooth
- Clock (y) smooth in BOTH → timekeeping preserved

### For Discussion Section

**Main points to make:**
1. Clock protein (y) maintains smooth oscillation despite ultradian light pulses
2. Masking pathway (M) directly tracks light and suppresses activity
3. Activity fragmentation occurs at output stage, not clock level
4. Period reduction (~0.3h) from parametric forcing via k2 modulation
5. Model demonstrates separation of timekeeping (clock) from time-telling (output)

**Comparison with data:**
- Day 1: Shows immediate response, no adaptation needed
- Day 20: Shows steady-state, period drift accumulated

---

## Script 2: light_intensity_period_analysis.py

### What It Does

Demonstrates that **light intensity parametrically controls circadian period** in LD 1:1 through k2 modulation.

Creates **two-panel figure**:
1. **Left**: Period vs light intensity (shows dose-dependent control)
2. **Right**: Activity traces at selected intensities (shows mechanism in action)

### Key Findings

**Linear period reduction:**
- Light = 0.0 (DD): Period = 23.40 hours
- Light = 2.0: Period = 22.65 hours
- Total reduction: 0.75 hours (45 minutes)
- Rate: 0.375 hours per unit light intensity

**Mechanism:**
```
k2_effective = k2_baseline × (1 + light_sensitivity × light)
k2_effective = 0.18 × (1 + 0.09 × light)
```
Higher k2 → faster protein degradation → shorter period

**Robustness:**
- All intensities (0-2.0) remain rhythmic
- No arrhythmic threshold within tested range
- Oscillations persist even with strong light

### Main Control Parameters

Located at **top of script**:

```python
# Light intensity range
LIGHT_MIN = 0.0
LIGHT_MAX = 2.0
LIGHT_STEP = 0.1

# Which intensities to plot traces for (MODIFY THIS LIST)
EXAMPLE_LIGHT_INTENSITIES = [0.0, 0.5, 1.0, 1.5, 2.0]

# Simulation duration
EQUILIBRATION_DAYS = 10   # Days to reach steady state
ANALYSIS_DAYS = 10        # Days analyzed for period
```

### How to Use

1. **Default run**:
   ```bash
   python light_intensity_period_analysis.py
   ```

2. **Change which traces to plot**:
   ```python
   # Show only extremes
   EXAMPLE_LIGHT_INTENSITIES = [0.0, 2.0]
   
   # Show fine gradation
   EXAMPLE_LIGHT_INTENSITIES = [0.0, 0.4, 0.8, 1.2, 1.6, 2.0]
   
   # Compare mid-range
   EXAMPLE_LIGHT_INTENSITIES = [0.8, 1.0, 1.2]
   ```

3. **Test higher intensities** (find arrhythmic threshold):
   ```python
   LIGHT_MAX = 3.0
   LIGHT_STEP = 0.2
   ```

### Output Files

- `light_intensity_period_control.png` - Two-panel figure

### Figure Interpretation

**Left Panel (Period vs Light):**
- Blue circles: Measured periods from activity (L)
- Connected line: Linear relationship
- Gray dashed: DD reference period
- Annotation box: Key statistics

**Key observation**: Period decreases linearly with light intensity

**Right Panel (48-hour traces):**
- Shows L (activity) over two complete cycles
- Different color = different light intensity
- Blue shading: 1-hour light pulses (LD 1:1 pattern)

**Key observations**:
- Purple (0.0): Two smooth peaks, natural period
- Blue/teal (0.5-1.0): Peaks shift earlier, fragmentation begins
- Green/yellow (1.5-2.0): Strong fragmentation, maximum phase shift
- Phase differences accumulate over 48h → period reduction visible

### For Discussion Section

**Main points to make:**
1. Light intensity parametrically controls period via k2 mechanism
2. Relationship is linear and dose-dependent (no threshold)
3. Period reduction: 9% k2 increase → ~0.36h period reduction
4. Clock remains robust across tested range (no rhythm failure)
5. 48-hour traces show phase accumulation from period differences
6. LD 1:1 tests parametric forcing vs entrainment (too fast to entrain)

**Biological interpretation:**
- Weak light coupling (sensitivity = 0.09) prevents inappropriate entrainment
- Maintains ~24h timing despite ultradian forcing
- Adaptive strategy for variable light environments
- Separation of timekeeping (stable period) from time-telling (flexible output)

---

## Model Parameters

Both scripts use the same **Larinioides parameter set**:

```python
PARAMS = {
    # Circadian oscillator (Goodwin)
    'v1': 0.84, 'v2': 0.84, 'v3': 0.84,      # Synthesis rates
    'k1': 0.18, 'k2': 0.18, 'k3': 0.18,      # Degradation rates
    'n': 12, 'K': 1.0,                        # Hill function
    
    # Light sensitivity
    'light_sensitivity': 0.09,                # k2 modulation (9%)
    
    # Masking pathway
    'masking_strength': 0.9,                  # 90% suppression
    'tau_m': 2.0,                             # Fast response (2h)
    
    # Locomotor gate
    'y_threshold': 0.5,                       # Gate opening threshold
    'h_steepness': 30,                        # Steep activation
    'tau_h': 2.0,                             # Gate time constant
    
    # Activity output
    'sigmoid_steepness': 12,                  # Sharp modulation
    'L_baseline': 0.0,
    'L_amplitude': 2.5
}
```

**Do not modify these** unless exploring different species strategies.

---

## Quick Reference

### For Appendix Figure 1 (Masking Mechanism):
```bash
python Larinioides_masking_experiment.py
```
- Use Day 1 to show immediate response
- Use Day 20 to show steady state
- Mechanism figure shows how clock and output separate

### For Appendix Figure 2 (Period Control):
```bash
python light_intensity_period_analysis.py
```
- Left panel for Discussion: "period is light-controlled"
- Right panel for Discussion: "shows mechanism in action"
- Modify EXAMPLE_LIGHT_INTENSITIES to show specific comparisons

---

## Common Modifications

### Show earlier vs later days:
```python
# In Larinioides_masking_experiment.py
MECHANISM_DAY = 1    # First day - aligned phases
MECHANISM_DAY = 20   # Last day - phase drift visible
```

### Test different light intensities:
```python
# In light_intensity_period_analysis.py
EXAMPLE_LIGHT_INTENSITIES = [0.0, 1.0, 2.0]  # Just three for clarity
```

### Find arrhythmic threshold:
```python
# In light_intensity_period_analysis.py
LIGHT_MAX = 5.0      # Test much higher intensities
LIGHT_STEP = 0.2     # Coarser steps for speed
```

### Change LD 1:1 phase:
```python
# In Larinioides_masking_experiment.py
FIRST_PULSE_ZT = 0   # Light at transition (default)
FIRST_PULSE_ZT = 12  # Light during subjective night
```

---

## Troubleshooting

**Figure looks wrong:**
- Check that you modified the right variable
- Run equilibration is 10-20 days (takes time but necessary)
- Verify light_sensitivity = 0.09 (if 0, no period effect)

**Period not changing:**
- Verify light_sensitivity > 0
- Check LIGHT_MAX > 0
- Ensure LD 1:1 condition (not DD)

**Can't see traces in right panel:**
- Verify light intensities in EXAMPLE_LIGHT_INTENSITIES are within tested range
- Check that values are in list: [0.0, 0.5, 1.0, 1.5, 2.0]

---

## File Checklist

Core scripts (both needed):
- ✓ `Larinioides_masking_experiment.py` - Mechanism demonstration
- ✓ `light_intensity_period_analysis.py` - Period control analysis

Output figures (for Appendix):
- ✓ `larinioides_masking_raster.png` - Actograms comparison
- ✓ `larinioides_mechanism_variables.png` - Internal variables
- ✓ `light_intensity_period_control.png` - Period vs light + traces

Documentation:
- ✓ `README.md` - This file
- ✓ `MECHANISM_DAY_GUIDE.md` - Details on mechanism figure
- ✓ `TWO_PANEL_GUIDE.md` - Details on period analysis figure

---

## Citation Suggestions

**For masking mechanism:**
"Supplementary Figure X shows that the Larinioides model maintains smooth circadian oscillations (clock protein y) while activity output (L) is fragmented by masking (M) in response to ultradian light pulses (LD 1:1). This demonstrates separation of timekeeping from time-telling, with period modestly reduced (~0.3h) through parametric forcing."

**For period control:**
"Supplementary Figure Y demonstrates parametric control of circadian period by light intensity in LD 1:1 conditions. Period decreases linearly from 23.4h (DD) to 22.65h (light=2.0), with activity traces showing progressive phase advancement and masking-induced fragmentation at higher intensities."

---

## Remember in 2 Weeks

**Script 1** = Mechanism figure (clock vs activity separation)
- Change MECHANISM_DAY to see different days
- Shows HOW masking works

**Script 2** = Period control figure (dose-response)  
- Change EXAMPLE_LIGHT_INTENSITIES to pick traces
- Shows THAT light controls period

Both use same model, same parameters. Just run them!
