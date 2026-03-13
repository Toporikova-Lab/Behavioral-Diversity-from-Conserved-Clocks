# QUICK REFERENCE - Larinioides Masking Model

## Two Scripts, Two Figures

### Script 1: Larinioides_masking_experiment.py
**Purpose**: Show masking mechanism separates clock from activity
**Figures**: Raster plots + 3-panel mechanism figure
**Key Control**: `MECHANISM_DAY = 1` (or 20)

### Script 2: light_intensity_period_analysis.py  
**Purpose**: Show light controls period dose-dependently
**Figures**: Period vs light + activity traces (48h)
**Key Control**: `EXAMPLE_LIGHT_INTENSITIES = [0.0, 0.5, 1.0, 1.5, 2.0]`

---

## Main Findings Summary

### Finding 1: Masking Separates Timekeeping from Time-Telling
- Clock (y): Smooth in both LD 1:1 and DD
- Activity (L): Fragmented in LD 1:1, smooth in DD
- Masking (M): Tracks light → suppresses activity
- **Conclusion**: Same clock, different outputs

### Finding 2: Light Intensity Controls Period Parametrically
- Light 0.0: 23.40h
- Light 2.0: 22.65h
- Reduction: 0.75h (45 min)
- Mechanism: k2_effective = 0.18 × (1 + 0.09 × light)
- **Conclusion**: Linear dose-response

---

## What to Change

### To see different days of mechanism:
```python
MECHANISM_DAY = 1    # First day (aligned phases)
MECHANISM_DAY = 20   # Last day (phase drift)
```

### To show different activity traces:
```python
EXAMPLE_LIGHT_INTENSITIES = [0.0, 2.0]              # Just extremes
EXAMPLE_LIGHT_INTENSITIES = [0.0, 0.5, 1.0, 1.5, 2.0]  # Full range
```

---

## Run Commands

```bash
# Generate masking mechanism figures
python Larinioides_masking_experiment.py

# Generate period control figures  
python light_intensity_period_analysis.py
```

---

## Output Files

From Script 1:
- `larinioides_masking_raster.png` - Actograms (LD 1:1 vs DD)
- `larinioides_mechanism_variables.png` - Clock/masking/activity panels

From Script 2:
- `light_intensity_period_control.png` - Period plot + 48h traces

---

## For Your Paper

### Appendix Figure 1 (Mechanism):
"Supplementary Figure shows LD 1:1 fragments activity via masking while clock oscillates smoothly, demonstrating separation of timekeeping from time-telling."

### Appendix Figure 2 (Period Control):  
"Supplementary Figure demonstrates linear period reduction with light intensity (23.4h→22.65h), supporting parametric forcing mechanism."

### Discussion Points:
1. Masking = output pathway, not clock disruption
2. Period controlled by k2 (light increases degradation)
3. Weak coupling (0.09) prevents entrainment to LD 1:1
4. Strong masking (0.9) provides immediate behavioral flexibility
5. Adaptive strategy: stable timing + flexible output

---

## Parameters (DON'T CHANGE)

```python
'light_sensitivity': 0.09    # 9% k2 increase per unit light
'masking_strength': 0.9      # 90% activity suppression
'y_threshold': 0.5           # Gate opens when y > 0.5
```

These define Larinioides strategy:
- Weak light coupling → stable period
- Strong masking → flexible behavior

---

## Remember

**Both scripts use SAME MODEL, different experiments:**
- Script 1: LD 1:1 vs DD (mechanism)
- Script 2: Light intensity sweep (control)

**If you forgot everything else:**
1. Open script
2. Change number at top (MECHANISM_DAY or EXAMPLE_LIGHT_INTENSITIES)
3. Run script
4. Get figure for Appendix
