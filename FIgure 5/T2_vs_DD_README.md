# T2 vs DD Period Comparison - User Guide

## Overview

This script compares circadian periods measured under two experimental conditions:

| Condition | Description |
|-----------|-------------|
| **DD** | Constant darkness - measures the free-running circadian period |
| **T2** | Ultradian light cycle (1h light / 1h dark) - tests how rapid light cycling affects period |

---

## Quick Start

1. Set file paths in CONFIGURATION section (lines 27-29)
2. Run script (F5 in Spyder)
3. Output files appear in the **same folder as the script**

---

## Output Files

All files saved in the same directory as the script:

| File | Contents |
|------|----------|
| `{prefix}_period_results.csv` | Per-spider T2 and DD periods |
| `{prefix}_statistics.csv` | Paired t-test, Wilcoxon test results |
| `{prefix}_figure.png` | Two-panel figure |

### period_results.csv columns

| Column | Description |
|--------|-------------|
| `Spider_ID` | Spider identifier |
| `T2_Period_hours` | Period detected in T2 |
| `T2_FAP` | False alarm probability (< 0.05 = significant) |
| `T2_Status` | "significant" or "non_significant" |
| `DD_Period_hours` | Matched DD period |

### statistics.csv columns

| Column | Description |
|--------|-------------|
| `test` | Test name (Paired t-test, Wilcoxon, Direction count) |
| `statistic` | Test statistic |
| `p_value` | P-value |
| `n` | Sample size |
| `mean_DD`, `mean_T2` | Group means (for t-test row) |

---

## Figure Description

**Left panel**: T2 period distribution
- Filled circles = significant (FAP < 0.05)
- Open circles = non-significant
- Box plot shows median and quartiles

**Right panel**: Paired T2 vs DD comparison (significant only)
- Dark green = T2, Light green = DD
- Black horizontal lines = means
- Gray lines connect same spider

---

## Common Modifications

| Change | Location |
|--------|----------|
| Figure size | Line ~175: `figsize=(10, 5)` |
| Point size | Lines ~185, ~197: `s=80` |
| Colors | Lines 40-41: `GREEN_DARK`, `GREEN_LIGHT` |
| Period range | Lines 34-35: `MIN_PERIOD_HOURS`, `MAX_PERIOD_HOURS` |

---

## Script Structure

```
CONFIGURATION (lines 27-41)
    ↓
load_t2_activity(), load_dd_data()     # Load input files
    ↓
analyze_all_spiders()                   # Run periodograms, match to DD
    ↓
run_statistics()                        # Paired tests
    ↓
create_figure()                         # Generate visualization
    ↓
run_analysis()                          # Orchestrates everything, saves files
```
