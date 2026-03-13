# T-cycle Masking Analysis - User Guide

## Overview

This script analyzes how spider locomotor activity differs between light and dark phases in T-cycle experiments (e.g., T2 = 1 hour light, 1 hour dark). It calculates a **Masking Index** for each spider and generates publication-ready figures and statistics.

---

## Quick Start

### Step 1: Prepare your data file

Your CSV file should look like this:

```
datetime,Light,LcF8,LcF23,AgM5
2025-01-16 11:00,1,0,5,2
2025-01-16 11:01,0,3,8,1
2025-01-16 11:02,1,2,4,0
```

- First column: `datetime` (date and time)
- Second column: `Light` (1 = light on, 0 = light off)
- Remaining columns: one per spider (column name = spider ID)

### Step 2: Set your file path

Open `Tcycle_masking_analysis.py` in Spyder. Find this line near the top (~line 72):

```python
DATA_FILE = r"C:\Users\...\your_file.csv"
```

Replace the path with your CSV file location. Keep the `r` before the quotes.

**Example:**
```python
DATA_FILE = r"C:\Users\StudentName\Documents\T2_experiment.csv"
```

### Step 3: Run the script

Press **F5** in Spyder (or click the green "Run" button).

---

## Output Files

The script creates three files in the same folder as your input CSV:

| File | Contents |
|------|----------|
| `{filename}_masking_results.csv` | Per-spider data (for detailed reporting) |
| `{filename}_statistics.csv` | Statistical test results (for methods/results sections) |
| `{filename}_masking_figure.png` | Two-panel figure (for publication) |

---

## Understanding the Masking Index

The Masking Index (MI) measures whether a spider is more active in light or dark:

```
MI = (Dark activity - Light activity) / (Dark activity + Light activity)
```

| MI Value | Interpretation |
|----------|----------------|
| MI > 0 | More active in **dark** (light suppresses activity) |
| MI < 0 | More active in **light** (light activates activity) |
| MI = 0 | Equal activity in light and dark |
| MI = +1 | Only active in dark |
| MI = -1 | Only active in light |

---

## Reading the Output Files

### masking_results.csv

| Column | Description |
|--------|-------------|
| `spider_id` | Spider identifier from your CSV column name |
| `mean_light` | Average activity (crossings/min) during light phases |
| `mean_dark` | Average activity (crossings/min) during dark phases |
| `masking_index` | Calculated MI for this spider |

**Use this file to:** Report individual spider values, calculate additional statistics, identify outliers.

### statistics.csv

| Column | Description |
|--------|-------------|
| `test` | Name of statistical test |
| `comparison` | What was compared |
| `statistic` | Test statistic value |
| `p_value` | P-value (significant if < 0.05) |
| `n` | Sample size |
| `mean` | Mean masking index (for t-test row) |
| `std` | Standard deviation (for t-test row) |
| `sem` | Standard error of mean (for t-test row) |

**Use this file to:** Write your results section and figure captions.

**Example results text:**
> "Spiders showed significantly higher activity during dark phases (Wilcoxon signed-rank test, p = 0.003, n = 15). The mean masking index (0.42 ± 0.08 SEM) was significantly greater than zero (one-sample t-test, p < 0.001), indicating light suppression of activity."

---

## Understanding the Figure

The script creates a two-panel figure:

### Panel A (Left): Light vs Dark Activity
- **Gold points/box**: Activity during light phases
- **Dark blue points/box**: Activity during dark phases
- **Gray lines**: Connect the same spider across conditions
- **Y-axis**: Mean crossings per minute

### Panel B (Right): Masking Index Distribution
- **Green points**: Spiders with positive MI (more active in dark)
- **Red points**: Spiders with negative MI (more active in light)
- **Red dashed line**: Zero reference line
- **Formula**: Shows how MI is calculated

---

## Troubleshooting

### "File not found" error
- Check that your file path is correct
- Make sure the `r` is before the opening quote: `r"C:\..."`
- Check for typos in the filename

### Empty or wrong results
- Verify your CSV has a `Light` column (capital L)
- Check that Light values are 0 and 1 (not "on"/"off")
- Make sure datetime column is first

### Figure looks wrong
- Check that your data has both light and dark periods
- Verify activity values are numbers, not text

### Script won't run
- Make sure you have the required packages installed:
  ```
  pip install pandas numpy matplotlib scipy
  ```

---

## Modifying the Figure

Common changes you might need (ask Claude or check the script comments):

| Change | What to modify |
|--------|----------------|
| Figure size | `figsize=(10, 5)` - change numbers |
| Point size | `s=50` - increase/decrease |
| Colors | Hex codes like `'#FFD700'` |
| Y-axis label | `ax1.set_ylabel('...')` |
| Formula position | `ax2.text(1, -0.75, ...)` - change -0.75 |

---

## Contact

If you encounter issues not covered here, contact Dr. Toporikova or bring your CSV file and error message to lab meeting.
