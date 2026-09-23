# Figure 6: Masking and period shortening under ultradian T2 cycles

Representative actogram under LD 1:1 (T2) cycles (A), light versus dark activity and the masking index for each spider (B, C), T2 periods (D), and T2 versus DD periods in the same individuals (E).

## Script

`Fig6.py`

## Inputs

* `Lc_Ag_1204-1223 2025 Monitor1.csv`: the whole recording, used for the raster in A (in this folder)
* `Lc_Ag_1204-1223 2025 Monitor1_T2.csv`: the T2 portion, used for B to D (in this folder)
* `Lc_Ag_1204-1223 2025 Monitor1_LD.csv`: the LD portion of the same run (in this folder; not read by the script, included for completeness)
* `Data/Larinioides/LC_spider_analysis_comprehensive_with_LD_split.csv`: DD periods of the same animals, for E

## How to run

1. Download or clone the whole repository, so that this folder and the `Data` folder sit side by side.
2. Open `Fig6.py` in Spyder. Spyder runs a file from its own folder by default, so outputs land next to the script.
3. Change `DATA_ROOT` (line 32), or edit `RASTER_CSV`, `T2_CSV` and `DD_CSV` (lines 34 to 36) directly so they point to the files listed under Inputs.
4. Run the whole file (F5). Outputs are saved in this folder.

## Outputs

* `Fig6.png` (1000 dpi)
* `Fig6_masking_results.csv`, `Fig6_statistics.csv`: masking index per spider and the Wilcoxon test
* `Fig6_T2_vs_DD_period_results.csv`, `Fig6_T2_vs_DD_statistics.csv`: T2 and DD periods and the paired t-test

## Requirements

Python 3.9 or newer with numpy, pandas, matplotlib, scipy, astropy.
