# Figure 2: Raster plots and free-running periods

Representative raster plots for one individual per species in DD, LD 12:12 and LL (A to C), and the distribution of Lomb-Scargle periods for all individuals in each condition, with the percentage of rhythmic animals (D to F).

## Script

`Fig2.py`

## Inputs

* Nine raw activity files in `Data/<species>/`, one per raster panel. The files and spider IDs are listed near the top of the script (lines 50 to 70).
* `Data/<species>/*_spider_analysis_comprehensive_with_LD_split.csv`: Lomb-Scargle period and p value for every spider (panels D to F)

## How to run

1. Download or clone the whole repository, so that this folder and the `Data` folder sit side by side.
2. Open `Fig2.py` in Spyder. Spyder runs a file from its own folder by default, so outputs land next to the script.
3. Change `BASE_FOLDER` (line 34) to the location of the `Data` folder.
4. Run the whole file (F5). Outputs are saved in this folder.

## Outputs

* `Fig2.png` (1000 dpi)

## Notes

* The three raster rows within a species come from different individuals, because no animal was recorded in all three conditions.

## Requirements

Python 3.9 or newer with numpy, pandas, matplotlib, scipy.
