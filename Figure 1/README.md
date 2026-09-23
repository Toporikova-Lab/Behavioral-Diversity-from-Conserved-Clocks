# Figure 1: LD activity profiles and phase of entrainment

Population activity profiles over the 24 h LD 12:12 cycle for each species (A to C), and circular plots of the activity phase and vector strength of every individual (D to F).

## Script

`Fig1.py`

## Inputs

* `Data/<species>/` raw LD activity files (every file with `LD` in its name)
* `Data/<species>/*_spider_analysis_comprehensive_with_LD_split.csv`: Lomb-Scargle period and p value for every spider, used to classify individuals as rhythmic and entrained
* `Data/<species>/LC.png`, `AP.png`, `SG.png`: spider drawings used as insets

## How to run

1. Download or clone the whole repository, so that this folder and the `Data` folder sit side by side.
2. Open `Fig1.py` in Spyder. Spyder runs a file from its own folder by default, so outputs land next to the script.
3. Change `BASE_FOLDER` (line 29) to the location of the `Data` folder, and `FIGURE_FOLDER` (line 35) to this folder.
4. Run the whole file (F5). Outputs are saved in this folder.

## Outputs

* `Fig1.png` (1000 dpi)
* `species_daily_activity_comparison_LD_v2.csv`: the hourly profiles plotted in A to C
* `Larinioides_circular_stats_v2.csv`, `Agelenopsis_circular_stats_v2.csv`, `Steatoda_circular_stats_v2.csv`: phase, vector strength and classification of every individual (the source of the entrainment percentages in D to F). These are written to Spyder's working directory, so set it to this folder first.

## Notes

* `Fig1-v2.py` is a reduced version that draws only the circular panels (D to F) and writes `Fig1-v2.png`.

## Requirements

Python 3.9 or newer with numpy, pandas, matplotlib.
