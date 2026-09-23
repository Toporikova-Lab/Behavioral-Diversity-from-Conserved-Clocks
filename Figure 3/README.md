# Figure 3: Mean hourly activity profiles

Mean hourly activity (± SEM) of rhythmic individuals for each species in DD, LD 12:12 and LL. DD and LL profiles are aligned to each animal's own free-running period with the peak at hour 12; LD profiles are aligned to Zeitgeber Time.

## Script

`Fig3.py`

## Inputs

* Raw DD, LD and LL activity files in `Data/<species>/`. The files used for each species and condition are listed in `CONFIGS` near the top of the script.
* `Data/<species>/*_spider_analysis_comprehensive_with_LD_split.csv`: the period of each spider, used for phase alignment and to select rhythmic individuals

## How to run

1. Download or clone the whole repository, so that this folder and the `Data` folder sit side by side.
2. Open `Fig3.py` in Spyder. Spyder runs a file from its own folder by default, so outputs land next to the script.
3. Change `BASE_FOLDER` (line 27) to the location of the `Data` folder.
4. Run the whole file (F5). Outputs are saved in this folder.

## Outputs

* `Fig3.png` (1000 dpi)
* `<species>_dd_aligned_activity.csv`, `_ld_`, `_ll_`: the profiles plotted in each panel
* `peak_detection_all_timepoints.csv`, `peak_detection_significant_only.csv`: hour-by-hour tests for activity peaks

## Notes

* The CSV outputs are written to Spyder's working directory, so set it to this folder first.

## Requirements

Python 3.9 or newer with numpy, pandas, matplotlib, scipy.
