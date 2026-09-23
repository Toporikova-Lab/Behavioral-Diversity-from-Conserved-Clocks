# Figure 5: Entrainment and masking during LD to DD transitions

Raster plots of two L. cornutus individuals across an LD to DD transition (A), and model simulations of the same protocol with both light pathways, without entrainment, and without masking (B).

## Script

`Fig5.py`

## Inputs

* `Data/Larinioides/pre-processing/LC LD-DD 01162025/LC 01162025 Monitor1.txt`: raw monitor output for the whole LD to DD run
* `Data/Larinioides/pre-processing/LC LD-DD 01162025/Monitors log 01162025.xlsx`: which spider was in which channel
* `Data/Larinioides/pre-processing/monitor channels info.csv`: channel layout

## How to run

1. Download or clone the whole repository, so that this folder and the `Data` folder sit side by side.
2. Open `Fig5.py` in Spyder. Spyder runs a file from its own folder by default, so outputs land next to the script.
3. Change `DATAFILES_PATH` (line 26) to the location of `Data/Larinioides/pre-processing`.
4. Run the whole file (F5). Outputs are saved in this folder.

## Outputs

* `Fig5.png` (1000 dpi)

## Notes

* Reading the Excel log requires the `openpyxl` package (included with Anaconda).

## Requirements

Python 3.9 or newer with numpy, pandas, matplotlib, scipy, openpyxl.
