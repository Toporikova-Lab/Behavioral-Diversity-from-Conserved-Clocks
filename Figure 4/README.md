# Figure 4: Dual-pathway model

Model schematic (A), LD/DD amplitude ratio across light sensitivity and masking strength (B), and simulated hourly activity profiles for the three species in DD, LD and LL (C).

## Script

`Fig4.py`

## Inputs

* None. This figure is produced entirely by the model. All parameter values are set at the top of the script and match Tables 1 and 2 of the paper.

## How to run

1. Download or clone the whole repository, so that this folder and the `Data` folder sit side by side.
2. Open `Fig4.py` in Spyder. Spyder runs a file from its own folder by default, so outputs land next to the script.
3. No edits needed.
4. Run the whole file (F5). Outputs are saved in this folder.

## Outputs

* `Fig4.png` (1000 dpi)

## Notes

* Panel A is drawn by the script with matplotlib shapes.
* `HEATMAP_RESOLUTION` (line 64) sets the grid for panel B. The published figure uses 15; larger values give a smoother map but take longer.

## Requirements

Python 3.9 or newer with numpy, scipy, matplotlib.
