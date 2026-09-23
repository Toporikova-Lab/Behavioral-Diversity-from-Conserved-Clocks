# Appendix Figure A1: Period shortening in ultradian light cycles (model)

Model actograms in LD 1:1 and DD (A, B), internal model variables over one cycle (C), the period of the LD 1:1 rhythm as a function of light intensity (D), and example activity traces at several intensities (E).

## Script

`Fig1-Appendix.py`

## Inputs

* None. This figure is produced entirely by the model, with the L. cornutus parameters from Table 2.

## How to run

1. Download or clone the whole repository, so that this folder and the `Data` folder sit side by side.
2. Open `Fig1-Appendix.py` in Spyder. Spyder runs a file from its own folder by default, so outputs land next to the script.
3. No edits needed.
4. Run the whole file (F5). Outputs are saved in this folder.

## Outputs

* `Fig1-Appendix.png` (1000 dpi), written to Spyder's working directory

## Notes

* The light-intensity sweep in D runs one simulation per intensity, so this is the slowest of the model figures.

## Requirements

Python 3.9 or newer with numpy, scipy, matplotlib, astropy.
