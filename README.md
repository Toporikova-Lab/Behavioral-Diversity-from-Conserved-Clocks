# Behavioral Diversity from Conserved Clocks

Data and code to reproduce every figure in

> Robb DT, Broda J, Petko J, Jones TC, Moore D, Ayoub NA, Toporikova N. *Behavioral Diversity from Conserved Clocks: Combined Modeling and Experimental Approaches Reveal Output Pathway Modulation in Nocturnal Spiders.* Journal of Theoretical Biology.

## Repository layout

| Folder | Reproduces | Uses data from | Main script |
|---|---|---|---|
| `Data/` | raw locomotor activity recordings (see below) | | |
| `Figure 1/` | Figure 1: LD activity profiles and circular phase plots | `Data/` | `Fig1.py` |
| `Figure 2/` | Figure 2: raster plots and free-running periods | `Data/` | `Fig2.py` |
| `Figure 3/` | Figure 3: mean hourly activity profiles | `Data/` | `Fig3.py` |
| `Figure 4/` | Figure 4: dual-pathway model | model only | `Fig4.py` |
| `FIgure 5/` | Figure 5: LD to DD transitions, data and model | `Data/Larinioides/pre-processing/` | `Fig5.py` |
| `Figure 6/` | Figure 6: T2 masking and period shortening | files in the folder, and `Data/` | `Fig6.py` |
| `Appendix/` | Appendix Figure A1: period shortening in the model | model only | `Fig1-Appendix.py` |
| `Appendix/` | Appendix Figures A2 and A3: population simulation of entrainment | files in the folder | `RUN_ALL.py` |
| `Supplementary Figure S2/` | Supplementary Figure S2 and Table S2: individual actogram atlas | `Data/` | `supplementary_figure_S2_atlas.py` |

Each folder has its own README listing the inputs, the one line to edit, and the outputs.

## Data

`Data/` holds one subfolder per species: `Agelenopsis` (*A. pennsylvanica*), `Larinioides` (*L. cornutus*) and `Steatoda` (*S. grossa*). Each contains:

* **Raw activity files** from TriKinetics LAM25 monitors, one CSV per monitor and lighting condition. The condition (`DD`, `LD` or `LL`) is in the file name. Each file has a timestamp column (`datetime` or `Date_Time`), one column per spider with beam crossings per minute, and a `Light` column (1 = lights on).
* **`*_spider_analysis_comprehensive_with_LD_split.csv`**: Lomb-Scargle period, amplitude and p value for every spider and condition, from the lab analysis pipeline. Figures 1, 2, 3 and 6 use these values to classify animals as rhythmic and entrained.
* A drawing of the species (`LC.png`, `AP.png`, `SG.png`), used as an inset in Figure 1.

`Data/Larinioides/pre-processing/` holds the unprocessed monitor output and channel log for the LD to DD experiment in Figure 5.

## Running the code

All scripts are plain Python, written to run in Spyder (Anaconda). Download or clone the whole repository so that the figure folders and `Data/` sit side by side. In each script, point the path near the top to your copy of `Data/` (the folder README gives the line), then run the file with F5.

Requirements: Python 3.9 or newer with numpy, pandas, scipy and matplotlib. Figure 6, Appendix Figure A1 and Supplementary Figure S2 also need astropy (`conda install astropy`). Figure 5 needs openpyxl, which is included with Anaconda.

## Contact

Natalia Toporikova, Department of Biology, Washington and Lee University, toporikovan@wlu.edu
