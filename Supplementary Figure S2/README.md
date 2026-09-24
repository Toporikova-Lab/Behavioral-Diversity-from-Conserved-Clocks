# Supplementary Figure S2 and Supplementary Table S2: individual actogram atlas

This folder reproduces Supplementary Figure S2 (the individual actogram and periodogram atlas) and Supplementary Table S2 (the recordings excluded by the quality filter) from

> Robb DT, Broda J, Petko J, Jones TC, Moore D, Ayoub NA, Toporikova N. *Disentangling Masking and Entrainment in Nocturnal Spiders Through Combined Modeling and Behavioral Experiment.* Journal of Theoretical Biology.

For every spider in every lighting condition (DD, LD 12:12 and LL), the atlas shows a double-plotted actogram next to its Lomb-Scargle periodogram. All actograms use one common activity scale, so you can compare amplitudes directly between individuals, species and conditions.

## What you need

* Python 3.9 or newer (Anaconda with Spyder works well)
* numpy, pandas and matplotlib (included with Anaconda)
* astropy, for the Lomb-Scargle periodogram. If it is missing, install it with `pip install astropy` or `conda install astropy`.

## Data

The script reads the raw TriKinetics activity files in the repository's `Data` folder (one subfolder per species: `Agelenopsis`, `Larinioides`, `Steatoda`). Each activity file has one timestamp column (`datetime` or `Date_Time`), one column per spider (beam crossings per minute) and a `Light` column (1 = lights on). The lighting condition is taken from `DD`, `LD` or `LL` in the file name.

The atlas is built from 16 of the 18 activity files in `Data`. They are listed in `ATLAS_FILES` near the top of the script. The other two files (`Ag 0918-0926 2025 Monitor2_DD.csv` and `StA 04032024_LD.csv`) are used by the main-text figures but were not part of the atlas. Set `ATLAS_FILES = None` to draw every activity file instead.

Set `SHOW_LAB_STATS = True` to print, in blue on each periodogram, the period and p value from the lab analysis files (`*_spider_analysis_comprehensive_with_LD_split.csv`) that Figure 2 uses. The default, `False`, reproduces the atlas exactly as submitted.

## How to run

1. Download this folder and the `Data` folder, or clone the whole repository.
2. Open `supplementary_figure_S2_atlas.py` in Spyder.
3. Near the top, set `DATA_DIR` to the location of the `Data` folder on your computer, for example
   `DATA_DIR = r"C:\Users\me\Downloads\Behavioral-Diversity-from-Conserved-Clocks\Data"`.
   If you cloned the repository and leave this line unchanged, the script finds `../Data` automatically.
4. Run the whole file (F5). You can also run it one cell at a time with Ctrl+Enter.

The full run takes about 5 to 10 minutes. Most of that time is spent writing one PNG per individual. To skip the PNGs, set `SAVE_PNGS = False`.

## Output

Everything is written to an `output` folder next to the script.

| File | Contents |
|---|---|
| `Supplementary_Figure_S2_atlas.pdf` | Title page with legend, a summary page (activity level and peak period of every individual), then one section per species and condition with four individuals per page |
| `Supplementary_Table_S2_excluded_recordings.csv` | Recordings removed by the quality filter, with the number of empty days for each |
| `periodogram_summary.csv` | One row per retained recording: mean activity, Lomb-Scargle peak period and p value |
| `individual_actograms/` | One PNG per retained recording, organized by species and condition |

## Expected result

With the 16 activity files in `ATLAS_FILES`, the script reads 281 individual recordings. It excludes 41 of them with the quality filter and keeps 240 in the atlas:

| Species | DD | LD | LL |
|---|---|---|---|
| *L. cornutus* | 33 | 26 | 16 |
| *A. pennsylvanica* | 23 | 42 | 19 |
| *S. grossa* | 35 | 13 | 33 |

The script prints these counts as it runs. If your numbers differ, check the list of files printed at the start of the run.

## Methods used in the script

* **Binning.** Actograms use 6-minute bins, and full row height equals 6 crossings per bin in every panel. The periodogram uses the unbinned 1-minute counts.
* **Time axis.** LD recordings are plotted in Zeitgeber Time (ZT0 = lights on, read from the `Light` column). DD and LL recordings are plotted in Circadian Time, with CT0 set to 08:00, the lights-on time of the housing LD cycle. You can change this for a whole batch in `CT0_BY_FILE`.
* **Periodogram.** The periodogram is a Lomb-Scargle periodogram (astropy) over 14 to 35 h. The dashed line marks the p = 0.05 false-alarm level (Baluev method). The periods agree closely with the lab pipeline used for Figure 2. The p values from this method are more permissive, so the rhythmicity percentages reported in the paper come from the lab pipeline, not from this script.
* **Quality filter.** A recording is excluded if the spider made zero beam crossings on 3 or more whole days, where a whole day has at least 20 h of recording. This happens when an animal dies, moults or stops moving part way through a run. The same rule is applied to DD, LD and LL. Set `MAX_EMPTY_DAYS = None` to keep every recording.
* **Excluded animals.** Channels named `At...` in the Agelenopsis LL file belong to a different species and are not included.

## Contact

Natalia Toporikova, Washington and Lee University, toporikovan@wlu.edu
