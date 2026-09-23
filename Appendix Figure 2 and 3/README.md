# Appendix Figures A2 and A3: Population simulation of entrainment

Simulations of 500 model oscillators per species, with free-running periods and rhythmicity in DD calibrated to the measured spiders. They show how the percentage of animals entrained to LD 12:12 depends on light sensitivity β (Figure A2), and which clocks entrain (Figure A3).

**Quick start:** open `RUN_ALL.py` in Spyder and press F5. It runs every stage in order and ends by checking the results against the published numbers (about 35 minutes; set `QUICK_TEST = True` for a 3-minute trial run).

Full instructions, including the stage-by-stage route, the inputs in `spider_data/` and `spider_profiles/`, and what each output file contains, are in [README_HOW_TO_REPRODUCE.md](README_HOW_TO_REPRODUCE.md).

The simulation uses fixed random seeds, so a rerun reproduces the published numbers exactly. `verify_reproduction.py` checks them and writes `verification_report.csv`.

Requirements: Python 3.9 or newer with numpy, scipy, pandas and matplotlib (astropy optional).
