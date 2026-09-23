# =============================================================================
# RUN_ALL.py  --  reproduce Simulation S1 end to end
#
# Press F5 in Spyder to run everything, or step through the cells with Ctrl+Enter.
# Each stage is the same script that produced the published results; this file only
# runs them in order and in the right folders, so nothing is reimplemented here.
#
# Inputs needed (already in this folder):
#   spider_data\      the three *_comprehensive_with_LD_split.csv and the three
#                     *_circular_stats_v2.csv files
#   spider_profiles\  the six *_aligned_activity.csv files
# Everything else is regenerated.
#
# Requires numpy, scipy, pandas, matplotlib. astropy is optional (one cross-check).
# =============================================================================

# %% configuration
import os
import sys
import time
import shutil
import subprocess

try:
    BASE = os.path.dirname(os.path.abspath(__file__))
except NameError:
    BASE = os.getcwd()
os.chdir(BASE)

# QUICK_TEST = True runs 40 individuals per species instead of 500: about 3 minutes
# instead of about 35, and the percentages will differ from the published ones by a
# few points. Use it to check the pipeline runs, then set it back to False.
QUICK_TEST = False
N_LIMIT = 40

# rebuilding the 20,000-clock library takes about 3 minutes; set False to keep the
# dd_library.csv already in the folder
REBUILD_LIBRARY = False

print(f"working in {BASE}")
print(f"python {sys.version.split()[0]}")
for mod in ["numpy", "scipy", "pandas", "matplotlib", "astropy"]:
    try:
        m = __import__(mod)
        print(f"  {mod:<12} {getattr(m, '__version__', '?')}")
    except ImportError:
        print(f"  {mod:<12} NOT INSTALLED" + ("  (optional)" if mod == "astropy" else "  <-- needed"))


# %% stage runner
def run(script, cwd=None, env_extra=None):
    """Run one stage as its own process, streaming its output live (nothing is
    buffered or truncated -- this is what makes stage 5's beta/oscillator
    progress lines visible while it runs), and report how long it took."""
    cwd = cwd or BASE
    env = dict(os.environ)
    if QUICK_TEST:
        env["S1_N_LIMIT"] = str(N_LIMIT)
    env.update(env_extra or {})
    env["PYTHONUNBUFFERED"] = "1"  # so the child's own prints aren't held back in a block buffer
    label = os.path.relpath(os.path.join(cwd, script), BASE)
    print(f"\n{'=' * 78}\n>>> {label}\n{'=' * 78}")
    t0 = time.time()
    proc = subprocess.Popen([sys.executable, script], cwd=cwd, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    for line in proc.stdout:
        print(line, end="")
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"{label} failed with exit code {proc.returncode}")
    print(f"--- {label} finished in {time.time() - t0:.0f} s")


T_START = time.time()

# %% stage 1  --  waveform yardstick from the recorded spiders   (seconds)
run("s1_spider_reference.py")

# %% stage 2  --  the 20,000-clock library in constant darkness   (about 3 min)
if REBUILD_LIBRARY:
    run("s1a_dd_library.py")
else:
    print("skipping s1a, using the dd_library.csv already in the folder")

# %% stage 3  --  fit the three population parameters per species   (about 20 s)
run("s1b_calibrate.py")

# %% stage 4  --  draw the population   (seconds)
run("s1c_population.py")

# %% stage 5  --  the beta sweep, inside track1_main/   (about 8 min)
TRACK = "track1_main"
os.makedirs(os.path.join(BASE, TRACK), exist_ok=True)
shutil.copy2(os.path.join(BASE, "s1d_beta_sweep.py"), os.path.join(BASE, TRACK, "s1d_beta_sweep.py"))
shutil.copy2(os.path.join(BASE, "population.csv"), os.path.join(BASE, TRACK, "population.csv"))
run("s1d_beta_sweep.py", cwd=os.path.join(BASE, TRACK))

# %% stage 6  --  quality control and results, inside track1_main/   (about 2 min)
for script in ["s1e_shape_qc.py", "s1f_results_figures.py"]:
    shutil.copy2(os.path.join(BASE, script), os.path.join(BASE, TRACK, script))
    run(script, cwd=os.path.join(BASE, TRACK))

# %% stage 7  --  the four supplement figures   (seconds)
run("s1g_supplement_figures.py")

# %% stage 7b  --  copy two supplement figures next to RUN_ALL.py under their
# appendix names, for convenience when assembling the manuscript appendix
shutil.copy2(os.path.join(BASE, "supplement_figures", "FigS1_1_entrainment_vs_beta.png"),
             os.path.join(BASE, "Fig2-Appendix.png"))
shutil.copy2(os.path.join(BASE, "supplement_figures", "FigS1_2_which_clocks_entrain.png"),
             os.path.join(BASE, "Fig3-Appendix.png"))
print("copied FigS1_1_entrainment_vs_beta.png -> Fig2-Appendix.png")
print("copied FigS1_2_which_clocks_entrain.png -> Fig3-Appendix.png")

# %% stage 8  --  did it reproduce?
run("verify_reproduction.py")
print(f"\nwhole pipeline: {(time.time() - T_START) / 60:.1f} minutes")
