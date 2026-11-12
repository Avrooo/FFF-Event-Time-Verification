NAMRC N1 — two toolpaths × 3 deposition speeds × 3 numerical time steps × 2 timing implementations

CURRENT STATUS: TOOLPATHS/CASES/INPUTS GENERATED AND CHECKED. NO NEW FULL THERMAL-SOLVER
RESULTS HAVE BEEN RUN IN THE CREATION ENVIRONMENT. DO NOT STATE THAT 36 SIMULATIONS ARE COMPLETE.

This package contains the ACTUAL R8.4 solver fork, recovered at git commit d812454,
not a replacement/simplified thermal solver. B01 is a 6-layer segment of the previously
executed G-code. The second toolpath is an independent 18-by-10-mm serpentine-raster design
with a perimeter, 6 layers at 0.20 mm, road spacing 0.55 mm; it is a numerical benchmark,
NOT a second physical print or experimental validation.

B01 sample: 190 extruding moves; 6662.8 mm deposited motion.
Serpentine: 126 extruding moves; 2172.0 mm deposited motion.
All feeds changed ONLY for extruding G1 motions. Travel speeds retained.

EXPERIMENT SETUP
- 20 / 40 / 80 mm/s for all extruding moves
- time_step 0.025 / 0.05 / 0.10 s
- legacy_timing 0 versus 1 (D1-D3 are grouped in the switch)
- Both independent toolpath patterns: 36 full thermal simulations / 18 matched pairs
- Same PLA material table, k=0.35 W/(m K), density 1271.1 kg/m3,
  nozzle/deposition model 215 C, bed 60 C, ambient 25 C, h=3 W/(m2 K), no fan,
  M2 bed contact Rc=0.005 m2 K/W, 0.55 x 0.2 mm bead, six layers,
  cooldown 60 s, output storage deadband 0.1 K.
- Fixed physical parameters used to isolate numerical timing; open-room B01 was 23-24 C,
  but this numerical factorial deliberately preserves the predeclared 25 C reference.

HOW TO RUN (WINDOWS ANACONDA PROMPT / POWERSHELL):
1. Install Rust/Cargo from https://rustup.rs/ (use official installation instructions),
   or open a terminal where 'cargo --version' works.
2. Unzip this package into a short folder path, e.g. C:\NAMRC_N1.
3. Open a terminal in that folder.
4. Optional smoke-test: python run_full_solver_N1.py --max-cases 2 --keep-outputs
5. Run all cases: python run_full_solver_N1.py
   To retain raw outputs and run independent energy audits, instead use:
   python run_full_solver_N1.py --keep-outputs
6. For one retained case, run:
   python N1_energy_audit.py runs/<exact_case_directory_name>
   The audit uses the real clock and true physical input laws. It is intentionally
   independent of the numerical update method, so a legacy timing mismatch may appear.
7. Review N1_results/paired_full_solver_results.csv and N1_results/failures.json.

CAUTION: This code has been syntax-checked and dry-tested for input preparation, but
Rust compilation and the full solver were not executable in the creation environment.
Expect to resolve platform-specific build or parsing issues on the first smoke test.

IMPORTANT SCIENTIFIC CONTROLS
* Do not use the original timing switch as a valid predictive model. It's a numerical
  negative control reproducing an identified implementation defect.
* Do not interpret 36 solver cases as 36 independent physical experiments.
* The second toolpath is intentionally different, but both benchmarks share the
  same numerical model and PLA assumptions.
* path_level_clock_diagnostics.csv is a G-code MOTION-level analytical diagnostic,
  not the Rust finite-element event clock and not thermal-solver output. It cannot
  substitute for full solver results, energy audits or temperature differences.
* The D1-D3 timing switch is GROUPED, so the experiment does NOT isolate D1 alone.
* If error is anomalous for one case, inspect full solver log and stability,
  compare geometry/number of active elements, and check stored temperatures.
* Complete an independent energy audit for at least one case at every speed/dt
  combination for the 2 toolpaths if computationally feasible; ideally all cases.
* After completion, archive source commit, raw outputs, run logs, and configuration files.

FILES
solver/                recovered original R8.4 Rust fork (git commit recorded; full git bundle included)
toolpaths/             B01 6-layer segment + independent serpentine G-code
build_n1.py            documented path generator and independent toolpath audits
run_full_solver_N1.py actual Cargo build, simulation matrix, paired result aggregation
N1_energy_audit.py     independent output-based physical energy audit (needs numpy)
experiment_matrix.csv  frozen 36-case design
toolpath_audit.json     independently checked G-code lengths/counts/layers
path_level_clock_diagnostics.csv  analytical PATH-level proxy ONLY
N1_results/            only populated on a machine that runs the solver

RECOMMENDED PAPER PRESENTATION
Call this a two-geometry numerical generalization test only once the Rust output
and numerical energy audits have actually been obtained and checked. The existing
NAMRC draft's R8.4 results are real historical full-solver runs; this package is
additional work that has NOT yet produced new thermal outcomes.
