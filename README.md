# Auditable Event-Time Verification for G-Code-Driven FFF Thermal Simulation

Research code, G-code inputs, and archived numerical results associated with the NAMRC 55 manuscript **“An Auditable Event-Time Verification Framework for G-Code-Driven Thermal Simulation of Fused Filament Fabrication.”**

## Contents

- `N0/`: 100-layer PLA benchmark, modified Rust solver, B01 toolpath, full-print runner, energy audit, eight archived case summaries and per-run manifests (A0, A2, A6, R2 × 0.10/0.25 K).
- `N1/`: modified Rust solver, two six-layer toolpaths, frozen 36-case design, Windows/Jupyter runners and complete archived case-level JSON, CSV summaries, input manifests, and stdout logs.
- `evidence/N1_representative/`: four selected postprocessing audit outputs, input configurations, raw-history hash manifest (the ~1 GB original raw histories are **not** distributed).
- `evidence/four_case_verifier/`: separate four-case re-audit script and notebook; requires underlying raw histories to actually rerun the audits.
- `scripts/check_archived_results.py`: fast dependency-free comparisons of archived numerical tables against the reported manuscript values; **does not run the thermal model**.
- `docs/ENERGY_LEDGER_CONVENTIONS.md`: explicit reconciliation of the archive and manuscript sign/reference conventions.
- `docs/PROVENANCE.md`: original archive names and technical scope.

## Quick start — check archived results (Python 3.10+)

From the repository root:

```powershell
python scripts/check_archived_results.py
```

This checks that the eight N0 entries and 36 N1 entries are present, verifies the archived values corresponding to the paper's Table 3 (complete-print energy accounting), Table 4 (operating-condition sensitivity), and the four selected N1 re-audit values. This is a consistency check of archived data, not fresh solver reproduction.



## Regenerate the numerical results (Windows 10/11)

Install Python 3.10+, NumPy and Rust/Cargo + MSVC C++ build tools. Refer to each folder's existing README for full prerequisites. Run from inside the respective directory, not the repository root:

```powershell
cd N0
python n0_runner.py --cases A0 A2 A6 R2 --thresholds 0.25
python n0_runner.py --cases A0 A2 A6 R2 --thresholds 0.1 --no-build
cd ..\N1
python run_N1_complete.py
```

N1 also has a `NAMRC_N1_Run_All_Windows.ipynb` notebook. Rerunning full N0 may require significant memory and disk. **Back up the archived result files before rerunning**, as scripts may overwrite outputs. New runs need to be compared to the frozen reference results; differences must be reported rather than silently replaced.

## Paper-result mapping

- Paper Table 3 / Fig. 3: `N0/N0_results/N0_full_print_summary.csv`, with detailed JSON audits for eight cases.
- Paper Table 4 / Figs. 5–7: `N1/N1_results/N1_energy_audit_summary.csv`, `N1/N1_results/paired_full_solver_results.csv`, and 36 case JSON audits.

- Paper Section 4.4: `evidence/N1_representative/N1_representative_verification/` (four verified audit summaries and raw-history hashes).
- Fig. 2 convergence values, and Fig. 1 schematic: this package does not assert complete automated figure regeneration. Underlying solver verification code is provided, but regenerate and review these figures separately before claiming exact plotting reproduction.

## Important scope and limitations

The data records are **archived results** generated in the earlier research workflow. This consolidation has verified their internal agreement with manuscript Tables 3-4 but **has not independently executed all 44 simulations on Windows**. Full cell-temperature CSV histories used in four historical audits are not included due to size; simulations can be rerun to generate new histories, but the original manifest hashes will only match the original bytes. Exact end-to-end reproduction remains unverified here.

Temperature differences are between numerical implementations, not experimental temperature-validation errors. This research models PLA with specified thermal assumptions and does not claim physical FLIR validation. Note that the original N1 runner documentation's `STATUS_NOT_RUN` statement applied to the runner-package creation stage; the separate subsequently generated `N1_results_for_review.zip` archive supplies completed results. That obsolete status marker was deliberately removed during consolidation.

## Attribution & license

This project adapts the published discrete-event thermal solver by S. Bhandari and R. A. Lopez-Anido (2020), *Materials* 13, 4985, https://doi.org/10.3390/ma13214985. The solver source contains its original license (`N0/solver/LICENSE.md`, `N1/solver/LICENSE.md`); retain all original notices and comply with those terms. Check the upstream repository's license and provenance before announcing your own license for non-solver scripts or uploading a public release. No additional blanket license is asserted here. The included `solver_R8_4_git.bundle` is source provenance, not a required execution dependency.

## Citation

Cite the associated NAMRC 55 manuscript once published. Update this README with the manuscript's final author list, DOI, and release DOI after publication. **No DOI has yet been assigned to this ZIP.**
