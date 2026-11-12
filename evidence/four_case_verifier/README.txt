# Final N1 raw-history verification (Windows)

1. Extract this ZIP; copy `verify_N1_four_cases.py` into your existing extracted `NAMRC_N1_FINAL` folder containing `run_N1_complete.py`, `N1_energy_audit.py`, `solver/`, `runs/` and `N1_results/`.
2. Open the included notebook or run in Anaconda Prompt: `python verify_N1_four_cases.py` from that folder.
3. The script reports which of the four raw temperature files survive and re-executes the archived independent-of-solver audit when present. It also checks aggregate energy identities, matches archived results, calculates SHA256 hashes, and packages the evidence.
4. If and only if any raw histories are missing, run `python verify_N1_four_cases.py --rerun-missing` to regenerate **only those missing cases** with the frozen solver. Original archived per-case JSON is backed up and restored.
5. Result: `N1_representative_evidence_for_review.zip`. Upload it for final comparison.

Cases: B01 segment and independent serpentine, both timing modes, 40 mm/s and dt=0.05 s.

Note: Re-executing the same audit is a reproducibility check, not a fully independent second audit algorithm. Raw temperature histories remain on your Windows machine unless `--include-raw` is explicitly used; the ZIP includes hashes.
