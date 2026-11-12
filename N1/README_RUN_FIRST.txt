NAMRC N1 — final Windows / Jupyter experiment

1. Extract the NAMRC_N1_FINAL folder into your Downloads folder.
2. Open NAMRC_N1_Run_All_Windows.ipynb and run cells in order, OR run:
      python run_N1_complete.py
   from inside NAMRC_N1_FINAL in Anaconda Prompt.
3. This compiles the actual R8.4 Rust solver, runs the 36-case frozen matrix, retains output histories, conducts energy audits for every case, and saves N1_results_for_review.zip.
4. Runs are restart-safe: previously successful cases whose raw histories are present and saved audits will be reused. If some cases fail, failures_final.json identifies them; they will not be represented as complete.
5. Upload N1_results_for_review.zip back to ChatGPT.

Scientific scope: two SIX-LAYER numerical toolpaths, not two experimental printed parts, speed 20/40/80 mm/s, dt 0.025/0.05/0.10 s, corrected vs grouped D1-D3 legacy timing. 36 real thermal solver simulations. Same frozen no-fan physical parameters at 25C and same output threshold of 0.1K; energy residuals refer to reconstructed physical-clock energy audit.

Storage: raw output histories are retained and may consume considerable disk space. Never edit the frozen matrix mid-run. Do not claim results before actual simulation and audit completion.
