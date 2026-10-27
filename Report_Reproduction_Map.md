# Final-report reproduction check

Checked against the 92-page final-report PDF dated 27 October 2026. This map distinguishes running extraction again from recalculating statistics on saved estimates. The ten supplied recordings are a demonstration subset. They cannot reproduce the full 142-walk validation by themselves.

## Commands

From the repository folder in MATLAB R2025a with Signal Processing Toolbox:

```matlab
addpath('Extraction Algorithms on Measured Data')
RealData_Summary_v7
RealData_SplitCheck_v1
Report_Comparisons
Make_Repeatability_MATLAB
Make_Figure_5_8

addpath('Extraction Algorithms on Simulated Data')
Render_Original_Comparisons
```

These commands read saved results. The measured demonstration runs extraction on the ten supplied WAV files:

```matlab
Run_Measured_Demo
```

## Results traced to supplied inputs

| Report result | Code and inputs | What can be repeated |
|---|---|---|
| Chapter 3 segment, phase and noise verification; Table 3.3 and Fig. A.3 | `Radar Simulator/Verify_Radar_Simulator.m` | Fresh simulator verification. All 36 segment comparisons passed the earlier package check. |
| Tables 3.1–3.2 and A.1 | Simulator and its joint-angle functions | Implemented constants, equations and reference quantities. |
| Chapter 4 CPI/SNR selection and upper-limit results, Table 4.5 | `Master_SimSweep_v6_TcEstimate.m`; saved Stage 1, 2 and 4 CSVs | Historical driver and saved numerical results. Full sweep was not repeated for this release. |
| Table 4.4; Figs. 4.5 and B.2 | `Master_SimVariants_v1.m`, E1 defaults; `E1_overlap.csv`; `Render_Original_Comparisons.m` | Original overlap/drop comparison and plots from the retained 30-draw results. |
| Fig. 4.7 and Appendix B.3 | `Render_Original_Comparisons.m`; Stage 1, 2 and 4 CSVs | Window/SNR and best-case accuracy/repeatability plots from saved results. The later report layout differs. |
| Slow-speed band comparison and Appendix B.4 | `Master_SlowSpeedFix_v1.m`; `Fix_Validation.csv` | Original experiment and saved results. |
| Step-window test, Appendix B.5 | `Master_StepWindowFineTune_v1.m` | Original experiment. Reduced check completed in the earlier package. |
| Table 5.1 and final scores in Table 5.4 | `RealData_Summary_v7.m`; full `Results_AllConfigs.csv` and `Report_References.csv` | Full-cohort statistics from saved estimates, including 140/142 exact counts. |
| Fig. 5.2/5.5 and Table C.2 stage values | `Report_Comparisons.m`; `Baseline_PerWalk.csv` and `Results_AllConfigs.csv` | Scores and count totals for the six stages, starting at 91/142. The baseline CSV came from the original fixed −45 dB run. |
| Table 5.5, strongest bin versus centroid | `Report_Comparisons.m`; `Sim_Trials.csv` and `Measured_Walks.csv` | Aggregate comparison from the **14 October** benchmark used in the report. The 25 October rerun has different values and is deliberately not substituted. |
| Fig. 5.8, P1_T47 and P2_T9 | `Make_Figure_5_8.m`; `Undercount_Comparison.mat` | Redraw stored spectrograms, detections and annotations. This does not reprocess those two recordings. |
| Fig. 5.10 and Table C.4 grouping | `Make_Repeatability_MATLAB.m`; final saved results | Participant/pace repeatability with sessions pooled. Overall pooled SDs: 3.5647, 4.2610 and 4.6150 percentage points. |
| 103/39 split and 500-split calibration analysis | `RealData_SplitCheck_v1.m`; saved estimates/references | Repeats the original seed-2026 analysis, including fitted k = 1.0440 and mean held-out length score 92.98. Development used these walks, so this is a stability check. |
| Table C.3, arm/leg/swing comparison | `Report_Comparisons.m`; `ArmLegAblation_Results.csv` | Aggregates the original **Octave** results used in the report. |
| MATLAB arm/leg mechanism test | `Master_ArmLegAblation_v1.m`, `Report Sources/RadarSimulator_Ablation.m`, retained selection and Pipeline v14 | Runs the corresponding model experiment. A one-speed smoke check completed; equality with the original Octave table has not been established. |

`Report_Comparisons.m` is a small reproduction adapter for the saved tables, rather than a historical extraction algorithm. It does not tune parameters or replace estimates. The original development and plotting scripts are retained where available, with local paths and shorter comments in the release copies.

## Participant-data permissions

The consent information sheet says recorded data will not be shared publicly. Do not expand the public recording dataset without clarification from the supervisor or ethics committee. The existing public recordings and participant-level result inputs also need review. Code and synthetic simulation inputs are separate from this permission question.

## Remaining gaps

This repository is **not yet a complete fresh-data reproduction of every final-report result**.

| Report item | Missing evidence or input |
|---|---|
| All 142 measured walks from raw signals | Only ten WAVs are supplied. Saved estimates allow the full statistics to be recalculated, but not independently extracted again for every walk. |
| Table 4.6 and its radar-height discussion | The exact original height-comparison run and the step-time ICC calculation have not been located or verified. The available simulator can accept changed geometry, but that alone does not establish the reported numbers. |
| Exact fresh simulation of Table C.3 | The report used Octave outputs. The MATLAB driver runs, but its full numerical agreement with that original run is unverified. |
| Table C.1, before/after boundary review | Requires the exact paired historical estimates and reference versions for the ten reviewed crops. Final v18 estimates alone are insufficient. |
| Measured SNR statements and the final SNR-to-performance comparison | The exact diagnostic input and matching run used for these statements still need to be pinned. The summary currently skips the missing SNR input. |
| Measured body-return diagnosis, Fig. 5.4 | Final per-walk estimates do not contain the body spectra. Exact saved plot inputs and matching diagnostic code still need to be packaged. |
| Appendix C.2–C.7 examples | The report uses twelve particular walks, which differ from the ten demonstration walks. Their raw recordings or explicitly labelled saved figure data are needed. |
| Exact report layouts and illustrative figures | Numerical analysis and a plot with the same data do not necessarily recreate the final report formatting. Figs. 3.1–3.3, 3.6, 4.1–4.4, 4.6, 5.1, 5.3, 5.6–5.7, 5.9, A.1–A.2, B.1 and C.1 are not all covered by a tested, exact-layout renderer. |
| Published comparison values in Tables 4.6 and 5.6 | These are external literature results. The project code cannot reproduce the other studies' experiments. Project-side aggregate statistics are available, but the complete Table 5.6 calculation has not been verified end to end. |

These gaps are recorded rather than filled with a different version or newly generated values. Original project files were left unchanged.
