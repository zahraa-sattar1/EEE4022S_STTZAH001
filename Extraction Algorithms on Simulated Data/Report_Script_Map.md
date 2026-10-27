# Chapter 4 experiment scripts

These are the historical experiment drivers, separate from the measured demonstration and the final adapted extraction check.

| Script | Report result | Retained run |
|---|---|---|
| `Master_SimSweep_v6_TcEstimate.m` | Estimated-cycle CPI/SNR selection, best-case parameter results and noise repeatability. | `SimSweep_TcEst_20261009_203704` |
| `Master_SlowSpeedFix_v1.m` | Torso-band comparison, including the appendix slow-speed figure. | `SlowSpeedFix_20261009_093327` |
| `Master_StepWindowFineTune_v1.m` | Step-window comparison and appendix `FineTune_Passes` figure. | `StepWindowFT_20261010_072324` |

The report links are recorded in the Chapter 4 performance section and appendix. The three files in `Report Sources` are the simulator and extractors needed by these drivers. They are kept separate because the drivers patch specific source expressions before running each experiment. Their numerical code is retained. The standalone demo extractors are not substituted for them.

Run a driver from this folder, or add this folder to the MATLAB path. `Zahraa_Spectrogram.m` is already supplied here. Outputs are written to new timestamped folders next to the drivers. Original project sources and saved report runs are not overwritten.

The full experiments can take substantially longer than the ten-recording demonstration. Reduced packaging checks verify that their source modifications and processing paths execute. They do not replace the archived full-run evidence or reproduce every reported number on their own.
