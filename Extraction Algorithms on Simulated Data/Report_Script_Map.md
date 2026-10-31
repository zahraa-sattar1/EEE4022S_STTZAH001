# Simulation experiments

| Script | Experiment | Saved run |
|---|---|---|
| `Master_SimSweep_v6_TcEstimate.m` | CPI/SNR selection with estimated cycle duration. | `SimSweep_TcEst_20261009_203704` |
| `Master_SlowSpeedFix_v1.m` | Torso-band comparison at slow speeds. | `SlowSpeedFix_20261009_093327` |
| `Master_StepWindowFineTune_v1.m` | Step-window comparison. | `StepWindowFT_20261010_072324` |

Run a script from this folder or add it to the MATLAB path. Each experiment writes a new output folder.

The drivers use the simulator and extractor versions in `Report Sources`. These are kept separate because the experiments modify specific source expressions in working copies. Full sweeps take longer than Quick mode.

## Optional checks

From the repository folder, run the simulator checks with:

```matlab
addpath('Radar Simulator')
Verify_Radar_Simulator
```

For a short simulated extraction check:

```matlab
addpath('Extraction Algorithms on Simulated Data')
Verify_Simulated_Extraction(struct('Quick',true))
```

These save verification tables and figures under `Results`. Quick mode uses two speeds and fewer noise draws. Run `Verify_Simulated_Extraction` without an argument for the full sweep.
