# SNR estimation

`EstimateSNR_RealData_v1.m` fits a peak-to-floor calibration to simulated radar signals, then uses it to estimate SNR for the recordings. It also estimates a noise-based threshold. The recording helpers are included in the file.

## Ten-recording example

From the repository folder, run:

```matlab
addpath('Extraction Algorithms on Measured Data')
refs = fullfile(pwd,'Participant Recordings','Demo_References.csv');
cfg = struct('Quick',true,'MaxWalks',10,'GroundTruthFile',refs);
EstimateSNR_RealData_v1(cfg)
```

This uses a short simulated calibration and processes the ten supplied recordings. For the full calibration, set `cfg.Quick = false` before calling the function.

To check the simulated calibration without processing recordings:

```matlab
EstimateSNR_RealData_v1(struct('Quick',true,'RunRealData',false))
```

## Other recordings

Set `RecordingFolder` to the WAV folder and `GroundTruthFile` to the recording list. Cropped files are preferred. The loader checks `CroppedRecordings`, participant folders such as `P1 Test Walks`, and files directly in the WAV folder.

The default recording folder is `Participant Recordings`. The default reference table, `Report Results/Report_References.csv`, lists the full dataset. Use `Demo_References.csv` when only the ten examples are available. Every listed walk is attempted, including double walks and limp recordings if present.

## Outputs

Each run saves a new `SNREstimate_...` folder beside the script.

| File | Contents |
|---|---|
| `Sim_Calibration.csv`, `Sim_Validation.csv` | Known SNRs and measured peak-to-floor ratios. |
| `Calibration_Summary.csv` | Fitted offsets and calibration range. |
| `RealData_SNR_Pth.csv` | Estimated SNR, range flags and thresholds for each walk. |
| `Self_Test.txt`, `SNREstimate_Report.txt` | Checks and results. |

The estimate is on the simulator's SNR scale. Values outside the calibration range are flagged and should not be read as precise measurements. This script does not change Pipeline v14's final threshold.

Tested in MATLAB R2025a with Signal Processing Toolbox.
