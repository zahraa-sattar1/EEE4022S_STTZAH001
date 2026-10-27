# SNR estimation

`EstimateSNR_RealData_v1.m` calibrates the peak-to-floor ratio using simulated radar returns, then estimates SNR for local IPM-365 recordings. It also estimates a noise-based threshold. The recording helpers are included in the script.

To run a short simulated calibration without recordings:

```matlab
addpath('Extraction Algorithms on Measured Data')
EstimateSNR_RealData_v1(struct('Quick',true,'RunRealData',false))
```

For measured data, set `RecordingFolder` to the folder containing the local WAV files. Cropped recordings are preferred. The folder can contain `CroppedRecordings`, participant folders such as `P1 Test Walks`, or the WAV files directly.

`GroundTruthFile` selects the recording list. Its default is `Report Results/Report_References.csv`. Every listed recording is attempted, including double walks and limp recordings if present. Omit Quick mode for the full calibration.

Each run writes a new `SNREstimate_...` folder beside the script:

| Output | Contents |
|---|---|
| `Sim_Calibration.csv`, `Sim_Validation.csv` | Known SNRs and measured peak-to-floor ratios. |
| `Calibration_Summary.csv` | Fitted offsets and calibration range. |
| `RealData_SNR_Pth.csv` | Estimated SNR, range flags and thresholds for each recording. |
| `Self_Test.txt`, `SNREstimate_Report.txt` | Checks and results. |

The SNR estimate is on the simulator's scale. Values outside the calibrated range are flagged and should not be treated as precise measurements. This diagnostic does not change the final threshold in Pipeline v14.

Tested in MATLAB R2025a with Signal Processing Toolbox. The WAV files are not included under the study's consent restriction.
