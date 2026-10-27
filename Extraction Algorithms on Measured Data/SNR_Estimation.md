# SNR estimation on measured data

`EstimateSNR_RealData_v1.m` calibrates a peak-to-floor ratio on simulated radar returns, then applies the calibration to local IPM-365 recordings. It also estimates a noise-based threshold and compares it with the fixed thresholds. The recording loader and path resolver are included in the file.

From the repository folder, add the measured-data code folder to the MATLAB path:

```matlab
addpath('Extraction Algorithms on Measured Data')
EstimateSNR_RealData_v1(struct('Quick',true,'RunRealData',false))
```

This checks calibration and validation on simulated signals without needing any recordings. Quick mode uses fewer speeds and noise draws, so its calibration is not the full experiment.

For measured data, set `RecordingFolder` to the folder containing authorised local WAV files and call the function without Quick mode. It searches for cropped files first. The folder can contain `CroppedRecordings`, participant folders such as `P1 Test Walks`, or the WAVs directly. `GroundTruthFile` selects the recording list. Its default is the supplied `Report Results/Report_References.csv`. Every listed walk is attempted; the script does not restrict the list to single passes.

The original script used a v4 reference list. The release defaults to the retained final reference table and includes a `RecordingFolder` option for the repository layout. The calibration settings, seeds, measurement equations and range flags are unchanged.

Each run creates a `SNREstimate_...` folder beside the script. Outputs include:

- `Sim_Calibration.csv` and `Sim_Validation.csv`: known SNRs and measured ratios.
- `Calibration_Summary.csv`: fitted offsets and calibration range.
- `RealData_SNR_Pth.csv`: each recording's equivalent SNR, range flag and threshold estimates.
- `Self_Test.txt` and `SNREstimate_Report.txt`: checks and results.

The SNR estimate is on the simulator's scale. Values outside the calibrated range are flagged and should not be treated as precise SNR measurements. This script is a separate diagnostic; it does not change Pipeline v14's final threshold settings.

The WAV recordings are not included in the public repository under the study's consent restriction. The original project files were left unchanged. MATLAB R2025a and Signal Processing Toolbox were used for the release check.
