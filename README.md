# Radar gait analysis

Zahraa Sattar | EEE4022S | University of Cape Town

This project estimates walking speed, step count, stride duration and stride length from 24GHz CW radar signals. Simulated signals are used to verify the methods, and recorded walks are used to test measured-data performance.

## Try the measured-data example

1. Select **Code → Download ZIP** and extract the repository.
2. Open the extracted folder as the Current Folder in MATLAB.
3. Run:

```matlab
Run_Measured_Demo
```

This processes the ten WAV files in `Participant Recordings` using Pipeline v14. It saves the estimates, plots and comparison with `Expected_Results.csv` in a new folder under `Results`.

The example uses four walks from P1 and three each from P2 and P3. The full measured analysis uses 142 single-pass walks, so the ten examples are only a sample of the dataset.

Use MATLAB with Signal Processing Toolbox. The measured example was tested in MATLAB R2025a. Keep the folder structure unchanged.

## Folders

| Folder | Contents |
|---|---|
| Radar Simulator | Walking model, radar signal generation, joint-angle functions and verification. |
| Extraction Algorithms on Simulated Data | Extraction algorithms, CPI/SNR sweeps and design comparisons. |
| Extraction Algorithms on Measured Data | Pipeline v14, SNR estimation and analysis of saved results. |
| Participant Recordings | Ten radar recordings, reference values and expected results. |

The standalone `Finding_...` scripts and individual sweeps were used for separate tests. They are not the sequence called by `Run_Measured_Demo`.

## Saved measured results

```matlab
addpath('Extraction Algorithms on Measured Data')
RealData_Summary_v7
```

This reads the saved estimates for the full dataset and writes the summary tables. It does not need the remaining WAV files. The results cover 142 single-pass walks. Some diagnostic plots need additional saved inputs.

The same folder also contains:

- `Report_Comparisons`: development-stage and centroid comparisons.
- `Make_Repeatability_MATLAB`: repeated-walk error variation.
- `RealData_SplitCheck_v1`: split and calibration analysis.

## SNR estimation

`EstimateSNR_RealData_v1` estimates SNR from each recording using a calibration fitted to simulated signals. See [SNR_Estimation.md](<Extraction Algorithms on Measured Data/SNR_Estimation.md>) for the commands and inputs.

## Parameter settings

`Measured_Final_Parameters.csv` and `Simulated_Final_Parameters.csv` are in their respective algorithm folders. The experiment drivers use the input versions in `Report Sources`. `Zahraa_Spectrogram.m` is included in each code folder.

## Other experiments

The simulator and extraction checks are described in [Simulation experiments](<Extraction Algorithms on Simulated Data/Report_Script_Map.md>). These are separate from the ten-recording example.
