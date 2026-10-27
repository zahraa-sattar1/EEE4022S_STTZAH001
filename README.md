# Radar gait analysis

Zahraa Sattar | EEE4022S | University of Cape Town

This project uses simulated and measured 24GHz CW radar signals to estimate walking speed, step count, stride duration and stride length. The simulator is used to verify the extraction methods. Participant recordings are used to evaluate measured-data performance.

## Folders

| Folder | Contents |
|---|---|
| Radar Simulator | Walking model, radar returns, joint-angle functions and simulator verification. |
| Extraction Algorithms on Simulated Data | Extraction algorithms, CPI/SNR sweeps and simulated comparisons. |
| Extraction Algorithms on Measured Data | Measured-data pipeline, SNR estimation and analysis of saved results. |
| Participant Recordings | Reference values and saved estimates for the ten-walk example. |

## Recordings

Ten example radar WAV recordings are included in `Participant Recordings`, with their reference values and saved estimates. These are a sample of the measured dataset. The full analysis uses saved results for 142 single-pass walks.

## MATLAB

The scripts use MATLAB R2025a and Signal Processing Toolbox. Run the commands below from the repository folder. Keep the folder structure unchanged.

### Simulator

```matlab
addpath('Radar Simulator')
Verify_Radar_Simulator
```

This checks segment Doppler, torso phase and noise power. Results are saved under `Results`.

### Simulated extraction

```matlab
addpath('Extraction Algorithms on Simulated Data')
Verify_Simulated_Extraction(struct('Quick',true))
```

Quick mode uses two speeds and fewer noise draws. Omit the argument for the full sweep. The CPI/SNR and design-comparison scripts are separate experiments. Their inputs are kept in `Report Sources` and `Report Results`.

### Saved measured results

```matlab
addpath('Extraction Algorithms on Measured Data')
RealData_Summary_v7
```

The main results use 142 single-pass walks. Six double walks are analysed separately. `Report_Comparisons` gives the development-stage and centroid comparisons. `Make_Repeatability_MATLAB` plots repeated-walk variation, and `RealData_SplitCheck_v1` performs the split and calibration analysis. These read saved estimates rather than reprocessing the recordings. Some diagnostic plots need additional saved inputs.

### SNR estimation

`EstimateSNR_RealData_v1` calibrates against simulated signals and estimates an equivalent SNR for each local recording. [SNR_Estimation.md](<Extraction Algorithms on Measured Data/SNR_Estimation.md>) gives the inputs and commands.

### Measured-data example

The ten example WAV files are in `Participant Recordings`. Run:

```matlab
Run_Measured_Demo
```

This processes all ten recordings and compares the four parameter estimates with `Expected_Results.csv`. The results, plots and comparison checks are saved in a new folder under `Results`.

## Settings

`Measured_Final_Parameters.csv` and `Simulated_Final_Parameters.csv` contain the final settings for each dataset. The standalone `Finding_...` scripts and individual sweeps are separate from the final measured pipeline. `Zahraa_Spectrogram.m` is included in each code folder.
