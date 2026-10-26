# Radar gait analysis

Zahraa Sattar · EEE4022S · University of Cape Town

This project investigates walking speed, step recognition, stride duration and stride length from 24GHz continuous-wave radar. Simulated walking is used for verification. Participant recordings are used for measured-data validation.

## Run the measured demonstration

The demonstration was tested in MATLAB R2025a with Signal Processing Toolbox. Download the repository using **Code → Download ZIP**, extract it, and open the extracted repository folder as the MATLAB Current Folder. Keep the folder structure intact.

Run:

```matlab
Run_Measured_Demo
```

The demonstration processes the ten supplied WAV files with Pipeline v14 and the retained Final settings. No recordings from outside this repository are needed. The results are saved in a new `Results/Measured_...` folder:

| Output | Contents |
|---|---|
| `Results.csv` | Reference and estimated values for each recording. |
| `Summary.csv` | MAE, RMSE and percentage error for the four parameters. |
| `Parameter_Comparison.png` | Reference and estimated values across the ten walks. |
| `WalkPlots/` | Leg detections and the torso-speed track for each walk. |
| `Reproduction_Check.csv` | Differences from the retained results, checked within 1e-6. |
| `References_Used.csv`, `Parameters_Used.csv` | Inputs used for that run. |

`Run_Measured_Recordings.m` processes the recordings. `Summarise_Measured_Results.m` creates the tables and comparison figure. These are the smaller demonstration versions of the measured batch and reporting workflow. They run the Final configuration only.

## Recordings and references

`Participant Recordings` contains four walks from P1 and three each from P2 and P3. Within each participant, walks were sorted by reference speed and selected at evenly spaced ranks, including the slowest and fastest eligible walks. Only single-pass, non-limp recordings were considered. Extraction errors were not used for selection. This is a demonstration subset, not the full validation dataset.

The exact WAV versions used in the completed analysis are supplied. Files ending in `_cropped.wav` are the retained crops. `Demo_References.csv` records each file, its sample count, sample rate and SHA-256 checksum. `Expected_Results.csv` contains the Final estimates retained from the completed v18 analysis.

Reference speed is walked distance divided by WAV duration. Reference stride duration is twice the WAV duration divided by the retained step count. Reference stride length is twice the walked distance divided by that count. These are walk-level references, rather than independently timed individual strides.

## Retained parameter settings

- `Extraction Algorithms on Measured Data/Measured_Final_Parameters.csv` records the measured Final configuration and processing settings.
- `Extraction Algorithms on Simulated Data/Simulated_Final_Parameters.csv` records the corresponding extraction configuration and the simulator-verification settings.

Both use the retained simulator-selected base settings, proportional spacing and the recording-edge prominence rule. The torso and leg CPI requests are 0.2s and 0.1s. The achieved windows depend on sample rate and the power-of-two window calculation. The configured leg threshold of -40dB overrides the -45dB option default. The simulator uses 3e8m/s for wave speed. The measured pipeline uses 299792458m/s.

## Check extraction on simulated data

From the repository folder, run:

```matlab
addpath('Extraction Algorithms on Simulated Data')
Verify_Simulated_Extraction(struct('Quick', true))
```

Quick mode uses two walking speeds, four SNR values and two noise draws per cell. It compares proportional spacing with and without the recording-edge rule on the same signals. The recording edges are excluded from scoring, so this test checks whether the edge rule harms the regular simulated gait. It does not establish its benefit at measured recording edges.

For the full speed and noise sweep, run `Verify_Simulated_Extraction` without an argument. Each run saves its tables and report in `Results/Simulated_...`.

## Check the radar simulator

```matlab
addpath('Radar Simulator')
Verify_Radar_Simulator
```

This checks segment Doppler behaviour, torso phase and added noise. Verification figures and tables are saved under `Results`. The verification functions use a temporary simulator copy with workspace resets removed so their test inputs are retained. The simulator equations are unchanged.

## Packaged checks

The measured demonstration reproduced all ten retained Final estimates within 1e-6. The quick simulated comparison passed its three gates. The physical verification completed, with all 36 segment Doppler comparisons passing. The full simulated extraction sweep was not rerun for this package.

## Other project scripts

The four standalone extraction scripts and the CPI/SNR sweeps preserve earlier design tests. Some use shared workspace variables or recording-specific references. They are separate from the ten-recording demonstration. The simulator includes its joint-angle helpers locally, with separate copies available for inspection. `Zahraa_Spectrogram.m` is supplied in each code folder.
