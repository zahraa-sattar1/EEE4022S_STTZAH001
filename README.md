# Radar gait analysis

Zahraa Sattar · EEE4022S · University of Cape Town

This project investigates walking speed, step recognition, stride duration and stride length from 24GHz continuous-wave radar. Simulated walking is used for verification. Participant recordings are used for measured-data validation.

## Data availability

The consent forms stated that recorded data would not be shared publicly. The participant WAV recordings are therefore not included in the current version of this repository. The repository remains public for the code, synthetic simulation inputs and saved report results.

The full measured extraction cannot be rerun from this repository alone. The saved estimates and matching references allow the report statistics to be recalculated without the recordings.

## Measured-data code

The ten-walk demonstration was tested in MATLAB R2025a with Signal Processing Toolbox before the WAV files were removed. It reproduced the retained Final estimates within 1e-6. It now requires the matching recordings to be supplied locally through an arrangement permitted by the consent. Downloading this repository alone is insufficient to run it.

With authorised local copies in `Participant Recordings`, run:

```matlab
Run_Measured_Demo
```

The folder retains `Demo_References.csv` and `Expected_Results.csv`. These record the ten selected files, their checksums, reference values and saved estimates. Four walks were selected from P1 and three each from P2 and P3 at evenly spaced speed ranks. Extraction errors were not used for selection. This was a demonstration subset, not the full validation dataset.

`Run_Measured_Recordings.m` processes the local recordings. `Summarise_Measured_Results.m` writes the parameter tables and plots. The wrapper compares the estimates against the retained results. Reference speed is distance divided by WAV duration. Reference stride duration is twice the duration divided by the retained step count, and reference stride length is twice the distance divided by that count.

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

## Report result scripts

The original Chapter 4 experiment drivers are also included in `Extraction Algorithms on Simulated Data`:

| Script | Report evidence |
|---|---|
| `Master_SimSweep_v6_TcEstimate.m` | Estimated-cycle CPI/SNR selection and best-case accuracy and repeatability. |
| `Master_SlowSpeedFix_v1.m` | Torso-band comparison and slow-speed appendix figure. |
| `Master_StepWindowFineTune_v1.m` | Step-window comparison and the appendix pass-count figure. |

`Report_Script_Map.md` identifies the retained result folders. These drivers use their own historical inputs from `Report Sources`, rather than the packaged demo extractors. Reduced MATLAB checks completed for all three drivers. Their generated extractor code matches the saved report-run copies. The full sweeps were not repeated for this upload.

The original full measured-data summary is included with saved v18 estimates and matching v17 references. From the repository folder, run:

```matlab
addpath('Extraction Algorithms on Measured Data')
RealData_Summary_v7
```

This regenerates the full-dataset report tables and figures under `Extraction Algorithms on Measured Data/Report Results/Summary_v7_Report_References`. It reads the saved estimates, so the remaining WAV recordings are not required. The main comparison covers 142 single-pass walks. Six double walks are reported separately. Extra diagnostics in the historical summary are not all claims made in the final report. SNR-specific figures need `RealData_SNR_Pth.csv` from a matching diagnostic run. The estimator is included, but that saved input is not.

## SNR estimation on measured data

`EstimateSNR_RealData_v1.m` is in the measured-data folder, with its recording helpers included. It calibrates against simulated signals, then estimates an equivalent SNR and noise-based threshold for each local recording. [SNR_Estimation.md](<Extraction Algorithms on Measured Data/SNR_Estimation.md>) gives the inputs and commands. The WAVs remain excluded under the consent restriction.

## Other project scripts

The four standalone extraction scripts and the CPI/SNR sweeps preserve earlier design tests. Some use shared workspace variables or recording-specific references. They are separate from the ten-recording demonstration. The simulator includes its joint-angle helpers locally, with separate copies available for inspection. `Zahraa_Spectrogram.m` is supplied in each code folder.

## Final-report reproduction check

The latest check uses the complete final PDF dated 27 October 2026. [Report_Reproduction_Map.md](Report_Reproduction_Map.md) lists the commands, matching result versions and remaining gaps. The repository supports simulated tests and recalculation of many full-report statistics from saved estimates. Measured extraction additionally requires authorised local recordings. It does not yet reproduce every report result from raw data.
