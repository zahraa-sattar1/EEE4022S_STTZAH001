# Radar gait analysis

Zahraa Sattar · EEE4022S · University of Cape Town

Human walking was modelled to generate synthetic 24GHz continuous-wave radar returns. Walking speed, step events, stride duration and stride length were then extracted from the Doppler signal.

This repository contains the standalone simulation and extraction scripts used during the project. The later measured-data batch pipeline is separate.

## Files

| Script | Purpose |
|---|---|
| `RadarSimulator_v11_torso_legs_arms_head_buffered_noplots.m` | Generates body-segment motion, radar returns and simulation references. |
| `Finding_StepRecognition_v5_noplots.m` | Detects leg-swing peaks and calculates step intervals. |
| `Finding_WalkingSpeed_v4_noplots.m` | Estimates walking speed from the torso Doppler track. |
| `Finding_StrideDuration_v1_noplots.m` | Compares detected stride intervals with the simulated gait cycle. |
| `Finding_StrideLength_v2_noplots_torsovelocity.m` | Integrates the torso-velocity track between detected events. |

## Running the code

MATLAB R2025a and Signal Processing Toolbox were used. Save all five files in one folder and make it the MATLAB Current Folder. Run them in this order:

1. Radar simulator
2. Step recognition
3. Walking-speed extraction
4. Stride duration
5. Stride length

The extraction scripts share the simulator workspace. Keep that workspace between steps. The spectrogram helper is included locally where needed, so no additional helper file is required.

## Scope

The simulator uses regular walking motion and additive Gaussian noise. Step detections correspond to peaks in the leg Doppler track, not independently measured foot-contact times. The standalone detector uses the known simulation cycle and noise settings. These scripts should therefore be interpreted as the simulation workflow rather than the final adapted measured-data method.
