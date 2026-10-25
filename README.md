# Radar gait analysis

Zahraa Sattar · EEE4022S · University of Cape Town

This project uses simulated 24GHz CW radar returns and measured walking recordings to investigate walking speed, step recognition, stride duration and stride length.

## Folders

| Folder | Contents |
|---|---|
| Radar Simulator | Walking model and simulated radar returns. |
| Extraction Algorithms on Simulated Data | Step recognition, walking speed, stride duration and stride length scripts. |
| Extraction Algorithms on Measured Data | Reserved for the measured-data pipeline and selected tests. |

The scripts were used in separate tests as the algorithms developed. CPI and SNR settings varied between experiments. Some scripts share workspace variables.

MATLAB R2025a and Signal Processing Toolbox were used. Spectrogram helpers are included where needed. The measured-data scripts and additional sweeps will be added after review.

The current standalone detector uses known simulation settings. Its leg-swing peaks are simulation references, rather than measured foot-contact times.
