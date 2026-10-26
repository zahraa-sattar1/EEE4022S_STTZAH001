# Radar gait analysis

Zahraa Sattar · EEE4022S · University of Cape Town

This project uses simulated 24GHz CW radar returns and measured walking recordings to investigate walking speed, step recognition, stride duration and stride length.

## Folders

| Folder | Contents |
|---|---|
| Radar Simulator | Walking radar simulator and hip, knee, ankle and elbow angle functions. |
| Extraction Algorithms on Simulated Data | Four parameter-extraction scripts and their CPI/SNR sweeps. |
| Extraction Algorithms on Measured Data | Pipeline v14 and four CPI sweeps on individual recordings. |
| Participant Recordings | Reserved for participant WAV recordings. |

The scripts were used in separate tests as the algorithms developed. CPI and SNR settings varied between experiments. Some scripts share workspace variables, while the sweep files include their processing helpers.

MATLAB R2025a and Signal Processing Toolbox were used. The simulator retains its local joint-angle functions. The separate joint-angle files expose the same calculations for inspection.

The standalone simulated-data detector uses known simulation settings. Its leg-swing references are distinct from measured foot-contact times. The measured CPI sweeps retain recording-specific reference values and require the corresponding WAV files. They are earlier design tests, separate from Pipeline v14.
