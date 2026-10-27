# Participant recordings

This folder contains ten example radar WAV recordings: four from P1 and three each from P2 and P3. The files ending in `_cropped.wav` are the retained crops used in the analysis.

`Demo_References.csv` lists each file, its reference values, sample count and checksum. `Expected_Results.csv` contains the saved estimates for comparison.

From the repository folder in MATLAB, run:

```matlab
Run_Measured_Demo
```

The script processes the ten recordings, plots the results and checks the estimates against the saved values. These ten walks are an example, rather than the full measured validation dataset.
