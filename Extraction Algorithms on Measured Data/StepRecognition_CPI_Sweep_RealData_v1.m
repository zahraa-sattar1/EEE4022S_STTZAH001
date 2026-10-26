%% Step recognition CPI sweep
% Zahraa Sattar | EEE4022S
%
% Compare CPI settings on one measured walking recording.
% Report step-count error against the recording reference.

groundTruthStepCount = 8;   % Actual number of steps taken in this recording

CPIValues = [0.025, 0.1, 0.125, 0.15, 0.2, 0.30, 0.4];   % CPI values
stepCountError = zeros(size(CPIValues));    % one result per CPI value
detectedCounts = zeros(size(CPIValues));    % raw detected counts, useful to look back at
CPI_actual_sweep_s = zeros(size(CPIValues));

%% Load the recording and get Tc once
% LoadRealRecording
wavFile_CW_All = {'HumanWalking.wav'; 'HumanRunning.wav'; 'walking_test.wav'; 'walking_test2.wav'; ...
    'MeTesting.wav'; ...  % 5th recording
    'Subject1_TestWalk.wav'; 'Subject1_TestWalk2.wav'; 'Subject1_TestWalk3.wav'; ...
    'P1_T1.wav'; 'P1_T2.wav'; 'P1_T3.wav'; 'P1_T4.wav'; 'P1_T5.wav'; 'P1_T6.wav'; 'P1_T7.wav'; 'P1_T8.wav'; 'P1_T9.wav'; 'P1_T10.wav'; 'P1_T11.wav';   % same list as the cantenna script
    'P2_T18.wav'; 'P2_T19.wav'; 'P2_T20.wav'; 'P2_T21.wav'; 'P2_T22.wav'; 'P2_T23.wav'; 'P2_T24.wav'; 'P2_T25.wav'; 'P2_T26.wav'; 'P2_T27.wav'; 'P1_T11.wav';   % same list as the cantenna script
    'P3_T2.wav'; 'P1_t6.wav'};

RecordingNo2Process = 32;   % Pick which recording to load

wavFile_CW = wavFile_CW_All{RecordingNo2Process};

c = 299792458;
fc = 24e9;    % radar centre frequency, matches the cantenna hardware

[Y, fs] = audioread(wavFile_CW, 'native');   % read the raw WAV samples at native (integer) precision
RdrRxSignal = -double(Y(:,1));
RdrRxSignal = RdrRxSignal(:).';

fprintf('Loaded %s: %d samples at fs = %d Hz (%.4f s long)\n', wavFile_CW, length(RdrRxSignal), fs, length(RdrRxSignal)/fs);

%% Spectrogram of the raw loaded recording
CPI_preview_s = 0.05;   
N_preview = 2^(nextpow2(round(CPI_preview_s*fs)));
win_preview = hamming(N_preview);
noverlap_preview = floor(0.8*N_preview);
nfft_preview = N_preview;

[S_preview, freq_preview, time_preview] = Zahraa_Spectrogram(RdrRxSignal, win_preview, noverlap_preview, nfft_preview, fs);

speed_preview_mps = freq_preview * (c/fc) / 2;             % converts Doppler frequency to radial speed for a more intuitive y-axis
speedIdx_preview = find(speed_preview_mps >= 0 & speed_preview_mps <= 5);   % only show 0-5 m/s, matches maxSpeed_m_s from cantenna_dop_v3_Zahraa
S_preview_norm = abs(S_preview(speedIdx_preview,:)) / max(max(abs(S_preview(speedIdx_preview,:))));   % normalise so 0dB is the loudest bin
figure('Name','Loaded recording - raw spectrogram preview');
imagesc(time_preview, speed_preview_mps(speedIdx_preview), 20*log10(S_preview_norm + eps), [-60 0]);
axis xy;
xlabel('Time (s)');
ylabel('Speed (m/s)');
title(sprintf('%s - raw spectrogram preview', wavFile_CW));
colorbar; colormap jet;
% Finding_WalkingSpeed_v3_noplots_realdata
%% Finding_walkingSpeed_v1_(NO PLOTS): STEPS 1-8:

%% 1. Initial passband (Step 1's K-scaling, using the real fc)
fc_paper = 5.8e9;              % Yardibi et al. (2011)'s radar carrier, Hz

if ~exist('CPI_s','var')
    CPI_s = 0.2;                % Introduce CPI input variable
end

K = fc / fc_paper;              % fc comes from the simulator workspace
f_low0  = 5 * K;    % multiply low and high frequencies by K scaler
f_high_min = 50*K;
maxExpectedSpeed_mps = 4.0;             % covers brisk walking/light jogging
f_high0 = maxExpectedSpeed_mps * 2 * fc / c;
fprintf('K = %.6f, initial band = [%.4f, %.4f] Hz\n', K, f_low0, f_high0);

%% 2. Steps 1-6
numIterations = 3;
resultsIter = d2_RunGaitBandRefinement(RdrRxSignal, fs, f_low0, f_high_min, f_high0, numIterations, CPI_s); % f_high_min is also an input
for iter = 1:numIterations
    fprintf('Iteration %d: band = [%.2f, %.2f] Hz, centre = %.4f Hz\n', ...
        iter, resultsIter(iter).freq_low, resultsIter(iter).freq_high, resultsIter(iter).centreFreq_Hz);
end
final = resultsIter(end);   % most refined pass

%% 3. Step 7: suppress no-motion frames
finalValidFrequencyIndex = find( abs(final.freq_axis) >= final.freq_low &  abs(final.freq_axis) <= final.freq_high);
energyPerFrame = sum(final.Ssmooth(finalValidFrequencyIndex, :), 1);
energyBaseline = median(energyPerFrame,'omitnan');
energySpread = median(abs(energyPerFrame-energyBaseline),'omitnan');
energyThreshold = max(0,energyBaseline-2*energySpread);
isNoMotion = energyPerFrame < energyThreshold;                          % identify frames with no motion (ie below threshold)
torsoFreq_suppressed = final.torsoFreq_refined_Hz;
torsoFreq_suppressed(isNoMotion) = 0;

%% 4. Step 8: compute walking speed
c = 3e8;   % speed of light, m/s
velocityTrack_mps = torsoFreq_suppressed * c / (2 * fc);
velocityTrack_mps(isNoMotion) = NaN;                              % instantaneous velocity
movingFrames = torsoFreq_suppressed ~= 0;
meanSpeed_mps = mean(velocityTrack_mps(movingFrames), 'omitnan');   % averages speed that are not equal to 0
fprintf('\nExtracted mean walking speed: %.4f m/s\n', meanSpeed_mps);

%% Ground-truth check
if exist('v0','var')
    percentError = 100 * (meanSpeed_mps - v0) / v0;
    fprintf('Ground-truth mean speed v0        : %.3f m/s\n', v0);
    fprintf('Extracted mean speed              : %.3f m/s\n', meanSpeed_mps);
    fprintf('Error vs. ground truth            : %+.1f %%\n', percentError);
end
if all(torsoFreq_suppressed == 0)
    warning('Every frame was suppressed as no-motion, energyThreshold is set too high.');
end

%% FUNCTIONS

%% Step 1: Ideal Bandpass Filter Function
minPeakDistance_s = 0.3;
minPeakProminence_Hz = 0.25*(max(final.torsoFreq_refined_Hz) - min(final.torsoFreq_refined_Hz));
[~, peakIdx] = findpeaks(final.torsoFreq_refined_Hz, ...
    'MinPeakDistance', round(minPeakDistance_s / (final.time_axis(2)-final.time_axis(1))), ...
    'MinPeakProminence', minPeakProminence_Hz);
peakTimes_s = final.time_axis(peakIdx);
stepIntervals_s = diff(peakTimes_s);
Tc = 2 * mean(stepIntervals_s);

%% Sweep CPI, re-running step recognition on the same recording each time
set(0, 'DefaultFigureVisible', 'off');
for cpiIdx = 1:length(CPIValues)
    CPI_s = CPIValues(cpiIdx);
    % Finding_StepRecognition_v5_noplots_realData
    %% Step RECOGNITION FOR REAL DATA  - STEPS 1-8:

    %% Step 0 - define input parameters
    freq_high_Yardibi = 350*24/5.8;
    c = 299792458;
    fc = 24e9;
    lambda = c/fc;

    %% Step 1: average filter on raw radar data
    marginFactor = 1.5;
    avgSamples = floor(fs / (marginFactor * freq_high_Yardibi));
    xAvg = d3_AverageFilterRadarData(RdrRxSignal, avgSamples);
    sampleWindow = 1:500;   % zoom into a short segment to actually see the smoothing
    fprintf('Std of raw signal (this segment): %.4f\n', std(real(RdrRxSignal(sampleWindow))));
    fprintf('Std of averaged signal (this segment): %.4f\n', std(real(xAvg(sampleWindow))));

    %% Step 2: compute spectrogram
    if ~exist('CPI_s','var')
        CPI_s = 0.2;   % Window-length parameter, matching the walking-speed algorithm
    end

    [S_step, freq_axis_step_full, time_axis_step, fsProc_step] = d3_ComputeGaitSpectrogram(xAvg, fs, CPI_s); % change freq_axis_step to freq_axis_step_full

    CPI_actual_s = length(freq_axis_step_full) / fsProc_step;     % CPI_actual_s = N*(1/fs)
    fprintf('Requested CPI %.4f:\n', CPI_s);
    fprintf('Actual CPI %.4f:\n', CPI_actual_s);

    %% Step 3: background noise subtraction
    nonzeroIndex = find( (freq_axis_step_full > 0) & (freq_axis_step_full <= freq_high_Yardibi));
    freq_axis_step = freq_axis_step_full(nonzeroIndex);
    S_step_leg = S_step(nonzeroIndex,:);
    P_step2 = abs(S_step_leg).^2;
    Pm = d3_SubtractBackgroundNoise(S_step_leg);
    fprintf('Check: Pm never negative -> %d\n', all(Pm(:) >= 0));
    fprintf('Check: Pm never exceeds original power -> %d\n', all(Pm(:) <= P_step2(:) + eps));

    %% Step 4: median filter on spectral data
    frame_period_s = time_axis_step(2) - time_axis_step(1);   % actual spacing between spectrogram frames
    medianWindow_s = Tc/8;
    medianFrames = max(1, round(medianWindow_s / frame_period_s));
    Pmedian = d3_MedianFilterSpectrogram(Pm, medianFrames);
    noiseBinIdx = 1;
    fprintf('Std at edge frequency bin BEFORE median filter: %.4f\n', std(Pm(noiseBinIdx,:)));
    fprintf('Std at edge frequency bin AFTER  median filter: %.4f\n', std(Pmedian(noiseBinIdx,:)));
    fprintf('Pmedian range (dB): [%.1f, %.1f]\n', ...
        10*log10(min(Pmedian(Pmedian>0)) + eps), 10*log10(max(Pmedian(:)) + eps));

    figure('Name','Step 4 verification - median filter');
    subplot(2,1,1);
    imagesc(time_axis_step, freq_axis_step, 10*log10(Pm + eps));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Before median filter (Step 3 output)');
    colorbar; colormap jet;
    subplot(2,1,2);
    imagesc(time_axis_step, freq_axis_step, 10*log10(Pmedian + eps));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('After median filter');
    colorbar; colormap jet;

    %% Step 5: Select Leg Frequency

    if exist('PthOverride','var')
        Pth = PthOverride;
    elseif exist('SNR_dB','var')
        Pth = 10 - SNR_dB;
        Pth = max(Pth, -40);
    else
        Pth = -30;
    end

    [legFreq_Hz, isDetected] = d3_SelectLegFrequency(Pmedian, freq_axis_step, Pth);
    fprintf('Frames with a leg-frequency detection: %d / %d\n', sum(isDetected), numel(isDetected));

    %% Step 6
    [legFreq_refined, isOutlier] = d3_RefineLegFrequencyOutliers(legFreq_Hz, isDetected);
    fprintf('Frames flagged as outliers: %d\n', sum(isOutlier));

    figure('Name','Step 6 verification - remove outliers');
    plot(time_axis_step, legFreq_Hz, 'b.');
    hold on;
    plot(time_axis_step, legFreq_refined, 'r-', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Leg frequency (Hz)');
    legend('Raw Step 5 detections', 'Step 6 refined (outliers removed)');
    title('Steps 5-6: leg frequency selection and outlier rejection');
    grid on;

    %% Step 7
    smoothWindow_s = Tc/8;
    avgFrames7 = max(1, round(smoothWindow_s / frame_period_s));
    F_smooth = d3_SmoothLegFrequency(legFreq_refined, avgFrames7);

    %% Step 8
    frequencyNormalization_Hz = max(F_smooth,[],'omitnan');
    normalizedLegFrequency =  F_smooth/frequencyNormalization_Hz;   % normalise leg frequency vector so maximum value is 1
    DopplerRes_Hz = 1/CPI_s;
    cPeak =  DopplerRes_Hz/freq_high_Yardibi;
    stepIdx8 = find( (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(1:end-2)) & ...
                      (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(3:end) + cPeak) ) + 1;
    fprintf('Step 8 raw peaks: %d\n', numel(stepIdx8));

    %% Step 9
    minInterval_s = 0.35;   % Wang et al., 2014
    stepIdx9 = d3_RemoveCloseSteps(stepIdx8, time_axis_step, F_smooth, minInterval_s);
    detectedStepTimes_s = time_axis_step(stepIdx9);
    fprintf('Step 9 final step count: %d\n', numel(stepIdx9));

    %% Step 10: Calculate step and stride parameters
    stepSpeeds_m_s =  F_smooth(stepIdx9)*lambda/2;
    stepDurations_s = diff(detectedStepTimes_s);
    if numel(detectedStepTimes_s) >= 3
        strideDurations_s = detectedStepTimes_s(3:end) - detectedStepTimes_s(1:end-2);
    else
        strideDurations_s = [];
    end

    %% Plot motion-only leg spectrogram and detected step events
    figure('Name','Wang Step Recognition from Leg Micro-Doppler');
    legPlot = Pmedian;
    legNorm = max(legPlot(:));
    legPlot = legPlot/legNorm;
    speed_axis_m_s = freq_axis_step*lambda/2;
    imagesc(time_axis_step, speed_axis_m_s,  10*log10(legPlot + eps), [-70 0]);
    hold on;
    if ~isempty(stepIdx9)
        detectedStepSpeeds_m_s =  F_smooth(stepIdx9)*lambda/2;
        plot(detectedStepTimes_s, detectedStepSpeeds_m_s, 'wo', 'MarkerFaceColor','m', 'MarkerSize', 10);
    end
    xlabel('Time (s)');
    ylabel('Radial speed (m/s)');
    title('Wang Algorithm III: Detected Leg-Swing / Step Instances');
    colorbar;
    colormap('jet');
    axis xy;
    grid on;

    %% FUNCTIONS

    %% Step 1: Average Filter
    detectedCounts(cpiIdx) = numel(detectedStepTimes_s);   % how many steps got detected at this CPI
    stepCountError(cpiIdx) = 100 * abs(detectedCounts(cpiIdx) - groundTruthStepCount) / groundTruthStepCount;   % error vs the known true count
    N_actual = 2^(nextpow2(round(CPI_s*fs)));
    CPI_actual_sweep_s(cpiIdx) = N_actual / fs;
    fprintf('CPI = %.3f s -> detected %d steps (true = %d), error = %.1f%%\n', CPI_s, detectedCounts(cpiIdx), groundTruthStepCount, stepCountError(cpiIdx));
end
set(0, 'DefaultFigureVisible', 'on');
close all;

%% Results
figure('Name','CPI vs step count error (real data)');
plot(CPI_actual_sweep_s, stepCountError, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor','b');
xlabel('Actual CPI (s)');
ylabel('Step count error (%)');
title('Step recognition step-count error vs. CPI (real data)');
grid on;

[bestError, bestIdx] = min(stepCountError);
bestCPI = CPI_actual_sweep_s(bestIdx);
hold on;
plot(bestCPI, bestError, 'rp', 'MarkerSize', 14, 'MarkerFaceColor','r');
fprintf('Best CPI = %.3f s, with step count error = %.1f%%\n', bestCPI, bestError);

fprintf('\n--- Results: CPI sweep on real data ---\n');
for i = 1:length(CPIValues)
    fprintf('CPI_s requested = %.3f s (actual = %.3f s) -> detected %d steps, error = %.1f%%\n', ...
        CPIValues(i), CPI_actual_sweep_s(i), detectedCounts(i), stepCountError(i));
end

function [y, f, mask] = d2_IdealBandpassFilter(x, fs, f_low, f_high)
    isRowInput = isrow(x);
    x = x(:).';
    N = length(x);
    X  = fft(x);
    Xs = fftshift(X);
    f = (-floor(N/2):ceil(N/2)-1) * (fs / N);
    mask = (abs(f) >= f_low) & (abs(f) <= f_high);
    Xf = Xs .* mask;
    Xu = ifftshift(Xf);
    y  = ifft(Xu);
    if ~isRowInput
        y = y.';        % restores the caller's original orientation
    end
end

%% Step 2: Compute Spectogram
function [S, freq_axis, time_axis, fsProc] = d2_ComputeGaitSpectrogram(x, fs, CPI_s, freq_high)  % included CPI_s
Seifert_approach = 0;  % Choose Seifert's approach of computing spectrogram
if (Seifert_approach == 1)
    processingFs = 2560;   % Hz
    if fs > processingFs
        [p, q]  = rat(processingFs / fs, 1e-6);
        xProc   = resample(x, p, q);
        fsProc  = fs * p / q;
    else
        xProc  = x;
        fsProc = fs;
    end
    N = 2^(nextpow2(round(CPI_s*fsProc)));
    win      = hamming(N);
    noverlap = N-1;
    nfft     = N*8;
else
        xProc  = x;
        fsProc = fs;
        N = 2^(nextpow2(round(CPI_s*fsProc)));
        win      = hamming(N);
        noverlap = N*3/4;                           % Wang uses 75% overlap
        nfft     = N;
end
    [S_temp, freq_axis_temp, time_axis] = d2_Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);
    validFrequencyIndex = find( (abs(freq_axis_temp) <= freq_high*3) );
    S = S_temp(validFrequencyIndex,:);                  % extract subset of spectogram
    freq_axis = freq_axis_temp(validFrequencyIndex);
end

%% Step 3: Reduce spectral noise with average filter
function Ssmooth = d2_SmoothSpectrogramNoise(S, avgFrames)
      Ssmooth = movmean(abs(S), avgFrames, 1);          % Apply moving average across Doppler bins (i.e dimension 1)
end

%% Step 4: torso Doppler frequency = frequency of peak
function [torsoFreq_Hz, peakVal] = d2_TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high)
    validFrequencyIndex = find( (abs(freq_axis) >= freq_low) & (abs(freq_axis) <= freq_high) );
    Subset_Ssmooth = Ssmooth(validFrequencyIndex,:);
    [peakVal, peakIdx] = max(Subset_Ssmooth, [], 1);
    torsoFreq_Hz = freq_axis(validFrequencyIndex(peakIdx));
end

%% Step 5: Outlier rejection
function [torsoFreq_refined, isOutlier] = d2_RefineOutliers(torsoFreq_Hz)
    mu    = mean(torsoFreq_Hz, 'omitnan');
    sigma = std(torsoFreq_Hz, 'omitnan');
    isOutlier = abs(torsoFreq_Hz - mu) > 3*sigma;
    torsoFreq_refined = torsoFreq_Hz;
    torsoFreq_refined(isOutlier) = NaN;
    torsoFreq_refined = fillmissing(torsoFreq_refined,'linear','EndValues','nearest'); % Fill in NaNs with 'linear' interpolated values
end

%% Step 6: Iterative band refinement (repeats Steps 1-5)
function results = d2_RunGaitBandRefinement(x, fs, freq_low0, f_high_min, freq_high0, numIterations, CPI_s)
    freq_low  = freq_low0;
    freq_high = freq_high0;
    results = struct('freq_low', {}, 'freq_high', {}, 'torsoFreq_Hz', {}, ...
        'torsoFreq_refined_Hz', {}, 'centreFreq_Hz', {});
    for iter = 1:numIterations
        xFilt = d2_IdealBandpassFilter(x, fs, freq_low, freq_high);
        [S, freq_axis, time_axis, fsProc] = d2_ComputeGaitSpectrogram(xFilt, fs, CPI_s, freq_high);
        Power_S = abs(S).^2;                                    % Compute Power of Spectrogram (linear units)
        Ssmooth = d2_SmoothSpectrogramNoise(Power_S, 21);          % Smooth the Power of Spetrogram (linear units)
        torsoFreq_Hz = d2_TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high);
        torsoFreq_refined_Hz = d2_RefineOutliers(torsoFreq_Hz);
        centreFreq_Hz = mean(torsoFreq_refined_Hz, 'omitnan');
        results(iter).freq_low             = freq_low;
        results(iter).freq_high            = freq_high;
        results(iter).torsoFreq_Hz         = torsoFreq_Hz;
        results(iter).torsoFreq_refined_Hz = torsoFreq_refined_Hz;
        results(iter).centreFreq_Hz        = centreFreq_Hz;
        results(iter).Ssmooth              = Ssmooth;
        results(iter).freq_axis            = freq_axis;
        results(iter).time_axis            = time_axis;
        results(iter).fsProc               = fsProc;
        if iter < numIterations
            freq_low = max(0.5*centreFreq_Hz,freq_low0);       % 0.5*MeanGaitFrequency (see Yardibi paper)
            freq_high = max(2.0*centreFreq_Hz,f_high_min);  % 2*MeanGaitFrequency (see Yardibi paper)
        end
    end
end

%% Step 7: Suppress no-motion frames
function torsoFreq_suppressed = d2_SuppressNoMotion(torsoFreq_Hz, Ssmooth, energyThreshold)
    energyPerFrame = sum(Ssmooth, 1);      % total energy across all bins, per frame
    isNoMotion = energyPerFrame < energyThreshold;
    torsoFreq_suppressed = torsoFreq_Hz;
    torsoFreq_suppressed(isNoMotion) = 0;   % wherever isNoMotion is true, that position is set to 0
end

%% Step 8: Compute walking speed
function [velocityTrack_mps, meanSpeed_mps] = d2_ComputeWalkingSpeed(torsoFreq_Hz, ft_Hz, c)
    velocityTrack_mps = torsoFreq_Hz * c / (2 * ft_Hz);
    movingFrames = torsoFreq_Hz ~= 0;
    meanSpeed_mps = mean(velocityTrack_mps(movingFrames), 'omitnan');   % averages speed that are not equal to 0
end

%% Zahraa_Spectrogram   Manually computed spectrogram (no built-in spectrogram).
function [S, freq_axis, time_axis] = d2_Zahraa_Spectrogram(x, win, noverlap, nfft, fs)
    x   = x(:).';
    win = win(:).';
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts
    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);
    frame_starts = zeros(1, M);      % keep track of where each frame starts (for time axis later)
    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % grab the chunk and window it
        S(:, m) = fft(segment, nfft).';
    end
    time_axis = (frame_starts - 1) / fs;
    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end

function xAvg = d3_AverageFilterRadarData(x, avgSamples)
    isRowInput = isrow(x);
    x = x(:).';
    xAvg = movmean(x, avgSamples);
    if ~isRowInput
        xAvg = xAvg.';
    end
end

%% Step 2: Compute Spectogram
function [S, freq_axis, time_axis, fsProc] = d3_ComputeGaitSpectrogram(x, fs, CPI_s)  % CPI request, s
    Seifert_approach = 0;
    if (Seifert_approach == 1)
        processingFs = 2560;   % Hz
        if fs > processingFs
            [p, q]  = rat(processingFs / fs, 1e-6);
            xProc   = resample(x, p, q);
            fsProc  = fs * p / q;
        else
            xProc  = x;
            fsProc = fs;
        end
        win      = hamming(256);   % Seifert et al. (2020): M = 256
        noverlap = 255;              % Seifert et al. (2020): M-1 samples
        nfft     = 2048;             % zero-padded
    else
        xProc  = x;
        fsProc = fs;
        N = 2^(nextpow2(round(CPI_s*fsProc)));
        win      = hamming(N);
        noverlap = N*3/4;                           % Wang uses 75% overlap
        nfft     = N;
    end
    [S, freq_axis, time_axis] = d3_Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);
end

%% Step 3: Subtract background noise
function Pm = d3_SubtractBackgroundNoise(S)
    P = abs(S).^2;                  % power, not magnitude (Wang's P(t,f) is a power quantity)
    Pback = zeros(size(P,1),1);                                 % intialise variable containing background power per Doppler bin
    numberBackgroundSamples = max(1,round(0.20*size(P,2)));
    for frequencyIndex = 1:size(P,1)
        sortedPower = sort(P(frequencyIndex,:),'ascend');                        % sort values in ascending order
        Pback(frequencyIndex) =  mean(sortedPower(1:numberBackgroundSamples));
    end
    backgroundMatrix =  repmat(Pback,1,size(P,2));
    Pm = P - backgroundMatrix;                                                   % removes background noise from signal to get Pm
    Pm(Pm < 0) = eps;
end

%% Step 4: Median Filter Spectogram
function Pmedian = d3_MedianFilterSpectrogram(Pm, medianFrames)
    Pmedian = movmedian(Pm, medianFrames, 1);  % median applied along Doppler dimension
end

%% Step 5: Select highest frequency for Leg
function [legFreq_Hz, isDetected] = d3_SelectLegFrequency(Pmedian, freq_axis, Pth)
    PmedianNorm =  Pmedian/(max(Pmedian(:)) + eps);    % normalise Spectrogram
    PmedianDb =  10*log10(PmedianNorm + eps);          % convert to dB, where maximum value is 0 dB
    aboveThresh = PmedianDb > Pth;         % logical mask, same size as Pmedian
    legFreq_Hz = NaN(1, size(Pmedian,2));    % initialise to NaN
    isDetected = false(1, size(Pmedian,2));
    for t = 1:size(Pmedian,2)
        idx = find(aboveThresh(:,t), 1, 'last');   % freq_axis is ascending, so 'last' = highest frequency
        if ~isempty(idx)
            legFreq_Hz(t) = freq_axis(idx);
            isDetected(t) = true;
        end
    end
end

%% Step 6: Finding frequency Outliers
function [legFreq_interp, isOutlier] = d3_RefineLegFrequencyOutliers(legFreq_Hz, isDetected)   % returns function legFreq_interp and not legFreq_refined
    mu    = mean(legFreq_Hz(isDetected));   % Finding average of leg frequencies where Step 5 found actual detections
    sigma = std(legFreq_Hz(isDetected));    % Finding the standard deviation
    isOutlier = false(size(legFreq_Hz));    % Array size of legFreq_Hz set to false
    isOutlier(isDetected) = abs(legFreq_Hz(isDetected) - mu) > 3*sigma; % Among detected frames, checks how far value is from mean
    legFreq_refined = legFreq_Hz;
    legFreq_refined(isOutlier) = NaN;   % Overwrites frequencies that are outliers
    legFreq_interp = fillmissing(legFreq_refined, 'linear', 'EndValues','nearest');  % fills the missing values using linear interpolation
end

%% Step 7: Smooth F(t) with an average filter
function F_smooth = d3_SmoothLegFrequency(F, avgFrames)
    F_smooth = movmean(F, avgFrames, 'omitnan');
end

%% Step 8: Obtain step instances through time domain peak selection for F (t),
function stepIdx = d3_DetectStepPeaks(F, c)
    N = length(F);                                  % Number of time frames in signal
    isPeak = false(1, N);                           % logical flag for each frame, all start as false
    for t = 2:N-1                                   % skip the first/last frame, both need a neighbour on each side
        if F(t) >= F(t-1) && F(t) >= F(t+1) + c     % rising/flat before t, and drops by at least c after t
            isPeak(t) = true;                       % if t satisfies both conditions, mark it a step
        end
    end
    stepIdx = find(isPeak);                         % convert the logical flags into a list of frame indices
end

%% Step 9
function stepIdx_refined = d3_RemoveCloseSteps(stepIdx, time_axis, F_smooth, minInterval_s)
    stepTimes = time_axis(stepIdx);              % look up the actual time of each candidate step
    keep = true(size(stepIdx));                  % logical flag for each candidate, all start as true
    lastKeptIdx = 1;
    for k = 2:length(stepIdx)                    % start from second candidate
        if (stepTimes(k) - stepTimes(lastKeptIdx)) <= minInterval_s    % too close to the last kept step
            if F_smooth(stepIdx(k)) > F_smooth(stepIdx(lastKeptIdx))
                 keep(lastKeptIdx) = false;
                 lastKeptIdx = k;
            else
                 keep(k) = false;
            end
        else
             lastKeptIdx = k;
        end
    end
    stepIdx_refined = stepIdx(keep);
end
function [S, freq_axis, time_axis] = d3_Zahraa_Spectrogram(x, win, noverlap, nfft, fs)
    x   = x(:).';
    win = win(:).';
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts
    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);
    frame_starts = zeros(1, M);      % keep track of where each frame starts (for time axis later)
    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % grab the chunk and window it
        S(:, m) = fft(segment, nfft).';
    end
    time_axis = (frame_starts - 1) / fs;
    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end

function [S, freq_axis, time_axis] = Zahraa_Spectrogram(x, win, noverlap, nfft, fs)

    x   = x(:).';
    win = win(:).';
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts

    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);

    frame_starts = zeros(1, M);      % keep track of where each frame starts (for time axis later)

    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % grab the chunk and window it
        S(:, m) = fft(segment, nfft).';
    end

    time_axis = (frame_starts - 1 + L/2) / fs;  % labelling the time_axis as the center of the window

    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end
