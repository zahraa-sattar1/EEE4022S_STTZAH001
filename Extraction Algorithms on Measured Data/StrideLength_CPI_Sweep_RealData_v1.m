%% Stride length CPI sweep
% Zahraa Sattar | EEE4022S
%
% Compare CPI settings on one measured walking recording.
% Report stride length estimates and RMSE against the recording reference.

clear all;
close all;

%% Load the recording once, this doesn't change between CPI values
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

if ~exist('fc','var')
    fc = 24e9;   % fallback only if LoadRealRecording didn't set it
end
if ~exist('c','var')
    c = 299792458;   % fallback only if LoadRealRecording didn't set it
end

%% Ground truth (SET THIS for whichever recording you're running)
groundTruthStrideLength_m = 1.3;   % known stride length, in metres

%% CPI values to sweep
CPI_values_s = 0.05:0.05:0.7;   % same set you used for the simulator version
numCPI = length(CPI_values_s);

rmse_m = NaN(1, numCPI);
meanStrideLength_results = NaN(1, numCPI);
percentError_results = NaN(1, numCPI);
numStrides_results = NaN(1, numCPI);
CPI_actual_s = NaN(1, numCPI);

for k = 1:numCPI
    CPI_s = CPI_values_s(k);   % this is what changes each loop, everything below reads it

%% Walking speed extraction (same as RealData_FullAnalysis_v1.m Section 2)
    fc_paper = 5.8e9;
    K = fc / fc_paper;
    f_low0  = 5 * K;
    f_high_min = 50*K;
    maxExpectedSpeed_mps = 4.0;
    f_high0 = maxExpectedSpeed_mps * 2 * fc / c;

    numIterations = 3;
    resultsIter = RunGaitBandRefinement(RdrRxSignal, fs, f_low0, f_high_min, f_high0, numIterations, CPI_s);
    final = resultsIter(end);   % most refined pass

    finalValidFrequencyIndex = find( abs(final.freq_axis) >= final.freq_low &  abs(final.freq_axis) <= final.freq_high);
    energyPerFrame = sum(final.Ssmooth(finalValidFrequencyIndex, :), 1);
    energyBaseline = median(energyPerFrame,'omitnan');
    energySpread = median(abs(energyPerFrame-energyBaseline),'omitnan');
    if exist('thresholdMultiplierOverride','var')
        thresholdMultiplier = thresholdMultiplierOverride;   % lets a sweep script force a specific multiplier
    else
        thresholdMultiplier = 2;   % 's original default
    end
    energyThreshold = max(0,energyBaseline-thresholdMultiplier*energySpread);
    isNoMotion = energyPerFrame < energyThreshold;
    torsoFreq_suppressed = final.torsoFreq_refined_Hz;
    torsoFreq_suppressed(isNoMotion) = 0;

    velocityTrack_mps = torsoFreq_suppressed * c / (2 * fc);
    velocityTrack_mps(isNoMotion) = NaN;

%% Tc / stride duration estimate via torso-peak findpeaks logic (same as Section 3)
    minPeakDistance_s = 0.3;
    minPeakProminence_Hz = 0.25*(max(final.torsoFreq_refined_Hz) - min(final.torsoFreq_refined_Hz));
    [~, peakIdx] = findpeaks(final.torsoFreq_refined_Hz, ...
        'MinPeakDistance', round(minPeakDistance_s / (final.time_axis(2)-final.time_axis(1))), ...
        'MinPeakProminence', minPeakProminence_Hz);
    peakTimes_s = final.time_axis(peakIdx);
    stepIntervals_s = diff(peakTimes_s);
    Tc = 2 * mean(stepIntervals_s);

%% Step recognition extraction (same as Section 4), this is what the step times come from
    freq_high_Yardibi = 350*(fc/1e9)/5.8;

    marginFactor = 1.5;
    avgSamples = floor(fs / (marginFactor * freq_high_Yardibi));
    xAvg = AverageFilterRadarData(RdrRxSignal, avgSamples);

    [S_step, freq_axis_step_full, time_axis_step, fsProc_step] = ComputeGaitSpectrogram_Legs(xAvg, fs, CPI_s);

    nonzeroIndex = find( (freq_axis_step_full > 0) & (freq_axis_step_full <= freq_high_Yardibi));
    freq_axis_step = freq_axis_step_full(nonzeroIndex);
    S_step_leg = S_step(nonzeroIndex,:);
    Pm = SubtractBackgroundNoise(S_step_leg);

    frame_period_s = time_axis_step(2) - time_axis_step(1);
    medianWindow_s = Tc/8;
    medianFrames = max(1, round(medianWindow_s / frame_period_s));
    Pmedian = MedianFilterSpectrogram(Pm, medianFrames);

    if exist('PthOverride','var')
        Pth = PthOverride;
    elseif exist('SNR_dB','var')
        Pth = 10 - SNR_dB;
        Pth = max(Pth, -40);
    else
        Pth = -40;
    end
    [legFreq_Hz, isDetected] = SelectLegFrequency(Pmedian, freq_axis_step, Pth);
    [legFreq_refined, isOutlier] = RefineLegFrequencyOutliers(legFreq_Hz, isDetected);

    smoothWindow_s = Tc/8;
    avgFrames7 = max(1, round(smoothWindow_s / frame_period_s));
    F_smooth = SmoothLegFrequency(legFreq_refined, avgFrames7);

    frequencyNormalization_Hz = max(F_smooth,[],'omitnan');
    normalizedLegFrequency =  F_smooth/frequencyNormalization_Hz;
    DopplerRes_Hz = 1/CPI_s;
    cPeak =  DopplerRes_Hz/freq_high_Yardibi;
    stepIdx8 = find( (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(1:end-2)) & ...
                      (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(3:end) + cPeak) ) + 1;

    minInterval_s = 0.35;
    stepIdx9 = RemoveCloseSteps(stepIdx8, time_axis_step, F_smooth, minInterval_s);
    detectedStepTimes_s = time_axis_step(stepIdx9);   % this is what stride length reads displacement at

%% Stride length (same as Section 5)
    torsoVelocity_m_s = velocityTrack_mps;
    torsoVelocity_m_s(isnan(torsoVelocity_m_s)) = 0;   % zero out the NaN gaps so integrating it doesn't break
    torsoVelocityTime_s = final.time_axis;
    [estimated_stepLengths_m, estimated_strideLengths_m, displacement_m] = EstimateStrideLength(torsoVelocity_m_s, torsoVelocityTime_s, detectedStepTimes_s);

%% Metrics for this CPI
    rmse_m(k) = sqrt(mean((estimated_strideLengths_m - groundTruthStrideLength_m).^2, 'omitnan'));   % pooled across every individual detected stride
    meanStrideLength_results(k) = mean(estimated_strideLengths_m, 'omitnan');
    percentError_results(k) = 100 * (meanStrideLength_results(k) - groundTruthStrideLength_m) / groundTruthStrideLength_m;
    numStrides_results(k) = numel(estimated_strideLengths_m);
    CPI_actual_s(k) = 2^(nextpow2(round(CPI_s*fs))) / fs;

    fprintf('CPI = %.3f s (actual %.3f s) -> mean stride length = %.3f m, RMSE = %.4f m, error = %+.1f%%, n = %d strides\n', ...
        CPI_s, CPI_actual_s(k), meanStrideLength_results(k), rmse_m(k), percentError_results(k), numStrides_results(k));
end

%% Results plot
figure('Name','Stride length CPI sweep (real data)');
subplot(2,1,1);
plot(CPI_actual_s, rmse_m, 'o-');
xlabel('Actual CPI (s)');
ylabel('RMSE (m)');
title('Stride length RMSE vs ground truth');
grid on;

subplot(2,1,2);
plot(CPI_actual_s, meanStrideLength_results, 'o-');
hold on;
yline(groundTruthStrideLength_m, 'r--');
xlabel('Actual CPI (s)');
ylabel('Mean stride length (m)');
legend('Estimated mean stride length','Ground truth','Location','best');
title('Mean stride length vs CPI');
grid on;

[bestRMSE, bestIdx] = min(rmse_m);
bestCPI = CPI_actual_s(bestIdx);

%% Results table
fprintf('\n--- Stride length CPI sweep summary ---\n');
fprintf('%-14s %-14s %-16s %-12s %-10s %-10s\n','CPI req.','CPI actual','Mean length','RMSE','Error %','n strides');
for k = 1:numCPI
    fprintf('%-14.3f %-14.3f %-16.3f %-12.4f %+.1f%%     %d\n', CPI_values_s(k), CPI_actual_s(k), meanStrideLength_results(k), rmse_m(k), percentError_results(k), numStrides_results(k));
end
fprintf('Best CPI = %.3f s (actual), with RMSE = %.4f m\n', bestCPI, bestRMSE);

%% LOCAL FUNCTIONS (everything the walking speed + Tc + step recognition + stride length chain needs)

function [y, f, mask] = IdealBandpassFilter(x, fs, f_low, f_high)
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
        y = y.';
    end
end

function [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram_Torso(x, fs, CPI_s, freq_high)
    xProc  = x;
    fsProc = fs;
    N = 2^(nextpow2(round(CPI_s*fsProc)));
    win      = hamming(N);
    noverlap = N*3/4;
    nfft     = N;
    [S_temp, freq_axis_temp, time_axis] = Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);
    validFrequencyIndex = find( (abs(freq_axis_temp) <= freq_high*3) );
    S = S_temp(validFrequencyIndex,:);
    freq_axis = freq_axis_temp(validFrequencyIndex);
end

function Ssmooth = SmoothSpectrogramNoise(S, avgFrames)
    Ssmooth = movmean(abs(S), avgFrames, 1);
end

function [torsoFreq_Hz, peakVal] = TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high)
    validFrequencyIndex = find( (abs(freq_axis) >= freq_low) & (abs(freq_axis) <= freq_high) );
    Subset_Ssmooth = Ssmooth(validFrequencyIndex,:);
    [peakVal, peakIdx] = max(Subset_Ssmooth, [], 1);
    torsoFreq_Hz = freq_axis(validFrequencyIndex(peakIdx));
end

function [torsoFreq_refined, isOutlier] = RefineOutliers(torsoFreq_Hz)
    mu    = mean(torsoFreq_Hz, 'omitnan');
    sigma = std(torsoFreq_Hz, 'omitnan');
    isOutlier = abs(torsoFreq_Hz - mu) > 3*sigma;
    torsoFreq_refined = torsoFreq_Hz;
    torsoFreq_refined(isOutlier) = NaN;
    torsoFreq_refined = fillmissing(torsoFreq_refined,'linear','EndValues','nearest');
end

function results = RunGaitBandRefinement(x, fs, freq_low0, f_high_min, freq_high0, numIterations, CPI_s)
    freq_low  = freq_low0;
    freq_high = freq_high0;
    results = struct('freq_low', {}, 'freq_high', {}, 'torsoFreq_Hz', {}, ...
        'torsoFreq_refined_Hz', {}, 'centreFreq_Hz', {});
    for iter = 1:numIterations
        xFilt = IdealBandpassFilter(x, fs, freq_low, freq_high);
        [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram_Torso(xFilt, fs, CPI_s, freq_high);
        Power_S = abs(S).^2;
        Ssmooth = SmoothSpectrogramNoise(Power_S, 21);
        torsoFreq_Hz = TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high);
        torsoFreq_refined_Hz = RefineOutliers(torsoFreq_Hz);
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
            freq_low = max(0.5*centreFreq_Hz,freq_low0);
            freq_high = max(2.0*centreFreq_Hz,f_high_min);
        end
    end
end

function xAvg = AverageFilterRadarData(x, avgSamples)
    isRowInput = isrow(x);
    x = x(:).';
    xAvg = movmean(x, avgSamples);
    if ~isRowInput
        xAvg = xAvg.';
    end
end

function [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram_Legs(x, fs, CPI_s)
    xProc  = x;
    fsProc = fs;
    N = 2^(nextpow2(round(CPI_s*fsProc)));
    win      = hamming(N);
    noverlap = N*3/4;
    nfft     = N;
    [S, freq_axis, time_axis] = Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);
end

function Pm = SubtractBackgroundNoise(S)
    P = abs(S).^2;
    Pback = zeros(size(P,1),1);
    numberBackgroundSamples = max(1,round(0.20*size(P,2)));
    for frequencyIndex = 1:size(P,1)
        sortedPower = sort(P(frequencyIndex,:),'ascend');
        Pback(frequencyIndex) =  mean(sortedPower(1:numberBackgroundSamples));
    end
    backgroundMatrix =  repmat(Pback,1,size(P,2));
    Pm = P - backgroundMatrix;
    Pm(Pm < 0) = eps;
end

function Pmedian = MedianFilterSpectrogram(Pm, medianFrames)
    Pmedian = movmedian(Pm, medianFrames, 1);
end

function [legFreq_Hz, isDetected] = SelectLegFrequency(Pmedian, freq_axis, Pth)
    PmedianNorm =  Pmedian/(max(Pmedian(:)) + eps);
    PmedianDb =  10*log10(PmedianNorm + eps);
    aboveThresh = PmedianDb > Pth;
    legFreq_Hz = NaN(1, size(Pmedian,2));
    isDetected = false(1, size(Pmedian,2));
    for t = 1:size(Pmedian,2)
        idx = find(aboveThresh(:,t), 1, 'last');
        if ~isempty(idx)
            legFreq_Hz(t) = freq_axis(idx);
            isDetected(t) = true;
        end
    end
end

function [legFreq_interp, isOutlier] = RefineLegFrequencyOutliers(legFreq_Hz, isDetected)
    mu    = mean(legFreq_Hz(isDetected));
    sigma = std(legFreq_Hz(isDetected));
    isOutlier = false(size(legFreq_Hz));
    isOutlier(isDetected) = abs(legFreq_Hz(isDetected) - mu) > 3*sigma;
    legFreq_refined = legFreq_Hz;
    legFreq_refined(isOutlier) = NaN;
    legFreq_interp = fillmissing(legFreq_refined, 'linear', 'EndValues','nearest');
end

function F_smooth = SmoothLegFrequency(F, avgFrames)
    F_smooth = movmean(F, avgFrames, 'omitnan');
end

function stepIdx_refined = RemoveCloseSteps(stepIdx, time_axis, F_smooth, minInterval_s)
    stepTimes = time_axis(stepIdx);
    keep = true(size(stepIdx));
    lastKeptIdx = 1;
    for k = 2:length(stepIdx)
        if (stepTimes(k) - stepTimes(lastKeptIdx)) <= minInterval_s
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

function [stepLengths_m, strideLengths_m, displacement_m] = EstimateStrideLength(velocity_m_s, velocityTime_s, stepTimes_s)
    displacement_m = cumtrapz(velocityTime_s, velocity_m_s);
    displacement_atSteps = interp1(velocityTime_s, displacement_m, stepTimes_s, 'linear');
    stepLengths_m = abs(diff(displacement_atSteps));
    strideLengths_m = abs(displacement_atSteps(3:end) - displacement_atSteps(1:end-2));
end

function [S, freq_axis, time_axis] = Zahraa_Spectrogram(x, win, noverlap, nfft, fs)
    x   = x(:).';
    win = win(:).';
    L = length(win);
    N = length(x);
    hop = L - noverlap;
    M = floor((N - L) / hop) + 1;
    S = zeros(nfft, M);
    frame_starts = zeros(1, M);
    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;
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
