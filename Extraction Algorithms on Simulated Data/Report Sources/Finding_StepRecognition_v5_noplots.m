%% Historical step detector
% Zahraa Sattar | EEE4022S
%
% Input source for the Chapter 4 experiment drivers.
% The drivers create and modify working copies for each comparison.

RadarSimulator_v11_torso_legs_arms_head_buffered_noplots;
freq_high_Yardibi = 350*24/5.8;
c = 299792458;
fc = 24e9;
lambda = c/fc;
marginFactor = 1.5;
avgSamples = floor(fs / (marginFactor * freq_high_Yardibi));
xAvg = AverageFilterRadarData(RdrRxSignal, avgSamples);
sampleWindow = 1:500;
fprintf('Std of raw signal (this segment): %.4f\n', std(real(RdrRxSignal(sampleWindow))));
fprintf('Std of averaged signal (this segment): %.4f\n', std(real(xAvg(sampleWindow))));
if ~exist('CPI_s','var')
    CPI_s = 0.1;
end
[S_step, freq_axis_step_full, time_axis_step, fsProc_step] = ComputeGaitSpectrogram(xAvg, fs, CPI_s);
CPI_actual_s = length(freq_axis_step_full) / fsProc_step;
fprintf('Requested CPI %.4f:\n', CPI_s);
fprintf('Actual CPI %.4f:\n', CPI_actual_s);
nonzeroIndex = find( (freq_axis_step_full > 0) & (freq_axis_step_full <= freq_high_Yardibi));
freq_axis_step = freq_axis_step_full(nonzeroIndex);
S_step_leg = S_step(nonzeroIndex,:);
P_step2 = abs(S_step_leg).^2;
Pm = SubtractBackgroundNoise(S_step_leg);
fprintf('Check: Pm never negative -> %d\n', all(Pm(:) >= 0));
fprintf('Check: Pm never exceeds original power -> %d\n', all(Pm(:) <= P_step2(:) + eps));
frame_period_s = time_axis_step(2) - time_axis_step(1);
medianWindow_s = Tc/8;
medianFrames = max(1, round(medianWindow_s / frame_period_s));
Pmedian = MedianFilterSpectrogram(Pm, medianFrames);
noiseBinIdx = 1;
fprintf('Std at edge frequency bin BEFORE median filter: %.4f\n', std(Pm(noiseBinIdx,:)));
fprintf('Std at edge frequency bin AFTER  median filter: %.4f\n', std(Pmedian(noiseBinIdx,:)));
fprintf('Pmedian range (dB): [%.1f, %.1f]\n', ...
    10*log10(min(Pmedian(Pmedian>0)) + eps), 10*log10(max(Pmedian(:)) + eps));
if exist('PthOverride','var')
    Pth = PthOverride;
elseif exist('SNR_dB','var')
    Pth = 10 - SNR_dB;
    Pth = max(Pth, -40);
else
    Pth = -25;
end
[legFreq_Hz, isDetected] = SelectLegFrequency(Pmedian, freq_axis_step, Pth);
fprintf('Frames with a leg-frequency detection: %d / %d\n', sum(isDetected), numel(isDetected));
[legFreq_refined, isOutlier] = RefineLegFrequencyOutliers(legFreq_Hz, isDetected);
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
F_smooth = SmoothLegFrequency(legFreq_refined, avgFrames7);
frequencyNormalization_Hz = max(F_smooth,[],'omitnan');
normalizedLegFrequency =  F_smooth/frequencyNormalization_Hz;
DopplerRes_Hz = 1/CPI_actual_s;
cPeak =  DopplerRes_Hz/freq_high_Yardibi;
if ~exist('cPeakScale','var')
    cPeakScale = 0.5;
end
cPeak = cPeak*cPeakScale;
stepIdx8 = find( (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(1:end-2)) & ...
                  (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(3:end) + cPeak) ) + 1;
fprintf('Step 8 raw peaks: %d\n', numel(stepIdx8));
minInterval_s = 0.35;
stepIdx9 = RemoveCloseSteps(stepIdx8, time_axis_step, F_smooth, minInterval_s);
detectedStepTimes_s = time_axis_step(stepIdx9);
fprintf('Step 9 final step count: %d\n', numel(stepIdx9));
if ~exist('stepGapFill','var')
    stepGapFill = 1;
end
stepIdx9_raw = stepIdx9;
detectedStepTimes_raw_s = time_axis_step(stepIdx9_raw);
if stepGapFill == 1
    stepIdx9 = RecoverMissedSteps(stepIdx9, time_axis_step, F_smooth, minInterval_s);
end
fprintf('Step 9b steps put back: %d\n', numel(stepIdx9) - numel(stepIdx9_raw));
if ~exist('stepRefineMode','var')
    stepRefineMode = 'centroid';
end
stepFrameTimes_s = time_axis_step(stepIdx9);
detectedStepTimes_s = RefineStepTimes(time_axis_step, F_smooth, stepIdx9, Tc, stepRefineMode);
stepSpeeds_m_s =  F_smooth(stepIdx9)*lambda/2;
stepDurations_s = diff(detectedStepTimes_s);
if numel(detectedStepTimes_s) >= 3
    strideDurations_s = detectedStepTimes_s(3:end) - detectedStepTimes_s(1:end-2);
else
    strideDurations_s = [];
end
function xAvg = AverageFilterRadarData(x, avgSamples)
    isRowInput = isrow(x);
    x = x(:).';
    xAvg = movmean(x, avgSamples);
    if ~isRowInput
        xAvg = xAvg.';
    end
end
function [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram(x, fs, CPI_s)
    Seifert_approach = 0;
    if (Seifert_approach == 1)
        processingFs = 2560;
        if fs > processingFs
            [p, q]  = rat(processingFs / fs, 1e-6);
            xProc   = resample(x, p, q);
            fsProc  = fs * p / q;
        else
            xProc  = x;
            fsProc = fs;
        end
        win      = hamming(256);
        noverlap = 255;
        nfft     = 2048;
    else
        xProc  = x;
        fsProc = fs;
        N = 2^(nextpow2(round(CPI_s*fsProc)));
        win      = hamming(N);
        noverlap = N*3/4;
        nfft     = N;
    end
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
function stepIdx = DetectStepPeaks(F, c)
    N = length(F);
    isPeak = false(1, N);
    for t = 2:N-1
        if F(t) >= F(t-1) && F(t) >= F(t+1) + c
            isPeak(t) = true;
        end
    end
    stepIdx = find(isPeak);
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
function stepIdx_out = RecoverMissedSteps(stepIdx, time_axis, F_smooth, minInterval_s)
    stepIdx_out = stepIdx;
    for pass = 1:3
        t = time_axis(stepIdx_out);
        iv = diff(t);
        if numel(iv) < 3
            break;
        end
        g = find(iv > 1.6*median(iv), 1);
        if isempty(g)
            break;
        end
        a = stepIdx_out(g);
        b = stepIdx_out(g+1);
        cand = [];
        for q = a+1:b-1
            isTop = F_smooth(q) >= F_smooth(q-1) && F_smooth(q) >= F_smooth(q+1);
            farEnough = (time_axis(q)-time_axis(a)) > minInterval_s && (time_axis(b)-time_axis(q)) > minInterval_s;
            if isTop && farEnough
                cand(end+1) = q;                                        %#ok<AGROW>
            end
        end
        if isempty(cand)
            break;
        end
        [~, m] = max(F_smooth(cand));
        stepIdx_out = sort([stepIdx_out, cand(m)]);
    end
end
function t = RefineStepTimes(time_axis, F, stepIdx, Tc, mode)
    t = time_axis(stepIdx);
    if strcmp(mode, 'none') || isempty(stepIdx) || numel(time_axis) < 2
        return;
    end
    dt = time_axis(2) - time_axis(1);
    N = numel(F);
    wAll = round(Tc/4/dt);
    for k = 1:numel(stepIdx)
        i = stepIdx(k);
        if strcmp(mode, 'parabola')
            if i > 1 && i < N
                y1 = F(i-1); y2 = F(i); y3 = F(i+1);
                den = y1 - 2*y2 + y3;
                if den < 0
                    d = max(-1, min(1, 0.5*(y1-y3)/den));
                    t(k) = time_axis(i) + d*dt;
                end
            end
        else
            wBack = wAll;
            wFwd = wAll;
            if k > 1
                wBack = min(wBack, floor((i - stepIdx(k-1))/2));
            end
            if k < numel(stepIdx)
                wFwd = min(wFwd, floor((stepIdx(k+1) - i)/2));
            end
            a = max(1, i - wBack);
            b = min(N, i + wFwd);
            seg = F(a:b);
            pk = max(seg);
            base = min(F(max(1,i-wAll):min(N,i+wAll)));
            wgt = max(seg - (pk - 0.5*(pk - base)), 0);
            wgt(~isfinite(wgt)) = 0;
            if sum(wgt) > 0
                t(k) = sum(time_axis(a:b).*wgt) / sum(wgt);
            end
        end
    end
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
