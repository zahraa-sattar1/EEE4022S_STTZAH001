%% Step detector input for the variant comparison
% Zahraa Sattar | EEE4022S
%
% Retain the original frame-based step rule for the report comparison.
% The experiment driver creates and patches its own working copy.

RadarSimulator_v11_torso_legs_arms_head_buffered_noplots;

%% Step 0 - define input parameters
freq_high_Yardibi = 300*24/5.8;

c = 299792458;
fc = 24e9;
lambda = c/fc;

%% Step 1: average filter on raw radar data
marginFactor = 1.5;
avgSamples = floor(fs / (marginFactor * freq_high_Yardibi));
xAvg = AverageFilterRadarData(RdrRxSignal, avgSamples);
sampleWindow = 1:500;   % zoom into a short segment to actually see the smoothing
fprintf('Std of raw signal (this segment): %.4f\n', std(real(RdrRxSignal(sampleWindow))));
fprintf('Std of averaged signal (this segment): %.4f\n', std(real(xAvg(sampleWindow))));

%% Step 2: compute spectrogram
if ~exist('CPI_s','var')
    CPI_s = 0.1;   % Window-length parameter, matching the walking-speed algorithm
end
[S_step, freq_axis_step_full, time_axis_step, fsProc_step] = ComputeGaitSpectrogram(xAvg, fs, CPI_s); % change freq_axis_step to freq_axis_step_full

CPI_actual_s = length(freq_axis_step_full) / fsProc_step;     % CPI_actual_s = N*(1/fs)
fprintf('Requested CPI %.4f:\n', CPI_s);
fprintf('Actual CPI %.4f:\n', CPI_actual_s);

%% Step 3: background noise subtraction
nonzeroIndex = find( (freq_axis_step_full > 0) & (freq_axis_step_full <= freq_high_Yardibi));
freq_axis_step = freq_axis_step_full(nonzeroIndex);
S_step_leg = S_step(nonzeroIndex,:);
P_step2 = abs(S_step_leg).^2;
Pm = SubtractBackgroundNoise(S_step_leg);
fprintf('Check: Pm never negative -> %d\n', all(Pm(:) >= 0));
fprintf('Check: Pm never exceeds original power -> %d\n', all(Pm(:) <= P_step2(:) + eps));

%% Step 4: median filter on spectral data
frame_period_s = time_axis_step(2) - time_axis_step(1);   % actual spacing between spectrogram frames
medianWindow_s = Tc/8;
medianFrames = max(1, round(medianWindow_s / frame_period_s));
Pmedian = MedianFilterSpectrogram(Pm, medianFrames);
noiseBinIdx = 1;
fprintf('Std at edge frequency bin BEFORE median filter: %.4f\n', std(Pm(noiseBinIdx,:)));
fprintf('Std at edge frequency bin AFTER  median filter: %.4f\n', std(Pmedian(noiseBinIdx,:)));
fprintf('Pmedian range (dB): [%.1f, %.1f]\n', ...
    10*log10(min(Pmedian(Pmedian>0)) + eps), 10*log10(max(Pmedian(:)) + eps));

%% Step 5: Select Leg Frequency

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

%% Step 6
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

%% Step 8
frequencyNormalization_Hz = max(F_smooth,[],'omitnan');
normalizedLegFrequency =  F_smooth/frequencyNormalization_Hz;   % normalise leg frequency vector so maximum value is 1

DopplerRes_Hz = 1/CPI_actual_s;

cPeak =  DopplerRes_Hz/freq_high_Yardibi;
if ~exist('cPeakScale','var')
    cPeakScale = 0.5;
end
cPeak = cPeak*cPeakScale;
stepIdx8 = find( (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(1:end-2)) & ...
                  (normalizedLegFrequency(2:end-1) >= normalizedLegFrequency(3:end) + cPeak) ) + 1;
fprintf('Step 8 raw peaks: %d\n', numel(stepIdx8));

%% Step 9
minInterval_s = 0.35;   % Wang et al., 2014
stepIdx9 = RemoveCloseSteps(stepIdx8, time_axis_step, F_smooth, minInterval_s);
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


%% FUNCTIONS

%% Step 1: Average Filter
function xAvg = AverageFilterRadarData(x, avgSamples)
    isRowInput = isrow(x);
    x = x(:).';
    xAvg = movmean(x, avgSamples);
    if ~isRowInput
        xAvg = xAvg.';
    end
end

%% Step 2: Compute Spectogram
function [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram(x, fs, CPI_s)  % added CPI_s
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
    [S, freq_axis, time_axis] = Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);
end

%% Step 3: Subtract background noise
function Pm = SubtractBackgroundNoise(S)
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
function Pmedian = MedianFilterSpectrogram(Pm, medianFrames)
    Pmedian = movmedian(Pm, medianFrames, 1);  % median applied along Doppler dimension
end

%% Step 5: Select highest frequency for Leg
function [legFreq_Hz, isDetected] = SelectLegFrequency(Pmedian, freq_axis, Pth)
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
function [legFreq_interp, isOutlier] = RefineLegFrequencyOutliers(legFreq_Hz, isDetected)   % returns function legFreq_interp and not legFreq_refined
    mu    = mean(legFreq_Hz(isDetected));   % Finding average of leg frequencies where Step 5 found actual detections
    sigma = std(legFreq_Hz(isDetected));    % Finding the standard deviation
    isOutlier = false(size(legFreq_Hz));    % Array size of legFreq_Hz set to false
    isOutlier(isDetected) = abs(legFreq_Hz(isDetected) - mu) > 3*sigma; % Among detected frames, checks how far value is from mean
    legFreq_refined = legFreq_Hz;
    legFreq_refined(isOutlier) = NaN;   % Overwrites frequencies that are outliers
    legFreq_interp = fillmissing(legFreq_refined, 'linear', 'EndValues','nearest');  % fills the missing values using linear interpolation
end

%% Average filter
function F_smooth = SmoothLegFrequency(F, avgFrames)
    F_smooth = movmean(F, avgFrames, 'omitnan');
end

%% Step 8: Obtain step instances through time domain peak selection for F (t),
function stepIdx = DetectStepPeaks(F, c)
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
function stepIdx_refined = RemoveCloseSteps(stepIdx, time_axis, F_smooth, minInterval_s)
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
    time_axis = (frame_starts - 1) / fs;
    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end
