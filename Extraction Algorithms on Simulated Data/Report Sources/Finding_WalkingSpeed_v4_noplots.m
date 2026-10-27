%% Historical walking-speed extractor
% Zahraa Sattar | EEE4022S
%
% Input source for the Chapter 4 experiment drivers.
% The drivers create and modify working copies for each comparison.

clear all;
close all;
RadarSimulator_v11_torso_legs_arms_head_buffered_noplots;
fc_paper = 5.8e9;
if ~exist('CPI_s','var')
    CPI_s = 0.1;
end
N_actual = 2^(nextpow2(round(CPI_s*fs)));
CPI_actual_s_WS = N_actual / fs;
fprintf('Requested CPI %.4f:\n', CPI_s);
fprintf('Actual CPI %.4f:\n', CPI_actual_s_WS);
K = fc / fc_paper;
f_low0  = 5 * K;
f_high_min = 50*K;
maxExpectedSpeed_mps = 4.0;
f_high0 = maxExpectedSpeed_mps * 2 * fc / c;
fprintf('K = %.6f, initial band = [%.4f, %.4f] Hz\n', K, f_low0, f_high0);
numIterations = 3;
resultsIter = RunGaitBandRefinement(RdrRxSignal, fs, f_low0, f_high_min, f_high0, numIterations, CPI_s);
for iter = 1:numIterations
    fprintf('Iteration %d: band = [%.2f, %.2f] Hz, centre = %.4f Hz\n', ...
        iter, resultsIter(iter).freq_low, resultsIter(iter).freq_high, resultsIter(iter).centreFreq_Hz);
end
final = resultsIter(end);
finalValidFrequencyIndex = find( abs(final.freq_axis) >= final.freq_low &  abs(final.freq_axis) <= final.freq_high);
energyPerFrame = sum(final.Ssmooth(finalValidFrequencyIndex, :), 1);
energyBaseline = median(energyPerFrame,'omitnan');
energySpread = median(abs(energyPerFrame-energyBaseline),'omitnan');
energyThreshold = max(0,energyBaseline-2*energySpread);
isNoMotion = energyPerFrame < energyThreshold;
torsoFreq_suppressed = final.torsoFreq_refined_Hz;
torsoFreq_suppressed(isNoMotion) = 0;
c = 3e8;
velocityTrack_mps = torsoFreq_suppressed * c / (2 * fc);
velocityTrack_mps(isNoMotion) = NaN;
movingFrames = torsoFreq_suppressed ~= 0;
meanSpeed_mps = mean(velocityTrack_mps(movingFrames), 'omitnan');
fprintf('\nExtracted mean walking speed: %.4f m/s\n', meanSpeed_mps);
if exist('v0','var')
    percentError = 100 * (meanSpeed_mps - v0) / v0;
    fprintf('Ground-truth mean speed v0        : %.3f m/s\n', v0);
    fprintf('Extracted mean speed              : %.3f m/s\n', meanSpeed_mps);
    fprintf('Error vs. ground truth            : %+.1f %%\n', percentError);
end
if all(torsoFreq_suppressed == 0)
    warning('Every frame was suppressed as no-motion, energyThreshold is set too high.');
end
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
function [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram(x, fs, CPI_s, freq_high)
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
    N = 2^(nextpow2(round(CPI_s*fsProc)));
    win      = hamming(N);
    noverlap = N-1;
    nfft     = N*8;
else
        xProc  = x;
        fsProc = fs;
        N = 2^(nextpow2(round(CPI_s*fsProc)));
        win      = hamming(N);
        noverlap = N*3/4;
        nfft     = N;
end
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
        [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram(xFilt, fs, CPI_s, freq_high);
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
function torsoFreq_suppressed = SuppressNoMotion(torsoFreq_Hz, Ssmooth, energyThreshold)
    energyPerFrame = sum(Ssmooth, 1);
    isNoMotion = energyPerFrame < energyThreshold;
    torsoFreq_suppressed = torsoFreq_Hz;
    torsoFreq_suppressed(isNoMotion) = 0;
end
function [velocityTrack_mps, meanSpeed_mps] = ComputeWalkingSpeed(torsoFreq_Hz, ft_Hz, c)
    velocityTrack_mps = torsoFreq_Hz * c / (2 * ft_Hz);
    movingFrames = torsoFreq_Hz ~= 0;
    meanSpeed_mps = mean(velocityTrack_mps(movingFrames), 'omitnan');
end
