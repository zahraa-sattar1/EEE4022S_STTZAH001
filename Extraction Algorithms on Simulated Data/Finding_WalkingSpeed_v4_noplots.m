%% Walking speed on simulated data
% Zahraa Sattar | EEE4022S
%
% Estimate walking speed from the torso Doppler track.
% The existing simulator workspace supplies the signal, and the velocity track is retained for stride length.

fc_paper = 5.8e9;              % Yardibi et al. (2011)'s radar carrier, Hz

if ~exist('CPI_s','var')       % retain a CPI supplied by the caller
    CPI_s = 0.1;                % Introduce CPI input variable
end

N_actual = 2^(nextpow2(round(CPI_s*fs)));
CPI_actual_s_WS = N_actual / fs;        % CPI_actual = N * (1/fs)
fprintf('Requested CPI %.4f:\n', CPI_s);
fprintf('Actual CPI %.4f:\n', CPI_actual_s_WS);

K = fc / fc_paper;              % fc comes from the simulator workspace
f_low0  = 5 * K;    % multiply low and high frequencies by K scaler
f_high_min = 50*K;
maxExpectedSpeed_mps = 4.0;             % covers brisk walking/light jogging
f_high0 = maxExpectedSpeed_mps * 2 * fc / c;
fprintf('K = %.6f, initial band = [%.4f, %.4f] Hz\n', K, f_low0, f_high0);

%% 2. Steps 1-6
numIterations = 3;
resultsIter = RunGaitBandRefinement(RdrRxSignal, fs, f_low0, f_high_min, f_high0, numIterations, CPI_s); % f_high_min is also an input
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
        y = y.';        % restores the caller's original orientation
    end
end

%% Step 2: Compute Spectrogram
function [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram(x, fs, CPI_s, freq_high)  % included CPI_s
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
        noverlap = N*3/4;                           % 75% overlap in this implementation
        nfft     = N;
end
    [S_temp, freq_axis_temp, time_axis] = Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);
    validFrequencyIndex = find( (abs(freq_axis_temp) <= freq_high*3) );
    S = S_temp(validFrequencyIndex,:);                  % extract subset of spectogram
    freq_axis = freq_axis_temp(validFrequencyIndex);
end

%% Step 3: Reduce spectral noise with average filter
function Ssmooth = SmoothSpectrogramNoise(S, avgFrames)
      Ssmooth = movmean(abs(S), avgFrames, 1);          % Apply moving average across Doppler bins (i.e dimension 1)
end

%% Step 4: torso Doppler frequency = frequency of peak
function [torsoFreq_Hz, peakVal] = TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high)
    validFrequencyIndex = find( (abs(freq_axis) >= freq_low) & (abs(freq_axis) <= freq_high) );
    Subset_Ssmooth = Ssmooth(validFrequencyIndex,:);
    [peakVal, peakIdx] = max(Subset_Ssmooth, [], 1);
    torsoFreq_Hz = freq_axis(validFrequencyIndex(peakIdx));
end

%% Step 5: Outlier rejection
function [torsoFreq_refined, isOutlier] = RefineOutliers(torsoFreq_Hz)
    mu    = mean(torsoFreq_Hz, 'omitnan');
    sigma = std(torsoFreq_Hz, 'omitnan');
    isOutlier = abs(torsoFreq_Hz - mu) > 3*sigma;
    torsoFreq_refined = torsoFreq_Hz;
    torsoFreq_refined(isOutlier) = NaN;
    torsoFreq_refined = fillmissing(torsoFreq_refined,'linear','EndValues','nearest'); % Fill in NaNs with 'linear' interpolated values
end

%% Step 6: Iterative band refinement (repeats Steps 1-5)
function results = RunGaitBandRefinement(x, fs, freq_low0, f_high_min, freq_high0, numIterations, CPI_s)
    freq_low  = freq_low0;
    freq_high = freq_high0;
    results = struct('freq_low', {}, 'freq_high', {}, 'torsoFreq_Hz', {}, ...
        'torsoFreq_refined_Hz', {}, 'centreFreq_Hz', {});
    for iter = 1:numIterations
        xFilt = IdealBandpassFilter(x, fs, freq_low, freq_high);
        [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram(xFilt, fs, CPI_s, freq_high);
        Power_S = abs(S).^2;                                    % Compute Power of Spectrogram (linear units)
        Ssmooth = SmoothSpectrogramNoise(Power_S, 21);          % Smooth the Power of Spetrogram (linear units)
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
            freq_low = max(0.5*centreFreq_Hz,freq_low0);       % 0.5*MeanGaitFrequency (see Yardibi paper)
            freq_high = max(2.0*centreFreq_Hz,f_high_min);  % 2*MeanGaitFrequency (see Yardibi paper)
        end
    end
end

%% Step 7: Suppress no-motion frames
function torsoFreq_suppressed = SuppressNoMotion(torsoFreq_Hz, Ssmooth, energyThreshold)
    energyPerFrame = sum(Ssmooth, 1);      % total energy across all bins, per frame
    isNoMotion = energyPerFrame < energyThreshold;
    torsoFreq_suppressed = torsoFreq_Hz;
    torsoFreq_suppressed(isNoMotion) = 0;   % wherever isNoMotion is true, that position is set to 0
end

%% Step 8: Compute walking speed
function [velocityTrack_mps, meanSpeed_mps] = ComputeWalkingSpeed(torsoFreq_Hz, ft_Hz, c)
    velocityTrack_mps = torsoFreq_Hz * c / (2 * ft_Hz);
    movingFrames = torsoFreq_Hz ~= 0;
    meanSpeed_mps = mean(velocityTrack_mps(movingFrames), 'omitnan');   % averages speed that are not equal to 0
end


function [S, freq_axis, time_axis] = Zahraa_Spectrogram(x, win, noverlap, nfft, fs)

    x   = x(:).';       % row vector without conjugation
    win = win(:).';     % row window without conjugation
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts

    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);

    frame_starts = zeros(1, M);      % frame start indices for the time axis

    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % window the current segment
        S(:, m) = fft(segment, nfft).';
    end

    time_axis = (frame_starts - 1 + L/2) / fs;  % frame-centre time, s

    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end
