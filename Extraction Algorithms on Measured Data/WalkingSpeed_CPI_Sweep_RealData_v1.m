%% Walking speed CPI sweep
% Zahraa Sattar | EEE4022S
%
% Compare CPI settings on one measured walking recording.
% Report walking speed estimates and RMSE against the recording reference.

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

%% Ground truth (SET THESE for whichever recording running)
walkedDistance_m = 5.2;   % distance actually walked, in metres
recordingDuration_s = length(RdrRxSignal)/fs;   % length of the recording in seconds
groundTruthSpeed_mps = walkedDistance_m / recordingDuration_s;   % the single ground truth speed we're sweeping against

%% CPI values to sweep
CPI_values_s = 0.05:0.05:0.5;   % range of CPI values being swept
numCPI = length(CPI_values_s);

rmse_results = NaN(1, numCPI);
meanSpeed_results = NaN(1, numCPI);
percentError_results = NaN(1, numCPI);

for k = 1:numCPI
    CPI_s = CPI_values_s(k);

%% Walking speed extraction (same as RealData_FullAnalysis_v1.m Section 2)
    fc_paper = 5.8e9;              % Yardibi et al. (2011)'s radar carrier, Hz
    K = fc / fc_paper;              % fc comes from the loaded recording
    f_low0  = 5 * K;    % multiply low and high frequencies by K scaler
    f_high_min = 50*K;
    maxExpectedSpeed_mps = 4.0;             % covers brisk walking/light jogging
    f_high0 = maxExpectedSpeed_mps * 2 * fc / c;   % converts max expected speed into a Doppler frequency

    numIterations = 3;
    resultsIter = RunGaitBandRefinement(RdrRxSignal, fs, f_low0, f_high_min, f_high0, numIterations, CPI_s);
    final = resultsIter(end);   % most refined pass

    finalValidFrequencyIndex = find( abs(final.freq_axis) >= final.freq_low &  abs(final.freq_axis) <= final.freq_high);
    energyPerFrame = sum(final.Ssmooth(finalValidFrequencyIndex, :), 1);
    energyBaseline = median(energyPerFrame,'omitnan');
    energySpread = median(abs(energyPerFrame-energyBaseline),'omitnan');
    energyThreshold = max(0,energyBaseline-2*energySpread);
    isNoMotion = energyPerFrame < energyThreshold;
    torsoFreq_suppressed = final.torsoFreq_refined_Hz;
    torsoFreq_suppressed(isNoMotion) = 0;

    velocityTrack_mps = torsoFreq_suppressed * c / (2 * fc);
    velocityTrack_mps(isNoMotion) = NaN;
    movingFrames = torsoFreq_suppressed ~= 0;
    meanSpeed_mps = mean(velocityTrack_mps(movingFrames), 'omitnan');

%% Metrics for this CPI
    rmse_results(k) = sqrt(mean((velocityTrack_mps(movingFrames) - groundTruthSpeed_mps).^2, 'omitnan'));
    meanSpeed_results(k) = meanSpeed_mps;
    percentError_results(k) = 100 * (meanSpeed_mps - groundTruthSpeed_mps) / groundTruthSpeed_mps;

    fprintf('CPI = %.3f s -> mean speed = %.3f m/s, RMSE = %.3f m/s, error = %+.1f%%\n', ...
        CPI_s, meanSpeed_mps, rmse_results(k), percentError_results(k));
end

%% Results plot
figure('Name','Walking speed CPI sweep (real data)');
subplot(2,1,1);
plot(CPI_values_s, rmse_results, 'o-');
xlabel('CPI (s)');
ylabel('RMSE (m/s)');
title('Velocity track RMSE vs ground truth speed');
grid on;

subplot(2,1,2);
plot(CPI_values_s, meanSpeed_results, 'o-');
hold on;
yline(groundTruthSpeed_mps, 'r--');
xlabel('CPI (s)');
ylabel('Mean speed (m/s)');
legend('Estimated mean speed','Ground truth speed','Location','best');
title('Mean walking speed vs CPI');
grid on;

%% Results table
fprintf('\n--- Walking speed CPI sweep summary ---\n');
fprintf('%-10s %-15s %-12s %-10s\n','CPI (s)','Mean speed','RMSE','Error %');
for k = 1:numCPI
    fprintf('%-10.3f %-15.3f %-12.3f %+.1f%%\n', CPI_values_s(k), meanSpeed_results(k), rmse_results(k), percentError_results(k));
end

%% LOCAL FUNCTIONS (only the ones the walking speed pipeline actually needs)

function [y, f, mask] = IdealBandpassFilter(x, fs, f_low, f_high)
    isRowInput = isrow(x);   % remember the original shape so we can give it back at the end
    x = x(:).';   % force it into a row vector so the fft math below is consistent
    N = length(x);   % number of samples
    X  = fft(x);   % into the frequency domain
    Xs = fftshift(X);
    f = (-floor(N/2):ceil(N/2)-1) * (fs / N);   % frequency axis matching that shifted fft
    mask = (abs(f) >= f_low) & (abs(f) <= f_high);   % keep only the frequencies inside our band, zero out the rest
    Xf = Xs .* mask;   % apply the mask
    Xu = ifftshift(Xf);   % shift back to normal fft order
    y  = ifft(Xu);   % back to the time domain, this is the filtered signal
    if ~isRowInput
        y = y.';   % put it back in its original column shape if that's what came in
    end
end

function [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram_Torso(x, fs, CPI_s, freq_high)
    xProc  = x;
    fsProc = fs;
    N = 2^(nextpow2(round(CPI_s*fsProc)));
    win      = hamming(N);   % hamming window to reduce spectral leakage
    noverlap = N*3/4;   % 75% overlap between windows
    nfft     = N;   % fft length matches the window
    [S_temp, freq_axis_temp, time_axis] = Zahraa_Spectrogram(xProc, win, noverlap, nfft, fsProc);   % run the actual spectrogram
    validFrequencyIndex = find( (abs(freq_axis_temp) <= freq_high*3) );   % only keep frequencies up to 3x freq_high
    S = S_temp(validFrequencyIndex,:);   % cut down to that range
    freq_axis = freq_axis_temp(validFrequencyIndex);   % cut down to that range too
end

function Ssmooth = SmoothSpectrogramNoise(S, avgFrames)
    Ssmooth = movmean(abs(S), avgFrames, 1);   % moving average across time frames to smooth out noise in the spectrogram
end

function [torsoFreq_Hz, peakVal] = TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high)
    validFrequencyIndex = find( (abs(freq_axis) >= freq_low) & (abs(freq_axis) <= freq_high) );   % only look within the current torso band
    Subset_Ssmooth = Ssmooth(validFrequencyIndex,:);   % restrict the spectrogram to that band
    [peakVal, peakIdx] = max(Subset_Ssmooth, [], 1);
    torsoFreq_Hz = freq_axis(validFrequencyIndex(peakIdx));   % convert that index back into an actual frequency value
end

function [torsoFreq_refined, isOutlier] = RefineOutliers(torsoFreq_Hz)
    mu    = mean(torsoFreq_Hz, 'omitnan');   % average torso frequency across the whole track
    sigma = std(torsoFreq_Hz, 'omitnan');   % how spread out the values are
    isOutlier = abs(torsoFreq_Hz - mu) > 3*sigma;   % flag anything more than 3 standard deviations away as junk
    torsoFreq_refined = torsoFreq_Hz;   % start from the raw track
    torsoFreq_refined(isOutlier) = NaN;   % blank out the outliers
    torsoFreq_refined = fillmissing(torsoFreq_refined,'linear','EndValues','nearest');   % fill those gaps back in by interpolating between good neighbours
end

function results = RunGaitBandRefinement(x, fs, freq_low0, f_high_min, freq_high0, numIterations, CPI_s)
    freq_low  = freq_low0;   % start with the initial band before narrowing it down
    freq_high = freq_high0;
    results = struct('freq_low', {}, 'freq_high', {}, 'torsoFreq_Hz', {}, ...
        'torsoFreq_refined_Hz', {}, 'centreFreq_Hz', {});   % empty struct array, one entry per iteration
    for iter = 1:numIterations
        xFilt = IdealBandpassFilter(x, fs, freq_low, freq_high);   % bandpass filter the raw signal to the current band
        [S, freq_axis, time_axis, fsProc] = ComputeGaitSpectrogram_Torso(xFilt, fs, CPI_s, freq_high);   % spectrogram of the filtered signal
        Power_S = abs(S).^2;   % magnitude squared, power
        Ssmooth = SmoothSpectrogramNoise(Power_S, 21);   % smooth the power across time
        torsoFreq_Hz = TrackPeakFrequency(Ssmooth, freq_axis, freq_low, freq_high);   % pick out the torso frequency at each time frame
        torsoFreq_refined_Hz = RefineOutliers(torsoFreq_Hz);   % clean up outlier values in that track
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
            freq_low = max(0.5*centreFreq_Hz,freq_low0);   % narrow the band around the centre frequency just found
            freq_high = max(2.0*centreFreq_Hz,f_high_min);   % ready for the next pass
        end
    end
end

function [S, freq_axis, time_axis] = Zahraa_Spectrogram(x, win, noverlap, nfft, fs)
    x   = x(:).';
    win = win(:).';
    L = length(win);   % window length
    N = length(x);   % total signal length
    hop = L - noverlap;   % how far the window slides each step, based on the overlap
    M = floor((N - L) / hop) + 1;   % number of frames that fit in the signal
    S = zeros(nfft, M);   % will hold the fft result for every frame
    frame_starts = zeros(1, M);   % will hold where each frame starts, in samples
    for m = 1:M
        start_idx = (m - 1) * hop + 1;   % this frame's starting sample
        stop_idx  = start_idx + L - 1;   % this frame's ending sample
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;   % grab that chunk and taper its edges with the window
        S(:, m) = fft(segment, nfft).';   % fft of this windowed frame becomes one column of the spectrogram
    end
    time_axis = (frame_starts - 1) / fs;   % convert frame start samples into seconds
    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);   % real signal, keep the axis from 0 to fs
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end
