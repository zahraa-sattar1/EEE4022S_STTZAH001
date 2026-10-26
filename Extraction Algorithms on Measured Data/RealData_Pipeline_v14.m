%% Measured radar gait parameters
% Zahraa Sattar | EEE4022S
%
% Extract walking speed, step events, stride duration and stride length from a measured radar signal.
% Configuration settings select the detector rules, including the recording-edge prominence check.

function R = RealData_Pipeline_v14(RdrRxSignal, fs, fc, c, cfgList, opts)
if nargin == 1 && ischar(RdrRxSignal) && strcmp(RdrRxSignal, 'configs')
    R = makeConfigs();                                   % just the list of settings
    return;
end
if nargin < 6 || isempty(opts)
    opts = struct();
end
opts = setDefault(opts, 'CPI_torso_s', 0.2);            % torso window request, gives the 0.372 s window
opts = setDefault(opts, 'CPI_legs_s', 0.1);             % leg window request, gives the 0.186 s window
opts = setDefault(opts, 'Pth_dB', -45);                 % fixed leg threshold for the real data
opts = setDefault(opts, 'gapFactor', 1.5);              % old gap fill rule, gap must be this many times Tc/2
opts = setDefault(opts, 'TcNominal_s', 1.2);            % Tc used in pass 1 and as the fallback
opts = setDefault(opts, 'TcRange_s', [0.3 4]);          % a Tc estimate outside this is not believed
opts = setDefault(opts, 'minInterval_s', 0.35);         % steps closer than this are merged (Wang et al.)
opts = setDefault(opts, 'thresholdMultiplier', 2);      % no-motion threshold, same as the speed script
opts = setDefault(opts, 'numIterations', 3);            % band narrowing passes
opts = setDefault(opts, 'ReturnPlotData', false);       % true = keep the spectrograms for plotting
opts = setDefault(opts, 'frameRef_s', 1.0);             % length the body level is taken over for frameNorm, about one stride

cfgList = fillCfg(cfgList);                              % v7 settings get the v8 switches switched off
x = RdrRxSignal(:).';                                    % row vector, like the algorithms expect
recordingDuration_s = length(x) / fs;                    % length of the recording in seconds

%% Torso part, once for each band setting that is asked for
bandList = unique([cfgList.bandA]);                      % e.g. [0 0.25]
torso = cell(1, numel(bandList));
for b = 1:numel(bandList)
    torso{b} = runTorso(x, fs, fc, c, bandList(b), opts);   % walking speed, torso track, torso Tc
end

%% Leg part, once for each leg window that is asked for
winList = unique([cfgList.legWin_s]);
legs = cell(1, numel(winList));
for w = 1:numel(winList)
    legs{w} = prepareLegs(x, fs, fc, opts, winList(w));
end

%% One result for every setting
for k = 1:numel(cfgList)
    cfg = cfgList(k);
    T = torso{find(bandList == cfg.bandA, 1)};           % the torso run that matches this band
    L = legs{find(winList == cfg.legWin_s, 1)};          % the leg spectrogram that matches this window
    TcTwo = NaN;                                         % only filled in for the two-pass method
    minIntervalUsed = opts.minInterval_s;                % spacing used in the final pass (0.35 s unless propSpacing)
    fallback = 0;                                        % 1 if the Tc estimate could not be used

    if strcmp(cfg.tcMethod, 'twopass')
        rec1 = cfg.recover;
        if strcmp(rec1, 'fill')
            rec1 = 'none';                               % old gap fill needs Tc/2, and Tc is not known yet
        end
        S1 = runSteps(L, opts.TcNominal_s, cfg, rec1, opts);   % pass 1 with Tc = 1.2 s
        if numel(S1.times) >= 3
            TcTwo = 2 * median(diff(S1.times));          % a stride is two steps, median so a stray step does not matter
        end
        isOk = isfinite(TcTwo) && TcTwo > opts.TcRange_s(1) && TcTwo < opts.TcRange_s(2);
        if isOk
            Tc = TcTwo;
            stepOpts = opts;
            if cfg.propSpacing == 1 && cfg.gapAlpha > 0
                stepPeriod = Tc / 2;                              % the walk's own step period, from pass 1
                stepOpts.minInterval_s = max(2 * L.frame_period, cfg.gapAlpha * stepPeriod);   % same fraction of a step at every pace
            end
            minIntervalUsed = stepOpts.minInterval_s;
            S = runSteps(L, Tc, cfg, cfg.recover, stepOpts); % pass 2 with the walk's own Tc
        else
            Tc = opts.TcNominal_s;                       % too few steps or a silly value, keep 1.2 s
            fallback = 1;
            S = S1;                                      % pass 2 would be the same as pass 1
        end
    else
        Tc = T.TcTorso;                                  % old method, Tc from the torso bobs
        if ~isfinite(Tc)
            Tc = opts.TcNominal_s;                       % torso found fewer than two bobs, keep 1.2 s
            fallback = 1;
        end
        S = runSteps(L, Tc, cfg, cfg.recover, opts);
    end

    times = S.times;                                     % final step times, s
    if numel(times) >= 3
        strideDurations_s = times(3:end) - times(1:end-2);   % same foot to same foot
    else
        strideDurations_s = [];
    end

    vel = fillmissing(T.vel, 'linear', 'EndValues', 'nearest');   % no-motion frames are NaN, fill them in
    if cfg.align == 1
        tVel = T.time;                                   % both time axes already label the frame centre
        tStep = times;
    else
        tVel = T.time + (T.N / 2) / fs;                  % old: centre added a second time
        tStep = times + (L.N / 2) / fs;
    end
    strideLengths_m = strideLengthFromVelocity(vel, tVel, tStep);

    iv = diff(times);
    if numel(iv) >= 3
        suspect = sum(abs(iv - median(iv)) > 0.5 * median(iv));
    else
        suspect = 0;
    end
    suspectTc = sum(abs(iv - Tc/2) > 0.5 * (Tc/2));      % the old definition, against Tc/2

    R(k).name = cfg.name;
    R(k).estSpeed_mps = T.meanSpeed;
    R(k).recordingDuration_s = recordingDuration_s;
    R(k).Tc_used_s = Tc;
    R(k).Tc_torso_s = T.TcTorso;
    R(k).Tc_twopass_s = TcTwo;
    R(k).TcFallback = fallback;
    R(k).CPI_legs_actual_s = L.CPI_actual_s;
    R(k).minInterval_used_s = minIntervalUsed;           % base spacing of the final pass
    R(k).stepTimes_s = times;
    R(k).stepCount = numel(times);
    R(k).strideDurations_s = strideDurations_s;
    R(k).strideLengths_m = strideLengths_m;
    R(k).suspectIntervals = suspect;
    R(k).suspectIntervalsTc = suspectTc;
    if opts.ReturnPlotData
        Pmed = medianFor(L, S.medianFrames);             % the leg spectrogram this setting used
        P.legTime = L.time_axis;
        P.legFreq = L.freq_axis;
        P.legPdB = 10 * log10(Pmed ./ S.ref + eps);      % in dB, measured from the same level the threshold uses
        P.Pth_dB = S.Pth_dB;
        P.F_smooth = S.F_smooth;                         % smoothed leg frequency track
        P.stepIdx = S.stepIdx;
        P.stepTimes = times;
        P.torsoTime = T.time;
        P.torsoFreq = T.torsoFreq;                       % tracked torso frequency
        P.vel = T.vel;                                   % torso speed track
        P.lambda = c / fc;
        R(k).plot = P;
    else
        R(k).plot = [];
    end
end
end

function cfgList = makeConfigs()
o = struct('name', '', 'bandA', 0, 'tcMethod', 'torso', 'halfC', 0, 'recover', 'fill', 'centroid', 0, 'align', 0, ...
    'regap', 0, 'promGamma', 0, 'gapAlpha', 0, 'legWin_s', 0, 'Pth', NaN, 'frameNorm_dB', 0, 'localGamma', 0, 'localReach_s', 0, 'medianHz', 0, 'frameCap_dB', 0, 'propSpacing', 0, 'propEdge', 0);
f8 = o;
f8.bandA = 0.25;  f8.tcMethod = 'twopass';  f8.halfC = 1;  f8.recover = 'median2';  f8.centroid = 1;  f8.align = 1;  f8.regap = 1;
f8.promGamma = 0.1;  f8.gapAlpha = 0.6;                  % Selection.csv of SimVerify_20261012_214009
f12 = f8;
f12.frameNorm_dB = 10;  f12.localGamma = 0.25;  f12.localReach_s = 0.6;  f12.promGamma = 0;  f12.Pth = -40;   % Selection.csv of SimVerify3_20261013_084256
cfgList = o;           cfgList(1).name = 'Old';
cfgList(2) = f8;       cfgList(2).name = 'F8';
cfgList(3) = f12;      cfgList(3).name = 'F12';
cfgList(4) = f12;      cfgList(4).name = 'FN15';    cfgList(4).frameNorm_dB = 15;
cfgList(5) = f12;      cfgList(5).name = 'FN20';    cfgList(5).frameNorm_dB = 20;
cfgList(6) = f12;      cfgList(6).name = 'Cap15';   cfgList(6).frameCap_dB = 15;
cfgList(7) = f12;      cfgList(7).name = 'MF16';    cfgList(7).medianHz = 16;
end

function out = fillCfg(in)
def = struct('name', '', 'bandA', 0, 'tcMethod', 'torso', 'halfC', 0, 'recover', 'fill', 'centroid', 0, 'align', 0, ...
    'regap', 0, 'promGamma', 0, 'gapAlpha', 0, 'legWin_s', 0, 'Pth', NaN, 'frameNorm_dB', 0, 'localGamma', 0, 'localReach_s', 0, 'medianHz', 0, 'frameCap_dB', 0, 'propSpacing', 0, 'propEdge', 0);
fn = fieldnames(def);
out = repmat(def, 1, numel(in));
for k = 1:numel(in)
    for q = 1:numel(fn)
        if isfield(in(k), fn{q}) && ~isempty(in(k).(fn{q}))
            out(k).(fn{q}) = in(k).(fn{q});
        end
    end
end
end

function T = runTorso(x, fs, fc, c, bandA, opts)
fc_paper = 5.8e9;                                        % Yardibi's carrier
K = fc / fc_paper;                                       % scales the paper's frequencies to our radar
f_low0 = 5 * K;                                          % lowest frequency the torso band can start from
f_high_min = 50 * K;                                     % smallest upper edge for the old band
f_high0 = 4.0 * 2 * fc / c;                              % 4 m/s as a Doppler frequency, the starting upper edge

res = RunGaitBandRefinement(x, fs, f_low0, f_high_min, f_high0, opts.numIterations, opts.CPI_torso_s, bandA);
final = res(end);                                        % the most narrowed pass

idxBand = find(abs(final.freq_axis) >= final.freq_low & abs(final.freq_axis) <= final.freq_high);
energyPerFrame = sum(final.Ssmooth(idxBand, :), 1);      % energy inside the torso band, per frame
energyBaseline = median(energyPerFrame, 'omitnan');      % typical frame energy
energySpread = median(abs(energyPerFrame - energyBaseline), 'omitnan');   % how much it normally varies
energyThreshold = max(0, energyBaseline - opts.thresholdMultiplier * energySpread);
isNoMotion = energyPerFrame < energyThreshold;           % frames where nothing is moving
torsoFreq_sup = final.torsoFreq_refined_Hz;
torsoFreq_sup(isNoMotion) = 0;
vel = torsoFreq_sup * c / (2 * fc);                      % Doppler frequency to speed
vel(isNoMotion) = NaN;                                   % NaN so they do not pull the mean down
moving = torsoFreq_sup ~= 0;
T.meanSpeed = mean(vel(moving), 'omitnan');              % walking speed, mean over the moving frames
if all(torsoFreq_sup == 0)
    warning('Every frame was suppressed as no motion.');
end
T.vel = vel;
T.time = final.time_axis;                                % frame centres
T.torsoFreq = final.torsoFreq_refined_Hz;
T.N = 2^(nextpow2(round(opts.CPI_torso_s * fs)));        % torso window length in samples

f = final.torsoFreq_refined_Hz;
dt = final.time_axis(2) - final.time_axis(1);
prom = 0.25 * (max(f) - min(f));                         % a bob must stand out by a quarter of the swing
try
    [~, pk] = findpeaks(f, 'MinPeakDistance', max(1, round(0.3 / dt)), 'MinPeakProminence', prom);
catch
    pk = [];
end
iv = diff(final.time_axis(pk));
if isempty(iv)
    T.TcTorso = NaN;
else
    T.TcTorso = 2 * mean(iv);                            % two steps per stride
end
end

function results = RunGaitBandRefinement(x, fs, freq_low0, f_high_min, freq_high0, numIterations, CPI_s, bandA)
freq_low = freq_low0;                                    % start with the wide band
freq_high = freq_high0;
results = struct('freq_low', {}, 'freq_high', {}, 'torsoFreq_refined_Hz', {}, 'centreFreq_Hz', {}, ...
    'Ssmooth', {}, 'freq_axis', {}, 'time_axis', {});
for iter = 1:numIterations
    xFilt = IdealBandpassFilter(x, fs, freq_low, freq_high);   % keep only the current band
    [S, freq_axis, time_axis] = ComputeGaitSpectrogram_Torso(xFilt, fs, CPI_s, freq_high);
    Power_S = abs(S).^2;                                 % power
    Ssmooth = movmean(abs(Power_S), 21, 1);              % smooth across the frequency bins
    validIdx = find(abs(freq_axis) >= freq_low & abs(freq_axis) <= freq_high);
    [~, peakIdx] = max(Ssmooth(validIdx, :), [], 1);     % strongest bin in the band at each frame
    torsoFreq_Hz = freq_axis(validIdx(peakIdx));         % that bin is the torso frequency
    mu = mean(torsoFreq_Hz, 'omitnan');
    sigma = std(torsoFreq_Hz, 'omitnan');
    isOutlier = abs(torsoFreq_Hz - mu) > 3 * sigma;      % more than 3 std away is junk
    torsoFreq_refined = torsoFreq_Hz;
    torsoFreq_refined(isOutlier) = NaN;
    torsoFreq_refined = fillmissing(torsoFreq_refined, 'linear', 'EndValues', 'nearest');
    centreFreq_Hz = mean(torsoFreq_refined, 'omitnan');  % average torso frequency, centres the next band
    results(iter).freq_low = freq_low;
    results(iter).freq_high = freq_high;
    results(iter).torsoFreq_refined_Hz = torsoFreq_refined;
    results(iter).centreFreq_Hz = centreFreq_Hz;
    results(iter).Ssmooth = Ssmooth;
    results(iter).freq_axis = freq_axis;
    results(iter).time_axis = time_axis;
    if iter < numIterations
        if bandA > 0
            freq_low = max((1 - bandA) * centreFreq_Hz, freq_low0);   % new: 0.75x to 1.25x
            freq_high = (1 + bandA) * centreFreq_Hz;
        else
            freq_low = max(0.5 * centreFreq_Hz, freq_low0);           % old: 0.5x to 2x
            freq_high = max(2.0 * centreFreq_Hz, f_high_min);
        end
    end
end
end

function [S, freq_axis, time_axis] = ComputeGaitSpectrogram_Torso(x, fs, CPI_s, freq_high)
N = 2^(nextpow2(round(CPI_s * fs)));                     % window length, power of two
win = hamming(N);
noverlap = N * 3 / 4;                                    % 75% overlap
[S_temp, freq_axis_temp, time_axis] = SpectrogramCentre(x, win, noverlap, N, fs);
keep = find(abs(freq_axis_temp) <= freq_high * 3);       % only the low frequencies, the torso lives there
S = S_temp(keep, :);
freq_axis = freq_axis_temp(keep);
end

function [y, f, mask] = IdealBandpassFilter(x, fs, f_low, f_high)
isRowInput = isrow(x);
x = x(:).';
N = length(x);
X = fft(x);
Xs = fftshift(X);
f = (-floor(N/2):ceil(N/2)-1) * (fs / N);
mask = (abs(f) >= f_low) & (abs(f) <= f_high);           % keep the band, zero the rest
Xu = ifftshift(Xs .* mask);
y = ifft(Xu);
if isreal(x)
    y = real(y);                                         % a real recording stays real (MATLAB already does this, Octave does not)
end
if ~isRowInput
    y = y.';
end
end

function L = prepareLegs(x, fs, fc, opts, legWin_s)
L.freq_high = 350 * (fc / 1e9) / 5.8;                    % leg Doppler limit (1448 Hz at 24 GHz)
avgSamples = floor(fs / (1.5 * L.freq_high));            % averaging length, keeps the cut-off above the leg band
xAvg = movmean(x, avgSamples);                           % smooths the raw signal before the spectrogram
if legWin_s > 0
    N = round(legWin_s * fs);                            % exact window length asked for
    nfft = 2^(nextpow2(N));                              % FFT zero-padded to a power of two
else
    N = 2^(nextpow2(round(opts.CPI_legs_s * fs)));       % leg window length, power of two (v7)
    nfft = N;
end
win = hamming(N);
[S, freq_axis_full, time_axis] = SpectrogramCentre(xAvg, win, round(N * 3 / 4), nfft, fs);   % 75% overlap
L.N = N;
L.CPI_actual_s = N / fs;                                 % the window actually used, 0.186 s
keep = find(freq_axis_full > 0 & freq_axis_full <= L.freq_high);
L.freq_axis = freq_axis_full(keep);
L.time_axis = time_axis;                                 % frame centres
L.frame_period = time_axis(2) - time_axis(1);
L.Pm = SubtractBackgroundNoise(S(keep, :));              % noise floor taken off
L.cache = containers.Map('KeyType', 'char', 'ValueType', 'any');   % median filter results, one per window length
end

function Pmedian = medianFor(L, frames)
key = sprintf('%d', frames);
if isKey(L.cache, key)
    Pmedian = L.cache(key);                              % already done for this window length
else
    Pmedian = movmedian(L.Pm, frames, 1);                % median filter across frequency
    L.cache(key) = Pmedian;
end
end

function S = runSteps(L, Tc, cfg, recover, opts)
Pth = cfg.Pth;
if isnan(Pth)
    Pth = opts.Pth_dB;                                        % the fixed threshold unless this setting has its own
end
if cfg.medianHz > 0
    medianFrames = max(1, round(cfg.medianHz / (L.freq_axis(2) - L.freq_axis(1))));   % width set in Hz, as a number of frequency bins
else
    medianFrames = max(1, round((Tc / 8) / L.frame_period));      % v9: Tc/8 counted in frames, but used along frequency
end
Pmedian = medianFor(L, medianFrames);
ref = referenceLevel(Pmedian, L.freq_axis, L.freq_high, Pth, cfg.frameNorm_dB, max(1, round(opts.frameRef_s / L.frame_period)), cfg.frameCap_dB);   % level the threshold is measured from
[legFreq, isDetected] = SelectLegFrequency(Pmedian, L.freq_axis, Pth, ref);   % highest bin above the threshold
legRefined = RefineLegFrequencyOutliers(legFreq, isDetected);
avgFrames = max(1, round((Tc / 8) / L.frame_period));         % smoothing window is also Tc/8
F_smooth = movmean(legRefined, avgFrames, 'omitnan');
normLeg = F_smooth / max(F_smooth, [], 'omitnan');            % scaled so the biggest value is 1

if cfg.halfC == 1
    cPeak = 0.5 * (1 / L.CPI_actual_s) / L.freq_high;         % half of c, from the actual window
else
    cPeak = (1 / opts.CPI_legs_s) / L.freq_high;              % old: from the requested 0.1 s
end
stepIdx = find((normLeg(2:end-1) >= normLeg(1:end-2)) & (normLeg(2:end-1) >= normLeg(3:end) + cPeak)) + 1;   % Wang's peak rule

if cfg.promGamma > 0
    Fok = F_smooth(isfinite(F_smooth));
    q = quantileSimple(Fok, [0.1 0.9]);
    span = q(2) - q(1);                                        % how far the leg track normally swings
    prom = peakProminence(F_smooth, stepIdx);                  % how far each candidate stands out
    stepIdx = stepIdx(prom >= cfg.promGamma * span);           % a small bump is not a step
end
if cfg.localGamma > 0 && numel(stepIdx) >= 2
    if cfg.propEdge == 1
        prom = peakProminenceEdge(F_smooth, stepIdx);          % a side cut off by the recording edge is not counted
    else
        prom = peakProminence(F_smooth, stepIdx);              % how far each candidate stands out
    end
    tS = L.time_axis(stepIdx);
    reach = max(0.75 * Tc / 2, cfg.localReach_s);              % three quarters of a step interval, at least localReach_s
    keepL = true(size(stepIdx));
    for k = 1:numel(stepIdx)
        nb = abs(tS - tS(k)) < reach;                          % candidates close to this one
        nb(k) = false;
        if any(nb)
            keepL(k) = prom(k) >= cfg.localGamma * max(prom(nb));   % a small bump next to a big swing is not a step
        end
    end
    stepIdx = stepIdx(keepL);
end
stepIdx = RemoveCloseSteps(stepIdx, L.time_axis, F_smooth, opts.minInterval_s);   % Wang's 0.35 s rule

minGap = opts.minInterval_s;
if cfg.gapAlpha > 0 && numel(stepIdx) >= 4
    m = median(diff(L.time_axis(stepIdx)));                    % the walk's own normal step gap
    minGap = max(opts.minInterval_s, cfg.gapAlpha * m);        % slow walkers get a longer minimum gap
    stepIdx = RemoveCloseSteps(stepIdx, L.time_axis, F_smooth, minGap);
end

switch recover
    case 'fill'
        stepIdx = FillMissedSteps_v1(stepIdx, L.time_axis, F_smooth, Tc / 2, opts.gapFactor);   % 29 Sep rule
    case 'median'
        stepIdx = RecoverMissedSteps(stepIdx, L.time_axis, F_smooth, opts.minInterval_s);       % v7 rule
    case 'median2'
        stepIdx = RecoverMissedSteps2(stepIdx, L.time_axis, F_smooth);                         % v8 rule
    otherwise
end

if cfg.centroid == 1
    times = RefineStepTimes(L.time_axis, F_smooth, stepIdx, Tc);   % between the frames
else
    times = L.time_axis(stepIdx);                                  % on the frame grid
end
times = times(:).';
if cfg.regap == 1 && numel(times) >= 2
    keep = mergeClose(times, F_smooth(stepIdx), minGap);           % refined times closer than the minimum gap are one step
    stepIdx = stepIdx(keep);
    times = times(keep);
end
S.stepIdx = stepIdx;
S.times = times;
S.F_smooth = F_smooth;
S.medianFrames = medianFrames;
S.Pth_dB = Pth;
S.ref = ref;
end

function ref = referenceLevel(Pmedian, freq_axis, freq_high, Pth, margin_dB, nRef, cap_dB)
gmax = max(Pmedian(:)) + eps;                            % strongest bin of the whole walk (v8)
if margin_dB <= 0
    ref = gmax * ones(1, size(Pmedian, 2));              % v8: one level for the whole walk
    return;
end
frameMax = max(Pmedian, [], 1);                          % strongest bin of each frame, the body return
ref = movmax(frameMax, nRef);                            % highest over about one stride, so one weak frame cannot set it
top = freq_axis >= 0.85 * freq_high;                     % top 15% of the leg band, above about 7.7 m/s
noiseDb = quantileSimple(reshape(10 * log10(Pmedian(top, :) / gmax + eps), 1, []), 0.95);   % noise peaks of this walk (95th percentile), dB below the strongest bin
refMin = gmax * 10^((noiseDb + margin_dB - Pth) / 10);   % lowest reference that keeps the threshold margin_dB above the noise
ref = min(gmax, max(ref, refMin));
if cap_dB > 0
    ref = max(ref, gmax * 10^(-cap_dB / 10));            % never more than cap_dB below the whole-walk level
end
end

function Pm = SubtractBackgroundNoise(S)
P = abs(S).^2;
Pback = zeros(size(P, 1), 1);
nBack = max(1, round(0.20 * size(P, 2)));                % quietest 20% of the frames
for i = 1:size(P, 1)
    sortedPower = sort(P(i, :), 'ascend');
    Pback(i) = mean(sortedPower(1:nBack));               % their average is the noise floor of this bin
end
Pm = P - repmat(Pback, 1, size(P, 2));
Pm(Pm < 0) = eps;                                        % power cannot be negative
end

function [legFreq_Hz, isDetected] = SelectLegFrequency(Pmedian, freq_axis, Pth, ref)
PmedianNorm = Pmedian ./ repmat(ref, size(Pmedian, 1), 1);   % reference level becomes 1 (v8: biggest value of the walk)
PmedianDb = 10 * log10(PmedianNorm + eps);
aboveThresh = PmedianDb > Pth;
legFreq_Hz = NaN(1, size(Pmedian, 2));
isDetected = false(1, size(Pmedian, 2));
for t = 1:size(Pmedian, 2)
    idx = find(aboveThresh(:, t), 1, 'last');            % highest frequency above the threshold
    if ~isempty(idx)
        legFreq_Hz(t) = freq_axis(idx);
        isDetected(t) = true;
    end
end
end

function legFreq_interp = RefineLegFrequencyOutliers(legFreq_Hz, isDetected)
mu = mean(legFreq_Hz(isDetected));
sigma = std(legFreq_Hz(isDetected));
isOutlier = false(size(legFreq_Hz));
isOutlier(isDetected) = abs(legFreq_Hz(isDetected) - mu) > 3 * sigma;   % 3 std rule
legFreq_refined = legFreq_Hz;
legFreq_refined(isOutlier) = NaN;
legFreq_interp = fillmissing(legFreq_refined, 'linear', 'EndValues', 'nearest');
end

function stepIdx_refined = RemoveCloseSteps(stepIdx, time_axis, F_smooth, minInterval_s)
keep = mergeClose(time_axis(stepIdx), F_smooth(stepIdx), minInterval_s);
stepIdx_refined = stepIdx(keep);
end

function keep = mergeClose(tms, strength, minGap)
keep = true(size(tms));
lastKept = 1;
for k = 2:numel(tms)
    if (tms(k) - tms(lastKept)) <= minGap                % too close to the last kept step
        if strength(k) > strength(lastKept)
            keep(lastKept) = false;                      % this one is stronger
            lastKept = k;
        else
            keep(k) = false;                             % the earlier one is stronger
        end
    else
        lastKept = k;
    end
end
end

function prom = peakProminence(F, idx)
n = numel(F);
prom = zeros(size(idx));
for k = 1:numel(idx)
    i = idx(k);
    pk = F(i);
    l = i;
    leftMin = pk;
    while l > 1 && ~(F(l-1) > pk)
        l = l - 1;
        leftMin = min(leftMin, F(l));                    % lowest point on the way left
    end
    r = i;
    rightMin = pk;
    while r < n && ~(F(r+1) > pk)
        r = r + 1;
        rightMin = min(rightMin, F(r));                  % lowest point on the way right
    end
    prom(k) = pk - max(leftMin, rightMin);
end
end

function prom = peakProminenceEdge(F, idx)
n = numel(F);
prom = zeros(size(idx));
for k = 1:numel(idx)
    i = idx(k);
    pk = F(i);
    l = i;
    leftMin = pk;
    while l > 1 && ~(F(l-1) > pk)
        l = l - 1;
        leftMin = min(leftMin, F(l));                    % lowest point on the way left
    end
    r = i;
    rightMin = pk;
    while r < n && ~(F(r+1) > pk)
        r = r + 1;
        rightMin = min(rightMin, F(r));                  % lowest point on the way right
    end
    leftCut = (l == 1) && (i > 1) && (F(1) <= leftMin);  % the lowest point on the left is the first frame
    rightCut = (r == n) && (i < n) && (F(n) <= rightMin);% the lowest point on the right is the last frame
    if leftCut && ~rightCut
        prom(k) = pk - rightMin;
    elseif rightCut && ~leftCut
        prom(k) = pk - leftMin;
    else
        prom(k) = pk - max(leftMin, rightMin);
    end
end
end

function stepIdx_out = RecoverMissedSteps2(stepIdx, time_axis, F)
stepIdx_out = stepIdx;
for added = 1:5                                          % at most five steps put back per walk
    t = time_axis(stepIdx_out);
    iv = diff(t);
    if numel(iv) < 3
        break;                                           % too few steps to know the normal gap
    end
    m = median(iv);
    gaps = find(iv > 1.6 * m);
    putBack = false;
    for g = gaps(:).'
        a = stepIdx_out(g);
        b = stepIdx_out(g + 1);
        qs = a+1:b-1;
        qs = qs((time_axis(qs) - time_axis(a)) > 0.5 * m & (time_axis(b) - time_axis(qs)) > 0.5 * m);
        if isempty(qs)
            continue;                                    % gap too short to hold a step, look at the next gap
        end
        tops = qs(F(qs) >= F(qs - 1) & F(qs) >= F(qs + 1));
        if isempty(tops)
            tops = qs;                                   % no clear top: the highest point in the gap
        end
        [~, mx] = max(F(tops));
        stepIdx_out = sort([stepIdx_out, tops(mx)]);
        putBack = true;
        break;                                           % gaps change after a step is added, so start again
    end
    if ~putBack
        break;
    end
end
end

function q = quantileSimple(v, p)
v = sort(v(:));
n = numel(v);
q = NaN(size(p));
if n == 0
    return;
end
for k = 1:numel(p)
    pos = 1 + (n - 1) * p(k);
    lo = floor(pos);
    hi = min(n, lo + 1);
    q(k) = v(lo) + (pos - lo) * (v(hi) - v(lo));
end
end

function stepIdx_out = RecoverMissedSteps(stepIdx, time_axis, F_smooth, minInterval_s)
stepIdx_out = stepIdx;
for pass = 1:3                                           % up to three steps can be put back
    t = time_axis(stepIdx_out);
    iv = diff(t);
    if numel(iv) < 3
        break;                                           % too few steps to know the normal gap
    end
    g = find(iv > 1.6 * median(iv), 1);                  % first gap far longer than normal
    if isempty(g)
        break;
    end
    a = stepIdx_out(g);                                  % step before the gap
    b = stepIdx_out(g + 1);                              % step after the gap
    cand = [];
    for q = a+1:b-1
        isTop = F_smooth(q) >= F_smooth(q-1) && F_smooth(q) >= F_smooth(q+1);
        farEnough = (time_axis(q) - time_axis(a)) > minInterval_s && (time_axis(b) - time_axis(q)) > minInterval_s;
        if isTop && farEnough
            cand(end+1) = q;                             %#ok<AGROW>
        end
    end
    if isempty(cand)
        break;
    end
    [~, m] = max(F_smooth(cand));                        % keep the highest top
    stepIdx_out = sort([stepIdx_out, cand(m)]);
end
end

function t = RefineStepTimes(time_axis, F, stepIdx, Tc)
t = time_axis(stepIdx);
if isempty(stepIdx) || numel(time_axis) < 2
    return;
end
dt = time_axis(2) - time_axis(1);
N = numel(F);
wAll = round(Tc / 4 / dt);                               % half a step interval, in frames
for k = 1:numel(stepIdx)
    i = stepIdx(k);
    wBack = wAll;
    wFwd = wAll;
    if k > 1
        wBack = min(wBack, floor((i - stepIdx(k-1)) / 2));      % never past the middle to the previous step
    end
    if k < numel(stepIdx)
        wFwd = min(wFwd, floor((stepIdx(k+1) - i) / 2));        % never past the middle to the next step
    end
    a = max(1, i - wBack);
    b = min(N, i + wFwd);
    seg = F(a:b);
    pk = max(seg);                                       % top of the peak
    base = min(F(max(1, i - wAll):min(N, i + wAll)));    % lowest leg frequency nearby
    wgt = max(seg - (pk - 0.5 * (pk - base)), 0);        % only the top half counts
    wgt(~isfinite(wgt)) = 0;
    if sum(wgt) > 0
        t(k) = sum(time_axis(a:b) .* wgt) / sum(wgt);    % weighted centre
    end
end
end

function strideLengths_m = strideLengthFromVelocity(velocity_m_s, velocityTime_s, stepTimes_s)
displacement_m = cumtrapz(velocityTime_s, velocity_m_s);                       % distance walked so far
displacement_atSteps = interp1(velocityTime_s, displacement_m, stepTimes_s, 'linear', 'extrap');   % distance at each step
strideLengths_m = abs(displacement_atSteps(3:end) - displacement_atSteps(1:end-2));   % two steps apart = one stride
end

function [S, freq_axis, t_centre] = SpectrogramCentre(x, win, noverlap, nfft, fs)
[S, freq_axis, ~] = Zahraa_Spectrogram(x, win, noverlap, nfft, fs);
Lw = numel(win);
hop = Lw - noverlap;
M = size(S, 2);
t_centre = ((0:M-1) * hop + Lw / 2) / fs;
end

function s = setDefault(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
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

function stepIdx = FillMissedSteps_v1(stepIdx, time_axis, F_smooth, expectedStep_s, gapFactor)

    if numel(stepIdx) < 2 || ~isfinite(expectedStep_s) || expectedStep_s <= 0
        return;   % nothing to fill between
    end
    frameStep_s  = time_axis(2) - time_axis(1);                           % time between spectrogram frames
    searchFrames = max(1, round((expectedStep_s/4) / frameStep_s));       % look a quarter of a step either side of each split for the leg peak
    newIdx = stepIdx(1);
    for k = 2:numel(stepIdx)
        gap_s    = time_axis(stepIdx(k)) - time_axis(stepIdx(k-1));      % time between these two detected steps
        numParts = round(gap_s / expectedStep_s);                         % how many steps this gap should really hold
        if gap_s > gapFactor * expectedStep_s && numParts >= 2
            for j = 1:numParts-1
                target_s = time_axis(stepIdx(k-1)) + j * gap_s / numParts;   % evenly spaced guess for the missing step
                [~, centre] = min(abs(time_axis - target_s));                % nearest frame to that guess
                lo = max(stepIdx(k-1) + 1, centre - searchFrames);           % search window, kept inside the gap
                hi = min(stepIdx(k) - 1, centre + searchFrames);
                if lo > hi
                    newIdx(end+1) = centre;                                  % window too small, use the guess itself
                else
                    [~, m] = max(F_smooth(lo:hi));                           % highest leg frequency near the guess, where the step most likely is
                    newIdx(end+1) = lo + m - 1;
                end
            end
        end
        newIdx(end+1) = stepIdx(k);   % keep the detected step
    end
    stepIdx = unique(newIdx);   % in time order, no repeats
end
