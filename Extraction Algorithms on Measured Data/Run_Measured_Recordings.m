%% Run measured recordings
% Zahraa Sattar | EEE4022S
%
% Process the ten locally held recordings using the retained Final settings.
% Save the four parameter estimates and a detection plot for each walk.

function outDir = Run_Measured_Recordings(userCfg)
root = fileparts(fileparts(mfilename('fullpath')));
cfg.ReferenceFile = fullfile(root, 'Participant Recordings', 'Demo_References.csv');
cfg.ParameterFile = fullfile(fileparts(mfilename('fullpath')), 'Measured_Final_Parameters.csv');
cfg.OutputFolder = fullfile(root, 'Results', ['Measured_' datestr(now, 'yyyymmdd_HHMMSSFFF')]);
cfg.SavePlots = true;
if nargin > 0
    fields = fieldnames(userCfg);
    for k = 1:numel(fields)
        assert(isfield(cfg, fields{k}), 'Unknown option: %s', fields{k});
        cfg.(fields{k}) = userCfg.(fields{k});
    end
end
G = readtable(cfg.ReferenceFile, 'TextType', 'string');
[settings, opts, fc, c] = readSettings(cfg.ParameterFile);
opts.ReturnPlotData = cfg.SavePlots;
outDir = cfg.OutputFolder;
mkdir(outDir);
copyfile(cfg.ReferenceFile, fullfile(outDir, 'References_Used.csv'));
copyfile(cfg.ParameterFile, fullfile(outDir, 'Parameters_Used.csv'));
rows = struct([]);
for k = 1:height(G)
    file = fullfile(root, 'Participant Recordings', G.recording_file(k));
    [samples, fs] = audioread(file, 'native');
    assert(size(samples, 1) == G.num_samples(k) && fs == G.sample_rate_Hz(k), ...
        'WAV samples or sample rate differ from the reference: %s', G.wav_file(k));
    duration = size(samples, 1) / fs;
    assert(abs(duration - G.duration_s(k)) < 1e-10, 'WAV duration differs from the reference.');
    assert(abs(G.walkedDistance_m(k)/duration - G.gtSpeed_mps(k)) < 1e-10, ...
        'Speed reference is inconsistent for %s.', G.wav_file(k));
    assert(abs(2*duration/G.gtStepCount(k) - G.gtStrideDuration_s(k)) < 1e-10 && ...
        abs(2*G.walkedDistance_m(k)/G.gtStepCount(k) - G.gtStrideLength_m(k)) < 1e-10, ...
        'Stride references are inconsistent for %s.', G.wav_file(k));
    x = -double(samples(:, 1)).';
    result = RealData_Pipeline_v14(x, fs, fc, c, settings, opts);
    rows(k).wav_file = G.wav_file(k);
    rows(k).participant = G.participant(k);
    rows(k).gtSpeed_mps = G.gtSpeed_mps(k);
    rows(k).estSpeed_mps = result.estSpeed_mps;
    rows(k).gtStepCount = G.gtStepCount(k);
    rows(k).estStepCount = result.stepCount;
    rows(k).gtStrideDuration_s = G.gtStrideDuration_s(k);
    rows(k).estStrideDuration_mean_s = mean(result.strideDurations_s, 'omitnan');
    rows(k).gtStrideLength_m = G.gtStrideLength_m(k);
    rows(k).estStrideLength_mean_m = mean(result.strideLengths_m, 'omitnan');
    rows(k).stepTimes_s = string(strtrim(sprintf('%.8f ', result.stepTimes_s)));
    rows(k).Tc_used_s = result.Tc_used_s;
    rows(k).TcFallback = result.TcFallback;
    if cfg.SavePlots
        saveWalkPlot(outDir, G.wav_file(k), result, G.gtSpeed_mps(k), G.gtStepCount(k));
    end
    fprintf('[%d/%d] %s: %d steps, %.4f m/s\n', k, height(G), G.wav_file(k), result.stepCount, result.estSpeed_mps);
end
writetable(struct2table(rows), fullfile(outDir, 'Results.csv'));
fprintf('Results saved in %s\n', outDir);
end

function [setting, opts, fc, c] = readSettings(file)
importOpts = detectImportOptions(file);
importOpts = setvartype(importOpts, importOpts.VariableNames, 'string');
T = readtable(file, importOpts);
configs = RealData_Pipeline_v14('configs');
setting = configs(strcmp({configs.name}, 'F12'));
assert(isscalar(setting), 'Expected one F12 base configuration.');
opts = struct();
fc = NaN;
c = NaN;
for k = 1:height(T)
    key = char(T.key(k));
    text = char(T.value(k));
    value = str2double(text);
    if startsWith(text, '[')
        value = sscanf(strrep(strrep(text, '[', ''), ']', ''), '%f').';
    elseif isnan(value) && ~strcmpi(text, 'NaN')
        value = text;
    end
    if isfield(setting, key)
        setting.(key) = value;
    elseif strcmp(key, 'fc_Hz')
        fc = value;
    elseif strcmp(key, 'c_mps')
        c = value;
    else
        opts.(key) = value;
    end
end
assert(isfinite(fc) && isfinite(c), 'Carrier frequency and wave speed are required.');
end

function saveWalkPlot(outDir, wavFile, result, referenceSpeed, referenceCount)
P = result.plot;
plotDir = fullfile(outDir, 'WalkPlots');
if ~isfolder(plotDir), mkdir(plotDir); end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 700]);
figCleanup = onCleanup(@() close(fig));
subplot(2, 1, 1);
imagesc(P.legTime, P.legFreq*P.lambda/2, P.legPdB);
axis xy;
caxis([P.Pth_dB-10, 0]);
ylim([0, P.legFreq(end)*P.lambda/2]);
colorbar;
hold on;
plot(P.legTime, P.F_smooth*P.lambda/2, 'w-', 'LineWidth', 1.2);
plot(P.stepTimes, interp1(P.legTime, P.F_smooth, P.stepTimes, 'linear', 'extrap')*P.lambda/2, ...
    'ro', 'MarkerFaceColor', 'r');
xlabel('Time (s)'); ylabel('Leg speed (m/s)');
title(sprintf('%s: %d detected steps, %d reference steps', strrep(wavFile, '_', ' '), result.stepCount, referenceCount));
subplot(2, 1, 2);
plot(P.torsoTime, P.vel, 'b-', 'LineWidth', 1.2);
hold on;
yline(referenceSpeed, 'r--', 'Reference');
xlabel('Time (s)'); ylabel('Walking speed (m/s)');
[~, stem] = fileparts(wavFile);
exportgraphics(fig, fullfile(plotDir, stem+'.png'), 'Resolution', 150);
end
