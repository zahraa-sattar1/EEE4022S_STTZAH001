%% Measured-data split check
% Zahraa Sattar | EEE4022S
%
% Repeat the report's 75/25 stability and calibration checks from saved estimates.
% These walks informed development, so this is not an unseen validation set.

function RealData_SplitCheck_v1(userCfg)

cfg.ResultsFolder = fullfile(fileparts(mfilename('fullpath')), 'Report Results');                                  % empty = reviewed current-data pointer
cfg.GroundTruthFile = fullfile(fileparts(mfilename('fullpath')), 'Report Results', 'Report_References.csv');                                % empty = current pointer reference file
cfg.nRandomSplits = 500;                                 % random 75/25 splits for the stability check
cfg.Seed = 2026;                                         % fixed seed, chosen before looking at any result
cfg.PaceEdges_mps = [0.8 1.0 1.2];                       % same pace regions as the summary
cfg.PassAccuracy_pct = 90;
if nargin > 0 && isstruct(userCfg)
    fn = fieldnames(userCfg);
    for q = 1:numel(fn)
        cfg.(fn{q}) = userCfg.(fn{q});
    end
end
[cfg.ResultsFolder, cfg.GroundTruthFile] = RealData_ResolveReportInputs_v1(cfg.ResultsFolder, cfg.GroundTruthFile);
T = readCsvStruct(fullfile(cfg.ResultsFolder, 'Results_AllConfigs.csv'));
T = applyGroundTruth(T, cfg.GroundTruthFile);
hasError = ~cellfun(@isempty, T.errorMsg);
singlePass = strcmpi(strtrim(T.walkType), 'Single pass');
configs = uniqueStable(T.config);
fin = configs{end};                                       % the Final setting is the last one the runner wrote
[~, gtName] = fileparts(cfg.GroundTruthFile);
outDir = fullfile(cfg.ResultsFolder, ['SplitCheck_v1_' strrep(gtName, 'RealData_', '')]);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

params = {'Walking speed', 'estSpeed_mps', 'gtSpeed_video_mps'; ...
          'Step count', 'estStepCount', 'gtStepCount'; ...
          'Stride duration', 'estStrideDuration_mean_s', 'gtStrideDuration_s'; ...
          'Stride length', 'estStrideLength_mean_m', 'gtStrideLength_m'};
nP = size(params, 1);

%% The walks that are split (one row per walk, taken from the Final rows)
rows = find(strcmp(T.config, fin) & ~hasError & singlePass);
wav = T.wav_file(rows);
part = T.participant(rows);
pace = 1 + (T.gtSpeed_video_mps(rows) >= cfg.PaceEdges_mps(1)) + (T.gtSpeed_video_mps(rows) >= cfg.PaceEdges_mps(2)) + (T.gtSpeed_video_mps(rows) >= cfg.PaceEdges_mps(3));
trialNo = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), T.trial_label(rows));
paceNames = {'Slow', 'Comfortable', 'Fast', 'Very fast'};
n = numel(rows);
parts = sort(unique(part));

%% 1. The fixed split: every 4th walk of each participant x pace group is held out
isTest = false(n, 1);
for q = 1:numel(parts)
    for b = 1:4
        idx = find(strcmp(part, parts{q}) & pace == b);
        [~, o] = sort(trialNo(idx));
        idx = idx(o);                                     % ordered by trial number
        nt = floor(0.25 * numel(idx) + 0.5);              % a quarter of the group, rounded
        if nt > 0
            isTest(idx(round(((1:nt) - 0.5) * numel(idx) / nt))) = true;   % evenly spaced through the trials
        end
    end
end
fid = fopen(fullfile(outDir, 'SplitCheck_Membership.csv'), 'w');
fprintf(fid, 'wav_file,participant,pace,set\n');
for i = 1:n
    fprintf(fid, '%s,%s,%s,%s\n', wav{i}, part{i}, paceNames{pace(i)}, ternary(isTest(i), 'held-out 25%', 'tuning 75%'));
end
fclose(fid);

R = {};
R{end+1} = sprintf('SPLIT CHECK v1   folder %s   made %s', outDir, datestr(now, 'yyyy-mm-dd HH:MM'));
R{end+1} = sprintf('Ground truth file: %s. Single-pass walks only (%d walks): %d in the tuning set (75%%), %d held out (25%%).', cfg.GroundTruthFile, n, sum(~isTest), sum(isTest));
R{end+1} = 'Split rule: a quarter of each participant x pace group (rounded), spread evenly through the trial order, is held out. No random numbers are used.';
R{end+1} = '  group count (tuning / held-out)   ' ;
line = sprintf('  %-6s', '');
for b = 1:4
    line = [line sprintf(' %12s', paceNames{b})]; %#ok<AGROW>
end
R{end+1} = line;
for q = 1:numel(parts)
    line = sprintf('  %-6s', parts{q});
    for b = 1:4
        m = strcmp(part, parts{q}) & pace == b;
        line = [line sprintf(' %12s', sprintf('%d / %d', sum(m & ~isTest), sum(m & isTest)))]; %#ok<AGROW>
    end
    R{end+1} = line; %#ok<AGROW>
end

%% 2. Every setting on the tuning set and on the held-out set
R{end+1} = '';
R{end+1} = 'ACCURACY (100 - RMSE%) OF EVERY SETTING ON THE TUNING WALKS AND ON THE HELD-OUT WALKS';
R{end+1} = sprintf('  %-8s | %-45s | %-45s | %s', 'setting', 'tuning 75%: speed steps strDur strLen', 'held-out 25%: speed steps strDur strLen', 'tuning mean');
tuneAcc = zeros(numel(configs), nP);
testAcc = zeros(numel(configs), nP);
for k = 1:numel(configs)
    [eT, gT, eH, gH] = deal([]);
    for p = 1:nP
        [estV, gtV] = valuesOfConfig(T, configs{k}, wav, params(p, :));
        tuneAcc(k, p) = accuracyOf(estV(~isTest), gtV(~isTest));
        testAcc(k, p) = accuracyOf(estV(isTest), gtV(isTest));
    end
    R{end+1} = sprintf('  %-8s | %10.2f %9.2f %9.2f %9.2f      | %10.2f %9.2f %9.2f %9.2f      | %.2f', configs{k}, tuneAcc(k, :), testAcc(k, :), mean(tuneAcc(k, :))); %#ok<AGROW>
end
[~, kb] = max(mean(tuneAcc, 2));
kf = find(strcmp(configs, fin), 1);
R{end+1} = sprintf('  Best setting on the tuning walks: %s. Its held-out accuracy: %.2f %.2f %.2f %.2f. Final (chosen on the simulator): %.2f %.2f %.2f %.2f.', ...
    configs{kb}, testAcc(kb, :), testAcc(kf, :));
R{end+1} = sprintf('  Final, tuning against held-out (accuracy points, held-out minus tuning): %+.2f %+.2f %+.2f %+.2f. A small difference means the result does not depend on which walks were scored.', testAcc(kf, :) - tuneAcc(kf, :));

%% 3. One speed scale fitted on the tuning walks and tested on the held-out walks
[estS, gtS] = valuesOfConfig(T, fin, wav, params(1, :));
[estL, gtL] = valuesOfConfig(T, fin, wav, params(4, :));
kScale = 1 / mean(estS(~isTest) ./ gtS(~isTest));       % video speed / radar speed over the tuning walks only
R{end+1} = '';
R{end+1} = 'ONE-NUMBER SPEED SCALE (calibration, reported separately from the main result)';
R{end+1} = sprintf('  Scale fitted on the tuning walks only: k = %.4f (the radar speed is multiplied by k, and so is the stride length, which is the integral of the speed).', kScale);
R{end+1} = sprintf('  Held-out walking speed accuracy:  %.2f -> %.2f', accuracyOf(estS(isTest), gtS(isTest)), accuracyOf(kScale * estS(isTest), gtS(isTest)));
R{end+1} = sprintf('  Held-out stride length accuracy:  %.2f -> %.2f', accuracyOf(estL(isTest), gtL(isTest)), accuracyOf(kScale * estL(isTest), gtL(isTest)));
for q = 1:numel(parts)
    m = isTest & strcmp(part, parts{q});
    R{end+1} = sprintf('    held-out %-4s (n = %2d): speed %.2f -> %.2f, stride length %.2f -> %.2f', parts{q}, sum(m), accuracyOf(estS(m), gtS(m)), accuracyOf(kScale * estS(m), gtS(m)), ...
        accuracyOf(estL(m), gtL(m)), accuracyOf(kScale * estL(m), gtL(m))); %#ok<AGROW>
end

%% 4. Many random 75/25 splits
rng(cfg.Seed, 'twister');
A0 = zeros(cfg.nRandomSplits, nP);                        % held-out accuracy without the scale
A1 = zeros(cfg.nRandomSplits, nP);                        % held-out accuracy with the scale (speed and stride length only)
kAll = zeros(cfg.nRandomSplits, 1);
estAll = cell(1, nP);
gtAll = cell(1, nP);
for p = 1:nP
    [estAll{p}, gtAll{p}] = valuesOfConfig(T, fin, wav, params(p, :));
end
for s = 1:cfg.nRandomSplits
    tst = false(n, 1);
    for q = 1:numel(parts)
        for b = 1:4
            idx = find(strcmp(part, parts{q}) & pace == b);
            nt = floor(0.25 * numel(idx) + 0.5);
            r = randperm(numel(idx));
            tst(idx(r(1:nt))) = true;
        end
    end
    kk = 1 / mean(estAll{1}(~tst) ./ gtAll{1}(~tst));
    kAll(s) = kk;
    for p = 1:nP
        A0(s, p) = accuracyOf(estAll{p}(tst), gtAll{p}(tst));
        if p == 1 || p == 4
            A1(s, p) = accuracyOf(kk * estAll{p}(tst), gtAll{p}(tst));
        else
            A1(s, p) = A0(s, p);
        end
    end
end
R{end+1} = '';
R{end+1} = sprintf('%d RANDOM 75/25 SPLITS (also split by participant and pace, seed %d), held-out accuracy of Final', cfg.nRandomSplits, cfg.Seed);
R{end+1} = sprintf('  %-26s %9s %9s %9s %9s', '', params{:, 1});
R{end+1} = sprintf('  %-26s %9.2f %9.2f %9.2f %9.2f', 'without the scale: mean', mean(A0, 1));
R{end+1} = sprintf('  %-26s %9.2f %9.2f %9.2f %9.2f', '   2.5th percentile', prctileSimple(A0, 2.5));
R{end+1} = sprintf('  %-26s %9.2f %9.2f %9.2f %9.2f', '   97.5th percentile', prctileSimple(A0, 97.5));
R{end+1} = sprintf('  %-26s %9.1f %9.1f %9.1f %9.1f', '   splits at or above 90%', 100 * mean(A0 >= cfg.PassAccuracy_pct, 1));
R{end+1} = sprintf('  %-26s %9.2f %9.2f %9.2f %9.2f', 'with the scale: mean', mean(A1, 1));
R{end+1} = sprintf('  %-26s %9.2f %9.2f %9.2f %9.2f', '   2.5th percentile', prctileSimple(A1, 2.5));
R{end+1} = sprintf('  %-26s %9.2f %9.2f %9.2f %9.2f', '   97.5th percentile', prctileSimple(A1, 97.5));
R{end+1} = sprintf('  %-26s %9.1f %9.1f %9.1f %9.1f', '   splits at or above 90%', 100 * mean(A1 >= cfg.PassAccuracy_pct, 1));
R{end+1} = sprintf('  Scale k over the splits: mean %.4f, standard deviation %.4f. A stable k means the bias is systematic, not a few walks.', mean(kAll), std(kAll));
R{end+1} = '';
R{end+1} = 'HOW TO REPORT THIS: the main result is Final without the scale. The scale is one fitted number, tested on walks that were not used to fit it,';
R{end+1} = 'and it corrects a bias that the video reference may also share (e.g. steps not landing exactly on the 0.65 m markers), so it is labelled a calibration.';
fid = fopen(fullfile(outDir, 'SplitCheck_Report.txt'), 'w');
fprintf(fid, '%s\n', R{:});
fclose(fid);
fprintf('%s\n', R{:});
fprintf('\nWritten to %s\n', outDir);
end

function [estV, gtV] = valuesOfConfig(T, cfgName, wav, par)
estV = NaN(numel(wav), 1);
gtV = NaN(numel(wav), 1);
inCfg = find(strcmp(T.config, cfgName));
for i = 1:numel(wav)
    j = inCfg(strcmpi(T.wav_file(inCfg), wav{i}));
    if ~isempty(j)
        estV(i) = T.(par{2})(j(1));
        gtV(i) = T.(par{3})(j(1));
    end
end
end

function v = prctileSimple(X, pc)
v = zeros(1, size(X, 2));
for c = 1:size(X, 2)
    x = sort(X(:, c));
    pos = 1 + (numel(x) - 1) * pc / 100;
    lo = floor(pos);
    hi = ceil(pos);
    v(c) = x(lo) + (pos - lo) * (x(hi) - x(lo));
end
end

function out = ternary(c, a, b)
if c
    out = a;
else
    out = b;
end
end

function T = readCsvStruct(file)
txt = fileread(file);
txt = strrep(txt, sprintf('\r\n'), sprintf('\n'));
txt = strrep(txt, sprintf('\r'), sprintf('\n'));
if ~isempty(txt) && double(txt(1)) == 65279
    txt = txt(2:end);                                    % byte order mark
end
if length(txt) >= 3 && all(double(txt(1:3)) == [239 187 191])
    txt = txt(4:end);
end
lines = strsplit(txt, sprintf('\n'));
lines = lines(~cellfun(@isempty, lines));
hdr = splitCsvLine(lines{1});
nC = numel(hdr);
nR = numel(lines) - 1;
raw = repmat({''}, nR, nC);
for r = 1:nR
    f = splitCsvLine(lines{r + 1});
    m = min(numel(f), nC);
    raw(r, 1:m) = f(1:m);
end
for j = 1:nC
    name = regexprep(strtrim(hdr{j}), '[^A-Za-z0-9_]', '_');
    col = raw(:, j);
    nonEmpty = ~cellfun(@isempty, col);
    v = str2double(col);
    isNum = any(nonEmpty) && all(~isnan(v(nonEmpty)) | strcmpi(strtrim(col(nonEmpty)), 'nan'));
    if isNum
        T.(name) = v;
    else
        T.(name) = col;
    end
end
end

function f = splitCsvLine(line)
f = {};
cur = '';
inQuote = false;
i = 1;
n = length(line);
while i <= n
    ch = line(i);
    if inQuote
        if ch == '"'
            if i < n && line(i + 1) == '"'
                cur(end + 1) = '"'; %#ok<AGROW>
                i = i + 1;
            else
                inQuote = false;
            end
        else
            cur(end + 1) = ch; %#ok<AGROW>
        end
    else
        if ch == '"'
            inQuote = true;
        elseif ch == ','
            f{end + 1} = cur; %#ok<AGROW>
            cur = '';
        else
            cur(end + 1) = ch; %#ok<AGROW>
        end
    end
    i = i + 1;
end
f{end + 1} = cur;
end

function T = applyGroundTruth(T, gtFile)
G = readCsvStruct(gtFile);
fn = fieldnames(G);
prevN = -1;
for q = 1:numel(fn)
    tok = regexp(fn{q}, '^v(\d+)_avgSpeed_logged_mps$', 'tokens');
    if ~isempty(tok)
        prevN = max(prevN, str2double(tok{1}{1}));
    end
end
n = numel(T.wav_file);
T.gtChanged = zeros(n, 1);
T.gtSpeed_prev = NaN(n, 1);
T.gtStepCount_prev = NaN(n, 1);
T.gtStrideDuration_prev = NaN(n, 1);
T.gtStrideLength_prev = NaN(n, 1);
T.gtStrideDuration_timed = NaN(n, 1);
if prevN >= 0
    p = sprintf('v%d_', prevN);
    T.prevPrefix = p;
    flag = sprintf('gtChanged_v%d', prevN + 1);
else
    T.prevPrefix = '(none)';
end
for i = 1:n
    j = find(strcmpi(G.wav_file, T.wav_file{i}), 1);
    assert(~isempty(j), '%s is not in %s.', T.wav_file{i}, gtFile);
    T.gtSpeed_video_mps(i) = G.avgSpeed_logged_mps(j);
    T.gtStepCount(i) = G.groundTruthStepCount(j);
    T.gtStrideDuration_s(i) = G.groundTruthStrideDuration_s(j);
    T.gtStrideLength_m(i) = G.groundTruthStrideLength_m(j);
    T.walkedDistance_m(i) = G.walkedDistance_m(j);
    T.walkType{i} = strtrim(G.walkType{j});
    if isfield(G, 'strideDuration_trialLog_s')
        T.gtStrideDuration_timed(i) = G.strideDuration_trialLog_s(j);   % the nine strides timed on the video (NaN for the rest)
    end
    if prevN >= 0
        T.gtSpeed_prev(i) = G.([p 'avgSpeed_logged_mps'])(j);
        T.gtStepCount_prev(i) = G.([p 'groundTruthStepCount'])(j);
        T.gtStrideDuration_prev(i) = G.([p 'groundTruthStrideDuration_s'])(j);
        T.gtStrideLength_prev(i) = 2 * G.([p 'walkedDistance_m'])(j) / G.([p 'groundTruthStepCount'])(j);   % 2 x distance / steps
        if isfield(G, flag)
            T.gtChanged(i) = G.(flag)(j);
        end
    else
        T.gtSpeed_prev(i) = T.gtSpeed_video_mps(i);          % no previous version, so "previous" = this file
        T.gtStepCount_prev(i) = T.gtStepCount(i);
        T.gtStrideDuration_prev(i) = T.gtStrideDuration_s(i);
        T.gtStrideLength_prev(i) = T.gtStrideLength_m(i);
    end
end
end

function a = accuracyOf(est, gt)
ok = isfinite(est) & isfinite(gt) & gt ~= 0;
e = 100 * (est(ok) - gt(ok)) ./ gt(ok);
a = 100 - sqrt(mean(e.^2));
end

function u = uniqueStable(c)
u = {};
for k = 1:numel(c)
    if ~any(strcmp(u, c{k}))
        u{end+1} = c{k}; %#ok<AGROW>
    end
end
end


function [resultsFolder, groundTruthFile] = RealData_ResolveReportInputs_v1(resultsFolder, groundTruthFile)
root = fileparts(mfilename('fullpath'));
if nargin < 1, resultsFolder = ''; end
if nargin < 2, groundTruthFile = ''; end
usingCurrent = isempty(resultsFolder);
if isempty(resultsFolder) || isempty(groundTruthFile)
    pointer = fullfile(root, 'Current_Measured_Results.json');
    assert(isfile(pointer), 'EEE4022S:MissingDataPointer', ...
        'Supply ResultsFolder and GroundTruthFile, or restore %s.', pointer);
    state = jsondecode(fileread(pointer));
    if isempty(resultsFolder), resultsFolder = state.results_folder; end
    if isempty(groundTruthFile), groundTruthFile = state.ground_truth_file; end
end
resultsFolder = char(resultsFolder);
groundTruthFile = char(groundTruthFile);
if ~isfolder(resultsFolder), resultsFolder = fullfile(root, resultsFolder); end
if ~isfile(groundTruthFile), groundTruthFile = fullfile(root, groundTruthFile); end
csvFile = fullfile(resultsFolder, 'Results_AllConfigs.csv');
assert(isfile(csvFile), 'EEE4022S:MissingResults', 'Results file not found: %s.', csvFile);
assert(isfile(groundTruthFile), 'EEE4022S:MissingReferences', 'Reference file not found: %s.', groundTruthFile);
importOpts = detectImportOptions(csvFile);
textFields = intersect({'wav_file','config','errorMsg','referenceWAV_file'}, importOpts.VariableNames);
importOpts = setvartype(importOpts, textFields, 'string');
T = readtable(csvFile, importOpts);
G = readtable(groundTruthFile, 'TextType', 'string');
assert(all(ismember({'wav_file','config','recordingDuration_s','errorMsg'}, T.Properties.VariableNames)), ...
    'EEE4022S:InvalidResults', 'Required results columns are missing in %s.', csvFile);
if usingCurrent
    final = T.config == string(state.configuration);
    assert(sum(final) == state.nNormal && numel(unique(lower(T.wav_file(final)))) == state.nNormal ...
        && all(ismissing(T.errorMsg(final)) | strlength(strtrim(T.errorMsg(final))) == 0), ...
        'EEE4022S:IncompleteCurrentRun', 'The current pointer does not select a complete %d-recording run: %s.', state.nNormal, resultsFolder);
end
if all(ismember({'wav_total_samples','wav_sample_rate_Hz'}, G.Properties.VariableNames))
    [found, j] = ismember(lower(T.wav_file), lower(G.wav_file));
    assert(all(found), 'EEE4022S:MissingTrialReference', 'A recording in %s has no reference in %s.', resultsFolder, groundTruthFile);
    valid = ismissing(T.errorMsg) | strlength(strtrim(T.errorMsg)) == 0;
    delta = abs(T.recordingDuration_s .* G.wav_sample_rate_Hz(j) - G.wav_total_samples(j));
    bad = find(valid & (~isfinite(delta) | delta > 0.5), 1);
    if ~isempty(bad)
        error('EEE4022S:WAVReferenceMismatch', ...
            '%s: retained WAV duration differs. Results: %s. References: %s. Select a matching completed run. Rerun the extractor only if the selected WAV has actually changed.', ...
            T.wav_file(bad), resultsFolder, groundTruthFile);
    end
    if all(ismember({'referenceWAV_file'}, T.Properties.VariableNames)) && ismember('wav_source_file', G.Properties.VariableNames)
        assert(all(~valid | strcmpi(T.referenceWAV_file, G.wav_source_file(j))), ...
            'EEE4022S:WAVReferenceMismatch', 'Retained WAV paths differ between %s and %s.', resultsFolder, groundTruthFile);
    end
end
fprintf('Report results: %s\nWAV references: %s\n', resultsFolder, groundTruthFile);
end
