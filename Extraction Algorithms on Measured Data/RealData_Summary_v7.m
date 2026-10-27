%% Full measured-data summary
% Zahraa Sattar | EEE4022S
%
% Regenerate report tables from the retained full-dataset results.
% This reads saved estimates and references without reprocessing the WAV files.

function RealData_Summary_v7(userCfg)

cfg.ResultsFolder = fullfile(fileparts(mfilename('fullpath')), 'Report Results');                                  % empty = reviewed current-data pointer
cfg.SNRFile = '';                                        % empty = the newest SNREstimate_ folder's RealData_SNR_Pth.csv
cfg.PassAccuracy_pct = 90;                               % pass mark for all four parameters
cfg.PaceEdges_mps = [0.8 1.0 1.2];                       % pace region boundaries
cfg.GroundTruthFile = fullfile(fileparts(mfilename('fullpath')), 'Report Results', 'Report_References.csv');                                % empty = current pointer reference file
cfg.ExcludeDouble = true;                                % leave the double walks out of the main results
cfg.IncludeArchivedVideoTiming = false;                  % optional historical check, excluded from Chapter 5
if nargin > 0 && isstruct(userCfg)
    fn = fieldnames(userCfg);
    for q = 1:numel(fn)
        cfg.(fn{q}) = userCfg.(fn{q});
    end
end
[cfg.ResultsFolder, cfg.GroundTruthFile] = RealData_ResolveReportInputs_v1(cfg.ResultsFolder, cfg.GroundTruthFile);
T = readCsvStruct(fullfile(cfg.ResultsFolder, 'Results_AllConfigs.csv'));
T = applyGroundTruth(T, cfg.GroundTruthFile);            % ground truth from the file, estimates from the run
isDouble = strcmpi(strtrim(T.walkType), 'Double');       % walks with a turn in the middle
if cfg.ExcludeDouble
    inMain = ~isDouble;
    mainName = 'Single-pass walks';
else
    inMain = true(size(isDouble));
    mainName = 'All walks';
end
[~, gtName] = fileparts(cfg.GroundTruthFile);
outDir = fullfile(cfg.ResultsFolder, ['Summary_v7_' strrep(gtName, 'RealData_', '')]);
if ~cfg.ExcludeDouble
    outDir = [outDir '_withDouble'];
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Rows that ran, and the pace region of every row
hasError = ~cellfun(@isempty, T.errorMsg);
paceNames = {'Slow', 'Comfortable', 'Fast', 'Very fast'};
paceIdx = 1 + (T.gtSpeed_video_mps >= cfg.PaceEdges_mps(1)) + (T.gtSpeed_video_mps >= cfg.PaceEdges_mps(2)) ...
    + (T.gtSpeed_video_mps >= cfg.PaceEdges_mps(3));    % 1 to 4
T.paceRegion = paceNames(paceIdx);
T.paceRegion = T.paceRegion(:);
configs = uniqueStable(T.config);                        % in the order the runner wrote them
participants = sort(unique(T.participant));

params = {'Walking speed',   'estSpeed_mps',             'gtSpeed_video_mps'; ...
          'Step count',      'estStepCount',             'gtStepCount'; ...
          'Stride duration', 'estStrideDuration_mean_s', 'gtStrideDuration_s'; ...
          'Stride length',   'estStrideLength_mean_m',   'gtStrideLength_m'};
nP = size(params, 1);

walkTypes = sort(unique(strtrim(T.walkType)));
groupType = [{'All'}, repmat({'Pace'}, 1, 4), repmat({'Participant'}, 1, numel(participants)), repmat({'Walk type'}, 1, numel(walkTypes))];
groupName = [{mainName}, paceNames, participants(:).', walkTypes(:).'];
groupMask = cell(1, numel(groupName));
groupMask{1} = inMain;
for g = 1:4
    groupMask{1 + g} = strcmp(T.paceRegion, paceNames{g}) & inMain;   % pace and participant groups use the main walks only
end
for g = 1:numel(participants)
    groupMask{5 + g} = strcmp(T.participant, participants{g}) & inMain;
end
for g = 1:numel(walkTypes)
    groupMask{5 + numel(participants) + g} = strcmp(strtrim(T.walkType), walkTypes{g});   % single pass or double walk with a turn
end

%% Accuracy table
acc = NaN(numel(configs), nP, numel(groupName));         % kept for the report
fidA = fopen(fullfile(outDir, 'Summary_Accuracy.csv'), 'w');
fprintf(fidA, 'config,groupType,group,parameter,n,nMissing,accuracy_pct,RMSE_pct,MAPE_pct,bias_pct,within10_pct,passMark_pct,pass\n');
for k = 1:numel(configs)
    isCfg = strcmp(T.config, configs{k}) & ~hasError;
    for g = 1:numel(groupName)
        m = isCfg & groupMask{g};
        for p = 1:nP
            st = errorStats(T.(params{p, 2})(m), T.(params{p, 3})(m));
            acc(k, p, g) = st.accuracy;
            fprintf(fidA, '%s,%s,%s,%s,%d,%d,%.2f,%.2f,%.2f,%.2f,%.1f,%g,%d\n', configs{k}, groupType{g}, groupName{g}, params{p, 1}, ...
                st.n, st.nMissing, st.accuracy, st.rmse, st.mape, st.bias, st.within10, cfg.PassAccuracy_pct, st.accuracy >= cfg.PassAccuracy_pct);
        end
    end
end
fclose(fidA);

%% Step count table
fidS = fopen(fullfile(outDir, 'Summary_StepCount.csv'), 'w');
fprintf(fidS, 'config,groupType,group,n,MAE_steps,meanSigned_steps,exact_pct,withinOne_pct,overCount_pct,underCount_pct\n');
stepTab = cell(numel(configs), numel(groupName));
for k = 1:numel(configs)
    isCfg = strcmp(T.config, configs{k}) & ~hasError;
    for g = 1:numel(groupName)
        m = isCfg & groupMask{g};
        d = T.estStepCount(m) - T.gtStepCount(m);        % extra steps (+) or missed steps (-)
        d = d(isfinite(d));
        s.n = numel(d);
        s.mae = mean(abs(d));
        s.signed = mean(d);
        s.exact = 100 * mean(d == 0);
        s.within1 = 100 * mean(abs(d) <= 1);
        s.over = 100 * mean(d > 0);
        s.under = 100 * mean(d < 0);
        stepTab{k, g} = s;
        fprintf(fidS, '%s,%s,%s,%d,%.2f,%.2f,%.1f,%.1f,%.1f,%.1f\n', configs{k}, groupType{g}, groupName{g}, s.n, s.mae, s.signed, s.exact, s.within1, s.over, s.under);
    end
end
fclose(fidS);

%% Tc against the trial log stride duration (a check of the Tc estimate, not a result)
fidT = fopen(fullfile(outDir, 'Summary_Tc.csv'), 'w');
fprintf(fidT, 'config,n,fallbacks,medianError_pct,meanAbsError_pct,within10_pct\n');
tcTab = zeros(numel(configs), 5);
for k = 1:numel(configs)
    m = strcmp(T.config, configs{k}) & ~hasError & inMain;
    e = 100 * (T.Tc_used_s(m) - T.gtStrideDuration_s(m)) ./ T.gtStrideDuration_s(m);
    e = e(isfinite(e));
    tcTab(k, :) = [numel(e), sum(T.TcFallback(m) == 1), median(e), mean(abs(e)), 100 * mean(abs(e) <= 10)];
    fprintf(fidT, '%s,%d,%d,%.2f,%.2f,%.1f\n', configs{k}, tcTab(k, :));
end
fclose(fidT);

%% SNR of every walk (from EstimateSNR_RealData_v1)
snr = loadSNR(cfg.SNRFile);
lastCfg = configs{end};                                   % 'All', the final pipeline
mFinal = strcmp(T.config, lastCfg) & ~hasError & inMain;
mOld = strcmp(T.config, configs{1}) & ~hasError & inMain;
wavFinal = T.wav_file(mFinal);
snrVal = NaN(numel(wavFinal), 1);
snrFlag = repmat({''}, numel(wavFinal), 1);
if ~isempty(snr)
    for a = 1:numel(wavFinal)
        j = find(strcmpi(snr.wav_file, wavFinal{a}), 1);
        if ~isempty(j)
            snrVal(a) = snr.SNRequiv_dB(j);
            snrFlag{a} = snr.SNRflag{j};
        end
    end
end
idxFinal = find(mFinal);
idxOld = find(mOld);
errFinal = NaN(numel(idxFinal), nP);
errOld = NaN(numel(idxFinal), nP);
for a = 1:numel(idxFinal)
    i = idxFinal(a);
    o = idxOld(strcmpi(T.wav_file(idxOld), T.wav_file{i}));
    for p = 1:nP
        errFinal(a, p) = 100 * (T.(params{p, 2})(i) - T.(params{p, 3})(i)) / T.(params{p, 3})(i);
        if ~isempty(o)
            errOld(a, p) = 100 * (T.(params{p, 2})(o(1)) - T.(params{p, 3})(o(1))) / T.(params{p, 3})(o(1));
        end
    end
end
fidW = fopen(fullfile(outDir, 'SNR_PerWalk.csv'), 'w');
fprintf(fidW, 'wav_file,participant,paceRegion,walkedDistance_m,SNRequiv_dB,SNRflag,errSpeed_pct_%s,errSteps_pct_%s,errStrideDur_pct_%s,errStrideLen_pct_%s,errSpeed_pct_Old,errSteps_pct_Old,errStrideDur_pct_Old,errStrideLen_pct_Old\n', ...
    lastCfg, lastCfg, lastCfg, lastCfg);
for a = 1:numel(idxFinal)
    i = idxFinal(a);
    fprintf(fidW, '"%s","%s","%s",%g,%.2f,"%s",%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n', T.wav_file{i}, T.participant{i}, T.paceRegion{i}, ...
        T.walkedDistance_m(i), snrVal(a), snrFlag{a}, errFinal(a, :), errOld(a, :));
end
fclose(fidW);

%% Report
R = {};
R{end+1} = sprintf('REAL DATA SUMMARY v7   folder %s   made %s', outDir, datestr(now, 'yyyy-mm-dd HH:MM'));
R{end+1} = sprintf('Ground truth file: %s (%d walks changed from the previous version).', cfg.GroundTruthFile, numel(unique(lower(T.wav_file(T.gtChanged == 1)))));
R{end+1} = sprintf('Main results: %s, %d walks (limp walks left out by the runner, %d double walks %s), %d rows failed.', mainName, ...
    numel(unique(lower(T.wav_file(inMain)))), numel(unique(lower(T.wav_file(isDouble)))), ternary(cfg.ExcludeDouble, 'reported on their own', 'included'), sum(hasError));
R{end+1} = sprintf('Accuracy = 100 - RMSE%%, RMSE%% = root mean square of the per-walk %% errors. Pass mark %g%% for all four parameters.', cfg.PassAccuracy_pct);
if T.referenceUsesWAV
    R{end+1} = 'References: retained trial-log distance/count. Speed = distance / selected WAV duration. Stride duration = 2 x WAV duration / count. Stride length = 2 x distance / count. Logged walking time was not used.';
else
    R{end+1} = 'References were read from the specified ground-truth file. Timing provenance follows that reference version.';
end
R{end+1} = '';
R{end+1} = sprintf('POOLED, %s (accuracy %%, * = passes)', upper(mainName));
R{end+1} = sprintf('  %-10s %16s %16s %16s %16s   passes', 'setting', params{:, 1});
for k = 1:numel(configs)
    line = sprintf('  %-10s', configs{k});
    nPass = 0;
    for p = 1:nP
        a = acc(k, p, 1);
        line = [line sprintf(' %15.2f%s', a, passMark(a, cfg.PassAccuracy_pct))]; %#ok<AGROW>
        nPass = nPass + (a >= cfg.PassAccuracy_pct);
    end
    R{end+1} = sprintf('%s   %d of %d', line, nPass, nP); %#ok<AGROW>
end
R{end+1} = '';
R{end+1} = 'CHANGE FROM Old (accuracy points, + = better)';
for k = 2:numel(configs)
    R{end+1} = sprintf('  %-10s %+16.2f %+16.2f %+16.2f %+16.2f', configs{k}, squeeze(acc(k, :, 1)) - squeeze(acc(1, :, 1))); %#ok<AGROW>
end
kb = find(strcmp(configs, 'F12'), 1);                    % round 4 reference
if isempty(kb)
    kb = find(strcmp(configs, 'F8'), 1);                 % round 3 reference
end
if ~isempty(kb)
    R{end+1} = '';
    R{end+1} = sprintf('CHANGE FROM %s (each candidate on its own, accuracy points, + = better)', configs{kb});
    for k = kb+1:numel(configs)
        R{end+1} = sprintf('  %-10s %+16.2f %+16.2f %+16.2f %+16.2f', configs{k}, squeeze(acc(k, :, 1)) - squeeze(acc(kb, :, 1))); %#ok<AGROW>
    end
end
for showK = unique([1, numel(configs)])
    R{end+1} = ''; %#ok<AGROW>
    R{end+1} = sprintf('BY GROUP, setting %s (accuracy %%, n walks)', configs{showK}); %#ok<AGROW>
    R{end+1} = sprintf('  %-14s %16s %16s %16s %16s   n', 'group', params{:, 1}); %#ok<AGROW>
    for g = 2:numel(groupName)
        nG = sum(strcmp(T.config, configs{showK}) & ~hasError & groupMask{g});
        R{end+1} = sprintf('  %-14s %15.2f%s %15.2f%s %15.2f%s %15.2f%s   %d', groupName{g}, ...
            acc(showK, 1, g), passMark(acc(showK, 1, g), cfg.PassAccuracy_pct), acc(showK, 2, g), passMark(acc(showK, 2, g), cfg.PassAccuracy_pct), ...
            acc(showK, 3, g), passMark(acc(showK, 3, g), cfg.PassAccuracy_pct), acc(showK, 4, g), passMark(acc(showK, 4, g), cfg.PassAccuracy_pct), nG); %#ok<AGROW>
    end
end
R{end+1} = '';
R{end+1} = sprintf('STEP COUNT, %s (steps: + = extra, - = missed)', upper(mainName));
R{end+1} = sprintf('  %-10s %6s %8s %8s %10s %8s %8s', 'setting', 'MAE', 'signed', 'exact%', 'within1%', 'over%', 'under%');
for k = 1:numel(configs)
    s = stepTab{k, 1};
    R{end+1} = sprintf('  %-10s %6.2f %8.2f %8.1f %10.1f %8.1f %8.1f', configs{k}, s.mae, s.signed, s.exact, s.within1, s.over, s.under); %#ok<AGROW>
end
R{end+1} = '';
R{end+1} = 'Tc USED BY EACH SETTING AGAINST THE TRIAL LOG STRIDE DURATION (check of the Tc estimate)';
R{end+1} = sprintf('  %-10s %5s %10s %10s %12s %10s', 'setting', 'n', 'fallbacks', 'median%', 'mean|err|%', 'within10%');
for k = 1:numel(configs)
    R{end+1} = sprintf('  %-10s %5d %10d %10.2f %12.2f %10.1f', configs{k}, tcTab(k, :)); %#ok<AGROW>
end
R{end+1} = '';
R{end+1} = 'SNR OF EVERY WALK (EstimateSNR_RealData_v1, spectrogram domain, 0.186 s window)';
if isempty(snr) || all(isnan(snrVal))
    R{end+1} = '  No SNR file found, so this part is skipped.';
else
    R{end+1} = sprintf('  SNR file: %s', snr.file);
    R{end+1} = '  The estimate is trustworthy from about 20 to 60 dB and saturates above about 70 dB, so values above 60 dB mean "at least 60 dB".';
    R{end+1} = snrLine('All walks', snrVal);
    for g = 1:numel(participants)
        R{end+1} = snrLine(participants{g}, snrVal(strcmp(T.participant(idxFinal), participants{g}))); %#ok<AGROW>
    end
    for g = 1:4
        R{end+1} = snrLine(paceNames{g}, snrVal(strcmp(T.paceRegion(idxFinal), paceNames{g}))); %#ok<AGROW>
    end
    dists = unique(T.walkedDistance_m(idxFinal));
    for g = 1:numel(dists)
        R{end+1} = snrLine(sprintf('%.2f m walks', dists(g)), snrVal(T.walkedDistance_m(idxFinal) == dists(g))); %#ok<AGROW>
    end
    R{end+1} = '';
    R{end+1} = sprintf('  Setting %s, accuracy by SNR third (lowest, middle, highest third of walks):', lastCfg);
    ok = isfinite(snrVal);
    edges = quantileSimple(snrVal(ok), [1/3 2/3]);
    bins = {snrVal < edges(1), snrVal >= edges(1) & snrVal < edges(2), snrVal >= edges(2)};
    binName = {sprintf('below %.1f dB', edges(1)), sprintf('%.1f to %.1f dB', edges(1), edges(2)), sprintf('above %.1f dB', edges(2))};
    for b = 1:3
        line = sprintf('    %-20s', binName{b});
        for p = 1:nP
            e = errFinal(bins{b} & ok, p);
            e = e(isfinite(e));
            line = [line sprintf(' %s %.2f', params{p, 1}, 100 - sqrt(mean(e.^2)))]; %#ok<AGROW>
        end
        R{end+1} = sprintf('%s   (n = %d)', line, sum(bins{b} & ok)); %#ok<AGROW>
    end
    for p = 1:nP
        e = abs(errFinal(:, p));
        k2 = ok & isfinite(e);
        r = corrSimple(snrVal(k2), e(k2));
        R{end+1} = sprintf('  Correlation between SNR and |error| of %s: r = %.2f (n = %d)', params{p, 1}, r, sum(k2)); %#ok<AGROW>
    end
end
R = [R, extraSections(T, configs, params, inMain, isDouble, hasError, cfg), doubleWalkSection(T, configs, params, isDouble, hasError), sensitivitySection(T, configs{end}, inMain, hasError), literatureSections(T, configs{end}, inMain, hasError, cfg)];
fid = fopen(fullfile(outDir, 'Summary_Report.txt'), 'w');
fprintf(fid, '%s\n', R{:});
fclose(fid);
fprintf('%s\n', R{:});

%% Figures (larger fonts, as asked for in the meeting)
try
    fig = figure('Visible', 'off', 'Position', [100 100 1000 500]);
    bar(squeeze(acc(:, :, 1)));                          % one group of bars per setting
    hold on;
    plot([0.5 numel(configs) + 0.5], cfg.PassAccuracy_pct * [1 1], 'k--', 'LineWidth', 1.5);   % pass mark
    set(gca, 'XTick', 1:numel(configs), 'XTickLabel', configs, 'FontSize', 14);
    ylabel('Accuracy (%)');
    ylim([max(0, min(min(acc(:, :, 1))) - 5), 100]);
    legend([params(:, 1); {sprintf('%g%% pass mark', cfg.PassAccuracy_pct)}], 'Location', 'southoutside', 'Orientation', 'horizontal');
    grid on;
    saveas(fig, fullfile(outDir, 'Fig_Accuracy_ByConfig.png'));
    close(fig);
    if ~isempty(snr) && any(isfinite(snrVal))
        fig = figure('Visible', 'off', 'Position', [100 100 800 450]);
        hold on;
        mk = {'o', 's', '^', 'd', 'v'};
        for g = 1:numel(participants)
            m = strcmp(T.participant(idxFinal), participants{g});
            plot(T.gtSpeed_video_mps(idxFinal(m)), snrVal(m), mk{1 + mod(g - 1, numel(mk))}, 'MarkerSize', 8, 'LineWidth', 1.2);
        end
        xlabel('Reference walking speed (m/s)');
        ylabel('Estimated SNR (dB)');
        legend(participants, 'Location', 'best');
        set(gca, 'FontSize', 14);
        grid on;
        saveas(fig, fullfile(outDir, 'Fig_SNR_PerWalk.png'));
        close(fig);
    end
    drawDoubleWalkFigure(T, configs, params, isDouble, hasError, cfg, outDir);
catch ME
    fprintf('Figures skipped: %s\n', ME.message);
end
fprintf('\nWritten to %s\n', outDir);
end

function st = errorStats(est, gt)
valid = isfinite(est) & isfinite(gt) & gt ~= 0;         % walks that gave a value
e = 100 * (est(valid) - gt(valid)) ./ gt(valid);         % per-walk % error
st.n = sum(valid);
st.nMissing = sum(~valid);                               % e.g. fewer than 3 steps, so no stride
st.rmse = sqrt(mean(e.^2));
st.accuracy = 100 - st.rmse;
st.mape = mean(abs(e));
st.bias = mean(e);
st.within10 = 100 * mean(abs(e) <= 10);
end

function R = extraSections(T, configs, params, inMain, isDouble, hasError, cfg)
nP = size(params, 1);
fin = configs{end};
R = {};

R{end+1} = '';
R{end+1} = 'EFFECT OF THE GROUND TRUTH VERSION AND OF THE DOUBLE WALKS (accuracy %)';
R{end+1} = sprintf('  "previous" = the values the ground truth file keeps from the version before (columns %s...).', T.prevPrefix);
R{end+1} = sprintf('  %-8s %-34s %8s %8s %8s %8s   n', 'setting', 'case', 'speed', 'steps', 'strDur', 'strLen');
cases = {'previous ground truth, all walks', 'prev', true(size(inMain)); ...
         'previous ground truth, single pass', 'prev', ~isDouble; ...
         'new ground truth, all walks', 'now', true(size(inMain)); ...
         'new ground truth, single pass', 'now', ~isDouble};
for k = [1, numel(configs)]
    for c = 1:size(cases, 1)
        m = strcmp(T.config, configs{k}) & ~hasError & cases{c, 3};
        a = zeros(1, nP);
        for p = 1:nP
            a(p) = accuracyOf(T.(params{p, 2})(m), gtOf(T, params{p, 3}, cases{c, 2}, m));
        end
        R{end+1} = sprintf('  %-8s %-34s %8.2f %8.2f %8.2f %8.2f   %d', configs{k}, cases{c, 1}, a, sum(m)); %#ok<AGROW>
    end
end

ch = find(strcmp(T.config, fin) & ~hasError & T.gtChanged == 1);
R{end+1} = '';
R{end+1} = sprintf('WALKS WHOSE GROUND TRUTH CHANGED, setting %s (%% error before -> after the change)', fin);
R{end+1} = sprintf('  %-11s %10s %14s %9s %9s %16s %16s %16s', 'walk', 'steps', 'speed (m/s)', 'est spd', 'est stp', 'speed err %', 'step err %', 'strDur err %');
for a = 1:numel(ch)
    i = ch(a);
    eS = 100 * (T.estSpeed_mps(i) ./ [T.gtSpeed_prev(i), T.gtSpeed_video_mps(i)] - 1);
    eN = 100 * (T.estStepCount(i) ./ [T.gtStepCount_prev(i), T.gtStepCount(i)] - 1);
    eD = 100 * (T.estStrideDuration_mean_s(i) ./ [T.gtStrideDuration_prev(i), T.gtStrideDuration_s(i)] - 1);
    R{end+1} = sprintf('  %-11s %4d -> %-3d %5.3f -> %5.3f %9.3f %9d %7.1f -> %6.1f %7.1f -> %6.1f %7.1f -> %6.1f', strrep(T.wav_file{i}, '.wav', ''), ...
        T.gtStepCount_prev(i), T.gtStepCount(i), T.gtSpeed_prev(i), T.gtSpeed_video_mps(i), T.estSpeed_mps(i), T.estStepCount(i), eS, eN, eD); %#ok<AGROW>
end

R{end+1} = '';
R{end+1} = sprintf('ERROR BREAKDOWN, setting %s, main walks (per-walk %% errors)', fin);
R{end+1} = '  accuracy = 100 - RMSE; RMSE^2 = bias^2 + spread^2, so a constant bias lowers the accuracy even when every walk is close.';
R{end+1} = sprintf('  %-16s %9s %8s %8s %8s %12s %10s %10s', 'parameter', 'accuracy', 'bias', 'spread', 'MAPE', 'median|e|', 'within10%', 'over 20%');
m = strcmp(T.config, fin) & ~hasError & inMain;
for p = 1:nP
    est = T.(params{p, 2})(m);
    gt = T.(params{p, 3})(m);
    ok = isfinite(est) & isfinite(gt) & gt ~= 0;
    e = 100 * (est(ok) - gt(ok)) ./ gt(ok);
    R{end+1} = sprintf('  %-16s %9.2f %8.2f %8.2f %8.2f %12.2f %10.1f %10d', params{p, 1}, 100 - sqrt(mean(e.^2)), mean(e), std(e, 1), ...
        mean(abs(e)), median(abs(e)), 100 * mean(abs(e) <= 10), sum(abs(e) > 20)); %#ok<AGROW>
end

parts = sort(unique(T.participant(inMain)));
R{end+1} = '';
R{end+1} = 'PARTICIPANT HOLD-OUT CHECK (main walks)';
R{end+1} = '  For each participant, the setting with the highest mean accuracy (four parameters) on the OTHER two';
R{end+1} = '  participants is chosen and then scored on the held-out participant. If this does not beat Final,';
R{end+1} = '  choosing the setting on the measured data would not have generalised better than choosing it on the simulator.';
poolEst = cell(1, nP);
poolGt = cell(1, nP);
for q = 1:numel(parts)
    train = inMain & ~strcmp(T.participant, parts{q}) & ~hasError;
    test = inMain & strcmp(T.participant, parts{q}) & ~hasError;
    best = 1;
    bestScore = -Inf;
    for k = 1:numel(configs)
        mk = train & strcmp(T.config, configs{k});
        s = 0;
        for p = 1:nP
            s = s + accuracyOf(T.(params{p, 2})(mk), T.(params{p, 3})(mk)) / nP;
        end
        if s > bestScore
            bestScore = s;
            best = k;
        end
    end
    mt = test & strcmp(T.config, configs{best});
    mf = test & strcmp(T.config, fin);
    aB = zeros(1, nP);
    aF = zeros(1, nP);
    for p = 1:nP
        aB(p) = accuracyOf(T.(params{p, 2})(mt), T.(params{p, 3})(mt));
        aF(p) = accuracyOf(T.(params{p, 2})(mf), T.(params{p, 3})(mf));
        poolEst{p} = [poolEst{p}; T.(params{p, 2})(mt)];
        poolGt{p} = [poolGt{p}; T.(params{p, 3})(mt)];
    end
    R{end+1} = sprintf('  held out %-4s chosen %-8s held-out accuracy %6.2f %6.2f %6.2f %6.2f   (%s: %6.2f %6.2f %6.2f %6.2f)', parts{q}, configs{best}, aB, fin, aF); %#ok<AGROW>
end
aP = zeros(1, nP);
aF = zeros(1, nP);
mf = inMain & ~hasError & strcmp(T.config, fin);
for p = 1:nP
    aP(p) = accuracyOf(poolEst{p}, poolGt{p});
    aF(p) = accuracyOf(T.(params{p, 2})(mf), T.(params{p, 3})(mf));
end
R{end+1} = sprintf('  pooled held-out          accuracy %6.2f %6.2f %6.2f %6.2f', aP);
R{end+1} = sprintf('  %-8s (simulator)      accuracy %6.2f %6.2f %6.2f %6.2f', fin, aF);
end

function R = doubleWalkSection(T, configs, params, isDouble, hasError, cfg)
nP = size(params, 1);
fin = configs{end};
R = {};
R{end+1} = '';
R{end+1} = 'DOUBLE WALKS INCLUDED OR LEFT OUT (accuracy %, 100 - RMSE%)';
R{end+1} = '  "single" = single-pass walks only. "all" = single-pass plus double walks. "double" = the double walks on their own.';
R{end+1} = '  change = all - single (negative: the double walks pull the accuracy down). The double walks have a turn in the middle,';
R{end+1} = '  which the algorithm does not model, so they are reported separately and are not part of the main result.';
R{end+1} = sprintf('  %-8s %-9s %9s %9s %9s %9s   n', 'setting', 'case', params{:, 1});
for k = unique([1, numel(configs)])
    cases = {'single', ~isDouble; 'all', true(size(isDouble)); 'double', isDouble};
    a = zeros(3, nP);
    nn = zeros(3, 1);
    for c = 1:3
        m = strcmp(T.config, configs{k}) & ~hasError & cases{c, 2};
        nn(c) = numel(unique(lower(T.wav_file(m))));
        for p = 1:nP
            a(c, p) = accuracyOf(T.(params{p, 2})(m), T.(params{p, 3})(m));
        end
        R{end+1} = sprintf('  %-8s %-9s %9.2f %9.2f %9.2f %9.2f   %d', configs{k}, cases{c, 1}, a(c, :), nn(c)); %#ok<AGROW>
    end
    R{end+1} = sprintf('  %-8s %-9s %+9.2f %+9.2f %+9.2f %+9.2f', configs{k}, 'change', a(2, :) - a(1, :)); %#ok<AGROW>
end
mf = strcmp(T.config, fin) & ~hasError;
idx = find(mf & isDouble);
R{end+1} = '';
R{end+1} = sprintf('  Every double walk, setting %s (per-walk %% error; steps are counts):', fin);
R{end+1} = sprintf('  %-10s %7s %7s %10s %14s %14s', 'walk', 'gt stp', 'est stp', 'speed %', 'strideDur %', 'strideLen %');
for a = 1:numel(idx)
    i = idx(a);
    R{end+1} = sprintf('  %-10s %7d %7d %+10.1f %+14.1f %+14.1f', strrep(T.wav_file{i}, '.wav', ''), T.gtStepCount(i), T.estStepCount(i), ...
        100 * (T.estSpeed_mps(i) / T.gtSpeed_video_mps(i) - 1), 100 * (T.estStrideDuration_mean_s(i) / T.gtStrideDuration_s(i) - 1), ...
        100 * (T.estStrideLength_mean_m(i) / T.gtStrideLength_m(i) - 1)); %#ok<AGROW>
end
aS = zeros(1, nP);
aA = zeros(1, nP);
for p = 1:nP
    aS(p) = accuracyOf(T.(params{p, 2})(mf & ~isDouble), T.(params{p, 3})(mf & ~isDouble));
    aA(p) = accuracyOf(T.(params{p, 2})(mf), T.(params{p, 3})(mf));
end
R{end+1} = '';
for p = 1:nP
    d = aA(p) - aS(p);
    if abs(d) < 0.3
        w = 'hardly changes';
    elseif d < 0
        w = 'lowers';
    else
        w = 'raises';
    end
    R{end+1} = sprintf('  %-16s including the %d double walks %s the accuracy (%.2f -> %.2f, %+.2f points).', params{p, 1}, sum(isDouble(mf)), w, aS(p), aA(p), d); %#ok<AGROW>
end
end

function drawDoubleWalkFigure(T, configs, params, isDouble, hasError, cfg, outDir)
fin = configs{end};
mf = strcmp(T.config, fin) & ~hasError;
nP = size(params, 1);
fig = figure('Visible', 'off', 'Position', [100 100 1200 800]);
for p = 1:nP
    subplot(2, 2, p);
    hold on;
    est = T.(params{p, 2});
    gt = T.(params{p, 3});
    e = 100 * (est - gt) ./ gt;
    x = T.gtSpeed_video_mps;
    if p == 2
        e = est - gt;                                    % step count: extra steps, not %
    end
    s = mf & ~isDouble & isfinite(e);
    d = mf & isDouble & isfinite(e);
    plot(x(s), e(s), 'o', 'Color', [0.45 0.45 0.45], 'MarkerSize', 6);
    plot(x(d), e(d), 's', 'MarkerFaceColor', [0.80 0.20 0.20], 'MarkerEdgeColor', [0.5 0.1 0.1], 'MarkerSize', 10);
    idd = find(d);
    for a = 1:numel(idd)
        text(x(idd(a)) + 0.01, e(idd(a)), strrep(strrep(T.wav_file{idd(a)}, '.wav', ''), '_', '\_'), 'FontSize', 10);
    end
    plot([min(x(mf)) max(x(mf))], [0 0], 'k-');
    if p ~= 2
        plot([min(x(mf)) max(x(mf))], [10 10], 'k--');
        plot([min(x(mf)) max(x(mf))], [-10 -10], 'k--');
        ylabel('Error (%)');
    else
        ylabel('Extra steps (+) or missed (-)');
    end
    xlabel('Reference walking speed (m/s)');
    aS = accuracyOf(est(mf & ~isDouble), gt(mf & ~isDouble));
    aA = accuracyOf(est(mf), gt(mf));
    title(sprintf('%s: single %.1f%%, with doubles %.1f%%', params{p, 1}, aS, aA), 'FontSize', 13);
    set(gca, 'FontSize', 12);
    grid on;
end
saveas(fig, fullfile(outDir, 'Fig_DoubleWalks.png'));
close(fig);
end

function R = sensitivitySection(T, fin, inMain, hasError)
R = {};
R{end+1} = '';
R{end+1} = sprintf('STRIDE LENGTH: SCORING AND CALIBRATION SENSITIVITY, setting %s, main walks', fin);
m = strcmp(T.config, fin) & ~hasError & inMain;
est = T.estStrideLength_mean_m(m);
gt = T.gtStrideLength_m(m);
ok = isfinite(est) & isfinite(gt);
e = 100 * (est(ok) - gt(ok)) ./ gt(ok);
R{end+1} = sprintf('  Same errors, three scores: 100 - RMSE%% = %.2f (the pass mark used in this thesis), 1 - MAPE = %.2f (Hadjipanayi definition), 100 - spread = %.2f (bias removed, i.e. repeatability only).', ...
    100 - sqrt(mean(e.^2)), 100 - mean(abs(e)), 100 - std(e, 1));
eS = 100 * (T.estSpeed_mps(m) - T.gtSpeed_video_mps(m)) ./ T.gtSpeed_video_mps(m);
eD = 100 * (T.estStrideDuration_mean_s(m) - T.gtStrideDuration_s(m)) ./ T.gtStrideDuration_s(m);
eL = 100 * (T.estStrideLength_mean_m(m) - T.gtStrideLength_m(m)) ./ T.gtStrideLength_m(m);
ok2 = isfinite(eS) & isfinite(eD) & isfinite(eL);
R{end+1} = sprintf('  Stride length error = speed error + stride duration error: mean %.2f%% against %.2f%% + %.2f%% = %.2f%% (correlation of the per-walk errors %.3f).', ...
    mean(eL(ok2)), mean(eS(ok2)), mean(eD(ok2)), mean(eS(ok2)) + mean(eD(ok2)), corrSimple(eL(ok2), eS(ok2) + eD(ok2)));
R{end+1} = '  A one-off speed scale k (applied to speed and to stride length) chosen on two participants and tested on the third:';
R{end+1} = sprintf('  %-12s %8s %12s %12s %12s %12s', 'held out', 'k', 'speed before', 'speed after', 'strLen before', 'strLen after');
parts = sort(unique(T.participant(m)));
pe = []; pl = [];
for q = 1:numel(parts)
    tr = m & ~strcmp(T.participant, parts{q});
    te = m & strcmp(T.participant, parts{q});
    ratio = T.estSpeed_mps(tr) ./ T.gtSpeed_video_mps(tr);
    k = 1 / mean(ratio(isfinite(ratio)));                % scale that removes the mean speed bias of the training participants
    sB = accuracyOf(T.estSpeed_mps(te), T.gtSpeed_video_mps(te));
    sA = accuracyOf(k * T.estSpeed_mps(te), T.gtSpeed_video_mps(te));
    lB = accuracyOf(T.estStrideLength_mean_m(te), T.gtStrideLength_m(te));
    lA = accuracyOf(k * T.estStrideLength_mean_m(te), T.gtStrideLength_m(te));
    R{end+1} = sprintf('  %-12s %8.3f %12.2f %12.2f %12.2f %12.2f', parts{q}, k, sB, sA, lB, lA); %#ok<AGROW>
    pe = [pe; 100 * (k * T.estStrideLength_mean_m(te) - T.gtStrideLength_m(te)) ./ T.gtStrideLength_m(te)]; %#ok<AGROW>
end
pe = pe(isfinite(pe));
R{end+1} = sprintf('  pooled held-out stride length accuracy with the scale: %.2f%% (without: %.2f%%).', 100 - sqrt(mean(pe.^2)), 100 - sqrt(mean(e.^2)));
R{end+1} = '  Mean radar stride length and mean speed bias per participant (a different value for each participant points to the reference, e.g. steps not landing on the 0.65 m markers):';
for q = 1:numel(parts)
    te = m & strcmp(T.participant, parts{q});
    R{end+1} = sprintf('    %-4s radar stride length %.3f m (reference 1.300 m), speed bias %+.2f%%, n = %d', parts{q}, mean(T.estStrideLength_mean_m(te), 'omitnan'), ...
        mean(100 * (T.estSpeed_mps(te) ./ T.gtSpeed_video_mps(te) - 1), 'omitnan'), sum(te)); %#ok<AGROW>
end
end

function a = accuracyOf(est, gt)
ok = isfinite(est) & isfinite(gt) & gt ~= 0;
e = 100 * (est(ok) - gt(ok)) ./ gt(ok);
a = 100 - sqrt(mean(e.^2));
end

function g = gtOf(T, field, which, m)
if strcmp(which, 'now')
    g = T.(field)(m);
    return;
end
switch field
    case 'gtSpeed_video_mps'
        g = T.gtSpeed_prev(m);
    case 'gtStepCount'
        g = T.gtStepCount_prev(m);
    case 'gtStrideDuration_s'
        g = T.gtStrideDuration_prev(m);
    otherwise
        g = T.gtStrideLength_prev(m);
end
end

function T = applyGroundTruth(T, gtFile)
G = readCsvStruct(gtFile);
T.referenceUsesWAV = isfield(G, 'wav_total_samples');
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
    if isfield(G, 'reference_version')
        flag = sprintf('gtChanged_v%d', G.reference_version(1));
    end
else
    T.prevPrefix = '(none)';
end
for i = 1:n
    j = find(strcmpi(G.wav_file, T.wav_file{i}), 1);
    assert(~isempty(j), '%s is not in %s.', T.wav_file{i}, gtFile);
    if isfield(G, 'wav_total_samples')
        deltaSamples = abs(T.recordingDuration_s(i) * G.wav_sample_rate_Hz(j) - G.wav_total_samples(j));
        assert(isfinite(deltaSamples) && deltaSamples <= 0.5, 'EEE4022S:WAVReferenceMismatch', ...
            '%s: saved run duration differs from the WAV reference. Select a matching ResultsFolder and GroundTruthFile. The current reviewed run is recorded in Current_Measured_Results.json.', T.wav_file{i});
        if isfield(T, 'referenceWAV_file')
            assert(strcmpi(T.referenceWAV_file{i}, G.wav_source_file{j}), 'EEE4022S:WAVReferenceMismatch', ...
                '%s: saved run used a different WAV file.', T.wav_file{i});
        end
    end
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

function R = literatureSections(T, fin, inMain, hasError, cfg)
m = find(strcmp(T.config, fin) & ~hasError & inMain);
R = {};
R{end+1} = '';
R{end+1} = sprintf('LITERATURE-COMPARABLE METRICS, setting %s, main walks', fin);
R{end+1} = '  MAPE = mean |error| / reference; 1 - MAPE is the accuracy definition of Hadjipanayi et al. (2024).';
R{end+1} = '  ICC(A,1) = absolute agreement, ICC(C,1) = consistency, two-way, single measure (McGraw and Wong), as reported by Wang et al. (2014).';
R{end+1} = '  Bland-Altman: bias = mean(radar - reference), limits of agreement = bias +- 1.96 SD.';
R{end+1} = sprintf('  %-16s %4s %9s %9s %7s %8s %7s %8s %8s %18s', 'parameter', 'n', 'MAE', 'RMSE', 'MAPE%', '1-MAPE%', 'r', 'ICC(A,1)', 'ICC(C,1)', 'BA bias [LoA]');
cad = @(sd) 120 ./ sd;                                   % steps per minute from the stride duration
list = {'Speed (m/s)', T.estSpeed_mps(m), T.gtSpeed_video_mps(m); ...
        'Stride dur. (s)', T.estStrideDuration_mean_s(m), T.gtStrideDuration_s(m); ...
        'Cadence (st/min)', cad(T.estStrideDuration_mean_s(m)), cad(T.gtStrideDuration_s(m)); ...
        'Stride length (m)', T.estStrideLength_mean_m(m), T.gtStrideLength_m(m)};
for k = 1:size(list, 1)
    e = list{k, 2};
    g = list{k, 3};
    ok = isfinite(e) & isfinite(g);
    e = e(ok);
    g = g(ok);
    d = e - g;
    mape = 100 * mean(abs(d) ./ g);
    if std(g) > 1e-9
        r = corrSimple(e, g);
        [iA, iC] = icc21(e, g);
        rs = sprintf('%7.3f %8.3f %8.3f', r, iA, iC);
    else
        rs = sprintf('%7s %8s %8s', 'n/a', 'n/a', 'n/a');
    end
    R{end+1} = sprintf('  %-16s %4d %9.4f %9.4f %7.2f %8.2f %s   %+.4f [%+.4f, %+.4f]', list{k, 1}, numel(e), mean(abs(d)), sqrt(mean(d.^2)), ...
        mape, 100 - mape, rs, mean(d), mean(d) - 1.96 * std(d), mean(d) + 1.96 * std(d)); %#ok<AGROW>
end
dN = T.estStepCount(m) - T.gtStepCount(m);
dN = dN(isfinite(dN));
R{end+1} = sprintf('  %-16s %4d  MAE %.2f steps, exact %.1f%%, within one step %.1f%%, MAPE %.2f%%', 'Step count', numel(dN), mean(abs(dN)), ...
    100 * mean(dN == 0), 100 * mean(abs(dN) <= 1), 100 * mean(abs(dN) ./ T.gtStepCount(m(isfinite(T.estStepCount(m) - T.gtStepCount(m))))));

R{end+1} = '';
R{end+1} = sprintf('95%% BOOTSTRAP CONFIDENCE INTERVALS of the accuracy (100 - RMSE%%), setting %s, main walks', fin);
R{end+1} = '  2000 resamples of the walks with replacement, fixed seed. Walks of one participant are not independent, so these intervals are a lower bound on the true uncertainty.';
pp = {'Walking speed', T.estSpeed_mps(m), T.gtSpeed_video_mps(m); 'Step count', T.estStepCount(m), T.gtStepCount(m); ...
      'Stride duration', T.estStrideDuration_mean_s(m), T.gtStrideDuration_s(m); 'Stride length', T.estStrideLength_mean_m(m), T.gtStrideLength_m(m)};
rngOld = rng;                                            % keep the random state of the session
rng(2026, 'twister');                                    % fixed seed, the same intervals every run
for k = 1:size(pp, 1)
    ePct = 100 * (pp{k, 2} - pp{k, 3}) ./ pp{k, 3};
    ePct = ePct(isfinite(ePct));
    nB = 2000;
    a = zeros(nB, 1);
    for b = 1:nB
        idx = randi(numel(ePct), numel(ePct), 1);          % resample the walks
        a(b) = 100 - sqrt(mean(ePct(idx).^2));
    end
    q = quantileSimple(a, [0.025 0.975]);
    R{end+1} = sprintf('  %-16s accuracy %6.2f%%   95%% CI %6.2f to %6.2f%%', pp{k, 1}, 100 - sqrt(mean(ePct.^2)), q(1), q(2)); %#ok<AGROW>
end
rng(rngOld);                                             % put the random state back

R{end+1} = '';
R{end+1} = 'REFERENCE CHECKS (how the stride duration reference itself behaves)';
mt = m(isfinite(T.gtStrideDuration_timed(m)) & T.gtStrideDuration_timed(m) > 0);
if cfg.IncludeArchivedVideoTiming && ~isempty(mt)
    est = T.estStrideDuration_mean_s(mt);
    der = T.gtStrideDuration_s(mt);
    tim = T.gtStrideDuration_timed(mt);
    R{end+1} = sprintf('  Walks with a stride duration timed on the video: %d (main walks only).', numel(mt));
    R{end+1} = sprintf('  %-12s %9s %9s %9s %12s %12s', 'walk', 'derived', 'timed', 'radar', 'radar-der%', 'radar-timed%');
    for a = 1:numel(mt)
        R{end+1} = sprintf('  %-12s %9.3f %9.3f %9.3f %12.1f %12.1f', strrep(T.wav_file{mt(a)}, '.wav', ''), der(a), tim(a), est(a), ...
            100 * (est(a) / der(a) - 1), 100 * (est(a) / tim(a) - 1)); %#ok<AGROW>
    end
    eD = 100 * (est ./ der - 1);
    eT = 100 * (est ./ tim - 1);
    R{end+1} = sprintf('  Against the derived reference: bias %+.2f%%, accuracy %.2f%%. Against the timed reference: bias %+.2f%%, accuracy %.2f%%.', ...
        mean(eD), 100 - sqrt(mean(eD.^2)), mean(eT), 100 - sqrt(mean(eT.^2)));
    R{end+1} = sprintf('  Derived reference is on average %+.1f%% longer than the timed one on these walks.', 100 * mean(der ./ tim - 1));
end
dists = unique(T.walkedDistance_m(m));
R{end+1} = '  Stride duration bias by walk distance (the shorter walks were cropped closely around the counted steps):';
for k = 1:numel(dists)
    mk = m(T.walkedDistance_m(m) == dists(k));
    e = 100 * (T.estStrideDuration_mean_s(mk) ./ T.gtStrideDuration_s(mk) - 1);
    e = e(isfinite(e));
    eS = 100 * (T.estSpeed_mps(mk) ./ T.gtSpeed_video_mps(mk) - 1);
    R{end+1} = sprintf('    %.2f m walks: n = %3d, stride duration bias %+.2f%%, walking speed bias %+.2f%%', dists(k), numel(mk), mean(e), mean(eS)); %#ok<AGROW>
end
end

function [iccA, iccC] = icc21(x, y)
X = [x(:), y(:)];
[n, k] = size(X);
gm = mean(X(:));
SSR = k * sum((mean(X, 2) - gm).^2);                     % between walks
SSC = n * sum((mean(X, 1) - gm).^2);                     % between radar and reference
SST = sum((X(:) - gm).^2);
SSE = SST - SSR - SSC;
MSR = SSR / (n - 1);
MSC = SSC / (k - 1);
MSE = SSE / ((n - 1) * (k - 1));
iccC = (MSR - MSE) / (MSR + (k - 1) * MSE);              % consistency, a constant offset does not count
iccA = (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n);   % absolute agreement, an offset does count
end

function f = newestGroundTruth()
d = dir('RealData_GroundTruth_v*.csv');
assert(~isempty(d), 'No RealData_GroundTruth_v*.csv found in the Current Folder.');
v = zeros(1, numel(d));
for k = 1:numel(d)
    tok = regexp(d(k).name, '_v(\d+)\.csv$', 'tokens');
    if ~isempty(tok)
        v(k) = str2double(tok{1}{1});
    end
end
[~, b] = max(v);
f = d(b).name;
end

function s = ternary(c, a, b)
if c
    s = a;
else
    s = b;
end
end

function s = passMark(a, mark)
if a >= mark
    s = '*';
else
    s = ' ';
end
end

function s = snrLine(label, v)
v = sort(v(isfinite(v)));
if isempty(v)
    s = sprintf('  %-16s no values', label);
    return;
end
q = quantileSimple(v, [0.25 0.5 0.75]);
s = sprintf('  %-16s n = %3d   median %.1f dB   middle half %.1f to %.1f dB   range %.1f to %.1f dB', label, numel(v), q(2), q(1), q(3), v(1), v(end));
end

function q = quantileSimple(v, p)
v = sort(v(:));
n = numel(v);
q = zeros(size(p));
for k = 1:numel(p)
    pos = 1 + (n - 1) * p(k);                            % linear interpolation between sorted values
    lo = floor(pos);
    hi = min(n, lo + 1);
    q(k) = v(lo) + (pos - lo) * (v(hi) - v(lo));
end
end

function r = corrSimple(a, b)
a = a(:) - mean(a);
b = b(:) - mean(b);
r = sum(a .* b) / sqrt(sum(a.^2) * sum(b.^2));
end

function u = uniqueStable(c)
u = {};
for k = 1:numel(c)
    if ~any(strcmp(u, c{k}))
        u{end+1} = c{k}; %#ok<AGROW>
    end
end
end

function snr = loadSNR(file)
snr = [];
if isempty(file)
    d = dir(fullfile('SNREstimate_*', 'RealData_SNR_Pth.csv'));
    if isempty(d)
        return;
    end
    folders = sort({d.folder});
    file = fullfile(folders{end}, 'RealData_SNR_Pth.csv');   % newest SNR run
end
if ~exist(file, 'file')
    return;
end
snr = readCsvStruct(file);
snr.file = file;
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
