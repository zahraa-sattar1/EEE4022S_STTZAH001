%% Verify simulated extraction
% Zahraa Sattar | EEE4022S
%
% Compare the retained detector with and without the recording-edge rule.
% Quick mode checks two walking speeds with two noise draws per SNR.

function Verify_Simulated_Extraction(userCfg)

cfg.Speeds = [0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];      % m/s, as round 4
cfg.SNRs = [50 60 70 80];                                    % dB, as round 4
cfg.SelectSNRs = [60 70 80];                                 % dB, the cells gates a and b use
cfg.nTrials = 10;                                            % noise draws per speed and SNR
cfg.R0_m = 10;                                               % same range as Chapter 4
cfg.AmplitudeMode = 'A';                                     % Chen amplitude model, as in Chapter 4
cfg.Seed = 20267001;                                         % simulator seed, same walks as round 4
cfg.SeedNoise = 20269002;                                    % new noise for this round
cfg.PassAcc_pct = 95;                                        % Chapter 4 pass mark
cfg.MaxExtraME_pct = 0.5;                                    % allowed rise in miss + extra (points)
cfg.SimName = 'RadarSimulator_v11_torso_legs_arms_head_buffered_noplots';
releaseRoot = fileparts(fileparts(mfilename('fullpath')));
cfg.SelectionFile = fullfile(fileparts(mfilename('fullpath')), 'Simulated_Final_Parameters.csv');
cfg.Quick = false;
if nargin > 0 && isstruct(userCfg)
    fn = fieldnames(userCfg);
    for q = 1:numel(fn)
        cfg.(fn{q}) = userCfg.(fn{q});
    end
end
if cfg.Quick
    cfg.Speeds = [0.5 1.4];
    cfg.nTrials = 2;
end

oldPath = path;
pathCleanup = onCleanup(@() path(oldPath));
addpath(fullfile(releaseRoot, 'Radar Simulator'), fullfile(releaseRoot, 'Extraction Algorithms on Measured Data'));                                               % keep this folder on the path
outDir = fullfile(releaseRoot, 'Results', ['Simulated_' datestr(now, 'yyyymmdd_HHMMSSFFF')]);
mkdir(outDir);
logFile = fullfile(outDir, 'Run_Log.txt');
simCopy = makeSimCopy(cfg.SimName, outDir);                  % copy without clear all / close all
final = readSelection(cfg.SelectionFile);                    % the frozen Final setting
final.name = 'Final';
final.propSpacing = 1;                                       % the adopted v13 setting (proportional spacing)
final.propEdge = 0;
prop = final;
prop.name = 'FinalEdge';
prop.propEdge = 1;                                           % the only change
configs = [final, prop];
opts = struct('CPI_torso_s', final.CPI_torso_s, 'CPI_legs_s', final.CPI_legs_s, 'Pth_dB', final.Pth_dB);
logMsg(logFile, 'Round 6: Final against FinalEdge, speeds %s m/s, SNR %s dB, %d draws each.', mat2str(cfg.Speeds), mat2str(cfg.SNRs), cfg.nTrials);

tStart = tic;
rows = runAll(cfg, configs, opts, simCopy, cfg.SeedNoise, logFile, tStart);
writeCells(fullfile(outDir, 'Round6_Cells.csv'), rows);

low = setdiff(cfg.SNRs, cfg.SelectSNRs);
[pF, meF] = cellScore(rows, 'Final', cfg.SelectSNRs, cfg);
[pP, meP] = cellScore(rows, 'FinalEdge', cfg.SelectSNRs, cfg);
[pFL, meFL] = cellScore(rows, 'Final', low, cfg);
[pPL, mePL] = cellScore(rows, 'FinalEdge', low, cfg);
ga = pP >= pF;
gb = meP <= meF + cfg.MaxExtraME_pct;
gc = pPL >= pFL && mePL <= meFL + cfg.MaxExtraME_pct;
L = {};
L{end + 1} = 'SIMULATOR CHECK, ROUND 6: edge-aware prominence (RealData_Pipeline_v14)';
L{end + 1} = sprintf('Speeds %s m/s, SNR %s dB, %d noise draws each. Pass = step time, stride duration and stride length accuracy all >= %g%%.', ...
    mat2str(cfg.Speeds), mat2str(cfg.SNRs), cfg.nTrials, cfg.PassAcc_pct);
L{end + 1} = '';
L{end + 1} = sprintf('Cells at %s dB', mat2str(cfg.SelectSNRs));
L{end + 1} = sprintf('  %-9s %8s %8s %8s %10s %9s %9s %9s', 'setting', 'passes', 'miss%', 'extra%', 'perfect%', 'step%', 'strDur%', 'strLen%');
L{end + 1} = summaryLine(rows, 'Final', cfg, cfg.SelectSNRs);
L{end + 1} = summaryLine(rows, 'FinalEdge', cfg, cfg.SelectSNRs);
L{end + 1} = sprintf('Cells at %s dB', mat2str(low));
L{end + 1} = summaryLine(rows, 'Final', cfg, low);
L{end + 1} = summaryLine(rows, 'FinalEdge', cfg, low);
L{end + 1} = '';
L{end + 1} = sprintf('Gate a (passes at 60-80 dB): FinalEdge %d, Final %d -> %s', pP, pF, yesNo(ga, 'PASS', 'FAIL'));
L{end + 1} = sprintf('Gate b (miss + extra at 60-80 dB): FinalEdge %.2f%%, Final %.2f%% -> %s', meP, meF, yesNo(gb, 'PASS', 'FAIL'));
L{end + 1} = sprintf('Gate c (50 dB): FinalEdge %d cells, %.2f%%; Final %d cells, %.2f%% -> %s', pPL, mePL, pFL, meFL, yesNo(gc, 'PASS', 'FAIL'));
L{end + 1} = sprintf('GATE 1: %s', yesNo(ga && gb && gc, 'PASSED, go on to the measured check (Gate 2)', 'FAILED, Final stays'));
L{end + 1} = '';
L{end + 1} = 'Per cell (speed, SNR): Final | FinalEdge  (miss%, extra%, step%, stride duration%, stride length%)';
rF = rows(strcmp({rows.config}, 'Final'));
rP = rows(strcmp({rows.config}, 'FinalEdge'));
for k = 1:numel(rF)
    L{end + 1} = sprintf('  %.1f m/s %2d dB: %5.1f %5.1f %6.2f %6.2f %6.2f | %5.1f %5.1f %6.2f %6.2f %6.2f', rF(k).speed, rF(k).snr, ...
        rF(k).miss_pct, rF(k).extra_pct, rF(k).stepAcc, rF(k).sdAcc, rF(k).slAcc, ...
        rP(k).miss_pct, rP(k).extra_pct, rP(k).stepAcc, rP(k).sdAcc, rP(k).slAcc); %#ok<AGROW>
end
fid = fopen(fullfile(outDir, 'SimVerify6_Report.txt'), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
logMsg(logFile, 'Finished in %.1f min. Results in %s', toc(tStart) / 60, outDir);
type(fullfile(outDir, 'SimVerify6_Report.txt'));
end

function s = readSelection(file)
fid = fopen(file, 'r');
fgetl(fid);                                                  % header line
s = struct();
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    p = strsplit(ln, ',');
    valueText = strjoin(p(2:end), ',');
    v = str2double(valueText);
    if startsWith(valueText, '[')
        s.(p{1}) = sscanf(strrep(strrep(valueText, '[', ''), ']', ''), '%f').';
        continue;
    end
    if isnan(v) && ~strcmpi(p{2}, 'NaN')
        s.(p{1}) = p{2};                                     % text value, e.g. twopass
    else
        s.(p{1}) = v;
    end
end
fclose(fid);
end

function rows = runAll(cfg, configs, opts, simCopy, seed, logFile, tStart)
rows = struct([]);
Nw = 2^nextpow2(round(0.1 * 44100));                        % SNR defined at the 0.186 s window (8192 samples), as in Chapter 4
for s = 1:numel(cfg.Speeds)
    D = simulate(cfg, cfg.Speeds(s), simCopy);              % clean signal and ground truth, once per speed
    for g = 1:numel(cfg.SNRs)
        acc = cell(1, numel(configs));                       % per-setting pooled scores for this cell
        for k = 1:cfg.nTrials
            rng(seed + 100000 * g + 1000 * round(10 * D.v0) + k, 'twister');   % fixed seed, so the run can be repeated
            z = randn(1, numel(D.clean)) + 1i * randn(1, numel(D.clean));
            noisePower = 10^(-(cfg.SNRs(g) - 10 * log10(Nw)) / 10);           % same noise equation as Chapter 4
            x = real(D.clean + sqrt(noisePower / 2) * z);    % one channel, like the IPM-365 recordings
            R = RealData_Pipeline_v14(x, D.fs, D.fc, D.c, configs, opts);
            for c = 1:numel(configs)
                acc{c} = addTrial(acc{c}, scoreTrial(R(c), D));
            end
        end
        for c = 1:numel(configs)
            row = summariseCell(acc{c});
            row.config = configs(c).name;
            row.speed = D.v0;
            row.snr = cfg.SNRs(g);
            if isempty(rows)
                rows = row;
            else
                rows(end + 1) = row; %#ok<AGROW>
            end
        end
        logMsg(logFile, '  %.1f m/s, %d dB done (%.1f min so far)', D.v0, cfg.SNRs(g), toc(tStart) / 60);
    end
end
end

function D = simulate(m_cfg, v0, m_simCopy)
SimulatorVerbose = false;                                    %#ok<NASGU>
SimulatorMakePlots = false;                                  %#ok<NASGU>
SimulatorR0_m = m_cfg.R0_m;                                  %#ok<NASGU>
AmplitudeMode = m_cfg.AmplitudeMode;                         %#ok<NASGU>
SNR_dB = 40;                                                 %#ok<NASGU>
CPI_s = 0.1;                                                 %#ok<NASGU>
rng(m_cfg.Seed + round(v0 * 100), 'twister');
evalc('run(m_simCopy);');
D.v0 = v0;
D.clean = RdrRxSignal_clean(:).';
D.fs = fs;
D.fc = fc;
D.c = c;
D.Tc = Tc;
D.t = t;
D.Torso_R = Torso_R;
D.buffer = bufferTime_s;
D.simT = SimulationTime_s;
D.truth = trueStepTimes_s_corrected(:).';                    % toe-speed peaks of the clean toe-only signal
end

function sc = scoreTrial(r, D)
edge = 0.1 + 0.1858 / 2;                                     % same walk interior for every setting
lo = D.buffer + edge;
hi = D.buffer + D.simT - edge;
t = r.stepTimes_s;
inside = t > lo & t < hi;
tIn = t(inside);
trIn = D.truth(D.truth > lo & D.truth < hi);
map = matchNew(tIn, trIn, min(0.2, D.Tc / 4));               % one-to-one matching in time order
sc.nTruth = numel(trIn);
sc.nMiss = numel(trIn) - nnz(map > 0);
sc.nExtra = nnz(map == 0);
sc.perfect = double(sc.nMiss == 0 && sc.nExtra == 0);
sc.stepErrPct = [];
if numel(tIn) >= 2
    ok = map(1:end-1) > 0 & diff(map) == 1;
    ii = find(ok);
    trueGap = trIn(map(ii) + 1) - trIn(map(ii));
    sc.stepErrPct = 100 * ((tIn(ii + 1) - tIn(ii)) - trueGap) ./ trueGap;
end
sc.sdErrPct = [];
if numel(tIn) >= 3
    sc.sdErrPct = 100 * ((tIn(3:end) - tIn(1:end-2)) - D.Tc) / D.Tc;
end
sc.slErrPct = [];
sl = r.strideLengths_m;
if numel(t) >= 3 && ~isempty(sl)
    trueTimes = D.buffer : D.Tc / 2 : D.buffer + D.simT;
    trueR = interp1(D.t, D.Torso_R, trueTimes, 'linear');
    trueLen = abs(trueR(3:end) - trueR(1:end-2));
    for k = 1:numel(sl)
        if t(k) > lo && t(k + 2) < hi
            [d, j] = min(abs(trueTimes(1:end-2) - t(k)));
            if d <= D.Tc / 4
                sc.slErrPct(end + 1) = 100 * (sl(k) - trueLen(j)) / trueLen(j);
            end
        end
    end
end
sc.wsErrPct = 100 * (r.estSpeed_mps - D.v0) / D.v0;
sc.tcErrPct = 100 * (r.Tc_used_s - D.Tc) / D.Tc;
end

function a = addTrial(a, sc)
if isempty(a)
    a = struct('nTruth', 0, 'nMiss', 0, 'nExtra', 0, 'perfect', [], 'step', [], 'sd', [], 'sl', [], 'ws', [], 'tc', []);
end
a.nTruth = a.nTruth + sc.nTruth;
a.nMiss = a.nMiss + sc.nMiss;
a.nExtra = a.nExtra + sc.nExtra;
a.perfect(end + 1) = sc.perfect;
a.step = [a.step, sc.stepErrPct];
a.sd = [a.sd, sc.sdErrPct];
a.sl = [a.sl, sc.slErrPct];
a.ws(end + 1) = sc.wsErrPct;
a.tc(end + 1) = sc.tcErrPct;
end

function row = summariseCell(a)
row.config = '';
row.speed = NaN;
row.snr = NaN;
row.miss_pct = 100 * a.nMiss / max(1, a.nTruth);
row.extra_pct = 100 * a.nExtra / max(1, a.nTruth);
row.perfect_pct = 100 * mean(a.perfect);
row.stepAcc = 100 - rmsOf(a.step);                           % 100 - RMSE%, Chapter 4 definition
row.sdAcc = 100 - rmsOf(a.sd);
row.slAcc = 100 - rmsOf(a.sl);
row.wsAcc = 100 - rmsOf(a.ws);
row.tcErr_pct = mean(a.tc);
row.nTruth = a.nTruth;
row.nMiss = a.nMiss;
row.nExtra = a.nExtra;
end

function v = rmsOf(x)
x = x(isfinite(x));
if isempty(x)
    v = 100;                                                 % nothing to score counts as 0% accuracy
else
    v = sqrt(mean(x.^2));
end
end

function [nPass, me] = cellScore(rows, name, snrs, cfg)
m = strcmp({rows.config}, name) & ismember([rows.snr], snrs);
r = rows(m);
nPass = sum([r.stepAcc] >= cfg.PassAcc_pct & [r.sdAcc] >= cfg.PassAcc_pct & [r.slAcc] >= cfg.PassAcc_pct);
me = 100 * (sum([r.nMiss]) + sum([r.nExtra])) / max(1, sum([r.nTruth]));
end

function s = yesNo(tf, a, b)
if tf
    s = a;
else
    s = b;
end
end

function writeCells(file, rows)
fid = fopen(file, 'w');
fprintf(fid, 'config,speed_mps,snr_dB,miss_pct,extra_pct,perfect_pct,stepAcc_pct,strideDurAcc_pct,strideLenAcc_pct,speedAcc_pct,TcErr_pct,nTrueSteps\n');
for k = 1:numel(rows)
    r = rows(k);
    fprintf(fid, '%s,%.1f,%d,%.2f,%.2f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n', r.config, r.speed, r.snr, r.miss_pct, r.extra_pct, ...
        r.perfect_pct, r.stepAcc, r.sdAcc, r.slAcc, r.wsAcc, r.tcErr_pct, r.nTruth);
end
fclose(fid);
end

function s = summaryLine(rows, name, cfg, snrs)
r = rows(strcmp({rows.config}, name) & ismember([rows.snr], snrs));
p = sum([r.stepAcc] >= cfg.PassAcc_pct & [r.sdAcc] >= cfg.PassAcc_pct & [r.slAcc] >= cfg.PassAcc_pct);
miss = 100 * sum([r.nMiss]) / max(1, sum([r.nTruth]));
extra = 100 * sum([r.nExtra]) / max(1, sum([r.nTruth]));
s = sprintf('  %-8s %5d/%-2d %8.2f %8.2f %10.1f %9.2f %9.2f %9.2f', name, p, numel(r), miss, extra, mean([r.perfect_pct]), ...
    mean([r.stepAcc]), mean([r.sdAcc]), mean([r.slAcc]));
end

function dest = makeSimCopy(simName, outDir)
src = [simName '.m'];
assert(exist(src, 'file') == 2, 'Simulator %s not found in the Current Folder.', src);
txt = fileread(src);
txt = regexprep(txt, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*(?:\r?\n|$)', '');   % no clear all inside a function
dest = fullfile(outDir, 'WC_Simulator.m');
fid = fopen(dest, 'w');
fprintf(fid, '%s', txt);
fclose(fid);
end

function logMsg(logFile, fmt, varargin)
msg = sprintf(fmt, varargin{:});
line = sprintf('[%s] %s', datestr(now, 'HH:MM:SS'), msg);
fprintf('%s\n', line);
fid = fopen(logFile, 'a');
if fid > 0
    fprintf(fid, '%s\n', line);
    fclose(fid);
end
end

function map = matchNew(det, truth, tolerance)
det = det(:).';
truth = truth(:).';
n = numel(det);
m = numel(truth);
map = zeros(1, n);
if n == 0 || m == 0
    return;
end
count = zeros(n + 1, m + 1);
cost = zeros(n + 1, m + 1);
action = zeros(n + 1, m + 1);
action(2:end, 1) = 1;
action(1, 2:end) = 2;
for i = 1:n
    for j = 1:m
        counts = [count(i, j+1), count(i+1, j)];
        costs = [cost(i, j+1), cost(i+1, j)];
        actions = [1, 2];
        distance = abs(det(i) - truth(j));
        if distance <= tolerance
            counts(end+1) = count(i, j) + 1; %#ok<AGROW>
            costs(end+1) = cost(i, j) + distance; %#ok<AGROW>
            actions(end+1) = 3; %#ok<AGROW>
        end
        keys = [-counts(:), costs(:), -actions(:)];
        [~, order] = sortrows(keys, [1, 2, 3]);
        best = order(1);
        count(i+1, j+1) = counts(best);
        cost(i+1, j+1) = costs(best);
        action(i+1, j+1) = actions(best);
    end
end
i = n;
j = m;
while i > 0 || j > 0
    a = action(i+1, j+1);
    if a == 3
        map(i) = j;
        i = i - 1;
        j = j - 1;
    elseif a == 1
        i = i - 1;
    elseif a == 2
        j = j - 1;
    else
        break;
    end
end
end
