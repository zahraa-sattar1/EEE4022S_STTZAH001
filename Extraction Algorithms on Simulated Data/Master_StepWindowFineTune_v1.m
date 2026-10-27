%% Step-window comparison
% Zahraa Sattar | EEE4022S
%
% Compare exact step windows with the power-of-two baseline.
% Keep noise draws paired when scoring step timing and stride estimates.

function Master_StepWindowFineTune_v1(userCfg)

if nargin < 1
    userCfg = struct();                                      % no overrides given
end
cfg = defaultConfig();                                       % default settings
if isfield(userCfg, 'ResumeFolder') && ~isempty(userCfg.ResumeFolder)
    cfgFile = fullfile(userCfg.ResumeFolder, 'Config.mat');
    if isfile(cfgFile)
        saved = load(cfgFile);                               % keep the settings the run started with
        cfg = saved.cfg;
    end
end
names = fieldnames(userCfg);
for k = 1:numel(names)
    cfg.(names{k}) = userCfg.(names{k});                     % apply overrides
end
if cfg.Quick                                                 % quick test, unless the caller set these themselves
    if ~isfield(userCfg, 'Speeds'), cfg.Speeds = [0.6 1.6]; end       % one slow speed, and the speed with the missed step
    if ~isfield(userCfg, 'nTrials'), cfg.nTrials = 3; end             % 3 noise draws per setting
end

codeFolder = cfg.CodeFolder;
if isempty(codeFolder)
    codeFolder = fileparts(mfilename('fullpath'));           % folder this script sits in
end
if isempty(codeFolder)
    codeFolder = pwd;
end
oldPath = path;
addpath(codeFolder);
pathCleanup = onCleanup(@() path(oldPath));                  % put the path back when finished

if isempty(cfg.ResumeFolder)
    outDir = fullfile(codeFolder, ['StepWindowFT_' datestr(now, 'yyyymmdd_HHMMSS')]);   % always a new folder
    mkdir(outDir);
else
    outDir = cfg.ResumeFolder;
    if ~exist(outDir, 'dir')
        outDir = fullfile(codeFolder, cfg.ResumeFolder);
    end
end
cfg.ResumeFolder = outDir;
cfg.CodeFolder = codeFolder;
if ~isfile(fullfile(outDir, 'Config.mat'))
    save(fullfile(outDir, 'Config.mat'), 'cfg');             % settings of this run, used on resume
end

global VX
VX = defaultVX();                                            % option block read by the patched copies

tRun = tic;                                                  % total run time
logMsg(outDir, 'Output folder: %s', outDir);
logMsg(outDir, 'Software: %s', softwareName());
logMsg(outDir, 'Seeds: A = %d, B = %d (trial seed = base + 7919*trial + 1000*speed, same draw for every window)', cfg.SeedA, cfg.SeedB);
logMsg(outDir, 'Step windows %s s (%.3f s = power-of-two baseline), walking speed window 0.372 s (requested %g s)', ...
    mat2str(cfg.StepWindows), cfg.BaselineWindow, cfg.WSCPIreq);
logMsg(outDir, 'Speeds %s, SNRs %s dB, trials %d, extra miss check %s m/s at %s dB, quick mode %d', ...
    mat2str(cfg.Speeds), mat2str(cfg.SNRs), cfg.nTrials, mat2str(cfg.MissCheckSpeeds), mat2str(cfg.MissCheckSNRs), cfg.Quick);

P = makeWorkingCopies(cfg, codeFolder, outDir);
logMsg(outDir, 'Working copies made. Patch log: Patch_Log.txt');
selfTest(P, cfg, outDir);                                    % stops the run if a patch did not take effect

tStart = tic;
for sIdx = 1:numel(cfg.Speeds)
    v0 = cfg.Speeds(sIdx);
    f = fullfile(outDir, sprintf('FT_v0_%.2f.mat', v0));     % results of one speed, so a crash loses at most one speed
    if isfile(f)
        logMsg(outDir, 'v0 = %.2f m/s already finished, skipping.', v0);
        continue;
    end
    logMsg(outDir, '===== v0 = %.2f m/s (%d of %d), elapsed %.1f min =====', v0, sIdx, numel(cfg.Speeds), toc(tStart)/60);
    try
        D = getSpeedData(P, cfg, v0);                        % one simulator run gives the clean signal for every trial
        snrs = cfg.SNRs;
        if any(abs(cfg.MissCheckSpeeds - v0) < 1e-9)
            snrs = [snrs, cfg.MissCheckSNRs];                % the extra 70 dB check for 1.6 and 2.0 m/s
        end
        rows = [];
        for snr = snrs
            rows = [rows, runCondition(P, cfg, D, snr)];     %#ok<AGROW>
            logMsg(outDir, '   %g dB done, elapsed %.1f min', snr, toc(tStart)/60);
        end
        save(f, 'rows');
    catch ME
        logMsg(outDir, 'Speed %.2f stopped early: %s', v0, ME.message);
        for q = 1:min(4, numel(ME.stack))
            logMsg(outDir, '   at %s line %d', ME.stack(q).name, ME.stack(q).line);
        end
    end
end

logMsg(outDir, 'All speeds done in %.1f min. Writing results.', toc(tStart)/60);
try
    assembleResults(cfg, outDir);
catch ME
    logMsg(outDir, 'Writing the results failed: %s', ME.message);
    for q = 1:min(5, numel(ME.stack))
        logMsg(outDir, '   at %s line %d', ME.stack(q).name, ME.stack(q).line);
    end
end
logMsg(outDir, 'Finished in %.1f min. Open FineTune_Report.txt in %s', toc(tRun)/60, outDir);
end

function cfg = defaultConfig()
cfg.Speeds = [0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];      % m/s
cfg.StepWindows = [0.12 0.14 0.16 0.186 0.21 0.25 0.30];     % s, step detection windows to test
cfg.BaselineWindow = 0.186;                                  % s, the power-of-two window used so far
cfg.BaselineCPIreq = 0.1;                                    % s, request that gives the baseline (N = 8192, 0.18576 s)
cfg.WSCPIreq = 0.2;                                          % s, walking speed request for stride length (power of two, 0.372 s)
cfg.SNRs = [30 40 60 80];                                    % dB
cfg.MissCheckSpeeds = [1.6 2.0];                             % m/s, speeds that missed a step in the final sweep
cfg.MissCheckSNRs = 70;                                      % dB, extra SNR for those speeds only (not counted in the passes)
cfg.nTrials = 30;                                            % noise draws per setting
cfg.Quick = false;                                           % true = 2 speeds, 3 trials
cfg.SeedA = 20261001;                                        % seed base of the main run (same as Master_SimVariants_v1)
cfg.SeedB = 20262002;                                        % seed base of the second baseline run (the noise estimate)
cfg.R0_m = 10;                                               % simulator range
cfg.AmplitudeMode = 'A';                                     % Chen mean sqrt(RCS) weights
cfg.Target_pct = 5;                                          % goal for RMSE %
cfg.RMSEFloor_pp = 0.05;                                     % smallest RMSE noise trusted (same floor as Master_SimVariants_v1)
cfg.StepRefine = 'centroid';                                 % v5 extractor default
cfg.StepGapFill = 1;                                         % v5 extractor default
cfg.MakePlots = true;
cfg.SimName = 'RadarSimulator_v11_torso_legs_arms_head_buffered_noplots';
cfg.WSName = 'Finding_WalkingSpeed_v4_noplots';
cfg.STName = 'Finding_StepRecognition_v5_noplots';
cfg.CodeFolder = '';
cfg.ResumeFolder = '';
cfg.CopyHook = [];                                           % only used for testing outside MATLAB
end

function VX = defaultVX()
VX = struct('WSexact', false, 'WSbandA', 0, 'WSpasses', 3, 'WSHop', 0.25, ...
    'STexact', false);                                       % STexact = step window not forced to a power of two
end

function N = exactN(cpi, fs)
N = 2*round(cpi*fs/2);
end

function Nw = noiseWindow(cpi, fs, exact)
if exact
    Nw = exactN(cpi, fs);                                    % same N as the patched step spectrogram
else
    Nw = 2^nextpow2(round(cpi*fs));                          % power of two, as before
end
end

function cpiAct = actualCPI(cpi, fs, exact)
if nargin < 3
    exact = false;
end
cpiAct = noiseWindow(cpi, fs, exact)/fs;
end

function w = windowSpec(cfg, win)
if abs(win - cfg.BaselineWindow) < 1e-9
    w = struct('win', win, 'req', cfg.BaselineCPIreq, 'exact', false);   % 0.1 s request -> N = 8192 -> 0.18576 s
else
    w = struct('win', win, 'req', win, 'exact', true);                   % exact window, N = even round(win*fs)
end
end

function P = makeWorkingCopies(cfg, codeFolder, outDir)
codeDir = fullfile(outDir, 'WorkingCopies');
if ~exist(codeDir, 'dir')
    mkdir(codeDir);
end
logFile = fullfile(outDir, 'Patch_Log.txt');
P.sim = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.SimName '.m']), fullfile(codeDir, 'WC_Simulator.m'), 'sim', cfg, logFile);
P.ws = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.WSName '.m']), fullfile(codeDir, 'WC_WalkingSpeed.m'), 'ws', cfg, logFile);
P.st = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.STName '.m']), fullfile(codeDir, 'WC_StepRecognition_v5.m'), 'st', cfg, logFile);
end

function dest = copyOne(src, dest, kind, cfg, logFile)
if isfile(dest)
    return;                                                  % never overwrite, reuse on resume
end
assert(isfile(src), 'Missing file: %s', src);
txt = fileread(src);
tag = [kind ':' src];
pat = '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*';
m = regexp(txt, pat, 'match');                               % lines about to be removed, for the log
for q = 1:numel(m)
    logPatch(logFile, tag, 'removed clear all / close all', m{q}, '(line removed)');
end
txt = regexprep(txt, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*(?:\r?\n|$)', '');   % no clear all / close all
if ~strcmp(kind, 'sim')
    pat = '(?m)^[ \t]*RadarSimulator_[A-Za-z0-9_]+;[^\r\n]*';
    m = regexp(txt, pat, 'match');                           % simulator call about to be removed, for the log
    for q = 1:numel(m)
        logPatch(logFile, tag, 'removed the simulator call', m{q}, '(line removed)');
    end
    txt = regexprep(txt, '(?m)^[ \t]*RadarSimulator_[A-Za-z0-9_]+;[^\r\n]*(?:\r?\n|$)', '');   % simulator is called separately
end
if strcmp(kind, 'st')
    plotMarker = 'figure(''Name'',''Step 6 verification - remove outliers'');';
    a = strfind(txt, plotMarker);
    b = strfind(txt, '%% Step 7');
    assert(numel(a) == 1 && numel(b) == 1 && b > a, 'Step 6 figure not found in %s.', src);
    txt = [txt(1:a-1), txt(b:end)];
    logPatch(logFile, tag, 'cut the Step 6 figure', plotMarker, '(removed up to %% Step 7)');
    txt = patchRegex(txt, '(?m)^[ \t]*freq_high_Yardibi[ \t]*=[ \t]*[0-9.]+\*24/5\.8;[^\r\n]*', ...
        'freq_high_Yardibi = 350*24/5.8;   % this copy: 350 Hz limit', 1, 'extractor frequency limit 350 Hz', tag, logFile);
    txt = patchRegex(txt, ['N = 2\^\(nextpow2\(round\(CPI_s\*fsProc\)\)\);[^\r\n]*\r?\n[ \t]*win[ \t]*=[ \t]*hamming\(N\);[^\r\n]*\r?\n' ...
        '[ \t]*noverlap = N\*3/4;[^\r\n]*\r?\n[ \t]*nfft[ \t]*=[ \t]*N;[^\r\n]*'], ...
        sprintf(['global VX\n' ...
        '        N = 2^(nextpow2(round(CPI_s*fsProc)));   %% original: window rounded up to a power of two\n' ...
        '        noverlap = N*3/4;                       %% original: 75 percent overlap\n' ...
        '        if VX.STexact\n' ...
        '            N = 2*round(CPI_s*fsProc/2);        %% exact window, not forced to a power of two (even, so the frequency axis stays exact)\n' ...
        '            noverlap = round(N*3/4);            %% 75 percent overlap in whole samples\n' ...
        '        end\n' ...
        '        win      = hamming(N);\n' ...
        '        nfft     = N;                           %% no zero padding, so the Doppler resolution is 1/window']), ...
        1, 'step spectrogram block (exact window switch)', tag, logFile);
end
if strcmp(kind, 'ws')
    txt = patchRegex(txt, 'N_actual = 2\^\(nextpow2\(round\(CPI_s\*fs\)\)\);[^\r\n]*', ...
        sprintf(['global VX\nN_actual = 2^(nextpow2(round(CPI_s*fs)));\nif VX.WSexact, N_actual = round(CPI_s*fs); end   %% window not forced to a power of two']), ...
        1, 'actual CPI line', tag, logFile);
    txt = patchRegex(txt, 'numIterations\s*=\s*3;[^\r\n]*', ...
        sprintf(['global VX\nnumIterations = VX.WSpasses;   %% number of passes (3 = original)']), ...
        1, 'number of passes', tag, logFile);
    txt = patchRegex(txt, ['N = 2\^\(nextpow2\(round\(CPI_s\*fsProc\)\)\);[^\r\n]*\r?\n[ \t]*win[ \t]*=[ \t]*hamming\(N\);[^\r\n]*\r?\n' ...
        '[ \t]*noverlap = N\*3/4;[^\r\n]*\r?\n[ \t]*nfft[ \t]*=[ \t]*N;[^\r\n]*'], ...
        sprintf(['global VX\nN = 2^(nextpow2(round(CPI_s*fsProc)));\nif VX.WSexact, N = round(CPI_s*fsProc); end   %% window not forced to a power of two\n' ...
        'win = hamming(N);\nnoverlap = N - max(1, round(N*VX.WSHop));   %% hop (0.25 = the original 75 percent overlap)\nnfft = 2^(nextpow2(N));   %% zero padded when N is not a power of two']), ...
        1, 'walking speed spectrogram block', tag, logFile);
end
if ~isempty(cfg.CopyHook)
    txt = cfg.CopyHook(txt, kind);                           % only used when testing outside MATLAB
end
fid = fopen(dest, 'w');
assert(fid >= 0, 'Could not create %s', dest);
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function m_o = runStep(m_P, m_x, m_D, m_w, m_snr)
global VX
VX = defaultVX();
VX.STexact = m_w.exact;                                      % exact window or the original power of two
m_o = struct('ok', false, 'det', [], 'key', [], 'sd', [], 'cpiAct', NaN, 'sc', [], 'framePeriod', NaN, 'nfft', NaN, 'err', '');
try
    RdrRxSignal = m_x;                                       %#ok<NASGU>
    fs = m_D.fs;                                             %#ok<NASGU>
    Tc = m_D.Tc;                                             %#ok<NASGU>
    CPI_s = m_w.req;                                         %#ok<NASGU>
    SNR_dB = m_snr;                                          %#ok<NASGU>
    stepRefineMode = 'centroid';                             %#ok<NASGU>
    stepGapFill = 1;                                         %#ok<NASGU>
    evalc('run(m_P.st);');
    m_o.det = detectedStepTimes_s(:).';                      % final step times
    m_o.key = stepFrameTimes_s(:).';                         % frame time stamps decide who is inside the walk
    m_o.cpiAct = CPI_actual_s;
    m_o.framePeriod = frame_period_s;
    m_o.nfft = numel(freq_axis_step_full);                   % FFT length really used
    edge = 0.1 + m_o.cpiAct/2;                               % same walk interior as the step time scoring
    tIn = m_o.det(m_o.key > m_D.buffer + edge & m_o.key < m_D.buffer + m_D.simT - edge);   % steps inside the walk only
    if numel(tIn) >= 3
        m_o.sd = (tIn(3:end) - tIn(1:end-2)) - m_D.Tc;       % stride durations, kept as (stride - Tc)
    end
    m_o.sc = scoreStepTimes(m_o.det, m_D, m_o.cpiAct, m_o.key);
    m_o.ok = true;
catch ME
    m_o.err = ME.message;
end
VX = defaultVX();
end

function rows = runCondition(P, cfg, D, snr)
nT = cfg.nTrials;
wins = cfg.StepWindows;
specs = cell(1, numel(wins));
for i = 1:numel(wins)
    specs{i} = windowSpec(cfg, wins(i));
end
ws = struct('cpi', cfg.WSCPIreq, 'exact', false, 'bandA', 0, 'passes', 3, 'hop', 0.25);   % walking speed, original settings
iBase = find(abs(wins - cfg.BaselineWindow) < 1e-9, 1);
rows = [];
for tag = {'A', 'B'}
    if strcmp(tag{1}, 'A')
        seed = cfg.SeedA;
        idx = 1:numel(wins);                                 % every window on seed A
    else
        seed = cfg.SeedB;
        idx = iBase;                                         % only the baseline on seed B
    end
    W = cell(nT, 1);
    S = cell(nT, numel(wins));
    for k = 1:nT
        z = trialNoise(D, seed, k);                          % the same draw for every window
        W{k} = runSpeed(P, addNoise(D.clean, noiseWindow(ws.cpi, D.fs, false), snr, z), D, ws);
        for i = idx
            x = addNoise(D.clean, noiseWindow(specs{i}.req, D.fs, specs{i}.exact), snr, z);   % SNR set on this window's N
            S{k, i} = runStep(P, x, D, specs{i}, snr);
        end
    end
    for i = idx
        for p = {'StepTime', 'StrideDuration', 'StrideLength'}
            rows = addRow(rows, makeRow(cfg, D, W, S(:, i), specs{i}, snr, p{1}, tag{1}, seed));
        end
    end
end
end

function row = makeRow(cfg, D, W, S, w, snr, param, tag, seed)
nT = numel(S);
E = cell(nT, 1);
G = cell(nT, 1);
miss = NaN(nT, 1);
extra = NaN(nT, 1);
perfect = NaN(nT, 1);
cpiV = NaN(nT, 1);
for k = 1:nT
    E{k} = NaN;
    G{k} = NaN;
    if ~S{k}.ok
        continue;
    end
    switch param
        case 'StepTime'
            E{k} = S{k}.sc.e;
            G{k} = S{k}.sc.g;
        case 'StrideDuration'
            E{k} = S{k}.sd;
            G{k} = D.Tc * ones(size(S{k}.sd));
        case 'StrideLength'
            [E{k}, G{k}] = strideLengthErrors(W{k}, S{k}, D);
    end
    miss(k) = S{k}.sc.miss;                                  % the step list behind all three parameters
    extra(k) = S{k}.sc.extra;
    perfect(k) = S{k}.sc.perfect;
    cpiV(k) = S{k}.cpiAct;
end
row = blankRow();
row.Parameter = param;
row.Speed_mps = D.v0;
row.SNR_dB = snr;
row.SeedTag = tag;
row.Seed = seed;
row.Window_s = w.win;
row.Window_req_s = w.req;
row.Window_act_s = actualCPI(w.req, D.fs, w.exact);
row.Window_samples = noiseWindow(w.req, D.fs, w.exact);
row.ExactWindow = double(w.exact);
row.IsBaseline = double(abs(w.win - cfg.BaselineWindow) < 1e-9);
row.ExtraCheck = double(~any(abs(cfg.SNRs - snr) < 1e-9));   % 1 = the 70 dB miss check, not counted in the passes
row.CPI_speed_act = actualCPI(cfg.WSCPIreq, D.fs, false);
row.CPI_used_check_s = meanFinite(cpiV);                     % what the extractor reported, should equal Window_act_s
row = fillStats(row, E, G, strcmp(param, 'StrideLength'), cfg);
row.MissRate_pct = meanFinite(miss);
row.ExtraRate_pct = meanFinite(extra);
row.PerfectTrials_pct = 100*meanFinite(perfect);
end

function row = blankRow()
row = struct('Parameter', '', 'Speed_mps', NaN, 'SNR_dB', NaN, 'SeedTag', '', 'Seed', NaN, ...
    'Window_s', NaN, 'Window_req_s', NaN, 'Window_act_s', NaN, 'Window_samples', NaN, 'ExactWindow', NaN, 'IsBaseline', NaN, 'ExtraCheck', NaN, ...
    'CPI_speed_act', NaN, 'CPI_used_check_s', NaN, 'Trials', NaN, 'FailedTrials', NaN, 'N_compared', NaN, 'TruthMean', NaN, ...
    'RMSE_abs', NaN, 'RMSE_pct', NaN, 'Accuracy_pct', NaN, 'MPE_pct', NaN, 'MeanErr_pct', NaN, 'TrialMean_pct', NaN, 'TrialStd_pct', NaN, ...
    'MissRate_pct', NaN, 'ExtraRate_pct', NaN, 'PerfectTrials_pct', NaN, 'Pass', NaN);
end

function row = fillStats(row, E, G, isSL, cfg)
nT = numel(E);
trialRMSE = NaN(nT, 1);
trialPct = NaN(nT, 1);
allE = [];
allG = [];
for k = 1:nT
    e = E{k}(:).';
    g = G{k}(:).';
    ok = isfinite(e) & isfinite(g) & g > 0;
    e = e(ok);
    g = g(ok);
    if ~isempty(e)
        trialRMSE(k) = rmsOf(e);
        trialPct(k) = 100*trialRMSE(k)/mean(g);
        allE = [allE, e];                                    %#ok<AGROW>
        allG = [allG, g];                                    %#ok<AGROW>
    end
end
row.Trials = nT;
row.FailedTrials = nnz(~isfinite(trialRMSE));
row.N_compared = numel(allE);
if isempty(allE)
    return;
end
row.TruthMean = mean(allG);
pooled = rmsOf(allE);
if isSL
    row.RMSE_abs = meanFinite(trialRMSE);                    % stride length averages the trial RMSE, as in the sweep
else
    row.RMSE_abs = pooled;
end
row.RMSE_pct = 100*row.RMSE_abs/row.TruthMean;
row.Accuracy_pct = 100 - row.RMSE_pct;                       % accuracy as asked for this test: 100 minus RMSE %
row.MPE_pct = 100*mean(allE ./ allG);                        % mean percentage error, signed, so it shows bias
row.MeanErr_pct = 100*mean(allE)/row.TruthMean;              % bias as the sweep defines it
row.TrialMean_pct = meanFinite(trialPct);
row.TrialStd_pct = stdFinite(trialPct);                      % spread of the per-trial RMSE % across the trials
row.Pass = double(row.RMSE_pct < cfg.Target_pct);
end

function selfTest(P, cfg, outDir)
global VX
res = {};
ok = true;
D = getSpeedData(P, cfg, 1.2);
setSeed(1);
z = randn(1, numel(D.clean)) + 1i*randn(1, numel(D.clean));
b = windowSpec(cfg, cfg.BaselineWindow);
o0 = runStep(P, addNoise(D.clean, noiseWindow(b.req, D.fs, b.exact), 60, z), D, b, 60);
[ok, res] = chk(ok, res, o0.ok && o0.nfft == 8192 && abs(o0.cpiAct - 8192/D.fs) < 1e-12 && abs(o0.framePeriod - 2048/D.fs) < 1e-12, ...
    sprintf('baseline window: N = %d, window %.5f s, frame period %.5f s (expected 8192, 0.18576 s, 0.04644 s)', o0.nfft, o0.cpiAct, o0.framePeriod), o0.err);
for win = [0.21 0.12]
    w = windowSpec(cfg, win);
    Nexp = exactN(win, D.fs);
    o1 = runStep(P, addNoise(D.clean, noiseWindow(w.req, D.fs, w.exact), 60, z), D, w, 60);
    hopExp = Nexp - round(Nexp*3/4);
    [ok, res] = chk(ok, res, o1.ok && o1.nfft == Nexp && abs(o1.cpiAct - Nexp/D.fs) < 1e-12 && abs(o1.framePeriod - hopExp/D.fs) < 1e-12 ...
        && noiseWindow(w.req, D.fs, true) == Nexp && abs(actualCPI(w.req, D.fs, true) - o1.cpiAct) < 1e-12, ...
        sprintf('exact window %.3f s: N = %d (expected %d), window %.5f s, hop %d samples, noise window and actualCPI agree', win, o1.nfft, Nexp, o1.cpiAct, hopExp), o1.err);
end
ws = struct('cpi', cfg.WSCPIreq, 'exact', false, 'bandA', 0, 'passes', 3, 'hop', 0.25);
w0 = runSpeed(P, addNoise(D.clean, noiseWindow(ws.cpi, D.fs, false), 40, z), D, ws);
[ok, res] = chk(ok, res, w0.ok && abs(w0.cpiAct - 16384/D.fs) < 1e-12, sprintf('walking speed window %.4f s (expected 0.3715 s), speed %.3f m/s', w0.cpiAct, w0.meanSpeed), w0.err);
VX = defaultVX();
fid = fopen(fullfile(outDir, 'Self_Test.txt'), 'w');
for k = 1:numel(res)
    fprintf(fid, '%s\n', res{k});
    logMsg(outDir, 'Self test: %s', res{k});
end
fclose(fid);
if ~ok
    error('Self test failed, see Self_Test.txt in %s. Nothing was run.', outDir);
end
end

function assembleResults(cfg, outDir)
ALL = [];
for v0 = cfg.Speeds
    f = fullfile(outDir, sprintf('FT_v0_%.2f.mat', v0));
    if isfile(f)
        s = load(f);
        ALL = [ALL, s.rows];                                 %#ok<AGROW>
    end
end
if isempty(ALL)
    logMsg(outDir, 'No results to write.');
    return;
end
writeRowsCSV(fullfile(outDir, 'FineTune_Results.csv'), ALL);   % every row, seed A and seed B
C = buildComparison(cfg, ALL);
writeRowsCSV(fullfile(outDir, 'FineTune_vs_Baseline.csv'), C);
PS = buildPassSummary(cfg, ALL, C);
writeRowsCSV(fullfile(outDir, 'Pass_Summary.csv'), PS);
if cfg.MakePlots
    try
        makePlots(cfg, outDir, ALL, PS);
    catch ME
        logMsg(outDir, 'Plots skipped: %s', ME.message);
    end
end
writeReport(cfg, outDir, ALL, C, PS);
end

function r = pickRow(ALL, param, v0, snr, win, tag)
r = ALL(strcmp({ALL.Parameter}, param) & abs([ALL.Speed_mps] - v0) < 1e-9 & abs([ALL.SNR_dB] - snr) < 1e-9 ...
    & abs([ALL.Window_s] - win) < 1e-9 & strcmp({ALL.SeedTag}, tag));
if numel(r) > 1
    r = r(1);
end
end

function C = buildComparison(cfg, ALL)
C = [];
A = ALL(strcmp({ALL.SeedTag}, 'A'));
for i = 1:numel(A)
    r = A(i);
    base = pickRow(ALL, r.Parameter, r.Speed_mps, r.SNR_dB, cfg.BaselineWindow, 'A');
    baseB = pickRow(ALL, r.Parameter, r.Speed_mps, r.SNR_dB, cfg.BaselineWindow, 'B');
    if isempty(base) || isempty(baseB)
        continue;
    end
    c = struct();
    c.Parameter = r.Parameter;
    c.Speed_mps = r.Speed_mps;
    c.SNR_dB = r.SNR_dB;
    c.ExtraCheck = r.ExtraCheck;
    c.Window_s = r.Window_s;
    c.Window_act_s = r.Window_act_s;
    c.RMSE_pct = r.RMSE_pct;
    c.RMSE_base = base.RMSE_pct;
    c.dRMSE_pp = r.RMSE_pct - base.RMSE_pct;                 % negative = better than the baseline
    c.NoiseRMSE_pp = max(abs(base.RMSE_pct - baseB.RMSE_pct), cfg.RMSEFloor_pp);   % seed A minus seed B of the baseline, with a floor
    c.RMSE_change = direction(-c.dRMSE_pp, c.NoiseRMSE_pp);  % lower RMSE is better
    c.Exact_pct = r.PerfectTrials_pct;
    c.Exact_base = base.PerfectTrials_pct;
    c.dExact_pp = r.PerfectTrials_pct - base.PerfectTrials_pct;
    c.NoiseExact_pp = max(abs(base.PerfectTrials_pct - baseB.PerfectTrials_pct), 100/cfg.nTrials);   % at least one trial
    c.Exact_change = direction(c.dExact_pp, c.NoiseExact_pp);   % more exact trials is better
    c.Miss_pct = r.MissRate_pct;
    c.Miss_base = base.MissRate_pct;
    c.dMiss_pp = r.MissRate_pct - base.MissRate_pct;
    c.NoiseMiss_pp = abs(base.MissRate_pct - baseB.MissRate_pct);
    c.Miss_change = direction(-c.dMiss_pp, c.NoiseMiss_pp);  % fewer misses is better
    c.Extra_pct = r.ExtraRate_pct;
    c.MPE_pct = r.MPE_pct;
    c.TrialStd_pct = r.TrialStd_pct;
    c.Pass = r.Pass;
    c.PassBase = base.Pass;
    C = addRow(C, c);
end
end

function s = direction(gain, noise)
if ~isfinite(gain) || ~isfinite(noise)
    s = 'no value';
elseif abs(gain) <= noise + 1e-9                         % equal to the noise also counts as within noise
    s = 'within noise';
elseif gain > 0
    s = 'better';
else
    s = 'worse';
end
end

function PS = buildPassSummary(cfg, ALL, C)
PS = [];
params = {'StepTime', 'StrideDuration', 'StrideLength'};
for snr = cfg.SNRs
    for win = [cfg.StepWindows, NaN]                         % NaN = the seed B baseline, the noise reference
        tag = 'A';
        w = win;
        if isnan(win)
            tag = 'B';
            w = cfg.BaselineWindow;
        end
        R = ALL(strcmp({ALL.SeedTag}, tag) & abs([ALL.SNR_dB] - snr) < 1e-9 & abs([ALL.Window_s] - w) < 1e-9);
        if isempty(R)
            continue;
        end
        p = struct();
        p.Window_s = w;
        p.SeedTag = tag;
        p.SNR_dB = snr;
        p.Window_act_s = R(1).Window_act_s;
        p.Passes = sum([R.Pass] == 1);
        p.Conditions = numel(R);
        for q = 1:3
            Rq = R(strcmp({R.Parameter}, params{q}));
            p.([params{q} '_passes']) = sum([Rq.Pass] == 1);
            p.([params{q} '_meanRMSE']) = meanFinite([Rq.RMSE_pct]);
        end
        st = R(strcmp({R.Parameter}, 'StepTime'));
        p.MeanRMSE_pct = meanFinite([R.RMSE_pct]);           % mean over the 30 conditions, used to rank the windows
        p.MeanExact_pct = meanFinite([st.PerfectTrials_pct]); % step recognition, mean over the speeds
        p.MeanMiss_pct = meanFinite([st.MissRate_pct]);
        p.MeanExtra_pct = meanFinite([st.ExtraRate_pct]);
        if strcmp(tag, 'A') && ~isempty(C)
            k = abs([C.SNR_dB] - snr) < 1e-9 & abs([C.Window_s] - w) < 1e-9;
            p.RMSE_better = sum(strcmp({C(k).RMSE_change}, 'better'));   % conditions better than the baseline beyond noise
            p.RMSE_worse = sum(strcmp({C(k).RMSE_change}, 'worse'));
            p.Exact_better = sum(strcmp({C(k).Exact_change}, 'better') & strcmp({C(k).Parameter}, 'StepTime'));
            p.Exact_worse = sum(strcmp({C(k).Exact_change}, 'worse') & strcmp({C(k).Parameter}, 'StepTime'));
        else
            p.RMSE_better = NaN;
            p.RMSE_worse = NaN;
            p.Exact_better = NaN;
            p.Exact_worse = NaN;
        end
        PS = addRow(PS, p);
    end
end
end

function v = verdictPerSNR(cfg, PS, snr)
A = PS(strcmp({PS.SeedTag}, 'A') & [PS.SNR_dB] == snr);
B = PS(strcmp({PS.SeedTag}, 'B') & [PS.SNR_dB] == snr);
v = struct('ok', false);
if isempty(A) || isempty(B)
    return;
end
[~, o] = sort([A.Window_s]);
A = A(o);
base = A(abs([A.Window_s] - cfg.BaselineWindow) < 1e-9);
v.passNoise = max(1, abs(base.Passes - B.Passes));           % passes can move by this much between seeds
v.meanNoise = max(cfg.RMSEFloor_pp, abs(base.MeanRMSE_pct - B.MeanRMSE_pct));   % mean RMSE can move by this much between seeds
v.exactNoise = max(100/cfg.nTrials, abs(base.MeanExact_pct - B.MeanExact_pct));
key = -[A.Passes]*1e6 + [A.MeanRMSE_pct];                    % most passes first, then the lowest mean RMSE
[~, ib] = min(key);
v.best = A(ib);
v.A = A;
v.base = base;
v.equiv = ([A.Passes] >= v.best.Passes - v.passNoise) & ([A.MeanRMSE_pct] <= v.best.MeanRMSE_pct + v.meanNoise);   % as good as the best within noise
lo = ib;
hi = ib;
while lo > 1 && v.equiv(lo - 1)
    lo = lo - 1;                                             % contiguous block of equivalent windows below the best
end
while hi < numel(A) && v.equiv(hi + 1)
    hi = hi + 1;                                             % and above the best
end
v.blockLo = A(lo).Window_s;
v.blockHi = A(hi).Window_s;
v.broad = (hi - lo + 1) >= 3;                                % broad = at least three neighbouring windows as good as the best
v.better = ([A.Passes] > base.Passes + v.passNoise) | (([A.Passes] >= base.Passes) & ([A.MeanRMSE_pct] < base.MeanRMSE_pct - v.meanNoise));
v.worse = ([A.Passes] < base.Passes - v.passNoise) | (([A.Passes] <= base.Passes) & ([A.MeanRMSE_pct] > base.MeanRMSE_pct + v.meanNoise));
v.ok = true;
end

function writeReport(cfg, outDir, ALL, C, PS)
L = {};
L{end+1} = 'STEP WINDOW FINE-TUNE REPORT (written automatically)';
L{end+1} = sprintf('Run folder: %s', outDir);
L{end+1} = sprintf('Software: %s. Simulated data only. R0 = %g m, Chen amplitude model %s, %d trials per setting, seeds A = %d and B = %d.', ...
    softwareName(), cfg.R0_m, cfg.AmplitudeMode, cfg.nTrials, cfg.SeedA, cfg.SeedB);
if cfg.Quick
    L{end+1} = 'QUICK TEST: 2 speeds and 3 trials. The numbers only show that the script runs. They are not results.';
end
L{end+1} = sprintf('Step windows: %s s. %.3f s is the power-of-two baseline (N = 8192). The others are exact windows (N even, 75 %% overlap, no zero padding).', ...
    sprintf('%g ', cfg.StepWindows), cfg.BaselineWindow);
L{end+1} = 'Walking speed window (stride length only): 0.372 s. Everything else as in Master_SimSweep_v5 (v5 extractor, Pth = max(10 - SNR, -40) dB, true Tc).';
L{end+1} = 'A change counts only when it is larger than the difference between seed A and seed B of the 0.186 s baseline at the same speed, SNR and parameter';
L{end+1} = sprintf('(RMSE noise floor %.2f pp, exact trials noise at least one trial = %.1f pp). Accuracy here = 100 - RMSE %%.', cfg.RMSEFloor_pp, 100/cfg.nTrials);
L{end+1} = '';

L{end+1} = sprintf('1. PASSES (RMSE below %g %%) OUT OF %d CONDITIONS PER WINDOW AND SNR (%d speeds x 3 parameters)', cfg.Target_pct, 3*numel(cfg.Speeds), numel(cfg.Speeds));
hdr = '  window s (actual)  ';
for snr = cfg.SNRs
    hdr = [hdr sprintf('| %2d dB: pass  meanRMSE%%  exact%%  better/worse ', snr)]; %#ok<AGROW>
end
L{end+1} = hdr;
for win = [cfg.StepWindows, NaN]
    tag = 'A';
    w = win;
    lab = '';
    if isnan(win)
        tag = 'B';
        w = cfg.BaselineWindow;
        lab = ' seed B';
    elseif abs(win - cfg.BaselineWindow) < 1e-9
        lab = ' base';
    end
    line = '';
    actW = NaN;                                              % actual window, filled from the first SNR that has a result
    for snr = cfg.SNRs
        p = PS(strcmp({PS.SeedTag}, tag) & [PS.SNR_dB] == snr & abs([PS.Window_s] - w) < 1e-9);
        if isempty(p)
            line = [line '|        -                                   ']; %#ok<AGROW>
            continue;
        end
        line = [line sprintf('|  %2d   %6.2f   %6.1f   %s ', p.Passes, p.MeanRMSE_pct, p.MeanExact_pct, ...
            ifElse(strcmp(tag, 'A'), sprintf('%2d/%-2d', p.RMSE_better, p.RMSE_worse), '  -  '))]; %#ok<AGROW>
        actW = p.Window_act_s;
    end
    L{end+1} = sprintf('  %.3f (%.4f)%-7s %s', w, actW, lab, line);
end
L{end+1} = sprintf('  better/worse = conditions (out of %d) whose RMSE', 3*numel(cfg.Speeds));
L{end} = [L{end} ' is better or worse than the baseline by more than the seed noise of that condition.'];
L{end+1} = '  The seed B row is the baseline again with other noise. The difference between "base" and "seed B" is the run-to-run noise.';
L{end+1} = '';

L{end+1} = '2. IS THE OPTIMUM BROAD OR NARROW, AND WHICH WINDOW IS BEST';
L{end+1} = '  Best = most passes, then the lowest mean RMSE over all speeds and parameters. "As good as the best" = passes within the seed noise of the passes (at least 1)';
L{end+1} = '  and mean RMSE within the seed noise of the mean RMSE. Broad = at least three neighbouring windows are as good as the best.';
nBetter = zeros(1, numel(cfg.StepWindows));
nWorse = zeros(1, numel(cfg.StepWindows));
bestList = [];
broadCount = 0;
nSNR = 0;
for snr = cfg.SNRs
    v = verdictPerSNR(cfg, PS, snr);
    if ~v.ok
        continue;
    end
    nSNR = nSNR + 1;
    eqW = [v.A(v.equiv).Window_s];
    L{end+1} = sprintf('  %2d dB: best %.3f s (%d passes, mean RMSE %.2f %%, exact %.1f %%). Baseline 0.186 s: %d passes, %.2f %%, %.1f %%. Seed noise: %d passes, %.2f pp.', ...
        snr, v.best.Window_s, v.best.Passes, v.best.MeanRMSE_pct, v.best.MeanExact_pct, v.base.Passes, v.base.MeanRMSE_pct, v.base.MeanExact_pct, v.passNoise, v.meanNoise);
    L{end+1} = sprintf('         as good as the best: %s s. Neighbouring block %.3f to %.3f s, so the optimum is %s.', sprintf('%.3f ', eqW), v.blockLo, v.blockHi, ifElse(v.broad, 'BROAD', 'NARROW'));
    if abs(v.best.Window_s - min(cfg.StepWindows)) < 1e-9 || abs(v.best.Window_s - max(cfg.StepWindows)) < 1e-9
        L{end+1} = '         The best window is at the end of the tested range, so the real optimum may lie outside it.';   % edge warning
    end
    bw = [v.A(v.better).Window_s];
    ww = [v.A(v.worse).Window_s];
    L{end+1} = sprintf('         better than the baseline beyond noise: %s. Worse: %s.', listOrNone(bw), listOrNone(ww));
    for i = 1:numel(cfg.StepWindows)
        j = find(abs([v.A.Window_s] - cfg.StepWindows(i)) < 1e-9, 1);
        if ~isempty(j)
            nBetter(i) = nBetter(i) + v.better(j);
            nWorse(i) = nWorse(i) + v.worse(j);
        end
    end
    bestList(end+1) = v.best.Window_s;                       %#ok<AGROW>
    broadCount = broadCount + v.broad;
end
if nSNR > 0
    cand = find(nWorse == 0 & nBetter > 0);                  % never worse than the baseline, better at one SNR or more
    if isempty(cand)
        rec = sprintf('Keep %.3f s. No window beats the baseline beyond the seed noise without losing at another SNR.', cfg.BaselineWindow);
    else
        [~, o] = sort(-nBetter(cand));
        c1 = cand(o(1));
        rec = sprintf('%.3f s is better than the baseline beyond noise at %d of %d SNRs and never worse. It is the best candidate.', cfg.StepWindows(c1), nBetter(c1), nSNR);
    end
    L{end+1} = sprintf('  Summary: the optimum is broad at %d of %d SNRs. Best window per SNR: %s s.', broadCount, nSNR, sprintf('%.3f ', bestList));
    L{end+1} = ['  Recommendation: ' rec];
end
L{end+1} = '';

L{end+1} = '3. DOES ANY WINDOW REMOVE THE 1.6 m/s MISSED STEP AT 70 TO 80 dB';
L{end+1} = '  step recognition at each window: miss % / extra % / exact trials %. "base" = 0.186 s seed A, "seed B" = the same with other noise.';
L{end+1} = '  70 dB is an extra check run only for these speeds. It is not counted in the passes.';
snrsMiss = unique([cfg.MissCheckSNRs, cfg.SNRs(cfg.SNRs >= 60)]);
for v0 = cfg.MissCheckSpeeds
    for snr = snrsMiss
        b = pickRow(ALL, 'StepTime', v0, snr, cfg.BaselineWindow, 'A');
        bB = pickRow(ALL, 'StepTime', v0, snr, cfg.BaselineWindow, 'B');
        if isempty(b)
            continue;
        end
        line = '';
        clean = [];
        for win = cfg.StepWindows
            r = pickRow(ALL, 'StepTime', v0, snr, win, 'A');
            if isempty(r)
                continue;
            end
            line = [line sprintf('%.3f: %.1f/%.1f/%.0f  ', win, r.MissRate_pct, r.ExtraRate_pct, r.PerfectTrials_pct)]; %#ok<AGROW>
            if r.MissRate_pct == 0 && r.ExtraRate_pct == 0
                clean(end+1) = win;                          %#ok<AGROW>
            end
        end
        L{end+1} = sprintf('  %.1f m/s, %g dB: %s', v0, snr, line);
        bText = sprintf('baseline miss %.1f %% (seed B %.1f %%)', b.MissRate_pct, bB.MissRate_pct);
        if b.MissRate_pct == 0 && bB.MissRate_pct == 0
            L{end+1} = sprintf('     The baseline missed no step with these seeds (%s), so the miss of the final sweep did not appear here.', bText);
        else
            noise = abs(b.MissRate_pct - bB.MissRate_pct);
            fixed = [];
            for win = clean
                r = pickRow(ALL, 'StepTime', v0, snr, win, 'A');
                if b.MissRate_pct - r.MissRate_pct > noise
                    fixed(end+1) = win;                      %#ok<AGROW>
                end
            end
            L{end+1} = sprintf('     %s. Windows with no miss and no extra: %s. Of these, the drop in misses is larger than the seed noise (%.1f pp) for: %s.', ...
                bText, listOrNone(clean), noise, listOrNone(fixed));
        end
    end
end
L{end+1} = '';

L{end+1} = '4. CONDITIONS THAT CHANGED BEYOND THE SEED NOISE (RMSE), per window';
L{end+1} = '  parameter speed SNR: window RMSE % (baseline RMSE %, noise pp)';
for win = cfg.StepWindows
    if abs(win - cfg.BaselineWindow) < 1e-9
        continue;
    end
    k = find(abs([C.Window_s] - win) < 1e-9 & ~strcmp({C.RMSE_change}, 'within noise') & [C.ExtraCheck] == 0);
    if isempty(k)
        L{end+1} = sprintf('  %.3f s: none', win);
        continue;
    end
    txtB = '';
    txtW = '';
    for q = k
        s = sprintf('%s %.1f %gdB: %.2f (%.2f, %.2f); ', shortName(C(q).Parameter), C(q).Speed_mps, C(q).SNR_dB, C(q).RMSE_pct, C(q).RMSE_base, C(q).NoiseRMSE_pp);
        if strcmp(C(q).RMSE_change, 'better')
            txtB = [txtB s]; %#ok<AGROW>
        else
            txtW = [txtW s]; %#ok<AGROW>
        end
    end
    L{end+1} = sprintf('  %.3f s better: %s', win, ifElse(isempty(txtB), 'none', txtB));
    L{end+1} = sprintf('  %.3f s worse : %s', win, ifElse(isempty(txtW), 'none', txtW));
end
L{end+1} = '';
L{end+1} = 'NOTES';
L{end+1} = '  The exact windows use an even number of samples. With an odd N the frequency axis of the extractor would be half a bin out after fftshift.';
L{end+1} = '  The walk interior used for scoring (0.1 s + half a window from each end) is the same rule as before, so longer windows score slightly fewer steps.';
L{end+1} = '  Every window sees the same noise draw (paired test). Seed B is only run for the baseline, so the noise is one difference per condition, not a full spread.';
L{end+1} = '  Simulated data only. The parameters are to be frozen before the IPM-365 recordings are processed.';
fid = fopen(fullfile(outDir, 'FineTune_Report.txt'), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
end

function s = listOrNone(x)
if isempty(x)
    s = 'none';
else
    s = strjoin(arrayfun(@(a) sprintf('%.3f s', a), x, 'UniformOutput', false), ', ');
end
end

function s = shortName(p)
switch p
    case 'StepTime', s = 'ST';
    case 'StrideDuration', s = 'SD';
    otherwise, s = 'SL';
end
end

function s = ifElse(tf, a, b)
if tf, s = a; else, s = b; end
end

function makePlots(cfg, outDir, ALL, PS)
speeds = cfg.Speeds;
n = numel(speeds);
cm = parula(n + 3);                                          % ordered colours, slow = dark blue to fast = green
cm = cm(1:n, :);                                             % the pale yellow end is left out, it is hard to see on white
mk = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x'};
panels = {'StepTime', 'PerfectTrials_pct', 'Step recognition: exact trials (%)'; ...
          'StepTime', 'RMSE_pct', 'Step time RMSE (%)'; ...
          'StrideDuration', 'RMSE_pct', 'Stride duration RMSE (%)'; ...
          'StrideLength', 'RMSE_pct', 'Stride length RMSE (%)'};
A = ALL(strcmp({ALL.SeedTag}, 'A'));
for snr = cfg.SNRs
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 850]);
    for q = 1:4
        ax = subplot(2, 2, q, 'Parent', fig);
        hold(ax, 'on');
        yMax = 0;
        h = [];
        lab = {};
        for i = 1:n
            R = A(strcmp({A.Parameter}, panels{q, 1}) & abs([A.Speed_mps] - speeds(i)) < 1e-9 & [A.SNR_dB] == snr);
            if isempty(R)
                continue;
            end
            [~, o] = sort([R.Window_act_s]);
            R = R(o);
            y = [R.(panels{q, 2})];
            yMax = max([yMax, y(isfinite(y))]);
            h(end+1) = plot(ax, [R.Window_act_s], y, ['-' mk{1 + mod(i - 1, numel(mk))}], 'Color', cm(i, :), ...
                'MarkerFaceColor', cm(i, :), 'LineWidth', 1.4, 'MarkerSize', 5);   %#ok<AGROW>
            lab{end+1} = sprintf('%.1f', speeds(i));         %#ok<AGROW>
        end
        if strcmp(panels{q, 2}, 'PerfectTrials_pct')
            ylim(ax, [0 105]);
        else
            yTop = min(30, max(8, 1.1*yMax));                % RMSE axis stops at 30 %
            plot(ax, [0.1 0.32], [cfg.Target_pct cfg.Target_pct], 'k--', 'LineWidth', 1.2);   % the 5 % target
            ylim(ax, [0 yTop]);
            if yMax > yTop
                text(ax, 0.99, 0.97, sprintf('Values above %g %% are clipped (see the CSV).', yTop), 'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8);
            end
        end
        yl = get(ax, 'YLim');
        plot(ax, 8192/44100*[1 1], yl, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);   % the 0.186 s baseline
        xlim(ax, [0.11 0.31]);
        set(ax, 'XTick', cfg.StepWindows, 'FontSize', 10, 'Box', 'on');
        grid(ax, 'on');
        xlabel(ax, 'Step detection window (s)', 'FontSize', 11);
        ylabel(ax, panels{q, 3}, 'FontSize', 11);
        title(ax, panels{q, 3}, 'FontSize', 12);
        if q == 2 && ~isempty(h)
            lg = legend(ax, h, lab, 'Location', 'eastoutside', 'FontSize', 9);
            try, title(lg, 'Speed (m/s)'); catch, end        % key heading
        end
    end
    annotation(fig, 'textbox', [0.02 0.95 0.96 0.04], 'String', sprintf(['Step window fine-tune at %d dB, %d trials per point, seed A. ' ...
        'Dotted line = 0.186 s power-of-two baseline, dashed line = 5 %% target.'], snr, cfg.nTrials), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11);
    saveLightFigure(fig, fullfile(outDir, sprintf('FineTune_SNR%d', snr)));
    close(fig);
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 560]);
ax = axes('Parent', fig);
hold(ax, 'on');
cs = [0.10 0.30 0.65; 0.20 0.55 0.75; 0.85 0.45 0.10; 0.55 0.15 0.30];   % one colour per SNR
h = [];
lab = {};
for j = 1:numel(cfg.SNRs)
    p = PS(strcmp({PS.SeedTag}, 'A') & [PS.SNR_dB] == cfg.SNRs(j));
    if isempty(p)
        continue;
    end
    [~, o] = sort([p.Window_act_s]);
    p = p(o);
    h(end+1) = plot(ax, [p.Window_act_s], [p.Passes], '-o', 'Color', cs(1 + mod(j - 1, 4), :), 'MarkerFaceColor', cs(1 + mod(j - 1, 4), :), 'LineWidth', 1.6);   %#ok<AGROW>
    lab{end+1} = sprintf('%d dB', cfg.SNRs(j));              %#ok<AGROW>
    pb = PS(strcmp({PS.SeedTag}, 'B') & [PS.SNR_dB] == cfg.SNRs(j));
    if ~isempty(pb)
        plot(ax, pb.Window_act_s, pb.Passes, 'x', 'Color', cs(1 + mod(j - 1, 4), :), 'MarkerSize', 10, 'LineWidth', 1.6);   % baseline on seed B
    end
end
xlim(ax, [0.11 0.31]);
ylim(ax, [0 3*numel(cfg.Speeds) + 1]);
set(ax, 'XTick', cfg.StepWindows, 'FontSize', 11, 'Box', 'on');
grid(ax, 'on');
xlabel(ax, 'Step detection window (s)', 'FontSize', 12);
ylabel(ax, sprintf('Conditions below %g %% RMSE (of %d)', cfg.Target_pct, 3*numel(cfg.Speeds)), 'FontSize', 12);
title(ax, {'Passes against step window (step time, stride duration and stride length)', 'x = 0.186 s baseline on seed B (run-to-run noise)'}, 'FontSize', 12);
if ~isempty(h)
    lg = legend(ax, h, lab, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);
    try, title(lg, 'SNR'); catch, end
end
saveLightFigure(fig, fullfile(outDir, 'FineTune_Passes'));
close(fig);
end

function logMsg(outDir, fmt, varargin)
msg = sprintf(fmt, varargin{:});
stamp = datestr(now, 'HH:MM:SS');
fprintf('[%s] %s\n', stamp, msg);
fid = fopen(fullfile(outDir, 'Run_Log.txt'), 'a');
if fid > 0
    fprintf(fid, '[%s] %s\n', stamp, msg);
    fclose(fid);
end
end

function s = softwareName()
if exist('OCTAVE_VERSION', 'builtin')
    s = ['Octave ' OCTAVE_VERSION ' (a stand-in, NOT MATLAB)'];
else
    s = ['MATLAB ' version];
end
end

function setSeed(s)
if exist('OCTAVE_VERSION', 'builtin')
    randn('state', s);                                       %#ok<RAND>
    rand('state', s);                                        %#ok<RAND>
else
    rng(s, 'twister');
end
end

function z = trialNoise(D, seedBase, k)
setSeed(seedBase + 7919*k + round(1000*D.v0));
n = numel(D.clean);
z = randn(1, n) + 1i*randn(1, n);
end

function txt = patchRegex(txt, pattern, replacement, expected, label, tag, logFile)
m = regexp(txt, pattern, 'match');
if numel(m) ~= expected
    error('Patch "%s" (%s): expected %d match, found %d. The extractor script has changed, nothing was patched.', label, tag, expected, numel(m));
end
rep = strrep(replacement, '\', '\\');                         % the replacement is plain text, not a regular expression
rep = strrep(rep, '$', '\$');
txt = regexprep(txt, pattern, rep);
logPatch(logFile, tag, label, m{1}, replacement);
end

function logPatch(logFile, tag, label, oldText, newText)
fid = fopen(logFile, 'a');
if fid > 0
    fprintf(fid, '[%s]\n  patch : %s\n  OLD   : %s\n  NEW   : %s\n\n', tag, label, strtrim(oldText), strrep(strtrim(newText), char(10), [char(10) '          ']));
    fclose(fid);
end
end

function o = runSpeed(m_P, m_x, m_D, m_vs)
global VX
VX = defaultVX();
VX.WSexact = m_vs.exact;
VX.WSbandA = m_vs.bandA;
VX.WSpasses = m_vs.passes;
VX.WSHop = m_vs.hop;
o = struct('ok', false, 'meanSpeed', NaN, 'vel', [], 'tv', [], 'cpiAct', NaN, 'nIter', NaN, ...
    'bandLow', NaN, 'bandHigh', NaN, 'centrePrev', NaN, 'framePeriod', NaN, 'err', '');
try
    RdrRxSignal = m_x;                                       %#ok<NASGU>
    fs = m_D.fs;                                             %#ok<NASGU>
    fc = m_D.fc;                                             %#ok<NASGU>
    c = m_D.c;                                               %#ok<NASGU>
    v0 = m_D.v0;                                             %#ok<NASGU>
    CPI_s = m_vs.cpi;                                        %#ok<NASGU>
    evalc('run(m_P.ws);');
    o.meanSpeed = meanSpeed_mps;
    o.vel = velocityTrack_mps;
    o.tv = final.time_axis;
    o.cpiAct = CPI_actual_s_WS;
    o.nIter = numel(resultsIter);
    o.bandLow = final.freq_low;
    o.bandHigh = final.freq_high;
    if numel(resultsIter) >= 2
        o.centrePrev = resultsIter(end-1).centreFreq_Hz;
    end
    o.framePeriod = final.time_axis(2) - final.time_axis(1);
    o.ok = isfinite(o.meanSpeed);
catch ME
    o.err = ME.message;
end
VX = defaultVX();
end

function [ok, res] = chk(ok, res, cond, text, errText)
if cond
    res{end+1} = ['PASS  ' text];
else
    res{end+1} = ['FAIL  ' text ' ' errText];
    ok = false;
end
end

function D = getSpeedData(m_P, m_cfg, v0)
SimulatorVerbose = false;                                    %#ok<NASGU>
SimulatorMakePlots = false;                                  %#ok<NASGU>
SimulatorR0_m = m_cfg.R0_m;                                  %#ok<NASGU>
AmplitudeMode = m_cfg.AmplitudeMode;                         %#ok<NASGU>
SNR_dB = 40;                                                 %#ok<NASGU>
CPI_s = 0.1;                                                 %#ok<NASGU>
evalc('run(m_P.sim);');
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
D.truth = trueStepTimes_s_corrected(:).';
D.Ngt = NumSamplesPerFrame_gt;
D.Nsim = NumSamplesPerFrame;
D.noiseSim = noise;
end

function x = addNoise(clean, Nw, snr_dB, z)
noisePower = 10^(-(snr_dB - 10*log10(Nw))/10);               % same equation as the simulator
x = clean + sqrt(noisePower/2) * z;                          % z is one unit complex noise draw
end

function sc = scoreStepTimes(det, D, cpiAct, key)
if nargin < 4
    key = det;                                               % key = the time stamps that decide which steps are inside the walk
end
tol = min(0.2, D.Tc/4);
sc = scoreEvents(det, D.truth, D, cpiAct, tol, key);         % original time stamps (window start)
b = scoreEvents(det + cpiAct/2, D.truth + (D.Ngt/D.fs)/2, D, cpiAct, tol, key + cpiAct/2);   % both moved to window centre
sc.missC = b.miss;
sc.extraC = b.extra;
end

function s = scoreEvents(det, truth, D, cpiAct, tol, key)
if nargin < 6
    key = det;
end
edge = 0.1 + cpiAct/2;                                       % ignore the walk start and end
lo = D.buffer + edge;
hi = D.buffer + D.simT - edge;
detN = det(key > lo & key < hi);                            % inside the walk is decided by the frame time stamp
trN = truth(truth > lo & truth < hi);
map = matchNew(detN, trN, tol);
nTruth = numel(trN);
nMatched = nnz(map > 0);
s.miss = pct(nTruth - nMatched, nTruth);
s.extra = pct(nnz(map == 0), nTruth);
s.perfect = double(nTruth == nMatched && nnz(map == 0) == 0);   % 1 if every true step was found and nothing extra was detected
s.e = [];
s.g = [];
s.offset = NaN;
if numel(detN) >= 2
    valid = map(1:end-1) > 0 & diff(map) == 1;               % both ends matched to neighbouring true steps
    ii = find(valid);
    est = diff(detN);
    s.g = trN(map(ii)+1) - trN(map(ii));
    s.e = est(ii) - s.g;
end
m = map > 0;
if any(m)
    s.offset = median(detN(m) - trN(map(m)));
end
end

function [e, g] = strideLengthErrors(W, S, D)
e = NaN;
g = NaN;
if ~W.ok || ~S.ok || numel(S.det) < 3
    return;
end
try
    edge = 0.1 + S.cpiAct/2;                                 % same walk interior as the step time and stride duration scoring
    inside = S.key > D.buffer + edge & S.key < D.buffer + D.simT - edge;
    detIn = S.det(inside);                                   % steps inside the walk only, decided by the frame time stamp
    if numel(detIn) < 3
        return;
    end
    vel = fillmissing(W.vel, 'linear', 'EndValues', 'nearest');
    dsp = cumtrapz(W.tv, vel);                               % torso position from velocity
    atDet = interp1(W.tv, dsp, detIn, 'linear');
    est = abs(atDet(3:end) - atDet(1:end-2));
    trueTimes = D.buffer : D.Tc/2 : D.buffer + D.simT;
    trueRange = interp1(D.t, D.Torso_R, trueTimes, 'linear');
    trueLen = abs(trueRange(3:end) - trueRange(1:end-2));
    [em, tm] = matchStridesExisting(est, detIn(1:end-2), trueLen, trueTimes(1:end-2), D.Tc/4);
    e = em - tm;
    g = tm;
catch
    e = NaN;
    g = NaN;
end
end

function [estimatedMatched, trueMatched] = matchStridesExisting(estimated, estimatedStarts, truth, trueStarts, tolerance)
estimatedMatched = [];
trueMatched = [];
for k = 1:numel(estimated)
    if isempty(trueStarts)
        continue;
    end
    [distance, index] = min(abs(trueStarts - estimatedStarts(k)));
    if distance <= tolerance
        estimatedMatched(end+1) = estimated(k);              %#ok<AGROW>
        trueMatched(end+1) = truth(index);                   %#ok<AGROW>
    end
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
            counts(end+1) = count(i, j) + 1;                 %#ok<AGROW>
            costs(end+1) = cost(i, j) + distance;            %#ok<AGROW>
            actions(end+1) = 3;                              %#ok<AGROW>
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

function v = pct(a, b)
if b > 0
    v = 100*a/b;
else
    v = NaN;
end
end

function v = rmsOf(x)
x = x(isfinite(x));
if isempty(x)
    v = NaN;
else
    v = sqrt(mean(x.^2));
end
end

function v = meanFinite(x)
x = x(isfinite(x));
if isempty(x)
    v = NaN;
else
    v = mean(x);
end
end

function v = stdFinite(x)
x = x(isfinite(x));
if numel(x) < 2
    v = NaN;
else
    v = std(x);
end
end

function writeRowsCSV(filename, rows)
if isempty(rows)
    return;
end
names = fieldnames(rows);
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', strjoin(names', ','));
for r = 1:numel(rows)
    parts = cell(1, numel(names));
    for k = 1:numel(names)
        val = rows(r).(names{k});
        if ischar(val)
            parts{k} = ['"' val '"'];
        elseif isempty(val)
            parts{k} = 'NaN';
        else
            parts{k} = sprintf('%.10g', val);   % 10 digits, so 8-digit seeds are written in full
        end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
fclose(fid);
end

function rows = addRow(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row;
end
end

function styleLightFigure(fig)
try, theme(fig, 'light'); catch, end                                 % R2025a and newer only, ignored on older versions
try
    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        set(ax(k), 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.35);
        try, ax(k).Toolbar.Visible = 'off'; catch, end               % hide the toolbar icons so they do not appear in the export
    end
    lgs = findall(fig, 'Type', 'legend');
    for k = 1:numel(lgs)
        set(lgs(k), 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);
    end
    txt = findall(fig, 'Type', 'text');
    for k = 1:numel(txt)
        set(txt(k), 'Color', [0.1 0.1 0.1]);
    end
catch
end
end

function saveLightFigure(fig, basePath)
styleLightFigure(fig);
try
    exportgraphics(fig, [basePath '.png'], 'Resolution', 200, 'BackgroundColor', 'white');
    exportgraphics(fig, [basePath '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
catch
    try
        print(fig, [basePath '.png'], '-dpng', '-r200');
    catch
        saveas(fig, [basePath '.png']);
    end
end
end

