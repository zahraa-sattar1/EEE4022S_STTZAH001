%% Step-detector variants
% Zahraa Sattar | EEE4022S
%
% Run the overlap and peak-drop comparison used in Figure 4.5 and Appendix B.2.
% The default runs E1, using historical extractor inputs and shared noise draws.

function Master_SimVariants_v1(userCfg)

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
if cfg.Quick
    if ~isfield(userCfg, 'Speeds'), cfg.Speeds = [0.6 1.4]; end
    if ~isfield(userCfg, 'nTrials'), cfg.nTrials = 3; end
    if ~isfield(userCfg, 'StepSNRs'), cfg.StepSNRs = 60; end
    if ~isfield(userCfg, 'SpeedSNRs'), cfg.SpeedSNRs = 30; end
end
cfg.AllSNRs = unique([cfg.StepSNRs(:).', cfg.SpeedSNRs(:).']);   % SNRs the baseline needs

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
    outDir = fullfile(codeFolder, ['SimVariants_' datestr(now, 'yyyymmdd_HHMMSS')]);
    mkdir(outDir);                                           % always a new folder, nothing is overwritten
else
    outDir = cfg.ResumeFolder;
    if ~exist(outDir, 'dir')
        outDir = fullfile(codeFolder, cfg.ResumeFolder);
    end
end
cfg.ResumeFolder = outDir;
cfg.OutDir = outDir;
cfg.CodeFolder = codeFolder;
if ~isfile(fullfile(outDir, 'Config.mat'))
    save(fullfile(outDir, 'Config.mat'), 'cfg');
end

global VX WCACHE PCACHE
VX = defaultVX();                                            % option block read by the patched copies
WCACHE = containers.Map();                                   % walking-speed results shared between variants
PCACHE = containers.Map();                                   % first-pass Tc estimates shared between variants

logMsg(outDir, 'Output folder: %s', outDir);
logMsg(outDir, 'Software: %s', softwareName());
logMsg(outDir, 'Seeds: A = %d, B = %d, validation = %d (trial seed = base + 7919*trial + 1000*speed)', cfg.SeedA, cfg.SeedB, cfg.SeedVal);
logMsg(outDir, 'Speeds %s, trials %d, step SNRs %s, speed SNRs %s, quick mode %d', mat2str(cfg.Speeds), cfg.nTrials, mat2str(cfg.StepSNRs), mat2str(cfg.SpeedSNRs), cfg.Quick);
logMsg(outDir, 'Extractor frequency limit that every copy must use: 350*24/5.8 = %.4f Hz', 350*24/5.8);

P = makeWorkingCopies(cfg, codeFolder, outDir);
logMsg(outDir, 'Working copies made. Patch log: Patch_Log.txt');
selfTest(P, cfg, outDir);

wantExp = cfg.Run;
tStart = tic;
tExp = struct();                                             % seconds spent per experiment
for sIdx = 1:numel(cfg.Speeds)
    v0 = cfg.Speeds(sIdx);
    logMsg(outDir, '===== v0 = %.2f m/s (%d of %d), elapsed %.1f min =====', v0, sIdx, numel(cfg.Speeds), toc(tStart)/60);
    todo = {};
    for q = {'E0', 'E1', 'E2', 'E3', 'E4a', 'E5', 'E6'}
        e = q{1};
        if any(strcmp(wantExp, e(1:2))) && ~isfile(expFile(outDir, e, v0))
            todo{end+1} = e;                                 %#ok<AGROW>
        end
    end
    if isempty(todo)
        logMsg(outDir, 'Phase 1 already finished for this speed, skipping.');
        continue;
    end
    try
        D = getSpeedData(P, cfg, v0);                        % one simulator run gives the clean signal for every trial
    catch ME
        logMsg(outDir, 'Simulator failed at v0 = %.2f: %s', v0, ME.message);
        continue;
    end
    WCACHE = containers.Map();                               % fresh caches for this speed
    PCACHE = containers.Map();
    for q = 1:numel(todo)
        e = todo{q};
        t0 = tic;
        try
            switch e
                case 'E0', rows = runE0(P, cfg, D);
                case 'E1', rows = runE1(P, cfg, D);
                case 'E2', rows = runE2(P, cfg, D);
                case 'E3', rows = runE3(P, cfg, D);
                case 'E4a', rows = runE4a(P, cfg, D);
                case 'E5', rows = runE5(P, cfg, D);
                case 'E6', rows = runE6(P, cfg, D);
            end
            save(expFile(outDir, e, v0), 'rows');            % saved after every experiment and speed
            dt = toc(t0);
            logMsg(outDir, '   %s finished in %.1f min (%d rows)', e, dt/60, numel(rows));
            tExp = addTime(tExp, e, dt);
        catch ME
            logMsg(outDir, '   %s stopped early at v0 = %.2f: %s', e, v0, ME.message);
            for s2 = 1:min(4, numel(ME.stack))
                logMsg(outDir, '      at %s line %d', ME.stack(s2).name, ME.stack(s2).line);
            end
        end
    end
end

if any(strcmp(wantExp, 'E4'))
    try
        sel4 = selectE4(cfg, outDir);                        % fixed and per-speed CPI, chosen from the sweep
        save(fullfile(outDir, 'E4_selection.mat'), 'sel4');
        logMsg(outDir, 'E4 selection: %s', sel4.text);
    catch ME
        sel4 = [];
        logMsg(outDir, 'E4 selection failed: %s', ME.message);
    end
else
    sel4 = [];
end
for sIdx = 1:numel(cfg.Speeds)
    v0 = cfg.Speeds(sIdx);
    todo = {};
    if any(strcmp(wantExp, 'E4')) && ~isempty(sel4) && ~isfile(expFile(outDir, 'E4v', v0))
        todo{end+1} = 'E4v';
    end
    if any(strcmp(wantExp, 'E7')) && ~isfile(expFile(outDir, 'E7', v0))
        todo{end+1} = 'E7';
    end
    if isempty(todo)
        continue;
    end
    logMsg(outDir, '===== (validation phase) v0 = %.2f m/s =====', v0);
    try
        D = getSpeedData(P, cfg, v0);
    catch ME
        logMsg(outDir, 'Simulator failed at v0 = %.2f: %s', v0, ME.message);
        continue;
    end
    WCACHE = containers.Map();
    PCACHE = containers.Map();
    for q = 1:numel(todo)
        e = todo{q};
        t0 = tic;
        try
            if strcmp(e, 'E4v')
                rows = runE4v(P, cfg, D, sel4);
            else
                rows = runE7(P, cfg, D);
            end
            save(expFile(outDir, e, v0), 'rows');
            dt = toc(t0);
            logMsg(outDir, '   %s finished in %.1f min (%d rows)', e, dt/60, numel(rows));
            tExp = addTime(tExp, e, dt);
        catch ME
            logMsg(outDir, '   %s stopped early at v0 = %.2f: %s', e, v0, ME.message);
            for s2 = 1:min(4, numel(ME.stack))
                logMsg(outDir, '      at %s line %d', ME.stack(s2).name, ME.stack(s2).line);
            end
        end
    end
end

logMsg(outDir, 'All experiments done in %.1f min. Writing results.', toc(tStart)/60);
try
    writeTiming(outDir, tExp);
    assembleResults(cfg, outDir);
catch ME
    logMsg(outDir, 'Writing the results failed: %s', ME.message);
    for s2 = 1:min(5, numel(ME.stack))
        logMsg(outDir, '      at %s line %d', ME.stack(s2).name, ME.stack(s2).line);
    end
end
logMsg(outDir, 'Finished. Open Variant_Report.txt in %s', outDir);
end

function cfg = defaultConfig()
cfg.Speeds = [0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];      % m/s
cfg.nTrials = 30;                                            % noise draws per condition (use 15 if the run is too slow, and say so)
cfg.Quick = false;                                           % true = 2 speeds, 3 trials, one SNR
cfg.SeedA = 20261001;                                        % seed base of the main runs
cfg.SeedB = 20262002;                                        % seed base of the second baseline run
cfg.SeedVal = 20263003;                                      % seed base of the fresh validation draws (E4 validation, E7)
cfg.Run = {'E1'};  % experiments to run
cfg.R0_m = 10;                                               % simulator range
cfg.AmplitudeMode = 'A';                                     % Chen mean sqrt(RCS) weights
cfg.Target_pct = 5;                                          % goal for RMSE %
cfg.MaxMissExtra_pct = 10;                                   % miss and extra limit for step time
cfg.StepCPIreq = 0.1;                                        % requested step window, actual 0.186 s
cfg.StepSNRs = [30 60];                                      % dB, step experiments (30 = noise limited, 60 = timing limited)
cfg.WSCPIreq = 0.2;
cfg.SpeedSNRs = [30 40];                                     % dB, walking speed experiments
cfg.E1Overlaps = [0.50 0.75 0.875];                          % overlap fractions
cfg.E1CScales = [1 0.5];                                     % multipliers on Wang's c
cfg.E1AlsoV5 = true;                                         % repeat E1 step detection on the v5 extractor
cfg.E2Hops = [0.25 0.125 0.0625];                            % hop as a fraction of the window
cfg.E2FixedLags = [0.10 0.15];                               % s, lags that need no Tc at all
cfg.E2LagK = 0.25;                                           % lag = K * Tc for the oracle and estimated lags
cfg.E2TcEstimator = 'twopass';                               % estimator used for the estimated lag
cfg.TcNominal = 1.2;                                         % s, nominal Tc of the first pass
cfg.E4CPI = [0.139 0.186 0.25 0.3 0.372 0.5 0.6 0.743 1.0 1.5];   % s, window lengths NOT forced to a power of two
cfg.E4K = [25 35 50];                                        % adaptive rule CPI = K / f, f = 160*v Hz
cfg.E4StepCPIreq = [0.025 0.075 0.1 0.2 0.4];                % requested step windows for the per-speed step choice
cfg.E4MinFrames = 8;
cfg.E5Bands = [0.15 0.25 0.4];                               % a: band edges (1-a) and (1+a) times the mean tracked frequency
cfg.E5Passes = [2 3 5];                                      % numbers of passes
cfg.E5CPIreq = 0.2;                                          % requested window used in E5, actual 0.372 s
cfg.E6Floors = [-30 -35 -40 -45];                            % dB
cfg.E6Margins = [6 10 14];                                   % dB above the estimated noise floor
cfg.E7TieTols = [0.25 0.5 1.0];                              % tie tolerances
cfg.E7RunDir = 'SimSweep_20261008_150337';                   % existing run whose Stage 1 and Stage 2 CSVs are re-used
cfg.StepRefine = 'centroid';                                 % used only by the v5 extractor
cfg.StepGapFill = 1;                                         % used only by the v5 extractor
cfg.MakePlots = true;
cfg.SimName = 'RadarSimulator_v11_torso_legs_arms_head_buffered_noplots';
cfg.WSName = 'Finding_WalkingSpeed_v4_noplots';
cfg.STName = 'Finding_StepRecognition_v4b_noplots';          % main extractor = the v3 settings (v4b + cPeakScale 0.5)
cfg.STName2 = 'Finding_StepRecognition_v5_noplots';          % improved extractor, used only for E1 and the v5 baseline
cfg.CodeFolder = fullfile(fileparts(mfilename('fullpath')), 'Report Sources');
cfg.ResumeFolder = '';
cfg.CopyHook = [];                                           % only used for testing outside MATLAB
end

function VX = defaultVX()
VX = struct('Hop', 0.25, 'CScale', 0.5, 'CHop', false, 'LagSec', [], 'PthMode', 'coupled', 'PthMargin', 10, ...
    'WSexact', false, 'WSbandA', 0, 'WSpasses', 3, 'WSHop', 0.25);
end

function f = expFile(outDir, e, v0)
f = fullfile(outDir, sprintf('%s_v0_%.2f.mat', e, v0));
end

function t = addTime(t, e, dt)
if isfield(t, e)
    t.(e) = t.(e) + dt;
else
    t.(e) = dt;
end
end

function writeTiming(outDir, t)
fid = fopen(fullfile(outDir, 'Run_Timing.csv'), 'w');
fprintf(fid, 'Experiment,Minutes\n');
n = fieldnames(t);
for k = 1:numel(n)
    fprintf(fid, '%s,%.2f\n', n{k}, t.(n{k})/60);
end
fclose(fid);
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

function P = makeWorkingCopies(cfg, codeFolder, outDir)
codeDir = fullfile(outDir, 'WorkingCopies');
if ~exist(codeDir, 'dir')
    mkdir(codeDir);
end
logFile = fullfile(outDir, 'Patch_Log.txt');
P.sim = copyOne(fullfile(codeFolder, [cfg.SimName '.m']), fullfile(codeDir, 'WC_Simulator.m'), 'sim', cfg, logFile);
P.ws = copyOne(fullfile(codeFolder, [cfg.WSName '.m']), fullfile(codeDir, 'WC_WalkingSpeed.m'), 'ws', cfg, logFile);
P.st = copyOne(fullfile(codeFolder, [cfg.STName '.m']), fullfile(codeDir, 'WC_StepRecognition.m'), 'st', cfg, logFile);
P.st2 = copyOne(fullfile(codeFolder, [cfg.STName2 '.m']), fullfile(codeDir, 'WC_StepRecognition_v5.m'), 'st', cfg, logFile);
end

function dest = copyOne(src, dest, kind, cfg, logFile)
if isfile(dest)
    return;                                                  % never overwrite, reuse on resume
end
assert(isfile(src), 'Missing file: %s', src);
txt = fileread(src);
tag = [kind ':' src];
txt = regexprep(txt, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*(?:\r?\n|$)', '');   % no clear all / close all
if ~strcmp(kind, 'sim')
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
        'freq_high_Yardibi = 350*24/5.8;   % this copy: 350 Hz limit for every variant', 1, 'extractor frequency limit 350 Hz', tag, logFile);
    txt = patchRegex(txt, 'noverlap = N\*3/4;[^\r\n]*', ...
        sprintf(['global VX\n        noverlap = N - max(1, round(N*VX.Hop));   %% hop set by the variant, 0.25 = the original 75 percent overlap']), ...
        1, 'overlap line (step spectrogram)', tag, logFile);
    txt = patchRegex(txt, 'cPeak\s*=\s*DopplerRes_Hz/freq_high_Yardibi;[^\r\n]*', ...
        sprintf(['global VX\ncPeak = DopplerRes_Hz/freq_high_Yardibi;\nif VX.CHop, cPeak = cPeak*VX.Hop/0.25; end   %% only when the variant asks for c to shrink with the hop']), ...
        1, 'cPeak line', tag, logFile);
    txt = patchRegex(txt, 'stepIdx8\s*=\s*find\([^;]*\)\s*\+\s*1;[^\r\n]*', ...
        sprintf(['global VX\nlagFrames_ = 1;   %% one frame = the original rule\n' ...
        'if ~isempty(VX.LagSec), lagFrames_ = max(1, round(VX.LagSec/frame_period_s)); end   %% lag in seconds turned into frames\n' ...
        'nF_ = numel(normalizedLegFrequency);\n' ...
        'ctr_ = (1+lagFrames_):(nF_-lagFrames_);   %% frames that have a neighbour at +/- the lag\n' ...
        'stepIdx8 = ctr_( (normalizedLegFrequency(ctr_) >= normalizedLegFrequency(ctr_-lagFrames_)) & ...\n' ...
        '                 (normalizedLegFrequency(ctr_) >= normalizedLegFrequency(ctr_+lagFrames_) + cPeak) );']), ...
        1, 'peak rule (stepIdx8)', tag, logFile);
    txt = patchRegex(txt, '\[legFreq_Hz, isDetected\] = SelectLegFrequency\(Pmedian, freq_axis_step, Pth\);[^\r\n]*', ...
        sprintf(['global VX\nif strcmp(VX.PthMode, ''noisefloor'')\n' ...
        '    nBack_ = max(1, round(0.20*size(P_step2, 2)));   %% lowest 20 percent of each bin, as in the background step\n' ...
        '    srt_ = sort(P_step2, 2, ''ascend'');\n' ...
        '    pb_ = mean(srt_(:, 1:nBack_), 2);\n' ...
        '    sig2_ = median(pb_)/0.1074;   %% 0.1074 = mean of the lowest 20 percent of an exponential cell\n' ...
        '    Pth = 10*log10(sig2_/(max(Pmedian(:)) + eps)) + VX.PthMargin;   %% noise floor relative to the peak, plus the margin\n' ...
        'end\nPthUsed_ = Pth;\n[legFreq_Hz, isDetected] = SelectLegFrequency(Pmedian, freq_axis_step, Pth);']), ...
        1, 'Pth noise-floor option', tag, logFile);
end
if strcmp(kind, 'ws')
    txt = patchRegex(txt, 'N_actual = 2\^\(nextpow2\(round\(CPI_s\*fs\)\)\);[^\r\n]*', ...
        sprintf(['global VX\nN_actual = 2^(nextpow2(round(CPI_s*fs)));\nif VX.WSexact, N_actual = round(CPI_s*fs); end   %% window not forced to a power of two']), ...
        1, 'actual CPI line', tag, logFile);
    txt = patchRegex(txt, 'numIterations\s*=\s*3;[^\r\n]*', ...
        sprintf(['global VX\nnumIterations = VX.WSpasses;   %% number of passes set by the variant (3 = original)']), ...
        1, 'number of passes', tag, logFile);
    txt = patchRegex(txt, ['N = 2\^\(nextpow2\(round\(CPI_s\*fsProc\)\)\);[^\r\n]*\r?\n[ \t]*win[ \t]*=[ \t]*hamming\(N\);[^\r\n]*\r?\n' ...
        '[ \t]*noverlap = N\*3/4;[^\r\n]*\r?\n[ \t]*nfft[ \t]*=[ \t]*N;[^\r\n]*'], ...
        sprintf(['global VX\nN = 2^(nextpow2(round(CPI_s*fsProc)));\nif VX.WSexact, N = round(CPI_s*fsProc); end   %% window not forced to a power of two\n' ...
        'win = hamming(N);\nnoverlap = N - max(1, round(N*VX.WSHop));   %% hop set by the variant\nnfft = 2^(nextpow2(N));   %% zero padded when N is not a power of two']), ...
        1, 'walking speed spectrogram block', tag, logFile);
    txt = patchRegex(txt, 'Ssmooth = SmoothSpectrogramNoise\(Power_S, 21\);[^\r\n]*', ...
        sprintf(['global VX\nsmoothBins_ = 21;\nif VX.WSexact\n    Nwin_ = round(CPI_s*fsProc);\n' ...
        '    smoothBins_ = 2*floor(21*(2^nextpow2(Nwin_))/Nwin_/2) + 1;   %% keeps the smoothing width in Hz when zero padding is used\nend\n' ...
        'Ssmooth = SmoothSpectrogramNoise(Power_S, smoothBins_);']), ...
        1, 'bin smoothing line', tag, logFile);
    txt = patchRegex(txt, 'freq_low\s*=\s*max\(0\.5\*centreFreq_Hz,freq_low0\);[^\r\n]*', ...
        sprintf(['if VX.WSbandA > 0\n    freq_low = max((1 - VX.WSbandA)*centreFreq_Hz, freq_low0);   %% narrower band below the mean tracked frequency\nelse\n' ...
        '    freq_low = max(0.5*centreFreq_Hz,freq_low0);\nend']), 1, 'band lower edge', tag, logFile);
    txt = patchRegex(txt, 'freq_high\s*=\s*max\(2\.0\*centreFreq_Hz,f_high_min\);[^\r\n]*', ...
        sprintf(['if VX.WSbandA > 0\n    freq_high = (1 + VX.WSbandA)*centreFreq_Hz;   %% narrower band above the mean tracked frequency\nelse\n' ...
        '    freq_high = max(2.0*centreFreq_Hz,f_high_min);\nend']), 1, 'band upper edge', tag, logFile);
end
if ~isempty(cfg.CopyHook)
    txt = cfg.CopyHook(txt, kind);
end
fid = fopen(dest, 'w');
assert(fid >= 0, 'Could not create %s', dest);
fprintf(fid, '%s\n', txt);
fclose(fid);
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

function m_o = runStep(m_P, m_x, m_D, m_v, m_snr, m_tc, m_lag)
global VX
VX = defaultVX();
VX.Hop = m_v.hop;
VX.CScale = m_v.cScale;
VX.CHop = m_v.cHop;
if ~isnan(m_lag)
    VX.LagSec = m_lag;
end
pthOverride = [];
switch m_v.pthMode
    case 'floor'
        pthOverride = max(10 - m_snr, m_v.floor);            % same rule as the extractor, with another floor
    case 'noisefloor'
        VX.PthMode = 'noisefloor';
        VX.PthMargin = m_v.margin;
end
m_o = struct('ok', false, 'det', [], 'key', [], 'sd', [], 'detRaw', [], 'sdAll', [], 'cpiAct', NaN, 'sc', [], 'scRaw', [], ...
    'Fs', [], 'tF', [], 'framePeriod', NaN, 'fhigh', NaN, 'lagFrames', NaN, 'Pth', NaN, 'TcUsed', NaN, 'cPeak', NaN, 'err', '');
try
    RdrRxSignal = m_x;                                       %#ok<NASGU>
    fs = m_D.fs;                                             %#ok<NASGU>
    Tc = m_tc;                                               %#ok<NASGU>
    CPI_s = m_v.cpiReq;                                      %#ok<NASGU>
    SNR_dB = m_snr;                                          %#ok<NASGU>
    cPeakScale = m_v.cScale;                                 %#ok<NASGU>
    stepRefineMode = m_v.refine;                             %#ok<NASGU>
    stepGapFill = m_v.gap;                                   %#ok<NASGU>
    if ~isempty(pthOverride)
        PthOverride = pthOverride;                           %#ok<NASGU>
    end
    if strcmp(m_v.ext, 'v5')
        evalc('run(m_P.st2);');
    else
        evalc('run(m_P.st);');
    end
    m_o.det = detectedStepTimes_s(:).';                      % final step times
    if exist('stepFrameTimes_s', 'var')
        m_o.key = stepFrameTimes_s(:).';                     % v5: frame time stamps decide who is inside the walk
    else
        m_o.key = m_o.det;
    end
    if exist('detectedStepTimes_raw_s', 'var')
        m_o.detRaw = detectedStepTimes_raw_s(:).';
    else
        m_o.detRaw = m_o.det;
    end
    m_o.cpiAct = CPI_actual_s;
    m_o.Fs = F_smooth(:).';
    m_o.tF = time_axis_step(:).';
    m_o.framePeriod = frame_period_s;
    m_o.fhigh = freq_high_Yardibi;
    m_o.lagFrames = lagFrames_;
    m_o.Pth = PthUsed_;
    m_o.TcUsed = Tc;
    m_o.cPeak = cPeak;
    edge = 0.1 + m_o.cpiAct/2;                               % same walk interior as the step time scoring
    lo = m_D.buffer + edge;
    hi = m_D.buffer + m_D.simT - edge;
    tIn = m_o.det(m_o.key > lo & m_o.key < hi);              % steps inside the walk only
    if numel(tIn) >= 3
        m_o.sd = (tIn(3:end) - tIn(1:end-2)) - m_D.Tc;       % stride durations inside the walk, kept as (stride - true Tc)
    end
    if numel(m_o.detRaw) >= 3
        m_o.sdAll = (m_o.detRaw(3:end) - m_o.detRaw(1:end-2)) - m_D.Tc;   % v3 rule: every detection, no trimming
    end
    m_o.sc = scoreStepTimes(m_o.det, m_D, m_o.cpiAct, m_o.key);
    m_o.scRaw = scoreStepTimes(m_o.detRaw, m_D, m_o.cpiAct, m_o.detRaw);
    m_o.ok = true;
catch ME
    m_o.err = ME.message;
end
VX = defaultVX();
end

function Nw = noiseWindow(cpi, fs, exact)
if exact
    Nw = round(cpi*fs);
else
    Nw = 2^nextpow2(round(cpi*fs));
end
end

function cpiAct = actualCPI(cpi, fs, exact)
if nargin < 3
    exact = false;
end
cpiAct = noiseWindow(cpi, fs, exact)/fs;
end

function selfTest(P, cfg, outDir)
global VX
res = {};
ok = true;
D = getSpeedData(P, cfg, 1.2);
nm = baseStepSpec(cfg);
snr = 60;
setSeed(1);
z = randn(1, numel(D.clean)) + 1i*randn(1, numel(D.clean));
x = addNoise(D.clean, noiseWindow(nm.cpiReq, D.fs, false), snr, z);
o0 = runStep(P, x, D, nm, snr, D.Tc, NaN);
[ok, res] = chk(ok, res, o0.ok, 'step extractor runs on the baseline spec', o0.err);
[ok, res] = chk(ok, res, abs(o0.fhigh - 350*24/5.8) < 1e-9, sprintf('frequency limit in the copy = %.4f Hz (350*24/5.8)', o0.fhigh), '');
[ok, res] = chk(ok, res, abs(o0.framePeriod - 0.25*o0.cpiAct) < 2/D.fs, sprintf('hop at the default = 0.25 of the window (%.4f s of %.4f s)', o0.framePeriod, o0.cpiAct), '');
v = nm; v.hop = 0.5;
o1 = runStep(P, x, D, v, snr, D.Tc, NaN);
[ok, res] = chk(ok, res, o1.ok && abs(o1.framePeriod - 0.5*o1.cpiAct) < 2/D.fs, sprintf('hop patch: 50 percent overlap gives frame period %.4f s (0.5 x %.4f)', o1.framePeriod, o1.cpiAct), o1.err);
v = nm; v.hop = 0.125;
o2 = runStep(P, x, D, v, snr, D.Tc, NaN);
[ok, res] = chk(ok, res, o2.ok && abs(o2.framePeriod - 0.125*o2.cpiAct) < 2/D.fs, sprintf('hop patch: 87.5 percent overlap gives frame period %.4f s', o2.framePeriod), o2.err);
o3 = runStep(P, x, D, nm, snr, D.Tc, o0.framePeriod);
[ok, res] = chk(ok, res, o3.ok && o3.lagFrames == 1 && isequal(o3.det, o0.det), 'lag patch: a lag of one frame reproduces the original peak rule exactly', o3.err);
o4 = runStep(P, x, D, nm, snr, D.Tc, 0.15);
[ok, res] = chk(ok, res, o4.ok && o4.lagFrames == max(1, round(0.15/o4.framePeriod)), sprintf('lag patch: 0.15 s = %d frames at this hop', o4.lagFrames), o4.err);
v = nm; v.hop = 0.125; v.cHop = true;
o5 = runStep(P, x, D, v, snr, D.Tc, NaN);
v.cHop = false;
o6 = runStep(P, x, D, v, snr, D.Tc, NaN);
[ok, res] = chk(ok, res, o5.ok && o6.ok && abs(o5.cPeak/o6.cPeak - 0.5) < 1e-9, 'c scaling patch: c halves at hop 0.125 only when asked', [o5.err o6.err]);
v = nm; v.pthMode = 'noisefloor'; v.margin = 10;
snr30 = 30;
x30 = addNoise(D.clean, noiseWindow(nm.cpiReq, D.fs, false), snr30, z);
o7 = runStep(P, x30, D, v, snr30, D.Tc, NaN);
o8 = runStep(P, x30, D, nm, snr30, D.Tc, NaN);
[ok, res] = chk(ok, res, o7.ok && o8.ok && isfinite(o7.Pth) && abs(o7.Pth - o8.Pth) > 1e-6, ...
    sprintf('Pth patch: noise-floor Pth = %.1f dB, coupled Pth = %.1f dB at Gamma = 30 dB', o7.Pth, o8.Pth), [o7.err o8.err]);
v = nm; v.pthMode = 'floor'; v.floor = -30;
o9 = runStep(P, x, D, v, snr, D.Tc, NaN);
[ok, res] = chk(ok, res, o9.ok && abs(o9.Pth - (-30)) < 1e-9, sprintf('Pth floor override: Pth = %.1f dB at Gamma = 60 dB with floor -30 (coupled would be -40)', o9.Pth), o9.err);
o10 = runStep(P, x, D, nm, snr, 1.2, NaN);
[ok, res] = chk(ok, res, o10.ok && abs(o10.TcUsed - 1.2) < 1e-12 && abs(D.Tc - 1.2) > 1e-6, sprintf('Tc handed to the extractor = %.3f s while the true Tc = %.3f s', o10.TcUsed, D.Tc), o10.err);
v = nm; v.ext = 'v5';
o11 = runStep(P, x, D, v, snr, D.Tc, NaN);
[ok, res] = chk(ok, res, o11.ok && abs(o11.fhigh - 350*24/5.8) < 1e-9, 'v5 copy runs and also uses 350 Hz', o11.err);
ws0 = struct('cpi', 0.2, 'exact', false, 'bandA', 0, 'passes', 3, 'hop', 0.25);
xs = addNoise(D.clean, noiseWindow(0.2, D.fs, false), 40, z);
w0 = runSpeed(P, xs, D, ws0);
[ok, res] = chk(ok, res, w0.ok && abs(w0.cpiAct - 16384/D.fs) < 1e-9 && w0.nIter == 3, sprintf('walking speed default: window %.4f s, %d passes, speed %.3f m/s', w0.cpiAct, w0.nIter, w0.meanSpeed), w0.err);
ws1 = ws0; ws1.cpi = 0.3; ws1.exact = true;
w1 = runSpeed(P, xs, D, ws1);
[ok, res] = chk(ok, res, w1.ok && abs(w1.cpiAct - round(0.3*D.fs)/D.fs) < 1e-9, sprintf('walking speed exact window: %.4f s for 0.3 s asked', w1.cpiAct), w1.err);
ws2 = ws0; ws2.bandA = 0.25;
w2 = runSpeed(P, xs, D, ws2);
[ok, res] = chk(ok, res, w2.ok && abs(w2.bandHigh - 1.25*w2.centrePrev) < 1e-9, sprintf('band patch: final upper edge %.2f Hz = 1.25 x %.2f Hz', w2.bandHigh, w2.centrePrev), w2.err);
ws3 = ws0; ws3.passes = 5;
w3 = runSpeed(P, xs, D, ws3);
[ok, res] = chk(ok, res, w3.ok && w3.nIter == 5, sprintf('passes patch: %d passes done', w3.nIter), w3.err);
ws4 = ws0; ws4.hop = 0.5;
w4 = runSpeed(P, xs, D, ws4);
[ok, res] = chk(ok, res, w4.ok && abs(w4.framePeriod - 0.5*w4.cpiAct) < 2/D.fs, sprintf('walking speed hop patch: frame period %.4f s', w4.framePeriod), w4.err);
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

function [ok, res] = chk(ok, res, cond, text, errText)
if cond
    res{end+1} = ['PASS  ' text];
else
    res{end+1} = ['FAIL  ' text ' ' errText];
    ok = false;
end
end

function v = baseStepSpec(cfg)
v = struct('name', 'BASE', 'exp', 'E0', 'ext', 'main', 'cpiReq', cfg.StepCPIreq, 'hop', 0.25, 'cScale', 0.5, 'cHop', false, ...
    'lagMode', 'frame', 'lagK', cfg.E2LagK, 'lagSec', NaN, 'tcMode', 'true', 'pthMode', 'coupled', 'floor', -40, 'margin', 10, ...
    'refine', cfg.StepRefine, 'gap', cfg.StepGapFill, 'seedTag', 'A');
end

function w = baseSpeedSpec(cfg)
w = struct('name', 'BASE', 'exp', 'E0', 'cpi', cfg.WSCPIreq, 'exact', false, 'bandA', 0, 'passes', 3, 'hop', 0.25, 'kAdapt', 0, 'seedTag', 'A');
end

function s = seedOf(cfg, tag)
switch tag
    case 'A', s = cfg.SeedA;
    case 'B', s = cfg.SeedB;
    otherwise, s = cfg.SeedVal;
end
end

function key = wKey(w, snr, seed, nT)
key = sprintf('%.5g|%d|%g|%d|%g|%g|%g|%d|%d', w.cpi, w.exact, w.bandA, w.passes, w.hop, w.kAdapt, snr, seed, nT);
end

function W = getW(P, cfg, D, w, snr, nT)
global WCACHE
seed = seedOf(cfg, w.seedTag);
key = wKey(w, snr, seed, nT);
if isKey(WCACHE, key)
    W = WCACHE(key);
    return;
end
W = cell(nT, 1);
for k = 1:nT
    z = trialNoise(D, seed, k);
    if w.kAdapt > 0
        Nw = noiseWindow(0.372, D.fs, true);
        x = addNoise(D.clean, Nw, snr, z);
        coarse = struct('cpi', 0.372, 'exact', true, 'bandA', 0, 'passes', 3, 'hop', 0.25);
        o0 = runSpeed(P, x, D, coarse);
        if o0.ok && o0.meanSpeed > 0.05
            cpiK = w.kAdapt/(160*o0.meanSpeed);
        else
            cpiK = 0.372;
        end
        cpiK = min(max(cpiK, min(cfg.E4CPI)), max(cfg.E4CPI));   % kept inside the tested range
        W{k} = runSpeed(P, x, D, struct('cpi', cpiK, 'exact', true, 'bandA', 0, 'passes', 3, 'hop', 0.25));
    else
        Nw = noiseWindow(w.cpi, D.fs, w.exact);
        x = addNoise(D.clean, Nw, snr, z);
        W{k} = runSpeed(P, x, D, struct('cpi', w.cpi, 'exact', w.exact, 'bandA', w.bandA, 'passes', w.passes, 'hop', w.hop));
    end
end
WCACHE(key) = W;
end

function prep = getPrep(P, cfg, D, snr, nT, seedTag)
global PCACHE
seed = seedOf(cfg, seedTag);
key = sprintf('%g|%d|%d|%g', snr, seed, nT, cfg.StepCPIreq);
if isKey(PCACHE, key)
    prep = PCACHE(key);
    return;
end
v = baseStepSpec(cfg);
prep.acorr = NaN(nT, 1);
prep.acorr2 = NaN(nT, 1);
prep.twopass = NaN(nT, 1);
for k = 1:nT
    z = trialNoise(D, seed, k);
    x = addNoise(D.clean, noiseWindow(v.cpiReq, D.fs, false), snr, z);
    o = runStep(P, x, D, v, snr, cfg.TcNominal, NaN);
    if ~o.ok
        continue;
    end
    prep.acorr(k) = tcAutocorr(o.Fs, o.tF, D, 0.6, 3.0, false);
    prep.acorr2(k) = tcAutocorr(o.Fs, o.tF, D, 0.3, 1.5, true);
    prep.twopass(k) = tcTwoPass(o, D);
end
PCACHE(key) = prep;
end

function T = tcAutocorr(F, tF, D, lagLo, lagHi, doubleIt)
T = NaN;
inWalk = tF >= D.buffer & tF <= D.buffer + D.simT;           % frames inside the walk
y = F(inWalk);
y = y(:).';
if numel(y) < 8
    return;
end
y(isnan(y)) = mean(y(~isnan(y)));                            % fill gaps with the mean
fp = median(diff(tF));                                       % frame period
y = y - mean(y);
s = sort(y);
n = numel(y);
lo = s(max(1, round(0.05*n)));
hi = s(min(n, round(0.95*n)));
y = min(max(y, lo), hi);                                     % clip the extremes
y = movmean(y, 3);                                           % light smoothing
r = real(ifft(abs(fft(y, 2*n)).^2));                         % autocorrelation through the FFT
r = r(1:n) ./ (n:-1:1);                                      % unbiased: divide by the number of overlapping samples
if r(1) <= 0
    return;
end
r = r / r(1);
lags = (0:n-1)*fp;
idx = find(lags >= lagLo & lags <= min(lagHi, 0.5*n*fp));    % allowed lags, never more than half the record
cand = [];
for i = idx
    if i > 1 && i < n && r(i) >= r(i-1) && r(i) >= r(i+1)    % local maximum
        cand(end+1) = i;                                     %#ok<AGROW>
    end
end
if isempty(cand)
    return;
end
[~, m] = max(r(cand));
T = lags(cand(m));
if doubleIt
    T = 2*T;                                                 % the strongest peak was taken to be the step period
end
end

function T = tcTwoPass(o, D)
T = NaN;
edge = 0.1 + o.cpiAct/2;
inWalk = o.key > D.buffer + edge & o.key < D.buffer + D.simT - edge;
t = o.det(inWalk);
if numel(t) >= 3
    T = 2*median(diff(t));
end
end

function [tcX, tcEst, failed] = pickTc(v, D, cfg, prep, k)
tcEst = NaN;
failed = false;
switch v.tcMode
    case 'true'
        tcX = D.Tc;
    case 'nominal'
        tcX = cfg.TcNominal;
    otherwise
        tcEst = prep.(v.tcMode)(k);
        if isfinite(tcEst) && tcEst > 0.3 && tcEst < 4
            tcX = tcEst;
        else
            tcX = cfg.TcNominal;                             % estimator failed, fall back to the nominal value
            failed = true;
        end
end
end

function rows = evalStep(P, cfg, D, v, snr, nT, wSpec)
seed = seedOf(cfg, v.seedTag);
Wk = getW(P, cfg, D, wSpec, snr, nT);
needPrep = any(strcmp(v.tcMode, {'acorr', 'acorr2', 'twopass'})) || strcmp(v.lagMode, 'est');
prep = [];
if needPrep
    prep = getPrep(P, cfg, D, snr, nT, v.seedTag);
end
S = cell(nT, 1);
tcEstV = NaN(nT, 1);
tcFail = 0;
lagEff = NaN(nT, 1);
pthV = NaN(nT, 1);
for k = 1:nT
    z = trialNoise(D, seed, k);
    x = addNoise(D.clean, noiseWindow(v.cpiReq, D.fs, false), snr, z);
    [tcX, tcEstV(k), failed] = pickTc(v, D, cfg, prep, k);
    tcFail = tcFail + failed;
    lagX = NaN;
    switch v.lagMode
        case 'oracle', lagX = v.lagK*D.Tc;
        case 'est'
            te = prep.(cfg.E2TcEstimator)(k);
            if ~(isfinite(te) && te > 0.3 && te < 4), te = cfg.TcNominal; end
            lagX = v.lagK*te;
        case 'fixed', lagX = v.lagSec;
    end
    S{k} = runStep(P, x, D, v, snr, tcX, lagX);
    if S{k}.ok
        lagEff(k) = S{k}.lagFrames*S{k}.framePeriod;         % the lag the rule really used, in seconds
        pthV(k) = S{k}.Pth;
    end
end
rows = [];
params = {'StepTime', 'StrideDuration', 'StrideDuration_v3rule', 'StrideLength'};
for p = 1:numel(params)
    E = cell(nT, 1);
    G = cell(nT, 1);
    miss = NaN(nT, 1);
    extra = NaN(nT, 1);
    missC = NaN(nT, 1);
    extraC = NaN(nT, 1);
    perfect = NaN(nT, 1);
    off = NaN(nT, 1);
    for k = 1:nT
        E{k} = NaN;
        G{k} = NaN;
        if ~S{k}.ok
            continue;
        end
        scK = S{k}.sc;
        switch params{p}
            case 'StepTime'
                E{k} = scK.e;
                G{k} = scK.g;
            case 'StrideDuration'
                E{k} = S{k}.sd;
                G{k} = D.Tc * ones(size(S{k}.sd));
            case 'StrideDuration_v3rule'
                E{k} = S{k}.sdAll;
                G{k} = D.Tc * ones(size(S{k}.sdAll));
            case 'StrideLength'
                [E{k}, G{k}] = strideLengthErrors(Wk{k}, S{k}, D);
        end
        miss(k) = scK.miss;
        extra(k) = scK.extra;
        missC(k) = scK.missC;
        extraC(k) = scK.extraC;
        perfect(k) = scK.perfect;
        off(k) = scK.offset;
    end
    row = newRow(cfg, D, v.exp, v.name, v.ext, params{p}, snr, seed);
    row.CPI_step_req = v.cpiReq;
    row.CPI_step_act = actualCPI(v.cpiReq, D.fs);
    if strcmp(params{p}, 'StrideLength')
        row.CPI_speed_req = wSpec.cpi;
        row.CPI_speed_act = actualCPI(wSpec.cpi, D.fs, wSpec.exact);
    end
    row.HopFraction = v.hop;
    row.Overlap = 1 - v.hop;
    row.CPeakScale = v.cScale;
    row.LagSeconds = meanFinite(lagEff);
    row.LagSource = v.lagMode;
    row.TcSource = v.tcMode;
    row.PthMode = v.pthMode;
    row.Pth_dB = meanFinite(pthV);
    if strcmp(v.pthMode, 'floor') || strcmp(v.pthMode, 'coupled')
        row.Floor_dB = v.floor;
    end
    if strcmp(v.pthMode, 'noisefloor')
        row.Margin_dB = v.margin;
    end
    row = fillStats(row, E, G, strcmp(params{p}, 'StrideLength'), cfg);
    if ~strcmp(params{p}, 'StrideLength')
        row.MissRate_pct = meanFinite(miss);
        row.ExtraRate_pct = meanFinite(extra);
        row.MissCentred_pct = meanFinite(missC);
        row.ExtraCentred_pct = meanFinite(extraC);
        row.PerfectTrials_pct = 100*meanFinite(perfect);
        row.StepOffset_ms = 1000*meanFinite(off);
    end
    if any(isfinite(tcEstV))
        row.TcTrue_s = D.Tc;
        row.TcEst_s = meanFinite(tcEstV);
        row.TcErr_pct = 100*(meanFinite(tcEstV) - D.Tc)/D.Tc;
        row.TcAbsErr_pct = 100*meanFinite(abs(tcEstV - D.Tc))/D.Tc;
        row.TcEstStd_s = stdFinite(tcEstV);
        row.TcFailed = tcFail;
    end
    rows = addRow(rows, row);
end
end

function [row, W] = evalSpeed(P, cfg, D, w, snr, nT)
seed = seedOf(cfg, w.seedTag);
W = getW(P, cfg, D, w, snr, nT);
E = cell(nT, 1);
G = cell(nT, 1);
cpiV = NaN(nT, 1);
for k = 1:nT
    if W{k}.ok
        E{k} = W{k}.meanSpeed - D.v0;
        G{k} = D.v0;
        cpiV(k) = W{k}.cpiAct;
    else
        E{k} = NaN;
        G{k} = NaN;
    end
end
row = newRow(cfg, D, w.exp, w.name, 'ws', 'WalkingSpeed', snr, seed);
if w.kAdapt > 0
    row.CPI_speed_req = NaN;
    row.AdaptK = w.kAdapt;
else
    row.CPI_speed_req = w.cpi;
end
row.CPI_speed_act = meanFinite(cpiV);
row.CyclesInWindow = meanFinite(cpiV)*160*D.v0;               % cycles of the torso Doppler (f = 160*v Hz) inside one window
row.HopFraction = w.hop;
row.Overlap = 1 - w.hop;
row.BandA = w.bandA;
row.Passes = w.passes;
row = fillStats(row, E, G, false, cfg);
end

function row = newRow(cfg, D, exp, variant, ext, param, snr, seed)
row = blankRow();
row.Stage = exp;
row.Experiment = exp;
row.Variant = variant;
row.Extractor = ext;
row.Parameter = param;
row.Speed_mps = D.v0;
row.SNR_dB = snr;
row.NoiseMode = 'extractorN';
row.Seed = seed;
end

function rows = runE0(P, cfg, D)
nT = cfg.nTrials;
rows = [];
for tag = {'A', 'B'}
    w = baseSpeedSpec(cfg);
    w.seedTag = tag{1};
    if strcmp(tag{1}, 'B'), w.name = 'BASE_B'; end
    for snr = cfg.AllSNRs
        rows = addRow(rows, evalSpeed(P, cfg, D, w, snr, nT));
    end
    for ext = {'main', 'v5'}
        v = baseStepSpec(cfg);
        v.ext = ext{1};
        v.seedTag = tag{1};
        if strcmp(tag{1}, 'B'), v.name = 'BASE_B'; end
        for snr = cfg.AllSNRs
            rows = [rows, evalStep(P, cfg, D, v, snr, nT, baseSpeedSpec(cfg))];   %#ok<AGROW>
        end
    end
end
end

function rows = runE1(P, cfg, D)
nT = cfg.nTrials;
rows = [];
exts = {'main'};
if cfg.E1AlsoV5, exts{end+1} = 'v5'; end
for e = 1:numel(exts)
    for ov = cfg.E1Overlaps
        for cs = cfg.E1CScales
            v = baseStepSpec(cfg);
            v.exp = 'E1';
            v.ext = exts{e};
            v.hop = 1 - ov;
            v.cScale = cs;
            v.name = sprintf('ov%g_c%g', 100*ov, cs);
            for snr = cfg.StepSNRs
                rows = [rows, evalStep(P, cfg, D, v, snr, nT, baseSpeedSpec(cfg))];   %#ok<AGROW>
            end
        end
    end
end
for ov = cfg.E1Overlaps                                       % the walking-speed algorithm also uses 75 % overlap
    w = baseSpeedSpec(cfg);
    w.exp = 'E1';
    w.hop = 1 - ov;
    w.name = sprintf('ws_ov%g', 100*ov);
    for snr = cfg.SpeedSNRs
        rows = addRow(rows, evalSpeed(P, cfg, D, w, snr, nT));
    end
end
end

function rows = runE2(P, cfg, D)
nT = cfg.nTrials;
rows = [];
for hop = cfg.E2Hops
    rules = {'frame', 'frame_cscaled', 'oracle', 'est'};
    for q = 1:numel(cfg.E2FixedLags)
        rules{end+1} = sprintf('fix%g', cfg.E2FixedLags(q));  %#ok<AGROW>
    end
    for r = 1:numel(rules)
        v = baseStepSpec(cfg);
        v.exp = 'E2';
        v.hop = hop;
        v.name = sprintf('hop%g_%s', hop, rules{r});
        switch rules{r}
            case 'frame'
                v.lagMode = 'frame';
            case 'frame_cscaled'
                v.lagMode = 'frame';
                v.cHop = true;                                % c shrinks with the hop, as in Master_SimSweep_v3 Stage 3
            case 'oracle'
                v.lagMode = 'oracle';
            case 'est'
                v.lagMode = 'est';
            otherwise
                v.lagMode = 'fixed';
                v.lagSec = cfg.E2FixedLags(r - 4);
        end
        for snr = cfg.StepSNRs
            rows = [rows, evalStep(P, cfg, D, v, snr, nT, baseSpeedSpec(cfg))];   %#ok<AGROW>
        end
    end
end
end

function rows = runE3(P, cfg, D)
nT = cfg.nTrials;
rows = [];
modes = {'true', 'nominal', 'acorr', 'acorr2', 'twopass'};
for m = 1:numel(modes)
    v = baseStepSpec(cfg);
    v.exp = 'E3';
    v.tcMode = modes{m};
    v.name = ['Tc_' modes{m}];
    for snr = cfg.StepSNRs
        rows = [rows, evalStep(P, cfg, D, v, snr, nT, baseSpeedSpec(cfg))];   %#ok<AGROW>
    end
end
end

function rows = runE5(P, cfg, D)
nT = cfg.nTrials;
rows = [];
for a = [0, cfg.E5Bands]
    for np = cfg.E5Passes
        if a == 0 && np == 3
            continue;                                         % that is the baseline
        end
        w = baseSpeedSpec(cfg);
        w.exp = 'E5';
        w.cpi = cfg.E5CPIreq;
        w.bandA = a;
        w.passes = np;
        w.name = sprintf('a%g_p%d', a, np);
        for snr = cfg.SpeedSNRs
            rows = addRow(rows, evalSpeed(P, cfg, D, w, snr, nT));
        end
    end
end
end

function rows = runE6(P, cfg, D)
nT = cfg.nTrials;
rows = [];
for f = cfg.E6Floors
    v = baseStepSpec(cfg);
    v.exp = 'E6';
    v.pthMode = 'floor';
    v.floor = f;
    v.name = sprintf('floor%g', f);
    for snr = cfg.StepSNRs
        rows = [rows, evalStep(P, cfg, D, v, snr, nT, baseSpeedSpec(cfg))];   %#ok<AGROW>
    end
end
for m = cfg.E6Margins
    v = baseStepSpec(cfg);
    v.exp = 'E6';
    v.pthMode = 'noisefloor';
    v.margin = m;
    v.name = sprintf('noisefloor+%g', m);
    for snr = cfg.StepSNRs
        rows = [rows, evalStep(P, cfg, D, v, snr, nT, baseSpeedSpec(cfg))];   %#ok<AGROW>
    end
end
end

function cpis = validE4CPIs(cfg, D)
cpis = [];
for c = cfg.E4CPI
    nFrames = floor((D.simT - c)/(0.25*c)) + 1;               % frames at 75 percent overlap that fit in the walk
    if nFrames >= cfg.E4MinFrames
        cpis(end+1) = c;                                      %#ok<AGROW>
    end
end
end

function rows = runE4a(P, cfg, D)
nT = cfg.nTrials;
rows = [];
cpis = validE4CPIs(cfg, D);
if numel(cpis) < numel(cfg.E4CPI)
    logMsg(cfg.OutDir, '   E4a: windows dropped at %.1f m/s (fewer than %d frames in the %.1f s walk): %s', D.v0, cfg.E4MinFrames, D.simT, mat2str(setdiff(cfg.E4CPI, cpis)));
end
for snr = cfg.SpeedSNRs
    for c = cpis
        w = baseSpeedSpec(cfg);
        w.exp = 'E4a';
        w.cpi = c;
        w.exact = true;
        w.name = sprintf('cpi%.3f', c);
        rows = addRow(rows, evalSpeed(P, cfg, D, w, snr, nT));
    end
    for cr = cfg.E4StepCPIreq
        v = baseStepSpec(cfg);
        v.exp = 'E4a';
        v.cpiReq = cr;
        v.name = sprintf('stepcpi%g', cr);
        rows = [rows, evalStep(P, cfg, D, v, snr, nT, baseSpeedSpec(cfg))];   %#ok<AGROW>
    end
end
end

function sel = selectE4(cfg, outDir)
rows = loadExp(outDir, 'E4a', cfg.Speeds);
assert(~isempty(rows), 'No E4a results to select from.');
sel = struct();
txt = '';
for snr = cfg.SpeedSNRs
    key = sprintf('s%g', snr);
    wsRows = rows(strcmp({rows.Parameter}, 'WalkingSpeed') & [rows.SNR_dB] == snr);
    cp = unique(round([wsRows.CPI_speed_req]*1000)/1000);
    passes = zeros(size(cp));
    worst = zeros(size(cp));
    complete = true(size(cp));
    for q = 1:numel(cp)
        for s = 1:numel(cfg.Speeds)
            r = wsRows(abs([wsRows.Speed_mps] - cfg.Speeds(s)) < 1e-9 & abs(round([wsRows.CPI_speed_req]*1000)/1000 - cp(q)) < 1e-9);
            if isempty(r)
                complete(q) = false;                          % this window was dropped at this speed
                continue;
            end
            passes(q) = passes(q) + (r(1).RMSE_pct < cfg.Target_pct);
            if cfg.Speeds(s) >= 0.5
                worst(q) = max(worst(q), r(1).RMSE_pct);
            end
        end
    end
    passes(~complete) = -1;
    best = find(passes == max(passes));
    [~, j] = min(worst(best) + 1e-6*cp(best));                % fewest failures, then lowest worst case, then shortest
    sel.(key).fixedCPI = cp(best(j));
    sel.(key).fixedPasses = passes(best(j));
    sel.(key).wsCPI = NaN(1, numel(cfg.Speeds));
    sel.(key).stCPI = NaN(1, numel(cfg.Speeds));
    for s = 1:numel(cfg.Speeds)
        r = wsRows(abs([wsRows.Speed_mps] - cfg.Speeds(s)) < 1e-9);
        if ~isempty(r)
            sel.(key).wsCPI(s) = r(pickTie([r.RMSE_pct], [r.CPI_speed_req], 0.25)).CPI_speed_req;
        end
        st = rows(strcmp({rows.Parameter}, 'StepTime') & [rows.SNR_dB] == snr & abs([rows.Speed_mps] - cfg.Speeds(s)) < 1e-9);
        if ~isempty(st)
            sel.(key).stCPI(s) = st(pickTie([st.RMSE_pct], [st.CPI_step_req], 0.25)).CPI_step_req;
        end
    end
    txt = [txt sprintf('Gamma %g dB: fixed speed window %.3f s (passes at %d speeds in the sweep); ', snr, sel.(key).fixedCPI, sel.(key).fixedPasses)]; %#ok<AGROW>
end
sel.text = txt;
end

function rows = runE4v(P, cfg, D, sel4)
nT = cfg.nTrials;
rows = [];
si = find(abs(cfg.Speeds - D.v0) < 1e-9, 1);
for snr = cfg.SpeedSNRs
    S = sel4.(sprintf('s%g', snr));
    cfgNames = {'baseline', 'fixed', 'perspeed'};
    for K = cfg.E4K
        cfgNames{end+1} = sprintf('adaptK%g', K);              %#ok<AGROW>
    end
    for q = 1:numel(cfgNames)
        w = baseSpeedSpec(cfg);
        w.exp = 'E4v';
        w.seedTag = 'V';
        v = baseStepSpec(cfg);
        v.exp = 'E4v';
        v.seedTag = 'V';
        switch cfgNames{q}
            case 'baseline'
            case 'fixed'
                w.cpi = S.fixedCPI;
                w.exact = true;
            case 'perspeed'
                if isfinite(S.wsCPI(si))
                    w.cpi = S.wsCPI(si);
                    w.exact = true;
                end
                if isfinite(S.stCPI(si))
                    v.cpiReq = S.stCPI(si);
                end
            otherwise
                w.kAdapt = cfg.E4K(q - 3);
        end
        w.name = cfgNames{q};
        v.name = cfgNames{q};
        rows = addRow(rows, evalSpeed(P, cfg, D, w, snr, nT));
        rows = [rows, evalStep(P, cfg, D, v, snr, nT, w)];    %#ok<AGROW>
    end
end
end

function rows = loadExp(outDir, e, speeds)
rows = [];
for s = 1:numel(speeds)
    f = expFile(outDir, e, speeds(s));
    if isfile(f)
        m = load(f);
        rows = [rows, m.rows];                                %#ok<AGROW>
    end
end
end

function idx = pickTie(vals, keys, tol)
vals = vals(:).';
keys = keys(:).';
valid = isfinite(vals);
if ~any(valid)
    idx = 1;
    return;
end
m = min(vals(valid));
cand = find(valid & vals <= m + tol);
[~, j] = min(keys(cand));
idx = cand(j);
end

function T = readCsvStruct(file)
txt = fileread(file);
lines = regexp(txt, '\r?\n', 'split');
lines = lines(~cellfun(@isempty, lines));
hdr = strsplit(lines{1}, ',');
n = numel(lines) - 1;
C = cell(n, numel(hdr));
for i = 1:n
    parts = strsplit(strrep(lines{i+1}, '"', ''), ',');
    m = min(numel(parts), numel(hdr));
    C(i, 1:m) = parts(1:m);
end
textCols = {'Stage', 'Parameter', 'NoiseMode'};
for j = 1:numel(hdr)
    col = C(:, j);
    if any(strcmp(hdr{j}, textCols))
        T.(hdr{j}) = col;
    else
        T.(hdr{j}) = str2double(col);
    end
end
end

function tabs = e7Tables(cfg)
persistent cacheTabs cacheDir
if ~isempty(cacheTabs) && strcmp(cacheDir, cfg.E7RunDir)
    tabs = cacheTabs;
    return;
end
tabs = [];
d = cfg.E7RunDir;
if ~exist(d, 'dir')
    d = fullfile(cfg.CodeFolder, cfg.E7RunDir);
end
f1 = fullfile(d, 'Stage1_CPI_sweep.csv');
f2 = fullfile(d, 'Stage2_SNR_sweep.csv');
if ~isfile(f1) || ~isfile(f2)
    return;
end
tabs.T1 = readCsvStruct(f1);
tabs.T2 = readCsvStruct(f2);
cacheTabs = tabs;
cacheDir = cfg.E7RunDir;
end

function sel = e7Select(tabs, v0, rule, cfg)
params = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideLength'};
sel = struct('param', params, 'cpiS', NaN, 'cpiT', NaN, 'snr', NaN, 'selRMSE', NaN, 'fallback', false, 'noS2', false);
T1 = tabs.T1;
T2 = tabs.T2;
for p = 1:4
    m1 = strcmp(T1.Parameter, params{p}) & abs(T1.Speed_mps - v0) < 1e-9;
    if ~any(m1)
        continue;
    end
    r1 = T1.RMSE_pct(m1);
    cs = T1.CPI_speed_req(m1);
    ct = T1.CPI_step_req(m1);
    okMask = true(size(r1));
    if strcmp(params{p}, 'StepTime')
        me = T1.MissRate_pct(m1) + T1.ExtraRate_pct(m1);
        if any(me <= cfg.MaxMissExtra_pct)
            okMask = me <= cfg.MaxMissExtra_pct;              % only windows that keep misses and extras acceptable, as in the v3 rule
        end
    end
    useCS = strcmp(params{p}, 'WalkingSpeed') || strcmp(params{p}, 'StrideLength');
    useCT = ~strcmp(params{p}, 'WalkingSpeed');
    key1 = zeros(size(r1));
    if useCS, key1 = key1 + cs; end
    if useCT, key1 = key1 + ct; end
    i1 = NaN;
    if strcmp(rule.kind, 'plateau')
        i1 = plateauPick(r1, cs, ct, useCS, useCT, cfg.Target_pct);
        if isnan(i1)
            sel(p).fallback = true;
        end
    end
    if isnan(i1)
        cand = find(okMask);
        j = pickTie(r1(cand), key1(cand), rule.tol);
        if strcmp(rule.kind, 'plateau')
            j = pickTie(r1(cand), key1(cand), 0.25);          % no plateau found: fall back to the baseline rule
        end
        i1 = cand(j);
    end
    if useCS, sel(p).cpiS = cs(i1); end
    if useCT, sel(p).cpiT = ct(i1); end
    m2 = strcmp(T2.Parameter, params{p}) & abs(T2.Speed_mps - v0) < 1e-9;
    if ~any(m2)
        sel(p).noS2 = true;
        sel(p).snr = cfg.StepSNRs(1);
        sel(p).selRMSE = r1(i1);
        continue;
    end
    snrs = T2.SNR_dB(m2);
    r2 = T2.RMSE_pct(m2);
    cs2 = T2.CPI_speed_req(m2);
    ct2 = T2.CPI_step_req(m2);
    if (useCS && abs(cs2(1) - sel(p).cpiS) > 1e-9) || (useCT && abs(ct2(1) - sel(p).cpiT) > 1e-9)
        sel(p).noS2 = true;                                   % the curve in the CSV belongs to another window
    end
    [snrs, order] = sort(snrs);
    r2 = r2(order);
    i2 = NaN;
    if strcmp(rule.kind, 'plateau')
        i2 = plateauPick(r2, snrs, zeros(size(snrs)), true, false, cfg.Target_pct);
        if isnan(i2)
            sel(p).fallback = true;
        end
    end
    if isnan(i2)
        tol = rule.tol;
        if strcmp(rule.kind, 'plateau'), tol = 0.25; end
        i2 = pickTie(r2, snrs, tol);
    end
    sel(p).snr = snrs(i2);
    sel(p).selRMSE = r2(i2);
end
end

function i = plateauPick(r, a, b, useA, useB, target)
n = numel(r);
if useA, ua = unique(a); else, ua = 0; a = zeros(size(a)); end
if useB, ub = unique(b); else, ub = 0; b = zeros(size(b)); end
score = inf(n, 1);
for k = 1:n
    ia = find(abs(ua - a(k)) < 1e-12, 1);
    ib = find(abs(ub - b(k)) < 1e-12, 1);
    nb = [];
    if ia > 1, nb = [nb, find(abs(a - ua(ia-1)) < 1e-12 & abs(b - b(k)) < 1e-12, 1)]; end          % neighbour with the next shorter first setting
    if ia < numel(ua), nb = [nb, find(abs(a - ua(ia+1)) < 1e-12 & abs(b - b(k)) < 1e-12, 1)]; end  % neighbour with the next longer first setting
    if ib > 1, nb = [nb, find(abs(b - ub(ib-1)) < 1e-12 & abs(a - a(k)) < 1e-12, 1)]; end          % neighbour with the next shorter second setting
    if ib < numel(ub), nb = [nb, find(abs(b - ub(ib+1)) < 1e-12 & abs(a - a(k)) < 1e-12, 1)]; end  % neighbour with the next longer second setting
    vals = [r(k); r(nb)];
    if all(isfinite(vals)) && all(vals < target)
        score(k) = mean(vals);
    end
end
[m, i] = min(score);
if ~isfinite(m)
    i = NaN;
end
end

function rows = runE7(P, cfg, D)
rows = [];
tabs = e7Tables(cfg);
if isempty(tabs)
    logMsg(cfg.OutDir, '   E7 skipped: Stage1_CPI_sweep.csv and Stage2_SNR_sweep.csv not found in %s. Set cfg.E7RunDir.', cfg.E7RunDir);
    return;
end
nT = cfg.nTrials;
rules = {};
for tol = cfg.E7TieTols
    rules{end+1} = struct('name', sprintf('tie%g', tol), 'kind', 'tie', 'tol', tol);   %#ok<AGROW>
end
rules{end+1} = struct('name', 'plateau', 'kind', 'plateau', 'tol', 0.25);
done = containers.Map();
for q = 1:numel(rules)
    sel = e7Select(tabs, D.v0, rules{q}, cfg);
    for p = 1:4
        s = sel(p);
        if ~isfinite(s.snr)
            continue;
        end
        k = sprintf('%s|%g|%g|%g', s.param, s.cpiS, s.cpiT, s.snr);
        if isKey(done, k)
            row = done(k);                                    % two rules chose the same setting, the validation is identical
        else
            row = validateSetting(P, cfg, D, s, nT);
            done(k) = row;
        end
        if isempty(row), continue; end
        row.Experiment = 'E7';
        row.Stage = 'E7';
        row.Variant = rules{q}.name;
        row.SelRMSE_pct = s.selRMSE;
        row.SelCPI_s = max([s.cpiS, s.cpiT]);
        row.SelSNR_dB = s.snr;
        row.SelFallback = double(s.fallback);
        row.SelNoStage2 = double(s.noS2);
        row.GapPct = row.RMSE_pct - s.selRMSE;
        rows = addRow(rows, row);
    end
end
end

function row = validateSetting(P, cfg, D, s, nT)
row = [];
switch s.param
    case 'WalkingSpeed'
        w = baseSpeedSpec(cfg);
        w.seedTag = 'V';
        w.cpi = s.cpiS;
        row = evalSpeed(P, cfg, D, w, s.snr, nT);
    case {'StepTime', 'StrideDuration'}
        v = baseStepSpec(cfg);
        v.seedTag = 'V';
        v.cpiReq = s.cpiT;
        w = baseSpeedSpec(cfg);
        w.seedTag = 'V';
        r = evalStep(P, cfg, D, v, s.snr, nT, w);
        row = r(strcmp({r.Parameter}, s.param));
    case 'StrideLength'
        v = baseStepSpec(cfg);
        v.seedTag = 'V';
        v.cpiReq = s.cpiT;
        w = baseSpeedSpec(cfg);
        w.seedTag = 'V';
        w.cpi = s.cpiS;
        r = evalStep(P, cfg, D, v, s.snr, nT, w);
        row = r(strcmp({r.Parameter}, 'StrideLength'));
end
end

function row = blankRow()
row = struct('Stage', '', 'Parameter', '', 'Speed_mps', NaN, ...
    'CPI_speed_req', NaN, 'CPI_speed_act', NaN, 'CPI_step_req', NaN, 'CPI_step_act', NaN, ...
    'SNR_dB', NaN, 'Pth_dB', NaN, 'HopFraction', NaN, 'CPeakScale', NaN, 'NoiseMode', '', ...
    'Trials', NaN, 'FailedTrials', NaN, 'N_compared', NaN, 'TruthMean', NaN, ...
    'RMSE_abs', NaN, 'RMSE_pct', NaN, 'RMSEpooled_pct', NaN, 'MeanErr_pct', NaN, ...
    'TrialMean_pct', NaN, 'TrialStd_pct', NaN, ...
    'MissRate_pct', NaN, 'ExtraRate_pct', NaN, 'StepOffset_ms', NaN, ...
    'MissCentred_pct', NaN, 'ExtraCentred_pct', NaN, 'PerfectTrials_pct', NaN, 'Pass', NaN, ...
    'Experiment', '', 'Variant', '', 'Extractor', '', 'Overlap', NaN, 'LagSeconds', NaN, 'LagSource', '', 'TcSource', '', ...
    'BandA', NaN, 'Passes', NaN, 'PthMode', '', 'Floor_dB', NaN, 'Margin_dB', NaN, 'Seed', NaN, ...
    'CyclesInWindow', NaN, 'AdaptK', NaN, ...
    'TcTrue_s', NaN, 'TcEst_s', NaN, 'TcErr_pct', NaN, 'TcAbsErr_pct', NaN, 'TcEstStd_s', NaN, 'TcFailed', NaN, ...
    'SelRMSE_pct', NaN, 'SelCPI_s', NaN, 'SelSNR_dB', NaN, 'SelFallback', NaN, 'SelNoStage2', NaN, 'GapPct', NaN);
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
            parts{k} = sprintf('%.6g', val);
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

function assembleResults(cfg, outDir)
R = struct();
expNames = {'E0', 'E1', 'E2', 'E3', 'E4a', 'E4v', 'E5', 'E6', 'E7'};
csvNames = {'E0_baseline.csv', 'E1_overlap.csv', 'E2_peak_rule_lag.csv', 'E3_Tc_estimate.csv', 'E4_window_sweep.csv', ...
    'E4_fixed_vs_perspeed.csv', 'E5_band_refinement.csv', 'E6_Pth_floor.csv', 'E7_selection.csv'};
ALL = [];
for q = 1:numel(expNames)
    rows = loadExp(outDir, expNames{q}, cfg.Speeds);
    key = strrep(expNames{q}, 'E4a', 'E4a');
    R.(key) = rows;
    if ~isempty(rows)
        writeRowsCSV(fullfile(outDir, csvNames{q}), rows);
        ALL = [ALL, rows];                                    %#ok<AGROW>
    end
end
if isempty(ALL)
    logMsg(outDir, 'No results to assemble.');
    return;
end
IDX = buildIndex(ALL);
C = buildComparison(cfg, ALL, IDX);
writeStructCSV(fullfile(outDir, 'Variant_Comparison.csv'), C);
PS = buildPassSummary(cfg, ALL, IDX);
writeStructCSV(fullfile(outDir, 'Pass_Summary.csv'), PS);
DD = designDecisions(cfg, ALL, IDX, C, PS);
writeDesignTable(fullfile(outDir, 'Design_Decisions_Tested.csv'), DD);
save(fullfile(outDir, 'All_Results.mat'), 'ALL', 'C', 'PS', 'DD');
if cfg.MakePlots
    try
        makePlots(cfg, outDir, ALL, IDX, PS);
    catch ME
        logMsg(outDir, 'Plots failed: %s', ME.message);
    end
end
writeReport(cfg, outDir, ALL, IDX, C, PS, DD);
end

function IDX = buildIndex(ALL)
IDX = containers.Map();
for i = 1:numel(ALL)
    IDX(rowKey(ALL(i).Experiment, ALL(i).Variant, ALL(i).Extractor, ALL(i).Parameter, ALL(i).SNR_dB, ALL(i).Speed_mps)) = i;
end
end

function k = rowKey(e, v, x, p, snr, spd)
k = sprintf('%s|%s|%s|%s|%g|%.2f', e, v, x, p, snr, spd);
end

function r = getRow(ALL, IDX, e, v, x, p, snr, spd)
k = rowKey(e, v, x, p, snr, spd);
if isKey(IDX, k)
    r = ALL(IDX(k));
else
    r = [];
end
end

function vec = getVec(ALL, IDX, e, v, x, p, snr, speeds, field)
vec = NaN(1, numel(speeds));
for s = 1:numel(speeds)
    r = getRow(ALL, IDX, e, v, x, p, snr, speeds(s));
    if ~isempty(r)
        vec(s) = r.(field);
    end
end
end

function x = extKey(param, ext)
if strcmp(param, 'WalkingSpeed')
    x = 'ws';
else
    x = ext;
end
end

function C = buildComparison(cfg, ALL, IDX)
C = [];
floorPP = 0.05;
for i = 1:numel(ALL)
    r = ALL(i);
    if any(strcmp(r.Experiment, {'E0', 'E7'}))
        continue;
    end
    xk = extKey(r.Parameter, r.Extractor);
    A = getRow(ALL, IDX, 'E0', 'BASE', xk, r.Parameter, r.SNR_dB, r.Speed_mps);
    B = getRow(ALL, IDX, 'E0', 'BASE_B', xk, r.Parameter, r.SNR_dB, r.Speed_mps);
    ref = A;
    refName = 'E0 BASE';
    if strcmp(r.Experiment, 'E4v')
        ref = getRow(ALL, IDX, 'E4v', 'baseline', xk, r.Parameter, r.SNR_dB, r.Speed_mps);   % paired: same fresh noise
        refName = 'E4v baseline';
    end
    if isempty(ref)
        continue;
    end
    noise = NaN;
    if ~isempty(A) && ~isempty(B)
        noise = abs(A.RMSE_pct - B.RMSE_pct);
    end
    c.Experiment = r.Experiment;
    c.Variant = r.Variant;
    c.Extractor = r.Extractor;
    c.Parameter = r.Parameter;
    c.Speed_mps = r.Speed_mps;
    c.SNR_dB = r.SNR_dB;
    c.Reference = refName;
    c.CPI_act_s = max([r.CPI_speed_act, r.CPI_step_act]);
    c.RMSE_pct = r.RMSE_pct;
    c.RMSE_ref = ref.RMSE_pct;
    c.dRMSE = r.RMSE_pct - ref.RMSE_pct;
    c.Bias_pct = r.MeanErr_pct;
    c.dBias = r.MeanErr_pct - ref.MeanErr_pct;
    c.TrialStd_pct = r.TrialStd_pct;
    c.dTrialStd = r.TrialStd_pct - ref.TrialStd_pct;
    c.Miss_pct = r.MissRate_pct;
    c.dMiss = r.MissRate_pct - ref.MissRate_pct;
    c.Extra_pct = r.ExtraRate_pct;
    c.dExtra = r.ExtraRate_pct - ref.ExtraRate_pct;
    c.Exact_pct = r.PerfectTrials_pct;
    c.dExact = r.PerfectTrials_pct - ref.PerfectTrials_pct;
    c.NoiseRMSE_pp = noise;
    c.NoiseUsed_pp = max(noise, floorPP);
    c.PassVariant = double(r.RMSE_pct < cfg.Target_pct);
    c.PassRef = double(ref.RMSE_pct < cfg.Target_pct);
    if ~isfinite(c.dRMSE)
        c.Direction = 'no value';
    elseif isnan(noise)
        c.Direction = 'no noise estimate';
    elseif abs(c.dRMSE) <= c.NoiseUsed_pp
        c.Direction = 'within noise';
    elseif c.dRMSE < 0
        c.Direction = 'better';
    else
        c.Direction = 'worse';
    end
    C = addCmp(C, c);
    clear c
end
end

function C = addCmp(C, c)
if isempty(C)
    C = c;
else
    C(end+1) = c;
end
end

function writeStructCSV(filename, S)
if isempty(S)
    return;
end
writeRowsCSV(filename, S);                                    % same writer, any flat struct array works
end

function PS = buildPassSummary(cfg, ALL, IDX)
PS = [];
keys = containers.Map();
for i = 1:numel(ALL)
    r = ALL(i);
    if any(strcmp(r.Experiment, {'E7'}))
        continue;
    end
    gk = sprintf('%s|%s|%s|%g', r.Experiment, r.Variant, extGroup(r), r.SNR_dB);
    if ~isKey(keys, gk)
        keys(gk) = struct('e', r.Experiment, 'v', r.Variant, 'x', extGroup(r), 'snr', r.SNR_dB);
    end
end
ks = keys.keys();
params = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideDuration_v3rule', 'StrideLength'};
for q = 1:numel(ks)
    g = keys(ks{q});
    ps.Experiment = g.e;
    ps.Variant = g.v;
    ps.Extractor = g.x;
    ps.SNR_dB = g.snr;
    cnt = NaN(1, numel(params));
    nSp = 0;
    for p = 1:numel(params)
        xk = extKey(params{p}, g.x);
        if strcmp(params{p}, 'WalkingSpeed') && ~strcmp(g.x, 'ws'), xk = 'ws'; end
        vals = NaN(1, numel(cfg.Speeds));
        for s = 1:numel(cfg.Speeds)
            r = getRow(ALL, IDX, g.e, g.v, xk, params{p}, g.snr, cfg.Speeds(s));
            if ~isempty(r), vals(s) = r.RMSE_pct; end
        end
        if any(isfinite(vals))
            cnt(p) = sum(vals < cfg.Target_pct);
            nSp = max(nSp, sum(isfinite(vals)));
        end
    end
    ps.nSpeeds = nSp;
    ps.WalkingSpeed_pass = cnt(1);
    ps.WS_fromBaseline = 0;
    if isnan(cnt(1)) && any(isfinite(cnt(2:5)))
        vals = getVec(ALL, IDX, 'E0', 'BASE', 'ws', 'WalkingSpeed', g.snr, cfg.Speeds, 'RMSE_pct');
        if any(isfinite(vals))
            ps.WalkingSpeed_pass = sum(vals < cfg.Target_pct);
            ps.WS_fromBaseline = 1;
        end
    end
    ps.StepTime_pass = cnt(2);
    ps.StrideDuration_pass = cnt(3);
    ps.StrideDuration_v3rule_pass = cnt(4);
    ps.StrideLength_pass = cnt(5);
    tot = [ps.WalkingSpeed_pass, ps.StepTime_pass, ps.StrideDuration_pass, ps.StrideLength_pass];
    if all(isfinite(tot))
        ps.Total_of_4 = sum(tot);
    else
        ps.Total_of_4 = NaN;
    end
    ps.OutOf = 4*nSp;
    PS = addCmp(PS, ps);
    clear ps
end
end

function x = extGroup(r)
if strcmp(r.Parameter, 'WalkingSpeed') && ~strcmp(r.Experiment, 'E4v')
    x = 'ws';
elseif strcmp(r.Experiment, 'E4v')
    x = 'main';
else
    x = r.Extractor;
end
end

function ps = getPS(PS, e, v, x, snr)
ps = [];
for i = 1:numel(PS)
    if strcmp(PS(i).Experiment, e) && strcmp(PS(i).Variant, v) && strcmp(PS(i).Extractor, x) && PS(i).SNR_dB == snr
        ps = PS(i);
        return;
    end
end
end

function vs = variantSummary(cfg, C, PS, e, v, ext)
vs = struct('nB', 0, 'nW', 0, 'nS', 0, 'passDiff', 0, 'passVar', 0, 'passBase', 0, 'worst', 'none', 'cells', 0, 'snrs', []);
if isempty(C)
    return;
end
m = strcmp({C.Experiment}, e) & strcmp({C.Variant}, v) & strcmp({C.Extractor}, ext) & ~strcmp({C.Parameter}, 'StrideDuration_v3rule');
if ~any(m)
    return;
end
Cm = C(m);
dir = {Cm.Direction};
vs.nB = sum(strcmp(dir, 'better'));
vs.nW = sum(strcmp(dir, 'worse'));
vs.nS = sum(strcmp(dir, 'within noise'));
vs.cells = numel(Cm);
vs.snrs = unique([Cm.SNR_dB]);
grp = ext;
isWS = strcmp(ext, 'ws');
for snr = vs.snrs
    pv = getPS(PS, e, v, grp, snr);
    if strcmp(e, 'E4v')
        pb = getPS(PS, 'E4v', 'baseline', grp, snr);
    else
        pb = getPS(PS, 'E0', 'BASE', grp, snr);
    end
    if isempty(pv) || isempty(pb)
        continue;
    end
    if isWS
        a = pv.WalkingSpeed_pass;
        b = pb.WalkingSpeed_pass;
    else
        a = pv.Total_of_4;
        b = pb.Total_of_4;
    end
    if isfinite(a) && isfinite(b)
        vs.passVar = vs.passVar + a;
        vs.passBase = vs.passBase + b;
    end
end
vs.passDiff = vs.passVar - vs.passBase;
w = find(strcmp(dir, 'worse'));
if ~isempty(w)
    [~, j] = max([Cm(w).dRMSE]);
    c = Cm(w(j));
    vs.worst = sprintf('%s at %.1f m/s, %g dB: %+.2f pp (RMSE %.2f %% against %.2f %%)', c.Parameter, c.Speed_mps, c.SNR_dB, c.dRMSE, c.RMSE_pct, c.RMSE_ref);
end
end

function [label, why] = verdictFor(vs)
if vs.cells == 0
    label = 'INCONCLUSIVE';
    why = 'no comparison cells (the variant did not run)';
elseif (vs.passDiff >= 1 && vs.nW <= vs.nB/3) || (vs.passDiff >= 0 && vs.nW == 0 && vs.nB >= max(3, 0.25*vs.cells))
    label = 'IMPLEMENT';
    why = sprintf('%+d passing conditions, %d cells better and %d worse beyond the baseline noise', vs.passDiff, vs.nB, vs.nW);
elseif vs.passDiff <= -1 || (vs.nW > vs.nB && vs.nW >= 2)
    label = 'DO NOT IMPLEMENT';
    why = sprintf('%+d passing conditions, %d cells better and %d worse beyond the baseline noise', vs.passDiff, vs.nB, vs.nW);
else
    label = 'INCONCLUSIVE';
    why = sprintf('%+d passing conditions, %d cells better, %d worse, %d within noise', vs.passDiff, vs.nB, vs.nW, vs.nS);
end
end

function s = softLine(cfg)
s = softwareName();
if exist('OCTAVE_VERSION', 'builtin')
    s = [s ' - these numbers must be repeated in MATLAB before they are called MATLAB results'];
end
end

function noise = noiseStat(cfg, ALL, IDX, param, snr, ext, field)
a = getVec(ALL, IDX, 'E0', 'BASE', extKey(param, ext), param, snr, cfg.Speeds, field);
b = getVec(ALL, IDX, 'E0', 'BASE_B', extKey(param, ext), param, snr, cfg.Speeds, field);
noise = meanFinite(abs(a - b));
end

function t = fmtv(x, f)
if nargin < 2, f = '%.2f'; end
if isnan(x), t = 'n/a'; else, t = sprintf(f, x); end
end

function L = tableLine(L, varargin)
L{end+1} = sprintf(varargin{:});
end

function DD = designDecisions(cfg, ALL, IDX, C, PS)
DD = {};
evid = evidenceText(cfg);
[txt1, lab1] = bestOf(cfg, C, PS, 'E1', {'ov50_c1', 'ov50_c0.5', 'ov87.5_c1', 'ov87.5_c0.5'}, 'main');
DD(end+1, :) = {'Spectrogram overlap of the step detector (and of the walking speed algorithm)', ...
    'Wang et al. [TBC: confirm the value in the paper, the code comment says 75 %, the student recalls 50 %]', ...
    'Overlap 50, 75 and 87.5 % crossed with c scale 1 and 0.5, on the main (v3) extractor and on the v5 extractor; walking speed overlap on its own', ...
    txt1, evid, proposal(lab1, 'keep 75 %')};
[txt1b, lab1b] = bestOf(cfg, C, PS, 'E1', {'ov75_c1'}, 'main');
DD(end+1, :) = {'Peak drop c of the step rule', 'Wang: c = (1/CPI)/f_h (one Doppler bin)', ...
    'c multiplied by 1 and by 0.5 (cPeakScale) at every overlap', ...
    ['c x 1 at 75 % overlap: ' txt1b], evid, proposal(lab1b, 'keep cPeakScale = 0.5 as in the v3 settings')};
names2 = {};
for h = cfg.E2Hops
    names2{end+1} = sprintf('hop%g_oracle', h);        %#ok<AGROW>
    names2{end+1} = sprintf('hop%g_est', h);           %#ok<AGROW>
    for q = 1:numel(cfg.E2FixedLags)
        names2{end+1} = sprintf('hop%g_fix%g', h, cfg.E2FixedLags(q)); %#ok<AGROW>
    end
end
[txt2, lab2] = bestOf(cfg, C, PS, 'E2', names2, 'main');
DD(end+1, :) = {'Neighbouring-frame peak rule', 'Wang: F(t) compared with F(t-1) and F(t+1)', ...
    'Compare with the frames at +/- a lag in seconds (0.25 x true Tc, 0.25 x estimated Tc, fixed 0.10 s and 0.15 s) at hops 0.25, 0.125, 0.0625', ...
    txt2, evid, proposal(lab2, 'keep the one-frame rule')};
[txt3, lab3] = bestOf(cfg, C, PS, 'E3', {'Tc_acorr', 'Tc_acorr2', 'Tc_twopass', 'Tc_nominal'}, 'main');
DD(end+1, :) = {'Tc used inside the extractor (Tc/8 median and smoothing windows)', 'Wang: Tc from the torso bob. This project: the true Tc from the simulator', ...
    'Estimated Tc: autocorrelation of the leg-frequency track, the same doubled from the step period, two-pass; and a fixed nominal 1.2 s', ...
    txt3, evid, proposal(lab3, 'an estimated Tc is needed on measured data, choose after reading the table in the report')};
names4 = {};
for c = cfg.E4CPI
    names4{end+1} = sprintf('cpi%.3f', c);             %#ok<AGROW>
end
[txt4, lab4] = bestOf(cfg, C, PS, 'E4a', names4, 'ws');
DD(end+1, :) = {'Walking speed window length', 'Power-of-two windows only (0.046, 0.093, 0.186, 0.372, 0.743 s) in the v3 sweeps', ...
    'Windows from 0.139 s to 1.5 s that are not powers of two (zero padded FFT)', txt4, evid, proposal(lab4, 'keep 0.372 s')};
names4v = {'fixed', 'perspeed'};
for K = cfg.E4K
    names4v{end+1} = sprintf('adaptK%g', K);           %#ok<AGROW>
end
[txt4v, lab4v] = bestOf(cfg, C, PS, 'E4v', names4v, 'main');
DD(end+1, :) = {'Per-speed settings versus one fixed setting versus adaptive window', 'Settings chosen per speed in v3', ...
    'One fixed speed window, per-speed windows, CPI = K/f with K = 25, 35, 50; validated on fresh noise', txt4v, evid, proposal(lab4v, 'keep per-speed settings')};
names5 = {};
for a = [0, cfg.E5Bands]
    for np = cfg.E5Passes
        if a == 0 && np == 3, continue; end
        names5{end+1} = sprintf('a%g_p%d', a, np);     %#ok<AGROW>
    end
end
[txt5, lab5] = bestOf(cfg, C, PS, 'E5', names5, 'ws');
DD(end+1, :) = {'Torso band refinement', 'Yardibi: band 0.5x to 2x the mean frequency, repeated passes (3 here)', ...
    'Band (1-a) to (1+a) times the mean tracked frequency with a = 0.15, 0.25, 0.4; 2, 3 and 5 passes', txt5, evid, proposal(lab5, 'keep 0.5x to 2x with 3 passes')};
names6 = {};
for f = cfg.E6Floors
    names6{end+1} = sprintf('floor%g', f);             %#ok<AGROW>
end
for m = cfg.E6Margins
    names6{end+1} = sprintf('noisefloor+%g', m);       %#ok<AGROW>
end
[txt6, lab6] = bestOf(cfg, C, PS, 'E6', names6, 'main');
DD(end+1, :) = {'Pth rule', 'Wang: a fixed threshold in dB. This project: Pth = max(10 - Gamma, -40) dB (supervisor rule with a floor)', ...
    'Floors -30, -35, -40, -45 dB; Pth from an estimated noise floor plus 6, 10 or 14 dB (no Gamma needed)', txt6, evid, proposal(lab6, 'keep the floor at -40 dB')};
DD(end+1, :) = {'Rule for choosing the window and SNR', 'Lowest RMSE, ties within 0.25 % go to the shortest window and lowest SNR', ...
    'Tie tolerance 0.5 and 1.0 %, and a plateau rule that needs neighbouring windows and SNRs to pass too', e7Text(cfg, ALL), evid, ...
    '[decide after reading the E7 section of the report]'};
DD(end+1, :) = {'Leg-frequency limit', 'Wang: 250 Hz at 5.8 GHz', 'Not varied. Fixed at 350 Hz scaled to 24 GHz in every copy (the last run used 300 Hz by accident)', ...
    'Not tested in this run', 'n/a', 'keep 350 Hz'};
end

function s = evidenceText(cfg)
if exist('OCTAVE_VERSION', 'builtin')
    s = 'Octave stand-in run, NOT MATLAB. Repeat in MATLAB before use';
else
    s = sprintf('MATLAB run, folder %s, %d trials, seeds %d / %d / %d', cfg.OutDir, cfg.nTrials, cfg.SeedA, cfg.SeedB, cfg.SeedVal);
end
end

function s = proposal(label, keepText)
switch label
    case 'IMPLEMENT'
        s = 'Proposed: adopt the best variant listed (needs review)';
    case 'DO NOT IMPLEMENT'
        s = ['Proposed: ' keepText];
    otherwise
        s = ['Inconclusive. Proposed: ' keepText];
end
end

function [txt, label] = bestOf(cfg, C, PS, e, names, ext)
best = '';
bestScore = -inf;
bestVs = [];
for q = 1:numel(names)
    vs = variantSummary(cfg, C, PS, e, names{q}, ext);
    if vs.cells == 0, continue; end
    sc = 1000*vs.passDiff + (vs.nB - vs.nW);
    if sc > bestScore
        bestScore = sc;
        best = names{q};
        bestVs = vs;
    end
end
if isempty(best)
    txt = 'Not run';
    label = 'INCONCLUSIVE';
    return;
end
[label, why] = verdictFor(bestVs);
txt = sprintf('Best variant %s: %s (%s)', best, why, label);
end

function s = e7Text(cfg, ALL)
s = 'Not run';
m = strcmp({ALL.Experiment}, 'E7');
if ~any(m)
    return;
end
R7 = ALL(m);
rules = unique({R7.Variant});
parts = {};
for q = 1:numel(rules)
    g = [R7(strcmp({R7.Variant}, rules{q})).GapPct];
    parts{end+1} = sprintf('%s mean gap %s pp', rules{q}, fmtv(meanFinite(g))); %#ok<AGROW>
end
s = ['Validation minus selection RMSE: ' strjoin(parts, ', ')];
end

function writeDesignTable(filename, DD)
fid = fopen(filename, 'w');
fprintf(fid, 'Decision,Published value,What was tried,What happened,Evidence,What was chosen\n');
for r = 1:size(DD, 1)
    parts = cell(1, 6);
    for c = 1:6
        parts{c} = ['"' strrep(DD{r, c}, '"', '''') '"'];
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
fclose(fid);
end

function writeReport(cfg, outDir, ALL, IDX, C, PS, DD)
L = {};
L{end+1} = 'VARIANT REPORT (written automatically)';
L{end+1} = sprintf('Run folder: %s', outDir);
L{end+1} = sprintf('Software: %s', softLine(cfg));
L{end+1} = sprintf('Speeds: %s m/s. Trials per condition: %d. Seeds: A = %d, B = %d (baseline repeat), V = %d (fresh validation).', mat2str(cfg.Speeds), cfg.nTrials, cfg.SeedA, cfg.SeedB, cfg.SeedVal);
L{end+1} = sprintf('Step detection: requested window %.3f s (actual %.4f s), SNR (Gamma) %s dB. Walking speed experiments: Gamma %s dB. Extractor limit: 350*24/5.8 = %.2f Hz in every copy.', ...
    cfg.StepCPIreq, 8192/44100, mat2str(cfg.StepSNRs), mat2str(cfg.SpeedSNRs), 350*24/5.8);
if cfg.nTrials < 30
    L{end+1} = sprintf('NOTE: %d trials were used instead of the planned 30, so the differences between variants are less certain.', cfg.nTrials);
end
L{end+1} = 'Everything below is SIMULATED data (Chen model A, R0 = 10 m). Nothing here says anything about measured IPM-365 data.';
L{end+1} = 'A statement marked [measured here] was computed in this run. A statement marked [hypothesis] was NOT measured. [TBC] = not computed.';
L{end+1} = '';
L{end+1} = 'HOW THE VERDICTS WERE DECIDED';
L{end+1} = '  Each variant is compared with the baseline cell by cell (speed x SNR x parameter) using the SAME noise draws (paired).';
L{end+1} = '  The baseline was run twice with different seeds. The difference between those two runs is the noise level of that cell.';
L{end+1} = '  A change smaller than the noise level (at least 0.05 percentage points) is called "within noise". Only larger changes count as better or worse.';
L{end+1} = '  The noise level comes from ONE pair of baseline runs, so it is a rough guide and not a confidence interval.';
L{end+1} = '  IMPLEMENT = the number of passing conditions (RMSE below 5 %) rises AND the cells that got worse are at most a third of the cells that got better,';
L{end+1} = '              or the passing conditions stay the same, nothing got worse, and at least a quarter of the cells (and at least 3) got better.';
L{end+1} = '  DO NOT IMPLEMENT = passes fall, or more cells got worse than better (at least 2 worse). Anything else = INCONCLUSIVE.';
L{end+1} = '  Stride duration is scored on the steps inside the walk (as in Master_SimSweep_v4). The v3 rule (every detection) is kept as StrideDuration_v3rule in the CSV files.';
L{end+1} = '';
L{end+1} = 'RUN-TO-RUN NOISE OF THE BASELINE (mean absolute difference of RMSE between seed A and seed B, over the speeds, percentage points)';
for snr = cfg.AllSNRs
    s = sprintf('  Gamma %g dB:', snr);
    s = [s ' walking speed ' fmtv(noiseStat(cfg, ALL, IDX, 'WalkingSpeed', snr, 'ws', 'RMSE_pct'))];
    if true
        s = [s ', step time ' fmtv(noiseStat(cfg, ALL, IDX, 'StepTime', snr, 'main', 'RMSE_pct'))];
        s = [s ', stride duration ' fmtv(noiseStat(cfg, ALL, IDX, 'StrideDuration', snr, 'main', 'RMSE_pct'))];
        s = [s ', stride length ' fmtv(noiseStat(cfg, ALL, IDX, 'StrideLength', snr, 'main', 'RMSE_pct'))];
    end
    L{end+1} = s;                                            %#ok<AGROW>
end
L{end+1} = '';
L = secE1(L, cfg, ALL, IDX, C, PS);
L = secE2(L, cfg, ALL, IDX, C, PS);
L = secE3(L, cfg, ALL, IDX, C, PS);
L = secE4(L, cfg, ALL, IDX, C, PS);
L = secE5(L, cfg, ALL, IDX, C, PS);
L = secE6(L, cfg, ALL, IDX, C, PS);
L = secE7(L, cfg, ALL);
L = secRanking(L, cfg, C, PS);
L{end+1} = 'WHAT WAS NOT TESTED';
L{end+1} = '  Measured (IPM-365) data. Other body heights, amplitude modes, ranges or aspect angles. Other noise models (the noise is complex Gaussian with the simulator equation).';
L{end+1} = '  Window lengths of the STEP detector other than 0.186 s in E1, E2, E3 and E6 (E4 sweeps five power-of-two step windows). Gamma values other than those listed above.';
L{end+1} = '  Walking speed overlap other than 50, 75 and 87.5 %. Any combination of two changes (each experiment changes one thing at a time).';
L{end+1} = '  The v5 changes (put back a skipped step, refined step times) were varied only in E1 (v5 baseline) and are otherwise not part of these experiments.';
L{end+1} = '  Whether Wang really used 50 % or 75 % overlap: [TBC] this must be read from the paper. These experiments only show what each choice does in this simulator.';
L{end+1} = '';
L{end+1} = 'DESIGN DECISIONS TESTED (also in Design_Decisions_Tested.csv)';
for r = 1:size(DD, 1)
    L{end+1} = sprintf('  Decision        : %s', DD{r, 1});          %#ok<AGROW>
    L{end+1} = sprintf('  Published value : %s', DD{r, 2});          %#ok<AGROW>
    L{end+1} = sprintf('  What was tried  : %s', DD{r, 3});          %#ok<AGROW>
    L{end+1} = sprintf('  What happened   : %s', DD{r, 4});          %#ok<AGROW>
    L{end+1} = sprintf('  Evidence        : %s', DD{r, 5});          %#ok<AGROW>
    L{end+1} = sprintf('  What was chosen : %s', DD{r, 6});          %#ok<AGROW>
    L{end+1} = '';                                                   %#ok<AGROW>
end
fid = fopen(fullfile(outDir, 'Variant_Report.txt'), 'w');
for k = 1:numel(L)
    fprintf(fid, '%s\n', L{k});
end
fclose(fid);
end

function L = variantLines(L, cfg, C, PS, e, names, ext)
for q = 1:numel(names)
    vs = variantSummary(cfg, C, PS, e, names{q}, ext);
    [label, why] = verdictFor(vs);
    L{end+1} = sprintf('  %-22s %s', names{q}, label);                                           %#ok<AGROW>
    L{end+1} = sprintf('      passing conditions: %d against %d for the baseline (summed over the SNRs run); %s [measured here]', vs.passVar, vs.passBase, why); %#ok<AGROW>
    L{end+1} = sprintf('      worst deterioration beyond noise: %s', vs.worst);                  %#ok<AGROW>
end
end

function [mn, nz] = fastMean(ALL, IDX, cfg, e, v, x, p, snr, field, minSpeed)
sp = cfg.Speeds(cfg.Speeds >= minSpeed);
mn = meanFinite(getVec(ALL, IDX, e, v, x, p, snr, sp, field));
a = getVec(ALL, IDX, 'E0', 'BASE', extKey(p, x), p, snr, sp, field);
b = getVec(ALL, IDX, 'E0', 'BASE_B', extKey(p, x), p, snr, sp, field);
nz = meanFinite(abs(a - b));
end

function L = secE1(L, cfg, ALL, IDX, C, PS)
m = strcmp({ALL.Experiment}, 'E1');
L{end+1} = 'E1  OVERLAP AND HOP (Wang 50 % against 75 %, and 87.5 %), crossed with the peak drop c';
L{end+1} = '  Hypothesis (stated before the run): with c x 1, a 50 % overlap reduces the fast-speed misses (>= 1.4 m/s) relative to 75 %, because frames are farther';
L{end+1} = '  apart so the drop between neighbours is larger, while the step-time error rises because the time-grid floor hop/sqrt(12) doubles (13.4 ms to 26.9 ms at 0.186 s).';
if ~any(m)
    L{end+1} = '  [not run]';
    L{end+1} = '';
    return;
end
exts = {'main'};
if any(strcmp({ALL(m).Extractor}, 'v5')), exts{end+1} = 'v5'; end
for xi = 1:numel(exts)
    ext = exts{xi};
    for snr = cfg.StepSNRs
        L{end+1} = sprintf('  Step detection, extractor %s, Gamma %g dB (mean over all speeds unless stated; "fast" = speeds >= 1.4 m/s):', ext, snr);
        L{end+1} = '    variant        ST_RMSE%  SD_RMSE%  fast_miss%  fast_extra%  exact%   pass ST/SD/SL   timing floor ms (calculated, hop/sqrt(12))';
        for ov = cfg.E1Overlaps
            for cs = cfg.E1CScales
                nm = sprintf('ov%g_c%g', 100*ov, cs);
                st = meanFinite(getVec(ALL, IDX, 'E1', nm, ext, 'StepTime', snr, cfg.Speeds, 'RMSE_pct'));
                sd = meanFinite(getVec(ALL, IDX, 'E1', nm, ext, 'StrideDuration', snr, cfg.Speeds, 'RMSE_pct'));
                fm = fastMean(ALL, IDX, cfg, 'E1', nm, ext, 'StepTime', snr, 'MissRate_pct', 1.4);
                fe = fastMean(ALL, IDX, cfg, 'E1', nm, ext, 'StepTime', snr, 'ExtraRate_pct', 1.4);
                ex = meanFinite(getVec(ALL, IDX, 'E1', nm, ext, 'StepTime', snr, cfg.Speeds, 'PerfectTrials_pct'));
                ps = getPS(PS, 'E1', nm, ext, snr);
                if isempty(ps), pstr = 'n/a'; else, pstr = sprintf('%d/%d/%d', ps.StepTime_pass, ps.StrideDuration_pass, ps.StrideLength_pass); end
                fl = 1000*(1 - ov)*(8192/44100)/sqrt(12);
                L{end+1} = sprintf('    %-14s %8s  %8s  %10s  %11s  %6s   %-13s   %.1f', nm, fmtv(st), fmtv(sd), fmtv(fm), fmtv(fe), fmtv(ex, '%.0f'), pstr, fl); %#ok<AGROW>
            end
        end
        if xi == 1
            [m50, n50] = fastMean(ALL, IDX, cfg, 'E1', 'ov50_c1', 'main', 'StepTime', snr, 'MissRate_pct', 1.4);
            m75 = fastMean(ALL, IDX, cfg, 'E1', 'ov75_c1', 'main', 'StepTime', snr, 'MissRate_pct', 1.4);
            r50 = meanFinite(getVec(ALL, IDX, 'E1', 'ov50_c1', 'main', 'StepTime', snr, cfg.Speeds, 'RMSE_pct'));
            r75 = meanFinite(getVec(ALL, IDX, 'E1', 'ov75_c1', 'main', 'StepTime', snr, cfg.Speeds, 'RMSE_pct'));
            nR = noiseStat(cfg, ALL, IDX, 'StepTime', snr, 'main', 'RMSE_pct');
            missNoise = max(n50, 0.5);
            rmseNoise = max(nR, 0.05);
            missDown = (m50 - m75) < -missNoise;
            rmseUp = (r50 - r75) > rmseNoise;
            if missDown && rmseUp
                lab = 'CONFIRMED';
            elseif ~missDown && ~rmseUp
                lab = 'REFUTED';
            else
                lab = 'PARTLY CONFIRMED';
            end
            L{end+1} = sprintf('    Check at Gamma %g dB (main extractor, c x 1): fast-speed miss 50 %% = %s, 75 %% = %s (noise %s pp): %s. Step-time RMSE 50 %% = %s, 75 %% = %s (noise %s pp): %s.', ...
                snr, fmtv(m50), fmtv(m75), fmtv(missNoise), yesNoText(missDown, 'reduced', 'not reduced'), fmtv(r50), fmtv(r75), fmtv(rmseNoise), yesNoText(rmseUp, 'rose', 'did not rise'));
            L{end+1} = sprintf('    Hypothesis at Gamma %g dB: %s [measured here]', snr, lab);
        end
    end
    a = getVec(ALL, IDX, 'E1', 'ov75_c0.5', ext, 'StepTime', cfg.StepSNRs(1), cfg.Speeds, 'RMSE_pct');
    b = getVec(ALL, IDX, 'E0', 'BASE', ext, 'StepTime', cfg.StepSNRs(1), cfg.Speeds, 'RMSE_pct');
    L{end+1} = sprintf('  Sanity check: ov75_c0.5 reproduces the baseline exactly (same seed): %s', yesNoText(isequaln(a, b), 'yes', 'NO, check the run'));
end
L{end+1} = '  Per-variant verdicts (against the baseline of the same extractor):';
nm = {};
for ov = cfg.E1Overlaps
    for cs = cfg.E1CScales
        if ov == 0.75 && cs == 0.5, continue; end
        nm{end+1} = sprintf('ov%g_c%g', 100*ov, cs);           %#ok<AGROW>
    end
end
for xi = 1:numel(exts)
    L{end+1} = sprintf('  Extractor %s:', exts{xi});
    L = variantLines(L, cfg, C, PS, 'E1', nm, exts{xi});
end
L{end+1} = '  Walking speed algorithm overlap (Gamma over the walking speed SNRs; passes out of the speeds run, mean RMSE and mean bias over speeds):';
for snr = cfg.SpeedSNRs
    for ov = cfg.E1Overlaps
        v = sprintf('ws_ov%g', 100*ov);
        ps = getPS(PS, 'E1', v, 'ws', snr);
        rm = meanFinite(getVec(ALL, IDX, 'E1', v, 'ws', 'WalkingSpeed', snr, cfg.Speeds, 'RMSE_pct'));
        bi = meanFinite(getVec(ALL, IDX, 'E1', v, 'ws', 'WalkingSpeed', snr, cfg.Speeds, 'MeanErr_pct'));
        if isempty(ps), pc = NaN; else, pc = ps.WalkingSpeed_pass; end
        L{end+1} = sprintf('    Gamma %g dB, overlap %g %%: pass %s, mean RMSE %s %%, mean bias %s %%', snr, 100*ov, fmtv(pc, '%d'), fmtv(rm), fmtv(bi)); %#ok<AGROW>
    end
end
L = variantLines(L, cfg, C, PS, 'E1', {'ws_ov50', 'ws_ov87.5'}, 'ws');
L{end+1} = '  Reading: the grid floor in the last column is a calculation. The measured effect of the overlap on the timing is the StepTime RMSE column.';
L{end+1} = '  Note: the v5 extractor places step times between frames (centroid), so the grid floor should matter much less there. Compare the two extractor blocks.';
L{end+1} = '';
end

function s = yesNoText(tf, a, b)
if tf, s = a; else, s = b; end
end

function L = secE2(L, cfg, ALL, IDX, C, PS)
m = strcmp({ALL.Experiment}, 'E2');
L{end+1} = 'E2  HOP-INDEPENDENT PEAK RULE (compare with the frames at +/- a lag in seconds)';
L{end+1} = '  Hypothesis (stated before the run): the stride-duration failures at a finer hop appear because the drop between neighbouring frames shrinks with the hop, so real peaks fail the';
L{end+1} = '  half-bin test. If so, a drop defined over a fixed TIME lag removes the hop dependence and the finer hop lowers the timing error without breaking stride duration.';
L{end+1} = '  Extra control: "frame_cscaled" shrinks c with the hop exactly as Master_SimSweep_v3 Stage 3 did. "frame" keeps c fixed (the published formula).';
if ~any(m)
    L{end+1} = '  [not run]';
    L{end+1} = '';
    return;
end
rules = {'frame', 'frame_cscaled', 'oracle', 'est'};
for q = 1:numel(cfg.E2FixedLags)
    rules{end+1} = sprintf('fix%g', cfg.E2FixedLags(q));      %#ok<AGROW>
end
for snr = cfg.StepSNRs
    L{end+1} = sprintf('  Gamma %g dB, step window 0.186 s, c scale 0.5 (mean over speeds unless stated):', snr);
    L{end+1} = '    hop     rule            pass ST/SD/SL   ST_RMSE%  SD_RMSE%  worst SD% (speed)   miss+extra%  exact%';
    for hop = cfg.E2Hops
        for r = 1:numel(rules)
            nm = sprintf('hop%g_%s', hop, rules{r});
            ps = getPS(PS, 'E2', nm, 'main', snr);
            if isempty(ps), continue; end
            st = meanFinite(getVec(ALL, IDX, 'E2', nm, 'main', 'StepTime', snr, cfg.Speeds, 'RMSE_pct'));
            sdv = getVec(ALL, IDX, 'E2', nm, 'main', 'StrideDuration', snr, cfg.Speeds, 'RMSE_pct');
            [wv, wi] = max(sdv);
            me = meanFinite(getVec(ALL, IDX, 'E2', nm, 'main', 'StepTime', snr, cfg.Speeds, 'MissRate_pct') + getVec(ALL, IDX, 'E2', nm, 'main', 'StepTime', snr, cfg.Speeds, 'ExtraRate_pct'));
            ex = meanFinite(getVec(ALL, IDX, 'E2', nm, 'main', 'StepTime', snr, cfg.Speeds, 'PerfectTrials_pct'));
            if isempty(wi) || isnan(wv), ws = 'n/a'; else, ws = sprintf('%.2f (%.1f)', wv, cfg.Speeds(wi)); end
            L{end+1} = sprintf('    %-7g %-15s %d/%d/%d           %8s  %8s  %-18s  %10s  %6s', hop, rules{r}, ps.StepTime_pass, ps.StrideDuration_pass, ps.StrideLength_pass, ...
                fmtv(st), fmtv(mean(sdv(isfinite(sdv)))), ws, fmtv(me), fmtv(ex, '%.0f')); %#ok<AGROW>
        end
    end
    finest = min(cfg.E2Hops);
    baseSD = getPS(PS, 'E2', sprintf('hop%g_frame', 0.25), 'main', snr);
    frameFine = getPS(PS, 'E2', sprintf('hop%g_frame', finest), 'main', snr);
    scaledFine = getPS(PS, 'E2', sprintf('hop%g_frame_cscaled', finest), 'main', snr);
    if ~isempty(baseSD) && ~isempty(frameFine)
        L{end+1} = sprintf('    Stride-duration passes (of %d speeds): hop 0.25 frame rule %d; finest hop frame rule (c fixed) %d; finest hop frame rule (c shrinks with hop, as in v3) %s [measured here].', ...
            baseSD.nSpeeds, baseSD.StrideDuration_pass, frameFine.StrideDuration_pass, ifEmpty(scaledFine, 'StrideDuration_pass'));
    end
    stBase = meanFinite(getVec(ALL, IDX, 'E0', 'BASE', 'main', 'StepTime', snr, cfg.Speeds, 'RMSE_pct'));
    nR = max(noiseStat(cfg, ALL, IDX, 'StepTime', snr, 'main', 'RMSE_pct'), 0.05);
    okRules = {};
    for r = 3:numel(rules)
        nm = sprintf('hop%g_%s', finest, rules{r});
        ps = getPS(PS, 'E2', nm, 'main', snr);
        if isempty(ps) || isempty(baseSD), continue; end
        st = meanFinite(getVec(ALL, IDX, 'E2', nm, 'main', 'StepTime', snr, cfg.Speeds, 'RMSE_pct'));
        if ps.StrideDuration_pass >= baseSD.StrideDuration_pass && (stBase - st) > nR
            okRules{end+1} = rules{r};                         %#ok<AGROW>
        end
    end
    if isempty(okRules)
        lab = 'REFUTED or INCONCLUSIVE: no lag rule keeps stride duration and lowers the step-time error beyond noise at the finest hop';
    else
        lab = ['CONFIRMED for: ' strjoin(okRules, ', ')];
    end
    L{end+1} = sprintf('    Hypothesis at Gamma %g dB: %s [measured here]', snr, lab);
end
names = {};
for hop = cfg.E2Hops
    for r = 1:numel(rules)
        if hop == 0.25 && strcmp(rules{r}, 'frame'), continue; end
        names{end+1} = sprintf('hop%g_%s', hop, rules{r});    %#ok<AGROW>
    end
end
L{end+1} = '  Per-variant verdicts (against the baseline = hop 0.25, one-frame rule):';
L = variantLines(L, cfg, C, PS, 'E2', names, 'main');
L{end+1} = '';
end

function s = ifEmpty(ps, field)
if isempty(ps), s = 'n/a'; else, s = sprintf('%d', ps.(field)); end
end

function L = secE3(L, cfg, ALL, IDX, C, PS)
m = strcmp({ALL.Experiment}, 'E3');
L{end+1} = 'E3  ESTIMATED Tc INSIDE THE EXTRACTOR (the scorer keeps the true Tc)';
L{end+1} = '  Hypothesis (stated before the run): replacing the true Tc by an estimate computed from the signal changes the step metrics by less than the baseline noise,';
L{end+1} = '  because Tc only sets the Tc/8 median and smoothing windows. Estimators: acorr = strongest autocorrelation peak of the leg-frequency track in 0.6 to 3 s;';
L{end+1} = '  acorr2 = twice the strongest peak in 0.3 to 1.5 s (peak taken as the step period); twopass = twice the median step interval of a first pass run with Tc = 1.2 s.';
if ~any(m)
    L{end+1} = '  [not run]';
    L{end+1} = '';
    return;
end
ests = {'acorr', 'acorr2', 'twopass'};
snr = cfg.StepSNRs(end);
L{end+1} = sprintf('  Estimator error (Gamma %g dB; estimate mean +/- std over trials, error in %% of the true Tc, and trials where the estimate failed and 1.2 s was used):', snr);
L{end+1} = '    speed   true Tc   acorr (err%/fail)           acorr2 (err%/fail)          twopass (err%/fail)';
for s = 1:numel(cfg.Speeds)
    line = sprintf('    %-6.1f  %-8s', cfg.Speeds(s), '');
    tcTrue = NaN;
    for q = 1:numel(ests)
        r = getRow(ALL, IDX, 'E3', ['Tc_' ests{q}], 'main', 'StepTime', snr, cfg.Speeds(s));
        if isempty(r)
            line = [line sprintf('%-28s', 'n/a')];            %#ok<AGROW>
        else
            tcTrue = r.TcTrue_s;
            line = [line sprintf('%-28s', sprintf('%.2f+/-%.2f (%+.1f/%d)', r.TcEst_s, r.TcEstStd_s, r.TcErr_pct, r.TcFailed))]; %#ok<AGROW>
        end
    end
    line(11:18) = sprintf('%-8s', fmtv(tcTrue, '%.3f'));
    L{end+1} = line;                                           %#ok<AGROW>
end
for q = 1:numel(ests)
    errs = getVec(ALL, IDX, 'E3', ['Tc_' ests{q}], 'main', 'StepTime', snr, cfg.Speeds, 'TcAbsErr_pct');
    L{end+1} = sprintf('    %-8s mean absolute error over speeds %s %% [measured here]', ests{q}, fmtv(meanFinite(errs)));
end
a = getVec(ALL, IDX, 'E3', 'Tc_true', 'main', 'StepTime', snr, cfg.Speeds, 'RMSE_pct');
b = getVec(ALL, IDX, 'E0', 'BASE', 'main', 'StepTime', snr, cfg.Speeds, 'RMSE_pct');
L{end+1} = sprintf('  Sanity check: Tc_true reproduces the baseline exactly (same seed): %s', yesNoText(isequaln(a, b), 'yes', 'NO, check the run'));
L{end+1} = '  Effect on the step metrics against the true-Tc baseline (nominal = a fixed 1.2 s with no estimation at all):';
L = variantLines(L, cfg, C, PS, 'E3', {'Tc_nominal', 'Tc_acorr', 'Tc_acorr2', 'Tc_twopass'}, 'main');
for q = 1:numel(ests)
    vs = variantSummary(cfg, C, PS, 'E3', ['Tc_' ests{q}], 'main');
    if vs.cells > 0
        if vs.nW == 0 && vs.nB == 0
            res = 'CONFIRMED (no cell changed beyond the noise)';
        elseif vs.nW <= 2
            res = sprintf('MOSTLY CONFIRMED (%d cells worse, %d better beyond the noise)', vs.nW, vs.nB);
        else
            res = sprintf('REFUTED (%d cells worse beyond the noise)', vs.nW);
        end
        L{end+1} = sprintf('  Hypothesis for %s: %s [measured here]', ests{q}, res); %#ok<AGROW>
    end
end
L{end+1} = '';
end

function L = secE4(L, cfg, ALL, IDX, C, PS)
m = strcmp({ALL.Experiment}, 'E4a');
L{end+1} = 'E4  WALKING SPEED WINDOW LENGTH (not forced to a power of two) AND FIXED / PER-SPEED / ADAPTIVE CHOICES';
L{end+1} = '  Hypothesis (stated before the run): the walking-speed RMSE is large below about 25 cycles of the torso Doppler (T x f, f = 160 x v Hz) in the window and flat above it,';
L{end+1} = '  so a longer window than 0.372 s gains nothing, and one fixed window (or CPI = K/f) passes as many speeds as per-speed choices.';
if ~any(m)
    L{end+1} = '  [not run]';
    L{end+1} = '';
    return;
end
for snr = cfg.SpeedSNRs
    L{end+1} = sprintf('  Walking speed RMSE %% against the window (Gamma %g dB, 75 %% overlap, 3 passes, band 0.5x to 2x). "-" = window dropped (too few frames in the walk):', snr);
    hdr = '    window s ';
    for s = 1:numel(cfg.Speeds), hdr = [hdr sprintf('%7.1f', cfg.Speeds(s))]; end %#ok<AGROW>
    L{end+1} = [hdr '   mean(>=0.5)'];
    for c = cfg.E4CPI
        nm = sprintf('cpi%.3f', c);
        v = getVec(ALL, IDX, 'E4a', nm, 'ws', 'WalkingSpeed', snr, cfg.Speeds, 'RMSE_pct');
        line = sprintf('    %-8.3f ', c);
        for s = 1:numel(v)
            if isnan(v(s)), line = [line sprintf('%7s', '-')]; else, line = [line sprintf('%7.2f', v(s))]; end %#ok<AGROW>
        end
        L{end+1} = [line sprintf('   %s', fmtv(meanFinite(v(cfg.Speeds >= 0.5))))]; %#ok<AGROW>
    end
    L{end+1} = '  Mean bias % (same grid) at the slow speeds 0.4 to 0.6 m/s:';
    for c = cfg.E4CPI
        nm = sprintf('cpi%.3f', c);
        b = getVec(ALL, IDX, 'E4a', nm, 'ws', 'WalkingSpeed', snr, cfg.Speeds, 'MeanErr_pct');
        sd = getVec(ALL, IDX, 'E4a', nm, 'ws', 'WalkingSpeed', snr, cfg.Speeds, 'TrialStd_pct');
        sp = find(cfg.Speeds <= 0.6);
        line = sprintf('    %-8.3f ', c);
        for s = sp
            line = [line sprintf(' %.1f m/s: bias %s, spread %s;', cfg.Speeds(s), fmtv(b(s), '%+.1f'), fmtv(sd(s), '%.1f'))]; %#ok<AGROW>
        end
        L{end+1} = line;                                       %#ok<AGROW>
    end
    starts = NaN(1, numel(cfg.Speeds));
    for s = 1:numel(cfg.Speeds)
        r2 = NaN(1, numel(cfg.E4CPI));
        for q = 1:numel(cfg.E4CPI)
            rr = getRow(ALL, IDX, 'E4a', sprintf('cpi%.3f', cfg.E4CPI(q)), 'ws', 'WalkingSpeed', snr, cfg.Speeds(s));
            if ~isempty(rr), r2(q) = rr.RMSE_pct; end
        end
        for q = 1:numel(r2)
            tail = r2(q:end);
            tail = tail(isfinite(tail));
            if isfinite(r2(q)) && all(tail < 2.5)
                starts(s) = cfg.E4CPI(q);
                break;
            end
        end
    end
    L{end+1} = sprintf('  Plateau start = the shortest window from which every longer window stays under 2.5 %% RMSE (Gamma %g dB):', snr);
    line = '    ';
    cyc = NaN(1, numel(cfg.Speeds));
    for s = 1:numel(cfg.Speeds)
        cyc(s) = starts(s)*160*cfg.Speeds(s);
        line = [line sprintf('%.1f m/s: %s s = %s cycles;  ', cfg.Speeds(s), fmtv(starts(s), '%.3f'), fmtv(cyc(s), '%.0f'))]; %#ok<AGROW>
    end
    L{end+1} = line;
    fastIdx = cfg.Speeds >= 0.8;
    medC = median(cyc(fastIdx & isfinite(cyc)));
    if isempty(medC) || isnan(medC)
        lab = 'INCONCLUSIVE (no plateau under 2.5 % found)';
    elseif medC >= 15 && medC <= 35
        lab = sprintf('CONFIRMED (median plateau start %.0f cycles at speeds >= 0.8 m/s)', medC);
    else
        lab = sprintf('REFUTED (median plateau start %.0f cycles at speeds >= 0.8 m/s, not near 25)', medC);
    end
    L{end+1} = sprintf('  Hypothesis about the 25-cycle rule at Gamma %g dB: %s [measured here]', snr, lab);
    ref = meanFinite(getVec(ALL, IDX, 'E4a', 'cpi0.372', 'ws', 'WalkingSpeed', snr, cfg.Speeds(cfg.Speeds >= 0.5), 'RMSE_pct'));
    lng = NaN;
    for c = sort(cfg.E4CPI, 'descend')
        v = getVec(ALL, IDX, 'E4a', sprintf('cpi%.3f', c), 'ws', 'WalkingSpeed', snr, cfg.Speeds(cfg.Speeds >= 0.5), 'RMSE_pct');
        if all(isfinite(v))
            lng = meanFinite(v);
            lc = c;
            break;
        end
    end
    if isfinite(ref) && isfinite(lng)
        L{end+1} = sprintf('  Longest window complete at all speeds >= 0.5 m/s (%.3f s): mean RMSE %s %% against %s %% at 0.372 s. [measured here]', lc, fmtv(lng), fmtv(ref));
    end
end
mv = strcmp({ALL.Experiment}, 'E4v');
if any(mv)
    for snr = cfg.SpeedSNRs
        L{end+1} = sprintf('  Validation on fresh noise (seed V), Gamma %g dB: passes out of %d speeds for each way of choosing the windows:', snr, numel(cfg.Speeds));
        L{end+1} = '    configuration   speed  stepT  strideD  strideL  total (of 4 x speeds)';
        nms = {'baseline', 'fixed', 'perspeed'};
        for K = cfg.E4K, nms{end+1} = sprintf('adaptK%g', K); end %#ok<AGROW>
        for q = 1:numel(nms)
            ps = getPS(PS, 'E4v', nms{q}, 'main', snr);
            if isempty(ps), continue; end
            L{end+1} = sprintf('    %-14s  %5d  %5d  %7d  %7d  %8s', nms{q}, ps.WalkingSpeed_pass, ps.StepTime_pass, ps.StrideDuration_pass, ps.StrideLength_pass, fmtv(ps.Total_of_4, '%d')); %#ok<AGROW>
        end
    end
    L{end+1} = '  Verdicts against the baseline configuration of the same fresh noise:';
    nms = {'fixed', 'perspeed'};
    for K = cfg.E4K, nms{end+1} = sprintf('adaptK%g', K); end %#ok<AGROW>
    L = variantLines(L, cfg, C, PS, 'E4v', nms, 'main');
    L{end+1} = '  Note: the fixed and per-speed windows were chosen on seed A (E4a) and validated on seed V. The adaptive rule keeps the raw noise at the level of the 0.372 s coarse pass.';
end
L{end+1} = '  Per-window verdicts for the sweep (against the 0.372 s power-of-two baseline):';
nm = {};
for c = cfg.E4CPI, nm{end+1} = sprintf('cpi%.3f', c); end       %#ok<AGROW>
L = variantLines(L, cfg, C, PS, 'E4a', nm, 'ws');
L{end+1} = '';
end

function L = secE5(L, cfg, ALL, IDX, C, PS)
m = strcmp({ALL.Experiment}, 'E5');
L{end+1} = 'E5  TORSO BAND REFINEMENT (narrower bands, other numbers of passes), window 0.372 s';
L{end+1} = '  Hypothesis (stated before the run): a band of (1 -/+ a) x the mean tracked frequency removes the limb energy that pulls the torso peak at 0.4 to 0.8 m/s,';
L{end+1} = '  so the bias and the trial spread fall there, with no loss at the higher speeds.';
if ~any(m)
    L{end+1} = '  [not run]';
    L{end+1} = '';
    return;
end
slow = cfg.Speeds(cfg.Speeds <= 0.8);
for snr = cfg.SpeedSNRs
    L{end+1} = sprintf('  Gamma %g dB. Slow speeds = %s m/s; "all" = every speed run:', snr, mat2str(slow));
    L{end+1} = '    variant       pass(all)  RMSE slow%  |bias| slow%  spread slow%  RMSE all%';
    names = {'BASE'};
    for a = [0, cfg.E5Bands]
        for np = cfg.E5Passes
            if a == 0 && np == 3, continue; end
            names{end+1} = sprintf('a%g_p%d', a, np);          %#ok<AGROW>
        end
    end
    for q = 1:numel(names)
        if q == 1
            e = 'E0'; v = 'BASE';
        else
            e = 'E5'; v = names{q};
        end
        rm = getVec(ALL, IDX, e, v, 'ws', 'WalkingSpeed', snr, cfg.Speeds, 'RMSE_pct');
        bi = getVec(ALL, IDX, e, v, 'ws', 'WalkingSpeed', snr, cfg.Speeds, 'MeanErr_pct');
        sp = getVec(ALL, IDX, e, v, 'ws', 'WalkingSpeed', snr, cfg.Speeds, 'TrialStd_pct');
        sl = ismember(cfg.Speeds, slow);
        L{end+1} = sprintf('    %-12s  %8d  %10s  %12s  %12s  %9s', names{q}, sum(rm < cfg.Target_pct), fmtv(meanFinite(rm(sl))), fmtv(meanFinite(abs(bi(sl)))), fmtv(meanFinite(sp(sl))), fmtv(meanFinite(rm))); %#ok<AGROW>
    end
end
L{end+1} = '  Per-variant verdicts (against the baseline 3 passes, band 0.5x to 2x; the baseline window is the same 0.372 s):';
nm = {};
for a = [0, cfg.E5Bands]
    for np = cfg.E5Passes
        if a == 0 && np == 3, continue; end
        nm{end+1} = sprintf('a%g_p%d', a, np);                 %#ok<AGROW>
    end
end
L = variantLines(L, cfg, C, PS, 'E5', nm, 'ws');
L{end+1} = '';
end

function L = secE6(L, cfg, ALL, IDX, C, PS)
m = strcmp({ALL.Experiment}, 'E6');
L{end+1} = 'E6  Pth FLOOR AND Pth FROM AN ESTIMATED NOISE FLOOR';
L{end+1} = '  Hypothesis (stated before the run): the -40 dB floor only hurts at 0.4 m/s at high SNR, and a Pth set from the estimated noise floor (no Gamma needed) matches the coupled rule';
L{end+1} = '  at 30 dB and is at least as good at 60 dB.';
if ~any(m)
    L{end+1} = '  [not run]';
    L{end+1} = '';
    return;
end
names = {};
for f = cfg.E6Floors, names{end+1} = sprintf('floor%g', f); end   %#ok<AGROW>
for mg = cfg.E6Margins, names{end+1} = sprintf('noisefloor+%g', mg); end %#ok<AGROW>
for snr = cfg.StepSNRs
    L{end+1} = sprintf('  Gamma %g dB (pass counts out of %d speeds; 0.4 m/s columns show where the floor mattered):', snr, numel(cfg.Speeds));
    L{end+1} = '    variant           pass ST/SD/SL   mean Pth dB   ST_RMSE 0.4 m/s   miss+extra 0.4 m/s   SD_RMSE 0.4 m/s';
    for q = 0:numel(names)
        if q == 0
            e = 'E0'; v = 'BASE';
        else
            e = 'E6'; v = names{q};
        end
        ps = getPS(PS, e, v, 'main', snr);
        if isempty(ps), continue; end
        pth = meanFinite(getVec(ALL, IDX, e, v, 'main', 'StepTime', snr, cfg.Speeds, 'Pth_dB'));
        r4 = getRow(ALL, IDX, e, v, 'main', 'StepTime', snr, 0.4);
        d4 = getRow(ALL, IDX, e, v, 'main', 'StrideDuration', snr, 0.4);
        if isempty(r4), s1 = 'n/a'; s2 = 'n/a'; else, s1 = fmtv(r4.RMSE_pct); s2 = fmtv(r4.MissRate_pct + r4.ExtraRate_pct); end
        if isempty(d4), s3 = 'n/a'; else, s3 = fmtv(d4.RMSE_pct); end
        L{end+1} = sprintf('    %-16s  %d/%d/%d           %9s     %13s     %16s     %13s', v, ps.StepTime_pass, ps.StrideDuration_pass, ps.StrideLength_pass, fmtv(pth, '%.1f'), s1, s2, s3); %#ok<AGROW>
    end
end
L{end+1} = '  Per-variant verdicts (against the coupled rule with the -40 dB floor):';
L = variantLines(L, cfg, C, PS, 'E6', names, 'main');
L{end+1} = '';
end

function L = secE7(L, cfg, ALL)
m = strcmp({ALL.Experiment}, 'E7');
L{end+1} = 'E7  SELECTION ROBUSTNESS (re-select from the existing Stage 1 and Stage 2 CSVs, validate on 30 fresh draws)';
L{end+1} = sprintf('  Hypothesis (stated before the run): the lowest-RMSE rule picks lucky draws, so its validation RMSE is worse than its selection RMSE (a winner''s-curse gap), and a larger tie tolerance or a plateau rule shrinks the gap. Source run: %s', cfg.E7RunDir);
if ~any(m)
    L{end+1} = '  [not run: the Stage 1 and Stage 2 CSV files of the source run were not found. Set cfg.E7RunDir.]';
    L{end+1} = '';
    return;
end
R7 = ALL(m);
rules = {};
for tol = cfg.E7TieTols, rules{end+1} = sprintf('tie%g', tol); end %#ok<AGROW>
rules{end+1} = 'plateau';
L{end+1} = '  Gap = validation RMSE minus selection RMSE, in percentage points. A positive gap means the selection was optimistic.';
L{end+1} = '    rule      parameter        cells  mean sel RMSE%  mean val RMSE%  mean gap  cells passing sel -> val   fallback  S2 from other window';
for q = 1:numel(rules)
    for p = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideLength'}
        r = R7(strcmp({R7.Variant}, rules{q}) & strcmp({R7.Parameter}, p{1}));
        if isempty(r), continue; end
        L{end+1} = sprintf('    %-8s  %-15s  %5d  %14s  %14s  %8s  %9d -> %-3d          %8d  %d', rules{q}, p{1}, numel(r), fmtv(meanFinite([r.SelRMSE_pct])), ...
            fmtv(meanFinite([r.RMSE_pct])), fmtv(meanFinite([r.GapPct])), sum([r.SelRMSE_pct] < cfg.Target_pct), sum([r.RMSE_pct] < cfg.Target_pct), ...
            sumFinite([r.SelFallback]), sumFinite([r.SelNoStage2])); %#ok<AGROW>
    end
end
g0 = meanFinite([R7(strcmp({R7.Variant}, 'tie0.25')).GapPct]);
gP = meanFinite([R7(strcmp({R7.Variant}, 'plateau')).GapPct]);
if isfinite(g0) && isfinite(gP)
    if g0 > 0.25
        s1 = sprintf('CONFIRMED that the baseline rule is optimistic (mean gap %.2f pp)', g0);
    else
        s1 = sprintf('NOT CONFIRMED for the baseline rule (mean gap %.2f pp)', g0);
    end
    if gP < g0 - 0.25
        s2 = sprintf('plateau rule shrinks the gap (%.2f pp)', gP);
    else
        s2 = sprintf('plateau rule does not clearly shrink the gap (%.2f pp)', gP);
    end
    L{end+1} = sprintf('  Hypothesis: %s; %s [measured here]', s1, s2);
end
L{end+1} = '  Caveat: the Stage 1 and Stage 2 CSVs of the source run came from the extractor settings of that run (the v3 run used a 300 Hz limit by accident),';
L{end+1} = '  while the validation here uses 350 Hz. Part of the gap can therefore come from that change and not from the selection rule. [hypothesis, not separated]';
L{end+1} = '  Point cfg.E7RunDir at a run made with 350 Hz (for example a new Master_SimSweep_v4 folder) to remove this.';
L{end+1} = '';
end

function L = secRanking(L, cfg, C, PS)
L{end+1} = 'RANKED LIST OF CHANGES WORTH ADOPTING (verdict IMPLEMENT, largest gain first)';
if isempty(C)
    L{end+1} = '  [no comparisons]';
    L{end+1} = '';
    return;
end
keyList = containers.Map();
for i = 1:numel(C)
    k = sprintf('%s|%s|%s', C(i).Experiment, C(i).Variant, C(i).Extractor);
    keyList(k) = {C(i).Experiment, C(i).Variant, C(i).Extractor};
end
ks = keyList.keys();
items = {};
score = [];
for q = 1:numel(ks)
    t = keyList(ks{q});
    vs = variantSummary(cfg, C, PS, t{1}, t{2}, t{3});
    [label, ~] = verdictFor(vs);
    if strcmp(label, 'IMPLEMENT') && ~strcmp(t{1}, 'E4a')
        nS = max(1, numel(vs.snrs));
        items{end+1} = sprintf('  %s / %s (%s): about %+.1f passing conditions per SNR out of %d, %d cells better and %d worse beyond the noise', t{1}, t{2}, t{3}, vs.passDiff/nS, 4*numel(cfg.Speeds), vs.nB, vs.nW); %#ok<AGROW>
        score(end+1) = 1000*vs.passDiff/nS + vs.nB - vs.nW; %#ok<AGROW>
    end
end
if isempty(items)
    L{end+1} = '  None of the tested changes met the IMPLEMENT rule in this run.';
else
    [~, order] = sort(score, 'descend');
    for q = order
        L{end+1} = items{q};                                   %#ok<AGROW>
    end
    L{end+1} = '  The gains are [measured here] on simulated data and are not independent of each other. Do not add them up.';
end
L{end+1} = '';
end

function s = sumFinite(x)
x = x(isfinite(x));
s = sum(x);
end

function makePlots(cfg, outDir, ALL, IDX, PS)
plotDir = fullfile(outDir, 'Plots');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end
nT = cfg.nTrials;
if any(strcmp({ALL.Experiment}, 'E4a'))
    for snr = cfg.SpeedSNRs
        for fld = {'RMSE_pct', 'MeanErr_pct'}
            try
                plotE4(cfg, plotDir, ALL, IDX, snr, fld{1});
            catch ME
                logMsg(outDir, 'E4 plot failed: %s', ME.message);
            end
        end
    end
end
if any(strcmp({ALL.Experiment}, 'E2'))
    for snr = cfg.StepSNRs
        try
            plotE2(cfg, plotDir, ALL, IDX, snr);
        catch ME
            logMsg(outDir, 'E2 plot failed: %s', ME.message);
        end
    end
end
if any(strcmp({ALL.Experiment}, 'E3'))
    for snr = cfg.StepSNRs
        try
            plotE3(cfg, plotDir, ALL, IDX, snr);
        catch ME
            logMsg(outDir, 'E3 plot failed: %s', ME.message);
        end
    end
end
groups = {'E1', 'main', 'E1 step detection, main extractor'; 'E1', 'v5', 'E1 step detection, v5 extractor'; 'E2', 'main', 'E2 peak rule and hop'; ...
    'E3', 'main', 'E3 estimated Tc'; 'E6', 'main', 'E6 Pth floor'; 'E4v', 'main', 'E4 fixed, per-speed and adaptive windows'; 'E5', 'ws', 'E5 torso band refinement (walking speed only)'; ...
    'E1', 'ws', 'E1 walking speed overlap (walking speed only)'};
for g = 1:size(groups, 1)
    try
        plotPass(cfg, plotDir, PS, groups{g, 1}, groups{g, 2}, groups{g, 3}, nT);
    catch ME
        logMsg(outDir, 'Pass summary plot (%s) failed: %s', groups{g, 1}, ME.message);
    end
end
end

function fig = newFig()
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 560]);
end

function ax = styledAxes(fig, pos)
ax = axes('Parent', fig, 'Position', pos);
hold(ax, 'on');
grid(ax, 'on');
set(ax, 'FontSize', 10, 'LineWidth', 1, 'Box', 'on');
end

function plotE4(cfg, plotDir, ALL, IDX, snr, fld)
fig = newFig();
cols = lines(numel(cfg.Speeds));
isBias = strcmp(fld, 'MeanErr_pct');
for panel = 1:2
    if panel == 1
        ax = styledAxes(fig, [0.08 0.12 0.37 0.72]);
    else
        ax = styledAxes(fig, [0.52 0.12 0.37 0.72]);
    end
    for s = 1:numel(cfg.Speeds)
        x = NaN(1, numel(cfg.E4CPI));
        y = NaN(1, numel(cfg.E4CPI));
        for q = 1:numel(cfg.E4CPI)
            r = getRow(ALL, IDX, 'E4a', sprintf('cpi%.3f', cfg.E4CPI(q)), 'ws', 'WalkingSpeed', snr, cfg.Speeds(s));
            if ~isempty(r)
                x(q) = r.CPI_speed_act;
                y(q) = r.(fld);
            end
        end
        if panel == 2, x = x*160*cfg.Speeds(s); end
        if isBias, y = max(min(y, 30), -30); else, y = min(y, 20); end
        plot(ax, x, y, '-o', 'Color', cols(s, :), 'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', sprintf('%.1f m/s', cfg.Speeds(s)));
    end
    set(ax, 'XScale', 'log');
    if panel == 1
        xlabel(ax, 'Window length (s)');
    else
        xlabel(ax, 'Cycles of torso Doppler in the window (T x 160 v)');
        xline_safe(ax, 25);
    end
    if isBias
        ylabel(ax, 'Mean error (% of true speed)');
        yline_safe(ax, 0);
    else
        ylabel(ax, 'RMSE (% of true speed)');
        yline_safe(ax, cfg.Target_pct);
    end
    if panel == 2
        lg = legend(ax, 'show', 'Location', 'eastoutside', 'FontSize', 8);
        try, title(lg, 'Speed'); catch, end
        pos = get(ax, 'Position');
        set(ax, 'Position', [pos(1) pos(2) 0.30 pos(4)]);
    end
end
what = 'RMSE';
if isBias, what = 'Mean error (bias)'; end
annotation(fig, 'textbox', [0.02 0.90 0.96 0.09], 'String', sprintf('E4 walking speed %s against window. Gamma = %g dB, 75%% overlap, 3 passes, band 0.5x to 2x, %d trials, all %d speeds. Dashed line: %s.', ...
    what, snr, cfg.nTrials, numel(cfg.Speeds), ifElse(isBias, 'zero', '5 % target, and 25 cycles on the right')), 'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');
saveLightFigure(fig, fullfile(plotDir, sprintf('E4_%s_vs_window_G%g', ifElse(isBias, 'bias', 'rmse'), snr)));
close(fig);
end

function plotE2(cfg, plotDir, ALL, IDX, snr)
fig = newFig();
rules = {'frame', 'frame_cscaled', 'oracle', 'est'};
labels = {'one frame, c fixed', 'one frame, c shrinks with hop (v3)', 'lag 0.25 x true Tc', 'lag 0.25 x estimated Tc'};
for q = 1:numel(cfg.E2FixedLags)
    rules{end+1} = sprintf('fix%g', cfg.E2FixedLags(q));        %#ok<AGROW>
    labels{end+1} = sprintf('fixed lag %g s', cfg.E2FixedLags(q)); %#ok<AGROW>
end
cols = lines(numel(rules));
hops = cfg.E2Hops;
for panel = 1:2
    if panel == 1
        ax = styledAxes(fig, [0.08 0.14 0.38 0.68]);
    else
        ax = styledAxes(fig, [0.55 0.14 0.38 0.68]);
    end
    for r = 1:numel(rules)
        y = NaN(1, numel(hops));
        for h = 1:numel(hops)
            nm = sprintf('hop%g_%s', hops(h), rules{r});
            if panel == 1
                me = getVec(ALL, IDX, 'E2', nm, 'main', 'StepTime', snr, cfg.Speeds, 'MissRate_pct') + getVec(ALL, IDX, 'E2', nm, 'main', 'StepTime', snr, cfg.Speeds, 'ExtraRate_pct');
                y(h) = meanFinite(me);
            else
                y(h) = meanFinite(getVec(ALL, IDX, 'E2', nm, 'main', 'StrideDuration', snr, cfg.Speeds, 'RMSE_pct'));
            end
        end
        if panel == 2, y = min(y, 40); end
        plot(ax, 1:numel(hops), y, '-o', 'Color', cols(r, :), 'LineWidth', 1.3, 'MarkerSize', 5, 'DisplayName', labels{r});
    end
    set(ax, 'XTick', 1:numel(hops), 'XTickLabel', compose_hops(hops));
    xlim(ax, [0.7 numel(hops) + 0.3]);
    xlabel(ax, 'Hop as a fraction of the window (overlap)');
    if panel == 1
        ylabel(ax, 'Mean miss + extra (% of true steps)');
    else
        ylabel(ax, 'Mean stride-duration RMSE (%, clipped at 40)');
        yline_safe(ax, cfg.Target_pct);
    end
end
legend(ax, 'show', 'Location', 'southoutside', 'FontSize', 8, 'NumColumns', 3);
annotation(fig, 'textbox', [0.02 0.90 0.96 0.09], 'String', sprintf('E2 peak rule against hop. Mean over %d speeds, step window 0.186 s, Gamma = %g dB, c scale 0.5, %d trials.', ...
    numel(cfg.Speeds), snr, cfg.nTrials), 'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');
saveLightFigure(fig, fullfile(plotDir, sprintf('E2_peak_rule_vs_hop_G%g', snr)));
close(fig);
end

function c = compose_hops(hops)
c = cell(1, numel(hops));
for k = 1:numel(hops)
    c{k} = sprintf('%g (%g%%)', hops(k), 100*(1 - hops(k)));
end
end

function plotE3(cfg, plotDir, ALL, IDX, snr)
fig = newFig();
ax = styledAxes(fig, [0.09 0.13 0.62 0.70]);
ests = {'acorr', 'acorr2', 'twopass'};
labels = {'autocorrelation (0.6 to 3 s)', 'autocorrelation x2 (step period)', 'two-pass'};
cols = lines(3);
allT = [];
for q = 1:numel(ests)
    xT = NaN(1, numel(cfg.Speeds));
    yE = NaN(1, numel(cfg.Speeds));
    sE = NaN(1, numel(cfg.Speeds));
    for s = 1:numel(cfg.Speeds)
        r = getRow(ALL, IDX, 'E3', ['Tc_' ests{q}], 'main', 'StepTime', snr, cfg.Speeds(s));
        if ~isempty(r)
            xT(s) = r.TcTrue_s;
            yE(s) = r.TcEst_s;
            sE(s) = r.TcEstStd_s;
        end
    end
    allT = [allT, xT, yE];                                    %#ok<AGROW>
    errorbar(ax, xT, yE, sE, 'o', 'Color', cols(q, :), 'LineWidth', 1.2, 'MarkerSize', 5, 'DisplayName', labels{q});
end
allT = allT(isfinite(allT));
if isempty(allT), allT = [0.5 2.5]; end
lim = [min(allT)*0.9, max(allT)*1.1];
plot(ax, lim, lim, 'k--', 'LineWidth', 1, 'DisplayName', 'estimate = true Tc');
xlim(ax, lim);
ylim(ax, lim);
xlabel(ax, 'True Tc (s), one point per speed');
ylabel(ax, 'Estimated Tc (s), mean and std over trials');
lg = legend(ax, 'show', 'Location', 'eastoutside', 'FontSize', 8);
try, title(lg, 'Estimator'); catch, end
annotation(fig, 'textbox', [0.02 0.90 0.96 0.09], 'String', sprintf('E3 estimated Tc against true Tc. Step window 0.186 s, Gamma = %g dB, %d trials, speeds %.1f to %.1f m/s.', ...
    snr, cfg.nTrials, min(cfg.Speeds), max(cfg.Speeds)), 'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');
saveLightFigure(fig, fullfile(plotDir, sprintf('E3_Tc_estimate_G%g', snr)));
close(fig);
end

function plotPass(cfg, plotDir, PS, e, ext, titleText, nT)
isWS = strcmp(ext, 'ws');
names = {};
for i = 1:numel(PS)
    if strcmp(PS(i).Experiment, e) && strcmp(PS(i).Extractor, ext) && ~any(strcmp(names, PS(i).Variant))
        names{end+1} = PS(i).Variant;                          %#ok<AGROW>
    end
end
if isempty(names)
    return;
end
if isWS, snrs = cfg.SpeedSNRs; else, snrs = cfg.StepSNRs; end
if strcmp(e, 'E4v'), snrs = cfg.SpeedSNRs; end
names = [{'BASE'}, names];
vals = NaN(numel(names), numel(snrs));
for q = 1:numel(names)
    for k = 1:numel(snrs)
        if q == 1
            ps = getPS(PS, 'E0', 'BASE', ext, snrs(k));
        else
            ps = getPS(PS, e, names{q}, ext, snrs(k));
        end
        if isempty(ps), continue; end
        if isWS
            vals(q, k) = ps.WalkingSpeed_pass;
        else
            vals(q, k) = ps.Total_of_4;
        end
    end
end
fig = newFig();
ax = styledAxes(fig, [0.30 0.13 0.52 0.72]);
b = barh(ax, flipud(vals));
cols = lines(numel(snrs));
for k = 1:numel(b)
    try, set(b(k), 'FaceColor', cols(k, :), 'EdgeColor', 'none'); catch, end
    try, set(b(k), 'DisplayName', sprintf('Gamma %g dB', snrs(k))); catch, end
end
set(ax, 'YTick', 1:numel(names), 'YTickLabel', flipud(names(:)), 'TickLabelInterpreter', 'none', 'FontSize', 8);
nSp = numel(cfg.Speeds);
if isWS, topv = nSp; what = 'walking speed conditions passing'; else, topv = 4*nSp; what = 'conditions passing (all four parameters)'; end
xlim(ax, [0 topv]);
xlabel(ax, sprintf('Number of %s (out of %d, RMSE below %g %%)', what, topv, cfg.Target_pct));
legend(ax, 'show', 'Location', 'eastoutside', 'FontSize', 8);
annotation(fig, 'textbox', [0.02 0.90 0.96 0.09], 'String', sprintf('%s. Passing conditions per variant. %d speeds, %d trials, step window 0.186 s (speed window 0.372 s unless the variant changes it). First bar = baseline.', ...
    titleText, nSp, nT), 'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold');
saveLightFigure(fig, fullfile(plotDir, sprintf('Pass_%s_%s', e, ext)));
close(fig);
end

function xline_safe(ax, v)
yl = get(ax, 'YLim');
plot(ax, [v v], yl, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

function yline_safe(ax, v)
xl = get(ax, 'XLim');
plot(ax, xl, [v v], '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

function s = ifElse(tf, a, b)
if tf, s = a; else, s = b; end
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
    row.RMSE_abs = meanFinite(trialRMSE);                    % existing stride length sweeps average trial RMSE
else
    row.RMSE_abs = pooled;
end
row.RMSE_pct = 100*row.RMSE_abs/row.TruthMean;
row.RMSEpooled_pct = 100*pooled/row.TruthMean;
row.MeanErr_pct = 100*mean(allE)/row.TruthMean;
row.TrialMean_pct = meanFinite(trialPct);
row.TrialStd_pct = stdFinite(trialPct);                      % spread across trials, used for error bars
row.Pass = double(row.RMSE_pct < cfg.Target_pct);
end

