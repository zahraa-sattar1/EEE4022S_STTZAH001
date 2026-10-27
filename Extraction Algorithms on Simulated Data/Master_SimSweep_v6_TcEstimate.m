%% CPI and SNR design sweep
% Zahraa Sattar | EEE4022S
%
% Compare CPI and SNR choices using an estimated gait cycle.
% Retain the step, speed and stride results used in Chapter 4.

function Master_SimSweep_v6_TcEstimate(userCfg)

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
    if ~isfield(userCfg, 'Speeds'),         cfg.Speeds = [1.2 0.6]; end
    if ~isfield(userCfg, 'nTrials'),        cfg.nTrials = 3; end
    if ~isfield(userCfg, 'nValidate'),      cfg.nValidate = 3; end
    if ~isfield(userCfg, 'RunDiagnostics'), cfg.RunDiagnostics = false; end
end
tRun = tic;                                                  % total run time, written to the log at the end

global StepRefineMode StepGapFill TcCfg
StepRefineMode = cfg.StepRefine;                             % read by runStep and handed to the v5 extractor
StepGapFill = double(cfg.StepGapFill);
TcCfg = struct('Method', cfg.TcMethod, 'Nominal', cfg.TcNominal, 'Range', cfg.TcRange, 'TorsoCPI', cfg.TorsoCPIreq);   % read by runStep
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
if ~isempty(cfg.PrevRun) && ~exist(cfg.PrevRun, 'dir')
    cfg.PrevRun = fullfile(codeFolder, cfg.PrevRun);         % previous run folder, looked for next to the scripts
end

if isempty(cfg.ResumeFolder)
    outDir = fullfile(codeFolder, ['SimSweep_TcEst_' datestr(now, 'yyyymmdd_HHMMSS')]);   % new folder name for v6
    mkdir(outDir);                                           % always a new folder, nothing is overwritten
else
    outDir = cfg.ResumeFolder;
    if ~exist(outDir, 'dir')
        outDir = fullfile(codeFolder, cfg.ResumeFolder);
    end
end
cfg.ResumeFolder = outDir;
if ~isfile(fullfile(outDir, 'Config.mat'))
    save(fullfile(outDir, 'Config.mat'), 'cfg');
end

logMsg(outDir, 'Output folder: %s', outDir);
logMsg(outDir, 'R0 = %g m, noise mode = %s, CPI sweep SNR = %g dB, trials = %d', ...
    cfg.R0_m, cfg.NoiseMode, cfg.CPIfixedSNR, cfg.nTrials);
logMsg(outDir, 'Script: Master_SimSweep_v6_TcEstimate. MATLAB/Octave version: %s', version());   % version of the software that ran this
logMsg(outDir, 'Tc method = %s, nominal Tc = %g s, accepted range %g to %g s, quick mode = %d', cfg.TcMethod, cfg.TcNominal, cfg.TcRange(1), cfg.TcRange(2), cfg.Quick);
logMsg(outDir, 'Random seeds: seed = %d + 1000000*stageCode + 1000*round(10*speed) + trial, set with rng(seed,''twister''). Stage codes S1=1 S2=2 S3a=3 S3e=4 S3b=5 S3c=6 S4=7. Every seed is in Seeds_Used.csv. Simulator seed = %d + round(100*speed).', ...
    cfg.Seed, cfg.Seed);
logMsg(outDir, 'Previous run for the comparison: %s', cfg.PrevRun);

P = makeWorkingCopies(cfg, codeFolder, outDir);
patchLogFile = fullfile(outDir, 'Patch_Log.txt');
if isfile(patchLogFile)
    pl = strsplit(strtrim(fileread(patchLogFile)), sprintf('\n'));   % copy the patch log into the run log
    for q = 1:numel(pl)
        logMsg(outDir, 'PATCH %s', pl{q});
    end
end
writeAssumptions(outDir, cfg);
selfTest(P, cfg, outDir);
selfTestTc(P, cfg, outDir);                                  % stops the run if the Tc patch is not working

tStart = tic;
for sIdx = 1:numel(cfg.Speeds)
    v0 = cfg.Speeds(sIdx);
    logMsg(outDir, '===== v0 = %.2f m/s (%d of %d), elapsed %.1f min =====', ...
        v0, sIdx, numel(cfg.Speeds), toc(tStart)/60);
    f1 = fullfile(outDir, sprintf('S1_v0_%.2f.mat', v0));
    f2 = fullfile(outDir, sprintf('S2_v0_%.2f.mat', v0));
    f3 = fullfile(outDir, sprintf('S3_v0_%.2f.mat', v0));
    f4 = fullfile(outDir, sprintf('S4_v0_%.2f.mat', v0));
    if isfile(f1) && isfile(f2) && isfile(f4) && (isfile(f3) || ~cfg.RunDiagnostics)
        logMsg(outDir, 'Already finished, skipping.');
        continue;
    end
    try
        D = getSpeedData(P, cfg, v0);                        % one simulator run gives the clean signal for every trial
    catch ME
        logMsg(outDir, 'Simulator failed at v0 = %.2f: %s', v0, ME.message);
        continue;
    end

    try
        if isfile(f1)
            s = load(f1);
            rows1 = s.rows1;
        else
            rows1 = stage1(P, cfg, D, outDir);
            save(f1, 'rows1');
        end
        sel = selectSettings(rows1, [], cfg);                % best CPI per parameter at the fixed SNR

        if isfile(f2)
            s = load(f2);
            rows2 = s.rows2;
            sel = s.sel;
        else
            rows2full = stage2(P, cfg, D, sel, outDir);
            sel = selectSettings(rows1, rows2full, cfg);     % best CPI and SNR chosen together from the full grid
            rows2 = keepSelectedCPIs(rows2full, sel);        % the CSV, plots and report show the chosen CPIs only
            logMsg(outDir, 'Chosen: speed %.3f s / %g dB, step %.3f s / %g dB, stride dur %.3f s / %g dB, stride len %.3f + %.3f s / %g dB', ...
                sel.ws.cpi, sel.ws.snr, sel.st.cpi, sel.st.snr, sel.sd.cpi, sel.sd.snr, sel.sl.cs, sel.sl.ct, sel.sl.snr);
            save(f2, 'rows2', 'rows2full', 'sel');
        end

        if cfg.RunDiagnostics && ~isfile(f3)
            rows3 = stage3(P, cfg, D, sel, outDir);
            save(f3, 'rows3');
        end

        if ~isfile(f4)
            [rows4, rows4old] = stage4(P, cfg, D, sel, outDir);
            save(f4, 'rows4', 'rows4old');
        end
    catch ME
        logMsg(outDir, 'Speed %.2f stopped early: %s', v0, ME.message);
        for q = 1:min(3, numel(ME.stack))
            logMsg(outDir, '   at %s line %d', ME.stack(q).name, ME.stack(q).line);
        end
    end
end

logMsg(outDir, 'All speeds done in %.1f min. Writing results.', toc(tStart)/60);
try
    assembleResults(cfg, outDir);
catch ME
    logMsg(outDir, 'Writing the results failed: %s', ME.message);
end
logMsg(outDir, 'Finished. Open TcEstimate_Report.txt and Diagnosis_Report.txt in %s', outDir);
logMsg(outDir, 'Total run time %.1f min', toc(tRun)/60);    % runtime for the log
end

function cfg = defaultConfig()
cfg.Speeds = [0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];      % m/s
cfg.nTrials = 15;                                            % noise draws per condition
cfg.nValidate = 15;                                          % fresh draws for Stage 4
cfg.R0_m = 10;                                               % simulator range, v11 default (v6 Section 7 used 50)
cfg.AmplitudeMode = 'A';                                     % Chen mean sqrt(RCS) weights
cfg.CPIreq = [0.025 0.075 0.1 0.2 0.4];                      % actual windows 0.046 0.093 0.186 0.372 0.743 s
cfg.SNRlist = 20:10:90;                                      % dB
cfg.CPIfixedSNR = 30;                                        % dB used in the CPI sweep
cfg.Target_pct = 5;                                          % goal for RMSE %
cfg.TieTol_pct = 0.25;
cfg.MaxMissExtra_pct = 10;                                   % step time CPI must keep miss + extra below this if any CPI can
cfg.NoiseMode = 'extractorN';                                % 'extractorN' or 'simulatorN' (original: 4096 window below 1 m/s)
cfg.PthList = [-15 -20 -25 -30 -35 -40 -45];                 % dB, for the Pth decoupling test
cfg.DiagGammas = [30 60];                                    % dB, SNRs used for the Pth test
cfg.HopList = [0.25 0.125 0.0625];                           % step spectrogram hop as fraction of window
cfg.HopSNR = 60;                                             % dB, high so noise is small and timing steps dominate
cfg.CPeakList = [1 0.5 0.25];                                % multipliers on the step peak-drop rule c, for the Stage 3e test
cfg.RunDiagnostics = true;
cfg.MakePlots = true;
cfg.SimName = 'RadarSimulator_v11_torso_legs_arms_head_buffered_noplots';
cfg.WSName = 'Finding_WalkingSpeed_v4_noplots';
cfg.STName = 'Finding_StepRecognition_v5_noplots';
cfg.StepRefine = 'centroid';                                 % 'centroid', 'parabola' or 'none' (Step 9c of the v5 extractor)
cfg.StepGapFill = true;                                      % Step 9b of the v5 extractor
cfg.CompareOldRules = true;                                  % Stage 4 also scores the v3 rules on the same noise
cfg.TcMethod = 'twopass';                                    % v6: 'twopass', 'torso' or 'true'
cfg.TcNominal = 1.2;                                         % s, Tc used in pass 1 and as the fallback
cfg.TcRange = [0.3 4];
cfg.TorsoCPIreq = 0.2;
cfg.Seed = 20266001;                                         % base random seed
cfg.Quick = false;                                           % true = 2 speeds, 3 trials, no diagnostics
cfg.PrevRun = 'SimSweep_20261008_211723';                    % previous run (true Tc), folder name next to the scripts or a full path
cfg.CodeFolder = '';
cfg.ResumeFolder = '';
cfg.CopyHook = [];                                           % only used for testing outside MATLAB
end

function writeAssumptions(outDir, cfg)
f = fullfile(outDir, 'Run_Assumptions.txt');
if isfile(f)
    return;
end
L = {
    'Run assumptions'
    sprintf('Simulator: %s, amplitude mode %s, R0 = %g m', cfg.SimName, cfg.AmplitudeMode, cfg.R0_m)
    'The clean signal is made once per speed. Noise is added here with the same equation as the simulator.'
    'Noise power = Nw*10^(-SNR/10), complex Gaussian, Nw = FFT window length.'
    sprintf('Noise mode %s: extractorN uses the extractor window for Nw. simulatorN copies the original 4096 window below 1 m/s.', cfg.NoiseMode)
    'The Pth rule in the step script is unchanged: Pth = max(10 - SNR, -40) dB unless a Pth is forced.'
    'Step time, miss and extra scoring copies ChenAmp: edge exclusion, one-to-one match, tolerance min(0.2, Tc/4).'
    'Stride length and stride duration follow the existing definitions (pooled / trial-averaged as in ChenAmp).'
    sprintf('Random seeds: every noise draw is seeded (base %d, see Seeds_Used.csv). v4 and v5 used no seed.', cfg.Seed)
    sprintf('Tc: the step extractor is given Tc_est (method %s), not the true Tc. The scorer keeps the true Tc.', cfg.TcMethod)
    'Results are simulated data only.'};
fid = fopen(f, 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
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
P.sim = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.SimName '.m']), fullfile(codeDir, 'WC_Simulator.m'), 'sim', cfg);
P.ws = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.WSName '.m']), fullfile(codeDir, 'WC_WalkingSpeed.m'), 'ws', cfg);
P.st = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.STName '.m']), fullfile(codeDir, 'WC_StepRecognition.m'), 'st', cfg);
end

function dest = copyOne(src, dest, kind, cfg)
if isfile(dest)
    return;                                                  % never overwrite, reuse on resume
end
assert(isfile(src), 'Missing file: %s', src);
txt = fileread(src);
patchLog = fullfile(fileparts(fileparts(dest)), 'Patch_Log.txt');   % v6: every patched line is written here
fidP = fopen(patchLog, 'a');
if fidP > 0
    fprintf(fidP, '==== %s (%s) ====\n', dest, src);          % heading for this working copy
    fclose(fidP);
end
logPatchMatches(patchLog, 'clear/close all removed', txt, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*');   % log before removing
txt = regexprep(txt, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*(?:\r?\n|$)', '');   % no clear all / close all
if ~strcmp(kind, 'sim')
    logPatchMatches(patchLog, 'simulator call removed', txt, '(?m)^[ \t]*RadarSimulator_[A-Za-z0-9_]+;[^\r\n]*');   % log before removing
    txt = regexprep(txt, '(?m)^[ \t]*RadarSimulator_[A-Za-z0-9_]+;[^\r\n]*(?:\r?\n|$)', '');   % simulator is called separately
end
if strcmp(kind, 'st')
    plotMarker = 'figure(''Name'',''Step 6 verification - remove outliers'');';
    a = strfind(txt, plotMarker);
    b = strfind(txt, '%% Step 7');
    assert(numel(a) == 1 && numel(b) == 1 && b > a, 'Step 6 figure not found in the step script.');
    logPatch(patchLog, 'Step 6 figure removed', sprintf('%d characters of figure code between the Step 6 marker and %%%% Step 7', b - a), '(block removed)');
    txt = [txt(1:a-1), txt(b:end)];                          % cut the Step 6 figure
    hopOld = 'noverlap = N*3/4;';
    assert(numel(strfind(txt, hopOld)) == 1, 'Overlap line not found in the step script.');
    hopNew = sprintf(['global StepHopFraction\n' ...
        '        hopFrac_ = 0.25;\n' ...
        '        if ~isempty(StepHopFraction), hopFrac_ = StepHopFraction; end\n' ...
        '        noverlap = N - max(1, round(N*hopFrac_));']);   % default is the original 75 percent overlap
    logPatch(patchLog, 'hop size made a setting', hopOld, hopNew);   % v6: logged
    txt = strrep(txt, hopOld, hopNew);
    cpkPattern = 'cPeak\s*=\s*DopplerRes_Hz/freq_high_Yardibi;[^\r\n]*';
    assert(numel(regexp(txt, cpkPattern, 'match')) == 1, 'cPeak line not found in the step script.');
    cpkNew = sprintf(['global StepHopFraction\n' ...
        'cPeak = DopplerRes_Hz/freq_high_Yardibi;\n' ...
        'if ~isempty(StepHopFraction), cPeak = cPeak*StepHopFraction/0.25; end']);   % peak drop rule scales with the hop
    if isempty(strfind(txt, 'cPeakScale'))                       % v4 has no scale of its own, add one (v4b already has cPeakScale)
        cpkNew = [cpkNew sprintf('\nglobal StepCPeakScale\nif ~isempty(StepCPeakScale), cPeak = cPeak*StepCPeakScale; end')];
    end
    logPatch(patchLog, 'peak rule scale made a setting', regexp(txt, cpkPattern, 'match', 'once'), cpkNew);   % v6: logged
    txt = regexprep(txt, cpkPattern, cpkNew);

    tcPatches = {'medianWindow_s = Tc/8;',                    'medianWindow_s = Tc_est/8;'; ...            % Step 4 median filter window
                 'smoothWindow_s = Tc/8;',                    'smoothWindow_s = Tc_est/8;'; ...            % Step 7 smoothing window
                 'stepIdx9, Tc, stepRefineMode)',             'stepIdx9, Tc_est, stepRefineMode)'};        % Step 9c search window Tc/4
    for q = 1:size(tcPatches, 1)
        assert(numel(strfind(txt, tcPatches{q, 1})) == 1, ...
            'Tc pattern "%s" was not found exactly once in %s. The extractor has changed, so the Tc patch was not applied.', tcPatches{q, 1}, src);   % stop if not found
        oldLine = lineContaining(txt, tcPatches{q, 1});      % the whole line, for the log
        txt = strrep(txt, tcPatches{q, 1}, tcPatches{q, 2});
        logPatch(patchLog, 'Tc replaced by Tc_est', oldLine, lineContaining(txt, tcPatches{q, 2}));
    end
    cut = strfind(txt, sprintf('\nfunction '));              % script part ends where the local functions start
    assert(~isempty(cut), 'No local functions found in the step script.');
    bodyTxt = txt(1:cut(1));
    left = regexp(bodyTxt, '[^\n]*(?<![A-Za-z0-9_])Tc(?![A-Za-z0-9_])[^\n]*', 'match');   % any use of the true Tc name left in the script part
    if ~isempty(left)                                        % stop if the true Tc is still used anywhere in the script part
        error('The step script still uses Tc in %d line(s) after the patch: %s', numel(left), strtrim(left{1}));
    end
    assert(numel(regexp(bodyTxt, '(?<![A-Za-z0-9_])Tc_est(?![A-Za-z0-9_])', 'match')) == 3, 'Tc_est should appear exactly 3 times in the patched step script.');
end
if ~isempty(cfg.CopyHook)
    txt = cfg.CopyHook(txt, kind);
end
fid = fopen(dest, 'w');
assert(fid >= 0, 'Could not create %s', dest);
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function D = getSpeedData(m_P, m_cfg, v0)
SimulatorVerbose = false;                                    %#ok<NASGU>
SimulatorMakePlots = false;                                  %#ok<NASGU>
SimulatorR0_m = m_cfg.R0_m;                                  %#ok<NASGU>
AmplitudeMode = m_cfg.AmplitudeMode;                         %#ok<NASGU>
SNR_dB = 40;                                                 %#ok<NASGU>
CPI_s = 0.1;                                                 %#ok<NASGU>
rng(m_cfg.Seed + round(v0*100), 'twister');                 % v6: fixed seed for anything random inside the simulator
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

function Nw = noiseWindow(cfg, v0, cpiReq, fs)
Nw = 2^nextpow2(round(cpiReq*fs));                           % extractor window length
if strcmp(cfg.NoiseMode, 'simulatorN') && v0 < 1.0
    Nw = 2^nextpow2(round(0.075*fs));                        % original simulator forces 0.075 s below 1 m/s
end
end

function x = addNoise(clean, Nw, snr_dB, z)
noisePower = 10^(-(snr_dB - 10*log10(Nw))/10);               % same equation as the simulator
x = clean + sqrt(noisePower/2) * z;                          % z is one unit complex noise draw
end

function cpiAct = actualCPI(cpiReq, fs)
cpiAct = 2^nextpow2(round(cpiReq*fs)) / fs;
end

function selfTest(P, cfg, outDir)
try
    D = getSpeedData(P, cfg, 1.2);
    snr = 30;
    Nw = D.Nsim;
    z = randn(1, numel(D.clean)) + 1i*randn(1, numel(D.clean));
    x = addNoise(D.clean, Nw, snr, z);
    mine = mean(abs(x - D.clean).^2);
    expected = 10^(-(snr - 10*log10(Nw))/10);
    simNoise = mean(abs(D.noiseSim).^2);
    simExpected = 10^(-(40 - 10*log10(Nw))/10);              % the simulator run above used 40 dB
    ok1 = abs(mine/expected - 1) < 0.05;
    ok2 = abs(simNoise/simExpected - 1) < 0.05;
    logMsg(outDir, 'Self test: injected noise power ratio %.3f (%s), simulator noise power ratio %.3f (%s)', ...
        mine/expected, passText(ok1), simNoise/simExpected, passText(ok2));
catch ME
    logMsg(outDir, 'Self test could not run: %s', ME.message);
end
end

function s = passText(ok)
if ok
    s = 'PASS';
else
    s = 'CHECK';
end
end

function m_o = runSpeed(m_P, m_x, m_D, m_cpi)
m_o = struct('ok', false, 'meanSpeed', NaN, 'vel', [], 'tv', [], 'cpiAct', NaN, 'err', '');
try
    RdrRxSignal = m_x;                                       %#ok<NASGU>
    fs = m_D.fs;                                             %#ok<NASGU>
    fc = m_D.fc;                                             %#ok<NASGU>
    c = m_D.c;                                               %#ok<NASGU>
    v0 = m_D.v0;                                             %#ok<NASGU>
    CPI_s = m_cpi;                                           %#ok<NASGU>
    evalc('run(m_P.ws);');
    m_o.meanSpeed = meanSpeed_mps;
    m_o.vel = velocityTrack_mps;
    m_o.tv = final.time_axis;
    m_o.cpiAct = CPI_actual_s_WS;
    m_o.ok = isfinite(m_o.meanSpeed);
catch ME
    m_o.err = ME.message;
end
end

function m_o = runStepPass(m_P, m_x, m_D, m_cpi, m_snr, m_pth, m_hop, m_cps, m_Tcx)
global StepHopFraction StepCPeakScale StepRefineMode StepGapFill
if nargin < 8
    m_cps = [];
end
StepHopFraction = m_hop;
StepCPeakScale = m_cps;
m_o = struct('ok', false, 'det', [], 'key', [], 'sd', [], 'detRaw', [], 'sdAll', [], 'cpiAct', NaN, 'sc', [], 'scRaw', [], 'err', '', ...
    'TcEst', NaN, 'TcFallback', 0, 'medWin', NaN, 'smWin', NaN);   % v6: the Tc given to the extractor and the two windows it produced
try
    RdrRxSignal = m_x;                                       %#ok<NASGU>
    fs = m_D.fs;                                             %#ok<NASGU>
    Tc = NaN;                                                %#ok<NASGU>
    Tc_est = m_Tcx;                                          %#ok<NASGU>
    CPI_s = m_cpi;                                           %#ok<NASGU>
    SNR_dB = m_snr;                                          %#ok<NASGU>
    if ~isempty(m_pth)
        PthOverride = m_pth;                                 %#ok<NASGU>
    end
    if ~isempty(m_cps)
        cPeakScale = m_cps;                                  %#ok<NASGU>
    end
    if ~isempty(StepRefineMode)
        stepRefineMode = StepRefineMode;                     %#ok<NASGU>
    end
    if ~isempty(StepGapFill)
        stepGapFill = StepGapFill;                           %#ok<NASGU>
    end
    evalc('run(m_P.st);');
    m_o.medWin = medianWindow_s;                             % v6: kept for the self test (should be Tc_est/8)
    m_o.smWin = smoothWindow_s;                              % v6: kept for the self test (should be Tc_est/8)
    m_o.det = detectedStepTimes_s(:).';                      % final step times (refined when the v5 extractor is used)
    if exist('stepFrameTimes_s', 'var')
        m_o.key = stepFrameTimes_s(:).';                     % frame time stamps of the same steps, they decide who is inside the walk
    else
        m_o.key = m_o.det;                                   % v4 and v4b: the step times are the frame time stamps
    end
    if exist('detectedStepTimes_raw_s', 'var')
        m_o.detRaw = detectedStepTimes_raw_s(:).';           % the v4b step list, for the before and after comparison
    else
        m_o.detRaw = m_o.det;
    end
    m_o.cpiAct = CPI_actual_s;
    edge = 0.1 + m_o.cpiAct/2;                               % same walk interior as the step time scoring
    lo = m_D.buffer + edge;
    hi = m_D.buffer + m_D.simT - edge;
    tIn = m_o.det(m_o.key > lo & m_o.key < hi);              % steps inside the walk only
    if numel(tIn) >= 3
        m_o.sd = (tIn(3:end) - tIn(1:end-2)) - m_D.Tc;       % stride durations, kept as (stride - Tc), the same way as before
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
StepHopFraction = [];
StepCPeakScale = [];
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

function row = blankRow()
row = struct('Stage', '', 'Parameter', '', 'Speed_mps', NaN, ...
    'CPI_speed_req', NaN, 'CPI_speed_act', NaN, 'CPI_step_req', NaN, 'CPI_step_act', NaN, ...
    'SNR_dB', NaN, 'Pth_dB', NaN, 'HopFraction', NaN, 'CPeakScale', NaN, 'NoiseMode', '', ...
    'Trials', NaN, 'FailedTrials', NaN, 'N_compared', NaN, 'TruthMean', NaN, ...
    'RMSE_abs', NaN, 'RMSE_pct', NaN, 'RMSEpooled_pct', NaN, 'MeanErr_pct', NaN, ...
    'TrialMean_pct', NaN, 'TrialStd_pct', NaN, ...
    'MAE_abs', NaN, 'MAPE_pct', NaN, 'Accuracy_pct', NaN, ...
    'MissRate_pct', NaN, 'ExtraRate_pct', NaN, 'StepOffset_ms', NaN, ...
    'MissCentred_pct', NaN, 'ExtraCentred_pct', NaN, 'PerfectTrials_pct', NaN, 'Pass', NaN, ...
    'Tc_est_mean_s', NaN, 'Tc_true_s', NaN, 'Tc_err_pct', NaN, 'Tc_fallbacks', NaN, 'Tc_worst_abs_err_pct', NaN);   % v6: Tc columns at the end
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
row.MAE_abs = mean(abs(allE));                               % mean of the absolute errors over all draws
row.MAPE_pct = 100*row.MAE_abs/row.TruthMean;                % mean absolute error as a % of the mean true value
row.Accuracy_pct = 100 - row.MAPE_pct;                       % accuracy = 100 minus the MAPE
row.RMSEpooled_pct = 100*pooled/row.TruthMean;
row.MeanErr_pct = 100*mean(allE)/row.TruthMean;
row.TrialMean_pct = meanFinite(trialPct);
row.TrialStd_pct = stdFinite(trialPct);                      % spread across trials, used for error bars
row.Pass = double(row.RMSE_pct < cfg.Target_pct);
end

function row = rowWS(cfg, D, W, cpi, snr, stage)
nT = numel(W);
E = cell(nT, 1);
G = cell(nT, 1);
for k = 1:nT
    if W{k}.ok
        E{k} = W{k}.meanSpeed - D.v0;
        G{k} = D.v0;
    else
        E{k} = NaN;
        G{k} = NaN;
    end
end
row = blankRow();
row.Stage = stage;
row.Parameter = 'WalkingSpeed';
row.Speed_mps = D.v0;
row.CPI_speed_req = cpi;
row.CPI_speed_act = actualCPI(cpi, D.fs);
row.SNR_dB = snr;
row.NoiseMode = cfg.NoiseMode;
row = fillStats(row, E, G, false, cfg);
end

function row = rowST(cfg, D, S, cpi, snr, stage, param, useRaw)
if nargin < 8
    useRaw = false;                                          % true = score the v3 rules (v4b step list, no stride trimming)
end
nT = numel(S);
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
    sdK = S{k}.sd;
    if useRaw
        scK = S{k}.scRaw;
        sdK = S{k}.sdAll;
    end
    if strcmp(param, 'StepTime')
        E{k} = scK.e;
        G{k} = scK.g;
    else
        E{k} = sdK;
        G{k} = D.Tc * ones(size(sdK));
    end
    miss(k) = scK.miss;
    extra(k) = scK.extra;
    missC(k) = scK.missC;
    extraC(k) = scK.extraC;
    perfect(k) = scK.perfect;
    off(k) = scK.offset;
end
row = blankRow();
row.Stage = stage;
row.Parameter = param;
row.Speed_mps = D.v0;
row.CPI_step_req = cpi;
row.CPI_step_act = actualCPI(cpi, D.fs);
row.SNR_dB = snr;
row.NoiseMode = cfg.NoiseMode;
row = fillStats(row, E, G, false, cfg);
row.MissRate_pct = meanFinite(miss);
row.ExtraRate_pct = meanFinite(extra);
row.MissCentred_pct = meanFinite(missC);
row.ExtraCentred_pct = meanFinite(extraC);
row.PerfectTrials_pct = 100*meanFinite(perfect);
row.StepOffset_ms = 1000*meanFinite(off);
row = addTcStats(row, S, D);                                 % v6: Tc columns
end

function row = rowSL(cfg, D, W, S, cs, ct, snr, stage)
nT = numel(S);
E = cell(nT, 1);
G = cell(nT, 1);
for k = 1:nT
    [E{k}, G{k}] = strideLengthErrors(W{k}, S{k}, D);
end
row = blankRow();
row.Stage = stage;
row.Parameter = 'StrideLength';
row.Speed_mps = D.v0;
row.CPI_speed_req = cs;
row.CPI_speed_act = actualCPI(cs, D.fs);
row.CPI_step_req = ct;
row.CPI_step_act = actualCPI(ct, D.fs);
row.SNR_dB = snr;
row.NoiseMode = cfg.NoiseMode;
row = fillStats(row, E, G, true, cfg);
row = addTcStats(row, S, D);                                 % v6: Tc columns (the Tc given to the step extractor)
end

function rows = addRow(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row;
end
end

function logFailures(outDir, label, objs)
for k = 1:numel(objs)
    if ~objs{k}.ok && ~isempty(objs{k}.err)
        logMsg(outDir, '   %s failed: %s', label, objs{k}.err);
        return;
    end
end
end

function rows = stage1(P, cfg, D, outDir)
logMsg(outDir, 'Stage 1: CPI sweep at %g dB', cfg.CPIfixedSNR);
cpis = cfg.CPIreq;
nC = numel(cpis);
nT = cfg.nTrials;
W = cell(nT, nC);
S = cell(nT, nC);
for k = 1:nT
    z = trialNoise(cfg, D, 'S1', k);                     % v6: fixed seed for this stage, speed and trial
    for c = 1:nC
        Nw = noiseWindow(cfg, D.v0, cpis(c), D.fs);
        x = addNoise(D.clean, Nw, cfg.CPIfixedSNR, z);
        W{k,c} = runSpeed(P, x, D, cpis(c));
        S{k,c} = runStep(P, x, D, cpis(c), cfg.CPIfixedSNR, [], []);
    end
    logMsg(outDir, '   trial %d of %d done', k, nT);
end
logFailures(outDir, 'walking speed', W(:));
logFailures(outDir, 'step recognition', S(:));
rows = [];
for c = 1:nC
    rows = addRow(rows, rowWS(cfg, D, W(:,c), cpis(c), cfg.CPIfixedSNR, 'S1'));
    rows = addRow(rows, rowST(cfg, D, S(:,c), cpis(c), cfg.CPIfixedSNR, 'S1', 'StepTime'));
    rows = addRow(rows, rowST(cfg, D, S(:,c), cpis(c), cfg.CPIfixedSNR, 'S1', 'StrideDuration'));
    for c2 = 1:nC
        rows = addRow(rows, rowSL(cfg, D, W(:,c), S(:,c2), cpis(c), cpis(c2), cfg.CPIfixedSNR, 'S1'));
    end
end
end

function idx = pickBest(rows, tieKey, useReliability, cfg, useExact)
if nargin < 5
    useExact = false;                                        % true = among tied settings prefer the most exact trials
end
rm = [rows.RMSE_pct];
fl = [rows.FailedTrials];
key = tieKey(:).';
ok = isfinite(rm);
if ~any(ok)
    idx = 1;
    return;
end
ok = ok & (fl == min(fl(ok)));
if useReliability
    rel = [rows.MissRate_pct] + [rows.ExtraRate_pct];
    rel(~isfinite(rel)) = 0;
    good = ok & (rel <= cfg.MaxMissExtra_pct);
    if any(good)
        ok = good;
    else
        ok = ok & (rel == min(rel(ok)));
    end
end
best = min(rm(ok));
cand = find(ok & (rm <= best + cfg.TieTol_pct));
if useExact
    pf = [rows.PerfectTrials_pct];
    pf(~isfinite(pf)) = 0;
    cand = cand(pf(cand) >= max(pf(cand)) - 1e-9);           % keep only the tied settings with the most exact trials
end
[~, i] = min(key(cand));
idx = cand(i);
end

function idx = pickRobust(rows, keyMat, tieKey, useReliability, useExact, cfg)
nR = numel(rows);
nD = size(keyMat, 2);
pos = zeros(nR, nD);
for d = 1:nD
    [~, ~, pos(:, d)] = unique(keyMat(:, d));                % position of each value along dimension d
end
rm = [rows.RMSE_pct];
rm(~isfinite(rm)) = inf;
rel = [rows.MissRate_pct] + [rows.ExtraRate_pct];
rel(~isfinite(rel)) = 0;                                     % speed and stride length rows have no miss or extra
worst = rm;
worstRel = rel;
for r = 1:nR
    for d = 1:nD
        for sgn = [-1 1]
            p = pos(r, :);
            p(d) = p(d) + sgn;
            j = find(all(pos == repmat(p, nR, 1), 2), 1);
            if ~isempty(j)
                worst(r) = max(worst(r), rm(j));
                worstRel(r) = max(worstRel(r), rel(j));
            end
        end
    end
end
fl = [rows.FailedTrials];
ok = isfinite(rm) & (rm < inf);
if ~any(ok)
    idx = 1;
    return;
end
ok = ok & (fl == min(fl(ok)));
if useReliability
    good = ok & (worstRel <= cfg.MaxMissExtra_pct);
    if any(good)
        ok = good;
    else
        ok = ok & (worstRel == min(worstRel(ok)));
    end
end
best = min(worst(ok));
cand = find(ok & (worst <= best + cfg.TieTol_pct));
if useExact
    pf = [rows.PerfectTrials_pct];
    pf(~isfinite(pf)) = 0;
    cand = cand(pf(cand) >= max(pf(cand)) - 1e-9);
end
[~, i] = min(tieKey(cand));
idx = cand(i);
end

function sel = selectSettings(rows1, rows2, cfg)
params = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideLength'};
sel = struct();
for p = 1:4
    r1 = rows1(strcmp({rows1.Parameter}, params{p}));
    switch params{p}
        case 'WalkingSpeed'
            i = pickBest(r1, [r1.CPI_speed_act], false, cfg);
            sel.ws.cpi = r1(i).CPI_speed_req;
        case 'StepTime'
            i = pickBest(r1, [r1.CPI_step_act], true, cfg);
            sel.st.cpi = r1(i).CPI_step_req;
        case 'StrideDuration'
            i = pickBest(r1, [r1.CPI_step_act], false, cfg);
            sel.sd.cpi = r1(i).CPI_step_req;
        case 'StrideLength'
            i = pickBest(r1, [r1.CPI_speed_act] + [r1.CPI_step_act], false, cfg);
            sel.sl.cs = r1(i).CPI_speed_req;
            sel.sl.ct = r1(i).CPI_step_req;
    end
end
sel.ws.snr = cfg.CPIfixedSNR;
sel.st.snr = cfg.CPIfixedSNR;
sel.sd.snr = cfg.CPIfixedSNR;
sel.sl.snr = cfg.CPIfixedSNR;
if isempty(rows2)
    return;
end
for p = 1:4                                                  % rows2 holds every CPI at every SNR, so the pair is chosen together
    r2 = rows2(strcmp({rows2.Parameter}, params{p}));
    switch params{p}
        case 'WalkingSpeed'
            i = pickRobust(r2, [[r2.CPI_speed_act]', [r2.SNR_dB]'], 1000*[r2.CPI_speed_act] + [r2.SNR_dB], false, false, cfg);
            sel.ws.cpi = r2(i).CPI_speed_req;
            sel.ws.snr = r2(i).SNR_dB;
        case 'StepTime'
            i = pickRobust(r2, [[r2.CPI_step_act]', [r2.SNR_dB]'], 1000*[r2.CPI_step_act] + [r2.SNR_dB], true, true, cfg);
            sel.st.cpi = r2(i).CPI_step_req;
            sel.st.snr = r2(i).SNR_dB;
        case 'StrideDuration'
            i = pickRobust(r2, [[r2.CPI_step_act]', [r2.SNR_dB]'], 1000*[r2.CPI_step_act] + [r2.SNR_dB], true, true, cfg);
            sel.sd.cpi = r2(i).CPI_step_req;
            sel.sd.snr = r2(i).SNR_dB;
        case 'StrideLength'
            i = pickRobust(r2, [[r2.CPI_speed_act]', [r2.CPI_step_act]', [r2.SNR_dB]'], 1000*([r2.CPI_speed_act] + [r2.CPI_step_act]) + [r2.SNR_dB], false, false, cfg);
            sel.sl.cs = r2(i).CPI_speed_req;
            sel.sl.ct = r2(i).CPI_step_req;
            sel.sl.snr = r2(i).SNR_dB;
    end
end
end

function rows = keepSelectedCPIs(rows2, sel)
keep = false(1, numel(rows2));
for r = 1:numel(rows2)
    switch rows2(r).Parameter
        case 'WalkingSpeed'
            keep(r) = rows2(r).CPI_speed_req == sel.ws.cpi;
        case 'StepTime'
            keep(r) = rows2(r).CPI_step_req == sel.st.cpi;
        case 'StrideDuration'
            keep(r) = rows2(r).CPI_step_req == sel.sd.cpi;
        case 'StrideLength'
            keep(r) = rows2(r).CPI_speed_req == sel.sl.cs && rows2(r).CPI_step_req == sel.sl.ct;
    end
end
rows = rows2(keep);
end

function rows = stage2(P, cfg, D, sel, outDir)
logMsg(outDir, 'Stage 2: SNR sweep at every CPI (%d CPIs x %d SNRs). Stage 1 favourites at %g dB: speed %.3f, step %.3f, stride dur %.3f, stride len %.3f / %.3f', ...
    numel(cfg.CPIreq), numel(cfg.SNRlist), cfg.CPIfixedSNR, sel.ws.cpi, sel.st.cpi, sel.sd.cpi, sel.sl.cs, sel.sl.ct);
cpis = cfg.CPIreq;
snrs = cfg.SNRlist;
nT = cfg.nTrials;
W = cell(nT, numel(snrs), numel(cpis));
S = cell(nT, numel(snrs), numel(cpis));
for k = 1:nT
    z = trialNoise(cfg, D, 'S2', k);                     % v6: fixed seed for this stage, speed and trial
    for g = 1:numel(snrs)
        for a = 1:numel(cpis)
            Nw = noiseWindow(cfg, D.v0, cpis(a), D.fs);
            x = addNoise(D.clean, Nw, snrs(g), z);
            W{k,g,a} = runSpeed(P, x, D, cpis(a));
            S{k,g,a} = runStep(P, x, D, cpis(a), snrs(g), [], []);
        end
    end
    logMsg(outDir, '   trial %d of %d done', k, nT);
end
rows = [];
for g = 1:numel(snrs)
    for a = 1:numel(cpis)
        rows = addRow(rows, rowWS(cfg, D, squeeze(W(:,g,a)), cpis(a), snrs(g), 'S2'));
        rows = addRow(rows, rowST(cfg, D, squeeze(S(:,g,a)), cpis(a), snrs(g), 'S2', 'StepTime'));
        rows = addRow(rows, rowST(cfg, D, squeeze(S(:,g,a)), cpis(a), snrs(g), 'S2', 'StrideDuration'));
        for b = 1:numel(cpis)
            rows = addRow(rows, rowSL(cfg, D, squeeze(W(:,g,a)), squeeze(S(:,g,b)), cpis(a), cpis(b), snrs(g), 'S2'));
        end
    end
end
end

function rows = stage3(P, cfg, D, sel, outDir)
logMsg(outDir, 'Stage 3: diagnostics');
rows = [];
nT = cfg.nTrials;

pths = cfg.PthList;
gams = cfg.DiagGammas;
S = cell(nT, numel(gams), numel(pths) + 1);
for k = 1:nT
    z = trialNoise(cfg, D, 'S3a', k);                     % v6: fixed seed for this stage, speed and trial
    Nw = noiseWindow(cfg, D.v0, sel.st.cpi, D.fs);
    for g = 1:numel(gams)
        x = addNoise(D.clean, Nw, gams(g), z);
        S{k,g,1} = runStep(P, x, D, sel.st.cpi, gams(g), [], []);                 % original coupled Pth
        for p = 1:numel(pths)
            S{k,g,p+1} = runStep(P, x, D, sel.st.cpi, gams(g), pths(p), []);      % forced Pth
        end
    end
end
for g = 1:numel(gams)
    for p = 1:numel(pths) + 1
        for param = {'StepTime', 'StrideDuration'}
            r = rowST(cfg, D, squeeze(S(:,g,p)), sel.st.cpi, gams(g), 'D_Pth', param{1});
            if p == 1
                r.Stage = 'D_PthCoupled';
                r.Pth_dB = max(10 - gams(g), -40);
            else
                r.Pth_dB = pths(p-1);
            end
            rows = addRow(rows, r);
        end
    end
end

scs = cfg.CPeakList;
Sc = cell(nT, numel(gams), numel(scs));
for k = 1:nT
    z = trialNoise(cfg, D, 'S3e', k);                     % v6: fixed seed for this stage, speed and trial
    Nw = noiseWindow(cfg, D.v0, sel.st.cpi, D.fs);
    for g = 1:numel(gams)
        x = addNoise(D.clean, Nw, gams(g), z);
        for q = 1:numel(scs)
            Sc{k,g,q} = runStep(P, x, D, sel.st.cpi, gams(g), [], [], scs(q));
        end
    end
end
for g = 1:numel(gams)
    for q = 1:numel(scs)
        for param = {'StepTime', 'StrideDuration'}
            r = rowST(cfg, D, squeeze(Sc(:,g,q)), sel.st.cpi, gams(g), 'D_CPeak', param{1});
            r.CPeakScale = scs(q);
            rows = addRow(rows, r);
        end
    end
end

hops = cfg.HopList;
cpiSD = sel.sd.cpi;
wsA = sel.sl.cs;
ctA = sel.sl.ct;
Wd = cell(nT, 1);
Sst = cell(nT, numel(hops));
Ssd = cell(nT, numel(hops));
Ssl = cell(nT, numel(hops));
for k = 1:nT
    z = trialNoise(cfg, D, 'S3b', k);                     % v6: fixed seed for this stage, speed and trial
    x = addNoise(D.clean, noiseWindow(cfg, D.v0, wsA, D.fs), cfg.HopSNR, z);
    Wd{k} = runSpeed(P, x, D, wsA);
    for h = 1:numel(hops)
        xs = addNoise(D.clean, noiseWindow(cfg, D.v0, sel.st.cpi, D.fs), cfg.HopSNR, z);
        Sst{k,h} = runStep(P, xs, D, sel.st.cpi, cfg.HopSNR, [], hops(h));
        xd = addNoise(D.clean, noiseWindow(cfg, D.v0, cpiSD, D.fs), cfg.HopSNR, z);
        Ssd{k,h} = runStep(P, xd, D, cpiSD, cfg.HopSNR, [], hops(h));
        xl = addNoise(D.clean, noiseWindow(cfg, D.v0, ctA, D.fs), cfg.HopSNR, z);
        Ssl{k,h} = runStep(P, xl, D, ctA, cfg.HopSNR, [], hops(h));
    end
end
for h = 1:numel(hops)
    r1 = rowST(cfg, D, Sst(:,h), sel.st.cpi, cfg.HopSNR, 'D_Hop', 'StepTime');
    r2 = rowST(cfg, D, Ssd(:,h), cpiSD, cfg.HopSNR, 'D_Hop', 'StrideDuration');
    r3 = rowSL(cfg, D, Wd, Ssl(:,h), wsA, ctA, cfg.HopSNR, 'D_Hop');
    r1.HopFraction = hops(h);
    r2.HopFraction = hops(h);
    r3.HopFraction = hops(h);
    rows = addRow(rows, r1);
    rows = addRow(rows, r2);
    rows = addRow(rows, r3);
end

if D.v0 < 1.0
    modes = {'extractorN', 'simulatorN'};
    Wm = cell(nT, 2);
    Sm = cell(nT, 2);
    for k = 1:nT
        z = trialNoise(cfg, D, 'S3c', k);                     % v6: fixed seed for this stage, speed and trial
        for m = 1:2
            cm = cfg;
            cm.NoiseMode = modes{m};
            xw = addNoise(D.clean, noiseWindow(cm, D.v0, sel.ws.cpi, D.fs), cfg.CPIfixedSNR, z);
            Wm{k,m} = runSpeed(P, xw, D, sel.ws.cpi);
            xs = addNoise(D.clean, noiseWindow(cm, D.v0, sel.st.cpi, D.fs), cfg.CPIfixedSNR, z);
            Sm{k,m} = runStep(P, xs, D, sel.st.cpi, cfg.CPIfixedSNR, [], []);
        end
    end
    for m = 1:2
        cm = cfg;
        cm.NoiseMode = modes{m};
        rows = addRow(rows, rowWS(cm, D, Wm(:,m), sel.ws.cpi, cfg.CPIfixedSNR, 'D_Noise'));
        rows = addRow(rows, rowST(cm, D, Sm(:,m), sel.st.cpi, cfg.CPIfixedSNR, 'D_Noise', 'StepTime'));
    end
end
end

function [rows, rowsOld] = stage4(P, cfg, D, sel, outDir)
logMsg(outDir, 'Stage 4: validation. SNR: speed %g, step %g, stride dur %g, stride len %g dB', ...
    sel.ws.snr, sel.st.snr, sel.sd.snr, sel.sl.snr);
nT = cfg.nValidate;
W1 = cell(nT, 1);
W2 = cell(nT, 1);
S1 = cell(nT, 1);
S2 = cell(nT, 1);
S3 = cell(nT, 1);
for k = 1:nT
    z = trialNoise(cfg, D, 'S4', k);                     % v6: fixed seed for this stage, speed and trial
    W1{k} = runSpeed(P, addNoise(D.clean, noiseWindow(cfg, D.v0, sel.ws.cpi, D.fs), sel.ws.snr, z), D, sel.ws.cpi);
    S1{k} = runStep(P, addNoise(D.clean, noiseWindow(cfg, D.v0, sel.st.cpi, D.fs), sel.st.snr, z), D, sel.st.cpi, sel.st.snr, [], []);
    S2{k} = runStep(P, addNoise(D.clean, noiseWindow(cfg, D.v0, sel.sd.cpi, D.fs), sel.sd.snr, z), D, sel.sd.cpi, sel.sd.snr, [], []);
    W2{k} = runSpeed(P, addNoise(D.clean, noiseWindow(cfg, D.v0, sel.sl.cs, D.fs), sel.sl.snr, z), D, sel.sl.cs);
    S3{k} = runStep(P, addNoise(D.clean, noiseWindow(cfg, D.v0, sel.sl.ct, D.fs), sel.sl.snr, z), D, sel.sl.ct, sel.sl.snr, [], []);
end
rows = [];
rows = addRow(rows, rowWS(cfg, D, W1, sel.ws.cpi, sel.ws.snr, 'S4'));
rows = addRow(rows, rowST(cfg, D, S1, sel.st.cpi, sel.st.snr, 'S4', 'StepTime'));
rows = addRow(rows, rowST(cfg, D, S2, sel.sd.cpi, sel.sd.snr, 'S4', 'StrideDuration'));
rows = addRow(rows, rowSL(cfg, D, W2, S3, sel.sl.cs, sel.sl.ct, sel.sl.snr, 'S4'));
rowsOld = [];                                                % the same noise scored with the v3 rules, for the before and after table
if cfg.CompareOldRules
    rowsOld = addRow(rowsOld, rowST(cfg, D, S1, sel.st.cpi, sel.st.snr, 'S4old', 'StepTime', true));
    rowsOld = addRow(rowsOld, rowST(cfg, D, S2, sel.sd.cpi, sel.sd.snr, 'S4old', 'StrideDuration', true));
end
end

function assembleResults(cfg, outDir)
S1 = []; S2 = []; S3 = []; S4 = []; S4old = [];
for v0 = cfg.Speeds
    S1 = [S1, loadRows(outDir, 'S1', v0, 'rows1')];                      %#ok<AGROW>
    S2 = [S2, loadRows(outDir, 'S2', v0, 'rows2')];                      %#ok<AGROW>
    S3 = [S3, loadRows(outDir, 'S3', v0, 'rows3')];                      %#ok<AGROW>
    S4 = [S4, loadRows(outDir, 'S4', v0, 'rows4')];                      %#ok<AGROW>
    S4old = [S4old, loadRows(outDir, 'S4', v0, 'rows4old')];             %#ok<AGROW>
end
writeRowsCSV(fullfile(outDir, 'Stage1_CPI_sweep.csv'), S1);
writeRowsCSV(fullfile(outDir, 'Stage2_SNR_sweep.csv'), S2);
writeRowsCSV(fullfile(outDir, 'Stage3_diagnostics.csv'), S3);
writeRowsCSV(fullfile(outDir, 'Stage4_validation.csv'), S4);
writeRowsCSV(fullfile(outDir, 'Stage4_validation_v3rules.csv'), S4old);
if cfg.MakePlots
    try
        makePlots(cfg, outDir, S1, S2);
    catch ME
        logMsg(outDir, 'Plots skipped: %s', ME.message);
    end
end
writeReport(cfg, outDir, S1, S2, S3, S4, S4old);
try
    writeTcComparison(cfg, outDir, S1, S4);                  % v6: comparison with the previous run, report and figure
catch ME
    logMsg(outDir, 'Tc comparison failed: %s', ME.message);
end
end

function rows = loadRows(outDir, tag, v0, name)
f = fullfile(outDir, sprintf('%s_v0_%.2f.mat', tag, v0));
rows = [];
if isfile(f)
    s = load(f);
    if isfield(s, name)
        rows = s.(name);
    end
end
if isempty(rows)
    rows = repmat(blankRow(), 1, 0);
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
            parts{k} = sprintf('%.6g', val);
        end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
fclose(fid);
end

function sub = subRows(rows, param, speed)
if isempty(rows)
    sub = rows;
    return;
end
keep = strcmp({rows.Parameter}, param);
if ~isempty(speed)
    keep = keep & abs([rows.Speed_mps] - speed) < 1e-9;
end
sub = rows(keep);
end

function makePlots(cfg, outDir, S1, S2)
params = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideLength'};
niceName = {'Walking speed', 'Step time', 'Stride duration', 'Stride length'};
speeds = cfg.Speeds;
yCap = 30;
for p = 1:4
    fig = newLightFigure();
    ax = axes('Parent', fig);
    hold(ax, 'on');
    clipped = false;
    for v = 1:numel(speeds)
        r = subRows(S1, params{p}, speeds(v));
        if isempty(r)
            continue;
        end
        if strcmp(params{p}, 'StrideLength')
            r = r([r.CPI_speed_req] == [r.CPI_step_req]);        % same CPI for both, to give a single curve
        end
        if strcmp(params{p}, 'WalkingSpeed') || strcmp(params{p}, 'StrideLength')
            x = [r.CPI_speed_act];
        else
            x = [r.CPI_step_act];
        end
        clipped = clipped | any([r.TrialMean_pct] > yCap);
        errorbar(ax, x, [r.TrialMean_pct], [r.TrialStd_pct], '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', sprintf('%.1f m/s', speeds(v)));
    end
    set(ax, 'XScale', 'log');
    hT = plot(ax, xlim(ax), [cfg.Target_pct cfg.Target_pct], 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Target (%g %%)', cfg.Target_pct));
    ylim(ax, [0 yCap]);
    finishLightAxes(ax, 'Actual CPI (s, log scale)', 'RMSE (% of the true value)', ...
        {sprintf('%s: RMSE vs CPI', niceName{p}), ...
         sprintf('SNR fixed at %g dB, %d trials per point, mean \\pm std, one line per walking speed', cfg.CPIfixedSNR, cfg.nTrials)}, ...
        clipped, yCap, 'Walking speed');
    saveLightFigure(fig, fullfile(outDir, sprintf('Plot_CPI_%s', params{p})));
    close(fig);

    fig = newLightFigure();
    ax = axes('Parent', fig);
    hold(ax, 'on');
    clipped = false;
    for v = 1:numel(speeds)
        r = subRows(S2, params{p}, speeds(v));
        if isempty(r)
            continue;
        end
        clipped = clipped | any([r.TrialMean_pct] > yCap);
        errorbar(ax, [r.SNR_dB], [r.TrialMean_pct], [r.TrialStd_pct], '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', sprintf('%.1f m/s', speeds(v)));
    end
    plot(ax, xlim(ax), [cfg.Target_pct cfg.Target_pct], 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Target (%g %%)', cfg.Target_pct));
    ylim(ax, [0 yCap]);
    finishLightAxes(ax, 'SNR (dB)', 'RMSE (% of the true value)', ...
        {sprintf('%s: RMSE vs SNR', niceName{p}), ...
         sprintf('Each speed at its own best CPI from the CPI sweep, %d trials per point, mean \\pm std', cfg.nTrials)}, ...
        clipped, yCap, 'Walking speed');
    saveLightFigure(fig, fullfile(outDir, sprintf('Plot_SNR_%s', params{p})));
    close(fig);
end

for sweepType = 1:2
    fig = newLightFigure();
    axTop = subplot(2, 1, 1, 'Parent', fig);
    hold(axTop, 'on');
    axBot = subplot(2, 1, 2, 'Parent', fig);
    hold(axBot, 'on');
    for v = 1:numel(speeds)
        if sweepType == 1
            r = subRows(S1, 'StepTime', speeds(v));
            x = [r.CPI_step_act];
        else
            r = subRows(S2, 'StepTime', speeds(v));
            x = [r.SNR_dB];
        end
        if isempty(r)
            continue;
        end
        nm = sprintf('%.1f m/s', speeds(v));
        plot(axTop, x, [r.MissRate_pct] + [r.ExtraRate_pct], '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', nm);
        plot(axBot, x, [r.PerfectTrials_pct], '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', nm);
    end
    plot(axTop, xlim(axTop), [cfg.MaxMissExtra_pct cfg.MaxMissExtra_pct], 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Limit (%g %%)', cfg.MaxMissExtra_pct));
    if sweepType == 1
        set(axTop, 'XScale', 'log');
        set(axBot, 'XScale', 'log');
        xText = 'Actual CPI (s, log scale)';
        sweepText = sprintf('CPI sweep at %g dB SNR', cfg.CPIfixedSNR);
        tag = 'Plot_CPI_StepRecognition';
    else
        xText = 'SNR (dB)';
        sweepText = 'SNR sweep at each speed''s best CPI';
        tag = 'Plot_SNR_StepRecognition';
    end
    ylim(axTop, [0 50]);
    ylim(axBot, [0 105]);
    finishLightAxes(axTop, '', 'Misses + extras (% of true steps)', ...
        {sprintf('Step recognition: %s', sweepText), sprintf('%d trials per point. Top: missed plus extra steps, lower is better. Bottom: trials with exactly the right steps, higher is better.', cfg.nTrials)}, ...
        false, 50, 'Walking speed');
    finishLightAxes(axBot, xText, 'Exactly right trials (% of trials)', {}, false, 105, '');
    try                                                              % line the lower plot up with the upper one, which has the key beside it
        pT = get(axTop, 'Position');
        pB = get(axBot, 'Position');
        pB([1 3]) = pT([1 3]);
        set(axBot, 'Position', pB);
    catch
    end
    saveLightFigure(fig, fullfile(outDir, tag));
    close(fig);
end
end

function fig = newLightFigure()
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 720]);
end

function finishLightAxes(ax, xText, yText, titleLines, clipped, yCap, legendTitle)
grid(ax, 'on');
set(ax, 'FontSize', 12, 'LineWidth', 1, 'Box', 'on');
if ~isempty(xText)
    xlabel(ax, xText, 'FontSize', 13);
end
ylabel(ax, yText, 'FontSize', 13);
if ~isempty(titleLines)
    title(ax, titleLines, 'FontSize', 13, 'FontWeight', 'bold');
end
if ~isempty(legendTitle)
    lg = legend(ax, 'show', 'Location', 'eastoutside', 'FontSize', 11);
    try, title(lg, legendTitle); catch, end                          % key heading
end
if clipped
    text(ax, 0.99, 0.97, sprintf('Values above %g %% are clipped. They are in the CSV files.', yCap), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9, 'Color', [0.35 0.35 0.35]);
end
styleLightFigure(ancestor(ax, 'figure'));
end

function styleLightFigure(fig)
try, theme(fig, 'light'); catch, end                                 % R2025a and newer only, ignored on older versions
try
    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');                 % white figure, keep the colours exactly as drawn
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
        set(txt(k), 'Color', [0.1 0.1 0.1]);                         % titles, labels and notes all dark
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
    saveas(fig, [basePath '.png']);                                  % older MATLAB without exportgraphics
end
end

function writeReport(cfg, outDir, S1, S2, S3, S4, S4old)
L = {};
L{end+1} = 'DIAGNOSIS REPORT (written automatically)';
L{end+1} = sprintf('Run folder: %s', outDir);
L{end+1} = sprintf('Simulated data only. R0 = %g m, amplitude mode %s, noise mode %s, %d trials per condition.', ...
    cfg.R0_m, cfg.AmplitudeMode, cfg.NoiseMode, cfg.nTrials);
L{end+1} = sprintf('CPI sweep SNR = %g dB. Target = RMSE below %g %%.', cfg.CPIfixedSNR, cfg.Target_pct);
L{end+1} = 'RMSE % = 100 * RMSE / mean true value (same definitions as ChenAmp / v6 Section 6).';
L{end+1} = '';

L{end+1} = '1. CHOSEN SETTINGS AND RESULTS (Stage 4 uses fresh noise, so it is not the data the settings were chosen on)';
L{end+1} = '  miss % = true steps with no detection. extra % = detections with no true step (inside the walk only), both as a percentage of the true steps. exact trials % = trials with zero misses AND zero extras. NaN = not a step based parameter.';
params = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideLength'};
for p = 1:4
    L{end+1} = sprintf('%s', params{p});
    L{end+1} = '  speed | CPI speed | CPI step | SNR | selection RMSE% | validated RMSE% | bias % | miss % | extra % | exact trials % | below target';
    vals = [];
    for v0 = cfg.Speeds
        r4 = subRows(S4, params{p}, v0);
        r2 = subRows(S2, params{p}, v0);
        if isempty(r4)
            L{end+1} = sprintf('  %.1f | not finished', v0);
            continue;
        end
        selRM = NaN;
        if ~isempty(r2)
            hit = r2(abs([r2.SNR_dB] - r4(1).SNR_dB) < 1e-9);
            if ~isempty(hit)
                selRM = hit(1).RMSE_pct;
            end
        end
        L{end+1} = sprintf('  %.1f | %s | %s | %g | %.2f | %.2f | %+.2f | %.1f | %.1f | %.0f | %s', v0, ...
            numText(r4(1).CPI_speed_act), numText(r4(1).CPI_step_act), r4(1).SNR_dB, selRM, r4(1).RMSE_pct, ...
            r4(1).MeanErr_pct, r4(1).MissRate_pct, r4(1).ExtraRate_pct, r4(1).PerfectTrials_pct, yesNo(r4(1).RMSE_pct < cfg.Target_pct));
        vals(end+1) = r4(1).RMSE_pct;                                    %#ok<AGROW>
    end
    if ~isempty(vals)
        L{end+1} = sprintf('  Speeds below target: %d of %d. Mean validated RMSE %.2f %%, worst %.2f %%.', ...
            nnz(vals < cfg.Target_pct), numel(vals), meanFinite(vals), max(vals));
    end
    L{end+1} = '';
end

L{end+1} = '2. WHY THE BEST CPI CHANGES WITH SPEED (walking speed)';
L{end+1} = '  The torso Doppler frequency is f = 2*v*fc/c. The FFT resolution is 1/CPI. A long window is needed when f is small.';
L{end+1} = '  speed | torso Doppler Hz | best CPI s | resolution as % of f | RMSE% at best | RMSE% at 0.186 s';
bestC = [];
spd = [];
for v0 = cfg.Speeds
    r1 = subRows(S1, 'WalkingSpeed', v0);
    if isempty(r1)
        continue;
    end
    [~, i] = min([r1.RMSE_pct]);
    ref = r1(abs([r1.CPI_speed_act] - 8192/44100) < 1e-3);
    refText = 'n/a';
    if ~isempty(ref)
        refText = sprintf('%.2f', ref(1).RMSE_pct);
    end
    f = 2*v0*24e9/3e8;
    L{end+1} = sprintf('  %.1f | %.1f | %.3f | %.1f | %.2f | %s', v0, f, r1(i).CPI_speed_act, ...
        100*(1/r1(i).CPI_speed_act)/f, r1(i).RMSE_pct, refText);
    bestC(end+1) = r1(i).CPI_speed_act;                                  %#ok<AGROW>
    spd(end+1) = v0;                                                     %#ok<AGROW>
end
if numel(bestC) >= 4
    lowM = meanFinite(bestC(spd <= 0.6));
    highM = meanFinite(bestC(spd >= 1.0));
    L{end+1} = sprintf('  Mean best CPI: %.3f s at 0.6 m/s and below, %.3f s at 1.0 m/s and above.', lowM, highM);
    if lowM > highM
        L{end+1} = '  This shows that slower walking needs a longer window, which is why one shared CPI chosen at 1.5 m/s fails at 0.4 to 0.6 m/s.';
    else
        L{end+1} = '  The best CPI did not grow at low speed in this run. The shared-CPI explanation is therefore not supported here.';
    end
end
L{end+1} = '';

L{end+1} = '3. STEP TIME STAMP OFFSET (step recognition)';
L{end+1} = '  Step times come from the spectrogram frame START, not the frame centre. The offset below is detected minus true time.';
L{end+1} = '  Miss/extra are shown with the original stamps and with both moved to the window centre.';
L{end+1} = '  CPI s | mean offset ms | miss % | miss centred % | extra % | extra centred %  (averaged over speeds, Stage 1)';
stRows = subRows(S1, 'StepTime', []);
if ~isempty(stRows)
    cpiVals = unique([stRows.CPI_step_act]);
    for cv = cpiVals
        sub = stRows(abs([stRows.CPI_step_act] - cv) < 1e-9);
        L{end+1} = sprintf('  %.3f | %.0f | %.1f | %.1f | %.1f | %.1f', cv, meanFinite([sub.StepOffset_ms]), ...
            meanFinite([sub.MissRate_pct]), meanFinite([sub.MissCentred_pct]), ...
            meanFinite([sub.ExtraRate_pct]), meanFinite([sub.ExtraCentred_pct]));
    end
end
L{end+1} = '  A growing offset with CPI, and centred miss/extra below the original, means the miss rate partly measures the time stamp convention.';
L{end+1} = '  Step intervals and stride durations are differences of two times, so the offset cancels there.';
L{end+1} = '';

L{end+1} = '4. SNR CURVES THAT GET WORSE AT HIGH SNR, AND THE Pth RULE';
L{end+1} = '  The original rule Pth = max(10 - SNR, -40) dB ties the leg threshold to the SNR and stops moving at 50 dB.';
flagged = 0;
for p = 1:4
    for v0 = cfg.Speeds
        r = subRows(S2, params{p}, v0);
        if numel(r) < 3
            continue;
        end
        rm = [r.RMSE_pct];
        sn = [r.SNR_dB];
        worse = false;
        for a = 1:numel(rm)
            for b = a+1:numel(rm)
                if isfinite(rm(a)) && isfinite(rm(b)) && rm(b) > rm(a) + 1.0
                    worse = true;
                end
            end
        end
        if worse
            flagged = flagged + 1;
            L{end+1} = sprintf('  %s at %.1f m/s: RMSE is not steadily falling with SNR (%s)', params{p}, v0, ...
                sprintf('%g dB:%.1f ', [sn; rm]));
        end
    end
end
if flagged == 0
    L{end+1} = '  No curve got worse by more than 1 percentage point at a higher SNR.';
end
L{end+1} = '  Forced-Pth test (step time, at the Stage 1 best step CPI): coupled rule vs forced Pth, miss + extra % and step RMSE %';
for v0 = cfg.Speeds
    for gm = cfg.DiagGammas
        if isempty(S3)
            continue;
        end
        c = S3(strcmp({S3.Stage}, 'D_PthCoupled') & strcmp({S3.Parameter}, 'StepTime') & abs([S3.Speed_mps] - v0) < 1e-9 & [S3.SNR_dB] == gm);
        f = S3(strcmp({S3.Stage}, 'D_Pth') & strcmp({S3.Parameter}, 'StepTime') & abs([S3.Speed_mps] - v0) < 1e-9 & [S3.SNR_dB] == gm);
        if isempty(c) || isempty(f)
            continue;
        end
        [~, i] = min([f.RMSE_pct] + [f.MissRate_pct] + [f.ExtraRate_pct]);
        L{end+1} = sprintf('  %.1f m/s at %g dB: coupled Pth %g dB -> RMSE %.2f%%, miss+extra %.1f%% | best forced Pth %g dB -> RMSE %.2f%%, miss+extra %.1f%%', ...
            v0, gm, c(1).Pth_dB, c(1).RMSE_pct, c(1).MissRate_pct + c(1).ExtraRate_pct, ...
            f(i).Pth_dB, f(i).RMSE_pct, f(i).MissRate_pct + f(i).ExtraRate_pct);
    end
end
L{end+1} = '  If a forced Pth beats the coupled rule at the same SNR, the Pth rule (not the radar noise) is limiting the result.';
L{end+1} = '';
L{end+1} = '  Peak-rule test (Stage 3e): the step detector keeps a frame as a step only if the leg frequency falls by at least c (= Doppler resolution divided by the leg-frequency cut-off freq_high_Yardibi) to the next frame.';
L{end+1} = '  At fast speeds the top of each swing can be flat, so the fall is smaller than c and the step is skipped at every noise draw. Same noise, cPeak x each scale:';
for v0 = cfg.Speeds
    for gm = cfg.DiagGammas
        if isempty(S3)
            continue;
        end
        r = S3(strcmp({S3.Stage}, 'D_CPeak') & strcmp({S3.Parameter}, 'StepTime') & abs([S3.Speed_mps] - v0) < 1e-9 & [S3.SNR_dB] == gm);
        if isempty(r)
            continue;
        end
        [~, ord] = sort([r.CPeakScale], 'descend');
        r = r(ord);
        txtq = '';
        for q = 1:numel(r)
            txtq = [txtq sprintf('x%g: miss+extra %.1f%%, RMSE %.2f%%   ', r(q).CPeakScale, r(q).MissRate_pct + r(q).ExtraRate_pct, r(q).RMSE_pct)];   %#ok<AGROW>
        end
        L{end+1} = sprintf('  %.1f m/s at %g dB: %s', v0, gm, txtq);
    end
end
L{end+1} = '  If a smaller scale removes the misses without adding extras, the peak rule (not noise, and not the Pth rule) is what limits step recognition at that speed.';
L{end+1} = '  The MATLAB run of Master_SimSweep_v2 on 8 October 2026 showed the same pattern (for example 2.0 m/s at 60 dB: miss+extra 11.1 % at x1 and 0.0 % at x0.5).';
L{end+1} = '';

L{end+1} = '5. STEP SPECTROGRAM HOP SIZE (time quantisation)';
L{end+1} = '  The step script uses 75 % overlap, so step times can only land on a grid of N/4 samples.';
L{end+1} = '  The peak-drop rule (cPeak) is scaled with the hop in the working copy. Without that scaling a smaller hop finds almost no peaks.';
L{end+1} = sprintf('  Tested at %g dB so noise is small. RMSE %% for each hop (fraction of the window):', cfg.HopSNR);
for params2 = {'StepTime', 'StrideDuration', 'StrideLength'}
    for v0 = cfg.Speeds
        if isempty(S3)
            continue;
        end
        r = S3(strcmp({S3.Stage}, 'D_Hop') & strcmp({S3.Parameter}, params2{1}) & abs([S3.Speed_mps] - v0) < 1e-9);
        if isempty(r)
            continue;
        end
        L{end+1} = sprintf('  %s %.1f m/s: %s', params2{1}, v0, sprintf('hop %.4g: %.2f%%   ', [[r.HopFraction]; [r.RMSE_pct]]));
    end
end
L{end+1} = '  A small gain from a smaller hop points to time quantisation. A very large change (tens of percent to a few percent) points to the peak rule being sensitive to the frame spacing, with extra or missing peaks.';
L{end+1} = '  Check Stage 3 miss and extra columns for the same rows to tell the two apart.';
L{end+1} = '';

L{end+1} = '6. NOISE WINDOW BELOW 1 m/s';
L{end+1} = '  The simulator forces a 4096-sample window below 1 m/s, so its SNR is set on that window, not on the extractor window.';
L{end+1} = '  The effective SNR at the extractor is therefore SNR + 10*log10(N_extractor/4096).';
if ~isempty(S1)
    for cv = unique([S1.CPI_step_act])
        if ~isfinite(cv)
            continue;
        end
        L{end+1} = sprintf('  CPI %.3f s: effective SNR shift %+.1f dB', cv, 10*log10((cv*44100)/4096));
    end
end
for v0 = cfg.Speeds(cfg.Speeds < 1)
    for params2 = {'WalkingSpeed', 'StepTime'}
        if isempty(S3)
            continue;
        end
        r = S3(strcmp({S3.Stage}, 'D_Noise') & strcmp({S3.Parameter}, params2{1}) & abs([S3.Speed_mps] - v0) < 1e-9);
        if numel(r) == 2
            L{end+1} = sprintf('  %s %.1f m/s: %s %.2f%%, %s %.2f%%', params2{1}, v0, r(1).NoiseMode, r(1).RMSE_pct, r(2).NoiseMode, r(2).RMSE_pct);
        end
    end
end
L{end+1} = '';

L{end+1} = '7. STRIDE LENGTH ERROR BUDGET';
L{end+1} = '  Stride length is the integrated speed over one detected stride, so its error should be close to the speed and stride duration errors combined.';
L{end+1} = '  speed | walking speed % | stride duration % | sqrt of sum of squares | stride length %';
for v0 = cfg.Speeds
    a = subRows(S4, 'WalkingSpeed', v0);
    b = subRows(S4, 'StrideDuration', v0);
    c = subRows(S4, 'StrideLength', v0);
    if isempty(a) || isempty(b) || isempty(c)
        continue;
    end
    L{end+1} = sprintf('  %.1f | %.2f | %.2f | %.2f | %.2f', v0, a(1).RMSE_pct, b(1).RMSE_pct, ...
        sqrt(a(1).RMSE_pct^2 + b(1).RMSE_pct^2), c(1).RMSE_pct);
end
L{end+1} = '  If stride length is much larger than the combined value, the extra error comes from time axis mismatch or unmatched strides, not from speed or timing alone.';
L{end+1} = '  Stage 1 also tried different CPIs for the speed and step parts. A split pair that beats the same-CPI pair is shown in the Stage 1 CSV.';
L{end+1} = '';

L{end+1} = '8. CONDITIONS STILL ABOVE THE TARGET';
anyMiss = false;
for p = 1:4
    for v0 = cfg.Speeds
        r4 = subRows(S4, params{p}, v0);
        if ~isempty(r4) && r4(1).RMSE_pct >= cfg.Target_pct
            anyMiss = true;
            L{end+1} = sprintf('  %s at %.1f m/s: %.2f %% (CPI speed %s, CPI step %s, SNR %g dB, failed trials %d)', ...
                params{p}, v0, r4(1).RMSE_pct, numText(r4(1).CPI_speed_act), numText(r4(1).CPI_step_act), r4(1).SNR_dB, r4(1).FailedTrials);
        end
    end
end
if ~anyMiss
    L{end+1} = '  None.';
end
L{end+1} = '';

L{end+1} = '9. BEFORE AND AFTER (the same noise draws scored two ways)';
L{end+1} = '  before = v3 rules: v4b step list, step times on the frame grid, every detection counted in the stride list.';
L{end+1} = '  after  = v4 rules: skipped step put back, step times placed between frames, stride list limited to the steps inside the walk.';
L{end+1} = '  parameter | speed | before RMSE % | after RMSE % | before miss/extra/exact % | after miss/extra/exact %';
for p = 2:3
    for v0 = cfg.Speeds
        rb = subRows(S4old, params{p}, v0);
        ra = subRows(S4, params{p}, v0);
        if isempty(rb) || isempty(ra)
            continue;
        end
        L{end+1} = sprintf('  %s | %.1f | %.2f | %.2f | %.1f / %.1f / %.0f | %.1f / %.1f / %.0f', params{p}, v0, rb(1).RMSE_pct, ra(1).RMSE_pct, ...
            rb(1).MissRate_pct, rb(1).ExtraRate_pct, rb(1).PerfectTrials_pct, ra(1).MissRate_pct, ra(1).ExtraRate_pct, ra(1).PerfectTrials_pct);
    end
end
L{end+1} = '  Walking speed and stride length are not changed by these rules, so they are not repeated here.';
L{end+1} = '';
L{end+1} = 'NOT CHECKED BY THIS SCRIPT';
L{end+1} = '  Real data. These are simulator results only.';
L{end+1} = '  R0 other than the value above. The v6 Section 7 runs used 50 m and these use the R0 printed at the top.';
L{end+1} = '  Per-speed settings need the speed to be known. In use, a coarse speed estimate would have to pick the CPI.';
fid = fopen(fullfile(outDir, 'Diagnosis_Report.txt'), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
end

function s = numText(v)
if isfinite(v)
    s = sprintf('%.3f', v);
else
    s = '-';
end
end

function s = yesNo(tf)
if tf
    s = 'yes';
else
    s = 'no';
end
end

function m_o = runStep(m_P, m_x, m_D, m_cpi, m_snr, m_pth, m_hop, m_cps)
global TcCfg                                                 % Tc method and limits, set once in the main function
if nargin < 8
    m_cps = [];                                              % no peak rule scale, same default as v5
end
tcEst = NaN;                                                 % the Tc estimate for this trial
firstPass = [];                                              % result of pass 1 (two-pass method only)
switch TcCfg.Method
    case 'twopass'
        firstPass = runStepPass(m_P, m_x, m_D, m_cpi, m_snr, m_pth, m_hop, m_cps, TcCfg.Nominal);   % pass 1 with the nominal Tc (1.2 s)
        tcEst = tcTwoPass(firstPass, m_D);                   % twice the median step interval of pass 1, steps inside the walk only
    case 'torso'
        tcEst = tcFromTorso(m_P, m_x, m_D);                  % twice the mean torso bob interval, as in the real data scripts
    case 'true'
        tcEst = m_D.Tc;                                      % check mode only: gives the extractor the true Tc, as v5 did
    otherwise
        error('Unknown TcMethod "%s". Use twopass, torso or true.', TcCfg.Method);
end
fallback = 0;                                                % 1 if the estimate could not be used
if ~strcmp(TcCfg.Method, 'true') && ~(isfinite(tcEst) && tcEst > TcCfg.Range(1) && tcEst < TcCfg.Range(2))
    tcEst = TcCfg.Nominal;                                   % too few steps (or a silly value): keep the nominal 1.2 s
    fallback = 1;                                            % counted as a fallback
end
if ~isempty(firstPass) && fallback == 1
    m_o = firstPass;                                         % pass 2 would use the same Tc as pass 1, so pass 1 is the output
else
    m_o = runStepPass(m_P, m_x, m_D, m_cpi, m_snr, m_pth, m_hop, m_cps, tcEst);   % pass 2 (or the only pass) with Tc_est
end
m_o.TcEst = tcEst;                                           % remember what the extractor was given
m_o.TcFallback = fallback;                                   % and whether it was the fallback
end

function T = tcTwoPass(o, D)
T = NaN;
if ~o.ok
    return;                                                  % pass 1 failed, so there is no estimate
end
edge = 0.1 + o.cpiAct/2;                                     % same walk interior as the scorer
inWalk = o.key > D.buffer + edge & o.key < D.buffer + D.simT - edge;   % steps inside the walk only
t = o.det(inWalk);                                           % their times
if numel(t) >= 3
    T = 2*median(diff(t));                                   % a full gait cycle is two steps
end
end

function T = tcFromTorso(m_P, m_x, m_D)
global TcCfg
T = NaN;
try
    w = runTorsoTrack(m_P, m_x, m_D, TcCfg.TorsoCPI);        % torso track from the walking speed script (not changed)
    if ~w.ok
        return;
    end
    f = w.torsoHz(:).';                                      % torso frequency track in Hz
    tv = w.tv(:).';                                          % its time axis
    if any(~isfinite(f))
        f = fillmissing(f, 'linear', 'EndValues', 'nearest');   % findpeaks cannot use gaps
    end
    dt = tv(2) - tv(1);                                      % frame spacing
    prom = 0.25*(max(f) - min(f));
    if prom <= 0
        return;
    end
    [~, idx] = findpeaks(f, 'MinPeakDistance', max(1, round(0.3/dt)), 'MinPeakProminence', prom);   % steps cannot be closer than 0.3 s (real data setting)
    t = tv(idx);                                             % times of the bobs
    edge = 0.1 + w.cpiAct/2;                                 % same walk interior as the scorer
    t = t(t > m_D.buffer + edge & t < m_D.buffer + m_D.simT - edge);   % bobs inside the walk only
    if numel(t) >= 3
        T = 2*mean(diff(t));                                 % mean step interval doubled, as in the real data scripts
    end
catch
    T = NaN;                                                 % any failure counts as no estimate
end
end

function m_o = runTorsoTrack(m_P, m_x, m_D, m_cpi)
m_o = struct('ok', false, 'torsoHz', [], 'tv', [], 'cpiAct', NaN, 'err', '');
try
    RdrRxSignal = m_x;                                       %#ok<NASGU>
    fs = m_D.fs;                                             %#ok<NASGU>
    fc = m_D.fc;                                             %#ok<NASGU>
    c = m_D.c;                                               %#ok<NASGU>
    v0 = m_D.v0;                                             %#ok<NASGU>
    CPI_s = m_cpi;                                           %#ok<NASGU>
    evalc('run(m_P.ws);');
    m_o.torsoHz = final.torsoFreq_refined_Hz;                % the torso track the real data scripts find their bobs in
    m_o.tv = final.time_axis;
    m_o.cpiAct = CPI_actual_s_WS;
    m_o.ok = numel(m_o.torsoHz) >= 3;
catch ME
    m_o.err = ME.message;
end
end

function z = trialNoise(cfg, D, tag, k)
tags = {'S1', 'S2', 'S3a', 'S3e', 'S3b', 'S3c', 'S4'};      % one code per place a noise draw is made
code = find(strcmp(tags, tag), 1);
seed = cfg.Seed + 1000000*code + 1000*round(D.v0*10) + k;   % unique for every stage, speed and trial
rng(seed, 'twister');                                        % fixed seed, Mersenne twister
z = randn(1, numel(D.clean)) + 1i*randn(1, numel(D.clean));  % the draw (same form as v5)
f = fullfile(cfg.ResumeFolder, 'Seeds_Used.csv');            % list of every seed used in this run
newFile = ~isfile(f);
fid = fopen(f, 'a');
if fid > 0
    if newFile
        fprintf(fid, 'Stage,Speed_mps,Trial,Seed\n');        % header on the first line
    end
    fprintf(fid, '%s,%.2f,%d,%d\n', tag, D.v0, k, seed);
    fclose(fid);
end
end

function row = addTcStats(row, S, D)
est = NaN(numel(S), 1);                                      % the Tc given to the extractor in each trial
fb = 0;                                                      % number of fallbacks
for k = 1:numel(S)
    if isfield(S{k}, 'TcEst') && isfinite(S{k}.TcEst)
        est(k) = S{k}.TcEst;
        fb = fb + S{k}.TcFallback;
    end
end
row.Tc_true_s = D.Tc;                                        % the true Tc (used by the scorer only)
row.Tc_fallbacks = fb;
good = est(isfinite(est));
if ~isempty(good)
    row.Tc_est_mean_s = mean(good);                          % average estimate over the trials
    row.Tc_err_pct = 100*(row.Tc_est_mean_s - D.Tc)/D.Tc;    % signed error of the average estimate
    row.Tc_worst_abs_err_pct = max(abs(100*(good - D.Tc)/D.Tc));   % worst single trial
end
end

function selfTestTc(P, cfg, outDir)
global TcCfg
txt = fileread(P.st);                                        % the patched working copy
nEst = numel(regexp(txt, '(?<![A-Za-z0-9_])Tc_est(?![A-Za-z0-9_])', 'match'));   % uses of Tc_est (must be exactly the three patched places)
ok1 = (nEst == 3);
D = getSpeedData(P, cfg, cfg.Speeds(1));
z = trialNoise(cfg, D, 'S4', 0);                             % trial 0 is only used here
x = addNoise(D.clean, noiseWindow(cfg, D.v0, 0.1, D.fs), 40, z);
tcTry = 1.5*D.Tc;                                            % a deliberately wrong Tc
o1 = runStepPass(P, x, D, 0.1, 40, [], [], [], tcTry);       % the extractor must follow it
ok2 = o1.ok && abs(o1.medWin - tcTry/8) < 1e-12 && abs(o1.smWin - tcTry/8) < 1e-12;   % both smoothing windows are Tc_est/8
o2 = runStep(P, x, D, 0.1, 40, [], []);                      % the full method
ok3 = o2.ok && isfinite(o2.TcEst) && abs(o2.medWin - o2.TcEst/8) < 1e-12 && abs(o2.smWin - o2.TcEst/8) < 1e-12;   % windows follow the estimate
ok4 = strcmp(TcCfg.Method, 'true') || o2.TcFallback == 1 || abs(o2.TcEst - D.Tc) > 1e-9;   % the estimate is not just a copy of the true Tc
logMsg(outDir, 'Tc self test: patched copy has %d uses of Tc_est, expected 3 (%s). Leftover uses of the true Tc were checked when the copy was made.', nEst, passText(ok1));
logMsg(outDir, 'Tc self test: extractor given Tc = %.3f s used windows of %.4f s and %.4f s, expected %.4f s (%s)', tcTry, o1.medWin, o1.smWin, tcTry/8, passText(ok2));
logMsg(outDir, 'Tc self test: method %s at %.1f m/s gave Tc_est = %.3f s (true %.3f s, fallback %d), windows %.4f s (%s, %s)', ...
    TcCfg.Method, D.v0, o2.TcEst, D.Tc, o2.TcFallback, o2.medWin, passText(ok3), passText(ok4));
if ~(ok1 && ok2 && ok3 && ok4)
    error('Tc self test failed. The patch did not take effect or Tc_est is not used. See Run_Log.txt and Patch_Log.txt in %s', outDir);
end
end

function logPatch(logFile, tag, oldTxt, newTxt)
fid = fopen(logFile, 'a');
if fid > 0
    fprintf(fid, '[%s]\n  before: %s\n  after : %s\n', tag, strtrim(oldTxt), strrep(strtrim(newTxt), sprintf('\n'), sprintf('\n          ')));
    fclose(fid);
end
end

function logPatchMatches(logFile, tag, txt, pattern)
m = regexp(txt, pattern, 'match');
for q = 1:numel(m)
    logPatch(logFile, tag, m{q}, '(line removed)');
end
end

function ln = lineContaining(txt, pat)
ln = regexp(txt, ['[^\n]*' regexptranslate('escape', pat) '[^\n]*'], 'match', 'once');
end

function T = readCsvSimple(file)
fid = fopen(file, 'r');
assert(fid > 0, 'Could not open %s', file);
hdr = strsplit(strtrim(fgetl(fid)), ',');
data = {};
while true
    ln = fgetl(fid);
    if ~ischar(ln)
        break;
    end
    ln = strtrim(ln);
    if isempty(ln)
        continue;
    end
    parts = strsplit(ln, ',');
    parts(end+1:numel(hdr)) = {'NaN'};                       % short lines padded
    data(end+1, :) = strrep(parts(1:numel(hdr)), '"', '');   %#ok<AGROW>
end
fclose(fid);
T.hdr = hdr;
T.data = data;
end

function v = csvNum(T, name)
i = find(strcmp(T.hdr, name), 1);
if isempty(i)
    v = NaN(size(T.data, 1), 1);
else
    v = str2double(T.data(:, i));
end
end

function v = csvStr(T, name)
i = find(strcmp(T.hdr, name), 1);
v = T.data(:, i);
end

function thr = noiseLimit(param, snrA, snrB)
switch param
    case 'WalkingSpeed'
        thr = 0.09;
    case 'StepTime'
        thr = 0.31;
    case 'StrideDuration'
        thr = 0.78;
    case 'StrideLength'
        s = min([snrA snrB]);                                % the lower SNR gives the larger noise
        s = min(max(s, 30), 60);
        thr = 0.28 + (0.05 - 0.28)*(s - 30)/30;
    otherwise
        thr = NaN;
end
end

function writeTcComparison(cfg, outDir, S1, S4)
prevFile = fullfile(cfg.PrevRun, 'Stage4_validation.csv');
if ~isfile(prevFile)
    logMsg(outDir, 'Comparison skipped: %s was not found. Set PrevRun to the folder of the previous run.', prevFile);
    return;
end
Tp = readCsvSimple(prevFile);
pParam = csvStr(Tp, 'Parameter');
pSpeed = csvNum(Tp, 'Speed_mps');
metrics = {'RMSE_pct', 'Accuracy_pct', 'MAPE_pct', 'MeanErr_pct', 'TrialStd_pct', 'MissRate_pct', 'ExtraRate_pct', 'PerfectTrials_pct'};
PV = struct();                                               % previous run, one numeric column per field
for q = 1:numel(metrics)
    PV.(metrics{q}) = csvNum(Tp, metrics{q});
end
PV.SNR_dB = csvNum(Tp, 'SNR_dB');
PV.CPI_speed_req = csvNum(Tp, 'CPI_speed_req');
PV.CPI_step_req = csvNum(Tp, 'CPI_step_req');
params = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideLength'};
cmp = [];
for p = 1:4
    for v0 = cfg.Speeds
        rn = subRows(S4, params{p}, v0);                     % this run
        ri = find(strcmp(pParam, params{p}) & abs(pSpeed - v0) < 1e-9, 1);   % previous run
        if isempty(rn) || isempty(ri)
            continue;
        end
        c = struct();
        c.Parameter = params{p};
        c.Speed_mps = v0;
        c.SNR_prev = PV.SNR_dB(ri);
        c.SNR_new = rn(1).SNR_dB;
        c.CPIspeed_prev = PV.CPI_speed_req(ri);
        c.CPIspeed_new = rn(1).CPI_speed_req;
        c.CPIstep_prev = PV.CPI_step_req(ri);
        c.CPIstep_new = rn(1).CPI_step_req;
        for q = 1:numel(metrics)
            a = PV.(metrics{q})(ri);                         % previous value
            b = rn(1).(metrics{q});                          % value in this run
            c.([metrics{q} '_prev']) = a;
            c.([metrics{q} '_new']) = b;
            c.([metrics{q} '_diff']) = b - a;                % this run minus previous run
        end
        c.Pass_prev = double(c.RMSE_pct_prev < cfg.Target_pct);
        c.Pass_new = double(c.RMSE_pct_new < cfg.Target_pct);
        c.Noise_limit_pp = noiseLimit(params{p}, c.SNR_prev, c.SNR_new);
        c.Beyond_noise = double(abs(c.RMSE_pct_diff) > c.Noise_limit_pp);
        if c.Pass_prev == 1 && c.Pass_new == 0
            c.Change = 'worse';                              % lost a pass
        elseif c.Pass_prev == 0 && c.Pass_new == 1
            c.Change = 'better';                             % gained a pass
        elseif c.Beyond_noise == 1 && c.RMSE_pct_diff > 0
            c.Change = 'worse';
        elseif c.Beyond_noise == 1
            c.Change = 'better';
        else
            c.Change = 'none';
        end
        if isempty(cmp)
            cmp = c;
        else
            cmp(end+1) = c;                                  %#ok<AGROW>
        end
    end
end
if isempty(cmp)
    logMsg(outDir, 'Comparison skipped: no speed is in both runs.');
    return;
end
writeRowsCSV(fullfile(outDir, 'TcEstimate_Comparison.csv'), cmp);   % same writer as the other CSV files (text fields in quotes)
try
    plotTcComparison(cfg, outDir, cmp);
catch ME
    logMsg(outDir, 'Comparison figure skipped: %s', ME.message);
end
writeTcReport(cfg, outDir, S1, S4, cmp);
end

function writeTcReport(cfg, outDir, S1, S4, cmp)
params = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideLength'};
nice = {'Walking speed', 'Step time', 'Stride duration', 'Stride length'};
N = numel(cmp);
L = {};
L{end+1} = 'TC ESTIMATE REPORT (written automatically)';
L{end+1} = sprintf('This run     : %s', outDir);
L{end+1} = sprintf('Previous run : %s (true Tc given to the extractor)', cfg.PrevRun);
L{end+1} = sprintf('Tc method    : %s. Pass 1 nominal Tc = %g s. Simulated data only. Noise seeds start at %d (see Seeds_Used.csv).', cfg.TcMethod, cfg.TcNominal, cfg.Seed);
if cfg.Quick
    L{end+1} = 'QUICK TEST: 2 speeds, 3 trials, no diagnostics. The numbers below are only a check that the script runs, not results.';
end
L{end+1} = 'Only one thing differs from the previous run: the extractor is given Tc_est instead of the true Tc. The scorer (matching tolerance min(0.2, Tc/4) and the stride duration reference) still uses the true Tc.';
L{end+1} = '';

L{end+1} = '1. HOW GOOD IS THE Tc ESTIMATE';
L{end+1} = '  error % = 100*(Tc_est - Tc_true)/Tc_true. mean = error of the average estimate. worst = largest single-trial error ignoring sign. A fallback means the nominal 1.2 s was used because too few steps were found.';
L{end+1} = '  A. Validation runs (Stage 4, fresh noise, each parameter at its own chosen setting)';
L{end+1} = '  speed | true Tc s | mean Tc_est s | mean error % | worst error % | fallbacks / runs';
totalFb = 0;
totalRuns = 0;
for v0 = cfg.Speeds
    r = [subRows(S4, 'StepTime', v0), subRows(S4, 'StrideDuration', v0), subRows(S4, 'StrideLength', v0)];
    if isempty(r)
        continue;
    end
    nRuns = sum([r.Trials]);
    fb = sum([r.Tc_fallbacks]);
    totalFb = totalFb + fb;
    totalRuns = totalRuns + nRuns;
    L{end+1} = sprintf('  %.1f | %.3f | %.3f | %+.2f | %.2f | %d / %d', v0, r(1).Tc_true_s, meanFinite([r.Tc_est_mean_s]), ...
        meanFinite([r.Tc_err_pct]), max([r.Tc_worst_abs_err_pct]), fb, nRuns);
end
L{end+1} = sprintf('  Fallbacks in Stage 4: %d of %d runs.', totalFb, totalRuns);
L{end+1} = '  B. CPI sweep (Stage 1, 30 dB, every step CPI)';
L{end+1} = '  speed | true Tc s | mean Tc_est s | mean error % | worst error % | fallbacks / runs';
fb1All = 0;
runs1All = 0;
for v0 = cfg.Speeds
    r = subRows(S1, 'StepTime', v0);
    if isempty(r)
        continue;
    end
    fb = sum([r.Tc_fallbacks]);
    nRuns = sum([r.Trials]);
    fb1All = fb1All + fb;
    runs1All = runs1All + nRuns;
    L{end+1} = sprintf('  %.1f | %.3f | %.3f | %+.2f | %.2f | %d / %d', v0, r(1).Tc_true_s, meanFinite([r.Tc_est_mean_s]), ...
        meanFinite([r.Tc_err_pct]), max([r.Tc_worst_abs_err_pct]), fb, nRuns);
end
L{end+1} = sprintf('  Fallbacks in Stage 1: %d of %d runs.', fb1All, runs1All);
L{end+1} = '';

pPrev = [cmp.Pass_prev];
pNew = [cmp.Pass_new];
L{end+1} = sprintf('2. CONDITIONS BELOW %g %% RMSE (Stage 4, %d conditions compared)', cfg.Target_pct, N);
L{end+1} = sprintf('  Previous run (true Tc)     : %d of %d', sum(pPrev), N);
L{end+1} = sprintf('  This run (estimated Tc)    : %d of %d', sum(pNew), N);
for p = 1:4
    k = strcmp({cmp.Parameter}, params{p});
    L{end+1} = sprintf('    %s: %d of %d before, %d of %d now', nice{p}, sum(pPrev(k)), nnz(k), sum(pNew(k)), nnz(k));
end
for tag = {'prev', 'new'}
    f = find(~([cmp.(['Pass_' tag{1}])] == 1));
    if isempty(f)
        L{end+1} = sprintf('  Not below the target (%s run): none.', tag{1});
    else
        txt = '';
        for q = f
            txt = [txt sprintf('%s %.1f m/s (%.2f %%); ', cmp(q).Parameter, cmp(q).Speed_mps, cmp(q).(['RMSE_pct_' tag{1}]))]; %#ok<AGROW>
        end
        L{end+1} = sprintf('  Not below the target (%s run): %s', tag{1}, txt);
    end
end
L{end+1} = '';

L{end+1} = '3. CONDITIONS THAT CHANGED BY MORE THAN THE RUN-TO-RUN NOISE (or changed pass/fail)';
L{end+1} = '  Noise limits (pp of RMSE): walking speed 0.09, step time 0.31, stride duration 0.78, stride length 0.28 at 30 dB falling to 0.05 at 60 dB (straight line between, the lower of the two runs'' SNR is used).';
L{end+1} = '  parameter | speed | RMSE % before | RMSE % now | change pp | noise limit pp | better or worse | pass before -> now | settings same? (CPI, SNR)';
chg = find(~strcmp({cmp.Change}, 'none'));
if isempty(chg)
    L{end+1} = '  None. Every condition is within the noise and no pass/fail flag changed.';
end
for q = chg
    c = cmp(q);
    same = isequaln(c.CPIspeed_prev, c.CPIspeed_new) && isequaln(c.CPIstep_prev, c.CPIstep_new) && isequaln(c.SNR_prev, c.SNR_new);
    L{end+1} = sprintf('  %s | %.1f | %.2f | %.2f | %+.2f | %.2f | %s | %s -> %s | %s', c.Parameter, c.Speed_mps, c.RMSE_pct_prev, c.RMSE_pct_new, ...
        c.RMSE_pct_diff, c.Noise_limit_pp, c.Change, yesNo(c.Pass_prev == 1), yesNo(c.Pass_new == 1), yesNo(same));
end
nSetDiff = 0;
for q = 1:N
    c = cmp(q);
    if ~(isequaln(c.CPIspeed_prev, c.CPIspeed_new) && isequaln(c.CPIstep_prev, c.CPIstep_new) && isequaln(c.SNR_prev, c.SNR_new))
        nSetDiff = nSetDiff + 1;
    end
end
L{end+1} = sprintf('  The sweep picks its own CPI and SNR in every run. In %d of %d conditions the chosen setting differs from the previous run, which can move a result without Tc being the cause.', nSetDiff, N);
ws = strcmp({cmp.Parameter}, 'WalkingSpeed');
L{end+1} = sprintf('  Walking speed does not use Tc at all. Its largest change is %.2f pp (noise limit 0.09 pp), so any change there shows the noise and setting effect rather than a Tc effect.', max(abs([cmp(ws).RMSE_pct_diff])));
border = find(abs([cmp.RMSE_pct_prev] - cfg.Target_pct) < 1 | abs([cmp.RMSE_pct_new] - cfg.Target_pct) < 1);
if ~isempty(border)
    txt = '';
    for q = border
        txt = [txt sprintf('%s %.1f m/s (%.2f -> %.2f); ', cmp(q).Parameter, cmp(q).Speed_mps, cmp(q).RMSE_pct_prev, cmp(q).RMSE_pct_new)]; %#ok<AGROW>
    end
    L{end+1} = ['  Within 1 pp of the target line in at least one run, so a pass/fail flip here may be noise: ' txt];
end
L{end+1} = '';

lost = find(pPrev == 1 & pNew == 0);
gained = find(pPrev == 0 & pNew == 1);
worse = find(strcmp({cmp.Change}, 'worse'));
better = find(strcmp({cmp.Change}, 'better'));
changed = find(~strcmp({cmp.Change}, 'none'));
L{end+1} = '4. VERDICT';
if isempty(changed)
    verdict = 'SAME RESULT: no pass/fail changes and no change beyond the run-to-run noise.';
elseif ~isempty(worse) && numel(lost) <= 2 && numel(changed) <= 6
    verdict = ['SMALL COST: ' condList(cmp, worse) ' got worse' impNote(cmp, better) '.'];
else
    verdict = ['DIFFERENT: worse: ' condList(cmp, worse) '. Better: ' condList(cmp, better) '.'];
end
L{end+1} = ['  ' verdict];
L{end+1} = '  Rule used: SAME RESULT = nothing changed beyond noise and no pass/fail flag moved. SMALL COST = something got worse, at most 2 passes were lost and at most 6 conditions changed in total. Anything else, including a run that only got better, is called DIFFERENT.';
L{end+1} = '';
L{end+1} = 'NOTES ON THE COMPARISON';
L{end+1} = '  The previous run used no random seed, so its noise cannot be repeated. This run is seeded (Seeds_Used.csv). Differences inside the noise limits are therefore expected even if Tc made no difference.';
L{end+1} = '  The previous run was made with Master_SimSweep_v4, which did not save MAPE or Accuracy. Those two columns are blank for the previous run in TcEstimate_Comparison.csv. Set PrevRun to a Master_SimSweep_v5 folder to fill them.';
L{end+1} = '  Walking speed rows have no Tc columns because the walking speed script does not use Tc.';
if strcmp(cfg.TcMethod, 'torso')
    L{end+1} = sprintf('  The torso method asks the walking speed script for a track at a requested window of %g s, as RealData_FullAnalysis_v6_params does with CPI_torso_s.', cfg.TorsoCPIreq);
end
fid = fopen(fullfile(outDir, 'TcEstimate_Report.txt'), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
logMsg(outDir, 'Verdict: %s', verdict);
end

function s = condList(cmp, idx)
if isempty(idx)
    s = 'none';
    return;
end
s = '';
for q = idx
    s = [s sprintf('%s at %.1f m/s (%+.2f pp)', cmp(q).Parameter, cmp(q).Speed_mps, cmp(q).RMSE_pct_diff)]; %#ok<AGROW>
    if q ~= idx(end)
        s = [s ', ']; %#ok<AGROW>
    end
end
end

function s = impNote(cmp, idx)
if isempty(idx)
    s = '';
else
    s = ['. Better: ' condList(cmp, idx)];
end
end

function plotTcComparison(cfg, outDir, cmp)
params = {'WalkingSpeed', 'StepTime', 'StrideDuration', 'StrideLength'};
nice = {'Walking speed', 'Step time', 'Stride duration', 'Stride length'};
cPrev = [0.61 0.64 0.69];                                    % grey: previous run
cNew = [0.06 0.46 0.43];                                     % teal: this run
fig = newLightFigure();
set(fig, 'Position', [100 100 1200 800]);
for p = 1:4
    ax = subplot(2, 2, p, 'Parent', fig);
    hold(ax, 'on');
    sub = cmp(strcmp({cmp.Parameter}, params{p}));
    [~, o] = sort([sub.Speed_mps]);
    sub = sub(o);
    vals = [[sub.RMSE_pct_prev]', [sub.RMSE_pct_new]'];
    yTop = min(30, max(10, 1.15*max(vals(:))));              % axis stops at 30 %, anything higher is flagged
    shown = min(vals, yTop);                                 % clipped at the top of the axis
    if size(shown, 1) == 1
        shown(2, :) = NaN;                                   % bar needs at least two groups to draw two bars side by side
    end
    b = bar(ax, 1:size(shown, 1), shown, 'grouped');
    set(b(1), 'FaceColor', cPrev);
    set(b(2), 'FaceColor', cNew);
    hT = plot(ax, [0.4 numel(sub) + 0.6], [cfg.Target_pct cfg.Target_pct], 'k--', 'LineWidth', 1.5);
    set(ax, 'XTick', 1:numel(sub), 'XTickLabel', arrayfun(@(s) sprintf('%g', s), [sub.Speed_mps], 'UniformOutput', false), 'FontSize', 11);
    xlim(ax, [0.4 numel(sub) + 0.6]);
    ylim(ax, [0 yTop]);
    grid(ax, 'on');
    xlabel(ax, 'Walking speed (m/s)', 'FontSize', 12);
    ylabel(ax, 'RMSE (% of the true value)', 'FontSize', 12);
    title(ax, nice{p}, 'FontSize', 13, 'FontWeight', 'bold');
    if any(vals(:) > yTop)
        text(ax, 0.99, 0.97, sprintf('Values above %g %% are clipped. They are in the CSV file.', yTop), 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 9);
    end
    if p == 1
        legend(ax, [b(1) b(2) hT], {'Previous run (true Tc)', 'This run (estimated Tc)', sprintf('Target (%g %%)', cfg.Target_pct)}, 'Location', 'northeast', 'FontSize', 10);
    end
end
annotation(fig, 'textbox', [0.02 0.94 0.96 0.05], 'String', ...
    sprintf('Validation RMSE with the gait cycle time Tc taken from the simulator (previous run) and estimated from the radar signal (%s method). Simulated data, Stage 4, fresh noise.', cfg.TcMethod), ...
    'EdgeColor', 'none', 'FontSize', 11, 'HorizontalAlignment', 'center');
saveLightFigure(fig, fullfile(outDir, 'TcEstimate_Comparison'));
close(fig);
end
