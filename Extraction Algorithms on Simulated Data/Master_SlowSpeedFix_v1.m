%% Torso-band comparison
% Zahraa Sattar | EEE4022S
%
% Compare the original torso band with narrower bands and longer windows.
% Use separate noise seeds for selection and validation.

function Master_SlowSpeedFix_v1(userCfg)

cfg = defaultConfig();
if nargin > 0 && isstruct(userCfg)                           % overwrite any default the user passes in
    f = fieldnames(userCfg);
    for q = 1:numel(f)
        cfg.(f{q}) = userCfg.(f{q});
    end
end
codeFolder = fileparts(mfilename('fullpath'));               % folder this script sits in (the radar folder)
if isempty(codeFolder)
    codeFolder = pwd;
end
oldPath = path;
addpath(codeFolder);                                         % lets the working copies find Zahraa_Spectrogram.m
pathCleanup = onCleanup(@() path(oldPath));                  %#ok<NASGU>
assert(exist('Zahraa_Spectrogram', 'file') ~= 0, 'Zahraa_Spectrogram.m must sit in %s', codeFolder);

outDir = fullfile(codeFolder, ['SlowSpeedFix_' datestr(now, 'yyyymmdd_HHMMSS')]);
mkdir(outDir);
P = makeWorkingCopies(cfg, codeFolder, outDir);
logMsg(outDir, 'Started. Speeds %s m/s, %d draws, SNRs %s dB.', mat2str(cfg.Speeds), cfg.nTrials, mat2str(cfg.SNRs));

S = settingList(cfg);                                        % all the settings that get tried
tStart = tic;

selRows = runGrid(P, cfg, S, cfg.seedA, outDir, 'selection', tStart);
writeRows(fullfile(outDir, 'Fix_Selection.csv'), selRows);

score = zeros(1, numel(S));                                  % passes over speed x SNR x (walking speed, stride length)
errSum = zeros(1, numel(S));
for s = 1:numel(S)
    m = strcmp({selRows.Setting}, S(s).name);
    score(s) = sum([selRows(m).WS_pass]) + sum([selRows(m).SL_pass]);
    errSum(s) = mean([selRows(m).WS_rmse, selRows(m).SL_rmse], 'omitnan');
end
best = 1;
for s = 2:numel(S)                                           % most passes wins, lowest mean error breaks a tie
    if score(s) > score(best) || (score(s) == score(best) && errSum(s) < errSum(best))
        best = s;
    end
end
logMsg(outDir, 'Picked setting: %s (%d passes, baseline had %d).', S(best).name, score(best), score(1));

runList = unique([1, best]);                                 % baseline is always setting 1
valRows = runGrid(P, cfg, S(runList), cfg.seedV, outDir, 'validation', tStart);
writeRows(fullfile(outDir, 'Fix_Validation.csv'), valRows);
writeReport(outDir, cfg, S, selRows, valRows, score, best);
logMsg(outDir, 'Finished in %.1f min. Open Fix_Report.txt in %s', toc(tStart)/60, outDir);
end

function cfg = defaultConfig()
cfg.Speeds = [0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];      % m/s, all of them so any new problem at another speed shows up
cfg.nTrials = 30;                                            % noise draws per condition
cfg.SNRs = [30 40 60];                                       % dB
cfg.seedA = 20264001;                                        % seed for picking the setting
cfg.seedV = 20265002;                                        % different seed for the validation
cfg.R0_m = 10;                                               % same range as the Master
cfg.AmplitudeMode = 'A';                                     % same amplitude model as the Master
cfg.NoiseMode = 'extractorN';                                % noise power from the extractor window, as in the Master
cfg.StepCPIreq = 0.1;                                        % step window 0.186 s, the one fixed step setting
cfg.Target_pct = 5;                                          % goal for RMSE %
cfg.SimName = 'RadarSimulator_v11_torso_legs_arms_head_buffered_noplots';
cfg.WSName = 'Finding_WalkingSpeed_v4_noplots';
cfg.STName = 'Finding_StepRecognition_v5_noplots';
cfg.CopyHook = [];                                           % only used for testing outside MATLAB
end

function S = settingList(cfg) %#ok<INUSD>
list = { ...
    'base_0.372',        0.2, false, 0,    3; ...            % the current setting, always first
    'base_0.743',        0.4, false, 0,    3; ...
    'band25_0.372',      0.2, false, 0.25, 3; ...
    'band25_0.743',      0.4, false, 0.25, 3; ...
    'band25_1.0',        1.0, true,  0.25, 3; ...
    'band15_0.372',      0.2, false, 0.15, 3; ...
    'band15_0.743',      0.4, false, 0.15, 3; ...
    'band40_0.743',      0.4, false, 0.40, 3};
S = struct('name', {}, 'cpi', {}, 'exact', {}, 'bandA', {}, 'passes', {});
for k = 1:size(list, 1)
    S(k) = struct('name', list{k, 1}, 'cpi', list{k, 2}, 'exact', list{k, 3}, 'bandA', list{k, 4}, 'passes', list{k, 5});
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

function logMsg(outDir, fmt, varargin)
msg = sprintf(fmt, varargin{:});
line = sprintf('[%s] %s', datestr(now, 'HH:MM:SS'), msg);
fprintf('%s\n', line);
fid = fopen(fullfile(outDir, 'Run_Log.txt'), 'a');
if fid > 0
    fprintf(fid, '%s\n', line);
    fclose(fid);
end
end

function P = makeWorkingCopies(cfg, codeFolder, outDir)
codeDir = fullfile(outDir, 'WorkingCopies');
mkdir(codeDir);
P.sim = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.SimName '.m']), fullfile(codeDir, 'WC_Simulator.m'), 'sim', cfg);
P.ws = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.WSName '.m']), fullfile(codeDir, 'WC_WalkingSpeed.m'), 'ws', cfg);
P.st = copyOne(fullfile(codeFolder, 'Report Sources', [cfg.STName '.m']), fullfile(codeDir, 'WC_StepRecognition.m'), 'st', cfg);
end

function dest = copyOne(src, dest, kind, cfg)
assert(isfile(src), 'Missing file: %s', src);
txt = fileread(src);
txt = regexprep(txt, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*(?:\r?\n|$)', '');   % no clear all or close all inside a loop
if ~strcmp(kind, 'sim')
    txt = regexprep(txt, '(?m)^[ \t]*RadarSimulator_[A-Za-z0-9_]+;[^\r\n]*(?:\r?\n|$)', '');   % the simulator is run separately
end
if strcmp(kind, 'st')
    plotMarker = 'figure(''Name'',''Step 6 verification - remove outliers'');';
    a = strfind(txt, plotMarker);
    b = strfind(txt, '%% Step 7');
    assert(numel(a) == 1 && numel(b) == 1 && b > a, 'Step 6 figure not found in the step script.');
    txt = [txt(1:a-1), txt(b:end)];                          % cut the Step 6 figure
    hopOld = 'noverlap = N*3/4;';
    assert(numel(strfind(txt, hopOld)) == 1, 'Overlap line not found in the step script.');
    hopNew = sprintf(['global StepHopFraction\n' ...
        '        hopFrac_ = 0.25;\n' ...
        '        if ~isempty(StepHopFraction), hopFrac_ = StepHopFraction; end\n' ...
        '        noverlap = N - max(1, round(N*hopFrac_));']);   % default stays the original 75 percent overlap
    txt = strrep(txt, hopOld, hopNew);
    cpkPattern = 'cPeak\s*=\s*DopplerRes_Hz/freq_high_Yardibi;[^\r\n]*';
    assert(numel(regexp(txt, cpkPattern, 'match')) == 1, 'cPeak line not found in the step script.');
    cpkNew = sprintf(['global StepHopFraction\n' ...
        'cPeak = DopplerRes_Hz/freq_high_Yardibi;\n' ...
        'if ~isempty(StepHopFraction), cPeak = cPeak*StepHopFraction/0.25; end']);
    if isempty(strfind(txt, 'cPeakScale'))
        cpkNew = [cpkNew sprintf('\nglobal StepCPeakScale\nif ~isempty(StepCPeakScale), cPeak = cPeak*StepCPeakScale; end')];
    end
    txt = regexprep(txt, cpkPattern, cpkNew);
end
if strcmp(kind, 'ws')
    txt = patchRegex(txt, 'N_actual = 2\^\(nextpow2\(round\(CPI_s\*fs\)\)\);[^\r\n]*', ...
        sprintf(['global VX\nN_actual = 2^(nextpow2(round(CPI_s*fs)));\nif VX.WSexact, N_actual = round(CPI_s*fs); end']), 'actual CPI line');
    txt = patchRegex(txt, 'numIterations\s*=\s*3;[^\r\n]*', ...
        sprintf(['global VX\nnumIterations = VX.WSpasses;']), 'number of passes');
    txt = patchRegex(txt, ['N = 2\^\(nextpow2\(round\(CPI_s\*fsProc\)\)\);[^\r\n]*\r?\n[ \t]*win[ \t]*=[ \t]*hamming\(N\);[^\r\n]*\r?\n' ...
        '[ \t]*noverlap = N\*3/4;[^\r\n]*\r?\n[ \t]*nfft[ \t]*=[ \t]*N;[^\r\n]*'], ...
        sprintf(['global VX\nN = 2^(nextpow2(round(CPI_s*fsProc)));\nif VX.WSexact, N = round(CPI_s*fsProc); end\n' ...
        'win = hamming(N);\nnoverlap = N - max(1, round(N*0.25));\nnfft = 2^(nextpow2(N));']), 'spectrogram block');
    txt = patchRegex(txt, 'Ssmooth = SmoothSpectrogramNoise\(Power_S, 21\);[^\r\n]*', ...
        sprintf(['global VX\nsmoothBins_ = 21;\nif VX.WSexact\n    Nwin_ = round(CPI_s*fsProc);\n' ...
        '    smoothBins_ = 2*floor(21*(2^nextpow2(Nwin_))/Nwin_/2) + 1;\nend\n' ...
        'Ssmooth = SmoothSpectrogramNoise(Power_S, smoothBins_);']), 'bin smoothing line');
    txt = patchRegex(txt, 'freq_low\s*=\s*max\(0\.5\*centreFreq_Hz,freq_low0\);[^\r\n]*', ...
        sprintf(['global VX\nif VX.WSbandA > 0\n    freq_low = max((1 - VX.WSbandA)*centreFreq_Hz, freq_low0);\nelse\n' ...
        '    freq_low = max(0.5*centreFreq_Hz,freq_low0);\nend']), 'band lower edge');
    txt = patchRegex(txt, 'freq_high\s*=\s*max\(2\.0\*centreFreq_Hz,f_high_min\);[^\r\n]*', ...
        sprintf(['global VX\nif VX.WSbandA > 0\n    freq_high = (1 + VX.WSbandA)*centreFreq_Hz;\nelse\n' ...
        '    freq_high = max(2.0*centreFreq_Hz,f_high_min);\nend']), 'band upper edge');
end
if ~isempty(cfg.CopyHook)
    txt = cfg.CopyHook(txt, kind);
end
fid = fopen(dest, 'w');
assert(fid >= 0, 'Could not create %s', dest);
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function txt = patchRegex(txt, pattern, replacement, label)
m = regexp(txt, pattern, 'match');
if numel(m) ~= 1
    error('Patch "%s": expected 1 match, found %d. The extractor script has changed, nothing was patched.', label, numel(m));
end
rep = strrep(replacement, '\', '\\');
rep = strrep(rep, '$', '\$');
txt = regexprep(txt, pattern, rep);
end

function rows = runGrid(P, cfg, S, seedBase, outDir, tag, tStart)
rows = struct('Stage', {}, 'Speed', {}, 'Setting', {}, 'SNR', {}, 'WS_rmse', {}, 'WS_mean', {}, 'WS_spread', {}, 'WS_pass', {}, ...
    'SL_rmse', {}, 'SL_mean', {}, 'SL_pass', {}, 'Failed', {});
for v0 = cfg.Speeds
    logMsg(outDir, '%s, v0 = %.2f m/s, elapsed %.1f min', tag, v0, toc(tStart)/60);
    D = getSpeedData(P, cfg, v0);                            % the simulator runs once per speed
    setSeed(seedBase + round(1000*v0));
    z = cell(1, cfg.nTrials);                                % same unit noise draws for every setting and SNR
    for k = 1:cfg.nTrials
        z{k} = randn(1, numel(D.clean)) + 1i*randn(1, numel(D.clean));
    end
    for snr = cfg.SNRs
        Nst = noiseWindow(cfg, v0, cfg.StepCPIreq, D.fs);
        Sst = cell(1, cfg.nTrials);
        for k = 1:cfg.nTrials                                % step extractor runs once per trial, all settings share it
            Sst{k} = runStep(P, addNoise(D.clean, Nst, snr, z{k}), D, cfg.StepCPIreq, snr);
        end
        for s = 1:numel(S)
            Nw = noiseWindow(cfg, v0, S(s).cpi, D.fs);
            est = NaN(1, cfg.nTrials);
            E = cell(cfg.nTrials, 1);
            G = cell(cfg.nTrials, 1);
            nFail = 0;
            for k = 1:cfg.nTrials
                W = runSpeed(P, addNoise(D.clean, Nw, snr, z{k}), D, S(s));
                if W.ok
                    est(k) = W.meanSpeed;
                else
                    nFail = nFail + 1;
                end
                [E{k}, G{k}] = strideLengthErrors(W, Sst{k}, D);
            end
            r.Stage = tag;
            r.Speed = v0;
            r.Setting = S(s).name;
            r.SNR = snr;
            e = est(isfinite(est)) - v0;
            if isempty(e)
                r.WS_rmse = NaN; r.WS_mean = NaN; r.WS_spread = NaN;
            else
                r.WS_rmse = 100*sqrt(mean(e.^2))/v0;         % RMSE as % of the true speed
                r.WS_mean = 100*mean(e)/v0;                  % bias
                r.WS_spread = 100*std(e)/v0;                 % spread between draws
            end
            r.WS_pass = double(isfinite(r.WS_rmse) && r.WS_rmse < cfg.Target_pct);
            [r.SL_rmse, r.SL_mean] = slStats(E, G);
            r.SL_pass = double(isfinite(r.SL_rmse) && r.SL_rmse < cfg.Target_pct);
            r.Failed = nFail;
            rows(end+1) = r;                                 %#ok<AGROW>
        end
    end
end
end

function [rm, me] = slStats(E, G)
trialRMSE = [];
allE = [];
allG = [];
for k = 1:numel(E)
    e = E{k}(:).';
    g = G{k}(:).';
    ok = isfinite(e) & isfinite(g) & g > 0;
    e = e(ok);
    g = g(ok);
    if ~isempty(e)
        trialRMSE(end+1) = sqrt(mean(e.^2));                 %#ok<AGROW>
        allE = [allE, e];                                    %#ok<AGROW>
        allG = [allG, g];                                    %#ok<AGROW>
    end
end
if isempty(allE)
    rm = NaN;
    me = NaN;
else
    rm = 100*mean(trialRMSE)/mean(allG);
    me = 100*mean(allE)/mean(allG);
end
end

function writeRows(file, rows)
fid = fopen(file, 'w');
fprintf(fid, 'Stage,Speed_mps,Setting,SNR_dB,WS_RMSE_pct,WS_MeanErr_pct,WS_Spread_pct,WS_Pass,SL_RMSE_pct,SL_MeanErr_pct,SL_Pass,FailedTrials\n');
for q = 1:numel(rows)
    r = rows(q);
    fprintf(fid, '%s,%.2f,%s,%g,%.3f,%.3f,%.3f,%d,%.3f,%.3f,%d,%d\n', r.Stage, r.Speed, r.Setting, r.SNR, ...
        r.WS_rmse, r.WS_mean, r.WS_spread, r.WS_pass, r.SL_rmse, r.SL_mean, r.SL_pass, r.Failed);
end
fclose(fid);
end

function writeReport(outDir, cfg, S, selRows, valRows, score, best)
fid = fopen(fullfile(outDir, 'Fix_Report.txt'), 'w');
fprintf(fid, 'SLOW SPEED FIX TEST (written automatically)\n');
fprintf(fid, 'Speeds %s m/s, %d draws per condition, SNRs %s dB. Target: RMSE below %g %%.\n', mat2str(cfg.Speeds), cfg.nTrials, mat2str(cfg.SNRs), cfg.Target_pct);
fprintf(fid, 'Seed A = %d (picking), seed V = %d (validation, quote these). Simulated data only.\n\n', cfg.seedA, cfg.seedV);
fprintf(fid, 'STEP 1: passes out of %d (walking speed + stride length, all speeds and SNRs), noise seed A\n', 2*numel(cfg.Speeds)*numel(cfg.SNRs));
for s = 1:numel(S)
    m = strcmp({selRows.Setting}, S(s).name);
    m04 = m & abs([selRows.Speed] - 0.4) < 1e-6;
    fprintf(fid, '  %-16s total %2d | at 0.4 m/s: walking speed RMSE %s, stride length RMSE %s (30 / 40 / 60 dB)\n', S(s).name, score(s), ...
        numList([selRows(m04).WS_rmse]), numList([selRows(m04).SL_rmse]));
end
fprintf(fid, '  Picked: %s\n\n', S(best).name);
fprintf(fid, 'STEP 2: fresh noise (seed V), baseline against the picked setting\n');
names = {S(1).name, S(best).name};
names = unique(names, 'stable');
for n = 1:numel(names)
    m = strcmp({valRows.Setting}, names{n});
    fprintf(fid, '\n  %s  (walking speed passes %d of %d, stride length passes %d of %d)\n', names{n}, ...
        sum([valRows(m).WS_pass]), sum(m), sum([valRows(m).SL_pass]), sum(m));
    fprintf(fid, '    speed | SNR | WS RMSE | WS mean err | SL RMSE | SL mean err\n');
    q = find(m);
    for k = q
        r = valRows(k);
        fprintf(fid, '    %.1f   | %3g | %7.2f | %+11.2f | %7.2f | %+10.2f\n', r.Speed, r.SNR, r.WS_rmse, r.WS_mean, r.SL_rmse, r.SL_mean);
    end
end
if numel(names) == 2
    a = strcmp({valRows.Setting}, names{1});
    b = strcmp({valRows.Setting}, names{2});
    A = valRows(a);
    B = valRows(b);
    fprintf(fid, '\nCHANGES (picked setting against baseline, fresh noise):\n');
    for k = 1:numel(A)
        for f = {'WS', 'SL'}
            pa = A(k).([f{1} '_pass']);
            pb = B(k).([f{1} '_pass']);
            if pa == 1 && pb == 0
                fprintf(fid, '  WORSE  %s at %.1f m/s, %g dB: now %.2f %%\n', f{1}, A(k).Speed, A(k).SNR, B(k).([f{1} '_rmse']));
            elseif pa == 0 && pb == 1
                fprintf(fid, '  FIXED  %s at %.1f m/s, %g dB: now %.2f %%\n', f{1}, A(k).Speed, A(k).SNR, B(k).([f{1} '_rmse']));
            end
        end
    end
    fprintf(fid, 'Total passes: baseline %d, picked %d (out of %d).\n', sum([A.WS_pass]) + sum([A.SL_pass]), sum([B.WS_pass]) + sum([B.SL_pass]), 2*numel(A));
end
fclose(fid);
end

function t = numList(x)
t = '';
for k = 1:numel(x)
    t = [t sprintf('%.1f ', x(k))];                          %#ok<AGROW>
    if k < numel(x)
        t = [t '/ '];                                        %#ok<AGROW>
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
end

function Nw = noiseWindow(cfg, v0, cpiReq, fs)
Nw = 2^nextpow2(round(cpiReq*fs));                           % noise power uses the extractor window length
if strcmp(cfg.NoiseMode, 'simulatorN') && v0 < 1.0
    Nw = 2^nextpow2(round(0.075*fs));
end
end

function x = addNoise(clean, Nw, snr_dB, z)
noisePower = 10^(-(snr_dB - 10*log10(Nw))/10);               % same equation as the simulator
x = clean + sqrt(noisePower/2)*z;
end

function m_o = runSpeed(m_P, m_x, m_D, m_vs)
global VX
VX = struct('WSexact', m_vs.exact, 'WSbandA', m_vs.bandA, 'WSpasses', m_vs.passes);
m_o = struct('ok', false, 'meanSpeed', NaN, 'vel', [], 'tv', [], 'cpiAct', NaN, 'err', '');
try
    RdrRxSignal = m_x;                                       %#ok<NASGU>
    fs = m_D.fs;                                             %#ok<NASGU>
    fc = m_D.fc;                                             %#ok<NASGU>
    c = m_D.c;                                               %#ok<NASGU>
    v0 = m_D.v0;                                             %#ok<NASGU>
    CPI_s = m_vs.cpi;                                        %#ok<NASGU>
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

function m_o = runStep(m_P, m_x, m_D, m_cpi, m_snr)
global StepHopFraction StepCPeakScale
StepHopFraction = [];
StepCPeakScale = 0.5;                                        % the half-bin peak drop used in the final settings
m_o = struct('ok', false, 'det', [], 'key', [], 'cpiAct', NaN, 'err', '');
try
    RdrRxSignal = m_x;                                       %#ok<NASGU>
    fs = m_D.fs;                                             %#ok<NASGU>
    Tc = m_D.Tc;                                             %#ok<NASGU>
    CPI_s = m_cpi;                                           %#ok<NASGU>
    SNR_dB = m_snr;                                          %#ok<NASGU>
    cPeakScale = 0.5;                                        %#ok<NASGU>
    stepRefineMode = 'centroid';                             %#ok<NASGU>
    stepGapFill = true;                                      %#ok<NASGU>
    evalc('run(m_P.st);');
    m_o.det = detectedStepTimes_s(:).';
    if exist('stepFrameTimes_s', 'var')
        m_o.key = stepFrameTimes_s(:).';
    else
        m_o.key = m_o.det;
    end
    m_o.cpiAct = CPI_actual_s;
    m_o.ok = true;
catch ME
    m_o.err = ME.message;
end
StepCPeakScale = [];
end

function [e, g] = strideLengthErrors(W, S, D)
e = NaN;
g = NaN;
if ~W.ok || ~S.ok || numel(S.det) < 3
    return;
end
try
    edge = 0.1 + S.cpiAct/2;                                 % same walk interior as the Master
    inside = S.key > D.buffer + edge & S.key < D.buffer + D.simT - edge;
    detIn = S.det(inside);
    if numel(detIn) < 3
        return;
    end
    vel = fillmissing(W.vel, 'linear', 'EndValues', 'nearest');
    dsp = cumtrapz(W.tv, vel);                               % torso position from the speed track
    atDet = interp1(W.tv, dsp, detIn, 'linear');
    est = abs(atDet(3:end) - atDet(1:end-2));
    trueTimes = D.buffer : D.Tc/2 : D.buffer + D.simT;
    trueRange = interp1(D.t, D.Torso_R, trueTimes, 'linear');
    trueLen = abs(trueRange(3:end) - trueRange(1:end-2));
    [em, tm] = matchStrides(est, detIn(1:end-2), trueLen, trueTimes(1:end-2), D.Tc/4);
    e = em - tm;
    g = tm;
catch
    e = NaN;
    g = NaN;
end
end

function [estimatedMatched, trueMatched] = matchStrides(estimated, estimatedStarts, truth, trueStarts, tolerance)
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
