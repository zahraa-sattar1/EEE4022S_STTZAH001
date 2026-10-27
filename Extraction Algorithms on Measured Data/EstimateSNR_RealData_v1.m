%% SNR estimation from measured radar recordings
% Zahraa Sattar | EEE4022S
%
% Calibrate the peak-to-floor ratio on simulated signals, then estimate SNR for each local recording.
% Estimates are reported on the simulator's SNR scale and flagged outside the calibrated range.
% The measurement and recording helpers are included below.

function EstimateSNR_RealData_v1(userCfg)

if nargin < 1
    userCfg = struct();                                      % no settings given, use the defaults
end
cfg = defaultConfig();                                       % default settings
if isfield(userCfg, 'Quick') && userCfg.Quick
    cfg = quickConfig(cfg);                                  % small run just to check everything works
end
names = fieldnames(userCfg);
for k = 1:numel(names)
    cfg.(names{k}) = userCfg.(names{k});                     % anything passed in overrides the defaults
end

codeFolder = cfg.CodeFolder;
if isempty(codeFolder)
    codeFolder = fileparts(mfilename('fullpath'));           % folder this script sits in
end
if isempty(codeFolder)
    codeFolder = pwd;
end
oldPath = path;
addpath(codeFolder);                                         % so the simulator, Zahraa_Spectrogram and the loader are found
oldDir = pwd;
cleanup = onCleanup(@() restoreState(oldPath, oldDir));     % put the path and current folder back at the end, even after an error
cd(codeFolder);                                              % the real-data loader looks for the recordings relative to this folder

outDir = fullfile(codeFolder, ['SNREstimate_' datestr(now, 'yyyymmdd_HHMMSS')]);
mkdir(outDir);                                               % always a new folder, nothing gets overwritten
save(fullfile(outDir, 'Config.mat'), 'cfg');
logMsg(outDir, 'EstimateSNR_RealData_v1 started. Quick mode: %d', cfg.Quick);
tStart = tic;

simCopy = makeSimCopy(cfg, codeFolder, outDir);              % working copy of the simulator, original untouched

%% Self test
selfTest(cfg, simCopy, outDir);

%% Part A: calibration and validation on the simulator
[calRows, valRows] = runSimulator(cfg, simCopy, outDir);
writeTable(fullfile(outDir, 'Sim_Calibration.csv'), calRows);
writeTable(fullfile(outDir, 'Sim_Validation.csv'), valRows);
cal = fitCalibration(cfg, calRows);                          % offset, margin and the SNR range where the estimate can be trusted
val = checkValidation(cfg, cal, valRows);                    % how well the fitted estimate does on fresh draws
writeSummary(fullfile(outDir, 'Calibration_Summary.csv'), cal);
logMsg(outDir, 'Calibration done. Offset I/Q %.2f dB, I channel %.2f dB, margin %.2f dB.', cal.offsetC, cal.offsetR, cal.margin);

%% Part B: real walks
realRows = struct([]);
if cfg.RunRealData
    realRows = runRealData(cfg, cal, outDir);
    writeTable(fullfile(outDir, 'RealData_SNR_Pth.csv'), realRows);
end

%% Figures and report
if cfg.MakePlots
    makePlots(cfg, cal, calRows, valRows, realRows, outDir);
end
writeReport(cfg, cal, val, realRows, outDir, toc(tStart));
logMsg(outDir, 'Finished in %.1f min. Results in %s', toc(tStart)/60, outDir);
end

function cfg = defaultConfig()
cfg.Speeds = [0.4 0.5 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0];      % m/s, the same ten speeds as Chapter 4
cfg.SNRlist = 0:10:90;                                       % dB, simulator SNR settings, 0 dB = signal and noise equally strong
cfg.nCal = 5;                                                % noise draws per speed and SNR used to fit the estimate
cfg.nVal = 5;                                                % fresh draws (other seed) used to check it
cfg.FitRange = [20 50];
cfg.LinearTol = 2;
cfg.R0_m = 10;                                               % simulator range, same as Chapter 4
cfg.AmplitudeMode = 'A';                                     % Chen amplitude model A, same as Chapter 4
cfg.fc = 24e9;                                               % Hz, IPM-365 carrier
cfg.StepCPIreq = 0.1;                                        % s, step-detection CPI, rounds up to the 0.186 s window used everywhere
cfg.TcMedian = 1.2;
cfg.NoiseBand = [0.5 1.0];                                   % fraction of the leg band limit used as the noise region (724 to 1448 Hz)
cfg.PthCap = -40;                                            % dB, the simulator rule's floor
cfg.PthFixed = -45;                                          % dB, the fixed threshold the real-data runs use now
cfg.SeedBase = 20261010;                                     % fixed seed so reruns match
cfg.GroundTruthFile = fullfile(fileparts(mfilename('fullpath')), 'Report Results', 'Report_References.csv');         % list of real walks
cfg.LoaderName = 'LoadRealRecordingByName_v2';
cfg.RecordingFolder = fullfile(fileparts(mfilename('fullpath')), '..', 'Participant Recordings');               % same loader as the real-data batch runners
cfg.MaxWalks = Inf;                                          % all walks
cfg.RunRealData = true;
cfg.MakePlots = true;
cfg.Quick = false;
cfg.SimName = 'RadarSimulator_v11_torso_legs_arms_head_buffered_noplots';
cfg.CodeFolder = '';
end

function cfg = quickConfig(cfg)
cfg.Speeds = [0.6 1.4];                                      % two speeds only
cfg.SNRlist = [0 20 40 60];                                  % four SNRs only
cfg.nCal = 2;                                                % two draws to fit
cfg.nVal = 1;                                                % one draw to check
cfg.FitRange = [20 40];                                      % fit range has to sit inside the short SNR list
cfg.MaxWalks = 3;                                            % first three real walks only
cfg.Quick = true;
end

function restoreState(oldPath, oldDir)
path(oldPath);                                               % path back to how it was
cd(oldDir);                                                  % and the folder you were in
end

function logMsg(outDir, fmt, varargin)
msg = sprintf(fmt, varargin{:});
stamp = datestr(now, 'HH:MM:SS');
fprintf('[%s] %s\n', stamp, msg);                            % show it in the command window
fid = fopen(fullfile(outDir, 'Run_Log.txt'), 'a');
if fid > 0
    fprintf(fid, '[%s] %s\n', stamp, msg);                   % and keep a copy in the log
    fclose(fid);
end
end

function dest = makeSimCopy(cfg, codeFolder, outDir)
src = fullfile(codeFolder, [cfg.SimName '.m']);
if ~isfile(src)
    src = fullfile(codeFolder, '..', 'Radar Simulator', [cfg.SimName '.m']);
end
assert(isfile(src), 'Missing simulator file: %s', src);
txt = fileread(src);
txt = regexprep(txt, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*(?:\r?\n|$)', '');   % no clear all / close all, so the loop variables survive
codeDir = fullfile(outDir, 'WorkingCopies');
mkdir(codeDir);
dest = fullfile(codeDir, 'WC_Simulator.m');
fid = fopen(dest, 'w');
assert(fid >= 0, 'Could not create %s', dest);
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function D = simulateSpeed(m_sim, m_cfg, v0)
SimulatorVerbose = false;                                    %#ok<NASGU>
SimulatorMakePlots = false;                                  %#ok<NASGU>
SimulatorR0_m = m_cfg.R0_m;                                  %#ok<NASGU>
AmplitudeMode = m_cfg.AmplitudeMode;                         %#ok<NASGU>
SNR_dB = 40;                                                 %#ok<NASGU>
CPI_s = m_cfg.StepCPIreq;                                    %#ok<NASGU>
evalc('run(m_sim);');                                        % run the simulator quietly, v0 is picked up from this workspace
D.v0 = v0;
D.clean = RdrRxSignal_clean(:).';                            % noise-free radar return
D.fs = fs;                                                   % sample rate the simulator used
end

function x = addNoise(clean, Nw, snr_dB, z)
noisePower = 10^(-(snr_dB - 10*log10(Nw))/10);               % same equation as the simulator and the sweeps
x = clean + sqrt(noisePower/2) * z;                          % z is one unit complex noise draw
end

function [calRows, valRows] = runSimulator(cfg, simCopy, outDir)
calRows = struct([]);
valRows = struct([]);
for iv = 1:numel(cfg.Speeds)
    v0 = cfg.Speeds(iv);
    D = simulateSpeed(simCopy, cfg, v0);                     % clean signal for this speed
    Nw = 2^nextpow2(round(cfg.StepCPIreq*D.fs));             % noise scaled to the step window, as in the sweeps
    for setIdx = 1:2                                         % 1 = calibration draws, 2 = validation draws
        if setIdx == 1
            nT = cfg.nCal;
        else
            nT = cfg.nVal;
        end
        for tr = 1:nT
            rng(cfg.SeedBase + 1e6*setIdx + 1000*round(10*v0) + tr);   % own seed per set, speed and trial, so the two sets never share noise
            z = randn(1, numel(D.clean)) + 1i*randn(1, numel(D.clean));   % one noise draw, reused at every SNR (paired)
            for snr = cfg.SNRlist
                x = addNoise(D.clean, Nw, snr, z);
                mC = measureSignal(x, D.fs, cfg);            % complex I/Q, as in Chapter 4
                mR = measureSignal(real(x), D.fs, cfg);      % I channel only, like the IPM-365 recordings
                r.Speed_mps = v0;
                r.SNR_dB = snr;
                r.Trial = tr;
                r.PFR_IQ_dB = mC.PFR_dB;
                r.PFR_I_dB = mR.PFR_dB;
                r.NoiseMed_IQ_dB = mC.NoiseMed_dB;
                r.NoiseMed_I_dB = mR.NoiseMed_dB;
                r.PthRule_dB = max(10 - snr, cfg.PthCap);    % the threshold the simulator rule gives
                r.PthRuleNoFloor_dB = 10 - snr;              % the same rule without the -40 dB floor
                if setIdx == 1
                    calRows = appendRow(calRows, r);
                else
                    valRows = appendRow(valRows, r);
                end
            end
        end
    end
    logMsg(outDir, 'Simulator %.1f m/s done (%d of %d speeds).', v0, iv, numel(cfg.Speeds));
end
end

function m = measureSignal(x, fs, cfg)
legLimit = 350*(cfg.fc/1e9)/5.8;                             % leg band limit, 1448 Hz at 24 GHz, same as step detection
avgSamples = floor(fs/(1.5*legLimit));                       % same average filter length as step detection
xAvg = movmean(x(:).', avgSamples);                          % smooths the raw signal like Step 1
N = 2^nextpow2(round(cfg.StepCPIreq*fs));                    % window length, 8192 samples = 0.186 s at 44.1 kHz
[S, f, t] = Zahraa_Spectrogram(xAvg, hamming(N), N*3/4, N, fs);   % Hamming window, 75 % overlap, no zero padding, like Step 2
inBand = f > 0 & f <= legLimit;                              % positive frequencies up to the leg limit, like Step 3
fBand = f(inBand);
P = abs(S(inBand, :)).^2;                                    % power in the leg band
noiseRows = fBand >= cfg.NoiseBand(1)*legLimit & fBand <= cfg.NoiseBand(2)*legLimit;   % upper half of the band, mostly noise
noiseFloor = median(reshape(P(noiseRows, :), [], 1));        % typical noise power in one cell
m.PFR_dB = 10*log10(max(P(:))/noiseFloor);                   % strongest return compared with the noise floor
Pm = subtractBackground(P);                                  % background subtraction, like Step 3
framePeriod = t(2) - t(1);                                   % time between frames
medianFrames = max(1, round((cfg.TcMedian/8)/framePeriod));  % median filter length, like Step 4
if medianFrames > 1
    Pmedian = movmedian(Pm, medianFrames, 1);                % median filter along frequency, like Step 4
else
    Pmedian = Pm;
end
PmedDb = 10*log10(Pmedian/(max(Pmedian(:)) + eps) + eps);
m.NoiseMed_dB = median(reshape(PmedDb(noiseRows, :), [], 1));   % where the noise typically sits below the maximum
m.FrameMaxNoise_dB = max(PmedDb(noiseRows, :), [], 1);       % loudest noise-band cell in each frame, used for the crossing check
end

function Pm = subtractBackground(P)
nB = max(1, round(0.20*size(P, 2)));                         % quietest 20 % of frames, same as step detection
Ps = sort(P, 2, 'ascend');                                   % each frequency bin sorted over time
Pback = mean(Ps(:, 1:nB), 2);                                % background power per bin
Pm = P - repmat(Pback, 1, size(P, 2));                       % take the background off
Pm(Pm < 0) = eps;                                            % power cannot go negative
end

function cal = fitCalibration(cfg, rows)
snr = [rows.SNR_dB];
inFit = snr >= cfg.FitRange(1) & snr <= cfg.FitRange(2);     % only the straight-line part is used to fit
cal.offsetC = median([rows(inFit).PFR_IQ_dB] - snr(inFit));  % PFR minus SNR, I/Q
cal.offsetR = median([rows(inFit).PFR_I_dB] - snr(inFit));   % PFR minus SNR, I channel only
cal.margin = median([rows(inFit).PthRuleNoFloor_dB] - [rows(inFit).NoiseMed_IQ_dB]);   % how far above the noise the rule puts Pth
cal.SNRs = unique(snr);
for k = 1:numel(cal.SNRs)
    sel = snr == cal.SNRs(k);
    cal.meanPFR_C(k) = mean([rows(sel).PFR_IQ_dB]);          % average PFR at this SNR, I/Q
    cal.sdPFR_C(k) = std([rows(sel).PFR_IQ_dB]);
    cal.meanPFR_R(k) = mean([rows(sel).PFR_I_dB]);           % same for the I channel
    cal.sdPFR_R(k) = std([rows(sel).PFR_I_dB]);
    cal.devC(k) = cal.meanPFR_C(k) - cal.SNRs(k) - cal.offsetC;   % how far this SNR is off the straight line
    cal.devR(k) = cal.meanPFR_R(k) - cal.SNRs(k) - cal.offsetR;
    cal.meanNoise_C(k) = mean([rows(sel).NoiseMed_IQ_dB]);
    cal.meanNoise_R(k) = mean([rows(sel).NoiseMed_I_dB]);
    cal.marginAt(k) = mean(10 - cal.SNRs(k) - [rows(sel).NoiseMed_IQ_dB]);   % the margin at this SNR, should be flat where the noise dominates
end
[cal.loC, cal.hiC] = linearRange(cal.SNRs, cal.devC, cfg);   % SNR range where the I/Q estimate is trustworthy
[cal.loR, cal.hiR] = linearRange(cal.SNRs, cal.devR, cfg);   % same for the I channel
[cal.loM, cal.hiM] = linearRange(cal.SNRs, cal.marginAt - cal.margin, cfg);   % SNR range where the margin stays flat
end

function [lo, hi] = linearRange(snrs, dev, cfg)
mid = find(snrs >= cfg.FitRange(1) & snrs <= cfg.FitRange(2));   % start inside the fit range
if isempty(mid)
    lo = NaN; hi = NaN;
    return;
end
k1 = mid(1);
while k1 > 1 && abs(dev(k1-1)) <= cfg.LinearTol
    k1 = k1 - 1;                                             % walk down while it stays on the line
end
k2 = mid(end);
while k2 < numel(snrs) && abs(dev(k2+1)) <= cfg.LinearTol
    k2 = k2 + 1;                                             % walk up while it stays on the line
end
lo = snrs(k1);
hi = snrs(k2);
end

function val = checkValidation(cfg, cal, rows)
snr = [rows.SNR_dB];
estC = [rows.PFR_IQ_dB] - cal.offsetC;                       % SNR estimate on fresh draws, I/Q
estR = [rows.PFR_I_dB] - cal.offsetR;                        % same for the I channel
pthC = [rows.NoiseMed_IQ_dB] + cal.margin;                   % Pth estimate on fresh draws
inC = snr >= cal.loC & snr <= cal.hiC;                       % only where the estimate is trustworthy
inR = snr >= cal.loR & snr <= cal.hiR;
inM = snr >= cal.loM & snr <= cal.hiM;
val.snrErrC_mean = mean(estC(inC) - snr(inC));               % average SNR error, dB
val.snrErrC_max = max(abs(estC(inC) - snr(inC)));            % worst SNR error, dB
val.snrErrR_mean = mean(estR(inR) - snr(inR));
val.snrErrR_max = max(abs(estR(inR) - snr(inR)));
val.pthErr_mean = mean(pthC(inM) - [rows(inM).PthRuleNoFloor_dB]);   % Pth estimate minus the rule, dB
val.pthErr_max = max(abs(pthC(inM) - [rows(inM).PthRuleNoFloor_dB]));
val.n = numel(rows);
end

function rows = runRealData(cfg, cal, outDir)
rows = struct([]);
gtFile = cfg.GroundTruthFile;
if ~isfile(gtFile)
    logMsg(outDir, 'Ground truth file %s not found, real walks skipped.', gtFile);
    return;
end
gt = readtable(gtFile, 'TextType', 'char');                  % list of walks with participant, pace and direction
nW = min(height(gt), cfg.MaxWalks);
if strcmp(cfg.LoaderName, 'LoadRealRecordingByName_v2')
    loader = @(name) LoadRealRecordingByName_v2(name, cfg.RecordingFolder);
else
    loader = str2func(cfg.LoaderName);
end
for i = 1:nW
    r = struct();
    r.wav_file = gt.wav_file{i};
    r.participant = textCol(gt, 'participant', i);
    r.pace = textCol(gt, 'pace', i);
    r.direction = textCol(gt, 'direction', i);
    r.isLimp = double(contains(lower(textCol(gt, 'notes', i)), 'limp'));   % same limp flag as the batch runners
    try
        [x, fs] = loader(r.wav_file);                        % cropped recording if there is one, same as the batch runners
        m = measureSignal(x, fs, cfg);
        r.PFR_dB = m.PFR_dB;
        r.SNRequiv_dB = m.PFR_dB - cal.offsetR;              % I channel calibration, the recordings are one channel
        r.SNRflag = rangeFlag(r.SNRequiv_dB, cal.loR, cal.hiR);   % says if it is outside the trustworthy range
        r.NoiseMed_dB = m.NoiseMed_dB;
        r.Pth_est_dB = m.NoiseMed_dB + cal.margin;           % threshold from the noise level, the same margin as the simulator rule
        r.Pth_estFloor_dB = max(r.Pth_est_dB, cfg.PthCap);   % with the simulator rule's -40 dB floor
        r.Cross_PthEst_pct = 100*mean(m.FrameMaxNoise_dB > r.Pth_est_dB);   % frames where noise crosses the estimated Pth
        r.Cross_Pth40_pct = 100*mean(m.FrameMaxNoise_dB > cfg.PthCap);      % frames where noise crosses -40 dB
        r.Cross_Pth45_pct = 100*mean(m.FrameMaxNoise_dB > cfg.PthFixed);    % frames where noise crosses the fixed -45 dB
        r.errorMsg = '';
    catch err
        r.PFR_dB = NaN; r.SNRequiv_dB = NaN; r.SNRflag = 'error';
        r.NoiseMed_dB = NaN; r.Pth_est_dB = NaN; r.Pth_estFloor_dB = NaN;
        r.Cross_PthEst_pct = NaN; r.Cross_Pth40_pct = NaN; r.Cross_Pth45_pct = NaN;
        r.errorMsg = err.message;                            % keep going, note what went wrong
    end
    rows = appendRow(rows, r);
    if mod(i, 10) == 0 || i == nW
        logMsg(outDir, 'Real walks: %d of %d done.', i, nW);
    end
end
end

function s = textCol(T, name, i)
s = '';
if any(strcmp(T.Properties.VariableNames, name))
    v = T.(name);
    if iscell(v) && ischar(v{i})
        s = strtrim(v{i});
    end
end
end

function f = rangeFlag(v, lo, hi)
if isnan(v)
    f = 'error';
elseif v < lo
    f = sprintf('below %g dB, at most', lo);                 % under the straight-line part, the real SNR is this or lower
elseif v > hi
    f = sprintf('above %g dB, at least', hi);                % over the straight-line part, the real SNR is this or higher
else
    f = 'in range';
end
end

function selfTest(cfg, simCopy, outDir)
L = {};
D = simulateSpeed(simCopy, cfg, 1.2);                        % one simulated walk
Nw = 2^nextpow2(round(cfg.StepCPIreq*D.fs));
rng(cfg.SeedBase);
z = randn(1, numel(D.clean)) + 1i*randn(1, numel(D.clean));
x = addNoise(D.clean, Nw, 30, z);
ratio = mean(abs(x - D.clean).^2) / (Nw*10^(-30/10));       % injected noise power compared with the formula
L{end+1} = sprintf('1. Noise power at 30 dB matches Nw*10^(-SNR/10): ratio %.3f (%s)', ratio, passText(abs(ratio - 1) < 0.05));
m1 = measureSignal(x, D.fs, cfg);
m2 = measureSignal(x, D.fs, cfg);
L{end+1} = sprintf('2. Same input gives the same output: %s', passText(m1.PFR_dB == m2.PFR_dB && m1.NoiseMed_dB == m2.NoiseMed_dB));
x40 = addNoise(D.clean, Nw, 40, z);
m3 = measureSignal(x40, D.fs, cfg);
step = m3.PFR_dB - m1.PFR_dB;                                % 10 dB less noise should raise PFR by about 10 dB
L{end+1} = sprintf('3. PFR rises %.2f dB from 30 to 40 dB SNR, expected about 10 (%s)', step, passText(abs(step - 10) < 2));
mR = measureSignal(real(x), D.fs, cfg);
L{end+1} = sprintf('4. I channel measurement is finite: %s (PFR %.1f dB, I/Q %.1f dB)', passText(isfinite(mR.PFR_dB) && isfinite(mR.NoiseMed_dB)), mR.PFR_dB, m1.PFR_dB);
hasLoader = strcmp(cfg.LoaderName, 'LoadRealRecordingByName_v2') || exist(cfg.LoaderName, 'file') ~= 0;
hasGT = isfile(cfg.GroundTruthFile);
L{end+1} = sprintf('5. Real-data loader found: %s. Ground truth file found: %s', yesNo(hasLoader), yesNo(hasGT));
fid = fopen(fullfile(outDir, 'Self_Test.txt'), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
for k = 1:numel(L)
    logMsg(outDir, 'Self test %s', L{k});
end
end

function s = passText(ok)
if ok
    s = 'PASS';
else
    s = 'FAIL';
end
end

function s = yesNo(ok)
if ok
    s = 'yes';
else
    s = 'no';
end
end

function rows = appendRow(rows, r)
if isempty(rows)
    rows = r;                                                % first row sets the fields
else
    rows(end+1) = r;
end
end

function writeTable(filename, rows)
if isempty(rows)
    return;
end
f = fieldnames(rows);
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', strjoin(f', ','));                      % header line
for i = 1:numel(rows)
    parts = cell(1, numel(f));
    for j = 1:numel(f)
        v = rows(i).(f{j});
        if ischar(v)
            parts{j} = ['"' strrep(v, '"', '''') '"'];       % text in quotes
        else
            parts{j} = sprintf('%.6g', v);                   % numbers to six figures
        end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
fclose(fid);
end

function writeSummary(filename, cal)
fid = fopen(filename, 'w');
fprintf(fid, 'SNR_dB,meanPFR_IQ_dB,sdPFR_IQ_dB,devFromLine_IQ_dB,meanPFR_I_dB,sdPFR_I_dB,devFromLine_I_dB,meanNoiseMed_IQ_dB,meanNoiseMed_I_dB,ruleMargin_dB\n');
for k = 1:numel(cal.SNRs)
    fprintf(fid, '%g,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n', cal.SNRs(k), cal.meanPFR_C(k), cal.sdPFR_C(k), cal.devC(k), ...
        cal.meanPFR_R(k), cal.sdPFR_R(k), cal.devR(k), cal.meanNoise_C(k), cal.meanNoise_R(k), cal.marginAt(k));
end
fprintf(fid, 'fit,offset_IQ=%.3f,offset_I=%.3f,margin=%.3f,range_IQ=%g-%g,range_I=%g-%g,range_margin=%g-%g\n', ...
    cal.offsetC, cal.offsetR, cal.margin, cal.loC, cal.hiC, cal.loR, cal.hiR, cal.loM, cal.hiM);
fclose(fid);
end

function fig = newFigure()
fig = figure('Visible', 'off', 'Color', 'w');                % white figure, not shown while it runs
try
    theme(fig, 'light');                                     % stops MATLAB dark mode turning it black
catch
end
end

function finishAxes(ax, xText, yText, titleText)
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');         % white plot area, black axes
xlabel(ax, xText);
ylabel(ax, yText);
title(ax, titleText);
grid(ax, 'on');
end

function saveFigure(fig, basePath)
try
    exportgraphics(fig, [basePath '.png'], 'BackgroundColor', 'white', 'Resolution', 200);   % PNG for the report draft
    exportgraphics(fig, [basePath '.pdf'], 'BackgroundColor', 'white', 'ContentType', 'vector');   % PDF for LaTeX
catch
    saveas(fig, [basePath '.png']);                          % older MATLAB fallback
    print(fig, [basePath '.pdf'], '-dpdf');
end
close(fig);
end

function makePlots(cfg, cal, calRows, valRows, realRows, outDir)
fig = newFigure();
ax = axes(fig);
hold(ax, 'on');
errorbar(ax, cal.SNRs, cal.meanPFR_C, cal.sdPFR_C, 'o-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
errorbar(ax, cal.SNRs, cal.meanPFR_R, cal.sdPFR_R, 's-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
plot(ax, cal.SNRs, cal.SNRs + cal.offsetC, '--', 'Color', [0 0.45 0.74]);   % straight line through the I/Q fit
plot(ax, cal.SNRs, cal.SNRs + cal.offsetR, '--', 'Color', [0.85 0.33 0.10]);   % straight line through the I channel fit
legend(ax, {'I/Q (Chapter 4)', 'I channel only (like IPM-365)', 'I/Q straight line', 'I channel straight line'}, 'Location', 'northwest');
finishAxes(ax, 'Simulator SNR setting (dB)', 'Peak-to-floor ratio PFR (dB)', ...
    sprintf('SNR calibration, simulated walks, %d speeds, mean \\pm 1 SD', numel(cfg.Speeds)));
saveFigure(fig, fullfile(outDir, 'Fig_SNR_Calibration'));

fig = newFigure();
ax = axes(fig);
hold(ax, 'on');
snr = [valRows.SNR_dB];
pthEst = [valRows.NoiseMed_IQ_dB] + cal.margin;             % estimate from the noise level
jitter = 0.6*(mod(0:numel(snr)-1, 5)/4 - 0.5);              % small fixed sideways spread so the dots don't sit on top of each other
plot(ax, snr + jitter, pthEst, '.', 'Color', [0 0.45 0.74], 'MarkerSize', 12);
plot(ax, cal.SNRs, 10 - cal.SNRs, 'k-', 'LineWidth', 1.5);   % the rule without the floor
plot(ax, cal.SNRs, max(10 - cal.SNRs, cfg.PthCap), 'k--', 'LineWidth', 1);   % the rule with the -40 dB floor
legend(ax, {'P_{th} estimated from the noise level', 'Rule 10 - SNR', 'Rule with -40 dB floor'}, 'Location', 'northeast');
finishAxes(ax, 'Simulator SNR setting (dB)', 'P_{th} (dB below the maximum)', 'Estimated threshold against the simulator rule (fresh draws)');
saveFigure(fig, fullfile(outDir, 'Fig_Pth_Validation'));

if isempty(realRows)
    return;
end
ok = ~isnan([realRows.SNRequiv_dB]);
fig = newFigure();
ax = axes(fig);
hist(ax, [realRows(ok).SNRequiv_dB], 20);                    % how clean the real recordings are
finishAxes(ax, 'Estimated SNR, equivalent simulator setting (dB)', 'Number of walks', ...
    sprintf('Estimated SNR of the real walks (n = %d)', sum(ok)));
saveFigure(fig, fullfile(outDir, 'Fig_RealData_SNR'));

fig = newFigure();
ax = axes(fig);
hold(ax, 'on');
hist(ax, [realRows(ok).Pth_est_dB], 20);
yl = ylim(ax);
plot(ax, [cfg.PthFixed cfg.PthFixed], yl, 'r-', 'LineWidth', 1.5);   % the fixed threshold used now
plot(ax, [cfg.PthCap cfg.PthCap], yl, 'k--', 'LineWidth', 1.5);    % the simulator rule's floor
legend(ax, {'Walks', sprintf('Fixed P_{th} = %d dB', cfg.PthFixed), sprintf('Floor %d dB', cfg.PthCap)}, 'Location', 'northwest');
finishAxes(ax, 'Estimated P_{th} (dB below the maximum)', 'Number of walks', 'Threshold estimated from each real walk''s noise level');
saveFigure(fig, fullfile(outDir, 'Fig_RealData_Pth'));
end

function writeReport(cfg, cal, val, realRows, outDir, tRun)
L = {};
L{end+1} = 'SNR AND PTH ESTIMATE REPORT (written automatically)';
L{end+1} = sprintf('Run folder: %s', outDir);
L{end+1} = sprintf('Run time %.1f min. Quick mode: %d.', tRun/60, cfg.Quick);
L{end+1} = sprintf('Simulator: %d speeds, SNR %s dB, %d calibration and %d validation draws each, R0 = %g m, Chen model %s, fixed seed %d.', ...
    numel(cfg.Speeds), mat2str(cfg.SNRlist), cfg.nCal, cfg.nVal, cfg.R0_m, cfg.AmplitudeMode, cfg.SeedBase);
L{end+1} = sprintf('Processing: step-detection spectrogram, CPI request %.3f s (window %.3f s at 44.1 kHz), noise region %.0f to %.0f Hz.', ...
    cfg.StepCPIreq, 2^nextpow2(round(cfg.StepCPIreq*44100))/44100, cfg.NoiseBand(1)*350*24/5.8, cfg.NoiseBand(2)*350*24/5.8);
L{end+1} = '';
L{end+1} = '1. CALIBRATION ON THE SIMULATOR';
L{end+1} = sprintf('  PFR = SNR + offset. Offset I/Q %.2f dB, I channel only %.2f dB (fitted on %g to %g dB).', cal.offsetC, cal.offsetR, cfg.FitRange(1), cfg.FitRange(2));
L{end+1} = sprintf('  The I channel offset is %.2f dB lower than I/Q: one channel carries half the signal but its noise is spread over both sides of 0 Hz.', cal.offsetC - cal.offsetR);
L{end+1} = sprintf('  Trustworthy SNR range (within %g dB of the straight line): I/Q %g to %g dB, I channel %g to %g dB.', cfg.LinearTol, cal.loC, cal.hiC, cal.loR, cal.hiR);
L{end+1} = '  Below that range the strongest noise cell is as strong as the signal, so the estimate cannot go lower.';
L{end+1} = '  Above it the noise region is no longer pure noise (window sidelobes and weak body returns), so the estimate stops rising.';
L{end+1} = sprintf('  The rule Pth = 10 - SNR sits %.2f dB above the median noise level (flat from %g to %g dB).', cal.margin, cal.loM, cal.hiM);
L{end+1} = '';
L{end+1} = '  SNR (dB) | PFR I/Q mean (SD) | off line | PFR I mean (SD) | off line | rule margin';
for k = 1:numel(cal.SNRs)
    L{end+1} = sprintf('  %7g  | %7.2f (%5.2f)     | %6.2f   | %7.2f (%5.2f)   | %6.2f   | %6.2f', cal.SNRs(k), cal.meanPFR_C(k), cal.sdPFR_C(k), cal.devC(k), ...
        cal.meanPFR_R(k), cal.sdPFR_R(k), cal.devR(k), cal.marginAt(k)); %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = '2. CHECK ON FRESH SIMULATED DRAWS (inside the trustworthy range only)';
L{end+1} = sprintf('  SNR estimate error, I/Q: mean %+.2f dB, worst %.2f dB. I channel: mean %+.2f dB, worst %.2f dB.', val.snrErrC_mean, val.snrErrC_max, val.snrErrR_mean, val.snrErrR_max);
L{end+1} = sprintf('  Pth estimate minus the rule 10 - SNR: mean %+.2f dB, worst %.2f dB.', val.pthErr_mean, val.pthErr_max);
L{end+1} = '';
L{end+1} = '3. REAL WALKS';
if isempty(realRows)
    L{end+1} = '  Not run.';
else
    ok = ~isnan([realRows.SNRequiv_dB]);
    nErr = sum(~ok);
    snrv = [realRows(ok).SNRequiv_dB];
    pth = [realRows(ok).Pth_est_dB];
    flags = {realRows(ok).SNRflag};
    L{end+1} = sprintf('  %d walks measured, %d failed to load or process.', sum(ok), nErr);
    if any(ok)
        L{end+1} = sprintf('  Estimated SNR (equivalent simulator setting, I channel calibration): median %.1f dB, middle half %.1f to %.1f dB, range %.1f to %.1f dB.', ...
            median(snrv), prctileSimple(snrv, 25), prctileSimple(snrv, 75), min(snrv), max(snrv));
        L{end+1} = sprintf('  Walks inside the trustworthy range: %d, below it: %d, above it: %d.', sum(strcmp(flags, 'in range')), ...
            sum(strncmp(flags, 'below', 5)), sum(strncmp(flags, 'above', 5)));
        L{end+1} = sprintf('  Estimated Pth: median %.1f dB, middle half %.1f to %.1f dB, range %.1f to %.1f dB.', ...
            median(pth), prctileSimple(pth, 25), prctileSimple(pth, 75), min(pth), max(pth));
        L{end+1} = sprintf('  Walks whose estimated Pth is below the fixed %d dB: %d. Above it: %d.', cfg.PthFixed, sum(pth < cfg.PthFixed), sum(pth >= cfg.PthFixed));
        L{end+1} = sprintf('  Frames where the noise band crosses the threshold, median over walks: estimated Pth %.1f %%, %d dB %.1f %%, %d dB %.1f %%.', ...
            median([realRows(ok).Cross_PthEst_pct]), cfg.PthCap, median([realRows(ok).Cross_Pth40_pct]), cfg.PthFixed, median([realRows(ok).Cross_Pth45_pct]));
        L = addGroups(L, realRows(ok), 'direction', cfg.PthFixed);
        L = addGroups(L, realRows(ok), 'pace', cfg.PthFixed);
        L = addGroups(L, realRows(ok), 'participant', cfg.PthFixed);
    end
end
L{end+1} = '';
L{end+1} = 'HOW TO READ THIS';
L{end+1} = '  Higher SNR = less noise. Every 10 dB is ten times less noise power. 0 dB means the signal and the noise in a spectrogram cell are equally strong.';
L{end+1} = '  The estimated SNR of a real walk is the simulator setting that would give the same peak-to-floor ratio on a one-channel recording.';
L{end+1} = '  On real data the strongest cell can be clutter near 0 Hz rather than the torso. Pth is measured from the same maximum, so the estimated Pth is still on the right scale.';
L{end+1} = '  A walk flagged above the range is at least that clean. A walk flagged below the range is that noisy or worse.';
L{end+1} = '  The crossing check counts frames where the loudest cell in the upper half of the leg band is above Pth. At high SNR this also counts real fast leg returns, so read it as a comparison between thresholds, not as a pure false-alarm rate.';
L{end+1} = '  The calibration offset changes a little with walking speed (see the SD columns), so an SNR estimate is good to about the worst error in Section 2.';
L{end+1} = '  Compare the estimated Pth with the fixed -45 dB. If they agree, the fixed value is backed by the noise level. If not, the walks that differ are the ones to look at.';
L{end+1} = '  Simulated calibration only. The real walks are measured, not scored, so nothing here uses the video ground truth.';
fid = fopen(fullfile(outDir, 'SNREstimate_Report.txt'), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
end

function L = addGroups(L, rows, field, pthFixed)
vals = {rows.(field)};
u = unique(vals(~cellfun(@isempty, vals)));
if isempty(u)
    return;
end
L{end+1} = sprintf('  By %s:', field);
for k = 1:numel(u)
    sel = strcmp(vals, u{k});
    L{end+1} = sprintf('    %-14s n = %3d   SNR median %.1f dB   Pth median %.1f dB   noise crossing at %d dB %.1f %%', ...
        u{k}, sum(sel), median([rows(sel).SNRequiv_dB]), median([rows(sel).Pth_est_dB]), pthFixed, median([rows(sel).Cross_Pth45_pct])); %#ok<AGROW>
end
end

function p = prctileSimple(x, q)
x = sort(x(:));
if isempty(x)
    p = NaN;
    return;
end
pos = 1 + (numel(x) - 1)*q/100;                              % linear interpolation between sorted values
lo = floor(pos);
hi = ceil(pos);
p = x(lo) + (pos - lo)*(x(hi) - x(lo));
end

function [RdrRxSignal, fs, fc, c, resolvedPath] = LoadRealRecordingByName_v2(wavFile_CW, recordingFolder)

c = 299792458;    % speed of light, m/s
fc = 24e9;         % radar centre frequency, matches the cantenna hardware

resolvedPath = ResolveRealRecordingPath_v1(wavFile_CW, recordingFolder);

[Y, fs] = audioread(resolvedPath, 'native');   % read the raw WAV samples at native (integer) precision
RdrRxSignal = -double(Y(:,1));
RdrRxSignal = RdrRxSignal(:).';
fprintf('Loaded %s: %d samples at fs = %d Hz (%.4f s long)\n', resolvedPath, length(RdrRxSignal), fs, length(RdrRxSignal)/fs);
end

function resolvedPath = ResolveRealRecordingPath_v1(wavFile, recordingFolder)
if nargin < 2 || isempty(recordingFolder)
    recordingFolder = fileparts(mfilename('fullpath'));
end
[~, stem] = fileparts(char(wavFile));
tok = regexp(stem, '^(P\d+)_', 'tokens', 'once');
assert(~isempty(tok), 'EEE4022S:RecordingName', 'Unrecognised trial filename: %s', char(wavFile));
participant = tok{1};
candidates = {
    fullfile(recordingFolder, 'CroppedRecordings', [stem '_cropped.wav']);
    fullfile(recordingFolder, [participant ' Test Walks'], [stem '_cropped.wav']);
    fullfile(recordingFolder, [participant ' Test Walks'], [stem '.wav']);
    fullfile(recordingFolder, [stem '_cropped.wav']);
    fullfile(recordingFolder, [stem '.wav'])};
for k = 1:numel(candidates)
    if isfile(candidates{k})
        resolvedPath = candidates{k};
        return;
    end
end
error('EEE4022S:RecordingMissing', 'No original or cropped WAV found for %s.', char(wavFile));
end
