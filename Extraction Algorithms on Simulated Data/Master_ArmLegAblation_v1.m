%% Arm and leg return comparison
% Zahraa Sattar | EEE4022S
%
% Test how limb amplitudes and shoulder swing affect torso-speed bias.
% This tests a possible mechanism in the model, rather than identifying it in the recordings.

function Master_ArmLegAblation_v1(userCfg)

cfg.Speeds = [0.8 1.0 1.2 1.4];                              % m/s, the range of the measured walks
cfg.SNR = 70;                                                % dB, the median of the measured walks
cfg.nTrials = 4;                                             % noise draws per case and speed
cfg.R0_m = 10;                                               % same start range as Chapter 4
cfg.AmplitudeMode = 'A';                                     % Chen amplitude model, as in Chapter 4
cfg.Seed = 20267001;                                         % simulator seed, same walks as round 4
cfg.SeedNoise = 20269101;                                    % new noise seed
base=fileparts(mfilename('fullpath'));
addpath(fullfile(base,'..','Radar Simulator'));
addpath(fullfile(base,'..','Extraction Algorithms on Measured Data'));
cfg.SimName = fullfile(base,'Report Sources','RadarSimulator_Ablation');
cfg.PipelineName = 'RealData_Pipeline_v14';                   % Final extractor with proportional spacing and the edge rule
cfg.SelectionFile = fullfile(base,'Report Results','Ablation_Selection.csv');   % the frozen Final setting
cfg.Quick = false;
if nargin > 0 && isstruct(userCfg)
    fn = fieldnames(userCfg);
    for q = 1:numel(fn)
        cfg.(fn{q}) = userCfg.(fn{q});
    end
end
if cfg.Quick
    cfg.Speeds = 1.0;
    cfg.nTrials = 1;
end
if exist('OCTAVE_VERSION', 'builtin'), pkg load signal; end
warning('off', 'all');
addpath(pwd);

cases = {
    'Baseline',            1, 1, 1
    'No arms',             0, 1, 1
    'No legs',             1, 0, 1
    'Arms x2',             2, 1, 1
    'Arms x4',             4, 1, 1
    'Arms x8',             8, 1, 1
    'Swing x1.5',          1, 1, 1.5
    'Swing x2',            1, 1, 2
    'Arms x4 + swing x1.5', 4, 1, 1.5
    'Arms x4 + swing x2',   4, 1, 2
    'Arms x8 + swing x2',   8, 1, 2
    };

final = readSelection(cfg.SelectionFile);                    % the frozen Final setting
final.name = 'Final';
final.propSpacing = 1;
final.propEdge = 1;
opts = struct('CPI_torso_s', 0.2, 'CPI_legs_s', 0.1, 'Pth_dB', -45);   % same fixed settings as the measured data

outDir = ['ArmLegAblation_' datestr(now, 'yyyymmdd_HHMMSS')];
mkdir(outDir);
logFile = fullfile(outDir, 'Run_Log.txt');
simCopy = makeSimCopy(cfg.SimName, outDir);
Nw = 2^nextpow2(round(0.1 * 44100));                         % SNR defined at the 0.186 s window, as in Chapter 4
lam = 299792458 / 24e9;
rows = [];
curves = [];
vAxis = [];
for k = 1:size(cases, 1)
    for s = 1:numel(cfg.Speeds)
        D = simulate(cfg, cfg.Speeds(s), simCopy, cases{k, 2}, cases{k, 3}, cases{k, 4});
        bias = zeros(1, cfg.nTrials);
        for q = 1:cfg.nTrials
            rng(cfg.SeedNoise + 1000 * k + 100 * s + q, 'twister');
            z = randn(1, numel(D.clean)) + 1i * randn(1, numel(D.clean));
            noisePower = 10^(-(cfg.SNR - 10 * log10(Nw)) / 10);
            x = real(D.clean + sqrt(noisePower / 2) * z);    % one channel, like the IPM-365
            R = feval(cfg.PipelineName, x, D.fs, D.fc, D.c, final, opts);
            bias(q) = 100 * (R(1).estSpeed_mps - D.v0) / D.v0;
            if q == 1
                walkSamples = x(round((D.buffer + 0.4) * D.fs) : round((D.buffer + D.simT - 0.4) * D.fs));
                [vAxis, A] = meanSpectrum(walkSamples, D.fs, lam);
                lo = lowerEdge(vAxis, A) / D.v0;
                curves = [curves; {cases{k, 1}, D.v0, vAxis(:)' / D.v0, A(:)'}]; %#ok<AGROW>
            end
        end
        rows = [rows; {cases{k, 1}, D.v0, mean(bias), std(bias), lo}]; %#ok<AGROW>
        logMsg(logFile, '%-22s %.1f m/s: bias %+.2f%%  edge %.2f', cases{k, 1}, D.v0, mean(bias), lo);
    end
end

fid = fopen(fullfile(outDir, 'ArmLegAblation_Results.csv'), 'w');
fprintf(fid, 'case,speed_mps,bias_pct,bias_sd_pct,lowerEdge_rel\n');
for r = 1:size(rows, 1)
    fprintf(fid, '%s,%.1f,%.3f,%.3f,%.3f\n', rows{r, 1}, rows{r, 2}, rows{r, 3}, rows{r, 4}, rows{r, 5});
end
fclose(fid);
fid = fopen(fullfile(outDir, 'ArmLegAblation_Spectra.csv'), 'w');
fprintf(fid, 'case,speed_mps,velocity_over_v0,meanNormalisedPower\n');
for r = 1:size(curves, 1)
    for j = 1:numel(curves{r, 3})
        fprintf(fid, '%s,%.1f,%.4f,%.5f\n', curves{r, 1}, curves{r, 2}, curves{r, 3}(j), curves{r, 4}(j));
    end
end
fclose(fid);
L = {};
L{end + 1} = 'ARM AND LEG ABLATION: torso-speed bias of the Final extractor in the simulator';
L{end + 1} = sprintf('Speeds %s m/s, SNR %d dB, %d noise draws. Bias = mean (estimate - true) / true. Edge = lower half-power edge / true speed.', mat2str(cfg.Speeds), cfg.SNR, cfg.nTrials);
names = unique(rows(:, 1), 'stable');
for k = 1:numel(names)
    m = strcmp(rows(:, 1), names{k});
    L{end + 1} = sprintf('  %-22s mean bias %+6.2f%%  (per speed %s)   mean edge %.2f', names{k}, mean(cell2mat(rows(m, 3))), ...
        mat2str(round(100 * cell2mat(rows(m, 3))') / 100), mean(cell2mat(rows(m, 5)))); %#ok<AGROW>
end
fid = fopen(fullfile(outDir, 'ArmLegAblation_Report.txt'), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
fprintf('%s\n', L{:});
end

function D = simulate(m_cfg, v0, m_simCopy, AblArm, AblLeg, AblSwing)
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
D.buffer = bufferTime_s;
D.simT = SimulationTime_s;
end

function [v, A] = meanSpectrum(x, fs, lam)
N = 16384;
hop = N / 4;
w = hamming(N);
nF = floor((numel(x) - N) / hop) + 1;
f = (0:N/2-1) * fs / N;
v = f * lam / 2;
sel = v > 0.3 & v < 2.4;
acc = zeros(sum(sel), 1);
for q = 1:nF
    seg = x((q - 1) * hop + (1:N));
    seg = seg(:) - mean(seg);
    P = abs(fft(seg .* w)).^2;
    P = P(1:N/2);
    P = P(sel);
    acc = acc + P(:) / max(P);
end
v = v(sel);
A = acc / nF;
A = conv(A, ones(3, 1) / 3, 'same');
A = A / max(A);
end

function e = lowerEdge(v, A)
j = find(A >= 0.5, 1, 'first');
e = v(j);
end

function s = readSelection(file)
fid = fopen(file, 'r');
fgetl(fid);
s = struct();
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    p = strsplit(ln, ',');
    v = str2double(p{2});
    if isnan(v) && ~strcmpi(p{2}, 'NaN')
        s.(p{1}) = p{2};
    else
        s.(p{1}) = v;
    end
end
fclose(fid);
end

function dest = makeSimCopy(simName, outDir)
src = [simName '.m'];
assert(exist(src, 'file') == 2, 'Simulator %s not found in the Current Folder.', src);
txt = fileread(src);
txt = regexprep(txt, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*(?:\r?\n|$)', '');
txt = strrep(txt, 'As = 9.88*Vr;', 'As = 9.88*Vr*AblSwing;');
gainCode = sprintf(['%%%% Ablation gains\nfor kk = [6 7 11 12], AppliedAmplitude{kk} = AppliedAmplitude{kk}*AblArm; end\n' ...
    'for kk = [3 4 5 8 9 10], AppliedAmplitude{kk} = AppliedAmplitude{kk}*AblLeg; end\n']);
txt = strrep(txt, '%% Amplitude variables for compatibility', [gainCode '%% Amplitude variables for compatibility']);
dest = fullfile(pwd, outDir, 'WC_Simulator.m');
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
