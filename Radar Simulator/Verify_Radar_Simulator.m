%% Verify the radar simulator
% Zahraa Sattar | EEE4022S
%
% Check body-segment Doppler, torso phase and added-noise behaviour.
% Save the verification tables and figures in the Results folder.

function Verify_Radar_Simulator

folder = fileparts(mfilename('fullpath'));
simFile = fullfile(folder, ...
    'RadarSimulator_v11_torso_legs_arms_head_buffered_noplots.m');
assert(exist(simFile,'file') ~= 0, 'Could not find RadarSimulator_v11.');

raw = fileread(simFile);
raw = regexprep(raw, '(?m)^[ \t]*(?:clear[ \t]+all|close[ \t]+all)[^\r\n]*(?:\r?\n|$)', '');
tempFolder = tempname;
mkdir(tempFolder);
tempCleanup = onCleanup(@() rmdir(tempFolder, 's'));
simFile = fullfile(tempFolder, 'Verification_Simulator.m');
fid = fopen(simFile, 'w');
fwrite(fid, raw);
fclose(fid);


oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(folder);

out = fullfile(fileparts(folder), 'Results', ...
    ['Verification_AllTests_V11_' datestr(now,'yyyymmdd_HHMMSS')]);
mkdir(out);

R0 = 10;
AmpMode = "A";

regimes = struct('name',{'a','b','c'},'speed',{0.3,1.0,1.6});
segTypes = ["torso","head","thigh","shank","foot","uarm","farm"];
weights = zeros(3,7);
for k = 1:3
    d = runV11(simFile,regimes(k).speed,R0,AmpMode,[]);
    weights(k,:) = d.summary.AppliedWeightMean(1:7).';
end
weightTable = array2table(weights,'VariableNames',cellstr(segTypes));
weightTable = addvars(weightTable,[regimes.speed].', ...
    'Before',1,'NewVariableNames','Speed_ms');
disp(weightTable);
writetable(weightTable,fullfile(out,'ModeA_WeightsBySpeed.csv'));

names = {'Torso','Head','Left upper arm','Left forearm', ...
         'Right upper arm','Right forearm','Left thigh','Left shank', ...
         'Left foot','Right thigh','Right shank','Right foot'};
rows = cell(0,9);

for k = 1:3
    d = runV11(simFile,regimes(k).speed,R0,AmpMode,[]);
    assert(d.regime == regimes(k).name, ...
        'The simulator chose a different gait regime from the one requested.');

    lambda = d.lambda;
    win = hamming(d.Nw);
    noverlap = floor(0.8*d.Nw);

    S = cell(12,1);
    for i = 1:12
        [S{i},F,T] = Zahraa_Spectrogram(d.clean(i,:),win,noverlap,d.Nw,d.fs);
    end
    [Scomposite,~,~] = Zahraa_Spectrogram(sum(d.clean,1),win,noverlap,d.Nw,d.fs);

    velocityAxis = F(:)*lambda/2;
    T = T(:).';
    binWidth = abs(mean(diff(velocityAxis)));

    halfWindow = d.Nw/(2*d.fs);
    inside = T-halfWindow >= d.buffer & T+halfWindow <= d.buffer+d.duration;
    frameTime = T(inside)-d.buffer;

    denseTime = linspace(0,d.duration,6001);
    [referenceDense,referenceRegime] = chenReference( ...
        denseTime,d.v0,d.height,d.R0,d.radar);
    assert(referenceRegime == d.regime, ...
        'Reference and simulator disagree on the gait regime.');
    referenceFrames = chenReference(frameTime,d.v0,d.height,d.R0,d.radar);

    velocityLimits = [floor(min(referenceDense(:))-1), ...
                      ceil(max(referenceDense(:))+1)];
    visible = velocityAxis >= velocityLimits(1) & ...
              velocityAxis <= velocityLimits(2);

    fig = lightFigure([100 100 1400 800]);
    layout = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');

    for i = 1:12
        ax = nexttile(layout);
        bandPower = abs(S{i}(visible,inside)).^2;
        bandPower = bandPower/max(bandPower(:));
        imagesc(ax,frameTime,velocityAxis(visible), ...
            10*log10(max(bandPower,1e-8)),[-50 0]);
        axis(ax,'xy');
        hold(ax,'on');
        colormap(ax,jet(256));
        referenceLine = plot(ax,denseTime,referenceDense(i,:), ...
            'w-','LineWidth',1.2);

        ridge = ridgeVelocity(S{i},inside,velocityAxis);
        ridgeError = ridge-referenceFrames(i,:);
        rmse = sqrt(mean(ridgeError.^2));
        maxError = max(abs(ridgeError));

        rows(end+1,:) = {d.regime,d.v0,names{i},rmse,maxError, ...
            rmse/binWidth,binWidth,1000*d.Nw/d.fs,rmse/binWidth < 0.5}; %#ok<AGROW>

        xlim(ax,[frameTime(1) frameTime(end)]);
        ylim(ax,velocityLimits);
        title(ax,sprintf('(%c) %s, RMSE %.3f m/s',char('a'+i-1),names{i},rmse), ...
            'FontWeight','normal','Color','k','FontSize',10);
        styleAxes(ax,10);
        grid(ax,'off');
    end

    spectrogramKey = patch(ax,NaN,NaN,[0.7 0 0],'EdgeColor','none');
    keyBox = legend(ax,[spectrogramKey referenceLine], ...
        {'Simulated spectrogram','Reference from motion equations'}, ...
        'Orientation','horizontal','Color',[0.15 0.15 0.15], ...
        'TextColor','w','EdgeColor',[0.15 0.15 0.15]);
    keyBox.Layout.Tile = 'north';

    colourBar = colorbar(ax);
    colourBar.Layout.Tile = 'east';
    colourBar.Label.String = 'Normalised power (dB)';
    colourBar.Color = 'k';

    xText = xlabel(layout,'Walking time (s)');
    yText = ylabel(layout,'Radial velocity (m/s)');
    set([xText yText],'Color','k','FontName','Times New Roman','FontSize',12);

    exportgraphics(fig,fullfile(out, ...
        sprintf('Verification_SegmentOverlays_Regime%s.png',d.regime)), ...
        'BackgroundColor','white','Resolution',300);
    close(fig);

    fig2 = lightFigure([100 100 950 520]);
    ax2 = axes(fig2);
    compositePower = abs(Scomposite(visible,inside)).^2;
    compositePower = compositePower/max(compositePower(:));
    imagesc(ax2,frameTime,velocityAxis(visible), ...
        10*log10(max(compositePower,1e-8)),[-50 0]);
    axis(ax2,'xy');
    hold(ax2,'on');
    colormap(ax2,jet(256));
    allLines = plot(ax2,denseTime,referenceDense.','w-','LineWidth',0.8);
    xlim(ax2,[frameTime(1) frameTime(end)]);
    ylim(ax2,velocityLimits);
    xlabel(ax2,'Walking time (s)');
    ylabel(ax2,'Radial velocity (m/s)');
    title(ax2,sprintf('Composite of all 12 segments, regime %s (v_0 = %.1f m/s)', ...
        d.regime,d.v0),'FontWeight','normal','Color','k');
    styleAxes(ax2,11);
    grid(ax2,'off');
    key2 = patch(ax2,NaN,NaN,[0.7 0 0],'EdgeColor','none');
    legend(ax2,[key2 allLines(1)], ...
        {'Simulated spectrogram','Reference (all 12 segments)'}, ...
        'Location','northoutside','Orientation','horizontal', ...
        'Color',[0.15 0.15 0.15],'TextColor','w','EdgeColor',[0.15 0.15 0.15]);
    colourBar2 = colorbar(ax2);
    colourBar2.Label.String = 'Normalised power (dB)';
    colourBar2.Color = 'k';
    exportgraphics(fig2,fullfile(out, ...
        sprintf('Verification_CompositeOverlay_Regime%s.png',d.regime)), ...
        'BackgroundColor','white','Resolution',300);
    close(fig2);
end

results1 = cell2table(rows,'VariableNames', ...
    {'Regime','Speed_ms','Segment','RMSE_ms','MaxError_ms','RMSE_bins', ...
     'VelocityBin_ms','Window_ms','Pass'});
disp(results1);
writetable(results1,fullfile(out,'Test1_SegmentOverlayErrors.csv'));
for k = 1:3
    r = results1(strcmp(results1.Regime,regimes(k).name),:);
    fprintf(['Regime %s (v0 = %.1f m/s): worst RMSE = %.4f m/s (%.3f bins), ' ...
        'largest frame error %.3f m/s, %d of %d pass\n'], ...
        regimes(k).name,regimes(k).speed,max(r.RMSE_ms),max(r.RMSE_bins), ...
        max(r.MaxError_ms),sum(r.Pass),height(r));
end

speeds = (0.4:0.2:2.0).';
positionError = zeros(size(speeds));
rangeError = zeros(size(speeds));
velocityError = zeros(size(speeds));
aspectGapMax = zeros(size(speeds));
aspectGapMean = zeros(size(speeds));

for n = 1:numel(speeds)
    d = runV11(simFile,speeds(n),R0,AmpMode,[]);
    lambda = d.lambda;

    walkingTime = d.t-d.buffer;
    inside = walkingTime > 0 & walkingTime < d.duration;

    [xRef,zRef] = torsoReference(walkingTime,speeds(n),d.height,d.R0);
    rangeRef = hypot(xRef-d.radar(2),zRef-d.radar(1));

    positionError(n) = max(hypot(d.torsoX(inside)-xRef(inside), ...
                                 d.torsoZ(inside)-zRef(inside)));
    rangeError(n) = max(abs(d.torsoR(inside)-rangeRef(inside)));

    pairs = inside(1:end-1) & inside(2:end);
    rangeStep = diff(rangeRef);
    assert(all(abs(4*pi*rangeStep(pairs)/lambda) < pi), ...
        'Phase increments exceed the unambiguous sampling limit.');

    phaseStep = angle(d.torsoSignal(2:end).*conj(d.torsoSignal(1:end-1)));
    fromPhase = lambda*d.fs/(4*pi)*phaseStep;
    fromMotion = -d.fs*rangeStep;
    velocityError(n) = max(abs(fromPhase(pairs)-fromMotion(pairs)));

    forwardStep = -d.fs*diff(xRef);
    gap = abs(fromMotion(pairs)-forwardStep(pairs));
    aspectGapMax(n) = 100*max(gap)/speeds(n);
    aspectGapMean(n) = 100*mean(gap)/speeds(n);

    if abs(speeds(n)-1.0) < 1e-9
        midTime = (walkingTime(1:end-1)+walkingTime(2:end))/2;
        exampleTime = midTime(pairs);
        exampleMotion = fromMotion(pairs);
        examplePhase = fromPhase(pairs);
    end
end

results2 = table(speeds,positionError,rangeError,velocityError, ...
    aspectGapMax,aspectGapMean,'VariableNames',{'Speed_ms', ...
    'MaxPositionError_m','MaxRangeError_m','MaxPhaseVelocityError_ms', ...
    'MaxRadialForwardGap_pct','MeanRadialForwardGap_pct'});
disp(results2);
writetable(results2,fullfile(out,'Test2_TorsoResults.csv'));

fig3 = lightFigure([100 100 1500 430]);
layout3 = tiledlayout(fig3,1,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout3);
plot(ax,exampleTime,exampleMotion,'k-','LineWidth',1.2);
hold(ax,'on');
markers = 1:1000:numel(exampleTime);
plot(ax,exampleTime(markers),examplePhase(markers),'ro','MarkerSize',4);
xlabel(ax,'Walking time (s)');
ylabel(ax,'Torso radial velocity (m/s)');
title(ax,'(a) Torso at 1.0 m/s','FontWeight','normal','Color','k');
lightLegend(legend(ax,{'From motion equations','From signal phase'},'Location','best'));
styleAxes(ax,11);

ax = nexttile(layout3);
semilogy(ax,speeds,positionError,'bo-',speeds,rangeError,'ms-','LineWidth',1.2);
xlabel(ax,'Prescribed walking speed (m/s)');
ylabel(ax,'Largest difference (m)');
title(ax,'(b) Simulator against motion equations','FontWeight','normal','Color','k');
lightLegend(legend(ax,{'Position','Range'},'Location','best'));
styleAxes(ax,11);

ax = nexttile(layout3);
plot(ax,speeds,velocityError,'bo-','LineWidth',1.2);
xlabel(ax,'Prescribed walking speed (m/s)');
ylabel(ax,'Largest difference (m/s)');
title(ax,'(c) Phase velocity against range velocity','FontWeight','normal','Color','k');
styleAxes(ax,11);

exportgraphics(fig3,fullfile(out,'Verification_Test2_TorsoPhase.png'), ...
    'BackgroundColor','white','Resolution',300);
close(fig3);

fig4 = lightFigure([100 100 700 450]);
ax = axes(fig4);
plot(ax,speeds,aspectGapMax,'mo-',speeds,aspectGapMean,'k^-','LineWidth',1.2);
xlabel(ax,'Prescribed walking speed (m/s)');
ylabel(ax,'|radial - forward| / v_0 (%)');
title(ax,sprintf('Viewing-angle effect on the torso, R_0 = %g m',R0), ...
    'FontWeight','normal','Color','k');
lightLegend(legend(ax,{'Largest over the walk','Mean over the walk'},'Location','best'));
styleAxes(ax,11);
exportgraphics(fig4,fullfile(out,'Verification_ViewingAngleGap.png'), ...
    'BackgroundColor','white','Resolution',300);
close(fig4);

gammas = (10:10:60).';
expectedPower = zeros(size(gammas));
measuredPower = zeros(size(gammas));
expectedRaw = zeros(size(gammas));
measuredRaw = zeros(size(gammas));
builtInRaw = zeros(size(gammas));

for k = 1:numel(gammas)
    d = runV11(simFile,1.0,R0,AmpMode,gammas(k));
    walking = d.t > d.buffer & d.t < d.buffer+d.duration;
    noise = d.received-d.composite;

    expectedPower(k) = d.Nw*10^(-gammas(k)/10);
    measuredPower(k) = mean(abs(noise(walking)).^2);
    signalPower = mean(abs(d.composite(walking)).^2);
    expectedRaw(k) = 10*log10(signalPower/expectedPower(k));
    measuredRaw(k) = 10*log10(signalPower/measuredPower(k));
    builtInRaw(k) = d.builtInSNR;   % v11's own value, whole record incl. buffers
end
differencePct = 100*(measuredPower-expectedPower)./expectedPower;

results3 = table(gammas,expectedPower,measuredPower,differencePct, ...
    expectedRaw,measuredRaw,builtInRaw,'VariableNames',{'Gamma_dB', ...
    'ExpectedNoisePower','MeasuredNoisePower','Difference_pct', ...
    'ExpectedRawSNR_dB','MeasuredRawSNR_dB','V11BuiltInRawSNR_dB'});
disp(results3);
writetable(results3,fullfile(out,'Test3_NoiseResults.csv'));

fig5 = lightFigure([100 100 1100 430]);
layout5 = tiledlayout(fig5,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout5);
semilogy(ax,gammas,expectedPower,'k--',gammas,measuredPower,'ro','LineWidth',1.2);
xlabel(ax,'SNR setting \Gamma (dB)');
ylabel(ax,'Noise power per sample');
title(ax,'(a) Noise power','FontWeight','normal','Color','k');
lightLegend(legend(ax,{'Expected','Measured'},'Location','best'));
styleAxes(ax,11);

ax = nexttile(layout5);
plot(ax,gammas,expectedRaw,'k--',gammas,measuredRaw,'bo','LineWidth',1.2);
xlabel(ax,'SNR setting \Gamma (dB)');
ylabel(ax,'Raw time-domain SNR (dB)');
title(ax,'(b) Raw SNR of the composite return','FontWeight','normal','Color','k');
lightLegend(legend(ax,{'Expected','Measured'},'Location','best'));
styleAxes(ax,11);

exportgraphics(fig5,fullfile(out,'Verification_Test3_NoisePower.png'), ...
    'BackgroundColor','white','Resolution',300);
close(fig5);

fprintf('\nAll results saved in:\n%s\n',out);
end

%% ================= helpers =================

function d = runV11(scriptFile,v0In,R0In,ampIn,snrIn)
v0 = v0In;
SimulatorR0_m = R0In;
AmplitudeMode = ampIn;
SimulatorVerbose = false;
SimulatorMakePlots = false;
CPI_s = 0.025;
SNR_dB = 40;
if ~isempty(snrIn)
    SNR_dB = snrIn;
end

figuresBefore = findall(0,'Type','figure');
evalc('run(scriptFile)');          % suppresses v11's console printing
newFigures = setdiff(findall(0,'Type','figure'),figuresBefore);
if ~isempty(newFigures)
    close(newFigures);
end

d.v0 = v0;
d.fs = fs;
d.lambda = lamda;
d.t = t;
d.buffer = bufferTime_s;
d.duration = SimulationTime_s;
d.Nw = NumSamplesPerFrame;
d.height = Height_m;
d.R0 = R0;
d.radar = RadarCord_ZX;
d.regime = gaitCharacteristic;
d.summary = AmplitudeSummary;

d.clean = [RdrRxSignal_Torso(:).'; RdrRxSignal_Head(:).'; ...
           RdrRxSignal_LUpperArm(:).'; RdrRxSignal_LForeArm(:).'; ...
           RdrRxSignal_RUpperArm(:).'; RdrRxSignal_RForeArm(:).'; ...
           RdrRxSignal_Thigh(:).'; RdrRxSignal_Shank(:).'; RdrRxSignal_Foot(:).'; ...
           RdrRxSignal_RThigh(:).'; RdrRxSignal_RShank(:).'; RdrRxSignal_RFoot(:).'];

d.torsoX = Torso_X(:).';
d.torsoZ = Torso_Z(:).';
d.torsoR = Torso_R(:).';
d.torsoSignal = RdrRxSignal_Torso(:).';
d.composite = RdrRxSignal_clean(:).';
d.received = RdrRxSignal(:).';
d.builtInSNR = MeasuredSampleSNR_dB;
end

function fig = lightFigure(position)
fig = figure('Visible','off','Color','w','Position',position);
try
    fig.Theme = 'light';            % R2025a and later, ignored otherwise
catch
end
end

function styleAxes(ax,fontSize)
set(ax,'FontName','Times New Roman','FontSize',fontSize,'Color','w', ...
    'XColor','k','YColor','k','GridColor',[0.6 0.6 0.6],'Box','on');
grid(ax,'on');
try
    ax.Toolbar.Visible = 'off';     % keeps the toolbar out of exported images
catch
end
end

function lightLegend(lg)
set(lg,'Color','w','TextColor','k','EdgeColor',[0.6 0.6 0.6]);
end

function velocity = ridgeVelocity(S,inside,velocityAxis)
P = 20*log10(abs(S(:,inside))+eps);
step = mean(diff(velocityAxis));
velocity = zeros(1,size(P,2));
for n = 1:size(P,2)
    [~,peak] = max(P(:,n));
    delta = 0;
    if peak > 1 && peak < size(P,1)
        left = P(peak-1,n);
        centre = P(peak,n);
        right = P(peak+1,n);
        delta = 0.5*(left-right)/(left-2*centre+right);
    end
    velocity(n) = velocityAxis(peak)+delta*step;
end
end

function [x,z] = torsoReference(tw,v,H,Rstart)
relativeSpeed = v/((0.245+0.246)*H);
cycleTime = 1.346/sqrt(relativeSpeed);
stanceFraction = (0.752*cycleTime-0.143)/cycleTime;
phase = tw/cycleTime;
verticalAmplitude = 0.015*relativeSpeed;
if relativeSpeed > 0.5
    foreAftAmplitude = -0.021;
    leanAmplitude = 2;
else
    foreAftAmplitude = -0.084*relativeSpeed^2+0.084*relativeSpeed;
    leanAmplitude = -8*relativeSpeed^2+8*relativeSpeed;
end
lean = (-leanAmplitude+leanAmplitude*sin(4*pi*phase-0.2*pi))*pi/180;
x = Rstart-v*tw ...
    +foreAftAmplitude*sin(4*pi*phase+4*pi*(0.625-stanceFraction)) ...
    +(0.288*H/2)*sin(lean);
z = 0.58*H-verticalAmplitude+verticalAmplitude*sin(4*pi*phase-0.7*pi) ...
    +(0.288*H/2)*cos(lean);
end

function [velocity,regime] = chenReference(tw,v,H,Rstart,radarZX)
thighLength = 0.245*H;
shankLength = 0.246*H;
footLength = 0.143*H;
torsoLength = 0.288*H;
upperArmLength = 0.188*H;
forearmLength = 0.145*H;
headLength = 0.13*H;

relativeSpeed = v/(thighLength+shankLength);
if relativeSpeed < 0.5
    regime = 'a';
elseif relativeSpeed < 1.3
    regime = 'b';
else
    regime = 'c';
end

cycleTime = 1.346/sqrt(relativeSpeed);
stanceFraction = (0.752*cycleTime-0.143)/cycleTime;
phase = tw/cycleTime;
omega = 4*pi/cycleTime;

verticalAmplitude = 0.015*relativeSpeed;
verticalOffset = -verticalAmplitude+verticalAmplitude*sin(4*pi*phase-0.7*pi);
verticalVelocity = verticalAmplitude*omega*cos(4*pi*phase-0.7*pi);

if relativeSpeed > 0.5
    foreAftAmplitude = -0.021;
    leanAmplitude = 2;
else
    foreAftAmplitude = -0.084*relativeSpeed^2+0.084*relativeSpeed;
    leanAmplitude = -8*relativeSpeed^2+8*relativeSpeed;
end
foreAftPhase = 0.625-stanceFraction;

baseX = Rstart-v*tw+foreAftAmplitude*sin(4*pi*phase+4*pi*foreAftPhase);
baseVX = -v+foreAftAmplitude*omega*cos(4*pi*phase+4*pi*foreAftPhase);
lean = (-leanAmplitude+leanAmplitude*sin(4*pi*phase-0.2*pi))*pi/180;
leanRate = leanAmplitude*omega*cos(4*pi*phase-0.2*pi)*pi/180;
spineZ = 0.58*H+verticalOffset;
velocity = zeros(12,numel(tw));

reach = torsoLength/2;
velocity(1,:) = projectRadial( ...
    baseX+reach*sin(lean),spineZ+reach*cos(lean), ...
    baseVX+reach*cos(lean).*leanRate, ...
    verticalVelocity-reach*sin(lean).*leanRate,radarZX);
reach = torsoLength+headLength;
velocity(2,:) = projectRadial( ...
    baseX+reach*sin(lean),spineZ+reach*cos(lean), ...
    baseVX+reach*cos(lean).*leanRate, ...
    verticalVelocity-reach*sin(lean).*leanRate,radarZX);

shoulderX = baseX+torsoLength*sin(lean);
shoulderZ = spineZ+torsoLength*cos(lean);
shoulderVX = baseVX+torsoLength*cos(lean).*leanRate;
shoulderVZ = verticalVelocity-torsoLength*sin(lean).*leanRate;

swingAmplitude = 9.88*relativeSpeed;
for side = 1:2
    if side == 1
        armPhase = phase+0.5;
    else
        armPhase = phase;
    end
    shoulderAngle = (3-swingAmplitude*cos(2*pi*armPhase))*pi/180;
    shoulderRate = swingAmplitude*2*pi*sin(2*pi*armPhase)/cycleTime*pi/180;
    upperAngle = shoulderAngle-lean;
    upperRate = shoulderRate-leanRate;

    [elbow,elbowRate] = angleAndRate('elbow',regime,relativeSpeed, ...
        stanceFraction,mod(armPhase,1),cycleTime);
    foreAngle = upperAngle+elbow;
    foreRate = upperRate+elbowRate;

    elbowX = shoulderX-upperArmLength*sin(upperAngle);
    elbowZ = shoulderZ-upperArmLength*cos(upperAngle);
    elbowVX = shoulderVX-upperArmLength*cos(upperAngle).*upperRate;
    elbowVZ = shoulderVZ+upperArmLength*sin(upperAngle).*upperRate;

    wristX = elbowX-forearmLength*sin(foreAngle);
    wristZ = elbowZ-forearmLength*cos(foreAngle);
    wristVX = elbowVX-forearmLength*cos(foreAngle).*foreRate;
    wristVZ = elbowVZ+forearmLength*sin(foreAngle).*foreRate;

    upperRow = 1+2*side;
    velocity(upperRow,:) = projectRadial( ...
        (shoulderX+elbowX)/2,(shoulderZ+elbowZ)/2, ...
        (shoulderVX+elbowVX)/2,(shoulderVZ+elbowVZ)/2,radarZX);
    velocity(upperRow+1,:) = projectRadial( ...
        (elbowX+wristX)/2,(elbowZ+wristZ)/2, ...
        (elbowVX+wristVX)/2,(elbowVZ+wristVZ)/2,radarZX);
end

hipZ = 0.53*H+verticalOffset;
for side = 1:2
    limbPhase = mod(phase+(side-1)/2,1);

    [hip,hipRate] = angleAndRate('hip',regime,relativeSpeed, ...
        stanceFraction,limbPhase,cycleTime);
    [knee,kneeRate] = angleAndRate('knee',regime,relativeSpeed, ...
        stanceFraction,limbPhase,cycleTime);
    [ankle,ankleRate] = angleAndRate('ankle',regime,relativeSpeed, ...
        stanceFraction,limbPhase,cycleTime);

    kneeX = baseX-thighLength*sin(hip);
    kneeZ = hipZ-thighLength*cos(hip);
    kneeVX = baseVX-thighLength*cos(hip).*hipRate;
    kneeVZ = verticalVelocity+thighLength*sin(hip).*hipRate;

    shankAngle = hip-knee;
    shankRate = hipRate-kneeRate;

    ankleX = kneeX-shankLength*sin(shankAngle);
    ankleZ = kneeZ-shankLength*cos(shankAngle);
    ankleVX = kneeVX-shankLength*cos(shankAngle).*shankRate;
    ankleVZ = kneeVZ+shankLength*sin(shankAngle).*shankRate;

    footAngle = hip-knee+ankle;
    footRate = hipRate-kneeRate+ankleRate;

    firstRow = 4+3*side;
    velocity(firstRow,:) = projectRadial( ...
        baseX-thighLength/2*sin(hip),hipZ-thighLength/2*cos(hip), ...
        baseVX-thighLength/2*cos(hip).*hipRate, ...
        verticalVelocity+thighLength/2*sin(hip).*hipRate,radarZX);
    velocity(firstRow+1,:) = projectRadial( ...
        kneeX-shankLength/2*sin(shankAngle), ...
        kneeZ-shankLength/2*cos(shankAngle), ...
        kneeVX-shankLength/2*cos(shankAngle).*shankRate, ...
        kneeVZ+shankLength/2*sin(shankAngle).*shankRate,radarZX);
    velocity(firstRow+2,:) = projectRadial( ...
        ankleX-footLength/2*cos(footAngle), ...
        ankleZ+footLength/2*sin(footAngle), ...
        ankleVX+footLength/2*sin(footAngle).*footRate, ...
        ankleVZ+footLength/2*cos(footAngle).*footRate,radarZX);
end
end

function [angle,rate] = angleAndRate(joint,regime,relativeSpeed, ...
        stanceFraction,phase,cycleTime)
[knots,values] = controlPoints(joint,regime,relativeSpeed,stanceFraction);
polynomial = pchip(knots,values);
[breaks,coefficients,~,order] = unmkpp(polynomial);
derivativeCoefficients = coefficients(:,1:end-1).*(order-1:-1:1);
derivative = mkpp(breaks,derivativeCoefficients);
angle = ppval(polynomial,phase)*pi/180;
rate = ppval(derivative,phase)*pi/(180*cycleTime);
end

function [knots,values] = controlPoints(joint,regime,Vr,stance)
s = (Vr-1.3)/1.7;
switch joint
    case 'hip'
        switch regime
            case 'a'
                x = [-0.1 0.5 0.9]; y = [50*Vr -30*Vr 50*Vr];
            case 'b'
                x = [-0.1 0.5 0.9]; y = [25 -15 25];
            otherwise
                x = [0.2*s-0.1 0.5 0.9]; y = [5*s+25 -15 6*s+25];
        end
        if abs((x(1)+1)-x(3)) < 1e-12
            knots = [x(1)-1 x(2)-1 x(1) x(2) x(3) x(2)+1 x(3)+1];
            values = [y(1) y(2) y(1) y(2) y(3) y(2) y(3)];
            return
        end
    case 'knee'
        switch regime
            case 'a'
                x = [0.17 0.4 0.75 1]; y = [3 3 140*Vr 3];
            case 'b'
                x = [0.17 0.4 0.75 1]; y = [3 3 70 3];
            otherwise
                x = [-0.05*s+0.17 0.4 -0.05*s+0.75 -0.03*s+1];
                y = [22*s+3 3 -5*s+70 3*s+3];
        end
    case 'ankle'
        switch regime
            case 'a'
                x = [0 0.08 0.5 stance 0.85];
                y = [-3 -30*Vr-3 22*Vr-3 -34*Vr-3 -3];
            case 'b'
                x = [0 0.08 0.5 stance 0.85];
                y = [-3 -18 8 -20 -3];
            otherwise
                x = [0 0.08 -0.1*s+0.5 stance 0.85];
                y = [5*s-3 4*s-18 -3*s+8 -8*s-20 5*s-3];
        end
    case 'elbow'
        r = (Vr-0.5)/0.8;
        switch regime
            case 'a'
                x = [0.05 0.5 0.9]; y = [6*Vr+3 34*Vr+3 10*Vr+3];
            case 'b'
                x = [0.05 0.01*r+0.5 0.9]; y = [8*r+6 24*r+20 9*r+8];
            otherwise
                x = [0.05 0.04*s+0.51 -0.1*s+0.9];
                y = [-6*s+14 26*s+44 -6*s+17];
        end
end
knots = [x-1 x x+1];
values = repmat(y,1,3);
end

function vr = projectRadial(x,z,vx,vz,radarZX)
dx = x-radarZX(2);
dz = z-radarZX(1);
vr = -(dx.*vx+dz.*vz)./hypot(dx,dz);
end
