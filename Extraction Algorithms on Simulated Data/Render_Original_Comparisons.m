%% Simulated algorithm comparisons
% Zahraa Sattar | EEE4022S
%
% Redraw the report comparisons from the original saved sweep results.
% The CSV files retain the original noise draws and selected settings.

%% Light mode, so the figures are never black even with MATLAB's dark theme on
try, settings().matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light'; catch, end
set(groot,'defaultFigureColor','w','defaultAxesColor','w','defaultTextColor','k', ...
    'defaultAxesXColor','k','defaultAxesYColor','k','defaultLegendTextColor','k','defaultLegendColor','w');

%% Where the results are and where the figures go
here = fullfile(fileparts(mfilename('fullpath')), 'Report Results');
E1  = readCsv(fullfile(here,'E1_overlap.csv'));         % published -> adapted, same noise
FIX = readCsv(fullfile(here,'Fix_Validation.csv'));    % original band against +-25% band
S4  = readCsv(fullfile(here,'Stage4_validation.csv')); % final validation, estimated Tc
outDir = fullfile(here,'Comparison Figures');
mkdir(outDir);

%% Sizes and colours
fontAx  = 12;    % tick labels
fontLab = 13;    % axis labels and panel titles
grey  = [0.60 0.63 0.67];   % published
amber = [0.89 0.66 0.34];   % + 75% overlap
mint  = [0.44 0.72 0.66];   % + half-bin drop
teal  = [0.06 0.44 0.39];   % + centering and gap recovery
blue  = [0.12 0.31 0.61];   % fully adapted, with the +-25% band
stageCols = [grey; amber; mint; teal];
stages = {'ov50_c1','main'; 'ov75_c1','main'; 'ov75_c0.5','main'; 'ov75_c0.5','v5'};
names  = {'Published (50% overlap, one-bin drop)','+ 75% overlap','+ half-bin drop', ...
          '+ centering and gap recovery','+ \pm25% band (fully adapted)'};

%% Fig 4.5: published algorithm against mine, one change at a time
for snr = [60 30]
    fig = figure('Units','centimeters','Position',[2 2 24 15],'Color','w','Visible','off');

    axes('Position',[0.065 0.57 0.27 0.34]); hold on;
    b = bar(1:10, [pick(FIX,'Setting','base_0.372',snr,'WS_RMSE_pct'), pick(FIX,'Setting','band25_0.372',snr,'WS_RMSE_pct')], 0.85, 'grouped', 'EdgeColor','none');
    b(1).FaceColor = grey; b(2).FaceColor = blue;
    yline5(); styleAxes('RMSE (%)','(a) Walking speed',fontAx,fontLab);
    speedHandles = b;

    axes('Position',[0.415 0.57 0.27 0.34]); hold on;
    colourBars(bar(1:10, stageMatrix(E1,stages,'StepTime',snr,'PerfectTrials_pct'), 0.85, 'grouped', 'EdgeColor','none'), stageCols);
    ylim([0 105]); styleAxes('Exact trials (%)','(b) Step detection: exact trials',fontAx,fontLab);

    axes('Position',[0.065 0.12 0.27 0.34]); hold on;
    colourBars(bar(1:10, stageMatrix(E1,stages,'StepTime',snr,'RMSE_pct'), 0.85, 'grouped', 'EdgeColor','none'), stageCols);
    yline5(); styleAxes('RMSE (%)','(c) Step time',fontAx,fontLab);
    xlabel('Walking speed (m/s)','FontSize',fontLab);

    axes('Position',[0.415 0.12 0.27 0.34]); hold on;
    M = stageMatrix(E1,stages,'StrideDuration',snr,'RMSE_pct');
    bb = bar(1:10, M, 0.85, 'grouped', 'EdgeColor','none');
    colourBars(bb, stageCols);
    set(gca,'YScale','log','YTick',[0.1 1 10 100],'YTickLabel',{'0.1','1','10','100'}); ylim([0.03 100]);
    yline5(); styleAxes('RMSE (%)','(d) Stride duration',fontAx,fontLab);
    xlabel('Walking speed (m/s)','FontSize',fontLab);
    lgSpeed = legend(speedHandles, {'Original band','\pm25% band (adapted)'}, ...
        'FontSize',11,'Box','off');
    lgSpeed.Title.String = 'Walking speed (a)';
    lgSpeed.Units = 'normalized'; lgSpeed.Position = [0.72 0.72 0.27 0.16];
    keyNames = {'Published: 50% overlap,\newline one-bin drop', ...
        '+ 75% overlap','+ half-bin drop', ...
        '+ centering and\newline gap recovery'};
    lg = legend(bb, keyNames, 'FontSize',11,'Box','off','NumColumns',1);
    lg.Title.String = 'Processing stages (b–d)';
    lg.Units = 'normalized'; lg.Position = [0.72 0.22 0.27 0.43];
    annotation(fig,'textbox',[0.73 0.10 0.25 0.07], ...
        'String','Dashed line: 5% limit','Color','k','EdgeColor','none','FontSize',11);
    savePair(fig, outDir, sprintf('Fig_PublishedVsAdapted_%ddB',snr));
end

S1 = readCsv(fullfile(here,'Stage1_CPI_sweep.csv'));
S2 = readCsv(fullfile(here,'Stage2_SNR_sweep.csv'));
lw=1.2; ms=4;

%% Figure 3: effect of the window (30 dB) and of the SNR (chosen window)
show = [0.4 0.6 1.0 1.4 2.0];                       % five speeds keep the plot readable
mk   = {'o-','s-','^-','d-','v-'};
fig = figure('Units','centimeters','Position',[2 2 23 15],'Color','w','Visible','off');
positions = [0.07 0.58 0.29 0.32;0.44 0.58 0.29 0.32;0.07 0.12 0.29 0.32;0.44 0.12 0.29 0.32];
ax = gobjects(1,4);
for p=1:4
    ax(p)=axes('Position',positions(p,:));
    colororder(ax(p),[0.12 0.31 0.61;0.75 0.35 0.08;0.12 0.52 0.40;0.48 0.25 0.60;0.30 0.30 0.30]);
end
for k = 1:numel(show)
    v = show(k);
    axes(ax(1)); hold on;   % walking speed against window at 30 dB
    [x,y] = curve(S1,'WalkingSpeed',v,'CPI_speed_act','RMSE_pct');
    semilogy(x,y,mk{k},'LineWidth',lw,'MarkerSize',ms);
    axes(ax(2)); hold on;   % exact step trials against window at 30 dB
    [x,y] = curve(S1,'StepTime',v,'CPI_step_act','PerfectTrials_pct');
    plot(x,y,mk{k},'LineWidth',lw,'MarkerSize',ms);
    axes(ax(3)); hold on;   % walking speed against SNR at the chosen window
    [x,y] = curve(S2,'WalkingSpeed',v,'SNR_dB','RMSE_pct');
    semilogy(x,y,mk{k},'LineWidth',lw,'MarkerSize',ms);
    axes(ax(4)); hold on;   % exact step trials against SNR at the chosen window
    [x,y] = curve(S2,'StepTime',v,'SNR_dB','PerfectTrials_pct');
    plot(x,y,mk{k},'LineWidth',lw,'MarkerSize',ms);
end
winTicks = [0.046 0.093 0.186 0.372 0.743];
for p = 1:2
    axes(ax(p)); set(gca,'XScale','log','XTick',winTicks,'XTickLabel',{'0.046','0.093','0.186','0.372','0.743'});
    xlim([0.04 0.85]); xlabel('Window (s)','FontSize',fontLab);
end
for p = 3:4
    axes(ax(p)); xlim([15 95]); xlabel('SNR \Gamma (dB)','FontSize',fontLab);
end
for p = [1 3]
    axes(ax(p)); set(gca,'YScale','log','YTick',[0.1 1 10 100],'YTickLabel',{'0.1','1','10','100'}); ylim([0.03 150]);
    yline5();
end
for p = [2 4]
    axes(ax(p)); ylim([-5 105]);
end
axes(ax(1)); styleSweep('Speed RMSE (%)','(a) Speed, 30dB',fontAx,fontLab);
axes(ax(2)); styleSweep('Exact trials (%)','(b) Steps, 30dB',fontAx,fontLab);
axes(ax(3)); styleSweep('Speed RMSE (%)','(c) Speed, chosen window',fontAx,fontLab);
axes(ax(4)); styleSweep('Exact trials (%)','(d) Steps, chosen window',fontAx,fontLab);
lg = legend(arrayfun(@(s) sprintf('%.1fm/s',s), show,'UniformOutput',false), ...
       'FontSize',fontAx,'Box','off');
lg.Units='normalized'; lg.Position=[0.77 0.36 0.21 0.32];
lg.Title.String='Walking speed';
savePair(fig, outDir, 'Fig_CPI_SNR_Effect');

%% Accuracy and repeatability of the final setting (bars = RMSE, error bars = trial SD)
params = {'WalkingSpeed','(a) Walking speed'; 'StepTime','(b) Step time'; ...
          'StrideDuration','(c) Stride duration'; 'StrideLength','(d) Stride length'};
cols = [0 0.45 0.70; 0.84 0.37 0; 0 0.62 0.45; 0.48 0.18 0.56];
fig = figure('Units','centimeters','Position',[2 2 26 7.5],'Color','w','Visible','off');
for p = 1:4
    axes('Position',[0.055+(p-1)*0.245 0.23 0.205 0.62]); hold on;
    r  = pick(S4,'Parameter',params{p,1},[],'RMSE_pct');      % error over all trials
    sd = pick(S4,'Parameter',params{p,1},[],'TrialStd_pct');  % spread from trial to trial (repeatability)
    bar(1:10, r, 0.65, 'FaceColor', cols(p,:), 'EdgeColor','none');
    errorbar(1:10, r, sd, 'k', 'LineStyle','none', 'LineWidth', 0.8, 'CapSize', 3);
    yline5(); ylim([0 10]);
    styleAxes('', params{p,2}, fontAx, fontLab);
    xlabel('Speed (m/s)','FontSize',fontLab);
    if p == 1, ylabel('RMSE (%)','FontSize',fontLab); end
end
savePair(fig, outDir, 'Fig_Accuracy_Repeatability');
fprintf('Figures saved in %s\n', outDir);

%% ---------------- helper functions ----------------
function T = readCsv(fileName)
    txt   = fileread(fileName);
    lines = regexp(txt,'\r?\n','split');
    lines = lines(~cellfun('isempty',lines));
    head  = strsplit(lines{1},',');
    rows  = cellfun(@(s) strsplit(s,','), lines(2:end),'UniformOutput',false);
    rows  = vertcat(rows{:});
    for c = 1:numel(head)
        name = matlab.lang.makeValidName(strrep(head{c},'"',''));
        col  = strrep(rows(:,c),'"','');
        num  = str2double(col);
        if all(isnan(num) == strcmpi(col,'NaN'))   % a numeric column (NaN only where written NaN)
            T.(name) = num;
        else
            T.(name) = col;
        end
    end
end

function y = pick(T, keyName, keyValue, snr, column)
    keep = strcmp(T.(keyName), keyValue);
    if ~isempty(snr), keep = keep & T.SNR_dB == snr; end
    [~,order] = sort(T.Speed_mps(keep));
    vals = T.(column)(keep);
    y = vals(order);
end

function M = stageMatrix(T, stages, parameter, snr, column)
    M = zeros(10, size(stages,1));
    for k = 1:size(stages,1)
        keep = strcmp(T.Variant,stages{k,1}) & strcmp(T.Extractor,stages{k,2}) & ...
               strcmp(T.Parameter,parameter) & T.SNR_dB == snr;
        [~,order] = sort(T.Speed_mps(keep));
        vals = T.(column)(keep);
        M(:,k) = vals(order);
    end
end

function [x,y] = curve(T, parameter, speed, xColumn, yColumn)
    keep = strcmp(T.Parameter,parameter) & abs(T.Speed_mps - speed) < 1e-6;
    [x,order] = sort(T.(xColumn)(keep));
    vals = T.(yColumn)(keep);
    y = vals(order);
end

function yline5()
    xl = xlim; plot(xl, [5 5], 'k--', 'LineWidth', 1); xlim(xl);
end

function styleAxes(yLabelText, titleText, fontAx, fontLab)
    grid on; box off;                                      % slim look: no box, light grid
    set(gca,'FontSize',fontAx,'LineWidth',0.6,'Color','w','XColor','k','YColor','k','GridColor',[0.85 0.85 0.85],'GridAlpha',1);
    if ~isempty(yLabelText), ylabel(yLabelText,'FontSize',fontLab); end
    title(titleText,'FontSize',fontLab,'FontWeight','normal');
    set(gca,'XTick',1:2:10,'XTickLabel',{'0.4','0.6','1','1.4','1.8'});   % every second speed, so the labels do not overlap
    xlim([0.4 10.6]);
end

function savePair(fig, outDir, name)
    savefig(fig,fullfile(outDir,[name '.fig']));
    axesHandles = findall(fig,'Type','axes');
    for n=1:numel(axesHandles), axesHandles(n).Toolbar.Visible='off'; end
    drawnow;
    set(fig,'PaperUnits','centimeters');
    pos = get(fig,'Position');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)],'PaperSize',[pos(3) pos(4)]);
    print(fig, fullfile(outDir,[name '.png']), '-dpng', '-r300');
    print(fig, fullfile(outDir,[name '.pdf']), '-dpdf');
end

function colourBars(b, cols)
    for k = 1:numel(b), b(k).FaceColor = cols(k,:); end
end

function styleSweep(yLabelText,titleText,fontAx,fontLab)
    grid on; box off;
    set(gca,'FontSize',fontAx,'LineWidth',0.6,'Color','w','XColor','k','YColor','k');
    ylabel(yLabelText,'FontSize',fontLab);
    title(titleText,'FontSize',fontLab,'FontWeight','normal');
end
